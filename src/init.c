/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's initialization routines
 *
 *  Copyright (C) 2010--2015 jeremy theler
 *
 *  This file is part of milonga.
 *
 *  milonga is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  milonga is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with wasora.  If not, see <http://www.gnu.org/licenses/>.
 *------------------- ------------  ----    --------  --     -       -         -
 */
#include <stdlib.h>
#include <signal.h>

#include "milonga.h"

// maneje por el cambio de PetscOptionsHasName en 3.7.0
#if PETSC_VERSION_LT(3,7,0)
 #define PetscOptionsHasNameWrapper(a, b, c) PetscOptionsHasName(a, b, c)
#else
 #define PetscOptionsHasNameWrapper(a, b, c) PetscOptionsHasName(PETSC_NULL, a, b, c)
#endif

// inicializa petsc y slepc, define variables internas y pone el esquema por default
#undef  __FUNCT__
#define __FUNCT__ "plugin_init_before_parser"
int plugin_init_before_parser(void) {

  PetscBool flag;
  char *dummy;
  int i;
  
  if (sizeof(PetscReal) != sizeof(double)) {
    wasora_push_error_message("PETSc should be compiled with double-precision real scalar types");
    return WASORA_PARSER_ERROR;
  }
  
  // amasamos la linea de comandos original (porque la que saca getopt puede tener un orden que no nos sirve)
  // el chiste es que por ejemplo "-log_summary" es atrapado por el getopt de wasora como "-l"
  // hay que re-escrbir eso como "--slepc_opt log_summary"
  // si alguna opcion tiene argumento hay que ponerlo como "--slepc_opt pc_type=sor"
  for (i = 0; i < wasora.argc_orig; i++) {
    if (strcmp(wasora.argv_orig[i], "--slepc_opt") == 0 || strcmp(wasora.argv_orig[i], "--petsc_opt") == 0) {
      if (i >= (wasora.argc_orig-1)) {
        wasora_push_error_message("commandline option --slepc_opt needs an argument");
        return WASORA_PARSER_ERROR;
      } else if (wasora.argv_orig[i+1][0] == '-') {
        wasora_push_error_message("the argument of commandline option --slepc_opt should not start with a dash (it is added automatically)");
        return WASORA_PARSER_ERROR;
      }
      
      if ((dummy = strchr(wasora.argv_orig[i+1], '=')) != NULL)  {
        char *tmp1, *tmp2;
        *dummy = '\0';
        tmp1 = strdup(wasora.argv_orig[i+1]);
        tmp2 = strdup(dummy+1);
        wasora.argv_orig[i]   = realloc(wasora.argv_orig[i],   strlen(wasora.argv_orig[i+1])+2);
        wasora.argv_orig[i+1] = realloc(wasora.argv_orig[i+1], strlen(dummy)+1);
        sprintf(wasora.argv_orig[i],  "-%s", tmp1);
        sprintf(wasora.argv_orig[i+1], "%s", tmp2);
        free(tmp1);
        free(tmp2);
        
      } else {
        char *tmp1;
        tmp1 = strdup(wasora.argv_orig[i+1]);
        wasora.argv_orig[i+1] = realloc(wasora.argv_orig[i+1], strlen(tmp1)+1);
        wasora.argv_orig[i][0] = '\0';
        sprintf(wasora.argv_orig[i+1],  "-%s", tmp1);
        free(tmp1);
      }
      i++;
    }
  }
  
  // inicializamos la slepc (que a su vez inicializa la petsc)
  // le pasamos la linea de comandos que acabamos de amasar
  petsc_call(SlepcInitialize(&wasora.argc_orig, &wasora.argv_orig, (char*)0, PETSC_NULL));
  // los segfaults son segfaults, no queremos que la petsc meta las narices
  signal(SIGSEGV, SIG_DFL);

  // esto lo vamos a usar despues cuando hagamos el chiste en paralelo
  MPI_Comm_rank(PETSC_COMM_WORLD, &milonga.rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &milonga.size);

  // instalamos nuestro error handler para errores la petsc 
  petsc_call(PetscPushErrorHandler(&milonga_handler, NULL));
  
  // registramos nuestros eventos
  petsc_call(PetscClassIdRegister("milonga", &milonga.petsc_classid));

  petsc_call(PetscLogStageRegister("Initialization", &milonga.petsc_stage_init));
  petsc_call(PetscLogStageRegister("Assembly", &milonga.petsc_stage_build));
  petsc_call(PetscLogStageRegister("Solution", &milonga.petsc_stage_solve));
  
  petsc_call(PetscLogEventRegister("milonga_init", milonga.petsc_classid, &milonga.petsc_event_init));
  petsc_call(PetscLogEventRegister("milonga_build", milonga.petsc_classid, &milonga.petsc_event_build));
  petsc_call(PetscLogEventRegister("milonga_solve", milonga.petsc_classid, &milonga.petsc_event_solve));
  

  // inicializamos mesh
  if (!wasora_mesh.initialized) {
    wasora_call(wasora_mesh_init_before_parser());
  }

  // variables internas
///va+keff+name keff
///va+keff+desc Efective multiplication factor as computed by solving the
///va+keff+desc multigroup neutron diffusion equation. It is equal to 1.0 until
///va+keff+desc `MILONGA_STEP` is executed.  
  milonga.vars.keff = wasora_define_variable("keff");
  wasora_var(milonga.vars.keff) = 1.0;
  
///va+power+name power
///va+power+desc Power setpoint used to normalize the computed fluxes.
///va+power+desc By default it is zero, meaning the fluxes are normalized so
///va+power+desc the mean value is equal to one. If this variable is set,
///va+power+desc then the property `eSigmaF` should be nonzero through the domain.
  milonga.vars.power = wasora_define_variable("power");

  
///va+sn_a_weight+name sn_alpha
///va+sn_a_weight+desc Factor used to control the upwinding of the discrete ordinates~$S_N$ formulation.
///va+sn_a_weight+desc For the finite volumes scheme, the difference factor $a = 1/2 (1+\alpha)$ is used to weight the adjacent cell fluxes,
///va+sn_a_weight+desc i.e. $\alpha = 0$ ($a=1/2$) corresponds to diamond difference
///va+sn_a_weight+desc  and $\alpha = 1$ ($a=1$) corresponds to full upwinding.
///va+sn_a_weight+desc For the finite elements scheme, a Streamline-Upwind Petrov-Galerkin stabilization term
///va+sn_a_weight+desc $1/2 \alpha \ell/\| \vec{\Omega} \cdot ( \vec{\Omega} \cdot \nabla h)\|$ is added to the weighting functions.
///va+sn_a_weight+desc A value of~$\alpha=0$ corresponds to the unstable original Galerkin scheme and
///va+sn_a_weight+desc            $\alpha=1$ to the cotangent-optimum upwind term.
///va+sn_a_weight+desc For~$\alpha < 1$ expect oscillations in the fluxes and even non-convergence for small values of~$\alpha$.
///va+sn_a_weight+desc Default is $\alpha=0.5$.  
  milonga.vars.sn_alpha = wasora_define_variable("sn_alpha");
  wasora_var(milonga.vars.sn_alpha) = 0.5;

///va+unknowns+name unknowns
///va+unknowns+desc Number of total unknowns (size) of the problem. It is equal
///va+unknowns+desc to the number of spatial unknowns times the number of energy groups.  
  milonga.vars.unknowns = wasora_define_variable("unknowns");

  
///va+eigen_rel_tolerance+name eigen_rel_tolerance
///va+eigen_rel_tolerance+desc Tolerance passed the iterative eigensolver, relative to the matrices norm.
///va+eigen_rel_tolerance+desc By default it is zero, meaning the default value of the SLEPc library.
  milonga.vars.rel_tolerance = wasora_define_variable("eigen_rel_tolerance");

  ///va+eigen_residual_norm+name eigen_residual_norm
///va+eigen_residual_norm+desc Residual norm $|R \phi - \lambda F \phi|_2$ (`|Rx-lFx|_2`)
///va+eigen_residual_norm+desc reached by the eigensolver.
  milonga.vars.residual_norm = wasora_define_variable("eigen_residual_norm");

///va+eigen_rel_error+name eigen_rel_error
///va+eigen_rel_error+desc Relative error $|R \phi - \lambda F \phi|_2/|\lambda \phi|_2$ `|Rx-lFx|_2/|kx|_2`
///va+eigen_rel_error+desc reached by the eigensolver.
  milonga.vars.rel_error = wasora_define_variable("eigen_rel_error");

///va+eigen_error_estimate+name eigen_error_estimate
///va+eigen_error_estimate+desc Absolute error estimate $|\lambda - \lambda_\text{real}|$ (`|lambda - lambda_real|`)
///va+eigen_error_estimate+desc of the obtained eigenvalue.
  milonga.vars.error_estimate = wasora_define_variable("eigen_error_estimate");

///va+time_wall_ini+name time_wall_ini
///va+time_wall_ini+desc Wall time insumed to initialize the problem, in seconds.
  milonga.vars.time_wall_ini  = wasora_define_variable("time_wall_ini");

///va+time_wall_build+name time_wall_build
///va+time_wall_build+desc Wall time insumed to build the problem matrices, in seconds.
  milonga.vars.time_wall_build = wasora_define_variable("time_wall_build");

///va+time_wall_solve+name time_wall_solve
///va+time_wall_solve+desc Wall time insumed to solve the eigen-problem, in seconds.
  milonga.vars.time_wall_solve = wasora_define_variable("time_wall_solve");

///va+time_wall_total+name time_wall_total
///va+time_wall_total+desc Wall time insumed to initialize, build and solve, in seconds.
  milonga.vars.time_wall_total = wasora_define_variable("time_wall_total");
  
///va+time_cpu_ini+name time_cpu_ini
///va+time_cpu_ini+desc CPU time insumed to initialize the problem, in seconds.
  milonga.vars.time_cpu_ini  = wasora_define_variable("time_cpu_ini");

///va+time_cpu_build+name time_cpu_build
///va+time_cpu_build+desc CPU time insumed to build the problem matrices, in seconds.
  milonga.vars.time_cpu_build = wasora_define_variable("time_cpu_build");

///va+time_cpu_solve+name time_cpu_solve
///va+time_cpu_solve+desc CPU time insumed to solve the eigen-problem, in seconds.
  milonga.vars.time_cpu_solve = wasora_define_variable("time_cpu_solve");

///va+time_wall_total+name time_cpu_total
///va+time_wall_total+desc CPU time insumed to initialize, build and solve, in seconds.
  milonga.vars.time_cpu_total = wasora_define_variable("time_cpu_total");
  
///va+time_petsc_ini+name time_petsc_ini
///va+time_petsc_ini+desc CPU time insumed by PETSc to initialize the problem, in seconds.
  milonga.vars.time_petsc_ini  = wasora_define_variable("time_petsc_ini");

///va+time_petsc_build+name time_petsc_build
///va+time_petsc_build+desc CPU time insumed by PETSc to build the problem matrices, in seconds.
  milonga.vars.time_petsc_build = wasora_define_variable("time_petsc_build");

///va+time_petsc_solve+name time_petsc_solve
///va+time_petsc_solve+desc CPU time insumed by PETSc to solve the eigen-problem, in seconds.
  milonga.vars.time_petsc_solve = wasora_define_variable("time_petsc_solve");

///va+time_wall_total+name time_wall_total
///va+time_wall_total+desc CPU time insumed by PETSc to initialize, build and solve, in seconds.
  milonga.vars.time_petsc_total = wasora_define_variable("time_petsc_total");

  ///va+petsc_flops+name petsc_flops
///va+petsc_flops+desc Number of floating point operations performed by PETSc/SLEPc.
  milonga.vars.flops_petsc = wasora_define_variable("flops_petsc");
         
///va+memory_use+name available_memory
///va+memory_use+desc Total available memory, in bytes.
  milonga.vars.available_memory = wasora_define_variable("available_memory");

///va+memory_usage_global+name global_memory_use
///va+memory_usage_global+desc Maximum resident set size (global memory used), in bytes.
  milonga.vars.memory_usage_global = wasora_define_variable("memory_usage_global");
  
///va+memory_usage_petsc+name petsc_memory_use
///va+memory_usage_petsc+desc Maximum resident set size (memory used by PETSc), in bytes.
  milonga.vars.memory_usage_petsc = wasora_define_variable("memory_usage_petsc");
  
  // por default ponemos un grupo de energia      
  milonga.groups = 1;
  
  // el chiste con la formulacion por la linea de comandos es que
  // hay que alocar cosas en tiempo de parseo dependiendo del valor de SN
  // asi que miramos esto aca e ignoramos keywords en parser.c si SN != 0
  // tomamos como caso particular diffusion como SN = 1
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--diffusion", &flag));
  if (flag) { milonga.formulation = formulation_diffusion; milonga.SN = 1; }
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--s2", &flag));
  if (flag) { milonga.formulation = formulation_sn; milonga.SN = 2; }
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--s4", &flag));
  if (flag) { milonga.formulation = formulation_sn; milonga.SN = 4; }
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--s6", &flag));
  if (flag) { milonga.formulation = formulation_sn; milonga.SN = 6; }
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--s8", &flag));
  if (flag) { milonga.formulation = formulation_sn; milonga.SN = 8; }

  // chequeos por si acaso  
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--s1", &flag));
  if (flag) { wasora_push_error_message("only s2, s4, s6 & s8 formulations are valid"); return WASORA_RUNTIME_ERROR; }
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--s3", &flag));
  if (flag) { wasora_push_error_message("only s2, s4, s6 & s8 formulations are valid"); return WASORA_RUNTIME_ERROR; }
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--s5", &flag));
  if (flag) { wasora_push_error_message("only s2, s4, s6 & s8 formulations are valid"); return WASORA_RUNTIME_ERROR; }
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--s7", &flag));
  if (flag) { wasora_push_error_message("only s2, s4, s6 & s8 formulations are valid"); return WASORA_RUNTIME_ERROR; }
  
  
  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "plugin_init_after_parser"
// verificamos la consistencia del input
int plugin_init_after_parser(void) {
  
  PetscBool flag;
  int g, g_prime;
  int n;
  function_t *function;
  material_t *material;
  mesh_post_t *mesh_post;
  mesh_post_dist_t *mesh_post_dist;
  mesh_fill_vector_t *mesh_fill_vector;
  mesh_find_max_t *mesh_find_max;
  mesh_integrate_t *mesh_integrate;  

  // las opciones en la linea de comando tienen precedencia sobre el input
  // por eso las miramos aca despues de leer el input
  
  // miramos si nos dijeron en la linea de comandos que esquema hay que usar
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--volumes", &flag));
  if (flag) milonga.scheme = scheme_volumes;
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--elements", &flag));
  if (flag) milonga.scheme = scheme_elements;

  // miramos si nos pidieron difusion, si pidieron sn ya lo hicimos en pre-parser
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--diffusion", &flag));
  if (flag) milonga.formulation = formulation_diffusion;

  // miramos si nos dicen que espectro tenemos que usar
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--largest", &flag));
  if (flag) milonga.spectrum = spectrum_largest_eigenvalue;
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--smallest", &flag));
  if (flag) milonga.spectrum = spectrum_smallest_eigenvalue;


  // en volumes somos cell centered (si no nos piedieron explicitamente otra cosa)
  if (milonga.scheme == scheme_volumes) {
    LL_FOREACH (wasora_mesh.posts, mesh_post) {
      LL_FOREACH(mesh_post->mesh_post_dists, mesh_post_dist) {
        if (mesh_post_dist->centering == centering_default) {
          mesh_post_dist->centering = centering_cells;
        }
      }
    }
    LL_FOREACH (wasora_mesh.fill_vectors, mesh_fill_vector) {
      if (mesh_fill_vector->centering == centering_default) {
        mesh_fill_vector->centering = centering_cells;
      }
    }
    LL_FOREACH (wasora_mesh.find_maxs, mesh_find_max) {
      if (mesh_find_max->centering == centering_default) {
        mesh_find_max->centering = centering_cells;
      }
    }
    LL_FOREACH (wasora_mesh.integrates, mesh_integrate) {
      if (mesh_integrate->centering == centering_default) {
        mesh_integrate->centering = centering_cells;
      }
    }
  }
  
  // resolvemos las XSs
  for (material = wasora_mesh.materials; material != NULL; material = material->hh.next) {
    xs_t *xs;

    xs = calloc(1, sizeof(xs_t));
    material->ext = (void *)xs;

    xs->D = calloc(milonga.groups, sizeof(expr_t *));
    xs->nuSigmaF = calloc(milonga.groups, sizeof(expr_t *));
    xs->SigmaT = calloc(milonga.groups, sizeof(expr_t *));
    xs->SigmaA = calloc(milonga.groups, sizeof(expr_t *));
    xs->eSigmaF = calloc(milonga.groups, sizeof(expr_t *));

    xs->SigmaS0 = calloc(milonga.groups, sizeof(expr_t **));
    xs->SigmaS1 = calloc(milonga.groups, sizeof(expr_t **));

    xs->S = calloc(milonga.groups, sizeof(expr_t *));
    
    for (g = 0; g < milonga.groups; g++) {
      milonga_resolve_xs_expr(material, "D",        &xs->D[g],        g, -1);
      milonga_resolve_xs_expr(material, "nuSigmaF", &xs->nuSigmaF[g], g, -1);
      milonga_resolve_xs_expr(material, "SigmaT",   &xs->SigmaT[g],   g, -1);
      milonga_resolve_xs_expr(material, "SigmaA",   &xs->SigmaA[g],   g, -1);
      milonga_resolve_xs_expr(material, "eSigmaF",  &xs->eSigmaF[g],  g, -1);

      milonga_resolve_xs_expr(material, "S",        &xs->S[g],        g, -1);

      xs->SigmaS0[g] = calloc(milonga.groups, sizeof(expr_t *));
      xs->SigmaS1[g] = calloc(milonga.groups, sizeof(expr_t *));
      for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
        milonga_resolve_xs_expr(material, "SigmaS",  &xs->SigmaS0[g][g_prime], g, g_prime);
        milonga_resolve_xs_expr(material, "SigmaSone",  &xs->SigmaS1[g][g_prime], g, g_prime);
      }
    }
  }
  
  // las funciones tipo propiedad van a parar a la malla del problema
  for (function = wasora.functions; function != NULL; function = function->hh.next) {
    if (function->type == type_pointwise_mesh_property) {
      function->mesh = milonga.mesh;
    }
  }
  
  // por default, el espectro de fision es 100% en el grupo rapido
  if ((milonga.vectors.chi = wasora_get_vector_ptr("chi")) == NULL) {
    milonga.vectors.chi = wasora_define_vector("chi", milonga.groups, NULL, NULL);
    wasora_call(wasora_vector_init(milonga.vectors.chi));
    gsl_vector_set(wasora_value_ptr(milonga.vectors.chi), 0, 1);
  } else if (milonga.vectors.chi->size != milonga.groups) {
    wasora_push_error_message("vector chi has size %d and problem has %d groups", milonga.vectors.chi->size, milonga.groups);
  }
  
  // procesamos los FLUX_POSTs
  LL_FOREACH (wasora_mesh.posts, mesh_post) {
    if (mesh_post->flags & POST_INCLUDE_FLUX) {
      
      // primero la potencia (si hay power)
      if (wasora_var(milonga.vars.power) != 0) {
        mesh_post_dist_t *mesh_post_dist = calloc(1, sizeof(mesh_post_dist_t));
        mesh_post_dist->scalar = milonga.functions.pow;
        LL_APPEND(mesh_post->mesh_post_dists, mesh_post_dist);
      }
      
      // despues los g flujos escalares
      for (g = 0; g < milonga.groups; g++) {
        mesh_post_dist_t *mesh_post_dist = calloc(1, sizeof(mesh_post_dist_t));
        mesh_post_dist->scalar = milonga.functions.phi[g];
        LL_APPEND(mesh_post->mesh_post_dists, mesh_post_dist);
      }
      
      // despues los flujos parciales
      if (milonga.SN > 1) {
        for (g = 0; g < milonga.groups; g++) {
          for (n = 0; n < milonga.directions; n++) {
            mesh_post_dist_t *mesh_post_dist = calloc(1, sizeof(mesh_post_dist_t));
            mesh_post_dist->scalar = milonga.functions.psi[n][g];
            LL_APPEND(mesh_post->mesh_post_dists, mesh_post_dist);
          }
        }
      }
    } 
    if (mesh_post->flags & POST_INCLUDE_XS) {
      property_data_t *property_data = NULL;
      function_t *function = NULL;
      
      for (property_data = wasora_mesh.materials->property_datums; property_data != NULL; property_data =  property_data->hh.next) {
        if ((function = wasora_get_function_ptr(property_data->property->name)) != NULL) {
          mesh_post_dist_t *mesh_post_dist = calloc(1, sizeof(mesh_post_dist_t));
          mesh_post_dist->scalar = function;
          mesh_post_dist->centering = centering_cells;       // las xs son siempre sobre celdas
          LL_APPEND(mesh_post->mesh_post_dists, mesh_post_dist);
        }
      }
    }
  }


  if (milonga.mesh != NULL) {
    wasora_call(milonga_set_entry_points());
  }
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "plugin_init_before_run"
int plugin_init_before_run(void) {
  wasora_var(milonga.vars.keff) = 1.0;

  milonga.problem_size = 0;
  milonga.spatial_unknowns = 0;

  if (milonga.problem_free != NULL) {
    wasora_call(milonga.problem_free());
  }
  wasora_call(milonga_free_global_objects());
  milonga.initialized = 0;

  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "plugin_finalize"
int plugin_finalize(void) {

  debug_t *debug;
  
  if (milonga.problem_free != NULL) {
    wasora_call(milonga.problem_free());
  }
  LL_FOREACH(milonga.debugs, debug) {
    PetscViewerDestroy(&debug->viewer);
  }
  
  wasora_call(milonga_free_global_objects());
  petsc_call(SlepcFinalize());  
  
  return WASORA_RUNTIME_OK;
}
