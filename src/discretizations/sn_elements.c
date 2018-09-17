/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's discrete ordinates with finite elements
 *
 *  Copyright (C) 2015--2016 jeremy theler
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <petscvec.h>
#include <petscmat.h>

#include "../milonga.h"
#include "sn_elements.h"

#define dof_index(n,g) ((n)*milonga.groups + (g))

#define BC_FACTOR 10

// nomenclatura como en la documentacion
int J;            // cantidad de nodos locales
int L;            // tamanio de la matriz elemental

// matrices de coeficientes
gsl_matrix *OMEGA;
gsl_matrix *A;
gsl_matrix *X;

// vector de fuentes
gsl_vector *S;

// matrices intermedias
gsl_matrix *P;   // esta es la H estabilizada con petrov-galerkin
gsl_matrix *OMEGAB;
gsl_matrix *AH;
gsl_matrix *XH;

// matriz elemental de rigidez
gsl_matrix *Ki;
// matriz elemental de scattering
gsl_matrix *Ai;
// matriz elemental de fision
gsl_matrix *Xi;

// vector elemental de fuentes
gsl_vector *Si;

// cosas para poner condiciones de contorno
int k_vacuum, current_size_vacuum, current_threshold_vacuum;
PetscInt *indexes_vacuum;

int k_mirror, current_size_mirror, current_threshold_mirror;
PetscInt *indexes_mirror;
PetscInt *indexes_mirror_reflected;


// caches para ensamblar mas rapido (en teoria)
double *Source;
double *SigmaA;
double *SigmaT;
double **SigmaS0;
double **SigmaS1;
double **chinuSigmaF;


#undef  __FUNCT__
#define __FUNCT__ "sn_elements_problem_init"
int sn_elements_problem_init(void) {

  int g, n;

  PetscFunctionBegin;
  milonga.spatial_unknowns = wasora_mesh.main_mesh->n_nodes;
  milonga_allocate_global_matrices(milonga.spatial_unknowns * milonga.directions * milonga.groups,
                                   wasora_mesh.main_mesh->max_first_neighbor_nodes * milonga.directions * milonga.groups,
                                   wasora_mesh.main_mesh->max_first_neighbor_nodes * milonga.directions * milonga.groups);
  milonga_allocate_global_vectors();
  

  // inicializamos los pesos de las ordenadas discretas
  wasora_call(sn_init_weights());
  
  wasora_var(wasora_mesh.vars.cells) = (double)wasora_mesh.main_mesh->n_cells;
  wasora_var(wasora_mesh.vars.nodes) = (double)wasora_mesh.main_mesh->n_nodes;
  wasora_var(wasora_mesh.vars.elements) = (double)wasora_mesh.main_mesh->n_elements;
  wasora_mesh.main_mesh->data_type = data_type_node;
  
  if (wasora_mesh.main_mesh->structured) {
    wasora_mesh_struct_init_rectangular_for_nodes(wasora_mesh.main_mesh);
  }
  
  for (n = 0; n < milonga.directions; n++) {
    for (g = 0; g < milonga.groups; g++) {
      wasora_call(sn_elements_results_fill_args(milonga.functions.psi[n][g]));
    }
  }
  
  for (g = 0; g < milonga.groups; g++) {
    wasora_call(sn_elements_results_fill_args(milonga.functions.phi[g]));
  }
  wasora_call(sn_elements_results_fill_args(milonga.functions.pow));

  wasora_call(mesh_node_indexes(wasora_mesh.main_mesh, milonga.groups * milonga.directions));

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


// esta rutina rellena datos administrativos de las funciones resultados
#undef  __FUNCT__
#define __FUNCT__ "sn_elements_results_fill_args"
int sn_elements_results_fill_args(function_t *function) {

  PetscFunctionBegin;
  
  // tenemos data
  function->data_size = milonga.spatial_unknowns;

  // pero tambien variables por si queremos hacer cuentitas
  function->var_argument = calloc(3, sizeof(var_t *));
  function->var_argument = wasora_mesh.vars.arr_x;

  function->data_argument = wasora_mesh.main_mesh->nodes_argument;
  function->data_value = calloc(function->data_size, sizeof(double));

  // y tipo milonga_status.mesh node en elementos
  function->type = type_pointwise_mesh_node;
  function->multidim_threshold = DEFAULT_MULTIDIM_INTERPOLATION_THRESHOLD;
  function->mesh = wasora_mesh.main_mesh;

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


#undef  __FUNCT__
#define __FUNCT__ "sn_elements_allocate_general_elemental_objects"
int sn_elements_allocate_general_elemental_objects(void) {
  
  int n, g, m;

  PetscFunctionBegin;
  
  OMEGA = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom, milonga.directions * milonga.groups * milonga.dimensions);
  
  for (m = 0; m < milonga.directions; m++) {
    for (g = 0; g < milonga.groups; g++) {
      for (n = 0; n < milonga.dimensions; n++) {
        gsl_matrix_set(OMEGA, dof_index(m,g), n*milonga.directions*milonga.groups + dof_index(m,g), Omega[m][n]);
      }
    }
  }
  
  A = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom, wasora_mesh.main_mesh->degrees_of_freedom);
  X = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom, wasora_mesh.main_mesh->degrees_of_freedom);
  S = gsl_vector_calloc(wasora_mesh.main_mesh->degrees_of_freedom);
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "sn_elements_allocate_particular_elemental_objects"
int sn_elements_allocate_particular_elemental_objects(element_t *element) {
  PetscFunctionBegin;
  
  J = element->type->nodes;
  L = milonga.directions * milonga.groups * element->type->nodes;

  // TODO: esta las tendria que alocar mesh
  gsl_matrix_free(wasora_mesh.main_mesh->fem.H);
  wasora_mesh.main_mesh->fem.H = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom, L);
  
  gsl_matrix_free(wasora_mesh.main_mesh->fem.B);
  wasora_mesh.main_mesh->fem.B = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom * wasora_mesh.main_mesh->bulk_dimensions, L);
  
  // esta no se, es la de petrov (h con un cacho de B)
  gsl_matrix_free(P);
  P = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom, L);

  gsl_matrix_free(OMEGAB);  
  OMEGAB = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom, L);
  gsl_matrix_free(AH);
  AH = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom, L);
  gsl_matrix_free(XH);
  XH = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom, L);

  gsl_matrix_free(Ki);
  Ki = gsl_matrix_calloc(L, L);
  gsl_matrix_free(Ai);
  Ai = gsl_matrix_calloc(L, L);
  gsl_matrix_free(Xi);
  Xi = gsl_matrix_calloc(L, L);
  gsl_vector_free(Si);
  Si = gsl_vector_calloc(L);

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "sn_elements_matrices_build"
int sn_elements_matrices_build(void) {
  
  int i;
  
  PetscFunctionBegin;

  wasora_call(sn_elements_allocate_general_elemental_objects());
  
  for (i = 0; i < wasora_mesh.main_mesh->n_elements; i++) {
    if (wasora_mesh.main_mesh->element[i].type->dim == wasora_mesh.main_mesh->bulk_dimensions) {

      // solo los elementos que tengan la dimension del problema
      // son los que usamos para las matrices elementales
      // las condiciones de contorno son todas de dirichlet (esenciales)
      // asi que van despues
      wasora_call(sn_elements_build_elemental_objects(&wasora_mesh.main_mesh->element[i]));
      
    }
  }
  
  wasora_call(sn_elements_set_essential_bc());
  
  // ensamblamos las matrices (para la boludina del mpi)
  wasora_call(milonga_assembly_objects(MAT_FINAL_ASSEMBLY));  
  
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


  
  
  
#undef  __FUNCT__
#define __FUNCT__ "sn_elements_build_elemental_objects"
int sn_elements_build_elemental_objects(element_t *element) {
  int v;           // punto de gauss
  int m, m_prime;  // indice de direccion
  int g, g_prime;  // grupo de energia
  int j;           // nodo local
  int n;           // dimension
  int this_element_has_fission = 0;
  int this_element_has_sources = 0;
  
  double tau;
  double w_gauss, xi;
  
  xs_t *material_xs;
  
  PetscFunctionBegin;
  
  if (Source == NULL) {
    Source = malloc(milonga.groups * sizeof(double));
    SigmaA = malloc(milonga.groups * sizeof(double));
    SigmaT = malloc(milonga.groups * sizeof(double));
    SigmaS0 = malloc(milonga.groups * sizeof(double *));
    SigmaS1 = malloc(milonga.groups * sizeof(double *));
    chinuSigmaF = malloc(milonga.groups * sizeof(double *));
    for (g = 0; g < milonga.groups; g++) {
      SigmaS0[g] = malloc(milonga.groups * sizeof(double));
      SigmaS1[g] = malloc(milonga.groups * sizeof(double));
      chinuSigmaF[g] = malloc(milonga.groups * sizeof(double));
    }
  }
  

  if (element->physical_entity == NULL) {
    wasora_push_error_message("element %d needs a physical entity", element->tag);
    PetscFunctionReturn(WASORA_RUNTIME_ERROR);
  } else if (element->physical_entity->material == NULL) {
    wasora_push_error_message("physical entity '%s' needs a material", element->physical_entity->material);
    PetscFunctionReturn(WASORA_RUNTIME_ERROR);
  }

  material_xs = (xs_t *)(element->physical_entity->material->ext);

  if (J != element->type->nodes) {
    wasora_call(sn_elements_allocate_particular_elemental_objects(element));
  }

  // inicializar Ki Ai Xi Si <- 0
  gsl_matrix_set_zero(Ki);
  gsl_matrix_set_zero(Ai);
  gsl_matrix_set_zero(Xi);
  gsl_vector_set_zero(Si);

  // factor de estabilizacion
  tau = wasora_var(milonga.vars.sn_alpha) * 0.5 * gsl_hypot3(element->node[1]->x[0]-element->node[0]->x[0],
                                                             element->node[1]->x[1]-element->node[0]->x[1],
                                                             element->node[1]->x[2]-element->node[0]->x[2]);
  
  // para cada punto de gauss
  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {

    // para este punto de gauss, calculamos las matrices H y B
    w_gauss = mesh_compute_fem_objects_at_gauss(wasora_mesh.main_mesh, element, v);
    
    // la estabilizacion de petrov
    for (j = 0; j < element->type->nodes; j++) {
      xi = element->type->h(j, wasora_mesh.main_mesh->fem.r);
      for (m = 0; m < milonga.directions; m++) {
        for (g = 0; g < milonga.groups; g++) {
          // parte base, igual a las h
          gsl_matrix_set(P, dof_index(m,g), milonga.directions*milonga.groups*j + dof_index(m,g), xi);
          // correccion
          for (n = 0; n < milonga.dimensions; n++) {
            gsl_matrix_add_to_element(P, dof_index(m,g), milonga.directions*milonga.groups*j + dof_index(m,g), 
                               tau * Omega[m][n] * gsl_matrix_get(wasora_mesh.main_mesh->fem.dhdx, j, n));
          }
        }
      }
    }
    
    
    // inicializamos las matrices con las XS (estas si dependen de la formulacion)
    gsl_matrix_set_zero(A);
    gsl_matrix_set_zero(X);
    gsl_vector_set_zero(S);

    for (g = 0; g < milonga.groups; g++) {
     
      if (material_xs->S[g]->n_tokens != 0) {
        if ((Source[g] = wasora_evaluate_expression(material_xs->S[g])) != 0) {
          milonga.has_sources = 1;
          this_element_has_sources = 1;
        } else {
          Source[g] = 0;
        }
      }
      if (material_xs->SigmaA[g]->n_tokens != 0) {
        SigmaA[g] = wasora_evaluate_expression(material_xs->SigmaA[g]);
      } else {
        SigmaA[g] = 0;
      }
      if (material_xs->SigmaT[g]->n_tokens != 0) {
        SigmaT[g] = wasora_evaluate_expression(material_xs->SigmaT[g]);
      } else {
        SigmaT[g] = 0;
      }
      
      for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
        if (material_xs->SigmaS0[g][g_prime]->n_tokens != 0) {
          SigmaS0[g][g_prime] = wasora_evaluate_expression(material_xs->SigmaS0[g][g_prime]);
        } else {
          SigmaS0[g][g_prime] = 0;
        }
        if (material_xs->SigmaS1[g][g_prime]->n_tokens != 0) {
          SigmaS1[g][g_prime] = wasora_evaluate_expression(material_xs->SigmaS1[g][g_prime]);
        } else {
          SigmaS1[g][g_prime] = 0;
        }
        if (material_xs->nuSigmaF[g]->n_tokens != 0 && gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), g_prime) != 0) {
          chinuSigmaF[g][g_prime] = gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), g_prime) * wasora_evaluate_expression(material_xs->nuSigmaF[g]);
          this_element_has_fission = 1;
          milonga.has_fission = 1;
        } else {
          chinuSigmaF[g][g_prime] = 0;
        }
      }
    }
       
      
    for (m = 0; m < milonga.directions; m++) {
      for (g = 0; g < milonga.groups; g++) {
        gsl_vector_set(S, dof_index(m,g), Source[g]);
        
        // scattering y fision
        for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
          for (m_prime = 0; m_prime < milonga.directions; m_prime++) {
            // scattering
            xi = -w[m_prime] * SigmaS0[g_prime][g];
            // si tenemos scattering anisotropico, l = 1
            if (material_xs->SigmaS1[g_prime][g]->n_tokens != 0) {
              xi -= w[m_prime] * SigmaS1[g_prime][g] * 3.0 * mesh_dot(Omega[m], Omega[m_prime]);
            }
            gsl_matrix_set(A, dof_index(m,g), dof_index(m_prime,g_prime), xi);

            // fision
            gsl_matrix_set(X, dof_index(m,g), dof_index(m_prime,g_prime), +w[m_prime] * chinuSigmaF[g_prime][g]);
          }
        }

        // absorcion
        xi = gsl_matrix_get(A, dof_index(m,g), dof_index(m,g));
        if (SigmaT[g] != 0) {
          xi += SigmaT[g];
        } else {
          xi += SigmaA[g];
          for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
            xi += SigmaS0[g][g_prime];
          }
        }
        gsl_matrix_set(A, dof_index(m,g), dof_index(m,g), xi);
      }
    }
    
    // armamos la matriz elemental del termino de fugas (estabilizada con petrov)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, OMEGA, wasora_mesh.main_mesh->fem.B, 0, OMEGAB);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, P, OMEGAB, 1, Ki);    

    // la matriz elemental de scattering (estabilizada con petrov)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, wasora_mesh.main_mesh->fem.H, 0, AH);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, P, AH, 1, Ai);
    
    // la matriz elemental de fision (estabilizada con petrov)
    if (this_element_has_fission) {
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, X, wasora_mesh.main_mesh->fem.H, 0, XH);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, P, XH, 1, Xi);
    }
    // el vector elemental de fuentes (estabilizado con petrov)
    if (this_element_has_sources) {
      gsl_blas_dgemv(CblasTrans, w_gauss, P, S, 1, Si);
    }
  }
  
  petsc_call(MatSetValues(milonga.R, L, wasora_mesh.main_mesh->fem.l, L, wasora_mesh.main_mesh->fem.l, gsl_matrix_ptr(Ki, 0, 0), ADD_VALUES));
  petsc_call(MatSetValues(milonga.R, L, wasora_mesh.main_mesh->fem.l, L, wasora_mesh.main_mesh->fem.l, gsl_matrix_ptr(Ai, 0, 0), ADD_VALUES));  
  if (this_element_has_fission) {
    petsc_call(MatSetValues(milonga.F, L, wasora_mesh.main_mesh->fem.l, L, wasora_mesh.main_mesh->fem.l, gsl_matrix_ptr(Xi, 0, 0), ADD_VALUES));
  }
  if (this_element_has_sources) {
    petsc_call(VecSetValues(milonga.S, L, wasora_mesh.main_mesh->fem.l, gsl_vector_ptr(Si, 0), ADD_VALUES));
  }
    

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "sn_elements_free_elemental_objects"
int sn_elements_free_elemental_objects(void) {

  PetscFunctionBegin;
  
  J = 0;
  L = 0;

  // matrices de interpolacion
  gsl_matrix_free(wasora_mesh.main_mesh->fem.H);
  wasora_mesh.main_mesh->fem.H = NULL;
  gsl_matrix_free(wasora_mesh.main_mesh->fem.B);
  wasora_mesh.main_mesh->fem.B = NULL;

  // matrices de coeficientes
  gsl_matrix_free(OMEGA);
  OMEGA = NULL;
  gsl_matrix_free(A);
  A = NULL;
  gsl_matrix_free(X);
  X = NULL;

  // matrices intermedias
  gsl_matrix_free(OMEGAB);
  OMEGAB = NULL;
  gsl_matrix_free(AH);
  AH = NULL;
  gsl_matrix_free(XH);
  XH = NULL;

  // matriz elemental de rigidez
  gsl_matrix_free(Ki);
  Ki = NULL;
  // matriz elemental de scattering
  gsl_matrix_free(Ai);
  Ai = NULL;
  // matriz elemental de fision
  gsl_matrix_free(Xi);
  Xi = NULL;

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


#undef  __FUNCT__
#define __FUNCT__ "sn_elements_set_essential_bc"
int sn_elements_set_essential_bc(void) {
  double outward_normal[3] = {0, 0, 0};
  element_t *surface_element;
  Vec zeros;
  int i, j;
  int bc_type;
  
  PetscFunctionBegin;

  k_vacuum = 0;
  current_size_vacuum = milonga.problem_size/BC_FACTOR;
  current_threshold_vacuum = current_size_vacuum - 2*milonga.groups;
  indexes_vacuum = malloc(current_size_vacuum * sizeof(int));

  k_mirror = 0;
  current_size_mirror = milonga.problem_size/BC_FACTOR;
  current_threshold_mirror = current_size_mirror - 2*milonga.groups;
  indexes_mirror = malloc(current_size_mirror * sizeof(int));
  indexes_mirror_reflected = malloc(current_size_mirror * sizeof(int));


/* 
  TODO: elegir barrer nodos o elementos 
  for (k = 0; k < wasora_mesh.main_mesh->n_nodes; k++) {
    if ((surface_element = mesh_find_node_neighbor_of_dim(&wasora_mesh.main_mesh->node[k], milonga.dimensions-1)) != NULL) {
      if (surface_element->physical_entity == NULL || surface_element->physical_entity->bc_type_int != BC_UNDEFINED) {
 */
  for (i = 0; i < wasora_mesh.main_mesh->n_elements; i++) {
    if (wasora_mesh.main_mesh->element[i].type->dim == milonga.dimensions-1) {
      surface_element = &wasora_mesh.main_mesh->element[i];
      if (surface_element->physical_entity == NULL || surface_element->physical_entity->bc_type_phys == BC_NULL
                                                   || surface_element->physical_entity->bc_type_phys == BC_UNDEFINED) {
        bc_type = BC_VACUUM;
      } else {
        bc_type = surface_element->physical_entity->bc_type_phys;
      }
        
      sn_elements_compute_outward_normal(surface_element, outward_normal);
      if (bc_type == BC_VACUUM) {
        for (j = 0; j < surface_element->type->nodes; j++) {
          sn_elements_add_vacuum_bc(surface_element->node[j], outward_normal);
        }
      } else if (bc_type == BC_MIRROR) {
        for (j = 0; j < surface_element->type->nodes; j++) {
          sn_elements_add_mirror_bc(surface_element->node[j], outward_normal);
        }
      }
    }
  }

  // hay que ensamblar porque creo que el chiste es que MatZeroRows usa INSERT_VALUES
  wasora_call(milonga_assembly_objects(MAT_FINAL_ASSEMBLY));

  // psi_n_g = 0
  if (k_vacuum > 0) {
    if (milonga.has_fission) {
      petsc_call(MatZeroRows(milonga.F, k_vacuum, indexes_vacuum, 0.0, PETSC_NULL, PETSC_NULL));
    }
    if (milonga.has_sources) {
      petsc_call(VecCreateSeq(PETSC_COMM_SELF, k_vacuum, &zeros));
      petsc_call(MatZeroRows(milonga.R, k_vacuum, indexes_vacuum, 1.0, zeros, milonga.S));
      VecDestroy(&zeros);
    } else {
      petsc_call(MatZeroRows(milonga.R, k_vacuum, indexes_vacuum, 1.0, PETSC_NULL, PETSC_NULL));
    }
  }
  
  // psi_n_g = psi_reflectedn_g
  if (k_mirror > 0) {
    if (milonga.has_fission) {
      petsc_call(MatZeroRows(milonga.F, k_mirror, indexes_mirror, 0.0, PETSC_NULL, PETSC_NULL));
    }

    if (milonga.has_sources) {
      petsc_call(VecCreateSeq(PETSC_COMM_SELF, k_mirror, &zeros));
      petsc_call(MatZeroRows(milonga.R, k_mirror, indexes_mirror, 1.0, zeros, milonga.S));
      petsc_call(VecDestroy(&zeros));
    } else {
      petsc_call(MatZeroRows(milonga.R, k_mirror, indexes_mirror, 1.0, PETSC_NULL, PETSC_NULL));
    }
    // TODO: improve!
    for (i = 0; i < k_mirror; i++) {
      petsc_call(MatSetValue(milonga.R, indexes_mirror[i], indexes_mirror_reflected[i], -1.0, INSERT_VALUES));
    }
  }
  
  
  free(indexes_vacuum);
  free(indexes_mirror);
  free(indexes_mirror_reflected);

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "sn_elements_compute_outward_normal"
int sn_elements_compute_outward_normal(element_t *surface_element, double *outward_normal) {
  double a[3], b[3];
  double surface_center[3];
  double volumetric_neighbor_center[3];
  element_t *volumetric_neighbor;

  PetscFunctionBegin;

  if ((volumetric_neighbor = mesh_find_element_volumetric_neighbor(surface_element)) == NULL) {
    wasora_push_error_message("cannot find any volumetric neighbor for surface element %d", surface_element->tag);
    PetscFunctionReturn(WASORA_RUNTIME_ERROR);
  }

  if (milonga.dimensions == 1) {
    // para una dimension el chiste es straightforward
    // OJO que esto no camina si el segment no esta en el eje x
    if (surface_element->node[0]->x[0] < 0.5*(volumetric_neighbor->node[0]->x[0] + volumetric_neighbor->node[1]->x[0])) {
      outward_normal[0] = -1;
    } else {
      outward_normal[0] = +1;
    }
  } else {
    // la normal exterior depende de si el problema es 2d o 3d
    // por la cantidad de nodos del elemento de la superficie
    if (surface_element->type->nodes == 2) {
      // esto es una linea
      // OJO que no camina con lineas que no estan en el plano xy!!
      double module = mesh_subtract_module(surface_element->node[1]->x, surface_element->node[0]->x);
      outward_normal[0] = -(surface_element->node[1]->x[1] - surface_element->node[0]->x[1])/module;
      outward_normal[1] = +(surface_element->node[1]->x[0] - surface_element->node[0]->x[0])/module;
      outward_normal[2] = 0;
    } else {
      // este puede ser un triangulo o un cuadrangulo, en cualquier caso usamos
      // los tres primeros nodos y ya
      a[0] = surface_element->node[1]->x[0] - surface_element->node[0]->x[0];
      a[1] = surface_element->node[1]->x[1] - surface_element->node[0]->x[1];
      a[2] = surface_element->node[1]->x[2] - surface_element->node[0]->x[2];
      b[0] = surface_element->node[2]->x[0] - surface_element->node[0]->x[0];
      b[1] = surface_element->node[2]->x[1] - surface_element->node[0]->x[1];
      b[2] = surface_element->node[2]->x[2] - surface_element->node[0]->x[2];
      mesh_normalized_cross(a, b, outward_normal);
    }

    // ahora tenemos que ver si la normal que elegimos es efectivamente la outward
    // para eso primero calculamos el centro del elemento de superficie
    wasora_call(mesh_compute_element_barycenter(surface_element, surface_center));

    // y despues el centro del elemento de volumen
    wasora_call(mesh_compute_element_barycenter(volumetric_neighbor, volumetric_neighbor_center));

    // calculamos el producto entre la normal propuesta y la resta de estos dos vectores
    // si elegimos la otra direccion, la damos tavuel
    if (mesh_subtract_dot(volumetric_neighbor_center, surface_center, outward_normal) > 0) {
      outward_normal[0] = -outward_normal[0];
      outward_normal[1] = -outward_normal[1];
      outward_normal[2] = -outward_normal[2];
    }
  }
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "sn_elements_add_vacuum_bc"
int sn_elements_add_vacuum_bc(node_t *node, double *outward_normal) {
  
  int n;
  
  PetscFunctionBegin;
 
  for (n = 0; n < milonga.directions; n++) {
    if (mesh_dot(Omega[n], outward_normal) < 0) {
      // si el producto interno de Omega con la normal es negativo entonces hacemos cero el flujo
      if (k_vacuum > current_threshold_vacuum) {
        current_size_vacuum += milonga.problem_size/BC_FACTOR;
        current_threshold_vacuum = current_size_vacuum - 2*milonga.groups;
        indexes_vacuum = realloc(indexes_vacuum, current_size_vacuum * sizeof(int));
      }

      petsc_call(PetscMemcpy(indexes_vacuum+k_vacuum, node->index_dof+dof_index(n,0), milonga.groups * sizeof(int)));
      k_vacuum += milonga.groups;
    }
  }
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "sn_elements_add_mirror_bc"
int sn_elements_add_mirror_bc(node_t *node, double *outward_normal) {
  
  double reflected[3] = {0, 0, 0};
  double Omega_dot_outward;
  int n, n_refl;
  
  PetscFunctionBegin;
  for (n = 0; n < milonga.directions; n++) {
    if ((Omega_dot_outward = mesh_dot(Omega[n], outward_normal)) < 0) {

      // si el producto interno de Omega con la normal es positivo entonces tenemos que reflejar
      // si Omega es la direccion de incidencia, la direccion reflejada con respecto a la normal outward_normal es 
      // reflected = Omega - 2*(Omega dot outward_normal) * outward_normal
      reflected[0] = Omega[n][0] - 2*Omega_dot_outward * outward_normal[0];
      reflected[1] = Omega[n][1] - 2*Omega_dot_outward * outward_normal[1];
      reflected[2] = Omega[n][2] - 2*Omega_dot_outward * outward_normal[2];
      for (n_refl = 0; n_refl < milonga.directions; n_refl++) {
        if (fabs(reflected[0]-Omega[n_refl][0]) < wasora_var(wasora_mesh.vars.eps) &&
            fabs(reflected[1]-Omega[n_refl][1]) < wasora_var(wasora_mesh.vars.eps) &&
            fabs(reflected[2]-Omega[n_refl][2]) < wasora_var(wasora_mesh.vars.eps)) {
          break;
        }
      }
      if (n_refl == milonga.directions) {
        wasora_push_error_message("cannot find a reflected direction for n=%d in node %d", n, node->tag);
        PetscFunctionReturn(WASORA_RUNTIME_ERROR);
      }

      if (k_mirror > current_threshold_mirror) {
        current_size_mirror += milonga.problem_size/BC_FACTOR;
        current_threshold_mirror = current_size_mirror - 2*milonga.groups;
        indexes_mirror = realloc(indexes_mirror, current_size_mirror * sizeof(int));
        indexes_mirror_reflected = realloc(indexes_mirror_reflected, current_size_mirror * sizeof(int));
      }
      petsc_call(PetscMemcpy(indexes_mirror          +k_mirror, node->index_dof+dof_index(n,0),      milonga.groups * sizeof(int)));
      petsc_call(PetscMemcpy(indexes_mirror_reflected+k_mirror, node->index_dof+dof_index(n_refl,0), milonga.groups * sizeof(int)));
      k_mirror += milonga.groups;
    }
  }
      
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}



#undef  __FUNCT__
#define __FUNCT__ "sn_elements_results_fill_flux"
int sn_elements_results_fill_flux(void) {

  int k, g, n;
  
  PetscFunctionBegin;
  

  for (k = 0; k < wasora_mesh.main_mesh->n_nodes; k++) {
    for (g = 0; g < milonga.groups; g++) {
      milonga.functions.phi[g]->data_value[k] = 0;
      for (n = 0; n < milonga.directions; n++) {
        VecGetValues(milonga.phi, 1, &wasora_mesh.main_mesh->node[k].index_dof[dof_index(n,g)], &milonga.functions.psi[n][g]->data_value[k]);
        milonga.functions.phi[g]->data_value[k] += w[n] * milonga.functions.psi[n][g]->data_value[k];
      }
    }
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "sn_elements_normalize_flux"
int sn_elements_normalize_flux(void) {

  int i, k, g, n;
  double factor;
  double num = 0;
  double den = 0;
  
  PetscFunctionBegin;

  if (wasora_var(milonga.vars.power) == 0) {

    // calculamos el factor de normalizacion
    for (i = 0; i < wasora_mesh.main_mesh->n_elements; i++) {
      if (wasora_mesh.main_mesh->element[i].type->dim == wasora_mesh.main_mesh->bulk_dimensions) {
        num += wasora_mesh.main_mesh->element[i].type->element_volume(&wasora_mesh.main_mesh->element[i]);
        for (g = 0; g < milonga.groups; g++) {
          den += mesh_integral_over_element(milonga.functions.phi[g], &wasora_mesh.main_mesh->element[i], NULL);
        }
      }
    }

  } else {

    xs_t *xs;

    num = wasora_var(milonga.vars.power);
    for (i = 0; i < wasora_mesh.main_mesh->n_elements; i++) {
      if (wasora_mesh.main_mesh->element[i].type->dim == wasora_mesh.main_mesh->bulk_dimensions) {
        if ((xs = (xs_t *)wasora_mesh.main_mesh->element[i].physical_entity->material->ext) == NULL) {
          wasora_push_error_message("physical entity '%s' needs a material", wasora_mesh.main_mesh->cell[i].element->physical_entity->name);
          PetscFunctionReturn(WASORA_RUNTIME_ERROR);
        }

        for (g = 0; g < milonga.groups; g++) {
          den += mesh_integral_over_element(milonga.functions.phi[g], &wasora_mesh.main_mesh->element[i], xs->eSigmaF[g]);
        }
      }
    }

    if (den == 0) {
      wasora_push_error_message("power setpoint was given but eSigmaF is identically zero");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    }
  }

  factor = num/den;

  // normalizamos los valores de las funciones flujo
  for (k = 0; k < wasora_mesh.main_mesh->n_nodes; k++) {
    for (g = 0; g < milonga.groups; g++) {
      milonga.functions.phi[g]->data_value[k] *= factor;
    }
  }
  
  for (n = 0; n < milonga.directions; n++) {
    for (g = 0; g < milonga.groups; g++) {
      for (i = 0; i < milonga.spatial_unknowns; i++) {
        milonga.functions.psi[n][g]->data_value[i] *= factor;
      }
    }
  }
  

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "sn_elements_results_fill_power"
int sn_elements_results_fill_power(void) {
  int g, k;
  double pow;
  xs_t *xs;
  element_t *element = NULL;
  element_list_item_t *associated_element = NULL;

  PetscFunctionBegin;

  for (k = 0; k < wasora_mesh.main_mesh->n_nodes; k++) {
    LL_FOREACH(wasora_mesh.main_mesh->node[k].associated_elements, associated_element) {
      element = associated_element->element;
    }
    if (element != NULL && element->physical_entity != NULL && element->physical_entity->material != NULL) {
      xs = (xs_t *)element->physical_entity->material->ext;
      if (xs != NULL && xs->eSigmaF != NULL) {
        wasora_var(wasora_mesh.vars.x) = wasora_mesh.main_mesh->node[k].x[0];
        wasora_var(wasora_mesh.vars.y) = wasora_mesh.main_mesh->node[k].x[1];
        wasora_var(wasora_mesh.vars.z) = wasora_mesh.main_mesh->node[k].x[2];

        pow = 0;
        for (g = 0; g < milonga.groups; g++) {
          pow += wasora_evaluate_expression(xs->eSigmaF[g]) * milonga.functions.phi[g]->data_value[k];
        }
        milonga.functions.pow->data_value[k] = pow;
      }
    }
  }

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


#undef  __FUNCT__
#define __FUNCT__ "sn_elements_problem_free"
int sn_elements_problem_free(void) {
  
  int g, n;
    
  PetscFunctionBegin;

  if (wasora_mesh.main_mesh != NULL && wasora_mesh.main_mesh->n_elements != 0) {
    if (milonga.functions.phi != NULL) {
      for (g = 0; g < milonga.groups; g++) {
        free(milonga.functions.phi[g]->data_value);
        milonga.functions.phi[g]->data_argument = NULL;
        milonga.functions.phi[g]->data_value = NULL;
        milonga.functions.phi[g]->var_argument = NULL;
      }
      free(milonga.functions.pow->data_value);
      milonga.functions.pow->data_argument = NULL;
      milonga.functions.pow->data_value = NULL;
      milonga.functions.pow->var_argument = NULL;
    }
    if (milonga.functions.psi != NULL) {
      // psi is allocated in parser.c end never de-allocated.
      // it is a ** function [n][g] where n is the number of directions (SN method)
      // and g the number of groups
      for (n = 0; n < milonga.directions; n++) {
	for (g = 0; g < milonga.groups; g++) {
	  free(milonga.functions.psi[n][g]->data_value);
	}
	free(milonga.functions.psi[n]);
      }
    }
  }
  
  sn_elements_free_elemental_objects();
  mesh_free(wasora_mesh.main_mesh);
  
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}
