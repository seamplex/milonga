/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's diffusion with finite elements
 *
 *  Copyright (C) 2012--2015 jeremy theler
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

#include "../milonga.h"
#include "diffusion_elements.h"

// nomenclatura como en la documentacion
int J;            // cantidad de nodos locales
int L;            // tamanio de la matriz elemental

// matrices de coeficientes
gsl_matrix *D;
gsl_matrix *A;
gsl_matrix *X;

// vector de fuentes
gsl_vector *S;

// matrices intermedias
gsl_matrix *DB;
gsl_matrix *AH;
gsl_matrix *XH;

// matriz elemental de rigidez
gsl_matrix *Ki;
// matriz elemental de scattering
gsl_matrix *Ai;
// matriz elemental de fision
gsl_matrix *Xi;

// matriz elemental de superficies de contorno
gsl_matrix *Ni;

// vector elemental de fuentes
gsl_vector *Si;

#undef  __FUNCT__
#define __FUNCT__ "diffusion_elements_problem_init"
int diffusion_elements_problem_init(void) {

  int g;

  milonga.spatial_unknowns = wasora_mesh.main_mesh->n_nodes;
  milonga_allocate_global_matrices(milonga.spatial_unknowns * milonga.groups,
                                   wasora_mesh.main_mesh->max_first_neighbor_nodes * milonga.groups,
                                   wasora_mesh.main_mesh->max_first_neighbor_nodes * milonga.groups);
  milonga_allocate_global_vectors();

  wasora_var(wasora_mesh.vars.cells) = (double)wasora_mesh.main_mesh->n_cells;
  wasora_var(wasora_mesh.vars.nodes) = (double)wasora_mesh.main_mesh->n_nodes;
  wasora_var(wasora_mesh.vars.elements) = (double)wasora_mesh.main_mesh->n_elements;
  wasora_mesh.main_mesh->data_type = data_type_node;

  if (wasora_mesh.main_mesh->structured) {
    wasora_mesh_struct_init_rectangular_for_nodes(wasora_mesh.main_mesh);
  }
  
  for (g = 0; g < milonga.groups; g++) {
    wasora_call(diffusion_elements_results_fill_args(milonga.functions.phi[g]));
  }
  wasora_call(diffusion_elements_results_fill_args(milonga.functions.pow));

  wasora_call(mesh_node_indexes(wasora_mesh.main_mesh, milonga.groups));

  return WASORA_RUNTIME_OK;
}


// esta rutina rellena datos administrativos de las funciones resultados
#undef  __FUNCT__
#define __FUNCT__ "results_fill_args_elements"
int diffusion_elements_results_fill_args(function_t *function) {

  // tenemos data
  function->data_size = milonga.spatial_unknowns;
  if (wasora_mesh.main_mesh->structured) {
    function->rectangular_mesh = 1;
    function->x_increases_first = 1;
    function->rectangular_mesh_size = wasora_mesh.main_mesh->rectangular_mesh_size;
    function->rectangular_mesh_point = wasora_mesh.main_mesh->rectangular_mesh_point;    
  }

  // pero tambien variables por si queremos hacer cuentitas
  function->var_argument = calloc(3, sizeof(var_t *));
  function->var_argument = wasora_mesh.vars.arr_x;

  function->data_argument = wasora_mesh.main_mesh->nodes_argument;
  function->data_value = calloc(function->data_size, sizeof(double));


  // y tipo milonga_status.mesh node en elementos
  function->type = type_pointwise_mesh_node;
  function->multidim_threshold = DEFAULT_MULTIDIM_INTERPOLATION_THRESHOLD;
  function->mesh = wasora_mesh.main_mesh;

  return WASORA_RUNTIME_OK;;
}


#undef  __FUNCT__
#define __FUNCT__ "diffusion_elements_allocate_general_elemental_objects"
int diffusion_elements_allocate_general_elemental_objects(void) {
  D = gsl_matrix_calloc(milonga.groups * milonga.dimensions, milonga.groups * milonga.dimensions);
  
  A = gsl_matrix_calloc(milonga.groups, milonga.groups);
  X = gsl_matrix_calloc(milonga.groups, milonga.groups);
  S = gsl_vector_calloc(milonga.groups);
  
  return WASORA_RUNTIME_OK;
  
}

#undef  __FUNCT__
#define __FUNCT__ "diffusion_elements_allocate_particular_elemental_objects"
int diffusion_elements_allocate_particular_elemental_objects(element_t *element) {

  J = element->type->nodes;
  L = milonga.groups * element->type->nodes;

  // TODO: esta las tendria que alocar mesh
  gsl_matrix_free(wasora_mesh.main_mesh->fem.H);
  wasora_mesh.main_mesh->fem.H = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom, L);
  
  gsl_matrix_free(wasora_mesh.main_mesh->fem.B);
  wasora_mesh.main_mesh->fem.B = gsl_matrix_calloc(wasora_mesh.main_mesh->degrees_of_freedom * wasora_mesh.main_mesh->bulk_dimensions, L);

  gsl_matrix_free(DB); 
  DB = gsl_matrix_calloc(milonga.groups * milonga.dimensions, L);
  gsl_matrix_free(AH); 
  AH = gsl_matrix_calloc(milonga.groups, L);
  gsl_matrix_free(XH); 
  XH = gsl_matrix_calloc(milonga.groups, L);
  

  gsl_matrix_free(Ki);
  Ki = gsl_matrix_calloc(L, L);
  gsl_matrix_free(Ai);
  Ai = gsl_matrix_calloc(L, L);
  gsl_matrix_free(Xi);
  Xi = gsl_matrix_calloc(L, L);
  gsl_matrix_free(Ni);
  Ni = gsl_matrix_calloc(L, L);
  gsl_vector_free(Si);
  Si = gsl_vector_calloc(L);
  
  return WASORA_RUNTIME_OK;

}


#undef  __FUNCT__
#define __FUNCT__ "matrices_build_elements"
int diffusion_elements_matrices_build(void) {
  
  int i;
  
  wasora_call(diffusion_elements_allocate_general_elemental_objects());

  for (i = 0; i < wasora_mesh.main_mesh->n_elements; i++) {
    
    if (wasora_mesh.main_mesh->element[i].type != NULL && wasora_mesh.main_mesh->element[i].type->dim == wasora_mesh.main_mesh->bulk_dimensions) {

      // solo los elementos que tengan la dimension del problema
      // son los que usamos para las matrices elementales
      wasora_call(diffusion_elements_build_volume_objects(&wasora_mesh.main_mesh->element[i]));
      
    } else if (wasora_mesh.main_mesh->element[i].type != NULL && wasora_mesh.main_mesh->element[i].type->dim == wasora_mesh.main_mesh->bulk_dimensions-1) {
      
      // si tienen dimension dim-1 entonces son candidatos a condiciones de contorno
      // pero aca miramos solo las de neumann y de robin porque las de dirichlet van
      // una vez que ensamblamos las matrizotas
      
      if (wasora_mesh.main_mesh->element[i].physical_entity == NULL ||
          wasora_mesh.main_mesh->element[i].physical_entity->bcs == NULL ||
          wasora_mesh.main_mesh->element[i].physical_entity->bcs->type_phys == BC_VACUUM ||
          wasora_mesh.main_mesh->element[i].physical_entity->bcs->type_phys == BC_UNDEFINED) {
          wasora_call(diffusion_elements_build_robin_objects(&wasora_mesh.main_mesh->element[i]));
      } else if (wasora_mesh.main_mesh->element[i].physical_entity->bcs->type_phys == BC_MIRROR) {
          // TODO: que se puedan poner corrientes no nulas (no debe haber o fision o fuentes)
          ; // no hay que hacer naranja!
      }
    }
  }
  
  wasora_call(diffusion_elements_set_essential_bc());
  
  // ensamblamos las matrices (para la boludina del mpi)
  wasora_call(milonga_assembly_objects(MAT_FINAL_ASSEMBLY));  
  
  
  return WASORA_RUNTIME_OK;
}





#undef  __FUNCT__
#define __FUNCT__ "diffusion_elements_build_volume_objects"
int diffusion_elements_build_volume_objects(element_t *element) {
  int v;           // punto de gauss
  int m;           // dimension
  int g, g_prime;  // grupo de energia
  int this_element_has_fission = 0;
  int this_element_has_sources = 0;

  double w_gauss, xi;
  
  xs_t *material_xs;

  if (element->physical_entity == NULL) {
    // esto pasa solo en malla estructuradas
    element_list_item_t *associated_element;
    int j;
    int nullflux;
  
    for (j = 0; j < element->type->nodes; j++) {
      // suponemos que tenemos que hacer cero el flujo
      nullflux = 1;
      LL_FOREACH(element->node[j]->associated_elements, associated_element) {
        if (associated_element->element->physical_entity != NULL) {
          // si alguno de los elementos asociados al nodo tiene entidad fisica
          // entonces no tocamos nada
          nullflux = 0;
        }
      }
      if (nullflux) {
        for (g = 0; g < milonga.groups; g++) {
          // phi = 0
          MatSetValue(milonga.R, element->node[j]->index_dof[g], element->node[j]->index_dof[g], 1.0, ADD_VALUES);
        }
      }
    }
  } else {
    
    if (element->physical_entity->material == NULL) {
      wasora_push_error_message("physical entity '%s' needs a material", element->physical_entity->name);
      return WASORA_RUNTIME_ERROR;
    }

    material_xs = (xs_t *)(element->physical_entity->material->ext);

    if (J != element->type->nodes) {
      wasora_call(diffusion_elements_allocate_particular_elemental_objects(element));
    }  

    // inicializar Ki Ai Xi Si <- 0
    gsl_matrix_set_zero(Ki);
    gsl_matrix_set_zero(Ai);
    gsl_matrix_set_zero(Xi);
    gsl_vector_set_zero(Si);

    // para cada punto de gauss
    for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {

      // para este punto de gauss, calculamos las matrices H y B
      w_gauss = mesh_compute_fem_objects_at_gauss(wasora_mesh.main_mesh, element, v);    

      // inicializamos las matrices con las XS (estas si dependen de la formulacion)
      gsl_matrix_set_zero(A);
      gsl_matrix_set_zero(X);
      gsl_vector_set_zero(S);
      gsl_matrix_set_zero(D);

      for (g = 0; g < milonga.groups; g++) {

        // fuentes
        if (material_xs->S[g]->n_tokens != 0) {
          if ((xi = wasora_evaluate_expression(material_xs->S[g])) != 0) {
            milonga.has_sources = 1;
            this_element_has_sources = 1;
            gsl_vector_set(S, g, xi);
          }
        }

        // scattering y fision
        for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
          if ((xi = -wasora_evaluate_expression(material_xs->SigmaS0[g_prime][g])) != 0) {
            gsl_matrix_set(A, g, g_prime, xi);
          }
          if ((xi = gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), g) * wasora_evaluate_expression(material_xs->nuSigmaF[g_prime])) != 0) {
            milonga.has_fission = 1;
            this_element_has_fission = 1;
            gsl_matrix_set(X, g, g_prime, xi);
          }
        }

        // absorcion
        xi = gsl_matrix_get(A, g, g);
        if (material_xs->SigmaT[g] != NULL && material_xs->SigmaT[g]->n_tokens != 0) {
          xi += wasora_evaluate_expression(material_xs->SigmaT[g]);
        } else {
          xi += wasora_evaluate_expression(material_xs->SigmaA[g]);
          for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
            xi += wasora_evaluate_expression(material_xs->SigmaS0[g][g_prime]);
          }
        }
        gsl_matrix_set(A, g, g, xi);

        // fugas
        for (m = 0; m < milonga.dimensions; m++) {
          if (material_xs->D[g]->n_tokens != 0) {
            xi = wasora_evaluate_expression(material_xs->D[g]);
          } else if (material_xs->SigmaT[g]->n_tokens != 0) {
            xi = 1.0/(3.0 * wasora_evaluate_expression(material_xs->SigmaT[g]));
          } else {
            wasora_push_error_message("neither D nor SigmaT given for material '%s' (diffusion does not handle void)", element->physical_entity->material->name);
            return WASORA_RUNTIME_ERROR;
          }
          if (xi != 0) {
            milonga.has_diffusion = 1;
          }
          gsl_matrix_set(D, g+milonga.groups*m,  g+milonga.groups*m, xi);
        }
      }

      // armamos la matriz elemental del termino de difusion
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, D, wasora_mesh.main_mesh->fem.B, 0, DB);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, wasora_mesh.main_mesh->fem.B, DB, 1, Ki);

      // la matriz elemental de scattering
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, wasora_mesh.main_mesh->fem.H, 0, AH);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, wasora_mesh.main_mesh->fem.H, AH, 1, Ai);

      // la matriz elemental de fision
      if (this_element_has_fission) {
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, X, wasora_mesh.main_mesh->fem.H, 0, XH);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, wasora_mesh.main_mesh->fem.H, XH, 1, Xi);
      }
      // el vector elemental de fuentes
      if (this_element_has_sources) {
        gsl_blas_dgemv(CblasTrans, w_gauss, wasora_mesh.main_mesh->fem.H, S, 1, Si);
      }
    }

    MatSetValues(milonga.R, L, wasora_mesh.main_mesh->fem.l, L, wasora_mesh.main_mesh->fem.l, gsl_matrix_ptr(Ki, 0, 0), ADD_VALUES);
    MatSetValues(milonga.R, L, wasora_mesh.main_mesh->fem.l, L, wasora_mesh.main_mesh->fem.l, gsl_matrix_ptr(Ai, 0, 0), ADD_VALUES);  
    if (this_element_has_fission) {
      MatSetValues(milonga.F, L, wasora_mesh.main_mesh->fem.l, L, wasora_mesh.main_mesh->fem.l, gsl_matrix_ptr(Xi, 0, 0), ADD_VALUES);
    }
    VecSetValues(milonga.S, L, wasora_mesh.main_mesh->fem.l, gsl_vector_ptr(Si, 0), ADD_VALUES);
    
  }
  
  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "diffusion_elements_build_robin_objects"
int diffusion_elements_build_robin_objects(element_t *element) {
  int v;
  double w_gauss;
  double a;
  
  if (J != element->type->nodes) {
    wasora_call(diffusion_elements_allocate_particular_elemental_objects(element));
  }  

  gsl_matrix_set_zero(Ni);

  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    w_gauss = mesh_compute_fem_objects_at_gauss(wasora_mesh.main_mesh, element, v);    
    // por default ponemos a = 1/2
    if (element->physical_entity == NULL ||
        element->physical_entity->bcs == NULL ||
        element->physical_entity->bcs->expr == NULL ||
        (a = fabs(wasora_evaluate_expression(element->physical_entity->bcs->expr))) == 0) {
      a = 0.5;
    }
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss * a, wasora_mesh.main_mesh->fem.H, wasora_mesh.main_mesh->fem.H, 1, Ni);
  }

  MatSetValues(milonga.R, L, wasora_mesh.main_mesh->fem.l, L, wasora_mesh.main_mesh->fem.l, gsl_matrix_ptr(Ni, 0, 0), ADD_VALUES);
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "elements_free_elemental_objects"
int diffusion_elements_free_elemental_objects(void) {

  J = 0;
  L = 0;

  // matrices de interpolacion
  gsl_matrix_free(wasora_mesh.main_mesh->fem.H);
  wasora_mesh.main_mesh->fem.H = NULL;
  
  gsl_matrix_free(wasora_mesh.main_mesh->fem.B);
  wasora_mesh.main_mesh->fem.B = NULL;

  // matrices de coeficientes
  gsl_matrix_free(D);
  D = NULL;
  gsl_matrix_free(A);
  A = NULL;
  gsl_matrix_free(X);
  X = NULL;

  // matrices intermedias
  gsl_matrix_free(DB);
  DB = NULL;
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

  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "elements_set_essential_bc"
int diffusion_elements_set_essential_bc(void) {

  int i, j;
  double *gzeros = calloc(milonga.groups, sizeof(double));

  // hay que ensamblar porque creo que el chiste es que MatZeroRows usa INSERT_VALUES
  wasora_call(milonga_assembly_objects(MAT_FINAL_ASSEMBLY));  

  for (i = 0; i < wasora_mesh.main_mesh->n_elements; i++) {
    if (wasora_mesh.main_mesh->element[i].type != NULL && wasora_mesh.main_mesh->element[i].type->dim == wasora_mesh.main_mesh->bulk_dimensions-1) {
      if (wasora_mesh.main_mesh->element[i].physical_entity == NULL ||
          wasora_mesh.main_mesh->element[i].physical_entity->bcs == NULL ||    
          wasora_mesh.main_mesh->element[i].physical_entity->bcs->type_phys == BC_NULL) {
        for (j = 0; j < wasora_mesh.main_mesh->element[i].type->nodes; j++) {
          if (milonga.has_fission) {
            MatZeroRows(milonga.F, milonga.groups, wasora_mesh.main_mesh->element[i].node[j]->index_dof, 0.0, PETSC_NULL, PETSC_NULL);
          }
          MatZeroRows(milonga.R, milonga.groups, wasora_mesh.main_mesh->element[i].node[j]->index_dof, 1.0, PETSC_NULL, PETSC_NULL);
          VecSetValues(milonga.S, milonga.groups, wasora_mesh.main_mesh->element[i].node[j]->index_dof, gzeros, INSERT_VALUES);
        }
      }
    }
  }
  
  free(gzeros);

  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "diffusion_elements_results_fill_flux"
int diffusion_elements_results_fill_flux(void) {

  int k, g;

  // rellenamos las funciones de los flujos con lo que dio PETSc
  for (k = 0; k < wasora_mesh.main_mesh->n_nodes; k++) {
    for (g = 0; g < milonga.groups; g++) {
      VecGetValues(milonga.phi, 1, &wasora_mesh.main_mesh->node[k].index_dof[g], &milonga.functions.phi[g]->data_value[k]);
    }
  }
  
  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "diffusion_elements_normalize_flux"
int diffusion_elements_normalize_flux(void) {

  int i, k, g;
  double factor;
  double num = 0;
  double den = 0;

  if (wasora_var(milonga.vars.power) == 0) {

    // calculamos el factor de normalizacion
    for (i = 0; i < wasora_mesh.main_mesh->n_elements; i++) {
      if (wasora_mesh.main_mesh->element[i].type != NULL && wasora_mesh.main_mesh->element[i].type->dim == wasora_mesh.main_mesh->bulk_dimensions) {
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
      if (wasora_mesh.main_mesh->element[i].type != NULL && wasora_mesh.main_mesh->element[i].type->dim == wasora_mesh.main_mesh->bulk_dimensions && wasora_mesh.main_mesh->element[i].physical_entity != NULL) {
        if ((xs = (xs_t *)wasora_mesh.main_mesh->element[i].physical_entity->material->ext) == NULL) {
          wasora_push_error_message("physical entity '%s' needs a material", wasora_mesh.main_mesh->cell[i].element->physical_entity->name);
          return WASORA_RUNTIME_ERROR;
        }

        for (g = 0; g < milonga.groups; g++) {
          den += mesh_integral_over_element(milonga.functions.phi[g], &wasora_mesh.main_mesh->element[i], xs->eSigmaF[g]);
        }
      }
    }

    if (den == 0) {
      wasora_push_error_message("power setpoint was given but eSigmaF is identically zero");
      return WASORA_RUNTIME_ERROR;
    }
  }

  factor = num/den;

  // normalizamos los valores de las funciones flujo
  for (k = 0; k < wasora_mesh.main_mesh->n_nodes; k++) {
    for (g = 0; g < milonga.groups; g++) {
      milonga.functions.phi[g]->data_value[k] *= factor;
    }
  }

  return WASORA_PARSER_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "diffusion_elements_results_fill_power"
int diffusion_elements_results_fill_power(void) {
  int g, k;
  double den;
  double pow;
  xs_t *xs;
  element_t *element = NULL;
  element_list_item_t *associated_element;


  for (k = 0; k < wasora_mesh.main_mesh->n_nodes; k++) {
    den = 0;
    LL_FOREACH(wasora_mesh.main_mesh->node[k].associated_elements, associated_element) {
      element = associated_element->element;
      pow = 0;
      if (element != NULL && element->physical_entity != NULL && element->physical_entity->material != NULL) {
        den += element->type->element_volume(element);
        xs = (xs_t *)element->physical_entity->material->ext;
        if (xs != NULL && xs->eSigmaF != NULL) {
          wasora_var(wasora_mesh.vars.x) = wasora_mesh.main_mesh->node[k].x[0];
          wasora_var(wasora_mesh.vars.y) = wasora_mesh.main_mesh->node[k].x[1];
          wasora_var(wasora_mesh.vars.z) = wasora_mesh.main_mesh->node[k].x[2];

          for (g = 0; g < milonga.groups; g++) {
            pow += wasora_evaluate_expression(xs->eSigmaF[g]) * milonga.functions.phi[g]->data_value[k];
          }
        }
        milonga.functions.pow->data_value[k] += element->type->element_volume(element) * pow;
      }
    }
    if (den != 0) {
      milonga.functions.pow->data_value[k] /= den;
    }
  }

  return WASORA_RUNTIME_OK;

}


#undef  __FUNCT__
#define __FUNCT__ "problem_free_elements"
int diffusion_elements_problem_free(void) {
  
  int g;

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
 
    diffusion_elements_free_elemental_objects();
    mesh_free(wasora_mesh.main_mesh);
  }
   
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_print_gsl_matrix"
int fino_print_gsl_matrix(gsl_matrix *A, FILE *file) {

  double xi;
  int i, j;

  for (i = 0; i < A->size1; i++) {
    for (j = 0; j < A->size2; j++) {
      xi = gsl_matrix_get(A, i, j);
      if (xi != 0) {
        fprintf(file, "% .1e ", xi);
      } else {
        fprintf(file, "    0    ");
      }
    }
    fprintf(file, "\n");
  }
  
  return WASORA_RUNTIME_OK;

}
