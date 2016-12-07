/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's diffusion with finite volumes routines
 *
 *  Copyright (C) 2012--2016 jeremy theler
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

#include "../milonga.h"
#include "diffusion_volumes.h"


#define flat_index(i,j,k) ((i) + (j)*milonga.mesh->ncells_x + (k)*milonga.mesh->ncells_x*milonga.mesh->ncells_y) 


#undef  __FUNCT__
#define __FUNCT__ "diffusion_volumes_problem_init"
int diffusion_volumes_problem_init(void) {

  int g;

  if (milonga.mesh == NULL) {
    wasora_push_error_message("no mesh found");
    return WASORA_RUNTIME_ERROR;
  }

  if (milonga.mesh->structured == 0) {
    if (milonga.mesh->cell == NULL) {
      wasora_call(mesh_element2cell(milonga.mesh));
    }
    if (milonga.mesh->cell[0].ifaces == NULL) {
      wasora_call(mesh_find_neighbors(milonga.mesh));
    }
  
    wasora_call(mesh_compute_coords(milonga.mesh));
    wasora_call(mesh_fill_neighbors(milonga.mesh));
  } else {
    wasora_mesh_struct_init_rectangular_for_cells(milonga.mesh);
  }


  milonga.spatial_unknowns = milonga.mesh->n_cells;
  milonga_allocate_global_objects(milonga.spatial_unknowns * milonga.groups,
                                  milonga.mesh->max_faces_per_element + milonga.groups,
                                  milonga.groups);
  
  wasora_var(wasora_mesh.vars.cells) = (double)milonga.mesh->n_cells;
  wasora_var(wasora_mesh.vars.nodes) = (double)milonga.mesh->n_nodes;
  wasora_var(wasora_mesh.vars.elements) = (double)milonga.mesh->n_elements;
  milonga.mesh->data_type = data_type_element;
  
  for (g = 0; g < milonga.groups; g++) {
    wasora_call(diffusion_volumes_results_fill_args(milonga.functions.phi[g]));
  }
  wasora_call(diffusion_volumes_results_fill_args(milonga.functions.pow));

  wasora_call(mesh_cell_indexes(milonga.mesh, milonga.groups));
  
  return WASORA_RUNTIME_OK;  
}

#undef  __FUNCT__
#define __FUNCT__ "diffusion_volumes_results_fill_args"
// esta rutina rellena datos administrativos de las funciones resultados
int diffusion_volumes_results_fill_args(function_t *function) {

  // estructurado o no, tenemos data
  function->data_size = milonga.spatial_unknowns;
  if (milonga.mesh->structured) {
    function->rectangular_mesh = 1;
    function->x_increases_first = 1;
    function->rectangular_mesh_size = milonga.mesh->rectangular_mesh_size;
    function->rectangular_mesh_point = milonga.mesh->rectangular_mesh_point;    
  }
  
  // pero tambien variables por si queremos hacer cuentitas
  function->var_argument = calloc(3, sizeof(var_t *));
  function->var_argument = wasora_mesh.vars.arr_x;

  function->data_argument = milonga.mesh->cells_argument;
  function->data_value = calloc(function->data_size, sizeof(double));
  
  // en volumes finitos decimos que la funcion es tipo mesh cell
  function->type = type_pointwise_mesh_cell;
  function->multidim_threshold = DEFAULT_MULTIDIM_INTERPOLATION_THRESHOLD;
  function->mesh = milonga.mesh;

  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "diffusion_volumes_matrices_build"
int diffusion_volumes_matrices_build(void) {

  int i, j;
  int g, g_prime;
  double xi;
  double w_ij, w_ji;

  double D_i;
  double D_j;
  xs_t *material_xs;
  xs_t *neighbor_xs;
  cell_t *cell;


  for (i = 0; i < milonga.mesh->n_cells; i++) {
    
    cell = &milonga.mesh->cell[i];
      
    if (cell->element->physical_entity == NULL || cell->element->physical_entity->material == NULL) {
      // esto pasa solo en malla estructuradas
      for (g = 0; g < milonga.groups; g++) {
        // phi = 0
        // TODO: no arruinar el radio espectral, aunque la petsc las escalea a traves
        // del precondicionador asi que no deberia haber mucho problema con esto
        MatSetValue(milonga.R, cell->index[g], cell->index[g], 1.0, ADD_VALUES);
      }    
    } else {
          
      material_xs = (xs_t *)(cell->element->physical_entity->material->ext);

      for (g = 0; g < milonga.groups; g++) {

        // ----- fuentes independientes  ----------------------------
        if (material_xs->S[g]->n_tokens != 0) {
          if ((xi = diffusion_volumes_cell_integral(cell, material_xs->S[g])) != 0) {
            VecSetValue(milonga.S, cell->index[g], xi, ADD_VALUES);
            milonga.has_sources = 1;
          }
        }

          // ----- absorcion total ----------------------------
        if (material_xs->SigmaT[g]->n_tokens != 0) {
          xi = diffusion_volumes_cell_integral(cell, material_xs->SigmaT[g]);
        } else {
          // podemos hacer esto porque la integracion es lineal thanks god!
          xi = diffusion_volumes_cell_integral(cell, material_xs->SigmaA[g]);
          for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
            xi += diffusion_volumes_cell_integral(cell, material_xs->SigmaS0[g][g_prime]);
          }
        }
        if (xi != 0) {
          MatSetValue(milonga.R, cell->index[g], cell->index[g], xi, ADD_VALUES);
        }


        // ----- fision ----------------------------
        if (gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), g) != 0) {
          for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
            xi = gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), g) * diffusion_volumes_cell_integral(cell, material_xs->nuSigmaF[g_prime]);
            if (isnan(xi)) {
              wasora_push_error_message("NaN found when computing fission term for group %d at cell %d (element %d)", g+1, cell->id, cell->element->id);
              return WASORA_RUNTIME_ERROR;
            } else if (xi < 0) {
              wasora_push_error_message("negative fission term for group %d at cell %d (element %d)", g+1, cell->id, cell->element->id);
              return WASORA_RUNTIME_ERROR;
            }
            if (xi != 0) {
              milonga.has_fission = 1;
              MatSetValue(milonga.F, cell->index[g], cell->index[g_prime], xi, INSERT_VALUES);  
            }
          }
        }

        // ----- scattering ----------------------------
        for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
          if ((xi = -diffusion_volumes_cell_integral(cell, material_xs->SigmaS0[g_prime][g])) != 0) {
            MatSetValue(milonga.R, cell->index[g], cell->index[g_prime], xi, ADD_VALUES);
          }
        }


        // ----- fugas ----------------------------
        wasora_value(wasora_mesh.vars.x) = cell->x[0];
        wasora_value(wasora_mesh.vars.y) = cell->x[1];
        wasora_value(wasora_mesh.vars.z) = cell->x[2];
        if (material_xs->D[g]->n_tokens != 0) {
          D_i = wasora_evaluate_expression(material_xs->D[g]);
        } else if (material_xs->SigmaT[g]->n_tokens != 0) {
          // nos fabricamos un coeficiente de difusion como 1/3*SigmaT (scattering isotropico)
          D_i = 1.0/(3.0 * wasora_evaluate_expression(material_xs->SigmaT[g]));
        } else {
          wasora_push_error_message("neither D nor SigmaT given for material '%s' (diffusion does not handle void)", cell->element->physical_entity->material->name);
          return WASORA_RUNTIME_ERROR;
        }

        for (j = 0; j < cell->n_neighbors; j++) {
          if (cell->neighbor[j].cell != NULL && cell->neighbor[j].cell->element->physical_entity != NULL) {
            // si el vecino tiene celda entonces representamos la ecuacion

            if (milonga.volhom == 0) {
              if (cell->element->physical_entity == cell->neighbor[j].cell->element->physical_entity) {

                // si son del mismo material se hace la cuenta a primeros vecinos
                // TODO: en verdad aca habria que poner coeficientes de continuidad
                // y siempre deberiamos hacer esta cuenta
                xi = D_i*cell->neighbor[j].S_ij
                       * mesh_subtract_dot(cell->neighbor[j].cell->x, cell->x, cell->neighbor[j].n_ij)
                       / mesh_subtract_squared_module(cell->neighbor[j].cell->x, cell->x);

              } else {

                // si no son del mismo material hay que hacer las cuentas para conservar la corriente
                if (cell->neighbor[j].cell->element->physical_entity == NULL) {
                  wasora_push_error_message("missing a physical entity for element %d'", cell->neighbor[j].cell->element->id);
                  return WASORA_RUNTIME_ERROR;
                }
                if (cell->neighbor[j].cell->element->physical_entity->material == NULL) {
                  wasora_push_error_message("missing a material for element %d", cell->neighbor[j].cell->element->id);
                  return WASORA_RUNTIME_ERROR;
                }
                if ((neighbor_xs = (xs_t *)cell->neighbor[j].cell->element->physical_entity->material->ext) == NULL) {
                  wasora_push_error_message("wrong XS set for element %d", cell->neighbor[j].cell->element->id);
                  return WASORA_RUNTIME_ERROR;
                }
                
                if (neighbor_xs->D[g]->n_tokens != 0) {
                  D_j = wasora_evaluate_expression(neighbor_xs->D[g]);
                } else if (neighbor_xs->SigmaT[g]->n_tokens != 0) {
                  // nos fabricamos un coeficiente de difusion como 1/3*SigmaT (scattering isotropico)
                  D_j = 1.0/(3.0 * wasora_evaluate_expression(neighbor_xs->SigmaT[g]));
                } else {
                  wasora_push_error_message("neither D nor SigmaT given for material '%s' (diffusion does not handle void)", cell->neighbor[j].cell->element->physical_entity->material->name);
                  return WASORA_RUNTIME_ERROR;
                }

                w_ij =  D_i * mesh_subtract_dot(cell->neighbor[j].x_ij, cell->x, cell->neighbor[j].n_ij)
                            / mesh_subtract_squared_module(cell->neighbor[j].x_ij, cell->x);

                w_ji = -D_j * mesh_subtract_dot(cell->neighbor[j].x_ij, cell->neighbor[j].cell->x, cell->neighbor[j].n_ij)
                            / mesh_subtract_squared_module(cell->neighbor[j].x_ij, cell->neighbor[j].cell->x);

                xi = D_i*cell->neighbor[j].S_ij
                       * (w_ji/(w_ij + w_ji))
                       * mesh_subtract_dot(cell->neighbor[j].x_ij, cell->x, cell->neighbor[j].n_ij)
                       / mesh_subtract_squared_module(cell->neighbor[j].x_ij, cell->x);
              }

            } else {
              // payasada para mostrar que pasa si haces homogeneo
              xi = D_i*cell->neighbor[j].S_ij
                      * mesh_subtract_dot(cell->neighbor[j].cell->x, cell->x, cell->neighbor[j].n_ij)
                       / mesh_subtract_squared_module(cell->neighbor[j].cell->x, cell->x);
            }

            if (isnan(xi)) {
              wasora_push_error_message("NaN found when computing leakage term for group %d at cell %d (element %d)", g+1, cell->id, cell->element->id);
              return WASORA_RUNTIME_ERROR;
            }
            if (xi != 0) {
              milonga.has_diffusion = 1;
              MatSetValue(milonga.R, cell->index[g], cell->index[g], +xi, ADD_VALUES);
              MatSetValue(milonga.R, cell->index[g], cell->neighbor[j].cell->index[g], -xi, ADD_VALUES);
            }

          } else {
            // sino es una condicion de contorno
            physical_entity_t *physical_entity;
            int bc_type = BC_UNDEFINED;
            expr_t *bc_args = NULL;
              
            if (milonga.mesh->structured == 0) {
              if (cell->neighbor[j].element == NULL ||
                  cell->neighbor[j].element->physical_entity == NULL ||
                  cell->neighbor[j].element->physical_entity->bc_type_int == BC_NULL) {
                bc_type = BC_NULL;
                bc_args = cell->neighbor[j].element->physical_entity->bc_args;
              } else  {
                bc_type = cell->neighbor[j].element->physical_entity->bc_type_int;
                bc_args = cell->neighbor[j].element->physical_entity->bc_args;
              }
            } else {
              LL_FOREACH(wasora_mesh.physical_entities, physical_entity) {
                if ((physical_entity->struct_bc_direction-1) == j) {
                  bc_type = physical_entity->bc_type_int;
                  bc_args = physical_entity->bc_args;
                }
              }
            }
            
            // default boundary condition
            if (bc_type == BC_UNDEFINED) {
              bc_type = BC_VACUUM;
            }

            
            switch (bc_type) {
              case BC_UNDEFINED:
                wasora_push_error_message("undefined boundary condition");
                return WASORA_RUNTIME_ERROR;
              case BC_NULL:
                // si es NULL entonces flujo a traves de la superficie es esto
                xi = D_i*cell->neighbor[j].S_ij
                       * mesh_subtract_dot(cell->neighbor[j].x_ij, cell->x, cell->neighbor[j].n_ij)
                       / mesh_subtract_squared_module(cell->neighbor[j].x_ij, cell->x);

                if (isnan(xi)) {
                  wasora_push_error_message("NaN found when computing boundary term for group %d at cell %d (element %d)", g+1, cell->id, cell->element->id);
                  return WASORA_RUNTIME_ERROR;
                }
                MatSetValue(milonga.R, cell->index[g], cell->index[g], +xi, ADD_VALUES);
              break;
              case BC_MIRROR: 
                // si es MIRROR, el gradiente es cero y no hacemos nada
                ;
              break;
              case BC_VACUUM:
                // si es vacio (robin) entonces el flujo es 
                // (el factor D_i aparece dividiendo en la expresion de la condicion de contorno)
                if (bc_args != NULL) {
                  xi = fabs(cell->neighbor[j].S_ij * wasora_evaluate_expression(bc_args));
                } else {
                  xi = cell->neighbor[j].S_ij * 0.5;
                }

                MatSetValue(milonga.R, cell->index[g], cell->index[g], xi, ADD_VALUES);
              break;
            }
          }
        }
      }
    }
  }
  
  if (milonga.has_fission == 0 && milonga.has_sources == 0) {
    wasora_push_error_message("sources (independent & fission) are identically zero through the domain");
    return WASORA_RUNTIME_ERROR;
  }
  if (milonga.has_diffusion == 0) {
    wasora_push_error_message("the domain does not contain any diffusive material");
    return WASORA_RUNTIME_ERROR;
  }
  

  // ensamblamos las matrices (para la boludina del mpi)
  wasora_call(milonga_assembly_objects(MAT_FINAL_ASSEMBLY));  

  return WASORA_RUNTIME_OK;
}



#undef  __FUNCT__
#define __FUNCT__ "diffusion_volumes_cell_integral"
double diffusion_volumes_cell_integral(cell_t *cell, expr_t *f) {

  if (f == NULL) {
    return 0;
  }

  // si todavia no calculamos el volumen de la celda, lo hacemos ahora
  if (cell->volume == 0) {
    cell->volume = cell->element->type->element_volume(cell->element);
  }

  wasora_value(wasora_mesh.vars.x) = cell->x[0];
  wasora_value(wasora_mesh.vars.y) = cell->x[1];
  wasora_value(wasora_mesh.vars.z) = cell->x[2];
  
  return wasora_evaluate_expression(f) * cell->volume;
  
}

#undef  __FUNCT__
#define __FUNCT__ "diffusion_volumes_normalize_flux"
int diffusion_volumes_normalize_flux(void) {

  int i, g;
  double factor;
  double num = 0;
  double den = 0;

  if (wasora_var(milonga.vars.power) == 0) {
    // calculamos el factor de normalizacion factor = num/den
    for (i = 0; i < milonga.mesh->n_cells; i++) {
      if (milonga.mesh->cell[i].volume == 0) {
        milonga.mesh->cell[i].volume = milonga.mesh->cell[i].element->type->element_volume(milonga.mesh->cell[i].element);
      }
      num += milonga.mesh->cell[i].volume;
      for (g = 0; g < milonga.groups; g++) {
        VecGetValues(milonga.phi, 1, &milonga.mesh->cell[i].index[g], &milonga.functions.phi[g]->data_value[i]);
        den += milonga.mesh->cell[i].volume * milonga.functions.phi[g]->data_value[i];
      }
    }

  } else {

    xs_t *xs;

    num = wasora_var(milonga.vars.power);
    for (i = 0; i < milonga.mesh->n_cells; i++) {
      if (milonga.mesh->cell[i].element->physical_entity != NULL &&
          milonga.mesh->cell[i].element->physical_entity->material != NULL &&
          (xs = (xs_t *)milonga.mesh->cell[i].element->physical_entity->material->ext) != NULL) {

        for (g = 0; g < milonga.groups; g++) {
          VecGetValues(milonga.phi, 1, &milonga.mesh->cell[i].index[g], &milonga.functions.phi[g]->data_value[i]);
          den += diffusion_volumes_cell_integral(&milonga.mesh->cell[i], xs->eSigmaF[g]) * milonga.functions.phi[g]->data_value[i];
        }
      }
    }

    if (den == 0) {
      wasora_push_error_message("power setpoint was given but eSigmaF is identically zero");
      return WASORA_RUNTIME_ERROR;
    }
  }


  factor = num/den;

  for (g = 0; g < milonga.groups; g++) {
    for (i = 0; i < milonga.spatial_unknowns; i++) {
      milonga.functions.phi[g]->data_value[i] *= factor;
    }
  }

  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "diffusion_volumes_results_fill_flux"
int diffusion_volumes_results_fill_flux(void) {

  int i, g;

  for (i = 0; i < milonga.mesh->n_cells; i++) {
    for (g = 0; g < milonga.groups; g++) {
      VecGetValues(milonga.phi, 1, &milonga.mesh->cell[i].index[g], &milonga.functions.phi[g]->data_value[i]);
    }
  }

  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "diffusion_volumes_results_fill_power"
// rellena los data_value de los resultados para el calculo no estructurado
int diffusion_volumes_results_fill_power(void) {

  int i, g;
  xs_t *xs;
  
  for (i = 0; i < milonga.mesh->n_cells; i++) {

    wasora_var(wasora_mesh.vars.x) = milonga.mesh->cell[i].x[0];
    wasora_var(wasora_mesh.vars.y) = milonga.mesh->cell[i].x[1];
    wasora_var(wasora_mesh.vars.z) = milonga.mesh->cell[i].x[2];
    
    milonga.functions.pow->data_value[i] = 0;

    if (milonga.mesh->cell[i].element->physical_entity != NULL &&
        milonga.mesh->cell[i].element->physical_entity->material != NULL &&
        (xs = (xs_t *)milonga.mesh->cell[i].element->physical_entity->material->ext) != NULL) {    

      for (g = 0; g < milonga.groups; g++) {
        milonga.functions.pow->data_value[i] += wasora_evaluate_expression(xs->eSigmaF[g]) * milonga.functions.phi[g]->data_value[i];
      }
    }
  }

  return WASORA_RUNTIME_OK;

}


#undef  __FUNCT__
#define __FUNCT__ "diffusion_volumes_problem_free"
int diffusion_volumes_problem_free(void) {
  
  int g;
  
  if (milonga.mesh != NULL && milonga.mesh->n_cells != 0) {
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

    mesh_free(milonga.mesh);
  }
   
  return WASORA_RUNTIME_OK;
}
