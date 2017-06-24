/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's discrete ordinates with finite volumes
 *
 *  Copyright (C) 2014--2016 jeremy theler
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
#include "sn_volumes.h"

#define flat_index(i,j,k) ((i) + (j)*milonga.mesh->ncells_x + (k)*milonga.mesh->ncells_x*milonga.mesh->ncells_y) 
#define dof_index(m,g) ((m)*milonga.groups + (g))
//#define dof_index(m,g) ((m) + (g)*milonga.directions)

#undef  __FUNCT__
#define __FUNCT__ "sn_volumes_problem_init"
int sn_volumes_problem_init(void) {

  int g, n;

  if (milonga.mesh == NULL) {
    wasora_push_error_message("no mesh found");
    return WASORA_RUNTIME_ERROR;
  }

//Check that the mesh is order 1 in finite volumes method.
  if(milonga.scheme == scheme_volumes && milonga.mesh->order > 1)
    {
    wasora_push_error_message("The finite volumes methods only accepts 1 order elements and your mesh has at least one element with order %d. Please remesh with order 1.\n", milonga.mesh->order);
    return WASORA_PARSER_ERROR;
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
  milonga_allocate_global_matrices(milonga.spatial_unknowns * milonga.directions * milonga.groups,
                                   milonga.mesh->max_faces_per_element + milonga.directions * milonga.groups,
                                   milonga.directions * milonga.groups);
  milonga_allocate_global_vectors();
  
  
  // inicializamos los pesos de las ordenadas discretas
  wasora_call(sn_init_weights());

  wasora_var(wasora_mesh.vars.cells) = (double)milonga.mesh->n_cells;
  wasora_var(wasora_mesh.vars.nodes) = (double)milonga.mesh->n_nodes;
  wasora_var(wasora_mesh.vars.elements) = (double)milonga.mesh->n_elements;
  milonga.mesh->data_type = data_type_element;

  for (n = 0; n < milonga.directions; n++) {
    for (g = 0; g < milonga.groups; g++) {
      wasora_call(sn_volumes_results_fill_args(milonga.functions.psi[n][g]));
    }
  }
  
  for (g = 0; g < milonga.groups; g++) {
    wasora_call(sn_volumes_results_fill_args(milonga.functions.phi[g]));
  }
  wasora_call(sn_volumes_results_fill_args(milonga.functions.pow));

  wasora_call(mesh_cell_indexes(milonga.mesh, milonga.groups * milonga.directions));
  
  return WASORA_RUNTIME_OK;  
}


#undef  __FUNCT__
#define __FUNCT__ "sn_volumes_results_fill_args"
// esta rutina rellena datos administrativos de las funciones resultados
int sn_volumes_results_fill_args(function_t *function) {

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

  // los argumentos los calculamos en results_init_structured_data_arg
  // y aca apuntamos cada funcion a ese array
  function->data_argument = milonga.mesh->cells_argument;
  function->data_value = calloc(function->data_size, sizeof(double));
  
  // en volumes finitos decimos que la funcion es tipo mesh cell
  function->type = type_pointwise_mesh_cell;
  function->multidim_threshold = DEFAULT_MULTIDIM_INTERPOLATION_THRESHOLD;
  function->mesh = milonga.mesh;

  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "sn_volumes_matrices_build"
int sn_volumes_matrices_build(void) {

  int i, j;
  int g, g_prime;
  int m, m_prime, m_refl;
  int p, p_prime;
  double reflected[3] = {0, 0, 0};
  double xi;
  xs_t *material_xs;
  cell_t *cell;
  neighbor_t *neighbor;

  double w_ij;
  double Omega_dot_outward;
  
  for (i = 0; i < milonga.mesh->n_cells; i++) {
    cell = &milonga.mesh->cell[i];
    
    if (cell->element->physical_entity != NULL && cell->element->physical_entity->material != NULL) {

      material_xs = (xs_t *)(cell->element->physical_entity->material->ext);

      for (m = 0; m < milonga.directions; m++) {
        for (g = 0; g < milonga.groups; g++) {

          p = dof_index(m, g);
          
          // ----- fuentes independientes  ----------------------------
          if (material_xs->S[g]->n_tokens != 0) {
            xi = sn_volumes_cell_integral(cell, material_xs->S[g]);
            if (xi != 0) {
              VecSetValue(milonga.S, cell->index[p], xi, ADD_VALUES);
              milonga.has_sources = 1;
            }
          }

          // ----- absorcion total ----------------------------
          if (material_xs->SigmaT[g]->n_tokens != 0) {
            xi = sn_volumes_cell_integral(cell, material_xs->SigmaT[g]);
          } else {
            // podemos hacer esto porque la integracion es lineal thanks god!
            xi = sn_volumes_cell_integral(cell, material_xs->SigmaA[g]);
            for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
              xi += sn_volumes_cell_integral(cell, material_xs->SigmaS0[g][g_prime]);
            }
          }
          if (xi != 0) {
            MatSetValue(milonga.R, cell->index[p], cell->index[p], xi, ADD_VALUES);
          }

          for (m_prime = 0; m_prime < milonga.directions; m_prime++) {
            for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
              
              p_prime = dof_index(m_prime, g_prime);
              
              // ----- fision ----------------------------
              if (gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), g) != 0) {
                xi = gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), g) * w[m_prime] * sn_volumes_cell_integral(cell, material_xs->nuSigmaF[g_prime]);
                if (xi != 0) {
                  milonga.has_fission = 1;
                  MatSetValue(milonga.F, cell->index[p], cell->index[p_prime], xi, ADD_VALUES);
                }
              }

              // ----- scattering ----------------------------
              xi = -w[m_prime] * sn_volumes_cell_integral(cell, material_xs->SigmaS0[g_prime][g]);
              // si tenemos scattering anisotropico, l = 1
              if (material_xs->SigmaS1[g_prime][g]->n_tokens != 0) {
                xi -= w[m_prime] * sn_volumes_cell_integral(cell, material_xs->SigmaS1[g_prime][g]) * 3.0 * mesh_dot(Omega[m], Omega[m_prime]);
              }
              if (xi != 0) {
                MatSetValue(milonga.R, cell->index[p], cell->index[p_prime], xi, ADD_VALUES);
              }
            }
          }

      
          // ----- fugas ----------------------------
          for (j = 0; j < cell->n_neighbors; j++) {

            neighbor = &cell->neighbor[j];
            Omega_dot_outward = mesh_dot(Omega[m], neighbor->n_ij);
            xi = Omega_dot_outward * neighbor->S_ij;
            
            if (neighbor->cell != NULL && neighbor->cell->element->physical_entity != NULL) {
              // tiene un vecino volumetrico
              // peso geometrico
              w_ij = mesh_subtract_module(cell->neighbor[j].x_ij, cell->neighbor[j].cell->x) / (mesh_subtract_module(cell->neighbor[j].x_ij, cell->x) + mesh_subtract_module(cell->neighbor[j].x_ij, cell->neighbor[j].cell->x));
              // peso estabilizado
              if (Omega_dot_outward > 0) {
                w_ij += wasora_value(milonga.vars.sn_alpha)*(1-w_ij);
              } else {
                w_ij -= wasora_value(milonga.vars.sn_alpha)*w_ij;
              }
              
              MatSetValue(milonga.R, cell->index[p], cell->index[p], w_ij*xi, ADD_VALUES);
              MatSetValue(milonga.R, cell->index[p], neighbor->cell->index[dof_index(m,g)], (1-w_ij)*xi, ADD_VALUES);
              
            } else {
              // condicion de contorno
              physical_entity_t *physical_entity;
              int bc_type = BC_UNDEFINED;
                            
              if (milonga.mesh->structured == 0) {
                if (cell->neighbor[j].element == NULL ||
                    cell->neighbor[j].element->physical_entity == NULL ||
                    cell->neighbor[j].element->physical_entity->bc_type_phys == BC_VACUUM ||
                    cell->neighbor[j].element->physical_entity->bc_type_phys == BC_NULL) {
                  bc_type = BC_VACUUM;
                } else if (cell->neighbor[j].element->physical_entity->bc_type_phys == BC_MIRROR) {
                  bc_type = BC_MIRROR;
                }
              } else {
                LL_FOREACH(wasora_mesh.physical_entities, physical_entity) {
                  if ((physical_entity->struct_bc_direction-1) == j) {
                    bc_type = physical_entity->bc_type_phys;
                  }
                }
              }
              
              // default boundary condition
              if (bc_type == BC_UNDEFINED) {
                bc_type = BC_VACUUM;
              }
              
              if (bc_type == BC_VACUUM) {
                if (Omega_dot_outward < 0) {
                  // corriente entrante nula, el flujo en S_ij es cero
                  ;
                } else {
                  // corriente saliente, no la sabemos pero suponemos que es igual la de la celda
                  MatSetValue(milonga.R, cell->index[p], cell->index[p], xi, ADD_VALUES);
                }
              } else if (bc_type == BC_MIRROR) {

                if (Omega_dot_outward < 0) {
                  // si el producto interno de Omega con la normal es positivo entonces tenemos que reflejar
                  // si Omega es la direccion de incidencia, la direccion reflejada con respecto a la normal outward_normal es 
                  // reflected = Omega - 2*(Omega dot outward_normal) * outward_normal
                  reflected[0] = Omega[m][0] - 2*Omega_dot_outward * neighbor->n_ij[0];
                  reflected[1] = Omega[m][1] - 2*Omega_dot_outward * neighbor->n_ij[1];
                  reflected[2] = Omega[m][2] - 2*Omega_dot_outward * neighbor->n_ij[2];
                  for (m_refl = 0; m_refl < milonga.directions; m_refl++) {
                    if (fabs(reflected[0]-Omega[m_refl][0]) < wasora_var(wasora_mesh.vars.eps) &&
                        fabs(reflected[1]-Omega[m_refl][1]) < wasora_var(wasora_mesh.vars.eps) &&
                        fabs(reflected[2]-Omega[m_refl][2]) < wasora_var(wasora_mesh.vars.eps)) {
                      break;
                    }
                  }
                  if (m_refl == milonga.directions) {
                    wasora_push_error_message("cannot find a reflected direction for n=%d (%.3f %.3f %.3f) in cell %d (normal %.3f %.3f %.3f)",
                            m, milonga.mesh->cell[i].id, Omega[m][0], Omega[m][1], Omega[m][2], neighbor->n_ij[0], neighbor->n_ij[1], neighbor->n_ij[2]);
                    return WASORA_RUNTIME_ERROR;
                  }                  

                  // la corriente entrante es igual a la reflejada
                  MatSetValue(milonga.R, cell->index[p], cell->index[dof_index(m_refl,g)], xi, ADD_VALUES);
                   
                } else {
                  // corriente saliente, no la sabemos pero suponemos que es igual la de la celda
                  MatSetValue(milonga.R, cell->index[p], cell->index[p], xi, ADD_VALUES);
                }
              }
            }
          }
        }
      }
    }
  }
     

  // TODO: dejar que las fuentes sean cero y poner condiciones de contorno no homogeneas
  if (milonga.has_fission == 0 && milonga.has_sources == 0) {
    wasora_push_error_message("sources (independent & fission) are identically zero through the domain (do you want multigroup and forgot to give MILONGA_PROBLEM GROUPS?)");
    return WASORA_RUNTIME_ERROR;
  }

  // ensamblamos las matrices (para la boludina del mpi)
  wasora_call(milonga_assembly_objects(MAT_FINAL_ASSEMBLY));  

  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "sn_volumes_results_fill_flux"
int sn_volumes_results_fill_flux(void) {

  int i, g, n;

  for (i = 0; i < milonga.mesh->n_cells; i++) {
    for (g = 0; g < milonga.groups; g++) {
      milonga.functions.phi[g]->data_value[i] = 0;
      for (n = 0; n < milonga.directions; n++) {
        VecGetValues(milonga.phi, 1, &milonga.mesh->cell[i].index[dof_index(n,g)], &milonga.functions.psi[n][g]->data_value[i]);
        milonga.functions.phi[g]->data_value[i] += w[n] * milonga.functions.psi[n][g]->data_value[i];
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "sn_volumes_normalize_flux"
int sn_volumes_normalize_flux(void) {

  int i, g, n;
  double factor;
  double num = 0;
  double den = 0;

  if (wasora_var(milonga.vars.power) == 0) {
    // calculamos el factor de normalizacion factor = num/den
    for (i = 0; i < milonga.mesh->n_cells; i++) {
      num += milonga.mesh->cell[i].volume;
      for (g = 0; g < milonga.groups; g++) {
        den += milonga.mesh->cell[i].volume * milonga.functions.phi[g]->data_value[i];
      }
    }

  } else {

    xs_t *xs;

    num = wasora_var(milonga.vars.power);
    for (i = 0; i < milonga.mesh->n_cells; i++) {
      if (milonga.mesh->cell[i].element->physical_entity != NULL && (xs = (xs_t *)milonga.mesh->cell[i].element->physical_entity->material->ext) != NULL) {
        for (g = 0; g < milonga.groups; g++) {
          milonga.functions.phi[g]->data_value[i] = 0;
          for (n = 0; n < milonga.directions; n++) {
            VecGetValues(milonga.phi, 1, &milonga.mesh->cell[i].index[g], &milonga.functions.psi[n][g]->data_value[i]);
            milonga.functions.phi[g]->data_value[i] += w[n] * milonga.functions.psi[n][g]->data_value[i];
          }
          den += sn_volumes_cell_integral(&milonga.mesh->cell[i], xs->eSigmaF[g]) * milonga.functions.phi[g]->data_value[i];
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

  for (n = 0; n < milonga.directions; n++) {
    for (g = 0; g < milonga.groups; g++) {
      for (i = 0; i < milonga.spatial_unknowns; i++) {
        milonga.functions.psi[n][g]->data_value[i] *= factor;
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "sn_volumes_results_fill_power"
int sn_volumes_results_fill_power(void) {

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
#define __FUNCT__ "sn_volumes_problem_free"
int sn_volumes_problem_free(void) {
  
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

#undef  __FUNCT__
#define __FUNCT__ "sn_volumes_cell_integral"
double sn_volumes_cell_integral(cell_t *cell, expr_t *f) {

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
