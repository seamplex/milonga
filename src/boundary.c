/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's boundary conditions routines
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

#include "milonga.h"

#undef  __FUNCT__
#define __FUNCT__ "milonga_boundaries"

int milonga_read_boundaries(void) {

  physical_entity_t *physical_entity;
  bc_t *bc;
  
  // barremos los physical entities y mapeamos cadenas a valores enteros
  // si alguna physical entity se llama "mirror" o "robin" ya le ponemos ese tipo de CC
  // por default, ponemos "vacuum" (discretization-dependent!)
  for (physical_entity = wasora_mesh.main_mesh->physical_entities; physical_entity != NULL; physical_entity = physical_entity->hh.next) {
    
//    LL_FOREACH(physical_entity->bcs, bc) {
    if ((bc = physical_entity->bcs) != NULL) {
    
      if (physical_entity->material != NULL) {
        wasora_push_error_message("both MATERIAL and BC given for physical entity '%s", physical_entity->name);
        return WASORA_PARSER_ERROR;
      }

      if (strcasecmp(bc->string, "null") == 0 || strcasecmp(bc->string, "dirichlet") == 0) {
        bc->type_phys = BC_NULL;

      } else if (strcasecmp(bc->string, "vacuum") == 0 || strcasecmp(bc->string, "robin") == 0) {
        bc->type_phys = BC_VACUUM;

      } else if (strcasecmp(bc->string, "mirror") == 0 || strcasecmp(bc->string, "neumann") == 0) {
        bc->type_phys = BC_MIRROR;
        
        if (bc->next != NULL) {
          bc->expr = calloc(1, sizeof(expr_t));
          wasora_parse_expression(bc->next->string, bc->expr);
        }

      } else {
        if (bc->type_phys == BC_UNDEFINED) {
          if (milonga.implicit_bc_none) {
            wasora_push_error_message("unknown boundary condition '%s' for physical entity '%s' (IMPLICIT_BC is set to NONE)", bc->string, physical_entity->name);
            return WASORA_RUNTIME_ERROR;
          } else {
            bc->type_phys = BC_VACUUM;   // default is vacuum
          }
        }
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}
