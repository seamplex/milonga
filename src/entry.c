/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's definition of methods entry points
 *
 *  Copyright (C) 2014 jeremy theler
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

// function prototypes
#include "discretizations/diffusion_volumes.h"
#include "discretizations/diffusion_elements.h"
#include "discretizations/sn_volumes.h"
#include "discretizations/sn_elements.h"


#undef  __FUNCT__
#define __FUNCT__ "milonga_set_entry_points"
int milonga_set_entry_points(void) {
  
  if (milonga.formulation == formulation_diffusion && milonga.scheme == scheme_volumes) {
    // difusion con volumenes finitos
    milonga.problem_init =       diffusion_volumes_problem_init;
    milonga.results_fill_args =  diffusion_volumes_results_fill_args;
    milonga.matrices_build =     diffusion_volumes_matrices_build;
    milonga.results_fill_flux =  diffusion_volumes_results_fill_flux;
    milonga.normalize_flux =     diffusion_volumes_normalize_flux;
    milonga.results_fill_power = diffusion_volumes_results_fill_power;
    milonga.problem_free =       diffusion_volumes_problem_free;

  } else if (milonga.formulation == formulation_diffusion && milonga.scheme == scheme_elements) {
    // difusion con elementos finitos
    milonga.problem_init =       diffusion_elements_problem_init;
    milonga.results_fill_args =  diffusion_elements_results_fill_args;
    milonga.matrices_build =     diffusion_elements_matrices_build;
    milonga.results_fill_flux =  diffusion_elements_results_fill_flux;
    milonga.normalize_flux =     diffusion_elements_normalize_flux;
    milonga.results_fill_power = diffusion_elements_results_fill_power;
    milonga.problem_free =       diffusion_elements_problem_free;

  } else if (milonga.formulation == formulation_sn && milonga.scheme == scheme_volumes) {
    // ordenandas discretas con volumenes finitos
    milonga.problem_init =       sn_volumes_problem_init;
    milonga.results_fill_args =  sn_volumes_results_fill_args;
    milonga.matrices_build =     sn_volumes_matrices_build;
    milonga.results_fill_flux =  sn_volumes_results_fill_flux;
    milonga.normalize_flux =     sn_volumes_normalize_flux;
    milonga.results_fill_power = sn_volumes_results_fill_power;
    milonga.problem_free =       sn_volumes_problem_free;

  } else if (milonga.formulation == formulation_sn && milonga.scheme == scheme_elements) {
    // ordenandas discretas con volumenes finitos
    milonga.problem_init =       sn_elements_problem_init;
    milonga.results_fill_args =  sn_elements_results_fill_args;
    milonga.matrices_build =     sn_elements_matrices_build;
    milonga.results_fill_flux =  sn_elements_results_fill_flux;
    milonga.normalize_flux =     sn_elements_normalize_flux;
    milonga.results_fill_power = sn_elements_results_fill_power;
    milonga.problem_free =       sn_elements_problem_free;

  }

  if (milonga.problem_init == NULL) {
    wasora_push_error_message("no suitable method found to satisfy MILONGA_PROBLEM requests :( wanna code it?");
    return WASORA_RUNTIME_ERROR;
  }
  
  return WASORA_RUNTIME_OK;
}
