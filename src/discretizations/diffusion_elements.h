/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's diffusion with finite elements over unstructured grids routines
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

extern int diffusion_elements_problem_init(void);
extern int diffusion_elements_results_init_args(void);
extern int diffusion_elements_results_fill_args(function_t *);
extern int diffusion_elements_matrices_build(void);
extern int diffusion_elements_results_fill_flux(void);
extern int diffusion_elements_normalize_flux(void);
extern int diffusion_elements_results_fill_power(void);
extern int diffusion_elements_problem_free(void);


extern int diffusion_elements_allocate_elemental_objects(void);
extern int diffusion_elements_build_volume_objects(element_t *);
extern int diffusion_elements_build_robin_objects(element_t *, expr_t *);
extern int diffusion_elements_set_essential_bc(void);
extern int diffusion_elements_free_elemental_objects(void);
extern void diffusion_elements_compute_H(void);
extern void diffusion_elements_compute_B(void);

