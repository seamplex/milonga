/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's diffusion with finite volumes over unstructured grids protoypes
 *
 *  Copyright (C) 2015 jeremy theler
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
         
extern int diffusion_volumes_problem_init(void);
extern int diffusion_volumes_results_init_args(void);
extern int diffusion_volumes_results_fill_args(function_t *);

extern int diffusion_volumes_matrices_build(void);
extern int diffusion_volumes_structured_matrices_build(void);
extern int diffusion_volumes_unstructured_matrices_build(void);

extern int diffusion_volumes_results_fill_flux(void);
extern int diffusion_volumes_normalize_flux(void);
extern int diffusion_volumes_results_fill_power(void);
extern int diffusion_volumes_problem_free(void);

extern double diffusion_volumes_cell_integral(cell_t *, expr_t *);
extern void diffusion_volumes_allocate_functions(void);
extern void diffusion_volumes_free_functions(void);
