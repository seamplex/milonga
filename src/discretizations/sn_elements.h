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

// direcciones y pesos para ordenadas discretas
extern double **Omega;
extern double *w;

extern int sn_elements_problem_init(void);
extern int sn_elements_results_init_args(void);
extern int sn_elements_results_fill_args(function_t *);
extern int sn_elements_matrices_build(void);
extern int sn_elements_results_fill_flux(void);
extern int sn_elements_normalize_flux(void);
extern int sn_elements_results_fill_power(void);
extern int sn_elements_problem_free(void);


extern int sn_elements_allocate_general_elemental_objects(void);
extern int sn_elements_allocate_particular_elemental_objects(element_t *);
extern int sn_elements_build_elemental_objects(element_t *);
extern int sn_elements_set_essential_bc(void);
extern int sn_elements_free_elemental_objects(void);
extern void sn_elements_compute_H(void);
extern void sn_elements_compute_B(void);

extern int sn_elements_compute_outward_normal(element_t *, double *);
extern int sn_elements_add_vacuum_bc(node_t *, double *);
extern int sn_elements_add_mirror_bc(node_t *, double *);

extern int sn_init_weights(void);
