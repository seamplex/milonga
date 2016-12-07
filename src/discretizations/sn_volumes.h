/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's transport with finite volumes over structured grids routines
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

extern int sn_volumes_problem_init(void);
extern int sn_volumes_results_init_args(void);
extern int sn_volumes_results_fill_args(function_t *);

extern int sn_volumes_matrices_build(void);
extern int sn_volumes_structured_matrices_build(void);
extern int sn_volumes_unstructured_matrices_build(void);

extern int sn_volumes_results_fill_flux(void);
extern int sn_volumes_normalize_flux(void);
extern int sn_volumes_results_fill_power(void);
extern int sn_volumes_problem_free(void);

extern double sn_volumes_cell_integral(cell_t *, expr_t *);


extern int sn_init_weights(void);
