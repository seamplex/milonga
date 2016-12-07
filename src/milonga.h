/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga main header
 *
 *  Copyright (C) 2009--2015 jeremy theler
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

#ifndef _MILONGA_H_
#define _MILONGA_H_


#include <wasora.h>
#include "mesh/mesh.h"

#include <petscsys.h>
#include <petscpc.h>
#include <petsctime.h>

#include <slepceps.h>

#if defined(PETSC_USE_COMPLEX)
#error "PETSc should be compiled with real scalars to be used with milonga."
#endif

PetscErrorCode petsc_err;
#define petsc_call(s) {petsc_err = s; CHKERRQ(petsc_err);}


#define DEBUG_MATRICES_ASCII                       1
#define DEBUG_MATRICES_ASCII_STRUCT                2
#define DEBUG_MATRICES_PETSC_BINARY                4
#define DEBUG_MATRICES_PETSC_COMPRESSED_BINARY     8
#define DEBUG_MATRICES_PETSC_ASCII                16
#define DEBUG_MATRICES_PETSC_OCTAVE               32
#define DEBUG_MATRICES_PETSC_DENSE                64
#define DEBUG_MATRICES_X                         128
#define DEBUG_MATRICES_SNG                       256
#define DEBUG_MATRICES_SNG_STRUCT                512

#define DEFAULT_MATRICES_SIZE       1024

#define POST_INCLUDE_FLUX           16
#define POST_INCLUDE_XS             32

#define BC_UNDEFINED  0    //    diffusion                      sn
#define BC_NULL       1    //   null flux                 null outward flux
#define BC_VACUUM     2    // robin with 0.5D             null outward flux
#define BC_MIRROR     3    // zero normal derivative      equal reflectected fluxes

// forward definitions
typedef struct xs_t xs_t;
typedef struct debug_t debug_t;
typedef struct milonga_step_t milonga_step_t;
typedef struct milonga_times_t milonga_times_t;


struct {
  
  int initialized;
  
  int spatial_unknowns;  // cant de incognitas espaciales (= celdas o nodos)
  int groups;
  int dimensions;
  int SN;               // orden de la formulacion de transporte (SN)
  int directions;       // cantidad de direcciones 
  
  int implicit_bc_none; // flag que si es true hay que poner CC explicita a cada superficie externa

  
  int problem_size;
  
  mesh_t *mesh;         // la malla del problema
  
// las matrices del problema en formato de PETSc
  Mat F;
  Mat R;

// el vector de fuentes independientes
  Vec S;

// el flujo (la solucion con petsc!)
  Vec phi;

// vector guess inicial
  Vec guess;

// informacion para la petsc
  PetscInt rank;
  PetscInt size;
  
  EPS eps;      // contexto eigensolver (SLEPc)
  ST st;        // contexto de la transformacion espectral asociada
  KSP ksp;      // contexto del solver lineal asociado a la transformacion espectral
  PC pc;        // contexto del precondicionador

  struct {
    var_t *keff;
    var_t *rel_tolerance;
    var_t *power;
    
    var_t *sn_alpha;     // upwinding

    var_t *unknowns;

    var_t *residual_norm;
    var_t *error_estimate;
    var_t *rel_error;

    var_t *time_wall_ini;
    var_t *time_wall_build;
    var_t *time_wall_solve;
    var_t *time_wall_total;

    var_t *time_cpu_ini;
    var_t *time_cpu_build;
    var_t *time_cpu_solve;
    var_t *time_cpu_total;

    var_t *time_petsc_ini;
    var_t *time_petsc_build;
    var_t *time_petsc_solve;
    var_t *time_petsc_total;

    var_t *flops_petsc;
    
    var_t *available_memory;
    var_t *memory_usage_global;
    var_t *memory_usage_petsc;
    
  } vars;

  struct {
    // espectro de fision
    vector_t *chi;
  } vectors;

  struct {
    function_t **phi;          // flujos escalares
    function_t ***psi;      // flujos en las ordenadas discretas
    function_t *pow;
  } functions;
  
  expr_t xs_zero;   // expresion que tiene "0" para xs incrementales
  xs_t *material_xss;
  debug_t *debugs;
  
  // esto capaz que deberia ir en otro lado
  PetscClassId petsc_classid;

  PetscLogStage petsc_stage_init;
  PetscLogStage petsc_stage_build;
  PetscLogStage petsc_stage_solve;
  
  PetscLogEvent petsc_event_init;
  PetscLogEvent petsc_event_build;
  PetscLogEvent petsc_event_solve;
  
  PetscLogDouble petsc_flops_init;
  PetscLogDouble petsc_flops_build;
  PetscLogDouble petsc_flops_solve;
  
  

  enum {
    scheme_volumes,
    scheme_elements,
  } scheme;

  enum {
    formulation_diffusion,
    formulation_sn,
  } formulation;
  
  loadable_routine_t *user_provided_eigensolver;
  loadable_routine_t *user_provided_linearsolver;

  char *eps_type;
  char *st_type;
  char *ksp_type;
  char *pc_type;
  
  expr_t eps_ncv;
  expr_t st_shift;
  expr_t st_anti_shift;

  enum {
    spectrum_largest_eigenvalue,              // F phi =  k  F phi
    spectrum_smallest_eigenvalue,             // R phi = 1/k F phi
  } spectrum;

  // flag academico para ver que pasa si hacemos volumenes mal
  int volhom;
  
  // flags que indican que es lo que hay
  int has_fission;
  int has_sources;
  int has_diffusion;

  // apuntadores a las funciones adecuadas segun elementos/volumenes
  int (*problem_init)(void);
  int (*results_fill_args)(function_t *);
  int (*matrices_build)(void);
  int (*results_fill_flux)(void);
  int (*normalize_flux)(void);
  int (*results_fill_power)(void);
  int (*problem_free)(void);

} milonga;

struct milonga_step_t {
  int do_not_build;
  int do_not_solve;
};

// para medir tiempos (wall y cpu)
struct milonga_times_t {
  PetscLogDouble init_begin;
  PetscLogDouble init_end;
  PetscLogDouble build_begin;
  PetscLogDouble build_end;
  PetscLogDouble solve_begin;
  PetscLogDouble solve_end;
};

// propiedades nucleares de un material:
// las XS son expresiones que pueden depender de variables (escalares)
// o de parametros (distribuciones dentro de la geometria)
struct xs_t {
  // arreglos de apuntadores a expresiones
  expr_t **D;
  expr_t **nuSigmaF;
// la sigma total, la de absorcion y las de scattering no son independientes!
  expr_t **SigmaT;
  expr_t **SigmaA;
  expr_t ***SigmaS0;
  expr_t ***SigmaS1;
  expr_t **eSigmaF;
  
  // la fuente de neutrones independiente
  expr_t **S;

  xs_t *next;
};

struct debug_t {
  file_t *file;
  PetscViewer viewer;
  
  int matrices;   // mascara con flags de que hay que exportar
  expr_t matrices_size;
  expr_t matrices_stride;
  
  int include_input;

  int file_opened;
  
  debug_t *next;
};

struct {
  int n_processors;

  char **model_name;
  char **mhz;
  char **cache_size;
  char **bogomips;
} cpuinfo;

// milonga.c
extern void milonga_resolve_xs_expr(material_t *, char *, expr_t **, int, int);
extern int milonga_assembly_objects(MatAssemblyType);

// boundary.c
extern int milonga_read_boundaries(void);

// init.c
extern int plugin_init_before_parser(void);
extern int plugin_init_after_parser(void);
extern int plugin_init_before_run(void);
extern int plugin_finalize(void);

extern int milonga_problem_init(void);
extern int milonga_instruction_step(void *);

// parser.c
extern int milonga_parse_line(char *);
extern int milonga_define_result_functions(void);

// allocate.c
extern int milonga_allocate_global_objects(int, int, int);
extern int milonga_free_global_objects(void);

// entry.c
extern int milonga_set_entry_points(void);

// eigen_slepc.c
extern int milonga_solve_eigen_slepc(Mat, Mat, Vec, PetscScalar *);

// linear_petsc.c
extern int milonga_solve_linear_petsc(Mat, Vec, Vec);


// debug.c
extern int milonga_debug_n_processors(void);
extern int milonga_debug_cpu_info(void);
extern int milonga_debug_cpu_info_free(void);
extern int milonga_debug_open(debug_t *);
extern int milonga_debug_initial(debug_t *);
extern int milonga_instruction_debug(void *);
extern int milonga_debug_close(debug_t *);
extern int milonga_print_petsc_vector(Vec, PetscViewer);
extern int milonga_print_petsc_matrix(Mat, PetscViewer);
extern int milonga_print_petsc_matrix_struct(Mat, PetscViewer);

// petschandler.c
PetscErrorCode milonga_handler(MPI_Comm comm, int, const char *, const char *, PetscErrorCode, PetscErrorType, const char *, void *);

// output.c
extern int milonga_instruction_post(void *);

// version.c
extern void milonga_usage(char *);
extern void milonga_version(FILE *, int, int);
extern void milonga_license(FILE *);

// times.c
extern double milonga_get_cpu_time(void);

extern const char *plugin_name(void);
extern const char *plugin_version(void);
extern const char *plugin_description(void);
extern const char *plugin_longversion(void);
extern const char *plugin_copyright(void);

#endif  /* _MILONGA_H_ */
