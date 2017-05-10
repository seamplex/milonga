/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's allocation routines
 *
 *  Copyright (C) 2010--2015 jeremy theler
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

#include <math.h>
#include <gsl/gsl_math.h>

#include <slepceps.h>

#include "milonga.h"

// estas estan separadas asi en problemas iterativos en lugar de hacer un
// matsetzero, lo que rellena todos los espacios prealocados, la destruimos
// y la volvemos a alocar (capaz que no sea lo optimo, pero funciona)

#undef  __FUNCT__
#define __FUNCT__ "milonga_allocate_global_matrices"
int milonga_allocate_global_matrices(int problem_size, int widthR, int widthF) {

  if (problem_size != 0) {
    milonga.problem_size = problem_size;
  }
  if (widthR != 0) {
    milonga.widthR = widthR;
  }
  if (widthF != 0) {
    milonga.widthF = widthF;
  }

  // matrices de remociones y de fisiones
  petsc_call(MatCreate(PETSC_COMM_WORLD, &milonga.R));
  petsc_call(MatSetSizes(milonga.R, PETSC_DECIDE, PETSC_DECIDE, milonga.problem_size, milonga.problem_size));
  petsc_call(MatSetFromOptions(milonga.R));
  petsc_call(MatMPIAIJSetPreallocation(milonga.R, milonga.widthR, PETSC_NULL, milonga.widthR, PETSC_NULL));
  petsc_call(MatSeqAIJSetPreallocation(milonga.R, milonga.widthR, PETSC_NULL));
  petsc_call(MatSetOption(milonga.R, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
  petsc_call(PetscObjectSetName((PetscObject)milonga.R, "R"));
  
  petsc_call(MatCreate(PETSC_COMM_WORLD, &milonga.F));
  petsc_call(MatSetSizes(milonga.F, PETSC_DECIDE, PETSC_DECIDE, milonga.problem_size, milonga.problem_size));
  petsc_call(MatSetFromOptions(milonga.F));
  petsc_call(MatMPIAIJSetPreallocation(milonga.F, milonga.widthF, PETSC_NULL, milonga.widthF, PETSC_NULL));
  petsc_call(MatSeqAIJSetPreallocation(milonga.F, milonga.widthF, PETSC_NULL));
  petsc_call(MatSetOption(milonga.F, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
  petsc_call(PetscObjectSetName((PetscObject)milonga.F, "F"));
  
  return WASORA_RUNTIME_OK;
  
}

#undef  __FUNCT__
#define __FUNCT__ "milonga_allocate_global_vectors"
int milonga_allocate_global_vectors(void) {

  // el autovector que va a tener el flujo
  petsc_call(MatCreateVecs(milonga.R, NULL, &milonga.phi));
  petsc_call(PetscObjectSetName((PetscObject)milonga.phi, "phi"));

  // el vector con las fuentes independientes
  petsc_call(MatCreateVecs(milonga.R, NULL, &milonga.S));
  petsc_call(PetscObjectSetName((PetscObject)milonga.S, "S"));
  
  // el vector de guess inicial (todos 1 ahora pero despues deberia tener la ultima solucion encontrada)
  // creamos el vector de guess inicial de tamanio N
  petsc_call(MatCreateVecs(milonga.R, NULL, &milonga.guess));
  petsc_call(PetscObjectSetName((PetscObject)milonga.guess, "phi0"));

  // hacemos un guess con todos unos
  petsc_call(VecSet(milonga.phi, 1.0));

  return WASORA_RUNTIME_OK;
  
}


#undef  __FUNCT__
#define __FUNCT__ "milonga_free_global_matrices"
int milonga_free_global_matrices(void) {

  if (milonga.R != PETSC_NULL) {
    petsc_call(MatDestroy(&milonga.R));
  }
  if (milonga.F != PETSC_NULL) {
    petsc_call(MatDestroy(&milonga.F));
  }
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "milonga_free_global_vectors"
int milonga_free_global_vectors(void) {

  if (milonga.phi != PETSC_NULL) {
    petsc_call(VecDestroy(&milonga.phi));
  }
  if (milonga.guess != PETSC_NULL) {
    petsc_call(VecDestroy(&milonga.guess));
  }
  if (milonga.ksp != PETSC_NULL) {
    petsc_call(EPSDestroy(&milonga.eps));
  }
  if (milonga.eps != PETSC_NULL) {
    petsc_call(EPSDestroy(&milonga.eps));
  }
  
  return WASORA_RUNTIME_OK;

}
