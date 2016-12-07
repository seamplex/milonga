/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga linear problem solution routines
 *
 *  Copyright (C) 2014--2015 jeremy theler
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

#include <petscksp.h>

#include "milonga.h"


#undef  __FUNCT__
#define __FUNCT__ "milonga_solve_linear_petsc"
int milonga_solve_linear_petsc(Mat R, Vec S, Vec phi) {

  KSPConvergedReason reason;
  
  // creamos un solver lineal
  petsc_call(KSPCreate(PETSC_COMM_WORLD, &milonga.ksp));
  petsc_call(KSPSetOperators(milonga.ksp, R, R));
/*  
  TODO: elegir tolerancias
  KSPSetTolerances(fino.solver, wasora_var(fino.rtol)/fino.N,
  		                          wasora_var(fino.atol)/fino.N,
  		                          wasora_var(fino.divtol)/fino.N,
  		                          (PetscInt)wasora_var(fino.max_it));
 */

  // precondicionarior por defecto
  petsc_call(KSPGetPC(milonga.ksp, &milonga.pc));
petsc_call(  PCSetType(milonga.pc, PCJACOBI));
  
  // el KSP
  if (milonga.ksp_type != NULL) {
    petsc_call(KSPSetType(milonga.ksp, milonga.ksp_type));
  }

  // el precondicionador
  if (milonga.pc_type != NULL) {
    petsc_call(PCSetType(milonga.pc, milonga.pc_type));
  }

  // sobreescribimos con la linea de comandos
  petsc_call(KSPSetFromOptions(milonga.ksp));
  
  // do the work!
  petsc_call(KSPSolve(milonga.ksp, S, phi));

  // chequeamos que haya convergido
  petsc_call(KSPGetConvergedReason(milonga.ksp, &reason));
  if (reason < 0) {
    wasora_push_error_message("PETSc's linear solver did not converge with reason '%s' (%d)", KSPConvergedReasons[reason], reason);
    return WASORA_RUNTIME_ERROR;
  }

  return WASORA_RUNTIME_OK;

}

