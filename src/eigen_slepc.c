/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga eigenvalue problem solution routines
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

#include <slepceps.h>
#include <petscksp.h>

#include "milonga.h"


#undef  __FUNCT__
#define __FUNCT__ "milonga_solve_eigen_slepc"
int milonga_solve_eigen_slepc(Mat R, Mat F, Vec phi, PetscScalar *keff) {

  PetscReal xi = 1.0;
  PetscReal lambda = 1.0;
  PetscInt nconv;

  // creamos el contexto del eigensolver
  if (milonga.eps != NULL) {
    petsc_call(EPSDestroy(&milonga.eps));
  }
  petsc_call(EPSCreate(PETSC_COMM_WORLD, &milonga.eps));
  petsc_call(PetscObjectSetName((PetscObject)milonga.eps, "eigenvalue-problem_solver"));
 
  // y obtenemos los contextos asociados
  petsc_call(EPSGetST(milonga.eps, &milonga.st));
  petsc_call(PetscObjectSetName((PetscObject)milonga.st, "spectral_transformation"));
  petsc_call(STGetKSP(milonga.st, &milonga.ksp));
  petsc_call(PetscObjectSetName((PetscObject)milonga.ksp, "linear_solver"));
  petsc_call(KSPGetPC(milonga.ksp, &milonga.pc));
  petsc_call(PetscObjectSetName((PetscObject)milonga.pc, "preconditioner"));

  // le decimos que matrices tiene que resolver
  // depende de si resolvemos el problema que da k o 1/k
  switch (milonga.spectrum) {
    case spectrum_largest_eigenvalue:
      petsc_call(EPSSetOperators(milonga.eps, F, R));
      petsc_call(EPSSetWhichEigenpairs(milonga.eps, EPS_LARGEST_MAGNITUDE));
    break;
    case spectrum_smallest_eigenvalue:
      petsc_call(EPSSetOperators(milonga.eps, R, F));
      petsc_call(EPSSetWhichEigenpairs(milonga.eps, EPS_SMALLEST_MAGNITUDE));
    break;
  }

  // problema generalizado no hermitico
  petsc_call(EPSSetProblemType(milonga.eps, EPS_GNHEP));
  
  // TODO: ver bien esto del guess inicial
  petsc_call(EPSSetInitialSpace(milonga.eps, 1, &milonga.guess));

  // elegimos el metodo de solucion del eps
  if (milonga.eps_type != NULL) {
    petsc_call(EPSSetType(milonga.eps, milonga.eps_type));
  } else {
    petsc_call(EPSSetType(milonga.eps, "jd"));
  }
  
  // la transformada espectral
  if (milonga.st_type != NULL) {
    petsc_call(STSetType(milonga.st, milonga.st_type));
  }
  // si no esta seteado el tipo se queja la SLEPc
  if (milonga.st_shift.n_tokens != 0) {
    petsc_call(STSetShift(milonga.st, wasora_evaluate_expression(&milonga.st_shift)));
  }
  if (milonga.st_anti_shift.n_tokens != 0) {
    petsc_call(STCayleySetAntishift(milonga.st, wasora_evaluate_expression(&milonga.st_anti_shift)));
  }

  // el KSP
  if (milonga.ksp_type != NULL) {
    petsc_call(KSPSetType(milonga.ksp, milonga.ksp_type));
  }

  // el precondicionador
  if (milonga.pc_type != NULL) {
    petsc_call(PCSetType(milonga.pc, milonga.pc_type));
  } else {
    petsc_call(PCSetType(milonga.pc, "asm"));
  }

  // convergencia con respecto a la norma de las matrices
  petsc_call(EPSSetConvergenceTest(milonga.eps, EPS_CONV_NORM));
  
  // tolerancia
  if (wasora_var(milonga.vars.rel_tolerance) != 0) {
    petsc_call(EPSSetTolerances(milonga.eps, wasora_var(milonga.vars.rel_tolerance), PETSC_DECIDE));
  }

  // dimension del sub espacio
  if (milonga.eps_ncv.n_tokens != 0) {
    petsc_call(EPSSetDimensions(milonga.eps, 1, (PetscInt)(wasora_evaluate_expression(&milonga.eps_ncv)), PETSC_DEFAULT));
  } else {
    petsc_call(EPSSetDimensions(milonga.eps, 1, PETSC_DEFAULT, PETSC_DEFAULT));
  }
  
  // sobreescribimos con la linea de comandos
  petsc_call(EPSSetFromOptions(milonga.eps));

  // do the work!
  petsc_call(EPSSolve(milonga.eps));

  // chequeamos que haya convergido
  petsc_call(EPSGetConverged(milonga.eps, &nconv));
  if (nconv < 1) {
    wasora_push_error_message("eigen-solver did not converge");
    return WASORA_RUNTIME_ERROR;
  }
  
  // leemos la solucion
  petsc_call(EPSGetEigenpair(milonga.eps, 0, &lambda, &xi, phi, PETSC_NULL));
  
  // chequeamos que el autovalor sea real
  if (xi != 0) {
    wasora_push_error_message("eigen-solver found a complex eigenvalue (%g + i %g)", lambda, xi);
    return WASORA_RUNTIME_ERROR;
  }

  // obtenemos informacion auxiliar
  petsc_call(EPSGetErrorEstimate(milonga.eps, 0, &xi));
  wasora_var(milonga.vars.error_estimate) = (double)xi;

  petsc_call(EPSComputeError(milonga.eps, 0, EPS_ERROR_ABSOLUTE, &xi));
  wasora_var(milonga.vars.residual_norm) = (double)xi;
  
  petsc_call(EPSComputeError(milonga.eps, 0, EPS_ERROR_RELATIVE, &xi));
   wasora_var(milonga.vars.rel_error) = (double)xi;

  // depende de si resolvemos el problema que da k o 1/k
  if (milonga.spectrum == spectrum_smallest_eigenvalue) {
    *keff = 1/(lambda);
  } else {
    *keff = lambda;
  }
  
  return WASORA_RUNTIME_OK;

}

