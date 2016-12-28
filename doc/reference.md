% Milonga reference sheet
% Jeremy Theler

This reference sheet is for [Milonga](index.html) v0.5.5-g7b0313a
. 

~~~
$ milonga
milonga v0.5.3-g67880a5+Δ 
free nuclear reactor core analysis code
$
~~~

Note that Milonga works on top of [wasora](/wasora), so you should also check the [wasora reference sheet](/wasora/reference.html) also---not to mention the [wasora RealBook](/wasora/realbook).

# Keywords

##  `FLUX_POST`

Writes a post-processing file with total and partial fluxes and optionally XS distributions.

~~~wasora
FLUX_POST { FILE <name> | FILE_PATH <file_path> } [ XS ] [ NO_MESH ] [ FORMAT { gmsh | vtk } ]
~~~



##  `IMPLICIT_BC`


~~~wasora
IMPLICIT_BC { NONE | ALLOWED }
~~~



##  `MILONGA_DEBUG`

Generates debugging and benchmarking output and/or dumps the matrices into files or the screen.

~~~wasora
MILONGA_DEBUG [ FILE <file_id> | [ FILE_PATH <file_path> ] [ MATRICES_ASCII ] [ MATRICES_ASCII_STRUCTURE ] [ MATRICES_PETSC_BINARY ] [ MATRICES_PETSC_COMPRESSED_BINARY ] [ MATRICES_PETSC_ASCII ] [ MATRICES_PETSC_OCTAVE ] [ MATRICES_PETSC_DENSE ] [ MATRICES_X ] [ MATRICES_SNG ] [ MATRICES_SNG_STRUCT ] [ MATRICES_SIZE <expr> ] [ MATRICES_STRIDE <expr> ] [ INCLUDE_INPUT ]
~~~



##  `MILONGA_PROBLEM`

Defines the number of spatial dimensions and groups of neutron energies.      
It also selects the formulation of the neutronic problem to be solved (i.e. diffusion or tranport)
and the spatial discretization scheme (i.e. finite volumes or finite elements).
If several meshes are defined, it selects over which one it is that the neutronic problem is solved.      

~~~wasora
MILONGA_PROBLEM [ DIMENSIONS <expr> ] [ GROUPS <expr> ] [ MESH <identifier> ] [ SCHEME { volumes | elements } ] [ FORMULATION { diffusion | s2 | s4 | s6 | s8 } ] [ VOLHOM ]
~~~



##  `MILONGA_SOLVER`

Sets options related to the eigen-solver.

~~~wasora
MILONGA_SOLVER [ ROUTINE <loadable_routine> ] [ SPECTRUM { largest_eigenvalue | smallest_eigenvalue } ] [ EPS_TYPE { krylovschur | gd | jd | power | arnoldi | subspace | ... } ] [ ST_TYPE { sinvert | shift | cayley | precond } ] [ KSP_TYPE { gmres | bcgs | bicg | richardson | chebyshev | ... } ] [ PC_TYPE { lu | none | sor | bjacobi | cholesky | ... } ] [ SUBSPACE_DIM <expr> ] [ ST_SHIFT <expr> ] [ ST_ANTI_SHIFT <expr> ]
~~~



List of `EPS_TYPE`s http://www.grycap.upv.es/slepc/documentation/current/docs/manualpages/EPS/EPSType.html
         
List of `ST_TYPE`s http://www.grycap.upv.es/slepc/documentation/current/docs/manualpages/ST/STType.html
         
List of `KSP_TYPE`s http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html
         
List of `PC_TYPE`s http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html
         

##  `MILONGA_STEP`

Solves the linear eigenvalue problem.

~~~wasora
MILONGA_STEP [ JUST_BUILD | JUST_SOLVE ]
~~~






# Variables

##  `eigen_error_estimate`

Absolute error estimate $|\lambda - \lambda_\text{real}|$ (`|lambda - lambda_real|`)
of the obtained eigenvalue.



##  `eigen_rel_error`

Relative error $|R \phi - \lambda F \phi|_2/|\lambda \phi|_2$ `|Rx-lFx|_2/|kx|_2`
reached by the eigensolver.



##  `eigen_rel_tolerance`

Tolerance passed the iterative eigensolver, relative to the matrices norm.
By default it is zero, meaning the default value of the SLEPc library.



##  `eigen_residual_norm`

Residual norm $|R \phi - \lambda F \phi|_2$ (`|Rx-lFx|_2`)
reached by the eigensolver.



##  `keff`

Efective multiplication factor as computed by solving the
multigroup neutron diffusion equation. It is equal to 1.0 until
`MILONGA_STEP` is executed.  



##  `memory_usage_global`

Maximum resident set size (global memory used), in bytes.



##  `memory_usage_petsc`

Maximum resident set size (memory used by PETSc), in bytes.



##  `memory_use`

Total available memory, in bytes.



##  `petsc_flops`

Number of floating point operations performed by PETSc/SLEPc.



##  `power`

Power setpoint used to normalize the computed fluxes.
By default it is zero, meaning the fluxes are normalized so
the mean value is equal to one. If this variable is set,
then the property `eSigmaF` should be nonzero through the domain.



##  `sn_a_weight`

Factor used to control the upwinding of the discrete ordinates~$S_N$ formulation.
For the finite volumes scheme, the difference factor $a = 1/2 (1+\alpha)$ is used to weight the adjacent cell fluxes,
i.e. $\alpha = 0$ ($a=1/2$) corresponds to diamond difference
 and $\alpha = 1$ ($a=1$) corresponds to full upwinding.
For the finite elements scheme, a Streamline-Upwind Petrov-Galerkin stabilization term
$1/2 \alpha \ell/\| \vec{\Omega} \cdot ( \vec{\Omega} \cdot \nabla h)\|$ is added to the weighting functions.
A value of~$\alpha=0$ corresponds to the unstable original Galerkin scheme and
           $\alpha=1$ to the cotangent-optimum upwind term.
For~$\alpha < 1$ expect oscillations in the fluxes and even non-convergence for small values of~$\alpha$.
Default is $\alpha=0.5$.  



##  `time_cpu_build`

CPU time insumed to build the problem matrices, in seconds.



##  `time_cpu_ini`

CPU time insumed to initialize the problem, in seconds.



##  `time_cpu_solve`

CPU time insumed to solve the eigen-problem, in seconds.



##  `time_petsc_build`

CPU time insumed by PETSc to build the problem matrices, in seconds.



##  `time_petsc_ini`

CPU time insumed by PETSc to initialize the problem, in seconds.



##  `time_petsc_solve`

CPU time insumed by PETSc to solve the eigen-problem, in seconds.



##  `time_wall_build`

Wall time insumed to build the problem matrices, in seconds.



##  `time_wall_ini`

Wall time insumed to initialize the problem, in seconds.



##  `time_wall_solve`

Wall time insumed to solve the eigen-problem, in seconds.



##  `time_wall_total`

Wall time insumed to initialize, build and solve, in seconds.
CPU time insumed to initialize, build and solve, in seconds.
CPU time insumed by PETSc to initialize, build and solve, in seconds.



##  `unknowns`

Number of total unknowns (size) of the problem. It is equal
to the number of spatial unknowns times the number of energy groups.  






