# Run as standalone  {#label_run_as_standalone}
The standalone application requires an input file within the specific format outlined below. 
Run the code as
``` markdown
> ./src/ddx_driver *path-to-input-file*
```
**Example:** (using the default input-file)
``` markdown
> ./src/ddx_driver ../tests/Input.txt
```

<br />
**Description of input parameters in input-file**

Line number  | Type           | Description   
------------- | ------------- | -------------- 
1  | Integer | Printing flag, range=[0,1,2,3,4]  <br /> 0: Minimal verbosity (energy and forces), 1: + number of iterations of solvers, timers, 2: + solution vector, 3: + right hand side vector, 4: + all debug outputs   
2  | Integer  | Number of OpenMP cores to be used
3  | Integer  | Specification of the model <br />  1: COSMO, 2: PCM, 3: LPB
4  | Integer  | Maximal degree \f$\ell_{\max}\f$ of modeling spherical harmonics
5  | Integer  | Approximate number of Lebedev grid points
6  | Float    | Dielectric permittivity constant \f$\varepsilon_s\f$ of the bulk solvent
7  | Float    | Shift \f$s\f$ of the regularized characteristic function \f$\chi_\eta\f$, range=[-1,1] 
8  | Float    | Regularization parameter \f$\eta\f$ of the smoothing function \f$\chi_\eta\f$, range=[0,1]
9  | Float    | Debye Hückel parameter \f$\kappa_s\f$ of the bulk solvent
10 | Integer  | Iterative solver for linear systems <br /> 1: Jacobi/DIIS, 2: GMRESR
11 | Float    | The relative threshold \f$tol\f$ for the iterative solver
12 | Integer  | Maximum number of iterations of the iterative solver before stopping
13 | Integer  | Number of Jacobi DIIS extrapolation points (Default)
14 | Integer  | Number of last vectors GMRESR works with
15 | Integer  | Dimension of GMRESR
16 | Bool     | Whether to compute (1) or not (0) the forces
17 | Bool    `| Whether to use (1) or not (0) the FMM (only for PCM and LPB???)
18 | Integer  | Max degree of multipole spherical harmonics \f$\tilde\ell_{\max}\f$ for the FMM (default value \f$\tilde\ell_{\max}=\ell_{\max}\f$???)
19 | Integer  | Max degree of local spherical harmonics \f$p_{\max}\f$ for the FMM (default???)
20 | Integer  | Number of spheres of the atomic structure
21--END | Float[5] | The flollowing lines specify for each atom (should be consistent with entry on line 20) the values c,x,y,z,r where the partial charge (c), x-coordinate (x), y-coordinate (y), z-coordinate (z), radius (r) are listed (see example below). (units???)

**Remarks:**
1. The order of the lines is essential


**Example of an input file:**
```
1           ! Printing flag. The larger the value[integer], the more verbose the output
1           ! Number of OpenMP cores to be used
1           ! Specification of the model: 1 for COSMO, 2 for PCM and 3 for LPB
7           ! Maximal degree of modeling spherical harmonics
302         ! Approximate number of Lebedev grid points
78.3553     ! Dielectric permittivity constant
0.0         ! Shift of the regularized characteristic function
0.1         ! Regularization parameter
0.0         ! Debye Hückel parameter
2           ! Iterative solver for linear systems: 1: Jacobi/DIIS , 2: GMRESR solver
1d-8        ! The relative threshold for the iterative solver
200         ! Maximum number of iterations
25          ! Number of Jacobi DIIS extrapolation points
1           ! Number of last vectors GMRESR works with
0           ! Dimension of GMRESR
1           ! Whether to compute (1) or not (0) forces
1           ! Whether to use (1) or not (0) the FMM
20          ! Max degree of multipole spherical harmonics for the FMM
20          ! Max degree of local spherical harmonics for the FMM
12          ! Number of spheres of the atomic structure
-0.04192   0.00000   2.29035   1.32281   4.00253    ! partial charge, x, y, z, radius[UNIT?????]
-0.04192   0.00000   2.29035  -1.32281   4.00253
-0.04198   0.00000   0.00000  -2.64562   4.00253
-0.04192   0.00000  -2.29035  -1.32281   4.00253
-0.04192   0.00000  -2.29035   1.32281   4.00253
-0.04198   0.00000   0.00000   2.64562   4.00253
 0.04193   0.00103   4.05914   2.34326   2.99956
 0.04193   0.00103   4.05914  -2.34326   2.99956
 0.04197   0.00000   0.00000  -4.68652   2.99956
 0.04193  -0.00103  -4.05914  -2.34326   2.99956
 0.04193  -0.00103  -4.05914   2.34326   2.99956
 0.04197   0.00000   0.00000   4.68652   2.99956
```


