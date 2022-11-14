# Python interface  {#label_python_interface}

## Download the python package pyddx
You can download the python-interface with ddx from pypi.org following:
``` markdown
pip install pyddx
```
Please find [here](https://pypi.org/project/pyddx/) the PyPI-webpage.


<br />
## Example
It is simplest by starting with an example. 
The following example can be found in examples/run_ddx.py. 

In the terminal, run the example (from the root directory)
``` markdown
python examples/run_ddx.py 
```

or launch it interactively within python:
``` markdown
import pyddx
import numpy as np

tobohr = 1 / 0.52917721092
charges = np.array([
    -0.04192, -0.04192, -0.04198, -0.04192, -0.04192, -0.04198,
    0.04193, 0.04193,  0.04197,  0.04193,  0.04193,  0.04197
])
rvdw = tobohr * np.array([
    4.00253, 4.00253, 4.00253, 4.00253, 4.00253, 4.00253,
    2.99956, 2.99956, 2.99956, 2.99956, 2.99956, 2.99956
])
centres = tobohr * np.array([
    [ 0.00000,  2.29035,  1.32281],  
    [ 0.00000,  2.29035, -1.32281],  
    [ 0.00000,  0.00000, -2.64562],  
    [ 0.00000, -2.29035, -1.32281],  
    [ 0.00000, -2.29035,  1.32281],  
    [ 0.00000,  0.00000,  2.64562],  
    [ 0.00103,  4.05914,  2.34326],  
    [ 0.00103,  4.05914, -2.34326],  
    [ 0.00000,  0.00000, -4.68652],  
    [-0.00103, -4.05914, -2.34326],  
    [-0.00103, -4.05914,  2.34326],  
    [ 0.00000,  0.00000,  4.68652],  
]).T

model = pyddx.Model("pcm", charges, centres, rvdw, solvent_epsilon=78.3553)
nuclear = model.solute_nuclear_contribution()
state = model.initial_guess()
state = model.solve(state, nuclear["phi"])
state = model.adjoint_solve(state, nuclear["psi"])
force = model.force_terms(state, nuclear["phi"], nuclear["gradphi"], nuclear["psi"])

energy = 0.5 * np.sum(state.x * nuclear["psi"])
print(energy)
print(force)
```

<br />
## Use ddX in python
We now proceed with further informations.
Within python, import the pyddx-package
``` markdown
>>> import pyddx
```

<br />
### Model constructor
Define the model as follows
``` markdown
model = pyddx.Model(modelstr, charges, centres, rvdw, solvent_epsilon)
```
Here, the model is constructed with the following mandatory arguments:

Argument number | name | type        | description   
------------- | ------------- | -------------- | -------------- 
1  | model | string | String with the model name, i.e. "cosmo", "pcm" or "lpb"
2  | sphere_charges | array  | Numpy-array containing the charges (of dimension 1 x n_spheres)
3  | sphere_centres | array  | Numpy-array  containing the centres of the spheres/atoms (of dimension 3 x n_spheres)
4  | sphere_radii | array  | Numpy-array  containing the radii of the cavity spheres (of dimension 1 x n_spheres)
5  | solvent_epsilon | double | Value of the solvent dielectric permittivity

The following arguments are optional:

Name | type        | description   
------------- | ------------- | -------------- 
solvent_kappa  | double | Debye Hückel parameter \f$\kappa_s\f$ of the bulk solvent
eta  | double  | Regularization parameter \f$\eta\f$ of the smoothing function \f$\chi_\eta\f$, range=[0,1]
lmax  | int  | Maximal degree \f$\ell_{\max}\f$ of modeling spherical harmonics
n_lebedev  | int  | Approximate number of Lebedev grid points
incore | bool | Whether to compute and store sparse matrices (1) or apply the matrix-vector product on the fly (0), range={0,1} <br /> The sparse matrices are the solution matrix of ddCOSMO referred to as \f$L\f$ used in ddCOSMO and ddPCM, and the matrices \f$A\f$ and \f$B\f$ in ddLPB.
maxiter  | int  | Maximum number of iterations of the iterative solver before stopping
jacobi_n_diis | int | Number of Jacobi/DIIS extrapolation points
enable_fmm | bool | Whether to use (1) or not (0) the FMM, range={0,1} 
fmm_multipole_lmax | int | Max degree of multipole spherical harmonics \f$\tilde\ell_{\max}\f$ for the FMM (recommended value \f$\tilde\ell_{\max}=\ell_{\max}\f$). 
fmm_local_lmax | int | Max degree of local spherical harmonics \f$p_{\max}\f$ for the FMM (recommended value \f$p_{\max}=6\f$). 
n_proc | int | Number of OpenMP cores to be used

<br />
### Functions
The following functions are available 
Name | description   
------------- | ------------- 
initial_guess | Builds an initial guess for the state variable
solve         | Solves the linear system that defines the state variable
adjoint_solve | Solves the adjoint linear system that is required for the force computation and/or the contribution to the Fock-/DFT-operator
force_terms   | Computes the forces due to the presence of the solvent, i.e. minus the gradient of the solvation energy w.r.t. the nuclear coordinates
solute_nuclear_contribution | Returns the terms of the nuclear contribution to the solvation as a python dictionary
scaled_ylm    | Computes the scaled harmonics with reference to a atomic sphere "sphere" of radius "r" centred at "a" compute \f$\frac{4π}{2 \ell +1} \cdot \frac{\|x-a\|^\ell}{r^\ell} \cdot Y_\ell^m(\|x-a\|)\f$

For further informations about the arguments, see the documentation in python 
``` markdown
>>> help(pyddx)
```


