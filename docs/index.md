# General information
ddX is an open-source software package for continuum solvation models based on
the domain decomposition paradigm. It contains a common interface for the three
different methods ddCOSMO, ddPCM and ddLPB for the numerical solution to the
COSMO, PCM and LPB solvation models. 
ddX computes for these methods the electrostatic contribution to the solvation
energy, the corresponding contribution to the forces and/or the contribution to
the Fock- or DFT Kohn-Sham-matrix.
Using a common [Python interface](@ref label_python_interface), all three
methods are accessible from a host-code that is used to model the solute. The
host code can be on the level of QM, QM/MM, polarizable MM or MM.
In this fashion an [interface to Psi4](https://psicode.org/psi4manual/master/ddx.html)
has been recently realised.
Additionally a standalone application, the [*ddx_driver* is available](@ref label_run_as_standalone).

For referencing the library and methods, see [the list of relevant references](@ref label_references).



