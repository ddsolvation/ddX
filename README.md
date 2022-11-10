# ddX
## General information
ddX is an open-source software package for continuum solvation models based on the domain decomposition paradigm. It contains a common interface for the three different methods ddCOSMO, ddPCM and ddLPB for the numerical solution to the COSMO, PCM and LPB solvation models. 
ddX computes for these methods the electrostatic contribution to the solvation energy, the corresponding contribution to the forces and/or the contribution to the Fock- or DFT Kohn-Sham-matrix.
Using the common interface, all three methods are accessible from a host-code that is used to model the solute and which can be on the level of QM, QM/MM, polarizable MM or MM, or as a standalone application.

## Documentation
The code is documented through the Doxygen, visit the [ddX-documentation](https://ddsolvation.github.io/ddX/).

## License 
ddX: an open-source software package for continuum solvation models based on the domain decomposition paradigm.

Copyright (c) 2022 The ddX Developers.

ddX is free software; you can redistribute it and/or modify it under the terms of the LGPL-3.0 License.
