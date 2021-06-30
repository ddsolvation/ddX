# ddX
## General information
ddX is an open-source software package for continuum solvation models based on the domain decomposition paradigm. It contains a common interface for the three different methods ddCOSMO, ddPCM and ddLPB for the numerical solution to the COSMO, PCM and LPB solvation models. 
Using the common interface, all three methods are accessible from a host-code that is used to model the solute and which can be on the level of QM, QM/MM, polarizable MM or MM.

For the above methods, ddX computes the electrostatic contribution to the solvation energy, the corresponding contribution to the forces and/or the contribution to the Fock- or DFT Kohn-Sham-matrix.

## Documentaton
The code is documented through the Doxygen.

## License [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

ddX: an open-source software package for continuum solvation models based on the domain decomposition paradigm.

Copyright (c) 2021 The ddX Developers.

ddX is free software; you can redistribute it and/or modify it under the terms of the MIT Licence.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## References
The relevant literature for the ddCOSMO, ddPCM and ddLPB method is listed below. If you're using one of the methods through the ddX-library, please cite the corresponding paper(s).

#### ddCOSMO: 
[1] E. Cancès, Y. Maday, B. Stamm, Domain decomposition for implicit solvation models, Journal of Chemical Physics, Vol. 139, No. 5, pp. 054111 (2013)

[2] F. Lipparini, B. Stamm, E. Cancès, Y. Maday, B. Mennucci, A Fast Domain Decomposition Algorithm for Continuum Solvation Models: Energy and First Derivatives, J. Chem. Theory Comput., Vol. 9, No. 8, pp. 3637–3648 (2013)

[3] F. Lipparini, L. Lagardère, G. Scalmani, B. Stamm, E. Cancès, Y. Maday, J.-P. Piquemal, M. J. Frisch, B. Mennucci, Quantum calculations in solution for large to very large molecules: a new linear scaling QM/continuum approach, J. Phys. Chem. Lett., Vol. 5, No. 4, pp. 953-958 (2014)

[4] F. Lipparini, G. Scalmani, L. Lagardère, B. Stamm, E. Cancès, Y. Maday, J.-P. Piquemal, M. Frisch and B. Mennucci, Quantum, Classical and Hybrid QM/MM Calculations in Solution: General Implementation of the ddCOSMO Linear Scaling Strategy, Journal of Chemical Physics., Vol. 141, pp. 184108 (2014)

[5] S. Caprasecca, S. Jurinovich, L. Lagardère, B. Stamm, F. Lipparini, Achieving linear scaling in computational cost for a fully polarizable MM/Continuum embedding, J. Chem. Theory Comput., Vol. 11, No. 2, pp. 694-704 (2015) 

[6] F. Lipparini, L. Lagardère, Ch. Raynaud, B. Stamm, E. Cancès, B. Mennucci, M. Schnieders, P. Ren,Y. Maday, J.-P. Piquemal, Polarizable Molecular Dynamics in a Polarizable Continuum Solvent, J. Chem. Theory Comput., Vol. 11, No. 2, pp. 623-634 (2015) 

[7] B. Stamm, L. Lagardère, G. Scalmani, P. Gatto, E. Cances, J. Piquemal, Y. Maday, B. Mennucci, F. Lipparini, How to make continuum solvation incredibly fast in a few simple steps: a practical guide to the domain decomposition paradigm for the Conductor-like Screening Model, Int. J. Quantum Chem. (2019)

#### ddPCM: 
[8] B. Stamm, E. Cancès, F. Lipparini,Y. Maday, A new discretization for the Polarizable Continuum Model within the domain decomposition paradigm, J. Chem. Physics, Vol. 144, 054101 (2016)

[9] P. Gatto, F. Lipparini, and B. Stamm, Computation of Forces arising from the Polarizable Continuum Model within the Domain-Decomposition Paradigm,  J. Chem. Physics, Vol. 147, 224108 (2017) 

[10] M. Nottoli, B. Stamm, G. Scalmani, F. Lipparini, Quantum calculations in solution of energies, structures and properties with a domain decomposition polarizable continuum model, J. Chem. Theory Comput., Vol. 15, No. 11, pp. 6061–6073 (2019)

#### ddLPB: 
[11] Ch. Quan, B. Stamm, Y. Maday, A Domain Decomposition Method for the Poisson-Boltzmann Solvation Models, SIAM J. Sci. Comput., Vol. 41, No. 2, pp. B320-B350, (2019)

