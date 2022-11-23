# ddX

| **Documentation** | **Build Status** | **Installation** |
|:----------------- |:---------------- |:----------- |
| [![][docs-img]][docs-url] | [![][ci-img]][ci-url] | [![][pypi-img]][pypi-url] [![][conda-img]][conda-url] [![][license-img]][license-url]

[docs-img]: https://img.shields.io/badge/doc-latest-blue.svg
[docs-url]: https://ddsolvation.github.io/ddX
[ci-img]: https://github.com/ddsolvation/ddX/actions/workflows/CI.yaml/badge.svg
[ci-url]: https://github.com/ddsolvation/ddX/actions/workflows/CI.yaml
[cov-img]: https://coveralls.io/repos/ddsolvation/ddX/badge.svg?branch=main&service=github
[cov-url]: https://coveralls.io/github/ddsolvation/ddX?branch=main
[license-img]: https://img.shields.io/badge/License-LGPL%20v3-blue.svg
[license-url]: https://github.com/ddsolvation/ddx/blob/master/LICENSE
[pypi-img]: https://img.shields.io/pypi/v/pyddx
[pypi-url]: https://pypi.org/project/pyddx
[conda-img]: https://anaconda.org/conda-forge/pyddx/badges/version.svg
[conda-url]: https://anaconda.org/conda-forge/pyddx

ddX is an open-source software package for continuum solvation models based on
the domain decomposition paradigm. It contains a common interface for the three
different methods ddCOSMO, ddPCM and ddLPB for the numerical solution to the
COSMO, PCM and LPB solvation models.  ddX computes for these methods the
electrostatic contribution to the solvation energy, the corresponding
contribution to the forces and/or the contribution to the Fock- or DFT
Kohn-Sham-matrix.  Using the common interface, all three methods are accessible
from a host-code that is used to model the solute and which can be on the level
of QM, QM/MM, polarizable MM or MM, or as a standalone application.

## Documentation
The code is documented through the Doxygen, visit the [ddX-documentation](https://ddsolvation.github.io/ddX/).

## License 
ddX: an open-source software package for continuum solvation models based on the domain decomposition paradigm.

Copyright (c) 2022 The ddX Developers.

ddX is free software; you can redistribute it and/or modify it under the terms of the LGPL-3.0 License.
