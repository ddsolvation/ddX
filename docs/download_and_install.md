# Download and Installation  {#label_download_and_install}
**Download precompiled binaries** at:
1. Conda
2. Pypi

**Download the ddX source code** at: 
``` markdown
> git clone git@github.com:ACoM-Computational-Mathematics/ddX.git
```
**Download and install** ddX as follows:
``` markdown
> git clone git@github.com:ACoM-Computational-Mathematics/ddX.git
> cd ddX
> mkdir build
> cd build
> cmake .. 
> make
```
Per default, the library is located in /src.

#### Hints and hacks
1. Specify XXX by 
``` markdown 
cmake -D CMAKE_CXX_COMPILER=/usr/local/bin/g++-11 CMAKE_C_COMPILER=/usr/local/bin/gcc-11 ..
```
2. disactivate OpenMP-support: 
``` markdown 
cmake -D OPENMP=OFF ..
```
