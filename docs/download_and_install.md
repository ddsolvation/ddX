# Download and Installation  {#label_download_and_install}

## Python package

ddx is available in the form of the `pyddx` python package. This can be installed
directly from pypi:
```
pip install pyddx
```
A conda-forge package will be added soon.

If you want to build the python package from source, simply run:
```
setup.py test
```

<br />
## Source code

**Download the ddX source code** at: 
```
git@github.com:ddsolvation/ddX.git
```
**Download and install** ddX as follows:
```
> git clone git@github.com:ddsolvation/ddX.git
> cd ddX
> mkdir build
> cd build
> cmake ..
> make
```
Per default, the library is located in /src.

**Build the documentation** as follows (after you have done the above process):
```
> cd build
> make docs
```
**To see the documentation**
```
> cd ../doxygen
> pwd
```
Copy the link shown by pwd and add /index.html in a web browser
#### Hints and hacks
1. For specifying compilers use
```
cmake -D CMAKE_CXX_COMPILER=/usr/local/bin/g++-11 CMAKE_Fortran_COMPILER=/usr/local/bin/gfortran-11 ..
```
or
```
cmake -D CMAKE_CXX_COMPILER=icx CMAKE_Fortran_COMPILER=ifort ..
```
**NOTE**: Replace with the compilers you desire.
