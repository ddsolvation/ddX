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

## Source Code

**Download the ddX source code** from GitHub:

```
git clone git@github.com:ddsolvation/ddX.git
```

Then change into the cloned directory:

```
cd ddX
```

The main Fortran sources for ddX reside in the `src/` folder. You can build ddX using either **CMake** or **Meson**, as outlined below.

---

## Building with CMake

1. Create a build directory and enter it:
   ```bash
   mkdir build
   cd build
   ```
2. Run CMake to configure:
   ```bash
   cmake ..
   ```
3. Compile the library and executables:
   ```bash
   make
   ```
4. (Optional) Run the test suite:
   ```bash
   make test
   ```

By default, the compiled library and executables will appear in the `build` folder. The original sources remain in `src/`.

**Build the documentation** as follows (after you have done the above process):
```bash
cd build
make docs
```
**To see the documentation**
```bash
cd ../doxygen
pwd
```
Copy the link shown by pwd and add /index.html in a web browser

**Specifying compilers** can be done by passing the desired compilers to CMake:
```
cmake -D CMAKE_CXX_COMPILER=/usr/local/bin/g++-11 CMAKE_Fortran_COMPILER=/usr/local/bin/gfortran-11 ..
```
or
```
cmake -D CMAKE_CXX_COMPILER=icx CMAKE_Fortran_COMPILER=ifort ..
```
**NOTE**: Replace with the compilers you desire.


---

## Building with Meson

1. Set up a build directory:
   ```bash
   meson setup build
   ```
2. Compile:
   ```bash
   meson compile -C build
   ```
3. Run the test suite (and print any error logs):
   ```bash
   meson test -C build
   ```

Again, the built library and any executables will appear in the `build` folder, while the ddX source remains in `src/`.


