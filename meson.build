# Minimum version and project definition
project('ddX',
  'fortran',
  'cpp',
  version: '0.6.0',
  license: 'LGPL-3.0-or-later',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)

install = not (meson.is_subproject() and get_option('default_library') == 'static')
has_cc = add_languages('c', required: get_option('api') or get_option('python'), native: false)

# Retrieve Fortran and C++ compiler objects
fortran_compiler = meson.get_compiler('fortran')
cxx_compiler = meson.get_compiler('cpp')

# Find dependencies: BLAS, LAPACK, OpenMP
blas_dep = dependency('blas', required: true)
lapack_dep = dependency('lapack', required: true)
openmp_dep = dependency('openmp', required: true)

subdir('config')

# Create a ddX library
ddx_library = library(
  'ddX',
  sources: [
    'src/ddx.f90', 
    'src/cbessel.f90',
    'src/ddx_driver.f90',
    'src/ddx_lpb_core.f90',
    'src/ddx_workspace.f90',
    'src/ddx_errors.f90',
    'src/ddx_lpb.f90',
    'src/llgnew.f',
    'src/ddx_cinterface.f90',
    'src/ddx_multipolar_solutes.f90',
    'src/ddx_constants.f90',
    'src/ddx_gradients.f90',
    'src/ddx_operators.f90',
    'src/ddx_core.f90',
    'src/ddx_parameters.f90',
    'src/ddx_cosmo.f90',
    'src/ddx_harmonics.f90',
    'src/ddx_pcm.f90',
    'src/ddx_definitions.f90',
    'src/ddx_legacy.f90',
    'src/ddx_solvers.f90'
  ],
  dependencies: [blas_dep, lapack_dep, openmp_dep],
)

# Export dependency for other projects and test suite
ddx_dep = declare_dependency(link_with: ddx_library)

# Add executable target for ddX (if needed)
executable_target = executable(
  'ddX_exec',
  sources: ['src/ddx_driver.f90'],
  link_with: [ddx_library, ],
  dependencies: [blas_dep, lapack_dep, openmp_dep]
)

# add the testsuite
subdir('tests')