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
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

# Retrieve Fortran compiler object
fortran_compiler = fc

# Find dependencies: BLAS, LAPACK, OpenMP
lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
   if fc_id == 'intel'
      lapack_vendor = 'mkl'
   endif
endif

lib_deps = []
if lapack_vendor == 'mkl'
   mkl_dep = []
   if fc_id == 'intel'
      mkl_dep += fc.find_library('mkl_intel_lp64')
      if get_option('openmp')
         mkl_dep += fc.find_library('mkl_intel_thread')
      endif
   elif fc_id == 'gcc'
      mkl_dep += fc.find_library('mkl_gf_lp64')
      if get_option('openmp')
         mkl_dep += fc.find_library('mkl_gnu_thread')
      endif
   else
      error('MKL not supported for this compiler')
   endif
   if not get_option('openmp')
      mkl_dep += fc.find_library('mkl_tbb_thread')
   endif
   # Ensure we add the MKL core library and assign it as blas_dep since it provides BLAS functionality.
   mkl_core = fc.find_library('mkl_core')
   mkl_dep += mkl_core
   lib_deps += mkl_dep
   add_project_arguments(['-DWITH_MKL'], language: 'fortran')
   blas_dep = mkl_core
   # Define lapack_dep even though MKL supplies LAPACK functionality:
   lapack_dep = dependency('lapack', required: false)
   if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack')
   endif
elif lapack_vendor == 'mkl-rt'
   mkl_dep = fc.find_library('mkl_rt')
   lib_deps += mkl_dep
   add_project_arguments(['-DWITH_MKL'], language: 'fortran')
   # Assume mkl_rt provides both LAPACK and BLAS functionality:
   blas_dep = mkl_dep
   lapack_dep = dependency('lapack', required: false)
   if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack')
   endif
elif lapack_vendor == 'openblas'
   blas_dep = dependency('openblas', required: false)
   if not blas_dep.found()
      blas_dep = fc.find_library('openblas')
   endif
   lib_deps += blas_dep
   if not fc.links('external dsytrs; call dsytrs(); end', dependencies: blas_dep)
      lapack_dep = dependency('lapack', required: false)
      if not lapack_dep.found()
         lapack_dep = fc.find_library('lapack')
      endif
      lib_deps += lapack_dep
   else
      # Even if the test passes the link check, define lapack_dep for consistency.
      lapack_dep = dependency('lapack', required: false)
      if not lapack_dep.found()
         lapack_dep = fc.find_library('lapack')
      endif
   endif
elif lapack_vendor == 'custom'
   foreach lib: get_option('custom_libraries')
      lib_deps += fc.find_library(lib)
   endforeach
   # Provide a blas_dep even for custom libraries
   blas_dep = dependency('blas', required: false)
   if not blas_dep.found()
      blas_dep = fc.find_library('blas')
   endif
   # Also define lapack_dep for tests:
   lapack_dep = dependency('lapack', required: false)
   if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack')
   endif
else
   lapack_dep = dependency('lapack', required: false)
   if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack')
   endif
   lib_deps += lapack_dep
   blas_dep = dependency('blas', required: false)
   if not blas_dep.found()
      blas_dep = fc.find_library('blas')
   endif
   lib_deps += blas_dep
endif

openmp_dep = dependency('openmp', required: get_option('openmp'))
if openmp_dep.found()
   lib_deps += openmp_dep
endif

subdir('config')

# Create the ddX library
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
   dependencies: lib_deps
)

# Export dependency for other projects and test suite
ddx_dep = declare_dependency(link_with: ddx_library)

# Add executable target for ddX (if needed)
executable_target = executable(
   'ddX_exec',
   sources: ['src/ddx_driver.f90'],
   link_with: [ddx_library],
   dependencies: lib_deps
)

# Add the testsuite
subdir('tests')

# Add the standalone tests
subdir('tests/standalone_tests')
