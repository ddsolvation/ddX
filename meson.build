project('ddX',
  'fortran',
  'cpp',
  version: '0.8.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.57.2',
  default_options: [
    'default_library=both',
  ],
)

install = not (meson.is_subproject() and get_option('default_library') == 'static')
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()
fortran_compiler = fc

lib_deps = []
if get_option('openmp')
  omp_dep = dependency('openmp')
  lib_deps += omp_dep
endif

lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
  if fc_id == 'intel' or fc_id == 'intel-llvm'
    lapack_vendor = 'mkl'
  endif
endif

if lapack_vendor == 'mkl'
  mkl_dep = []
  if fc_id == 'intel' or fc_id == 'intel-llvm'
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
  mkl_dep += fc.find_library('mkl_core')
  lib_deps += mkl_dep
  add_project_arguments(['-DWITH_MKL'], language: 'fortran')
elif lapack_vendor == 'mkl-rt'
  mkl_dep = fc.find_library('mkl_rt')
  lib_deps += mkl_dep
  add_project_arguments(['-DWITH_MKL'], language: 'fortran')
elif lapack_vendor == 'openblas'
  openblas_dep = dependency('openblas', required: false)
  if not openblas_dep.found()
    openblas_dep = fc.find_library('openblas')
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack')
    endif
    lib_deps += lapack_dep
  endif

elif lapack_vendor == 'custom'
  foreach lib: get_option('custom_libraries')
    lib_deps += fc.find_library(lib)
  endforeach

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

subdir('config')

ddx_library = library(
  'ddx',
  sources: [
    'src/ddx.f90',
    'src/cbessel.f90',
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
  version: meson.project_version(),
  dependencies: lib_deps,
  install: install,
  link_language: 'fortran',
)

ddx_inc = [include_directories('src'), ddx_library.private_dir_include()]
ddx_dep = declare_dependency(
  link_with: ddx_library,
  include_directories: ddx_inc,
  dependencies: lib_deps,
  version: meson.project_version(),
)

executable_target = executable(
  'ddX_exec',
  sources: ['src/ddx_driver.f90'],
  dependencies: ddx_dep,
  install: install,
  link_language: 'fortran',
)

if install
  install_headers('src/ddx.h')

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('config'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    ddx_library,
    filebase: meson.project_name(),
    name: meson.project_name(),
    description: 'Domain decomposition continuum solvation library',
    subdirs: ['', module_id],
  )
endif

if get_option('tests')
  subdir('tests')
  subdir('tests/standalone_tests')
endif
