
project_dir = meson.source_root()

# Test sources with their arguments
# Test sources with their arguments
test_sources = [
  ['ddx_core.f90', []],
  ['ddx_operators.f90', []],
  ['bessel.f90', []],
  ['ddx_driver.f90', [
    join_paths(project_dir, 'tests/data/ddpcm_force_fmm.in') + ' ' + join_paths(project_dir, 'tests/data/ddpcm_force_fmm.out') + ' 1E-12',
    join_paths(project_dir, 'tests/data/ddcosmo_force_fmm.in') + ' ' + join_paths(project_dir, 'tests/data/ddcosmo_force_fmm.out') + ' 1E-12'
  ]],
  ['force_ddlpb.f90', [join_paths(project_dir, 'tests/data/ddlpb_force.txt')]],
  ['ddlpb_esolv.f90', [join_paths(project_dir, 'tests/data/ddlpb_force.txt')]],
  ['matrix_derivatives.f90', [join_paths(project_dir, 'tests/data/ddlpb_force.txt')]],
  ['matrix_adjoint.f90', [join_paths(project_dir, 'tests/data/ddlpb_force.txt')]],
  ['matrix_solvers.f90', [join_paths(project_dir, 'tests/data/ddlpb_force.txt')]],
  ['m2l.f90', []],
  ['multipolar_solutes.f90', []],
  ['error.f90', []]
]

# Create test executables and add tests
foreach source_entry : test_sources
  src = source_entry[0]
  args = source_entry[1]
  executable_name = src.split('.')[0]

  test_executable = executable(
    executable_name,
    src,
    dependencies: [blas_dep, lapack_dep, openmp_dep, ddx_dep],  # Link to ddx_dep
    install: false
  )

  # Add the test with ddX runtime dependency
  if args.length() == 0
    test(executable_name, test_executable)
  else
    idx = 0
    foreach arg_set : args
      test(
        executable_name + '_' + idx.to_string(),
        test_executable,
        args: arg_set.split(),
        timeout: 240
      )

      # print out the exact command with argument for each test
      message('test ' + idx.to_string() + ' ' + executable_name + ' ' + arg_set)

      idx += 1
    endforeach
  endif
endforeach