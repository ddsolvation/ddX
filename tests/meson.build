
ddx_prefix = ''
# Test sources with their arguments
test_sources = [
  ['ddx_core.f90', []],
  ['ddx_operators.f90', []],
  ['bessel.f90', []],
  ['force.f90', ['Input_force.txt']],
  ['ddx_driver.f90', [
    'data/ddpcm_force_fmm.in', 'data/ddpcm_force_fmm.out', '1E-12',
    'data/ddcosmo_force_fmm.in', 'data/ddcosmo_force_fmm.out', '1E-12'
  ]],
  ['force_ddlpb.f90', ['data/ddlpb_force.txt']],
  ['ddlpb_esolv.f90', ['data/ddlpb_force.txt']],
  ['matrix_derivatives.f90', ['data/ddlpb_force.txt']],
  ['matrix_adjoint.f90', ['data/ddlpb_force.txt']],
  ['matrix_solvers.f90', ['data/ddlpb_force.txt']],
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
    ddx_prefix + executable_name,
    src,
    dependencies: [blas_dep, lapack_dep, openmp_dep, ddx_dep],  # Link to ddx_dep
    install: false
  )

  # Add the test with ddX runtime dependency
  if args.length() == 0
    test(ddx_prefix + executable_name, test_executable)
  else
    idx = 0
    foreach arg_set : args
      test(ddx_prefix + executable_name + '_' + idx.to_string(), test_executable, args: arg_set.split())
      idx += 1
    endforeach
  endif
endforeach