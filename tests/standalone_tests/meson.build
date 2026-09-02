run_test_py = configure_file(
   input: 'run_test.py',
   output: 'run_test.py',
   copy: true
)

python3 = find_program('python3')

# Define a list of test basenames.
standalone_tests = [
   'cosmo',
   'cosmo_fmm',
   'cosmo_incore',
   'pcm',
   'pcm_fmm',
   'pcm_incore',
   'lpb',
   'lpb_fmm',
   'lpb_incore'
]

foreach t : standalone_tests
   # Copy each test’s .txt and .ref files into the build dir.
   configure_file(
      input: t + '.txt',
      output: t + '.txt',
      copy: true
   )
   configure_file(
      input: t + '.ref',
      output: t + '.ref',
      copy: true
   )

   test(t, python3, args: [run_test_py, t, '--exe', executable_target], timeout: 240)
endforeach
