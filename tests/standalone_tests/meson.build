# Copy run_test.py to the build directory
# (If run_test.py already has +x permission, this copy should preserve it on most OS.)
configure_file(
   input: 'run_test.py',
   output: 'run_test.py',
   copy: true
)

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

   # Generate a small shell wrapper from meson_test_wrapper.in
   # We replace the placeholders with the actual test name and build directory.
   meson_test_wrapper = configure_file(
      input: 'meson_test_wrapper.in',
      output: t + '_wrapper.sh',
      configuration: {
         'BUILD_DIR': meson.project_build_root(),
         'TESTNAME': t
      },
      copy: false
   )

   # Finally, register the wrapper script as a test with Meson.
   # This call has exactly two arguments: the test name, and the wrapper "File".
   test(t, meson_test_wrapper)
endforeach
