function(ddx_add_executable_test basename)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${basename}.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${basename}.ref DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_test(NAME ${basename}
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/run_test.py
                     ${basename} --exe $<TARGET_FILE:ddx_driver>)
endfunction(ddx_add_executable_test)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/run_test.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

ddx_add_executable_test("cosmo")
ddx_add_executable_test("cosmo_fmm")
ddx_add_executable_test("cosmo_incore")
ddx_add_executable_test("pcm")
ddx_add_executable_test("pcm_fmm")
ddx_add_executable_test("pcm_incore")
ddx_add_executable_test("lpb")
ddx_add_executable_test("lpb_fmm")
ddx_add_executable_test("lpb_incore")
