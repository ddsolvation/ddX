function(ddx_add_executable_test basename)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${basename}.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${basename}.ref DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_test(NAME ${basename} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.py ${basename})
endfunction(ddx_add_executable_test)

add_executable("ddx_driver_testing" "ddx_driver_testing.f90")
target_include_directories("ddx_driver_testing" PUBLIC ${PROJECT_BINARY_DIR}/src)
target_link_libraries("ddx_driver_testing" ddx)
add_custom_command(TARGET "ddx_driver_testing" POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/run_test.py
            $<TARGET_FILE_DIR:ddx_driver_testing>)
file(RELATIVE_PATH path ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

ddx_add_executable_test("cosmo")
ddx_add_executable_test("cosmo_fmm")
ddx_add_executable_test("cosmo_incore")
ddx_add_executable_test("pcm")
ddx_add_executable_test("pcm_fmm")
ddx_add_executable_test("pcm_incore")
ddx_add_executable_test("lpb")
ddx_add_executable_test("lpb_fmm")
ddx_add_executable_test("lpb_incore")
