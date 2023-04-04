if(APPLE)
    if(DEFINED ADD_RPATH)
        message(STATUS "RPATH of each test will be extended by ${ADD_RPATH}")
    else()
        message(STATUS "`make test` command might not work properly on "
            "Apple system due to System Integrity Protection (SIP). A "
            "walkaround is based on providing paths with dynamic libraries"
            " through -DADD_RPATH additional argument to the cmake "
            "configuration command.")
    endif()
endif()

# Add single test executable, presented by its single source, along with a list
# of combinations of parameters (stored as strings). One test executable = one
# source file. Number of tests = number of combinations. These tests are then
# available for ctest.
function(ddx_add_test src)
    # Parse path to get name of executable and its corresponding unique logical
    # name, based on its path and name
    get_filename_component(exec ${src} NAME_WE)
    file(RELATIVE_PATH path ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    string(REPLACE "/" "_" suffix ${path})
    set(logical ${suffix}_${exec})
    # Register executable with its unique logical name
    add_executable(${ddx_}${logical} ${src})
    target_include_directories(${ddx_}${logical} PUBLIC ${PROJECT_BINARY_DIR}/src)
    # Set name of executable file
    set_target_properties(${ddx_}${logical} PROPERTIES OUTPUT_NAME ${exec})
    # Add dependency to the global coverage target
    if(COVERAGE)
        add_dependencies(${ddx_}coverage ${ddx_}${logical})
    endif()
    # Try to set environment for a test on Apple MacOS
    if(APPLE AND DEFINED ADD_RPATH)
        set_target_properties(${ddx_}${logical} PROPERTIES BUILD_RPATH
            "${ADD_RPATH}")
    endif()
    # Link to ddx
    target_link_libraries(${ddx_}${logical} ddx)
    if(COVERAGE)
        # Add coverage for all tests together (they have ${logical}_${i} names)
        setup_target_for_coverage_lcov(NAME ${ddx_}coverage_${logical}
            EXECUTABLE ctest -R ${ddx_}${logical}
            DEPENDENCIES ${ddx_}${logical}
            LCOV_ARGS --no-external)
    endif()
    # If no arguments were provided, add test without parameters
    if(${ARGC} EQUAL 1)
        add_test(NAME ${ddx_}${logical} COMMAND
            $<TARGET_FILE:${ddx_}${logical}>)
        message(STATUS "Adding 1 test: ${path}/${exec}")
    # If only one test with arguments must be added
    elseif(${ARGC} EQUAL 2)
        # Convert arguments from string to list
        string(REPLACE " " ";" param ${ARGV1})
        string(REPLACE "\n" ";" param "${param}")
        add_test(NAME ${ddx_}${logical} COMMAND
            $<TARGET_FILE:${ddx_}${logical}> ${param})
        message(STATUS "Adding 1 test: ${path}/${exec}")
    # Add multiple tests and add indexing for better readability of output
    else()
        math(EXPR ntests ${ARGC}-1)
        message(STATUS "Adding ${ntests} tests: ${path}/${exec}")
        # For every string of parameters, register corresponding test
        foreach(i RANGE 1 ${ntests})
            # Convert arguments from string to list
            string(REPLACE " " ";" param ${ARGV${i}})
            string(REPLACE "\n" ";" param "${param}")
            add_test(NAME ${ddx_}${logical}_${i} COMMAND
                $<TARGET_FILE:${ddx_}${logical}> ${param})
        endforeach()
    endif()
endfunction(ddx_add_test)

#
# Actual tests
#
ddx_add_test("ddx_core.f90")
ddx_add_test("ddx_operators.f90")
ddx_add_test("bessel.f90")
ddx_add_test("force.f90" "${CMAKE_CURRENT_SOURCE_DIR}/Input_force.txt")
ddx_add_test("ddx_driver.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/ddpcm_force_fmm.in
    ${CMAKE_CURRENT_SOURCE_DIR}/data/ddpcm_force_fmm.out 1E-12"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/ddcosmo_force_fmm.in
    ${CMAKE_CURRENT_SOURCE_DIR}/data/ddcosmo_force_fmm.out 1E-12")
add_test(NAME ${ddx_}tests_driver_cosmo COMMAND ddx_driver
    "${CMAKE_CURRENT_SOURCE_DIR}/Input_cosmo.txt")

ddx_add_test("force_ddlpb.f90" "${CMAKE_CURRENT_SOURCE_DIR}/data/ddlpb_force.txt")
ddx_add_test("ddlpb_esolv.f90" "${CMAKE_CURRENT_SOURCE_DIR}/data/ddlpb_force.txt")
ddx_add_test("matrix_derivatives.f90" "${CMAKE_CURRENT_SOURCE_DIR}/data/ddlpb_force.txt")
ddx_add_test("matrix_adjoint.f90" "${CMAKE_CURRENT_SOURCE_DIR}/data/ddlpb_force.txt")
ddx_add_test("matrix_solvers.f90" "${CMAKE_CURRENT_SOURCE_DIR}/data/ddlpb_force.txt")
ddx_add_test("m2l.f90")
ddx_add_test("multipolar_solutes.f90")
