# @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2020-12-17


###############################################################################
#                       THIS IS A TOP-LEVEL CMAKELISTS.txt                    #
#                                                                             #
#        It is intended to find all dependencies (required or optional)       #
#                    and set up corresponding variables                       #
###############################################################################


###############################################################################
##                            PRELIMINARIES                                  ##
###############################################################################

# Need to identify lowest possible CMake version
cmake_minimum_required(VERSION 3.2.3)

# Restrict building in top-level directory
if("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed.\nPlease create a "
        "build directory first and execute cmake configuration from this "
        "directory. Example: mkdir build && cd build && cmake ..")
endif()

# Read version
file(STRINGS "VERSION.txt" VERSION LIMIT_COUNT 1)

# Notify user about project name and version
message(STATUS "Configuring ddX ${VERSION}")

# Create project and check Fortran compiler
# CXX compiler is only needed to properly find Intel MKL
project(ddX VERSION ${VERSION} LANGUAGES Fortran CXX)

# If ddX is a subproject, append ddx_ prefix to all targets.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(ddx_is_subproject false)
    set(ddx_ "")
else()
    set(ddx_is_subproject true)
    set(ddx_ "ddx_")
endif()

# Add flag in case of Intel compiler to avoid fast inaccurate math
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_compile_options(-fp-model=precise)
    message(STATUS "Extending Intel compiler flags by \"-fp-model=precise\"")
endif()

# Extend path for additional cmake modules
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

#add_compile_options(-qopenmp -i8 -DINT64)
#add_compile_options(-g -debug extended -traceback -check all -fp-stack-check)
#add_compile_options(-i8  -I"${MKLROOT}/include")
#add_link_options(-L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)
#add_link_options(-mkl=parallel)

#add_compile_options(-O3)
#add_compile_options(-fopenmp)
#add_compile_options(-fdefault-integer-8)
#add_compile_options(-fcheck=all -g -fbacktrace -Og -fbounds-check)
#add_compile_options(-Waliasing -Wampersand -Wconversion -Wsurprising -Wc-binding-type -Wintrinsics-std)
#add_compile_options(-Wtabs -Wintrinsic-shadow -Wline-truncation -Wtarget-lifetime -Winteger-division)
#add_compile_options(-Wreal-q-constant -Wundefined-do-loop -Wmaybe-uninitialized)


###############################################################################
##                       REQUIREMENTS: BLAS AND LAPACK                       ##
###############################################################################

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)

###############################################################################
##                          BUILDS FOR PYTHON                                ##
###############################################################################
if (PYTHON)
    message(STATUS "Switching to python-only build.")
    find_package(pybind11 CONFIG PATHS "${PYBIND11_DIR}")
    if (NOT pybind11_FOUND)
        message(FATAL_ERROR "Pybind11 not found. Please install it.")
    endif()
    include(src/define_sources.cmake)
    list(TRANSFORM SRC_PYDDX PREPEND src/)
    pybind11_add_module(pyddx ${SRC_PYDDX})
    target_compile_definitions(pyddx PRIVATE VERSION_INFO=${VERSION})
    target_link_libraries(pyddx PUBLIC
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_Fortran
    )
    return()  # Ignore the rest of the file
endif()

###############################################################################
##                           DEFINE OPTIONS                                  ##
###############################################################################

# Options for build type for a convenient ccmake use
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    None Debug Release RelWithDebInfo MinSizeRel)

# Build tests
option(TESTS "Generate testing binaries" ON)

# Build documentation
option(DOCS "Build documentation" ON)

# Check if a code coverage report is needed
option(COVERAGE "Generate code coverage report" OFF)

###############################################################################
##                               CHECK OPTIONS                               ##
###############################################################################

# Check if need to build documentation
if(DOCS)
    find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
        set(DOCS OFF)
    endif()
endif()

#add_compile_options(-O3 -xHost -g -check all -traceback)
#add_compile_options(-O0 -g -fcheck=all -fbacktrace)

# Update compiler flags in case of code coverage
if(COVERAGE)
    # Tell user what we are doing here
    message(STATUS "Code coverage report was requested, so option TESTS is ON."
        " The report itself can be generated by \"make coverage\" command.")
    # Enable tests even if they were manually disabled
    set(TESTS ON)
    # Use CodeCoverage.cmake from cmake_modules
    include(CodeCoverage)
    # Append coverage flags
    append_coverage_compiler_flags()
    # Set excluded coverage paths
    set(COVERAGE_EXCLUDES "tests/*" "src/llgnew.f" "src/bessel.f90" "src/CBESSEL.F90")
    # Setup a target for an overall coverage
    setup_target_for_coverage_lcov(NAME ${ddx_}coverage
        EXECUTABLE ctest -R ${ddx_}tests
        LCOV_ARGS --no-external)
endif()

# Check if testing is required
if(TESTS)
    include_directories(${PROJECT_BINARY_DIR}/src)
    enable_testing()
    if(APPLE)
        if(DEFINED ADD_RPATH)
            message(STATUS "RPATH of each test will be extended by "
                "${ADD_RPATH}")
        else()
            message(STATUS "`make test` command might not work properly on "
                "Apple system due to System Integrity Protection (SIP). A "
                "walkaround is based on providing paths with dynamic libraries"
                " through -DADD_RPATH additional argument to the cmake "
                "configuration command.")
        endif()
    endif()
endif()

###############################################################################
##          AUXILIARY FUNCTION TO COPE WITH UNIQUE LOGICAL NAMES             ##
###############################################################################

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

###############################################################################
##                 BUILD LIBRARY, EXAMPLES, TESTS AND DOCS                   ##
###############################################################################

add_subdirectory("src")

# Build examples
if(EXAMPLES)
    add_subdirectory("examples")
endif()

# Build tests
if(TESTS)
    add_subdirectory("tests")
endif()

# Build documentation
if(DOCS)
    add_subdirectory("docs")
endif()


#-------------------------------------------------------------------------------
# Install rules.
# GNU Filesystem Conventions
include( GNUInstallDirs )
set(install_configdir "${CMAKE_INSTALL_LIBDIR}/ddX")

# Install library and add to <package>Targets.cmake
install(
    TARGETS ddx
    EXPORT ddXTargets
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

# Install ddx_driver
install(
    TARGETS ddx_driver
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

# Install header files
install(
    # / copies contents, not directory itself
    DIRECTORY "${PROJECT_BINARY_DIR}/src/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ddX"
    #TODO change behaviour for different os's.
    FILES_MATCHING REGEX "\\.(mod|h)$"
)

#add include dirs to target ddx
target_include_directories(ddx INTERFACE
    $<INSTALL_INTERFACE:include/ddX>  # <prefix>/include/mylib
)


# Install <package>Targets.cmake
install(
    EXPORT ddXTargets
    DESTINATION "${install_configdir}"
)

# Also export <package>Targets.cmake in build directory
export(
    EXPORT ddXTargets
    FILE "ddXTargets.cmake"
)

# Install <package>Config.cmake and <package>ConfigVersion.cmake,
# to enable find_package( <package> ).
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "ddXConfig.cmake.in"
    "ddXConfig.cmake"
    INSTALL_DESTINATION "${install_configdir}"
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/ddXConfig.cmake"
          #"${CMAKE_CURRENT_BINARY_DIR}/ddXConfigVersion.cmake"
    DESTINATION "${install_configdir}"
)
