# Define common source files
set(SRC
    ddx_definitions.f90
    ddx_errors.f90
    ddx_parameters.f90
    ddx_constants.f90
    ddx_workspace.f90
    ddx_harmonics.f90
    ddx_core.f90
    ddx_lpb_core.f90
    ddx_operators.f90
    ddx_gradients.f90
    ddx_solvers.f90
    ddx_cosmo.f90
    ddx_pcm.f90
    ddx_lpb.f90
    ddx.f90
    ddx_legacy.f90
    ddx_multipolar_solutes.f90
    ddx_profiling.f90
    llgnew.f
    cbessel.f90
    ddx_cinterface.f90
)

# llgnew is not compliant with the f2008 standard
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(llgnew.f PROPERTIES COMPILE_FLAGS "-std=gnu")
endif()
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_source_files_properties(llgnew.f PROPERTIES COMPILE_FLAGS "-stand none")
endif()

# and this warning it is not real problem for cbesself.90
set_source_files_properties(cbessel.f90 PROPERTIES COMPILE_FLAGS "-Wno-do-subscript")

# Define source files to be build into the python version of ddx
set(SRC_PYDDX ${SRC} pyddx.cpp pyddx_classes.cpp pyddx_data.cpp)

# Define source files for the ddx.so shared library
set(SRC_DDX ${SRC})

#
# Python build
#
option(PYTHON "Build python module (Generally you want to use ./setup.py build_ext for this)" OFF)
if (PYTHON)
    find_package(pybind11 CONFIG PATHS "${PYBIND11_DIR}")
    if (NOT pybind11_FOUND)
        message(FATAL_ERROR "Pybind11 not found. Please install it.")
    endif()
    pybind11_add_module(pyddx ${SRC_PYDDX})
    target_compile_definitions(pyddx PRIVATE VERSION_INFO=${PROJECT_VERSION})
    target_link_libraries(pyddx PUBLIC
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_Fortran
    )
endif()

#
# Library build
#
option(DDX_LIBRARY "Build ddx library" ON)
if (DDX_LIBRARY)
    add_library(ddx SHARED ${SRC_DDX})
    target_link_libraries(ddx PUBLIC
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_Fortran)
    if(COVERAGE)
        target_link_libraries(ddx PRIVATE gcov)
    endif()

    add_executable(ddx_driver ddx_driver.f90)
    target_link_libraries(ddx_driver PUBLIC ddx)

    #
    # Installation rules
    #
    include(GNUInstallDirs)

    # Install library, binaries and headers
    install(TARGETS ddx EXPORT ddXTargets)
    install(TARGETS ddx_driver)
    install(FILES ddx.h TYPE INCLUDE)
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/" TYPE INCLUDE FILES_MATCHING PATTERN "*.mod")
    target_include_directories(ddx INTERFACE $<INSTALL_INTERFACE:include>)

    set(install_configdir "${CMAKE_INSTALL_LIBDIR}/ddX")
    install(EXPORT ddXTargets DESTINATION "${install_configdir}")
    export(EXPORT ddXTargets FILE ddXTargets.cmake)
    install(FILES "${PROJECT_SOURCE_DIR}/cmake/ddXConfig.cmake"
        DESTINATION "${install_configdir}"
    )
endif()
