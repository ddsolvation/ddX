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
if (DDX_USE_PYTHON)
    find_package(pybind11 CONFIG PATHS "${PYBIND11_DIR}")
    if (NOT pybind11_FOUND)
        message(FATAL_ERROR "Pybind11 not found. Please install it.")
    endif()
    pybind11_add_module(pyddx ${SRC_PYDDX})
    target_compile_definitions(pyddx PRIVATE VERSION_INFO=${PROJECT_VERSION} BUILD_PYBIND11)
    target_link_libraries(pyddx PUBLIC
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_Fortran
    )
    if(PYTHON_EXT_SUFFIX)
        set_target_properties(pyddx PROPERTIES SUFFIX "${PYTHON_EXT_SUFFIX}")
    endif()
endif()

#
# Library build
#
if (DDX_BUILD_LIBRARY)
    add_library(ddx ${SRC_DDX})
    add_library(ddx::ddx ALIAS ddx)
    set_target_properties(ddx PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}
	POSITION_INDEPENDENT_CODE ON
	Fortran_MODULE_DIRECTORY "modules")
    target_include_directories(ddx PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ddX>)
    target_link_libraries(ddx PUBLIC
        ${LAPACK_LIBRARIES}
        OpenMP::OpenMP_Fortran)
    if(DDX_COVERAGE)
        target_link_libraries(ddx PRIVATE gcov)
    endif()

    add_executable(ddx_driver ddx_driver.f90)
    set_target_properties(ddx_driver PROPERTIES
        Fortran_MODULE_DIRECTORY "modules")
    target_link_libraries(ddx_driver PUBLIC ddx)

    # Installation rules
    include(GNUInstallDirs)

    # Install library, binaries and headers
    install(TARGETS ddx EXPORT ddXTargets)
    install(TARGETS ddx_driver)
    install(FILES ddx.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ddX")

    set(DDX_INSTALL_MODULEDIR "${CMAKE_INSTALL_INCLUDEDIR}/ddX" CACHE PATH
        "Installation directory for Fortran module files")
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/"
      DESTINATION "${DDX_INSTALL_MODULEDIR}"
      FILES_MATCHING PATTERN "*.mod")

    include(CMakePackageConfigHelpers)
    configure_package_config_file(
      "${PROJECT_SOURCE_DIR}/cmake/ddXConfig.cmake.in"
      "${CMAKE_CURRENT_BINARY_DIR}/ddXConfig.cmake"
      INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    )
    write_basic_package_version_file(
      "${CMAKE_CURRENT_BINARY_DIR}/ddXConfigVersion.cmake"
      VERSION "${PROJECT_VERSION}"
      COMPATIBILITY SameMinorVersion
    )
    install(EXPORT ddXTargets
      NAMESPACE ddx::
      DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
    install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/ddXConfig.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/ddXConfigVersion.cmake"
      DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    )

    # Pkg-config
    if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
      set(DDX_PC_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
    else()
      set(DDX_PC_LIBDIR "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
    endif()
    if(IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
      set(DDX_PC_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}")
    else()
      set(DDX_PC_INCLUDEDIR "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
    endif()
    if(IS_ABSOLUTE "${DDX_INSTALL_MODULEDIR}")
      set(DDX_PC_FMODDIR "${DDX_INSTALL_MODULEDIR}")
    else()
      set(DDX_PC_FMODDIR "\${prefix}/${DDX_INSTALL_MODULEDIR}")
    endif()
    configure_file(
      ${PROJECT_SOURCE_DIR}/cmake/ddX.pc.in
      ${CMAKE_CURRENT_BINARY_DIR}/ddX.pc
      @ONLY)
    install(FILES
      ${CMAKE_CURRENT_BINARY_DIR}/ddX.pc
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
endif()
