include(define_sources.cmake)

# Define ddX library
add_library(ddx SHARED ${SRC_DDX})
target_link_libraries(ddx PUBLIC
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    OpenMP::OpenMP_Fortran)
if(COVERAGE)
    target_link_libraries(ddx PRIVATE gcov)
endif()

# Define ddx driver executable
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
