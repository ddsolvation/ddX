# @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file src/CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2020-02-17

include(define_sources.cmake)

set(COVERAGE_EXCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/llgnew.f" PARENT_SCOPE
    "${CMAKE_CURRENT_SOURCE_DIR}/cbessel.f90" PARENT_SCOPE)

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
