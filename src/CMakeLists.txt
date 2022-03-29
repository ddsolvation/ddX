# @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file src/CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2020-02-17

set(SRC
    "ddx_definitions.f90"
    "ddx_parameters.f90"
    "ddx_constants.f90"
    "ddx_workspace.f90"
    "ddx_harmonics.f90"
    "ddx_core.f90"
    "ddx_lpb_core.f90"
    "ddx_operators.f90"
    "ddx_solvers.f90"
    "ddx_cosmo.f90"
    "ddx_pcm.f90"
    "ddx_lpb.f90"
    "ddx.f90"
    "ddx_mkrhs.f90"
    "llgnew.f"
    "cbessel.f90"
    "matrix_debug.f90")

set(COVERAGE_EXCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/llgnew.f" PARENT_SCOPE
    "${CMAKE_CURRENT_SOURCE_DIR}/cbessel.f90" PARENT_SCOPE)

list(TRANSFORM SRC PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)

# Define ddX library
add_library(ddx SHARED ${SRC})
# Link against BLAS and LAPACK
target_link_libraries(ddx PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
# OpenMP linkage
target_link_libraries(ddx PUBLIC OpenMP::OpenMP_Fortran)

# Define ddx driver executable
add_executable(ddx_driver "ddx_driver.f90")
target_link_libraries(ddx_driver PRIVATE ddx)

# Code coverage
if(CODE_COVERAGE)
    target_link_libraries(ddx PRIVATE gcov)
endif()

# Move SRC to parent directory
set(SRC ${SRC} PARENT_SCOPE)

