# @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file src/CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2020-02-17

set(SRC
    "ddx_core.f90"
    "ddx_operators.f90"
    "ddx_solvers.f90"
    "ddx.f90"
    "ddx_lpb.f90"
    "bessel.f90"
    "llgnew.f")

# Preprocess OPENMP directive (disabled currently)
#set_source_files_properties("dd_core.f90" PROPERTIES COMPILE_FLAGS -cpp)

set(COVERAGE_EXCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/llgnew.f" PARENT_SCOPE
    "${CMAKE_CURRENT_SOURCE_DIR}/bessel.f90" PARENT_SCOPE)

list(TRANSFORM SRC PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)

# Define ddX library
add_library(ddx SHARED ${SRC})
# Link against BLAS and LAPACK
target_link_libraries(ddx PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
# OpenMP linkage
if(OPENMP)
    target_link_libraries(ddx PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Define ddx driver executable
add_executable(ddx_driver "ddx_driver.f90")
target_link_libraries(ddx_driver PRIVATE ddx)

# Code coverage
if(CODE_COVERAGE)
    target_link_libraries(ddx PRIVATE gcov)
endif()

# Move SRC to parent directory
set(SRC ${SRC} PARENT_SCOPE)

