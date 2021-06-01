# @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file tests/CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2021-02-25

set(TESTS
    "ddx_core.f90"
    "ddx_operators.f90"
    "force.f90"
    "main_ddLPB.f90")

ddx_add_test("ddx_core.f90")
ddx_add_test("ddx_operators.f90")
ddx_add_test("main_ddLPB.f90" "${CMAKE_SOURCE_DIR}/tests/Input_ddLPB.txt")
ddx_add_test("force.f90" "${CMAKE_SOURCE_DIR}/tests/Input_force.txt")
ddx_add_test("ddx_driver.f90"
    "${CMAKE_SOURCE_DIR}/tests/data/ddpcm_force_fmm.in
    ${CMAKE_SOURCE_DIR}/tests/data/ddpcm_force_fmm.out 1E-13"
    "${CMAKE_SOURCE_DIR}/tests/data/ddcosmo_force_fmm.in
    ${CMAKE_SOURCE_DIR}/tests/data/ddcosmo_force_fmm.out 1E-13")

add_test(NAME driver_cosmo COMMAND ddx_driver
    "${CMAKE_SOURCE_DIR}/tests/Input_cosmo.txt")

