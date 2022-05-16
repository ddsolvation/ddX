# @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file tests/CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2021-02-25

# Override compile options for accuracy
add_compile_options(-O1)

set(TESTS
    "ddx_core.f90"
    "ddx_operators.f90"
    "force.f90"
    "force_ddlpb.f90")

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

