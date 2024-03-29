!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/force.f90
!! Test of analytical forces against numerical
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

program main
use ddx_core
use ddx_operators
use ddx_solvers
use ddx
use ddx_legacy
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
type(ddx_error_type) :: ddx_error
type(ddx_state_type) :: state
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), &
    & hessianphi_cav(:, :, :), psi(:, :), &
    & force(:, :), force_num(:, :), charges(:)
real(dp) :: tol, esolv1, esolv2, step=0.0001, relerr
integer :: isph, i
real(dp), external :: dnrm2
character(len=255) :: dummy_file_name = ''

! Read input file name
call get_command_argument(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, ddx_data, tol, charges, ddx_error)
call check_error(ddx_error)
call allocate_state(ddx_data % params, ddx_data % constants, state, ddx_error)
call check_error(ddx_error)

allocate(phi_cav(ddx_data % constants % ncav), gradphi_cav(3, ddx_data % constants % ncav), &
    & hessianphi_cav(3, 3, ddx_data % constants % ncav), &
    & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), &
    & force(3, ddx_data % params % nsph), &
    & force_num(3, ddx_data % params % nsph))
call mkrhs(ddx_data % params, ddx_data % constants, ddx_data % workspace, 1, &
    & phi_cav, 1, gradphi_cav, 1, hessianphi_cav, psi, charges)
call ddsolve_legacy(ddx_data, state, phi_cav, -gradphi_cav, hessianphi_cav, psi, &
    & tol, esolv1, force, ddx_error)
call check_error(ddx_error)
call grad_phi_for_charges(ddx_data % params, ddx_data % constants, &
    & ddx_data % workspace, state, charges, force, ddx_error)
call check_error(ddx_error)

ddx_data % params % force = 0
do isph = 1, ddx_data % params % nsph
    do i = 1, 3
        ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) + step
        call test_solve(ddx_data, state, tol, esolv1, phi_cav, gradphi_cav, &
            & hessianphi_cav, psi, force, charges)
        ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) - two*step
        call test_solve(ddx_data, state, tol, esolv2, phi_cav, gradphi_cav, &
            & hessianphi_cav, psi, force, charges)
        ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) + step
        force_num(i, isph) = (esolv1-esolv2) / two / step
    end do
end do
relerr = dnrm2(3*ddx_data % params % nsph, force_num-force, 1) / &
    & dnrm2(3*ddx_data % params % nsph, force, 1)

write(6,'(2A60)') 'Analytical forces', 'Numerical forces'
do i = 1, ddx_data % params % nsph
  write(6,'(6E20.10)') force(1,i), force(2,i), force(3,i), force_num(1,i), &
      & force_num(2,i), force_num(3,i)
end do

deallocate(phi_cav, gradphi_cav, psi, force, force_num, charges)
call deallocate_state(state, ddx_error)
call deallocate_model(ddx_data, ddx_error)

write(*, *) "Rel.ddx_error of forces:", relerr
if (relerr .gt. 1d-5) stop 1
contains 

subroutine test_solve(ddx_data, state, tol, esolv, phi_cav, gradphi_cav, &
        & hessianphi_cav, psi, force, charges)
    type(ddx_type), intent(inout) :: ddx_data
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: charges(ddx_data % params % nsph)
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: esolv, phi_cav(ddx_data % constants % ncav), &
        & gradphi_cav(3, ddx_data % constants % ncav), &
        & hessianphi_cav(3, 3, ddx_data % constants % ncav), &
        & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), &
        & force(3, ddx_data % params % nsph)
    type(ddx_type) :: ddx_data2
    type(ddx_error_type) :: error2
    call allocate_model(ddx_data % params % nsph, ddx_data % params % csph(1, :), &
        & ddx_data % params % csph(2, :), ddx_data % params % csph(3, :), ddx_data % params % rsph, &
        & ddx_data % params % model, ddx_data % params % lmax, ddx_data % params % ngrid, 0, &
        & ddx_data % params % fmm, ddx_data % params % pm, ddx_data % params % pl, &
        & ddx_data % params % se, &
        & ddx_data % params % eta, ddx_data % params % eps, ddx_data % params % kappa, &
        & ddx_data % params % matvecmem, ddx_data % params % maxiter, &
        & ddx_data % params % jacobi_ndiis, &
        & ddx_data % params % nproc, dummy_file_name, ddx_data2, error2)
    call mkrhs(ddx_data2 % params, ddx_data2 % constants, ddx_data2 % workspace, &
        & 1, phi_cav, 1, gradphi_cav, 1, hessianphi_cav, psi, charges)
    call ddsolve_legacy(ddx_data2, state, phi_cav, -gradphi_cav, hessianphi_cav, psi, tol, esolv, &
        & force, error2)
    call check_error(error2)
    call deallocate_model(ddx_data2, error2)
end subroutine test_solve

end program main


