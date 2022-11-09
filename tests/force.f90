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
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
type(ddx_state_type) :: state
integer :: iprint, info, pmax=30
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), &
    & hessianphi_cav(:, :, :), psi(:, :), &
    & force(:, :), force_num(:, :)
real(dp) :: tol, esolv1, esolv2, start_time, finish_time, step=0.0001, relerr
integer :: isph, i
real(dp), external :: dnrm2

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, ddx_data, tol, iprint, info)
call ddx_init_state(ddx_data % params, ddx_data % constants, state)
if(info .ne. 0) stop "info != 0"
allocate(phi_cav(ddx_data % constants % ncav), gradphi_cav(3, ddx_data % constants % ncav), &
    & hessianphi_cav(3, 3, ddx_data % constants % ncav), &
    & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), &
    & force(3, ddx_data % params % nsph), &
    & force_num(3, ddx_data % params % nsph))
call mkrhs(ddx_data % params, ddx_data % constants, ddx_data % workspace, 1, &
    & phi_cav, 1, gradphi_cav, 1, hessianphi_cav, psi)
call ddsolve(ddx_data, state, phi_cav, gradphi_cav, hessianphi_cav, psi, &
    & tol, esolv1, force, info)
do isph = 1, ddx_data % params % nsph
    do i = 1, 3
        ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) + step
        call solve(ddx_data, state, tol, esolv1, phi_cav, gradphi_cav, &
            & hessianphi_cav, psi, force, info)
        ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) - two*step
        call solve(ddx_data, state, tol, esolv2, phi_cav, gradphi_cav, &
            & hessianphi_cav, psi, force, info)
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

deallocate(phi_cav, gradphi_cav, psi, force, force_num)
call ddx_free_state(state)
call ddfree(ddx_data)

write(*, *) "Rel.error of forces:", relerr
if (relerr .gt. 3d-6) stop 1
contains 

subroutine solve(ddx_data, state, tol, esolv, phi_cav, gradphi_cav, &
        & hessianphi_cav, psi, force, info)
    type(ddx_type), intent(inout) :: ddx_data
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: esolv, phi_cav(ddx_data % constants % ncav), &
        & gradphi_cav(3, ddx_data % constants % ncav), &
        & hessianphi_cav(3, 3, ddx_data % constants % ncav), &
        & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), &
        & force(3, ddx_data % params % nsph)
    integer, intent(out) :: info
    type(ddx_type) :: ddx_data2
    call ddinit(ddx_data % params % nsph, ddx_data % params % charge, ddx_data % params % csph(1, :), &
        & ddx_data % params % csph(2, :), ddx_data % params % csph(3, :), ddx_data % params % rsph, &
        & ddx_data % params % model, ddx_data % params % lmax, ddx_data % params % ngrid, 0, &
        & ddx_data % params % fmm, ddx_data % params % pm, ddx_data % params % pl, &
        & ddx_data % params % se, &
        & ddx_data % params % eta, ddx_data % params % eps, ddx_data % params % kappa, &
        & ddx_data % params % matvecmem, ddx_data % params % itersolver, ddx_data % params % maxiter, &
        & ddx_data % params % jacobi_ndiis, &
        & ddx_data % params % nproc, ddx_data2, info)
    call mkrhs(ddx_data2 % params, ddx_data2 % constants, ddx_data2 % workspace, &
        & 1, phi_cav, 1, gradphi_cav, 1, hessianphi_cav, psi)
    call ddsolve(ddx_data2, state, phi_cav, gradphi_cav, hessianphi_cav, psi, tol, esolv, &
        & force, info)
    call ddfree(ddx_data2)
end subroutine solve

end program main


