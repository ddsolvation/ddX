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
integer :: info, pmax=30
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), psi(:, :), &
    & force(:, :), force_num(:, :)
real(dp) :: esolv1, esolv2, start_time, finish_time, step=0.0001, relerr
integer :: isph, i
real(dp), external :: dnrm2

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, ddx_data, info)
if(info .ne. 0) stop "info != 0"
allocate(phi_cav(ddx_data % ncav), gradphi_cav(3, ddx_data % ncav), &
    & psi(ddx_data % nbasis, ddx_data % nsph), force(3, ddx_data % nsph), &
    & force_num(3, ddx_data % nsph))
call mkrhs(ddx_data, phi_cav, gradphi_cav, psi)
call ddsolve(ddx_data, phi_cav, gradphi_cav, psi, esolv1, force)
do isph = 1, ddx_data % nsph
    do i = 1, 3
        ddx_data % csph(i, isph) = ddx_data % csph(i, isph) + step
        call solve(ddx_data, esolv1, phi_cav, gradphi_cav, psi, force)
        ddx_data % csph(i, isph) = ddx_data % csph(i, isph) - two*step
        call solve(ddx_data, esolv2, phi_cav, gradphi_cav, psi, force)
        ddx_data % csph(i, isph) = ddx_data % csph(i, isph) + step
        force_num(i, isph) = (esolv1-esolv2) / two / step
    end do
end do
relerr = dnrm2(3*ddx_data % nsph, force_num-force, 1) / &
    & dnrm2(3*ddx_data % nsph, force, 1)

write(6,'(2A60)') 'Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') force(1,i), force(2,i), force(3,i), force_num(1,i), &
      & force_num(2,i), force_num(3,i)
end do

deallocate(phi_cav, gradphi_cav, psi, force, force_num)
call ddfree(ddx_data)

write(*, *) "Rel.error of forces:", relerr
if (relerr .gt. 1d-6) stop 1
contains 

subroutine solve(ddx_data, esolv, phi_cav, gradphi_cav, psi, force)
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(out) :: esolv, phi_cav(ddx_data % ncav), &
        & gradphi_cav(3, ddx_data % ncav), &
        & psi(ddx_data % nbasis, ddx_data % nsph), force(3, ddx_data % nsph)
    type(ddx_type) :: ddx_data2
    call ddinit(ddx_data % nsph, ddx_data % charge, ddx_data % csph(1, :), &
        & ddx_data % csph(2, :), ddx_data % csph(3, :), ddx_data % rsph, &
        & ddx_data % model, ddx_data % lmax, ddx_data % ngrid, 0, &
        & ddx_data % fmm, ddx_data % pm, ddx_data % pl, &
        & ddx_data % fmm_precompute, ddx_data % iprint, ddx_data % se, &
        & ddx_data % eta, ddx_data % eps, ddx_data % kappa, &
        & ddx_data % itersolver, ddx_data % tol, ddx_data % maxiter, &
        & ddx_data % ndiis, ddx_data % nproc, ddx_data2, info)
    call mkrhs(ddx_data2, phi_cav, gradphi_cav, psi)
    call ddsolve(ddx_data2, phi_cav, gradphi_cav, psi, esolv, force)
    call ddfree(ddx_data2)
end subroutine solve

end program main


