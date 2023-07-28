!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/force.f90
!! Test of analytical forces against numerical for ddLPB
!!
!! @version 1.0.0
!! @author Abhinav Jha
!! @date 2021-02-25

program main
use ddx_core
use ddx_operators
use ddx_solvers
use ddx
use ddx_lpb
use ddx_legacy
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
type(ddx_error_type) :: error
type(ddx_state_type) :: state

real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), &
                       & hessianphi_cav(:,:,:), psi(:, :), &
                       & force(:, :), force_num(:, :), charges(:)
real(dp) :: esolv, esolv_plus_h, esolv_minus_h, &
          & step = 2.d-5, relerr, tol
integer :: isph, i
real(dp), external :: dnrm2
character(len=255) :: dummy_file_name = ''

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, ddx_data, tol, charges, error)
call check_error(error)
call allocate_state(ddx_data % params, ddx_data % constants, state, error)
call check_error(error)

! Allocation for variable vectors
allocate(phi_cav(ddx_data % constants % ncav), gradphi_cav(3, ddx_data % constants % ncav), &
       & hessianphi_cav(3,3,ddx_data % constants % ncav), &
       & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), &
       & force(3, ddx_data % params % nsph), &
       & force_num(3, ddx_data % params % nsph))

gradphi_cav = zero
phi_cav = zero
hessianphi_cav = zero
psi = zero
esolv = zero
force = zero
force_num = zero
relerr = zero
esolv_plus_h = zero
esolv_minus_h = zero


call mkrhs(ddx_data % params, ddx_data % constants, ddx_data % workspace, &
    & 1, phi_cav, 1, gradphi_cav, 1, hessianphi_cav, psi, charges)
gradphi_cav = - gradphi_cav

call ddlpb(ddx_data % params, ddx_data % constants, ddx_data % workspace, &
    & state, phi_cav, gradphi_cav, psi, tol, esolv, hessianphi_cav, force, error)
call check_error(error)

! add the solute specific contributions to the forces
call grad_phi_for_charges(ddx_data % params, &
    & ddx_data % constants, ddx_data % workspace, state, &
    & charges, force, error)
call check_error(error)
call grad_e_for_charges(ddx_data % params, ddx_data % constants, &
    & ddx_data % workspace, state, charges, force, error)
call check_error(error)

do isph = 1, ddx_data % params % nsph
  do i = 1, 3
    esolv_plus_h = zero
    esolv_minus_h = zero
    ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) + step
    call test_solve(ddx_data, esolv_plus_h, tol, charges)
    ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) - step
    !call test_solve(ddx_data, esolv_minus_h)
    !write(*,*) 'esolv  :', esolv, 'esolv+h  : ', esolv_plus_h, ' , esolv : ', esolv_minus_h, ' , step :', step
    force_num(i, isph) = (esolv_plus_h - esolv) / step
    if(abs(force_num(i,isph)) .le. 1.d-8) force_num(i,isph) = zero
    if(abs(force(i,isph)) .le. 1.d-8) force(i,isph) = zero
  end do
end do
! force is minus the gradient
relerr = dnrm2(3*ddx_data % params % nsph, force_num-force, 1) / &
    & dnrm2(3*ddx_data % params % nsph, force, 1)

write(6,'(2A60)') 'Analytical forces', 'Numerical forces'
do i = 1, ddx_data % params % nsph
  write(6,'(6E20.10)') force(1,i), force(2,i), force(3,i), force_num(1,i), &
      & force_num(2,i), force_num(3,i)
end do

deallocate(phi_cav, gradphi_cav, hessianphi_cav, psi, force, force_num, charges)

call deallocate_state(state, error)
call deallocate_model(ddx_data, error)

write(*, *) "Rel.error of forces:", relerr
if (relerr .gt. 1.d-5) stop 1
contains

subroutine test_solve(ddx_data, esolv_in, tol, charges)
    implicit none
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(inout) :: esolv_in
    real(dp), intent(in) ::tol
    real(dp), intent(in) :: charges(ddx_data % params % nsph)

    type(ddx_type) :: ddx_data2
    type(ddx_error_type) :: error2
    type(ddx_state_type) :: state
    real(dp), allocatable :: phi_cav2(:)
    real(dp), allocatable :: gradphi_cav2(:,:)
    real(dp), allocatable :: hessianphi_cav2(:,:,:)
    real(dp), allocatable :: psi2(:,:)
    real(dp), allocatable :: force2(:,:)

    call allocate_model(ddx_data % params % nsph, ddx_data % params % csph(1, :), &
        & ddx_data % params % csph(2, :), ddx_data % params % csph(3, :), ddx_data % params % rsph, &
        & ddx_data % params % model, ddx_data % params % lmax, ddx_data % params % ngrid, 0, &
        & ddx_data % params % fmm, ddx_data % params % pm, ddx_data % params % pl, &
        & ddx_data % params % se, &
        & ddx_data % params % eta, ddx_data % params % eps, ddx_data % params % kappa, &
        & ddx_data % params % matvecmem, ddx_data % params % maxiter, &
        & ddx_data % params % jacobi_ndiis, &
        & ddx_data % params % nproc, dummy_file_name, ddx_data2, error2)
    call check_error(error2)

    call allocate_state(ddx_data2 % params, ddx_data2 % constants, state, error2)
    call check_error(error2)

    allocate(phi_cav2(ddx_data2 % constants % ncav), gradphi_cav2(3, ddx_data2 % constants % ncav), &
            & hessianphi_cav2(3, 3, ddx_data2 % constants % ncav), &
            & psi2(ddx_data2 % constants % nbasis, ddx_data2 % params % nsph), &
            & force2(3, ddx_data2 % params % nsph))

    gradphi_cav2 = zero; phi_cav2 = zero
    hessianphi_cav2 = zero; psi2 = zero; force2 = zero

    call mkrhs(ddx_data2 % params, ddx_data2 % constants, ddx_data2 % workspace, &
        & 1, phi_cav2, 1, gradphi_cav2, 1, hessianphi_cav2, psi2, charges)
    gradphi_cav2 = - gradphi_cav2
    call ddsolve_legacy(ddx_data2, state, phi_cav2, gradphi_cav2, hessianphi_cav2, &
        & psi2, tol, esolv_in, force2, error2)
    call check_error(error2)

    call deallocate_state(state, error2)
    call deallocate_model(ddx_data2, error2)
    deallocate(phi_cav2, gradphi_cav2, hessianphi_cav2, psi2, force2)
end subroutine test_solve

end program main


