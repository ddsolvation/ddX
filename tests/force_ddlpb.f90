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
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
integer :: info

real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), &
                       & hessian_cav(:,:,:), psi(:, :), &
                       & force(:, :), force_num(:, :)
real(dp) :: esolv, esolv_plus_h, esolv_minus_h, &
          & step = 0.00001, relerr
integer :: isph, i
real(dp), external :: dnrm2

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, ddx_data, info)
if(info .ne. 0) stop "info != 0"

! Allocation for variable vectors
allocate(phi_cav(ddx_data % ncav), gradphi_cav(3, ddx_data % ncav), &
       & hessian_cav(3,3,ddx_data % ncav), &
       & psi(ddx_data % nbasis, ddx_data % nsph), &
       & force(3, ddx_data % nsph), &
       & force_num(3, ddx_data % nsph))

gradphi_cav = zero
phi_cav = zero
hessian_cav = zero
psi = zero
esolv = zero
force = zero
force_num = zero
relerr = zero
esolv_plus_h = zero
esolv_minus_h = zero


call mkrhs(ddx_data, phi_cav, gradphi_cav, hessian_cav, psi)

call ddlpb(ddx_data, phi_cav, gradphi_cav, hessian_cav, psi, esolv, force)

isph = 1; i = 3
write(*,*) 'Sphere : ', isph, ' , Dimension : ', i

ddx_data % csph(i, isph) = ddx_data % csph(i, isph) + step

call solve(ddx_data, esolv_plus_h)
ddx_data % csph(i, isph) = ddx_data % csph(i, isph) - two*step

call solve(ddx_data, esolv_minus_h)
ddx_data % csph(i, isph) = ddx_data % csph(i, isph) + step

write(*,*) 'esolv+h  : ', esolv_plus_h, ' , esolv-h : ', esolv_minus_h, ' , step :', step
force_num(i, isph) = (esolv_plus_h - esolv_minus_h) / two / step

relerr = abs(force_num(i,isph) -force(i,isph))/ &
           & abs(force(i,isph))

write(6,'(2A60)') 'Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') force(1,i), force(2,i), force(3,i), force_num(1,i), &
      & force_num(2,i), force_num(3,i)
end do

deallocate(phi_cav, gradphi_cav, hessian_cav, psi, force, force_num)

call ddfree(ddx_data)

write(*, *) "Rel.error of forces:", relerr
if (relerr .gt. 1d-6) stop 1
contains

subroutine solve(ddx_data, esolv_in)
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(inout) :: esolv_in

    type(ddx_type) :: ddx_data2
    real(dp), allocatable :: phi_cav2(:)
    real(dp), allocatable :: gradphi_cav2(:,:)
    real(dp), allocatable :: hessian_cav2(:,:,:)
    real(dp), allocatable :: psi2(:,:)
    real(dp), allocatable :: force2(:,:)


    call ddinit(ddx_data % nsph, ddx_data % charge, ddx_data % csph(1, :), &
        & ddx_data % csph(2, :), ddx_data % csph(3, :), ddx_data % rsph, &
        & ddx_data % model, ddx_data % lmax, ddx_data % ngrid, 0, &
        & ddx_data % fmm, ddx_data % pm, ddx_data % pl, &
        & ddx_data % fmm_precompute, ddx_data % iprint, ddx_data % se, &
        & ddx_data % eta, ddx_data % eps, ddx_data % kappa, &
        & ddx_data % itersolver, ddx_data % tol, ddx_data % maxiter, &
        & ddx_data % ndiis, ddx_data % nproc, ddx_data2, info)

    allocate(phi_cav2(ddx_data2 % ncav), gradphi_cav2(3, ddx_data2 % ncav), &
            & hessian_cav2(3, 3, ddx_data2 % ncav), &
            & psi2(ddx_data2 % nbasis, ddx_data2 % nsph), &
            & force2(3, ddx_data2 % nsph))

    gradphi_cav2 = zero; phi_cav2 = zero
    hessian_cav2 = zero; psi2 = zero; force2 = zero

    call mkrhs(ddx_data2, phi_cav2, gradphi_cav2, hessian_cav2, psi2)
    call ddlpb(ddx_data2, phi_cav2, gradphi_cav2, hessian_cav2, psi2, esolv_in, force2)
    call ddfree(ddx_data2)
    deallocate(phi_cav2, gradphi_cav2, hessian_cav2, psi2, force2)
    return
end subroutine solve

end program main


