!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/ddx_operators.f90
!! Tests for ddx_operators module
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-24

program test_ddx_operators
use ddx_operators
implicit none

integer :: i, iprint=1, info, ngrid=590, nproc=1
! Some implementation do not work for alpha=1d-307, so range of test values is
! reduced. This can be fixed by enforcing proper order of multiplications like
! a*b*c into a*(b*c).
real(dp) :: alpha(4)=(/1d0, -1d0, 1d-100, 1d+100/)
type(ddx_type) :: ddx_data
integer, parameter :: nsph=10, lmax=7, force=1, itersolver=1, &
    & maxiter=1000, jacobi_ndiis=25, gmresr_j=1, gmresr_dim=10
real(dp), parameter :: se=0d0, eta=0.1d0, eps=78d0, kappa=0d0
real(dp) :: gcsph(3, nsph), csph(3, nsph), grsph(nsph), rsph(nsph), &
    & gcharge(nsph), charge(nsph)

gcsph(:, 1) = (/1d0, 1d0, 1d0/)
gcsph(:, 2) = (/2d0, 2d0, 2d0/)
gcsph(:, 3) = (/1d0, 1d0, 3d0/)
gcsph(:, 4) = (/2d0, 2d0, 4d0/)
gcsph(:, 5) = (/1d0, 1d0, 5d0/)
gcsph(:, 6) = (/2d0, 2d0, 6d0/)
gcsph(:, 7) = (/1d0, 1d0, 7d0/)
gcsph(:, 8) = (/2d0, 2d0, 8d0/)
gcsph(:, 9) = (/1d0, 1d0, 9d0/)
gcsph(:, 10) = (/2d0, 2d0, 10d0/)
grsph = (/1d-1, 2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, 1d0/)
gcharge = one

do i = 1, size(alpha)
    charge = abs(alpha(i)) * gcharge
    csph = alpha(i) * gcsph
    rsph = abs(alpha(i)) * grsph
    call ddinit(nsph, charge, csph(1, :), csph(2, :), csph(3, :), rsph, 2, &
        lmax, ngrid, force, 0, -1, -1, se, eta, eps, kappa, &
        & itersolver, maxiter, jacobi_ndiis, gmresr_j, gmresr_dim, &
        & nproc, ddx_data, info)
    if(info .ne. 0) stop 1
    call check_mkrhs(ddx_data, 0, 0, iprint, 1d-1)
    call check_mkrhs(ddx_data, 1, 1, iprint, 1d-2)
    call check_mkrhs(ddx_data, 3, 3, iprint, 1d-3)
    call check_mkrhs(ddx_data, 5, 5, iprint, 1d-4)
    call check_mkrhs(ddx_data, 20, 20, iprint, 1d-9)
    call check_mkrhs(ddx_data, 40, 40, iprint, 1d-15)
    call check_dx(ddx_data, lmax, lmax, iprint, 1d-4)
    call check_dx(ddx_data, 40, 40, iprint, 1d-15)
    call check_gradr(ddx_data, lmax, lmax, iprint, 1d-4)
    call check_gradr(ddx_data, lmax+1, lmax+1, iprint, 1d-4)
    call check_gradr(ddx_data, 40, 40, iprint, 1d-15)
    call ddfree(ddx_data)
end do

contains

subroutine check_mkrhs(ddx_data, pm, pl, iprint, threshold)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    integer, intent(in) :: pm, pl, iprint
    real(dp), intent(in) :: threshold
    ! Local variables
    type(ddx_type) :: ddx_data_fmm
    integer :: info
    real(dp), allocatable :: phi_cav(:), phi2_cav(:), gradphi_cav(:, :), &
        & gradphi2_cav(:, :), hessianphi_cav(:, :, :), &
        & hessianphi2_cav(:, :, :), psi(:, :), psi2(:, :), force(:, :)
    real(dp) :: fnorm, fdiff
    real(dp), external :: dnrm2
    ! Init FMM-related ddx_data
    call ddinit(ddx_data % params % nsph, ddx_data % params % charge, ddx_data % params % csph(1, :), &
        & ddx_data % params % csph(2, :), ddx_data % params % csph(3, :), ddx_data % params % rsph, &
        & ddx_data % params % model, ddx_data % params % lmax, ddx_data % params % ngrid, ddx_data % params % force, &
        & 1, pm, pl, ddx_data % params % se, ddx_data % params % eta, &
        & ddx_data % params % eps, ddx_data % params % kappa, ddx_data % params % itersolver, &
        & ddx_data % params % maxiter, ddx_data % params % jacobi_ndiis, &
        & ddx_data % params % gmresr_j, ddx_data % params % gmresr_dim, &
        & ddx_data % params % nproc, &
        & ddx_data_fmm, info)
    ! Allocate resources
    allocate(phi_cav(ddx_data % constants % ncav), &
        & phi2_cav(ddx_data % constants % ncav), &
        & gradphi_cav(3, ddx_data % constants % ncav), &
        & gradphi2_cav(3, ddx_data % constants % ncav), &
        & hessianphi_cav(3, 3, ddx_data % constants % ncav), &
        & hessianphi2_cav(3, 3, ddx_data % constants % ncav), &
        & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), &
        & psi2(ddx_data % constants % nbasis, ddx_data % params % nsph), &
        & force(3, ddx_data % params % nsph), &
        & stat=info)
    if(info .ne. 0) call error(-1, "Allocation failed")
    ! Dense operator mkrhs is trusted to have no errors, this must be somehow
    ! checked in the future.
    call mkrhs(ddx_data, 1, phi_cav, 1, gradphi_cav, 1, hessianphi_cav, psi)
    call mkrhs(ddx_data_fmm, 1, phi2_cav, 1, gradphi2_cav, 1, &
        & hessianphi2_cav, psi2)
    ! Compare potentials
    phi2_cav = phi2_cav - phi_cav
    gradphi2_cav = gradphi2_cav - gradphi_cav
    hessianphi2_cav = hessianphi2_cav - hessianphi_cav
    fnorm = dnrm2(ddx_data % constants % ncav, phi_cav, 1)
    fdiff = dnrm2(ddx_data % constants % ncav, phi2_cav, 1)
    print *, "Pot ", fdiff, fnorm, fdiff/fnorm
    if (fdiff .gt. threshold*fnorm) then
        call error(-1, "Potentials are different")
    end if
    fnorm = dnrm2(3*ddx_data % constants % ncav, gradphi_cav, 1)
    fdiff = dnrm2(3*ddx_data % constants % ncav, gradphi2_cav, 1)
    print *, "Grad", fdiff, fnorm, fdiff/fnorm
    if (fdiff .gt. threshold*fnorm) then
        call error(-1, "Gradients are different")
    end if
    fnorm = dnrm2(9*ddx_data % constants % ncav, hessianphi_cav, 1)
    fdiff = dnrm2(9*ddx_data % constants % ncav, hessianphi2_cav, 1)
    print *, "Hess", fdiff, fnorm, fdiff/fnorm
    if (fdiff .gt. threshold*fnorm) then
        call error(-1, "Hessians are different")
    end if
end subroutine check_mkrhs

subroutine check_dx(ddx_data, pm, pl, iprint, threshold)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    integer, intent(in) :: pm, pl, iprint
    real(dp), intent(in) :: threshold
    ! Local variables
    type(ddx_type) :: ddx_data_fmm
    integer :: info, irand, iseed(4)=(/0, 0, 0, 1/), do_diag
    integer, parameter :: nrand=10
    real(dp) :: x(ddx_data % constants % nbasis, ddx_data % params % nsph, nrand), &
        & y(ddx_data % constants % nbasis, ddx_data % params % nsph, nrand), &
        & z(ddx_data % constants % nbasis, ddx_data % params % nsph, nrand), &
        & xx(nrand, nrand), yy(nrand, nrand), full_norm, diff_norm, &
        & forces(3, ddx_data % params % nsph), forces2(3, ddx_data % params % nsph)
    real(dp), external :: dnrm2
    ! Init FMM-related ddx_data
    call ddinit(ddx_data % params % nsph, ddx_data % params % charge, ddx_data % params % csph(1, :), &
        & ddx_data % params % csph(2, :), ddx_data % params % csph(3, :), ddx_data % params % rsph, &
        & ddx_data % params % model, ddx_data % params % lmax, ddx_data % params % ngrid, ddx_data % params % force, &
        & 1, pm, pl, ddx_data % params % se, ddx_data % params % eta, &
        & ddx_data % params % eps, ddx_data % params % kappa, ddx_data % params % itersolver, &
        & ddx_data % params % maxiter, ddx_data % params % jacobi_ndiis, &
        & ddx_data % params % gmresr_j, ddx_data % params % gmresr_dim, &
        & ddx_data % params % nproc, &
        & ddx_data_fmm, info)
    ! Dense operator dx is trusted to have no errors, this must be somehow
    ! checked in the future.
    ! Get random x
    call dlarnv(2, iseed, ddx_data % constants % n * nrand, x)
    write(*, *) "pm=", pm, "pl=", pl
    do do_diag = 0, 1
        write(*, *) "do_diag=", do_diag
        ! Random check of FMM dx operator against dense dx operator
        do irand = 1, nrand
            call dx_dense(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, do_diag, x(:, :, irand), y(:, :, irand))
        end do
        full_norm = dnrm2(ddx_data % constants % n * nrand, y, 1)
        do irand = 1, nrand
            call dx_fmm(ddx_data_fmm % params, ddx_data_fmm % constants, &
                & ddx_data_fmm % workspace, do_diag, x(:, :, irand), z(:, :, irand))
        end do
        diff_norm = dnrm2(ddx_data % constants % n * nrand, y-z, 1)
        write(*, *) "dx_dense vs dx_fmm rel.error=", diff_norm/full_norm
        ! Check dense adjoint operator dstarx
        do irand = 1, nrand
            call dstarx_dense(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, do_diag, x(:, :, irand), y(:, :, irand))
        end do
        call dgemm('T', 'N', nrand, nrand, ddx_data % constants % n, one, y, ddx_data % constants % n, &
            & y, ddx_data % constants % n, zero, xx, nrand)
        full_norm = dnrm2(nrand**2, xx, 1)
        do irand = 1, nrand
            call dx_dense(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, do_diag, y(:, :, irand), z(:, :, irand))
        end do
        call dgemm('T', 'N', nrand, nrand, ddx_data % constants % n, one, z, ddx_data % constants % n, &
            & x, ddx_data % constants % n, zero, yy, nrand)
        diff_norm = dnrm2(nrand**2, xx-yy, 1)
        write(*, *) "dstarx_dense vs dx_dense rel.error=", diff_norm/full_norm
        ! Check FMM adjoint operator dstarx (without precomputed FMM matrices)
        do irand = 1, nrand
            call dstarx_fmm(ddx_data_fmm % params, ddx_data_fmm % constants, &
                & ddx_data_fmm % workspace, do_diag, x(:, :, irand), y(:, :, irand))
        end do
        call dgemm('T', 'N', nrand, nrand, ddx_data % constants % n, one, y, ddx_data % constants % n, &
            & y, ddx_data % constants % n, zero, xx, nrand)
        full_norm = dnrm2(nrand**2, xx, 1)
        do irand = 1, nrand
            call dx_fmm(ddx_data_fmm % params, ddx_data_fmm % constants, &
                & ddx_data_fmm % workspace, do_diag, y(:, :, irand), z(:, :, irand))
        end do
        call dgemm('T', 'N', nrand, nrand, ddx_data % constants % n, one, z, ddx_data % constants % n, &
            & x, ddx_data % constants % n, zero, yy, nrand)
        diff_norm = dnrm2(nrand**2, xx-yy, 1)
        write(*, *) "dstarx_fmm vs dx_fmm rel.error=", diff_norm/full_norm
    end do
    ! Free temporary objects
    call ddfree(ddx_data_fmm)
end subroutine check_dx

subroutine check_gradr(ddx_data, pm, pl, iprint, threshold)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    integer, intent(in) :: pm, pl, iprint
    real(dp), intent(in) :: threshold
    ! Local variables
    type(ddx_type) :: ddx_data_fmm
    integer :: info, irand, iseed(4)=(/0, 0, 0, 1/), do_diag
    integer, parameter :: nrand=10
    real(dp) :: x(ddx_data % constants % nbasis, ddx_data % params % nsph, nrand), &
        & y(ddx_data % constants % nbasis, ddx_data % params % nsph, nrand), &
        & z(ddx_data % constants % nbasis, ddx_data % params % nsph, nrand), &
        & xx(nrand, nrand), yy(nrand, nrand), full_norm, diff_norm, &
        & forces(3, ddx_data % params % nsph), forces2(3, ddx_data % params % nsph)
    real(dp), external :: dnrm2
    ! Init FMM-related ddx_data
    call ddinit(ddx_data % params % nsph, ddx_data % params % charge, ddx_data % params % csph(1, :), &
        & ddx_data % params % csph(2, :), ddx_data % params % csph(3, :), ddx_data % params % rsph, &
        & ddx_data % params % model, ddx_data % params % lmax, ddx_data % params % ngrid, ddx_data % params % force, &
        & 1, pm, pl, ddx_data % params % se, ddx_data % params % eta, &
        & ddx_data % params % eps, ddx_data % params % kappa, ddx_data % params % itersolver, &
        & ddx_data % params % maxiter, ddx_data % params % jacobi_ndiis, &
        & ddx_data % params % gmresr_j, ddx_data % params % gmresr_dim, &
        & ddx_data % params % nproc, &
        & ddx_data_fmm, info)
    ! Dense operator dx is trusted to have no errors, this must be somehow
    ! checked in the future.
    ! Get random ygrid and g
    call dlarnv(2, iseed, ddx_data % params % ngrid * ddx_data % params % nsph, ddx_data % ygrid)
    call dlarnv(2, iseed, ddx_data % constants % n, ddx_data % g)
    ! Check gradr
    call gradr_dense(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, ddx_data % g, ddx_data % ygrid, forces)
    full_norm = dnrm2(3*ddx_data % params % nsph, forces, 1)
    call gradr_fmm(ddx_data_fmm % params, ddx_data_fmm % constants, &
        & ddx_data_fmm % workspace, ddx_data % g, ddx_data % ygrid, forces2)
    diff_norm = dnrm2(3*ddx_data % params % nsph, forces-forces2, 1)
    write(*, *) "gradr dense vs fmm rel.error=", diff_norm / full_norm
    call ddfree(ddx_data_fmm)
end subroutine check_gradr

end program test_ddx_operators

