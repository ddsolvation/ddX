!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/ddx_core.f90
!! Tests for ddx_core module
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-11

program test_ddx_core
use ddx_core
implicit none

!character(len=255) :: testname, iprint_string
integer :: argc, iprint=1
integer :: p=10, i

! Input points and spheres
integer, parameter :: nsph = 10
real(dp) :: gsrc0(3), gsrc(3), gdst(3, 2), gcsph(3, nsph), grsph(nsph), &
    & gdst_csph(3, 2), gdst_rsph(2)
! Some implementation do not work for alpha=1d-307, so range of test values is
! reduced. This can be fixed by enforcing proper order of multiplications like
! a*b*c into a*(b*c).
real(dp) :: alpha(4)=(/1d0, -1d0, 1d-294, 1d+294/)

real(dp), external :: dnrm2
! Set inputs
gsrc0 = zero
gsrc = (/1d-1, 2d-1, 3d-1/)
! This test shall be aligned along OZ axis
gdst(:, 1) = (/-0.2d0, 0.1d0, 25.2d0/)
gdst_csph(:, 1) = (/0d0, 0d0, 25d0/)
gdst_rsph(1) = pt5
! This test shall NOT be aligned along OZ axis
gdst(:, 2) = (/-16.3d0, 20.2d0, 0.1d0/)
gdst_csph(:, 2) = (/-16d0, 20d0, 0d0/)
gdst_rsph(2) = pt5
gcsph(:, 1) = zero
grsph(1) = one
gcsph(:, 2) = (/0d0, 0d0, -1.1d0/)
grsph(2) = grsph(1) + dnrm2(3, gcsph(:, 2)-gcsph(:, 1), 1)
gcsph(:, 3) = (/0d0, 0d0, 7d-1/)
grsph(3) = grsph(1) + dnrm2(3, gcsph(:, 3)-gcsph(:, 1), 1)
gcsph(:, 4) = -pt5
grsph(4) = grsph(1) + dnrm2(3, gcsph(:, 4)-gcsph(:, 1), 1)
gcsph(:, 5) = pt5
grsph(5) = grsph(1) + dnrm2(3, gcsph(:, 4)-gcsph(:, 1), 1)
gcsph(:, 6) = zero
grsph(6) = three
gcsph(:, 7) = (/1d-1, 2d-1, 1.1d0/)
grsph(7) = grsph(1) + dnrm2(3, gcsph(:, 7)-gcsph(:, 1), 1)
gcsph(:, 8) = (/4d-1, 2d-9, 1.1d0/)
grsph(8) = grsph(1) + dnrm2(3, gcsph(:, 8)-gcsph(:, 1), 1)
gcsph(:, 9) = (/1.1d0, 0d0, 0d0/)
grsph(9) = grsph(1) + dnrm2(3, gcsph(:, 9)-gcsph(:, 1), 1)
gcsph(:, 10) = (/-4d-1, 1.1d0, 0d0/)
grsph(10) = grsph(1) + dnrm2(3, gcsph(:, 10)-gcsph(:, 1), 1)

! Check correctness of info for valid and invalid input parameters of ddinit
call check_ddinit_args()

! Check P2M and M2P operations of the FMM
do i = 1, size(alpha)
    call check_p2m_m2p(0, alpha(i), iprint, 6d-2)
    call check_p2m_m2p(1, alpha(i), iprint, 3d-3)
    call check_p2m_m2p(p, alpha(i), iprint, 120*epsilon(zero))
end do

! Check P2L and L2P operations of the FMM
do i = 1, size(alpha)
    call check_p2l_l2p(0, alpha(i), iprint, 6d-2)
    call check_p2l_l2p(1, alpha(i), iprint, 3d-3)
    call check_p2l_l2p(p, alpha(i), iprint, 120*epsilon(zero))
end do

! Check M2M operations of the FMM
do i = 1, size(alpha)
    call check_m2m(0, alpha(i), iprint, 10d0*epsilon(zero))
    call check_m2m(1, alpha(i), iprint, 10d0*epsilon(zero))
    call check_m2m(p, alpha(i), iprint, 20d0*epsilon(zero))
end do

! Check L2L operations of the FMM
do i = 1, size(alpha)
    call check_l2l(0, alpha(i), iprint, 10d0*epsilon(zero))
    call check_l2l(1, alpha(i), iprint, 10d0*epsilon(zero))
    call check_l2l(p, alpha(i), iprint, 40d0*epsilon(zero))
end do

! Check M2L operations of the FMM
do i = 1, size(alpha)
    call check_m2l(0, 0, alpha(i), iprint, 6d-2)
    call check_m2l(0, 1, alpha(i), iprint, 6d-2)
    call check_m2l(0, p, alpha(i), iprint, 6d-2)
    call check_m2l(1, 0, alpha(i), iprint, 2d-2)
    call check_m2l(1, 1, alpha(i), iprint, 3d-3)
    call check_m2l(1, p, alpha(i), iprint, 3d-3)
    call check_m2l(p, 0, alpha(i), iprint, 2d-2)
    call check_m2l(p, 1, alpha(i), iprint, 2d-4)
    call check_m2l(p, p, alpha(i), iprint, 20d0*epsilon(zero))
end do

! Check recursive inertial tree
do i = 1, size(alpha)
    call check_tree_rib(alpha(i))
end do

! Check tree M2M
do i = 1, size(alpha)
    call check_tree_m2m(0, alpha(i), iprint, 10d0*epsilon(one))
    call check_tree_m2m(1, alpha(i), iprint, 10d0*epsilon(one))
    call check_tree_m2m(10, alpha(i), iprint, 10d0*epsilon(one))
end do

! Check tree L2L
do i = 1, size(alpha)
    call check_tree_l2l(0, alpha(i), iprint, 10d0*epsilon(one))
    call check_tree_l2l(1, alpha(i), iprint, 10d0*epsilon(one))
    call check_tree_l2l(10, alpha(i), iprint, 10d0*epsilon(one))
end do

! Check tree M2L
do i = 1, size(alpha)
    call check_tree_m2l(0, 0, alpha(i), iprint, 6d-2)
    call check_tree_m2l(1, 0, alpha(i), iprint, 4d-2)
    call check_tree_m2l(10, 0, alpha(i), iprint, 4d-2)
    call check_tree_m2l(0, 1, alpha(i), iprint, 3d-2)
    call check_tree_m2l(1, 1, alpha(i), iprint, 4d-3)
    call check_tree_m2l(10, 1, alpha(i), iprint, 4d-3)
    call check_tree_m2l(0, 10, alpha(i), iprint, 3d-2)
    call check_tree_m2l(1, 10, alpha(i), iprint, 3d-3)
    call check_tree_m2l(10, 10, alpha(i), iprint, 4d-9)
    !call check_tree_m2l(20, 20, alpha(i), iprint, 1d-14)
end do

contains

subroutine check_ddinit_args()
    ! Example of correct args
    integer :: n=1, model=1, lmax=0, ngrid=1000, force=1, fmm=1, pm=0, pl=0, &
        & fmm_precompute=0, iprint=0, itersolver=1, maxiter=10, ndiis=10, &
        & nproc=2
    real(dp) :: charge(10), x(10), y(10), z(10), rvdw(10), se=zero, eta=1d-1, &
        & eps=zero, kappa=zero, tol=1d-6
    type(ddx_type) :: ddx_data
    integer :: info=0, i, j
    real(dp) :: tmp
    ! Generate coordinates and radii
    rvdw = 4d0
    charge = one
    do i = 1, 10
        x(i) = dble(2*i)
        y(i) = x(i)
        z(i) = x(i)
    end do
    ! Check correct input
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop "ddinit tests failed: correct input"
    call ddfree(ddx_data)
    ! Check different correct inputs with different n <= 10 (hardcoded value)
    do i = 1, 10
        call ddinit(i, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, &
            & pm, pl, fmm_precompute, iprint, se, eta, eps, kappa, &
            & itersolver, tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop "Correct input failed"
        call ddfree(ddx_data)
    end do
    ! Check incorrect input nsph = 0
    i = 0
    call ddinit(i, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -1) stop "nsph=0 failed"
    call ddfree(ddx_data)
    ! Check all possible models with other correct inputs
    do i = 1, 3
        call ddinit(n, charge, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
            & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, &
            & tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check incorrect models
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -7) stop 1
    call ddfree(ddx_data)
    i = 4
    call ddinit(n, charge, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -7) stop 1
    call ddfree(ddx_data)
    ! Check correct lmax
    do i = 1, 6
        call ddinit(n, charge, x, y, z, rvdw, model, i, ngrid, force, fmm, &
            & pm, pl, fmm_precompute, iprint, se, eta, eps, kappa, &
            & itersolver, tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check incorrect lmax < 0
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, i, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -8) stop 1
    call ddfree(ddx_data)
    ! Check correct ngrid
    do i = 0, 1000, 100
        j = i
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, j, force, fmm, pm, &
            & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, &
            & tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check incorrect ngrid < 0
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, i, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -9) stop 1
    call ddfree(ddx_data)
    ! Check correct force
    do i = 0, 1
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
            & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, &
            & tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check incorrect force
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -10) stop 1
    call ddfree(ddx_data)
    i = 2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -10) stop 1
    call ddfree(ddx_data)
    ! Check correct fmm
    do i = 0, 1
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, i, &
            & pm, pl, fmm_precompute, iprint, se, eta, eps, kappa, &
            & itersolver, tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check incorrect fmm
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -11) stop 1
    call ddfree(ddx_data)
    i = 2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -11) stop 1
    call ddfree(ddx_data)
    ! Check correct pm (ignored if fmm=0)
    j = 0
    do i = -2, 2
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
            & i, pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, &
            & tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check correct pm (fmm=1)
    j = 1
    do i = 0, 20, 5
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
            & i, pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, &
            & tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check incorrect pm (fmm=1)
    j = 1
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, i, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -12) stop 1
    call ddfree(ddx_data)
    ! Check correct pl (ignored if fmm=0)
    j = 0
    do i = -2, 2
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
            & pm, i, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, &
            & tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check correct pl (fmm=1)
    j = 1
    do i = 0, 20, 5
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
            & pm, i, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, &
            & tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check incorrect pl (fmm=1)
    j = 1
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, pm, &
        & i, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -13) stop 1
    call ddfree(ddx_data)
    ! Check correct fmm_precompute
    do i = 0, 1
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, &
            & pm, pl, i, iprint, se, eta, eps, kappa, itersolver, tol, &
            & maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check incorrect fmm_precompute
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, i, iprint, se, eta, eps, kappa, itersolver, tol, maxiter, &
        & ndiis, nproc, ddx_data, info)
    if (info .ne. -14) stop 1
    call ddfree(ddx_data)
    i = 2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, i, iprint, se, eta, eps, kappa, itersolver, tol, maxiter, &
        & ndiis, nproc, ddx_data, info)
    if (info .ne. -14) stop 1
    call ddfree(ddx_data)
    ! Check correct iprint
    do i = 0, 10
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, &
            & pm, pl, fmm_precompute, i, se, eta, eps, kappa, itersolver, &
            & tol, maxiter, ndiis, nproc, ddx_data, info)
        if (info .ne. 0) stop 1
        call ddfree(ddx_data)
    end do
    ! Check incorrect iprint
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, i, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -15) stop 1
    call ddfree(ddx_data)
    ! Check correct se
    tmp = -one
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, tmp, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    tmp = zero
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, tmp, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    tmp = one
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, tmp, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    ! Check incorrect se
    tmp = 1.01d0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, tmp, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -16) stop 1
    call ddfree(ddx_data)
    tmp = -1.01d0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, tmp, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -16) stop 1
    call ddfree(ddx_data)
    ! Check correct eta
    tmp = zero
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, tmp, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    tmp = pt5
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, tmp, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    tmp = one
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, tmp, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    ! Check incorrect eta
    tmp = 1.01d0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, tmp, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -17) stop 1
    call ddfree(ddx_data)
    tmp = -1d-2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, tmp, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -17) stop 1
    call ddfree(ddx_data)
    ! Check correct eps
    tmp = zero
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, tmp, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    tmp = pt5
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, tmp, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    tmp = one
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, tmp, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    tmp = dble(1000)
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, tmp, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    ! Check incorrect eps
    tmp = -1d-2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, tmp, kappa, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -18) stop 1
    call ddfree(ddx_data)
    ! Check correct kappa
    tmp = 1d-2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, tmp, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    ! Check incorrect kappa
    tmp = -1d-2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, tmp, itersolver, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -19) stop 1
    call ddfree(ddx_data)
    ! Check correct itersolver
    i = 1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, i, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    ! Check incorrect itersolver
    i = 0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, i, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -20) stop 1
    call ddfree(ddx_data)
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, i, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -20) stop 1
    call ddfree(ddx_data)
    i = 2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, i, tol, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -20) stop 1
    call ddfree(ddx_data)
    ! Check correct tol
    tmp = 1d-14
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tmp, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    tmp = 1d-7
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tmp, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    tmp = one
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tmp, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    ! Check incorrect tol
    tmp = 9d-15
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tmp, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -21) stop 1
    call ddfree(ddx_data)
    tmp = zero
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tmp, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -21) stop 1
    call ddfree(ddx_data)
    tmp = -1d-1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tmp, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -21) stop 1
    call ddfree(ddx_data)
    tmp = 1.1d0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tmp, &
        & maxiter, ndiis, nproc, ddx_data, info)
    if (info .ne. -21) stop 1
    call ddfree(ddx_data)
    ! Check correct maxiter
    i = 1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & i, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    i = 1000000
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & i, ndiis, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    ! Check incorrect maxiter
    i = 0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & i, ndiis, nproc, ddx_data, info)
    if (info .ne. -22) stop 1
    call ddfree(ddx_data)
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & i, ndiis, nproc, ddx_data, info)
    if (info .ne. -22) stop 1
    call ddfree(ddx_data)
    ! Check correct ndiis
    i = 0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, i, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    i = 1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, i, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    i = 1000
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, i, nproc, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    ! Check incorrect ndiis
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, i, nproc, ddx_data, info)
    if (info .ne. -23) stop 1
    call ddfree(ddx_data)
    ! Check correct nproc
    i = 1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, i, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    i = 100000
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, i, ddx_data, info)
    if (info .ne. 0) stop 1
    call ddfree(ddx_data)
    ! Check incorrect nproc
    i = 0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, i, ddx_data, info)
    if (info .ne. -24) stop 1
    call ddfree(ddx_data)
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, &
        & maxiter, ndiis, i, ddx_data, info)
    if (info .ne. -24) stop 1
    call ddfree(ddx_data)
end subroutine check_ddinit_args

! Check P2M and M2P for spherical harmonics
! List of explicitly checked functions:
!       fmm_p2m
!       fmm_m2p
! List of implicitly checked functions:
!       ylmscale
!       ylmbas
!       trgev
!       polleg
! TODO: check adjoint routines
subroutine check_p2m_m2p(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasis, i, j, istatus
    real(dp) :: vscales((p+1)**2), vscales_rel((p+1)**2), v4pi2lp1(p+1)
    real(dp) :: src0(3), src(3), dst(3), csph(3, nsph), rsph(nsph), v0, v(nsph)
    real(dp), dimension((p+1)**2, nsph) :: coef, coef2
    real(dp) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    real(dp) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(dp) :: full_mat((p+1)**2, (p+1)**2)
    real(dp) :: z, r1(3, 3), dst1(3), full_norm, diff_norm
    logical :: ok
    real(dp), external :: dnrm2
    !! Scale inputs
    src0 = gsrc0
    src = alpha * gsrc
    dst = alpha * gdst(:, 1)
    csph = alpha * gcsph
    rsph = abs(alpha) * grsph
    ! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    !! Check P2M parameters q and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES12.4E3)") "Check P2M q and beta params for p=", &
            & p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    ! Check beta=zero
    coef(:, 1) = zero
    call fmm_p2m(src, one, alpha, p, vscales, zero, coef(:, 1))
    coef2(:, 1) = coef(:, 1)
    call fmm_p2m(src, one, alpha, p, vscales, zero, coef2(:, 1))
    v(1) = dnrm2(nbasis, coef(:, 1)-coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1)
    if (iprint .gt. 0) then
        write(*, *) "P2M beta=zero param check: ", v(1) .eq. zero
    end if
    if (v(1) .ne. zero) stop 1
    coef2(:, 1) = coef(:, 1)
    ! Check non-zero beta
    call fmm_p2m(src, one, alpha, p, vscales, one, coef2(:, 1))
    v(1) = dnrm2(nbasis, two*coef(:, 1)-coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1) / two
    if (iprint .gt. 0) then
        write(*, *) "P2M q and beta params check:", v(1) .eq. zero
    end if
    if (v(1) .ne. zero) stop 1
    call fmm_p2m(src, -pt5, alpha, p, vscales, -two, coef2(:, 1))
    v(1) = dnrm2(nbasis, 4.5d0*coef(:, 1)+coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1) / 4.5d0
    if (iprint .gt. 0) then
        write(*, *) "P2M q and beta params check:", v(1) .eq. zero
    end if
    if (v(1) .ne. zero) stop 1
    !! Check M2P parameters alpha and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,A,I0,A,ES12.4E3)") "Check M2P alpha and beta params ", &
            & "for p=", p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    ! Check beta=zero
    v(1) = zero
    call fmm_m2p(dst, alpha, p, vscales_rel, one, coef(:, 1), zero, v(1))
    v0 = v(1)
    call fmm_m2p(dst, alpha, p, vscales_rel, one, coef(:, 1), zero, v(1))
    if (iprint .gt. 0) then
        write(*, *) "M2P beta=zero param check:", v0 .eq. v(1)
    end if
    if (v0 .ne. v(1)) stop 1
    ! Check alpha=zero
    call fmm_m2p(dst, alpha, p, vscales_rel, zero, coef(:, 1), -one, v(1))
    if (iprint .gt. 0) then
        write(*, *) "M2P alpha=zero param check:", v0 .eq. -v(1)
    end if
    if (v0 .ne. -v(1)) stop 1
    ! Check non-zero alpha and beta
    v(1) = v0
    call fmm_m2p(dst, alpha, p, vscales_rel, one, coef(:, 1), one, v(1))
    ok = abs(v(1)/v0-two) .le. 10d0*epsilon(zero)
    if (iprint .gt. 0) then
        write(*, *) "M2P alpha and beta params check:", ok
    end if
    if (.not. ok) stop 1
    call fmm_m2p(dst, alpha, p, vscales_rel, -two, coef(:, 1), pt5, v(1))
    ok = abs(v(1)/v0+one) .le. 10d0*epsilon(zero)
    if (iprint .gt. 0) then
        write(*, *) "M2P alpha and beta parameter check:", ok
    end if
    if (.not. ok) stop 1
    !! Check M2P with a location of P equal to the center of M
    call fmm_m2p(src0, rsph(1), p, vscales_rel, one, coef(:, 1), zero, v0)
    if (iprint .gt. 0) then
        write(*, *) "Check M2P with center of harmonics=particle location"
        write(*, *) "================================="
        write(*, *) "Required result must be 0"
        write(*, *) "Got", v0
    end if
    if (v0 .ne. zero) stop 1
    !! Check P2M+M2P with source particle and multipole harmonics are centered
    !! at the origin
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES12.4E3)") "Check P2M and M2P for p=", p, &
            & " alpha=", alpha
        write(*, "(A)") "================================"
        write(*, *) "threshold=" , threshold
    end if
    call fmm_p2m(src0, one, alpha, p, vscales, zero, coef(:, 1))
    call fmm_m2p(dst, alpha, p, vscales_rel, one, coef(:, 1), zero, v(1))
    v0 = one / dnrm2(3, dst, 1)
    v(1) = abs((v(1)-v0) / v0)
    if (iprint .gt. 0) then
        write(*, "(A,A,ES24.16E3)") "|| P2M(0) + M(0)2P - P2P ||  /  ", &
            & "|| P2P || =", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    call fmm_p2m(src0, pt5, alpha, p, vscales, pt5, coef(:, 1))
    call fmm_m2p(dst, alpha, p, vscales_rel, one, coef(:, 1), zero, v(1))
    v(1) = abs((v(1)-v0) / v0)
    if (iprint .gt. 0) then
        write(*, "(A,A,ES24.16E3)") "|| P2M(0) + M(0)2P - P2P ||  /  ", &
            & "|| P2P || =", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    !! Check P2M+M2P for predefined spheres
    ! Compute multipole coefficients for given spheres from source particle
    ! Get potential, spawned by each sphere with its multipole expansion
    v0 = one / dnrm2(3, src-dst, 1)
    do i = 1, nsph
        call fmm_p2m(src-csph(:, i), one, rsph(i), p, vscales, zero, &
            & coef(:, i))
        call fmm_m2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v(i))
        ! Finally check p2m+m2p
        v(i) = abs((v(i)-v0) / v0)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(", i, ") + M(", i, &
                & ")2P - P2P ||  /  || P2P || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check adjoint M2P fully
    do i = 1, nsph
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_m2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
                & coef(:, i), zero, coef2(j, i))
        end do
        call fmm_m2p_adj(dst-csph(:, i), one, rsph(i), p, vscales_rel, zero, &
            & coef(:, i))
        !!write(*, *) coef(:, i) - coef2(:, i)
        v(i) = dnrm2((p+1)**2, coef(:, i)-coef2(:, i), 1) / &
            & dnrm2((p+1)**2, coef2(:, i), 1)
        write(*, *) "adj=", v(i)
        !write(*, *) maxval(abs(coef(:, i)-coef2(:, i)) / coef2(:, i))
        if (v(i) .gt. 1d-15) stop 1
        call fmm_m2p_adj(dst-csph(:, i), -one, rsph(i), p, vscales_rel, two, &
            & coef(:, i))
        v(i) = dnrm2((p+1)**2, coef(:, i)-coef2(:, i), 1) / &
            & dnrm2((p+1)**2, coef2(:, i), 1)
        write(*, *) "adj(beta!=0)=", v(i)
        if (v(i) .gt. 1d-15) stop 1
    end do
end subroutine check_p2m_m2p

! Check P2L and L2P for spherical harmonics
! List of explicitly checked functions:
!       fmm_p2l
!       fmm_l2p
subroutine check_p2l_l2p(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasis, i, j
    real(dp) :: delta = 10d0 * epsilon(one)
    real(dp) :: vscales((p+1)**2), vscales_rel((p+1)**2), v4pi2lp1(p+1)
    real(dp) :: dst0(3), dst(3), src(3), csph(3, nsph), rsph(nsph), v0, v(nsph)
    real(dp), dimension((p+1)**2, nsph) :: coef, coef2
    logical :: ok
    real(dp), external :: dnrm2
    ! Scale inputs
    dst0 = gsrc0
    dst = alpha * gsrc
    src = alpha * gdst(:, 1)
    csph = alpha * gcsph
    rsph = abs(alpha) * grsph
    ! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    !! Check P2L parameters q and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES12.4E3)") "Check P2L q and beta params for p=", &
            & p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    ! Check beta=zero
    coef(:, 1) = zero
    call fmm_p2l(src, one, alpha, p, vscales, zero, coef(:, 1))
    coef2(:, 1) = coef(:, 1)
    call fmm_p2l(src, one, alpha, p, vscales, zero, coef2(:, 1))
    v(1) = dnrm2(nbasis, coef(:, 1)-coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1)
    ok = v(1) .lt. delta
    if (iprint .gt. 0) then
        write(*, *) "P2L beta=zero param check: ", ok
    end if
    if (.not. ok) stop 1
    ! Check non-zero beta
    coef2(:, 1) = coef(:, 1)
    call fmm_p2l(src, one, alpha, p, vscales, one, coef2(:, 1))
    v(1) = dnrm2(nbasis, two*coef(:, 1)-coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1) / two
    ok = v(1) .lt. delta
    if (iprint .gt. 0) then
        write(*, *) "P2L q and beta params check:", ok
    end if
    if (.not. ok) stop 1
    call fmm_p2l(src, -pt5, alpha, p, vscales, -two, coef2(:, 1))
    v(1) = dnrm2(nbasis, 4.5d0*coef(:, 1)+coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1) / 4.5d0
    ok = v(1) .lt. delta
    if (iprint .gt. 0) then
        write(*, *) "P2L q and beta params check:", ok
    end if
    if (.not. ok) stop 1
    !! Check L2P parameters alpha and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,A,I0,A,ES12.4E3)") "Check L2P alpha and beta params ", &
            & "for p=", p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    ! Check beta=zero
    v(1) = zero
    call fmm_l2p(dst, alpha, p, vscales_rel, one, coef(:, 1), zero, v(1))
    v0 = v(1)
    call fmm_l2p(dst, alpha, p, vscales_rel, one, coef(:, 1), zero, v(1))
    ok = v0 .eq. v(1)
    if (iprint .gt. 0) then
        write(*, *) "L2P beta=zero param check:", ok
    end if
    if (.not. ok) stop 1
    ! Check alpha=zero
    call fmm_l2p(dst, alpha, p, vscales_rel, zero, coef(:, 1), -one, v(1))
    ok = v0 .eq. -v(1)
    if (iprint .gt. 0) then
        write(*, *) "L2P alpha=zero param check:", ok
    end if
    if (.not. ok) stop 1
    ! Check non-zero alpha and beta
    v(1) = v0
    call fmm_l2p(dst, alpha, p, vscales_rel, one, coef(:, 1), one, v(1))
    ok = abs(v(1)/v0-two) .le. delta
    if (iprint .gt. 0) then
        write(*, *) "L2P alpha and beta params check:", ok
    end if
    if (.not. ok) stop 1
    call fmm_l2p(dst, alpha, p, vscales_rel, -two, coef(:, 1), pt5, v(1))
    ok = abs(v(1)/v0+one) .le. delta
    if (iprint .gt. 0) then
        write(*, *) "L2P alpha and beta parameter check:", ok
    end if
    if (.not. ok) stop 1
    !! Check P2L with a location of P equal to the center of L
    if (iprint .gt. 0) then
        write(*, *) "Check P2L with center of harmonics=particle location"
        write(*, *) "================================="
        write(*, *) "Required norm of resulting vector must be 0"
    end if
    call fmm_p2l(dst0, one, rsph(1), p, vscales, zero, coef(:, 1))
    v0 = dnrm2((p+1)**2, coef(:, 1), 1)
    ok = v0 .eq. zero
    if (iprint .gt. 0) then
        write(*, *) "|| diff || =", v0
    end if
    if (.not. ok) stop 1
    coef(:, 1) = one
    call fmm_p2l(dst0, one, rsph(1), p, vscales, -one, coef(:, 1))
    v0 = dnrm2((p+1)**2, coef(:, 1)+one, 1)
    ok = v0 .eq. zero
    if (iprint .gt. 0) then
        write(*, *) "|| diff || =", v0
    end if
    if (.not. ok) stop 1
    !! Check P2L+L2P with source particle at the origin and local harmonics
    !! and target particle are at the same point
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES12.4E3)") "Check P2L and L2P for p=", p, &
            & " alpha=", alpha
        write(*, "(A)") "================================"
        write(*, *) "threshold=" , threshold
    end if
    call fmm_p2l(src, one, alpha, p, vscales, zero, coef(:, 1))
    call fmm_l2p(dst0, alpha, p, vscales_rel, one, coef(:, 1), zero, v(1))
    v0 = one / dnrm2(3, src, 1)
    v(1) = abs((v(1)-v0) / v0)
    if (iprint .gt. 0) then
        write(*, "(A,A,ES24.16E3)") "|| P2L(0) + L(0)2P - P2P ||  /  ", &
            & "|| P2P || =", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    call fmm_p2l(src, one, alpha, p, vscales, zero, coef(:, 1))
    call fmm_l2p(dst0, alpha, p, vscales_rel, one, coef(:, 1), zero, v(1))
    v(1) = abs((v(1)-v0) / v0)
    if (iprint .gt. 0) then
        write(*, "(A,A,ES24.16E3)") "|| P2L(0) + L(0)2P - P2P ||  /  ", &
            & "|| P2P || =", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    !! Check P2L+L2P for predefined spheres
    ! Get potential, spawned by each sphere with its multipole expansion
    v0 = one / dnrm2(3, src-dst, 1)
    ! Compute multipole coefficients for given spheres from source particle
    do i = 1, nsph
        call fmm_p2l(src-csph(:, i), one, rsph(i), p, vscales, zero, &
            & coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v(i))
        ! Finally check p2m+m2p
        v(i) = abs((v(i)-v0) / v0)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2L(", i, ") + L(", i, &
                & ")2P - P2P ||  /  || P2P || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check adjoint L2P fully
    do i = 1, nsph
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
                & coef(:, i), zero, coef2(j, i))
        end do
        call fmm_l2p_adj(dst-csph(:, i), one, rsph(i), p, vscales_rel, zero, &
            & coef(:, i))
        !!write(*, *) coef(:, i) - coef2(:, i)
        v(i) = dnrm2((p+1)**2, coef(:, i)-coef2(:, i), 1) / &
            & dnrm2((p+1)**2, coef2(:, i), 1)
        write(*, *) "adj=", v(i)
        !write(*, *) maxval(abs(coef(:, i)-coef2(:, i)) / coef2(:, i))
        if (v(i) .gt. 1d-15) stop 1
        call fmm_l2p_adj(dst-csph(:, i), -one, rsph(i), p, vscales_rel, two, &
            & coef(:, i))
        v(i) = dnrm2((p+1)**2, coef(:, i)-coef2(:, i), 1) / &
            & dnrm2((p+1)**2, coef2(:, i), 1)
        write(*, *) "adj(beta!=0)=", v(i)
        if (v(i) .gt. 1d-15) stop 1
    end do
end subroutine check_p2l_l2p

! Check M2M for spherical harmonics
! List of explicitly checked functions:
!       fmm_m2m_ztranslate
!       fmm_m2m_ztranslate_adj
!       fmm_m2m_ztranslate_get_mat
!       fmm_m2m_ztranslate_use_mat
!       fmm_m2m_ztranslate_use_mat_adj
!       fmm_m2m_rotation
!       fmm_m2m_rotation_adj
!       fmm_m2m_reflection_get_mat
!       fmm_m2m_reflection_use_mat
!       fmm_m2m_reflection_use_mat_adj
subroutine check_m2m(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasis, i, indi, j
    real(dp) :: vscales((p+1)**2), vfact(2*p+1), vscales_rel((p+1)**2), &
        & v4pi2lp1(p+1), vcnk((2*p+1)*(p+1)), m2l_ztranslate_coef(p+1, 1, 1), &
        & m2l_ztranslate_adj_coef(1, 1, p+1)
    real(dp) :: src0(3), src(3), dst(3), csph(3, nsph), rsph(nsph), v0, v(nsph)
    real(dp), dimension((p+1)**2, nsph) :: coef, coef2
    real(dp) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    real(dp) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(dp) :: full_mat((p+1)**2, (p+1)**2)
    real(dp) :: z, r1(3, 3), dst1(3), full_norm, diff_norm
    logical :: ok
    real(dp), external :: dnrm2
    !! Scale inputs
    src0 = gsrc0
    src = alpha * gsrc
    dst = alpha * gdst(:, 1)
    csph = alpha * gcsph
    rsph = alpha * grsph
    !! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    vfact(1) = one
    do i = 2, 2*p+1
        vfact(i) = vfact(i-1) * sqrt(dble(i-1))
    end do
    !! Get multipole coefficients of all spheres by P2M
    do i = 1, nsph
        call fmm_p2m(src-csph(:, i), one, rsph(i), p, vscales, zero, &
            & coef(:, i))
    end do
    ! Compute FMM-related constants
    call fmm_constants(p, p, 0, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    !! Check M2M OZ translation by fmm_m2m_ztranslate. Spheres 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_ztranslate"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vcnk, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vcnk, -one, coef(:, 1), three, coef2(:, i))
        v(i) = dnrm2(nbasis, pt5*coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    ! Check M2M with zero OZ translation. Sphere 6 is explicitly set for this
    call fmm_m2m_ztranslate(zero, rsph(1), rsph(6), p, vscales, vcnk, one, &
        & coef(:, 1), zero, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
    end if
    if (.not. ok) stop 1
    call fmm_m2m_ztranslate(zero, rsph(1), rsph(6), p, vscales, vcnk, -one, &
        & coef(:, 1), three, coef2(:, 6))
    v(6) = dnrm2(nbasis, pt5*coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check M2M OZ translation by precomputed matrices. Spheres 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_ztranslate_get_mat"
        write(*, "(A)") "Check fmm_m2m_ztranslate_use_mat"
        write(*, "(A)") "Check fmm_m2m_scale"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, coef(:, 1), &
            & zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, -one, coef(:, 1), &
            & three, coef2(:, i))
        v(i) = dnrm2(nbasis, pt5*coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    ! M2M with 0 shift is simply scaling
    call fmm_m2m_scale(rsph(1), rsph(6), p, one, coef(:, 1), zero, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
    end if
    if (.not. ok) stop 1
    call fmm_m2m_scale(rsph(1), rsph(6), p, -one, coef(:, 1), three, &
        & coef2(:, 6))
    v(6) = dnrm2(nbasis, pt5*coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check adjoint M2M OZ translation by precomputed matrices of direct M2M
    !! OZ translation. Spheres 2 and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_ztranslate_use_mat_adj"
        write(*, "(A)") "Check fmm_m2m_scale_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || M(1)2M(i) - [adj M(1)2M(i)]^T || ", &
            & "/ || M(1)2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        ! Generate entire matrix of a direct M2M OZ translation
        call fmm_m2m_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        do j = 1, (p+1)**2
            coef2(:, i) = zero
            coef2(j, i) = one
            call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, &
                & coef2(:, i), zero, full_mat(:, j))
        end do
        full_norm = dnrm2((p+1)**4, full_mat, 1)
        ! Subtract transpose of an adjoint M2M OZ translation matrix
        do j = 1, (p+1)**2
            coef2(:, i) = zero
            coef2(j, i) = one
            call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
                & coef2(:, i), two, full_mat(j, :))
        end do
        diff_norm = dnrm2((p+1)**4, full_mat, 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, one, &
            & coef(:, 1), zero, full_mat(:, 1))
        full_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
            & coef(:, 1), two, full_mat(:, 1))
        diff_norm = pt5*dnrm2((p+1)**2, full_mat(:, 1), 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    ! Check fmm_m2m_scale_adj. Sphere 6 is intended for this. Result is the
    ! same as with direct fmm_m2m_scale.
    call fmm_m2m_scale_adj(rsph(6), rsph(1), p, one, coef(:, 1), zero, &
        & coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
    end if
    if (.not. ok) stop 1
    call fmm_m2m_scale_adj(rsph(6), rsph(1), p, -one, coef(:, 1), three, &
        & coef2(:, 6))
    v(6) = dnrm2(nbasis, pt5*coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check adjoint M2M OZ translation without precomputed matrices. Spheres 2
    !! and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_ztranslate_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj M(1)2M(i)] - [adj M(1)2M(i)]", &
            & " || / || [adj M(1)2M(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_m2m_ztranslate_adj(csph(3, i), rsph(i), rsph(1), p, &
            & vscales, vcnk, one, coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, -one, &
            & coef(:, 1), one, coef2(:, i))
        diff_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, one, &
            & coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        call fmm_m2m_ztranslate_adj(csph(3, i), rsph(i), rsph(1), p, &
            & vscales, vcnk, -two, coef(:, 1), two, coef2(:, i))
        diff_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    !! Check that attempt to use M2M matrix with 0 shift raises NaNs
    call fmm_m2m_ztranslate_get_mat(zero, rsph(1), rsph(6), p, vscales, &
        & vfact, ztranslate_mat)
    call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, coef(:, 1), zero, &
        & coef2(:, 6))
    ! Adjoint with zero translation is the same as direct with zero translation
    call fmm_m2m_ztranslate_adj(zero, rsph(6), rsph(1), p, vscales, vcnk, &
        & one, coef(:, 1), zero, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
    end if
    if (.not. ok) stop 1
    call fmm_m2m_ztranslate_adj(zero, rsph(6), rsph(1), p, vscales, vcnk, &
        & -one, coef(:, 1), three, coef2(:, 6))
    v(6) = dnrm2(nbasis, pt5*coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check M2M by rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_rotation"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        call fmm_m2m_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vcnk, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vcnk, -one, coef(:, 1), three, coef2(:, i))
        v(i) = dnrm2(nbasis, pt5*coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check M2M by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_m2m_reflection_get_mat"
        write(*, "(A,I0)") "Check fmm_m2m_reflection_use_mat"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, one, coef(:, 1), &
            & zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, -one, coef(:, 1), &
            & two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint M2M by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_m2m_reflection_use_mat_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || M(1)2M(i) - [adj M(1)2M(i)]^T || ", &
            & "/ || M(1)2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        ! Generate entire matrix of a direct M2M translation
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        do j = 1, (p+1)**2
            coef2(:, i) = zero
            coef2(j, i) = one
            call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
                & rsph(i), p, transform_mat, ztranslate_mat, one, &
                & coef2(:, i), zero, full_mat(:, j))
        end do
        full_norm = dnrm2((p+1)**4, full_mat, 1)
        ! Subtract transpose of an adjoint M2M translation matrix
        do j = 1, (p+1)**2
            coef2(:, i) = zero
            coef2(j, i) = one
            call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
                & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, -two, &
                & coef2(:, i), two, full_mat(j, :))
        end do
        diff_norm = dnrm2((p+1)**4, full_mat, 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, one, coef(:, i), &
            & zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, -two, coef(:, i), &
            & two, coef2(:, i))
        diff_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, one, coef(:, i), &
            & zero, coef2(:, i))
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint M2M by rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_m2m_rotation_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj M(1)2M(i)] - [adj M(1)2M(i)]", &
            & " || / || [adj M(1)2M(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        ! Check with beta=zero
        call fmm_m2m_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vcnk, one, coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        ! Subtract reference value
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, -one, coef(:, 1), &
            & one, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        ! Get reference value
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, one, coef(:, 1), &
            & zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        ! Subtract value to be checked
        call fmm_m2m_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vcnk, -two, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
end subroutine check_m2m

! Check L2L for spherical harmonics
! List of explicitly checked functions:
!       fmm_l2l_ztranslate
!       fmm_l2l_ztranslate_adj
!       fmm_l2l_ztranslate_get_mat
!       fmm_l2l_ztranslate_use_mat
!       fmm_l2l_ztranslate_use_mat_adj
!       fmm_l2l_rotation
!       fmm_l2l_rotation_adj
!       fmm_l2l_reflection_get_mat
!       fmm_l2l_reflection_use_mat
!       fmm_l2l_reflection_use_mat_adj
subroutine check_l2l(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasis, i, j
    real(dp) :: vscales((p+1)**2), vfact(2*p+1), vscales_rel((p+1)**2), &
        & v4pi2lp1(p+1)
    real(dp) :: dst0(3), dst(3), src(3), csph(3, nsph), rsph(nsph), v0, v1
    real(dp), dimension((p+1)**2, nsph) :: coef
    real(dp) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    real(dp) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(dp) :: full_mat((p+1)**2, (p+1)**2)
    real(dp) :: z, r1(3, 3), dst1(3), full_norm, diff_norm
    logical :: ok
    real(dp), external :: dnrm2
    !! Scale inputs
    dst0 = gsrc0
    dst = alpha * gsrc
    src = alpha * gdst(:, 1)
    csph = alpha * gcsph
    rsph = abs(alpha) * grsph
    !! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    vfact(1) = one
    do i = 2, 2*p+1
        vfact(i) = vfact(i-1) * sqrt(dble(i-1))
    end do
    !! Get local coefficients of the main sphere (source of L2L) and
    !! corresponding potential
    call fmm_p2l(src-csph(:, 1), one, rsph(1), p, vscales, zero, coef(:, 1))
    call fmm_l2p(dst-csph(:, 1), rsph(1), p, vscales_rel, one, coef(:, 1), &
        & zero, v0)
    !! Check L2L OZ translation by fmm_l2l_ztranslate. Spheres 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_ztranslate"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_l2l_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vfact, one, coef(:, 1), zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_l2l_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vfact, -one, coef(:, 1), two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    ! Check L2L with zero OZ translation. Sphere 6 is explicitly set for this
    call fmm_l2l_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, one, &
        & coef(:, 1), zero, coef(:, 6))
    call fmm_l2p(dst-csph(:, 6), rsph(6), p, vscales_rel, one, coef(:, 6), &
        & zero, v1)
    v1 = abs(v1/v0 - one)
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
    end if
    if (.not. ok) stop 1
    call fmm_l2l_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, -one, &
        & coef(:, 1), two, coef(:, 6))
    call fmm_l2p(dst-csph(:, 6), rsph(6), p, vscales_rel, one, coef(:, 6), &
        & zero, v1)
    v1 = abs(v1/v0 - one)
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check L2L OZ translation by precomputed matrices. Spheres 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_ztranslate_get_mat"
        write(*, "(A)") "Check fmm_l2l_ztranslate_use_mat"
        write(*, "(A)") "Check fmm_l2l_scale"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_l2l_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, one, coef(:, 1), &
            & zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, -one, coef(:, 1), &
            & two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    ! L2L with 0 shift is simply scaling
    call fmm_l2l_scale(rsph(1), rsph(6), p, one, coef(:, 1), zero, coef(:, 6))
    call fmm_l2p(dst-csph(:, 6), rsph(6), p, vscales_rel, one, coef(:, 6), &
        & zero, v1)
    v1 = abs(v1/v0 - one)
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
    end if
    if (.not. ok) stop 1
    call fmm_l2l_scale(rsph(1), rsph(6), p, -one, coef(:, 1), two, &
        & coef(:, 6))
    call fmm_l2p(dst-csph(:, 6), rsph(6), p, vscales_rel, one, coef(:, 6), &
        & zero, v1)
    v1 = abs(v1/v0 - one)
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check adjoint L2L OZ translation by precomputed matrices of direct L2L
    !! OZ translation. Spheres 2 and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_ztranslate_use_mat_adj"
        write(*, "(A)") "Check fmm_l2l_scale_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) - [adj L(1)2L(i)]^T || ", &
            & "/ || L(1)2L(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        ! Generate entire matrix of a direct L2L OZ translation
        call fmm_l2l_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, one, &
                & coef(:, i), zero, full_mat(:, j))
        end do
        full_norm = dnrm2((p+1)**4, full_mat, 1)
        ! Subtract transpose of an adjoint L2L OZ translation matrix
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
                & coef(:, i), two, full_mat(j, :))
        end do
        diff_norm = dnrm2((p+1)**4, full_mat, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, one, &
            & coef(:, 1), zero, full_mat(:, 1))
        full_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
            & coef(:, 1), two, full_mat(:, 1))
        diff_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    ! Check fmm_l2l_scale_adj. Sphere 6 is intended for this. Result is the
    ! same as with direct fmm_l2l_scale.
    call fmm_l2l_scale_adj(rsph(6), rsph(1), p, one, coef(:, 1), zero, &
        & coef(:, 6))
    full_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    call fmm_l2l_scale(rsph(1), rsph(6), p, -one, coef(:, 1), one, coef(:, 6))
    diff_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    v1 = diff_norm / full_norm
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
    end if
    call fmm_l2l_scale(rsph(1), rsph(6), p, one, coef(:, 1), zero, coef(:, 6))
    full_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    call fmm_l2l_scale_adj(rsph(6), rsph(1), p, one, coef(:, 1), -one, &
        & coef(:, 6))
    diff_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    v1 = diff_norm / full_norm
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check adjoint L2L OZ translation without precomputed matrices. Spheres 2
    !! and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_ztranslate_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj L(1)2L(i)] - [adj L(1)2L(i)]", &
            & " || / || [adj L(1)2L(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_l2l_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_l2l_ztranslate_adj(csph(3, i), rsph(i), rsph(1), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        full_norm = dnrm2((p+1)**2, coef(:, i), 1)
        call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
            & coef(:, 1), two, coef(:, i))
        diff_norm = dnrm2((p+1)**2, coef(:, i), 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, one, &
            & coef(:, 1), zero, coef(:, i))
        full_norm = dnrm2((p+1)**2, coef(:, i), 1)
        call fmm_l2l_ztranslate_adj(csph(3, i), rsph(i), rsph(1), p, &
            & vscales, vfact, -two, coef(:, 1), two, coef(:, i))
        diff_norm = dnrm2((p+1)**2, coef(:, i), 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    !! Check that attempt to use L2L matrix with 0 shift raises NaNs
    call fmm_l2l_ztranslate_get_mat(zero, rsph(1), rsph(6), p, vscales, &
        & vfact, ztranslate_mat)
    call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, one, coef(:, 1), zero, &
        & coef(:, 6))
    ! Adjoint with zero translation is the same as direct with zero translation
    call fmm_l2l_ztranslate_adj(zero, rsph(6), rsph(1), p, vscales, vfact, &
        & one, coef(:, 1), zero, coef(:, 6))
    full_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    call fmm_l2l_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, &
        & one, coef(:, 1), -one, coef(:, 6))
    diff_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    v1 = diff_norm / full_norm
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
    end if
    if (.not. ok) stop 1
    call fmm_l2l_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, &
        & one, coef(:, 1), zero, coef(:, 6))
    full_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    call fmm_l2l_ztranslate_adj(zero, rsph(6), rsph(1), p, vscales, vfact, &
        & one, coef(:, 1), -one, coef(:, 6))
    diff_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    v1 = diff_norm / full_norm
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check L2L by rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_rotation"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        call fmm_l2l_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, -one, coef(:, 1), two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check L2L by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_l2l_reflection_get_mat"
        write(*, "(A,I0)") "Check fmm_l2l_reflection_use_mat"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        call fmm_l2l_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        call fmm_l2l_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, one, coef(:, 1), &
            & zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, -one, coef(:, 1), &
            & two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales_rel, one, &
            & coef(:, i), zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint L2L by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_l2l_reflection_use_mat_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) - [adj L(1)2L(i)]^T || ", &
            & "/ || L(1)2L(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        ! Generate entire matrix of a direct L2L translation
        call fmm_l2l_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_l2l_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
                & rsph(i), p, transform_mat, ztranslate_mat, one, coef(:, i), &
                & zero, full_mat(:, j))
        end do
        full_norm = dnrm2((p+1)**4, full_mat, 1)
        ! Subtract transpose of an adjoint L2L translation matrix
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
                & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, -two, &
                & coef(:, i), two, full_mat(j, :))
        end do
        diff_norm = dnrm2((p+1)**4, full_mat, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, one, coef(:, 1), &
            & zero, full_mat(:, 1))
        full_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, -two, coef(:, 1), &
            & two, full_mat(:, 1))
        diff_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint L2L by rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_l2l_rotation_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj L(1)2L(i)] - [adj L(1)2L(i)]", &
            & " || / || [adj L(1)2L(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        ! Get reference value
        call fmm_l2l_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, one, coef(:, 1), &
            & zero, coef(:, i))
        full_norm = dnrm2((p+1)**2, coef(:, i), 1)
        ! Subtract value to be checked
        call fmm_l2l_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef(:, i))
        v1 = dnrm2((p+1)**2, coef(:, i), 1) / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        ! Check with beta=zero
        call fmm_l2l_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, -one, coef(:, 1), &
            & one, coef(:, i))
        ! Subtract reference value
        v1 = dnrm2((p+1)**2, coef(:, i), 1) / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
end subroutine check_l2l

subroutine check_m2l(pm, pl, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: pm, pl, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasism, nbasisl, i, j, k, p, indi
    real(dp) :: vscales((pm+pl+1)**2), vfact(2*(pm+pl)+1), dst_csph(3, 2), &
        & dst_rsph(2), vscales_rel((pm+pl+1)**2), v4pi2lp1(pm+pl+1), &
        & vcnk((2*(pm+pl)+1)*(pm+pl+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1)
    real(dp) :: src0(3), src(3), dst(3, 2), csph(3, nsph), rsph(nsph), v0(2), &
        & v1
    real(dp) :: coefm((pm+1)**2, nsph), coef2m((pm+1)**2), coefl((pl+1)**2)
    real(dp) :: ztranslate_mat((min(pm,pl)+1)*(min(pm,pl)+2)* &
        & (3*max(pm,pl)+3-min(pm,pl))/6)
    real(dp) :: transform_mat((max(pm,pl)+1)*(2*max(pm,pl)+1)* &
        & (2*max(pm,pl)+3)/3)
    real(dp) :: full_mat((pl+1)**2, (pm+1)**2)
    real(dp) :: z, r1(3, 3), dst1(3), full_norm, diff_norm
    logical :: ok
    real(dp), external :: dnrm2
    !! Scale inputs
    src0 = gsrc0
    src = alpha * gsrc
    dst = alpha * gdst
    csph = alpha * gcsph
    rsph = abs(alpha) * grsph
    dst_csph = alpha * gdst_csph
    dst_rsph = abs(alpha) * gdst_rsph
    !! Preliminaries
    nbasism = (pm+1)**2
    nbasisl = (pl+1)**2
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(pm+pl, vscales, v4pi2lp1, vscales_rel)
    vfact(1) = one
    do i = 2, 2*(pm+pl)+1
        vfact(i) = vfact(i-1) * sqrt(dble(i-1))
    end do
    !! Get multipole coefficients of all spheres by P2M
    do i = 1, nsph
        call fmm_p2m(src-csph(:, i), one, rsph(i), pm, vscales, zero, &
            & coefm(:, i))
    end do
    ! Compute FMM-related constants
    call fmm_constants(pm+pl, pm, pl, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    !! Get reference value of potentials
    v0(1) = one / dnrm2(3, src-dst(:, 1), 1)
    v0(2) = one / dnrm2(3, src-dst(:, 2), 1)
    !! Check M2L OZ translation by fmm_m2l_ztranslate. Spheres 1, 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2l_ztranslate"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(i) + M(i)2L(1) + L(1)2P - P2P", &
            & " || / || P2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, 3
        call fmm_m2l_ztranslate(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, m2l_ztranslate_coef, one, &
            & coefm(:, i), zero, coefl)
        call fmm_l2p(dst(:, 1)-dst_csph(:, 1), dst_rsph(1), pl, vscales_rel, &
            & one, coefl, zero, v1)
        v1 = abs(v1/v0(1) - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_m2l_ztranslate(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, m2l_ztranslate_coef, -one, &
            & coefm(:, i), two, coefl)
        call fmm_l2p(dst(:, 1)-dst_csph(:, 1), dst_rsph(1), pl, vscales_rel, &
            & one, coefl, zero, v1)
        v1 = abs(v1/v0(1) - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check M2L OZ translation by precomputed matrices. Spheres 1, 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2l_ztranslate_get_mat"
        write(*, "(A)") "Check fmm_m2l_ztranslate_use_mat"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || M(i)2L(1) - M(i)2L(1) || / ", &
            & "|| M(i)2L(1) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, 3
        call fmm_m2l_ztranslate_get_mat(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, ztranslate_mat)
        call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, &
            & coefm(:, i), zero, coefl)
        full_norm = dnrm2((pl+1)**2, coefl, 1)
        call fmm_m2l_ztranslate(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, m2l_ztranslate_coef, -one, &
            & coefm(:, i), one, coefl)
        diff_norm = dnrm2((pl+1)**2, coefl, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_m2l_ztranslate(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, m2l_ztranslate_coef, one, &
            & coefm(:, i), zero, coefl)
        full_norm = dnrm2((pl+1)**2, coefl, 1)
        call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, -two, &
            & coefm(:, i), two, coefl)
        diff_norm = dnrm2((pl+1)**2, coefl, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Since centers of M and L harmonics are not allowed to be in the same
    !! place, there is no fmm_m2l_scale function.
    !! Check adjoint M2L OZ translation by precomputed matrices of direct L2L
    !! OZ translation. Spheres 1, 2 and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2l_ztranslate_use_mat_adj"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || M(i)2L(1) - [adj M(i)2L(1)]^T || ", &
            & "/ || M(i)2L(1) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, 3
        ! Generate entire matrix of a direct M2L OZ translation
        call fmm_m2l_ztranslate_get_mat(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, ztranslate_mat)
        do j = 1, (pm+1)**2
            coef2m(:) = zero
            coef2m(j) = one
            call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, &
                & coef2m, zero, coefl)
            full_mat(:, j) = coefl
        end do
        full_norm = dnrm2(((pm+1)*(pl+1))**2, full_mat, 1)
        ! Subtract transpose of an adjoint M2L OZ translation matrix
        do j = 1, (pl+1)**2
            coefl(:) = zero
            coefl(j) = one
            call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, -one, &
                & coefl, one, full_mat(j, :))
        end do
        diff_norm = dnrm2(((pm+1)*(pl+1))**2, full_mat, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, &
            & coefm(:, i), zero, coefl)
        call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, one, &
            & coefl, zero, coef2m)
        full_norm = dnrm2((pm+1)**2, coef2m, 1)
        call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, -two, &
            & coefl, two, coef2m)
        diff_norm = dnrm2((pm+1)**2, coef2m, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint M2L OZ translation without precomputed matrices. Spheres
    !! 1, 2 and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2l_ztranslate_adj"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj M(i)2L(1)] - [adj M(i)2L(1)]", &
            & " || / || [adj M(i)2L(1)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, 3
        call fmm_m2l_ztranslate_get_mat(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, ztranslate_mat)
        call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, &
            & coefm(:, i), zero, coefl)
        call fmm_m2l_ztranslate_adj(dst_csph(3, 1)-csph(3, i), dst_rsph(1), &
            & rsph(i), pl, pm, vscales, m2l_ztranslate_adj_coef, one, coefl, &
            & zero, coef2m)
        full_norm = dnrm2((pm+1)**2, coef2m, 1)
        call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, -two, &
            & coefl, two, coef2m)
        diff_norm = dnrm2((pm+1)**2, coef2m, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, one, &
            & coefl, zero, coef2m)
        full_norm = dnrm2((pm+1)**2, coef2m, 1)
        call fmm_m2l_ztranslate_adj(dst_csph(3, 1)-csph(3, i), dst_rsph(1), &
            & rsph(i), pl, pm, vscales, m2l_ztranslate_adj_coef, -two, coefl, &
            & two, coef2m)
        diff_norm = dnrm2((pm+1)**2, coef2m, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check M2L by rotation
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_rotation"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || M(i)2L(", k, &
                & ") - M(i)2L(", k, ") || / || M(i)2L(", k, ") ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            call fmm_m2l_rotation(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, m2l_ztranslate_coef, one, &
                & coefm(:, i), zero, coefl)
            call fmm_l2p(dst(:, k)-dst_csph(:, k), dst_rsph(k), pl, &
                & vscales_rel, one, coefl, zero, v1)
            v1 = abs(v1/v0(k) - one)
            ok = v1 .le. threshold
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            call fmm_m2l_rotation(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, m2l_ztranslate_coef, -one, &
                & coefm(:, i), three, coefl)
            call fmm_l2p(dst(:, k)-dst_csph(:, k), dst_rsph(k), pl, &
                & vscales_rel, pt5, coefl, zero, v1)
            v1 = abs(v1/v0(k) - one)
            ok = v1 .le. threshold
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check M2L by reflection with help of matrices
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_reflection_get_mat"
            write(*, "(A)") "Check fmm_m2l_reflection_use_mat"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || M(i)2L(", k, &
                & ") - M(i)2L(", k, ") || / || M(i)2L(", k, ") ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            call fmm_m2l_reflection_get_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, vscales, vfact, &
                & transform_mat, ztranslate_mat)
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, transform_mat, &
                & ztranslate_mat, one, coefm(:, i), zero, coefl)
            full_norm = dnrm2((pl+1)**2, coefl, 1)
            call fmm_m2l_rotation(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, m2l_ztranslate_coef, -one, &
                & coefm(:, i), one, coefl)
            diff_norm = dnrm2((pl+1)**2, coefl, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            call fmm_m2l_rotation(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, m2l_ztranslate_coef, one, &
                & coefm(:, i), zero, coefl)
            full_norm = dnrm2((pl+1)**2, coefl, 1)
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(K), pm, pl, transform_mat, &
                & ztranslate_mat, -two, coefm(:, i), two, coefl)
            diff_norm = dnrm2((pl+1)**2, coefl, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check adjoint M2L by reflection with help of matrices
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_reflection_use_mat_adj"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || M(i)2L(", k, &
                & ") - [adj M(i)2L(", k, ")]^T || / || M(i)2L(", k, ") ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            ! Generate entire matrix of a direct M2L translation
            call fmm_m2l_reflection_get_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, vscales, vfact, &
                & transform_mat, ztranslate_mat)
            do j = 1, (pm+1)**2
                coef2m = zero
                coef2m(j) = one
                call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                    & rsph(i), dst_rsph(k), pm, pl, transform_mat, &
                    & ztranslate_mat, one, coef2m, zero, full_mat(:, j))
            end do
            full_norm = dnrm2(((pm+1)*(pl+1))**2, full_mat, 1)
            ! Subtract transpose of an adjoint L2L translation matrix
            do j = 1, (pl+1)**2
                coefl = zero
                coefl(j) = one
                call fmm_m2l_reflection_use_mat_adj( &
                    & dst_csph(:, k)-csph(:, i), dst_rsph(k), rsph(i), pl, &
                    & pm, transform_mat, ztranslate_mat, -two, coefl, two, &
                    & full_mat(j, :))
            end do
            diff_norm = dnrm2(((pm+1)*(pl+1))**2, full_mat, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(K), pm, pl, transform_mat, &
                & ztranslate_mat, one, coefm(:, i), zero, coefl)
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, one, coefl, zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, -two, coefl, two, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check adjoint M2L by rotation
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_rotation_adj"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || [adj M(i)2L(", k, &
                & ")] - [adj M(i)2L(", k, ")] || / || [adj M(i)2L(", k, ")] ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            ! Check with help of fmm_m2l_reflection_use_mat_adj
            call fmm_m2l_reflection_get_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, vscales, vfact, &
                & transform_mat, ztranslate_mat)
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, transform_mat, &
                & ztranslate_mat, one, coefm(:, i), zero, coefl)
            ! Get reference value
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, one, coefl, zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            ! Subtract value to be checked
            call fmm_m2l_rotation_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, vscales, &
                & m2l_ztranslate_adj_coef, -two, coefl, two, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            ! Check with beta=zero
            call fmm_m2l_rotation_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, vscales, &
                & m2l_ztranslate_adj_coef, one, coefl, zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            ! Subtract reference value
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, -one, coefl, one, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
end subroutine

subroutine check_tree_rib(alpha)
    real(dp), intent(in) :: alpha
    integer, parameter :: nsph = 10
    real(dp) :: csph(3, nsph), rsph(nsph), csph2(3, nsph), rsph2(nsph), &
        & cnode(3, 2*nsph-1), rnode(2*nsph-1)
    integer :: order(nsph), i, reorder(nsph), cluster(2, 2*nsph-1), &
        & children(2, 2*nsph-1), parent(2*nsph-1), snode(nsph)
    ! Scale inputs
    csph(:, 1) = alpha * (/1d0, 1d0, 1d0/)
    csph(:, 2) = alpha * (/2d0, 2d0, 2d0/)
    csph(:, 3) = alpha * (/1d0, 1d0, 3d0/)
    csph(:, 4) = alpha * (/2d0, 2d0, 4d0/)
    csph(:, 5) = alpha * (/1d0, 1d0, 5d0/)
    csph(:, 6) = alpha * (/2d0, 2d0, 6d0/)
    csph(:, 7) = alpha * (/1d0, 1d0, 7d0/)
    csph(:, 8) = alpha * (/2d0, 2d0, 8d0/)
    csph(:, 9) = alpha * (/1d0, 1d0, 9d0/)
    csph(:, 10) = alpha * (/2d0, 2d0, 10d0/)
    rsph = abs(alpha) * (/1d-1, 2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, &
        & 9d-1, 1d0/)
    ! Reorder inputs
    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
    do i = 1, nsph
        csph2(:, i) = csph(:, order(i))
        rsph2(i) = rsph(order(i))
    end do
    ! Build a recursive inertial binary tree
    call tree_rib_build(nsph, csph2, rsph2, reorder, cluster, children, &
        & parent, cnode, rnode, snode)
end subroutine check_tree_rib

subroutine check_tree_m2m(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer, parameter :: nsph = 10
    integer :: i, indi, j, k, order(nsph), istatus
    real(dp) :: csph(3, nsph), rsph(nsph), src(3, nsph), csph2(3, nsph), &
        & rsph2(nsph), sph_m((p+1)**2, nsph), node_m((p+1)**2, 2*nsph-1), &
        & node_m2((p+1)**2, 2*nsph-1), vscales((p+1)**2), vfact(2*p+1), &
        & rel_err, sph_m2((p+1)**2, nsph), full_norm, diff_norm, &
        & node_m3((p+1)**2, 2*nsph-1)
    type(ddx_type) :: ddx_data
    logical :: ok
    real(dp), external :: dnrm2
    real(dp) :: full_mat((p+1)**2, 2*nsph-1, (p+1)**2, nsph)
    ! Scale inputs
    csph(:, 1) = alpha * (/1d0, 1d0, 1d0/)
    csph(:, 2) = alpha * (/2d0, 2d0, 2d0/)
    csph(:, 3) = alpha * (/1d0, 1d0, 3d0/)
    csph(:, 4) = alpha * (/2d0, 2d0, 4d0/)
    csph(:, 5) = alpha * (/1d0, 1d0, 5d0/)
    csph(:, 6) = alpha * (/2d0, 2d0, 6d0/)
    csph(:, 7) = alpha * (/1d0, 1d0, 7d0/)
    csph(:, 8) = alpha * (/2d0, 2d0, 8d0/)
    csph(:, 9) = alpha * (/1d0, 1d0, 9d0/)
    csph(:, 10) = alpha * (/2d0, 2d0, 10d0/)
    rsph = abs(alpha) * (/2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, &
        & 1d0, 1.1d0/)
    src(:, 1) = alpha * (/1.1d0, 0.9d0, 1d0/)
    src(:, 2) = alpha * (/2.1d0, 1.9d0, 2.1d0/)
    src(:, 3) = alpha * (/1.1d0, 1.1d0, 3.1d0/)
    src(:, 4) = alpha * (/1.9d0, 1.9d0, 4.1d0/)
    src(:, 5) = alpha * (/1.1d0, 0.9d0, 5d0/)
    src(:, 6) = alpha * (/1.9d0, 2.1d0, 5.9d0/)
    src(:, 7) = alpha * (/1.1d0, 0.9d0, 7d0/)
    src(:, 8) = alpha * (/2.1d0, 2.1d0, 8d0/)
    src(:, 9) = alpha * (/1.1d0, 0.9d0, 9d0/)
    src(:, 10) = alpha * (/2.1d0, 1.9d0, 9.9d0/)
    ! Reorder inputs
    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
    do i = 1, nsph
        csph2(:, i) = csph(:, order(i))
        rsph2(i) = rsph(order(i))
    end do
    ddx_data % nclusters = 2*nsph-1
    ddx_data % pl = -1
    ddx_data % pm = p
    ddx_data % nsph = nsph
    ! Allocate space for a tree
    allocate(ddx_data % cluster(2, 2*nsph-1), ddx_data % children(2, 2*nsph-1), &
        & ddx_data % parent(2*nsph-1), ddx_data % cnode(3, 2*nsph-1), &
        & ddx_data % rnode(2*nsph-1), ddx_data % snode(nsph), &
        & ddx_data % order(nsph), ddx_data % vscales((p+1)**2), &
        & ddx_data % vfact(2*p+1), ddx_data % vscales_rel((p+1)**2), &
        & ddx_data % v4pi2lp1(p+1), ddx_data % vcnk((2*p+1)*(p+1)), &
        & ddx_data % m2l_ztranslate_coef(p+1, 0, 0), &
        & ddx_data % m2l_ztranslate_adj_coef(0, 0, p+1), stat=istatus)
    if(istatus .ne. 0) stop "ALLOC FAILED"
    ! Build a recursive inertial binary tree
    call tree_rib_build(nsph, csph2, rsph2, ddx_data % order, &
        & ddx_data % cluster, ddx_data % children, ddx_data % parent, &
        & ddx_data % cnode, ddx_data % rnode, ddx_data % snode)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, ddx_data % vscales, ddx_data % v4pi2lp1, &
        & ddx_data % vscales_rel)
    ddx_data % vfact(1) = one
    do i = 2, 2*p+1
        ddx_data % vfact(i) = ddx_data % vfact(i-1) * sqrt(dble(i-1))
    end do
    call fmm_constants(p, p, -1, ddx_data % vcnk, &
        & ddx_data % m2l_ztranslate_coef, &
        & ddx_data % m2l_ztranslate_adj_coef)
    ! Init input harmonics
    do i = 1, nsph
        call fmm_p2m(src(:, i)-csph(:, i), one, rsph(i), p, &
            & ddx_data % vscales, zero, sph_m(:, i))
    end do
    ! Get reference result of M2M operation
    do i = 1, 2*nsph-1
        node_m(:, i) = zero
        do j = ddx_data % cluster(1, i), ddx_data % cluster(2, i)
            k = ddx_data % order(j)
            call fmm_m2m_rotation(csph2(:, k)-ddx_data % cnode(:, i), &
                & rsph2(k), ddx_data % rnode(i), p, ddx_data % vscales, &
                & ddx_data % vcnk, one, sph_m(:, k), one, node_m(:, i))
        end do
    end do
    ! Check tree_m2m_rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_m2m_rotation"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree M2M] - [plain M2M] ||", &
            & " / || [plain M2M] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    node_m2 = one
    do i = 1, nsph
        node_m2(:, ddx_data % snode(i)) = sph_m(:, i)
    end do
    call tree_m2m_rotation(ddx_data, node_m2)
    rel_err = dnrm2(((p+1)**2)*(2*nsph-1), node_m-node_m2, 1) / &
        & dnrm2(((p+1)**2)*(2*nsph-1), node_m, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Allocate space for transfer matrices
    ddx_data % m2m_reflect_mat_size = (p+1)*(2*p+1)*(2*p+3)/3
    ddx_data % m2m_ztranslate_mat_size = (p+1)*(p+2)*(p+3)/6
    allocate(ddx_data % m2m_reflect_mat(ddx_data % m2m_reflect_mat_size, &
        & ddx_data % nclusters-1), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    allocate(ddx_data % m2m_ztranslate_mat(ddx_data % m2m_ztranslate_mat_size, &
        & ddx_data % nclusters-1), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    ! Compute transfer matrices
    call tree_m2m_reflection_get_mat(ddx_data)
    ! Check tree_m2m_reflection_use_mat
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_m2m_reflection_get_mat"
        write(*, "(A)") "Check tree_m2m_reflection_use_mat"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree M2M] - [plain M2M] ||", &
            & " / || [plain M2M] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    node_m2 = one
    do i = 1, nsph
        node_m2(:, ddx_data % snode(i)) = sph_m(:, i)
    end do
    call tree_m2m_reflection_use_mat(ddx_data, node_m2)
    rel_err = dnrm2(((p+1)**2)*(2*nsph-1), node_m-node_m2, 1) / &
        & dnrm2(((p+1)**2)*(2*nsph-1), node_m, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Check tree_m2m_reflection_use_mat_adj
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_m2m_reflection_use_mat_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree M2M] - [plain M2M] ||", &
            & " / || [plain M2M] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        do j = 1, (p+1)**2
            node_m2 = zero
            node_m2(j, ddx_data % snode(i)) = one
            call tree_m2m_reflection_use_mat(ddx_data, node_m2)
            full_mat(:, :, j, i) = node_m2
        end do
    end do
    full_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
    do i = 1, 2*nsph-1
        do j = 1, (p+1)**2
            node_m2 = zero
            node_m2(j, i) = one
            call tree_m2m_reflection_use_mat_adj(ddx_data, node_m2)
            do k = 1, nsph
                full_mat(j, i, :, k) = full_mat(j, i, :, k) - &
                    & node_m2(:, ddx_data % snode(k))
            end do
        end do
    end do
    diff_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
    rel_err = diff_norm / full_norm
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Check tree_m2m_rotation_adj
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_m2m_rotation_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree M2M] - [plain M2M] ||", &
            & " / || [plain M2M] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    node_m2 = node_m
    call tree_m2m_reflection_use_mat_adj(ddx_data, node_m2)
    node_m3 = node_m
    call tree_m2m_rotation_adj(ddx_data, node_m3)
    rel_err = dnrm2((2*nsph-1)*((p+1)**2), node_m2-node_m3, 1) / &
        dnrm2((2*nsph-1)*((p+1)**2), node_m2, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Deallocate tree
    call ddfree(ddx_data)
end subroutine check_tree_m2m

subroutine check_tree_l2l(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer, parameter :: nsph = 10
    integer :: i, j, k, order(nsph), istatus
    real(dp) :: csph(3, nsph), rsph(nsph), src(3, nsph), csph2(3, nsph), &
        & rsph2(nsph), sph_l((p+1)**2, nsph), node_l((p+1)**2, 2*nsph-1), &
        & node_l2((p+1)**2, 2*nsph-1), vscales((p+1)**2), vfact(2*p+1), &
        & rel_err, sph_l2((p+1)**2, nsph), full_norm, diff_norm, &
        & node_l3((p+1)**2, 2*nsph-1)
    type(ddx_type) :: ddx_data
    logical :: ok
    real(dp), external :: dnrm2
    real(dp) :: full_mat((p+1)**2, nsph, (p+1)**2, 2*nsph-1)
    ! Scale inputs
    csph(:, 1) = alpha * (/1d0, 1d0, 1d0/)
    csph(:, 2) = alpha * (/2d0, 2d0, 2d0/)
    csph(:, 3) = alpha * (/1d0, 1d0, 3d0/)
    csph(:, 4) = alpha * (/2d0, 2d0, 4d0/)
    csph(:, 5) = alpha * (/1d0, 1d0, 5d0/)
    csph(:, 6) = alpha * (/2d0, 2d0, 6d0/)
    csph(:, 7) = alpha * (/1d0, 1d0, 7d0/)
    csph(:, 8) = alpha * (/2d0, 2d0, 8d0/)
    csph(:, 9) = alpha * (/1d0, 1d0, 9d0/)
    csph(:, 10) = alpha * (/2d0, 2d0, 10d0/)
    rsph = abs(alpha) * (/2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, &
        & 1d0, 1.1d0/)
    src(:, 1) = alpha * (/1.1d0, 0.9d0, 1d0/)
    src(:, 2) = alpha * (/2.1d0, 1.9d0, 2.1d0/)
    src(:, 3) = alpha * (/1.1d0, 1.1d0, 3.1d0/)
    src(:, 4) = alpha * (/1.9d0, 1.9d0, 4.1d0/)
    src(:, 5) = alpha * (/1.1d0, 0.9d0, 5d0/)
    src(:, 6) = alpha * (/1.9d0, 2.1d0, 5.9d0/)
    src(:, 7) = alpha * (/1.1d0, 0.9d0, 7d0/)
    src(:, 8) = alpha * (/2.1d0, 2.1d0, 8d0/)
    src(:, 9) = alpha * (/1.1d0, 0.9d0, 9d0/)
    src(:, 10) = alpha * (/2.1d0, 1.9d0, 9.9d0/)
    ! Reorder inputs
    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
    do i = 1, nsph
        csph2(:, i) = csph(:, order(i))
        rsph2(i) = rsph(order(i))
    end do
    ddx_data % nclusters = 2*nsph-1
    ddx_data % pl = p
    ddx_data % pm = p
    ddx_data % nsph = nsph
    ! Allocate space for a tree
    allocate(ddx_data % cluster(2, 2*nsph-1), ddx_data % children(2, 2*nsph-1), &
        & ddx_data % parent(2*nsph-1), ddx_data % cnode(3, 2*nsph-1), &
        & ddx_data % rnode(2*nsph-1), ddx_data % snode(nsph), &
        & ddx_data % order(nsph), ddx_data % vscales((p+1)**2), &
        & ddx_data % vfact(2*p+1), ddx_data % vscales_rel((p+1)**2), &
        & ddx_data % v4pi2lp1(p+1), ddx_data % vcnk((2*p+1)*(p+1)), &
        & ddx_data % m2l_ztranslate_coef(0, p+1, p+1), &
        & ddx_data % m2l_ztranslate_adj_coef(p+1, p+1, 0), stat=istatus)
    ! Build a recursive inertial binary tree
    call tree_rib_build(nsph, csph2, rsph2, ddx_data % order, &
        & ddx_data % cluster, ddx_data % children, ddx_data % parent, &
        & ddx_data % cnode, ddx_data % rnode, ddx_data % snode)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, ddx_data % vscales, ddx_data % v4pi2lp1, &
        & ddx_data % vscales_rel)
    ddx_data % vfact(1) = one
    do i = 2, 2*p+1
        ddx_data % vfact(i) = ddx_data % vfact(i-1) * sqrt(dble(i-1))
    end do
    call fmm_constants(p, -1, p, ddx_data % vcnk, &
        & ddx_data % m2l_ztranslate_coef, &
        & ddx_data % m2l_ztranslate_adj_coef)
    ! Init input harmonics
    node_l = zero
    do i = 1, nsph
        call fmm_p2m(src(:, i)-csph(:, i), one, rsph(i), p, &
            & ddx_data % vscales, zero, node_l(:, ddx_data % snode(i)))
    end do
    call tree_m2m_rotation(ddx_data, node_l)
    ddx_data % pm = -1
    ! Get reference result of L2L operation
    sph_l = zero
    do i = 1, 2*nsph-1
        do j = ddx_data % cluster(1, i), ddx_data % cluster(2, i)
            k = ddx_data % order(j)
            call fmm_l2l_rotation(ddx_data % cnode(:, i)-csph2(:, k), &
                & ddx_data % rnode(i), rsph2(k), p, ddx_data % vscales, &
                & ddx_data % vfact, one, node_l(:, i), one, sph_l(:, k))
        end do
    end do
    ! Check tree_l2l_rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_l2l_rotation"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree L2L] - [plain L2L] ||", &
            & " / || [plain L2L] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    node_l2 = node_l
    call tree_l2l_rotation(ddx_data, node_l2)
    do i = 1, nsph
        sph_l2(:, i) = node_l2(:, ddx_data % snode(i))
    end do
    rel_err = dnrm2(((p+1)**2)*nsph, sph_l-sph_l2, 1) / &
        & dnrm2(((p+1)**2)*nsph, sph_l, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Allocate space for transfer matrices
    ddx_data % l2l_reflect_mat_size = (p+1)*(2*p+1)*(2*p+3)/3
    ddx_data % l2l_ztranslate_mat_size = (p+1)*(p+2)*(p+3)/6
    allocate(ddx_data % l2l_reflect_mat(ddx_data % l2l_reflect_mat_size, &
        & ddx_data % nclusters-1), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    allocate(ddx_data % l2l_ztranslate_mat(ddx_data % l2l_ztranslate_mat_size, &
        & ddx_data % nclusters-1), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    ! Compute transfer matrices
    call tree_l2l_reflection_get_mat(ddx_data)
    ! Check tree_l2l_reflection_use_mat
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_l2l_reflection_get_mat"
        write(*, "(A)") "Check tree_l2l_reflection_use_mat"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree L2L] - [plain L2L] ||", &
            & " / || [plain L2L] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    node_l2 = node_l
    call tree_l2l_reflection_use_mat(ddx_data, node_l2)
    do i = 1, nsph
        sph_l2(:, i) = node_l2(:, ddx_data % snode(i))
    end do
    rel_err = dnrm2(((p+1)**2)*nsph, sph_l-sph_l2, 1) / &
        & dnrm2(((p+1)**2)*nsph, sph_l, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Check tree_l2l_reflection_use_mat_adj
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_l2l_reflection_use_mat_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree L2L] - [plain L2L] ||", &
            & " / || [plain L2L] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, 2*nsph-1
        do j = 1, (p+1)**2
            node_l2 = zero
            node_l2(j, i) = one
            call tree_l2l_reflection_use_mat(ddx_data, node_l2)
            do k = 1, nsph
                full_mat(:, k, j, i) = node_l2(:, ddx_data % snode(k))
            end do
        end do
    end do
    full_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
    do i = 1, nsph
        do j = 1, (p+1)**2
            node_l2 = one
            do k = 1, nsph
                node_l2(:, ddx_data % snode(k)) = zero
            end do
            node_l2(j, ddx_data % snode(i)) = one
            call tree_l2l_reflection_use_mat_adj(ddx_data, node_l2)
            full_mat(j, i, :, :) = full_mat(j, i, :, :) - node_l2
        end do
    end do
    diff_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
    rel_err = diff_norm / full_norm
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Check tree_l2l_rotation_adj
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_l2l_rotation_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree L2L] - [plain L2L] ||", &
            & " / || [plain L2L] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    node_l2 = zero
    do i = 1, nsph
        node_l2(:, ddx_data % snode(i)) = node_l(:, ddx_data % snode(i))
    end do
    call tree_l2l_reflection_use_mat_adj(ddx_data, node_l2)
    node_l3 = zero
    do i = 1, nsph
        node_l3(:, ddx_data % snode(i)) = node_l(:, ddx_data % snode(i))
    end do
    call tree_l2l_rotation_adj(ddx_data, node_l3)
    rel_err = dnrm2((2*nsph-1)*((p+1)**2), node_l2-node_l3, 1) / &
        dnrm2((2*nsph-1)*((p+1)**2), node_l2, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Deallocate tree
    call ddfree(ddx_data)
end subroutine check_tree_l2l

subroutine check_tree_m2l(pm, pl, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: pm, pl, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer, parameter :: nsph = 10, lwork = 1000
    integer :: i, j, k, l, order(nsph), istatus, iwork, jwork, work(3, lwork)
    real(dp) :: csph(3, nsph), rsph(nsph), src(3, nsph), csph2(3, nsph), &
        & rsph2(nsph), node_m((pm+1)**2, 2*nsph-1), &
        & node_l((pl+1)**2, 2*nsph-1), vscales((pm+pl+1)**2), &
        & vfact(2*(pm+pl)+1), rel_err, full_norm, diff_norm, &
        & node_m2((pm+1)**2, 2*nsph-1), node_l2((pl+1)**2, 2*nsph-1)
    type(ddx_type) :: ddx_data
    logical :: ok
    real(dp), external :: dnrm2
    real(dp) :: full_mat((pl+1)**2, 2*nsph-1, (pm+1)**2, 2*nsph-1)
    real(dp) :: far_p2p(nsph), far_p2p2(nsph)
    ! Scale inputs
    csph(:, 1) = alpha * (/1d0, 1d0, 1d0/)
    csph(:, 2) = alpha * (/2d0, 2d0, 2d0/)
    csph(:, 3) = alpha * (/1d0, 1d0, 3d0/)
    csph(:, 4) = alpha * (/2d0, 2d0, 4d0/)
    csph(:, 5) = alpha * (/1d0, 1d0, 5d0/)
    csph(:, 6) = alpha * (/2d0, 2d0, 6d0/)
    csph(:, 7) = alpha * (/1d0, 1d0, 7d0/)
    csph(:, 8) = alpha * (/2d0, 2d0, 8d0/)
    csph(:, 9) = alpha * (/1d0, 1d0, 9d0/)
    csph(:, 10) = alpha * (/2d0, 2d0, 10d0/)
    rsph = abs(alpha) * (/2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, &
        & 1d0, 1.1d0/)
    src(:, 1) = alpha * (/1.1d0, 0.9d0, 1d0/)
    src(:, 2) = alpha * (/2.1d0, 1.9d0, 2.1d0/)
    src(:, 3) = alpha * (/1.1d0, 1.1d0, 3.1d0/)
    src(:, 4) = alpha * (/1.9d0, 1.9d0, 4.1d0/)
    src(:, 5) = alpha * (/1.1d0, 0.9d0, 5d0/)
    src(:, 6) = alpha * (/1.9d0, 2.1d0, 5.9d0/)
    src(:, 7) = alpha * (/1.1d0, 0.9d0, 7d0/)
    src(:, 8) = alpha * (/2.1d0, 2.1d0, 8d0/)
    src(:, 9) = alpha * (/1.1d0, 0.9d0, 9d0/)
    src(:, 10) = alpha * (/2.1d0, 1.9d0, 9.9d0/)
    ! Reorder inputs
    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
    do i = 1, nsph
        csph2(:, i) = csph(:, order(i))
        rsph2(i) = rsph(order(i))
    end do
    ddx_data % nclusters = 2*nsph-1
    ddx_data % pl = pl
    ddx_data % pm = pm
    ddx_data % nsph = nsph
    ! Allocate space for a tree
    allocate(ddx_data % cluster(2, 2*nsph-1), ddx_data % children(2, 2*nsph-1), &
        & ddx_data % parent(2*nsph-1), ddx_data % cnode(3, 2*nsph-1), &
        & ddx_data % rnode(2*nsph-1), ddx_data % snode(nsph), &
        & ddx_data % order(nsph), ddx_data % vscales((pm+pl+1)**2), &
        & ddx_data % vfact(2*(pm+pl)+1), ddx_data % vscales_rel((pm+pl+1)**2), &
        & ddx_data % v4pi2lp1(pm+pl+1), &
        & ddx_data % vcnk((2*(pm+pl)+1)*(pm+pl+1)), &
        & ddx_data % m2l_ztranslate_coef(pm+1, pl+1, pl+1), &
        & ddx_data % m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), &
        & stat=istatus)
    allocate(ddx_data % nfar(ddx_data % nclusters), &
        & ddx_data % nnear(ddx_data % nclusters))
    ! Build a recursive inertial binary tree
    call tree_rib_build(nsph, csph2, rsph2, ddx_data % order, &
        & ddx_data % cluster, ddx_data % children, ddx_data % parent, &
        & ddx_data % cnode, ddx_data % rnode, ddx_data % snode)
    ! Get list of neighbours for M2L operation
    iwork = 0
    call tree_get_farnear_work(ddx_data % nclusters, ddx_data % children, &
        & ddx_data % cnode, ddx_data % rnode, lwork, iwork, jwork, work, &
        & ddx_data % nnfar, ddx_data % nfar, ddx_data % nnnear, ddx_data % nnear)
    if (iwork .le. jwork) then
        write(*, "(A,A)") "Value of lwork, size of temporary buffer, ", &
            & "is too low, please increase it"
        stop 1
    end if
    allocate(ddx_data % far(ddx_data % nnfar), &
        & ddx_data % sfar(ddx_data % nclusters+1), &
        & ddx_data % near(ddx_data % nnnear), &
        & ddx_data % snear(ddx_data % nclusters+1))
    call tree_get_farnear(jwork, lwork, work, ddx_data % nclusters, &
        & ddx_data % nnfar, ddx_data % nfar, ddx_data % sfar, ddx_data % far, &
        & ddx_data % nnnear, ddx_data % nnear, ddx_data % snear, ddx_data % near)
    ! Get far-field P2P for a reference result
    do i = 1, nsph
        far_p2p(i) = zero
        do j = 1, nsph
            ok = .true.
            do k = ddx_data % snear(ddx_data % snode(i)), &
                & ddx_data % snear(ddx_data % snode(i)+1)-1
                if (ddx_data % near(k) .eq. ddx_data % snode(j)) ok = .false.
            end do
            if (ok) then
                far_p2p(i) = far_p2p(i) + &
                    & one/dnrm2(3, src(:, order(i))-src(:, order(j)), 1)
            end if
        end do
    end do
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(pm+pl, ddx_data % vscales, ddx_data % v4pi2lp1, &
        & ddx_data % vscales_rel)
    ddx_data % vfact(1) = one
    do i = 2, 2*(pm+pl)+1
        ddx_data % vfact(i) = ddx_data % vfact(i-1) * sqrt(dble(i-1))
    end do
    call fmm_constants(pm+pl, pm, pl, ddx_data % vcnk, &
        & ddx_data % m2l_ztranslate_coef, &
        & ddx_data % m2l_ztranslate_adj_coef)
    ! Init input harmonics
    do i = 1, nsph
        call fmm_p2m(src(:, order(i))-csph2(:, i), one, rsph2(i), pm, &
            & ddx_data % vscales, zero, node_m(:, ddx_data % snode(i)))
    end do
    ! Prepare M2M
    call tree_m2m_rotation(ddx_data, node_m)
    ! Check tree_m2l_rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_m2l_rotation"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || P2M + [tree M2L] + L2P - P2P ||", &
            & " / || P2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    call tree_m2l_rotation(ddx_data, node_m, node_l)
    call tree_l2l_rotation(ddx_data, node_l)
    do i = 1, nsph
        call fmm_l2p(src(:, order(i))-csph2(:, i), rsph2(i), pl, &
            & ddx_data % vscales_rel, one, node_l(:, ddx_data % snode(i)), &
            & zero, far_p2p2(i))
    end do
    rel_err = dnrm2(nsph, far_p2p-far_p2p2, 1) / dnrm2(nsph, far_p2p, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Reference value of M2L
    call tree_m2l_rotation(ddx_data, node_m, node_l)
    ! Allocate space for transfer matrices
    ddx_data % m2l_reflect_mat_size = (max(pm,pl)+1)*(2*max(pm,pl)+1)* &
        & (2*max(pm,pl)+3)/3
    ddx_data % m2l_ztranslate_mat_size = (min(pm,pl)+1)*(min(pm,pl)+2)* &
        & (3*max(pm,pl)+3-min(pm,pl))/6
    allocate(ddx_data % m2l_reflect_mat(ddx_data % m2l_reflect_mat_size, &
        & ddx_data % nnfar), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    allocate(ddx_data % m2l_ztranslate_mat(ddx_data % m2l_ztranslate_mat_size, &
        & ddx_data % nnfar), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    ! Check tree_m2l_reflection_use_mat
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_m2l_reflection_get_mat"
        write(*, "(A)") "Check tree_m2l_reflection_use_mat"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree M2L] - [tree_m2l_rotation] ||", &
            & " / || [tree_m2l_rotation] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", 1d-15
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    call tree_m2l_reflection_get_mat(ddx_data)
    call tree_m2l_reflection_use_mat(ddx_data, node_m, node_l2)
    rel_err = dnrm2((2*nsph-1)*((pl+1)**2), node_l2-node_l, 1) / &
        & dnrm2((2*nsph-1)*((pl+1)**2), node_l, 1)
    ok = rel_err .le. 1d-15
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Check tree_m2l_reflection_use_mat_adj
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_m2l_reflection_use_mat_adj"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [tree M2L] - [adj tree M2L]^T ||", &
            & " / || [tree M2L] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", 1d-15
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, ddx_data % nclusters
        do j = 1, (pm+1)**2
            node_m2 = zero
            node_m2(j, i) = one
            call tree_m2l_reflection_use_mat(ddx_data, node_m2, node_l2)
            full_mat(:, :, j, i) = node_l2
        end do
    end do
    full_norm = dnrm2((((2*nsph-1)*((pm+1)*(pl+1)))**2), full_mat, 1)
    do i = 1, ddx_data % nclusters
        do j = 1, (pl+1)**2
            node_l2 = zero
            node_l2(j, i) = one
            call tree_m2l_reflection_use_mat_adj(ddx_data, node_l2, node_m2)
            full_mat(j, i, :, :) = full_mat(j, i, :, :) - node_m2
        end do
    end do
    diff_norm = dnrm2((((2*nsph-1)*((pm+1)*(pl+1)))**2), full_mat, 1)
    rel_err = diff_norm / full_norm
    ok = rel_err .le. 1d-15
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Check tree_m2l_rotation_adj
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check tree_m2l_rotation_adj"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err = || [adj tree M2L] - [adj tree M2L] ||", &
            & " / || [adj tree M2L] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", 1d-15
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") " ok | err"
        write(*, "(A)") repeat("=", 40)
    end if
    call tree_m2l_reflection_use_mat_adj(ddx_data, node_l, node_m)
    call tree_m2l_rotation_adj(ddx_data, node_l, node_m2)
    rel_err = dnrm2((2*nsph-1)*((pm+1)**2), node_m2-node_m, 1) / &
        & dnrm2((2*nsph-1)*((pm+1)**2), node_m, 1)
    ok = rel_err .le. 1d-15
    if (iprint .gt. 0) then
        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    ! Deallocate tree
    call ddfree(ddx_data)
end subroutine check_tree_m2l

subroutine check_tree_l2p(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
end subroutine check_tree_l2p

subroutine check_tree_m2p
end subroutine check_tree_m2p

end program test_ddx_core

