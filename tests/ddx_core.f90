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

character(len=255) :: testname, p_string
integer :: argc
integer :: p=30, i, j

integer, parameter :: nx=13, nrand=10
real(dp) :: x(3, nx)
!real(dp) :: alpha(4)=(/1d0, -1d0, 1d-294, 1d+294/)
!real(dp) :: alpha(4)=(/1d0, -1d0, 1d-10, 1d+2/)
real(dp) :: alpha(2)=(/1d0, -1d0/)
real(dp), external :: dnrm2

! Read argument (which tests to run)
argc = command_argument_count()
if (argc .eq. 0) then
    testname = "all"
else
    ! Read test name
    call get_command_argument(1, testname)
    select case(testname)
        case ('all')
        case ('init')
        case ('polleg')
        case ('ylmbas')
        case ('m2p')
        case ('m2p_bessel')
        case ('m2p_m2l')
        case ('m2p_adj')
        case ('m2p_bessel_grad')
        case ('l2p')
        case ('l2p_bessel')
        case ('l2p_adj')
        case ('l2p_bessel_grad')
        case ('m2m')
        case ('m2m_bessel')
        case ('m2m_adj')
        case ('m2m_bessel_adj')
        case ('l2l')
        case ('l2l_bessel')
        case ('l2l_adj')
        case ('l2l_bessel_adj')
        case ('m2l')
        case ('m2l_bessel')
        case ('m2l_adj')
        case ('m2l_bessel_adj')
        case ('tree_init')
!        case ('tree_m2m')
!        case ('tree_l2l')
!        case ('tree_m2l')
        case default
            stop "Wrong testname value"
    end select
    ! Read maximal degree p
    if (argc .gt. 1) then
        call get_command_argument(2, p_string)
        read(p_string, "(I3)") p
    end if
end if

! Set points x
x(:, 1) = zero
x(:, 2) = (/6d-1, 0d0, 0d0/)
x(:, 3) = (/-4d-1, 0d0, 0d0/)
x(:, 4) = (/0d0, 5d-1, 0d0/)
x(:, 5) = (/0d0, -7d-1, 0d0/)
x(:, 6) = (/0d0, 0d0, 8d-1/)
x(:, 7) = (/0d0, 0d0, -4d-1/)
x(:, 8) = (/3d-1, -2d-1, 0d0/)
x(:, 9) = (/-2d-1, -3d-1, 0d0/)
x(:, 10) = (/3d-1, 0d0, 6d-1/)
x(:, 11) = (/-3d-1, 0d0, 7d-1/)
x(:, 12) = (/3d-1, 5d-1, 3d-1/)
x(:, 13) = (/-2d-1, 4d-1, 3d-1/)

! Check correctness of info for valid and invalid input parameters of ddinit
if ((testname .eq. 'all') .or. (testname .eq. 'init')) then
    call check_ddinit_args()
end if

! Check polleg (Legendre polynomials)
if ((testname .eq. 'all') .or. (testname .eq. 'polleg')) then
    do i = 0, 2*p
        call check_polleg(i)
    end do
end if

! Check ylmbas (spherical harmonics)
if ((testname .eq. 'all') .or. (testname .eq. 'ylmbas')) then
    do i = 0, 2*p
        call check_ylmbas(i)
    end do
end if

! Check M2P
if ((testname .eq. 'all') .or. (testname .eq. 'm2p')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_m2p(j, alpha(i))
        end do
    end do
end if

! Check M2P Bessel
if ((testname .eq. 'all') .or. (testname .eq. 'm2p_bessel')) then
    do i = 1, size(alpha)
        do j = 0, 20
            call check_m2p_bessel(j, alpha(i))
        end do
    end do
end if

! Check M2P against M2L with pl=0
if ((testname .eq. 'all') .or. (testname .eq. 'm2p_m2l')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_m2p_m2l(j, alpha(i))
        end do
    end do
end if

! Check adjoint M2P
if ((testname .eq. 'all') .or. (testname .eq. 'm2p_adj')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_m2p_adj(j, alpha(i))
        end do
    end do
end if

! Check M2P Bessel gradient
!if ((testname .eq. 'all') .or. (testname .eq. 'm2p_bessel')) then
if (testname .eq. 'm2p_bessel_grad') then
    !do i = 1, size(alpha)
        !do j = 0, 20
            j = 5
            call check_m2p_bessel_grad(j)
        !end do
    !end do
end if

! Check L2P
if ((testname .eq. 'all') .or. (testname .eq. 'l2p')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_l2p(j, alpha(i))
        end do
    end do
end if

! Check L2P Bessel
if ((testname .eq. 'all') .or. (testname .eq. 'l2p_bessel')) then
    do i = 1, size(alpha)
        do j = 0, 20
            call check_l2p_bessel(j, alpha(i))
        end do
    end do
end if

! Check adjoint L2P
if ((testname .eq. 'all') .or. (testname .eq. 'l2p_adj')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_l2p_adj(j, alpha(i))
        end do
    end do
end if

! Check L2P Bessel gradient
!if ((testname .eq. 'all') .or. (testname .eq. 'l2p_bessel_grad')) then
if (testname .eq. 'l2p_bessel_grad') then
    !do i = 1, size(alpha)
        !do j = 0, 20
            j = 5
            call check_l2p_bessel_grad(j)
        !end do
    !end do
end if

! Check M2M
if ((testname .eq. 'all') .or. (testname .eq. 'm2m')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_m2m(j, alpha(i))
        end do
    end do
end if

! Check M2M Bessel manually
!if ((testname .eq. 'all') .or. (testname .eq. 'm2m_bessel')) then
if (testname .eq. 'm2m_bessel') then
    !do i = 1, size(alpha)
    j = 20
        !do j = 0, p
            call check_m2m_bessel(j)
        !end do
    !end do
end if

! Check adjoint M2M
if ((testname .eq. 'all') .or. (testname .eq. 'm2m_adj')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_m2m_adj(j, alpha(i))
        end do
    end do
end if

! Check adjoint M2M Bessel
!if ((testname .eq. 'all') .or. (testname .eq. 'm2m_bessel_adj')) then
if (testname .eq. 'm2m_bessel_adj') then
    do i = 1, 1 !size(alpha)
        do j = 0, p
            call check_m2m_bessel_adj(j, alpha(i))
        end do
    end do
end if

! Check L2L
if ((testname .eq. 'all') .or. (testname .eq. 'l2l')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_l2l(j, alpha(i))
        end do
    end do
end if

! Check L2L Bessel manually
!if ((testname .eq. 'all') .or. (testname .eq. 'l2l_bessel')) then
if (testname .eq. 'l2l_bessel') then
    !do i = 1, size(alpha)
    j = 10
        !do j = 0, p
            call check_l2l_bessel(j)
        !end do
    !end do
end if

! Check adjoint L2L
if ((testname .eq. 'all') .or. (testname .eq. 'l2l_adj')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_l2l_adj(j, alpha(i))
        end do
    end do
end if

! Check adjoint L2L Bessel
!if ((testname .eq. 'all') .or. (testname .eq. 'l2l_bessel_adj')) then
if (testname .eq. 'l2l_bessel_adj') then
    do i = 1, 1 !size(alpha)
        do j = 0, p
            call check_l2l_bessel_adj(j, alpha(i))
        end do
    end do
end if

! Check M2L
if ((testname .eq. 'all') .or. (testname .eq. 'm2l')) then
    do i = 1, size(alpha)
        do j = 0, p
        !do j = 1, 1
            call check_m2l(j, j, alpha(i))
        end do
        do j = 1, p
            call check_m2l(0, j, alpha(i))
            call check_m2l(j, 0, alpha(i))
        end do
    end do
end if

! Check M2L Bessel
!if ((testname .eq. 'all') .or. (testname .eq. 'm2l_bessel')) then
if (testname .eq. 'm2l_bessel') then
    j = 10
    call check_m2l_bessel(j)
end if

! Check adjoint M2L
if ((testname .eq. 'all') .or. (testname .eq. 'm2l_adj')) then
    do i = 1, size(alpha)
        do j = 0, p
            call check_m2l_adj(j, j, alpha(i))
        end do
        do j = 1, p
            call check_m2l_adj(0, j, alpha(i))
            call check_m2l_adj(j, 0, alpha(i))
        end do
    end do
end if

! Check adjoint M2L Bessel
!if ((testname .eq. 'all') .or. (testname .eq. 'm2l_bessel_adj')) then
if (testname .eq. 'm2l_bessel_adj') then
    do i = 1, 1 !size(alpha)
        do j = 0, p
            call check_m2l_bessel_adj(j, alpha(i))
        end do
    end do
end if

! Check recursive inertial tree
if ((testname .eq. 'all') .or. (testname .eq. 'tree_init')) then
    do i = 1, size(alpha)
        call check_tree_rib(alpha(i))
    end do
end if

!! Check tree M2M
!if ((testname .eq. 'all') .or. (testname .eq. 'tree_m2m')) then
!    do i = 1, size(alpha)
!        do j = 0, p
!            call check_tree_m2m(j, alpha(i))
!        end do
!    end do
!end if
!
!! Check tree L2L
!if ((testname .eq. 'all') .or. (testname .eq. 'tree_l2l')) then
!    do i = 1, size(alpha)
!        do j = 0, p
!            call check_tree_l2l(j, alpha(i))
!        end do
!    end do
!end if
!
!! Check tree M2L
!if ((testname .eq. 'all') .or. (testname .eq. 'tree_m2l')) then
!    do i = 1, size(alpha)
!        call check_tree_m2l(0, 0, alpha(i), 6d-2)
!        call check_tree_m2l(1, 0, alpha(i), 4d-2)
!        call check_tree_m2l(10, 0, alpha(i), 4d-2)
!        call check_tree_m2l(0, 1, alpha(i), 3d-2)
!        call check_tree_m2l(1, 1, alpha(i), 4d-3)
!        call check_tree_m2l(10, 1, alpha(i), 4d-3)
!        call check_tree_m2l(0, 10, alpha(i), 3d-2)
!        call check_tree_m2l(1, 10, alpha(i), 3d-3)
!        call check_tree_m2l(10, 10, alpha(i), 4d-9)
!        !call check_tree_m2l(20, 20, alpha(i), 1d-14)
!    end do
!end if

contains

subroutine check_ddinit_args()
    ! Example of correct args
    integer :: n=1, model=1, lmax=1, ngrid=1202, force=1, fmm=1, pm=0, pl=0, &
        & fmm_precompute=0, iprint=0, matvecmem=0, maxiter=10, &
        & jacobi_ndiis=10, nproc=1
    real(dp) :: charge(10), x(10), y(10), z(10), rvdw(10), se=zero, eta=1d-1, &
        & eps=1.1d1, kappa=1d0
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
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "correct test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check different correct inputs with different n <= 10 (hardcoded value)
    do i = 1, 10
        call ddinit(i, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, &
            & pm, pl, se, eta, eps, kappa, matvecmem, &
            & maxiter, jacobi_ndiis, nproc, &
            & '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`nsph` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    ! Check incorrect input nsph <= 0
    i = 0
    call ddinit(i, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`nsph` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = -1
    call ddinit(i, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`nsph` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check all possible models (1, 2, 3) with other correct inputs
    do i = 1, 3
        write(*, *) "model=", i
        call ddinit(n, charge, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
            & pl, se, eta, eps, kappa, matvecmem, &
            & maxiter, jacobi_ndiis, nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`model` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    ! Check incorrect models
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`model` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = 4
    call ddinit(n, charge, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`model` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct lmax >= 0
    !do i = 0, 6
    do i = 1, 6
        call ddinit(n, charge, x, y, z, rvdw, model, i, ngrid, force, fmm, &
            & pm, pl, se, eta, eps, kappa, &
            & matvecmem, maxiter, jacobi_ndiis, &
            & nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`lmax` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    ! Check incorrect lmax < 0
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, i, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`lmax` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct ngrid >= 0
    do i = 1, nllg
        j = ng0(i)
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, j, force, fmm, pm, &
            & pl, se, eta, eps, kappa, matvecmem, &
            & maxiter, jacobi_ndiis, nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`ngrid` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    ! Check incorrect ngrid < 0
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, i, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`ngrid` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct force (0, 1)
    do i = 0, 1
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
            & pl, se, eta, eps, kappa, matvecmem, &
            & maxiter, jacobi_ndiis, nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`force` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    ! Check incorrect force
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`force` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = 2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`force` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct fmm (0, 1)
    do i = 0, 1
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, i, &
            & pm, pl, se, eta, eps, kappa, &
            & matvecmem, maxiter, jacobi_ndiis, nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`fmm` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    ! Check incorrect fmm
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`fmm` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = 2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`fmm` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct pm (ignored if fmm=0)
    j = 0
    do i = -2, 2
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
            & i, pl, se, eta, eps, kappa, matvecmem, &
            & maxiter, jacobi_ndiis, nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`pm` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    ! Check correct pm (fmm=1)
    j = 1
    do i = 0, 20, 5
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
            & i, pl, se, eta, eps, kappa, matvecmem, &
            & maxiter, jacobi_ndiis, nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`pm` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
        & i, pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`pm` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check incorrect pm (fmm=1)
    j = 1
    i = -2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, i, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`pm` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct pl (ignored if fmm=0)
    j = 0
    do i = -2, 2
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
            & pm, i, se, eta, eps, kappa, matvecmem, &
            & maxiter, jacobi_ndiis, nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`pl` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    ! Check correct pl (fmm=1)
    j = 1
    do i = 0, 20, 5
        call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
            & pm, i, se, eta, eps, kappa, matvecmem, &
            & maxiter, jacobi_ndiis, nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`pl` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, &
        & pm, i, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`pl` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check incorrect pl (fmm=1)
    j = 1
    i = -2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, j, pm, &
        & i, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`pl` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct se (interval [-1,1])
    tmp = -one
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, tmp, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`se` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = zero
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, tmp, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`se` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = one
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, tmp, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`se` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check incorrect se
    tmp = 1.01d0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, tmp, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`se` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = -1.01d0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, tmp, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`se` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct eta (interval [0,1])
    tmp = pt5
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, tmp, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`eta` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = one
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, tmp, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`eta` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check incorrect eta
    tmp = -0.0000005
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, tmp, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`eta` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = 1.01d0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, tmp, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`eta` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = -1d-2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, tmp, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`eta` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct eps
    tmp = 1.01d0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, tmp, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`eps` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = dble(1000)
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, tmp, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`eps` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check incorrect eps
    tmp = zero
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, tmp, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`eps` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = pt5
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, tmp, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`eps` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = one
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, tmp, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`eps` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = -1d-2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, tmp, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`eps` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct kappa
    tmp = 1d-2
    j = 3 ! only referenced in case of LPB model
    call ddinit(n, charge, x, y, z, rvdw, j, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, tmp, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`kappa` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    tmp = -1d-2 ! not referenced in case of COSMO and PCM models
    do j = 1, 2
        call ddinit(n, charge, x, y, z, rvdw, j, lmax, ngrid, force, fmm, pm, &
            & pl, se, eta, eps, tmp, matvecmem, &
            & maxiter, jacobi_ndiis, nproc, '', ddx_data)
        if (ddx_data % error_flag .ne. 0) call error(-1, "`kappa` test failed in " // &
            & "check_ddinit_args()")
        call ddfree(ddx_data)
    end do
    ! Check incorrect kappa
    j = 3
    tmp = -1d-2
    call ddinit(n, charge, x, y, z, rvdw, j, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, tmp, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`kappa` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct matvecmem
    i = 0
    call ddinit(n, charge, x, y, z, rvdw, j, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, i, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`matvecmem` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = 1
    call ddinit(n, charge, x, y, z, rvdw, j, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, i, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`matvecmem` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, j, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, i, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`matvecmem` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = 2
    call ddinit(n, charge, x, y, z, rvdw, j, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, i, &
        & maxiter, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`matvecmem` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct maxiter
    i = 1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & i, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`maxiter` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = 1000000
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & i, jacobi_ndiis, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`maxiter` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check incorrect maxiter
    i = 0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & i, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`maxiter` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & i, jacobi_ndiis, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`maxiter` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct ndiis
    i = 0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, i, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`jacobi_ndiis` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = 1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, i, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`jacobi_ndiis` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = 1000
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, i, nproc, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`jacobi_ndiis` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check incorrect jacobi_ndiis
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, i, nproc, '', ddx_data)
    if (info .eq. 0) call error(-1, "`jacobi_ndiis` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check correct nproc
    i = 0
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, i, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`nproc` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = 1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, i, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`nproc` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    ! Check incorrect nproc
    i = 2
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, i, '', ddx_data)
    if (ddx_data % error_flag .ne. 0) call error(-1, "`nproc` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
    i = -1
    call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, i, '', ddx_data)
    if (info .eq. 0) call error(-1, "`nproc` test failed in " // &
        & "check_ddinit_args()")
    call ddfree(ddx_data)
end subroutine check_ddinit_args

! Check Legendre polynomials
subroutine check_polleg(p)
    ! Input
    integer, intent(in) :: p
    ! Local variables
    real(dp) :: c(3), rho, ctheta, stheta, cphi, sphi, vplm((p+1)**2), &
        & vplm2((p+1)**2), err
    real(dp), parameter :: threshold=2d-15
    logical :: ok
    integer :: i, j
    real(dp), external :: dnrm2
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check polleg"
    print "(A,I0)", " p=", p
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    ! Init output, since P_l^m for negative m are not set in polleg procedure
    vplm = zero
    vplm2 = zero
    ! Check against baseline
    do i = 1, nx
        c = x(:, i)
        call carttosph(c, rho, ctheta, stheta, cphi, sphi)
        call polleg_baseline(ctheta, stheta, p, vplm)
        call polleg(ctheta, stheta, p, vplm2)
        err = dnrm2((p+1)**2, vplm-vplm2, 1) / dnrm2((p+1)**2, vplm, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
end subroutine check_polleg

! Check spherical harmonics
subroutine check_ylmbas(p)
    ! Input
    integer, intent(in) :: p
    ! Local variables
    real(dp) :: c(3), rho, ctheta, stheta, cphi, sphi, vplm((p+1)**2), &
        & vylm((p+1)**2), vylm2((p+1)**2), vcos(p+1), vsin(p+1), err, &
        & vscales((p+1)**2), v4pi2lp1(p+1), vscales_rel((p+1)**2), &
        & threshold
    logical :: ok
    integer :: i, j
    real(dp), external :: dnrm2
    ! Set scaling factors
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check ylmbas"
    print "(A,I0)", " p=", p
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = dble((p+2)**2) * 2d-16
    ! Check against baseline
    do i = 2, nx
        c = x(:, i)
        call ylmbas_baseline(c, p, vscales, vylm, vplm, vcos, vsin)
        call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm2, &
            & vplm, vcos, vsin)
        err = dnrm2((p+1)**2, vylm-vylm2, 1) / dnrm2((p+1)**2, vylm, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
end subroutine check_ylmbas

! Check M2P
subroutine check_m2p(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    real(dp) :: y(3, nx), r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), src_m((p+1)**2, nrand), dst_v(nrand), &
        & dst_v2(nrand), err, threshold
    logical :: ok
    integer :: i, j, iseed(4)
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = alpha / sqrt(13d-2) * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_m randomly
    call dlarnv(3, iseed, nrand*((p+1)**2), src_m)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check M2P"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = dble(p+3) * 3d-16
    ! Check against the baseline, i=1 is ignored since it y(:,1)=zero
    do i = 2, nx
        c = y(:, i)
        ! Check alpha=beta=zero
        do j = 1, nrand
            call fmm_m2p(c, r, p, vscales_rel, zero, src_m(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=one, beta=zero
        do j = 1, nrand
            call fmm_m2p_baseline(c, r, p, vscales, one, src_m(:, j), zero, &
                & dst_v(j))
            call fmm_m2p(c, r, p, vscales_rel, one, src_m(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta=zero
        do j = 1, nrand
            call fmm_m2p_baseline(c, r, p, vscales, -three, src_m(:, j), &
                & zero, dst_v(j))
            call fmm_m2p(c, r, p, vscales_rel, -three, src_m(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta=one
        do j = 1, nrand
            call fmm_m2p(c, r, p, vscales_rel, zero, src_m(:, j), one, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta!={zero,one}
        do j = 1, nrand
            dst_v(j) = -pt5 * dst_v(j)
            call fmm_m2p(c, r, p, vscales_rel, zero, src_m(:, j), -pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta!={zero,one}
        do j = 1, nrand
            call fmm_m2p_baseline(c, r, p, vscales, -three, src_m(:, j), pt5, &
                & dst_v(j))
            call fmm_m2p(c, r, p, vscales_rel, -three, src_m(:, j), pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2p

! Check M2P Bessel
subroutine check_m2p_bessel(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    real(dp) :: y(3, nx), r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), src_m((p+1)**2, nrand), dst_v(nrand), &
        & dst_v2(nrand), err, threshold
    logical :: ok
    integer :: i, j, iseed(4)
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = alpha / sqrt(13d-2) * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_m randomly
    call dlarnv(3, iseed, nrand*((p+1)**2), src_m)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check M2P"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = dble(p+3) * 3d-14
    ! Check against the baseline, i=1 is ignored since y(:,1)=zero
    do i = 2, nx
        c = y(:, i)
        ! Check alpha=beta=zero
        do j = 1, nrand
            call fmm_m2p_bessel(c, r, p, vscales, zero, src_m(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=one, beta=zero
        do j = 1, nrand
            call fmm_m2p_bessel_baseline(c, r, p, vscales, one, src_m(:, j), zero, &
                & dst_v(j))
            call fmm_m2p_bessel(c, r, p, vscales, one, src_m(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta=zero
        do j = 1, nrand
            call fmm_m2p_bessel_baseline(c, r, p, vscales, -three, src_m(:, j), &
                & zero, dst_v(j))
            call fmm_m2p_bessel(c, r, p, vscales, -three, src_m(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta=one
        do j = 1, nrand
            call fmm_m2p_bessel(c, r, p, vscales, zero, src_m(:, j), one, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta!={zero,one}
        do j = 1, nrand
            dst_v(j) = -pt5 * dst_v(j)
            call fmm_m2p_bessel(c, r, p, vscales, zero, src_m(:, j), -pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta!={zero,one}
        do j = 1, nrand
            call fmm_m2p_bessel_baseline(c, r, p, vscales, -three, src_m(:, j), pt5, &
                & dst_v(j))
            call fmm_m2p_bessel(c, r, p, vscales, -three, src_m(:, j), pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2p_bessel

! Check M2P by M2L with pl=0
subroutine check_m2p_m2l(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    real(dp) :: y(3, nx), r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), src_m((p+1)**2, nrand), dst_v(nrand), &
        & dst_v2(nrand), err, &
        & vcnk((2*p+1)*(p+1)), &
        & m2l_ztranslate_coef(p+1, 1, 1), &
        & m2l_ztranslate_adj_coef(1, 1, p+1)
    logical :: ok
    integer :: i, j, iseed(4)
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = alpha / sqrt(13d-2) * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    call fmm_constants(p, p, 0, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_m randomly
    call dlarnv(3, iseed, nrand*((p+1)**2), src_m)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check M2P by M2L with pl=0"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    ! Check against the baseline, i=1 is ignored since it y(:,1)=zero
    do i = 2, nx
        c = y(:, i)
        ! Check alpha=beta=zero
        do j = 1, nrand
            call fmm_m2p(c, r, p, vscales_rel, zero, src_m(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        !if (.not. ok) stop 1
        ! Check alpha=one, beta=zero
        do j = 1, nrand
            call fmm_m2l_rotation(-c, r, r, p, 0, vscales, &
                & m2l_ztranslate_coef, vscales_rel(1), src_m(:, j), zero, &
                & dst_v(j))
            call fmm_m2p(c, r, p, vscales_rel, one, src_m(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. 2d-15
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta=zero
        do j = 1, nrand
            call fmm_m2l_rotation(-c, r, r, p, 0, vscales, &
                & m2l_ztranslate_coef, -three*vscales_rel(1), src_m(:, j), &
                & zero, dst_v(j))
            call fmm_m2p(c, r, p, vscales_rel, -three, src_m(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. 4d-15
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta=one
        do j = 1, nrand
            call fmm_m2p(c, r, p, vscales_rel, zero, src_m(:, j), one, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. 4d-15
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta!={zero,one}
        do j = 1, nrand
            dst_v(j) = -pt5 * dst_v(j)
            call fmm_m2p(c, r, p, vscales_rel, zero, src_m(:, j), -pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. 4d-15
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta!={zero,one}
        do j = 1, nrand
            call fmm_m2l_rotation(-c, r, r, p, 0, vscales, &
                & m2l_ztranslate_coef, -three*vscales_rel(1), src_m(:, j), &
                & pt5, dst_v(j))
            call fmm_m2p(c, r, p, vscales_rel, -three, src_m(:, j), pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. 4d-15
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2p_m2l

! Check adjoint M2P
subroutine check_m2p_adj(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    real(dp) :: y(3, nx), r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), dst_m((p+1)**2), dst_m2((p+1)**2), err, &
        & threshold
    logical :: ok
    integer :: i
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = alpha / sqrt(13d-2) * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check adjoint M2P"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = dble(p+3) * 3d-16
    ! Check against the baseline, i=1 is ignored since it y(:,1)=zero
    do i = 2, nx
        c = y(:, i)
        ! Check src_q=beta=zero
        call fmm_m2p_adj(c, zero, r, p, vscales_rel, zero, dst_m2)
        err = dnrm2((p+1)**2, dst_m2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check src_q=one, beta=zero
        call fmm_m2p_adj_baseline(c, one, r, p, vscales, zero, dst_m)
        call fmm_m2p_adj(c, one, r, p, vscales_rel, zero, dst_m2)
        err = dnrm2((p+1)**2, dst_m-dst_m2, 1) / dnrm2((p+1)**2, dst_m, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check src_q!={zero,one}, beta=zero
        call fmm_m2p_adj_baseline(c, -three, r, p, vscales, zero, dst_m)
        call fmm_m2p_adj(c, -three, r, p, vscales_rel, zero, dst_m2)
        err = dnrm2((p+1)**2, dst_m-dst_m2, 1) / dnrm2((p+1)**2, dst_m, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta=one
        call fmm_m2p_adj(c, zero, r, p, vscales_rel, one, dst_m2)
        err = dnrm2((p+1)**2, dst_m-dst_m2, 1) / dnrm2((p+1)**2, dst_m, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta!={zero,one}
        call fmm_m2p_adj_baseline(c, zero, r, p, vscales, -pt5, dst_m)
        call fmm_m2p_adj(c, zero, r, p, vscales_rel, -pt5, dst_m2)
        err = dnrm2((p+1)**2, dst_m-dst_m2, 1) / dnrm2((p+1)**2, dst_m, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta!={zero,one}
        call fmm_m2p_adj_baseline(c, -three, r, p, vscales, pt5, dst_m)
        call fmm_m2p_adj(c, -three, r, p, vscales_rel, pt5, dst_m2)
        err = dnrm2((p+1)**2, dst_m-dst_m2, 1) / dnrm2((p+1)**2, dst_m, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2p_adj

! Check M2P Bessel gradient
subroutine check_m2p_bessel_grad(p)
    ! Inputs
    integer, intent(in) :: p
    ! Local variables
    real(dp) :: x(3), r, src_c(3), src_r, dst_r, dst_v, dst_v2, dst_v3
    real(dp) :: src_m((p+1)**2), dst_m((p+2)**2), vscales((p+2)**2), &
        & vscales_rel((p+2)**2), v4pi2lp1(p+2), vcnk((2*p+3)*(p+2)), &
        & tmp(p+2), work(p+2), src_sk(p+2), sk(p+1), dk(p+1), basloc((p+1)**2), &
        & dbsloc(3, (p+1)**2), vplm((p+1)**2), vcos(p+1), vsin(p+1), &
        & dst_g1(3), dst_g2(3), dst_m_grad((p+2)**2, 3)
    complex(dp) :: work_complex(p+2)
    real(dp), external :: dnrm2
    integer :: i, l, m, indm
    type(ddx_params_type) :: params
    type(ddx_constants_type) :: constants
    ! Set params
    params % lmax = p
    constants % nbasis = (p+1)**2
    ! Compute special FMM constants
    call ylmscale(p+1, vscales, v4pi2lp1, vscales_rel)
    allocate(constants % vscales((p+2)**2))
    constants % vscales = vscales
    x = 7.5d0 * (/1.1d0, -2d0, one/)
    r = dnrm2(3, x, 1)
    src_r = one
    src_m = zero
    do i = 0, p
        !src_m(i*i+i+1) = 10d0 ** (-i)
        src_m(i*i+1:i*i+2*i+1) = 10d0 ** (-i)
    end do
    !src_m = one
    call modified_spherical_bessel_second_kind(p, r, sk, dk, &
        & work_complex)
    call modified_spherical_bessel_second_kind(p+1, src_r, src_sk, work, &
        & work_complex)
    call dbasis(params, constants, x/r, basloc, dbsloc, vplm, vcos, vsin)
    dst_g1 = zero
    do l = 0, p
        do m = -l, l
            indm = l*l+l+1+m
            dst_g1 = dst_g1 + src_m(indm)/src_sk(l+1)/r*(dk(l+1)*x* &
                & basloc(indm)+sk(l+1)*dbsloc(:, indm))
        end do
    end do
    write(*, *) dst_g1
    call fmm_m2p_bessel_grad(x, src_r, p, vscales, one, src_m, zero, dst_g2)
    write(*, *) dst_g2
    deallocate(constants % vscales)
end subroutine check_m2p_bessel_grad

! Check L2P
subroutine check_l2p(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    real(dp) :: y(3, nx), r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), src_l((p+1)**2, nrand), dst_v(nrand), &
        & dst_v2(nrand), err, threshold
    logical :: ok
    integer :: i, j, iseed(4)
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = 1.25d0 * alpha * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_m randomly
    call dlarnv(3, iseed, nrand*((p+1)**2), src_l)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check L2P"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = dble(p+1) * 2d-15
    ! Check against the baseline
    do i = 1, nx
        c = y(:, i)
        ! Check alpha=beta=zero
        do j = 1, nrand
            call fmm_l2p(c, r, p, vscales_rel, zero, src_l(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=one, beta=zero
        do j = 1, nrand
            call fmm_l2p_baseline(c, r, p, vscales, one, src_l(:, j), zero, &
                & dst_v(j))
            call fmm_l2p(c, r, p, vscales_rel, one, src_l(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta=zero
        do j = 1, nrand
            call fmm_l2p_baseline(c, r, p, vscales, -three, src_l(:, j), &
                & zero, dst_v(j))
            call fmm_l2p(c, r, p, vscales_rel, -three, src_l(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta=one
        do j = 1, nrand
            call fmm_l2p(c, r, p, vscales_rel, zero, src_l(:, j), one, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta!={zero,one}
        do j = 1, nrand
            call fmm_l2p_baseline(c, r, p, vscales, zero, src_l(:, j), -pt5, &
                & dst_v(j))
            call fmm_l2p(c, r, p, vscales_rel, zero, src_l(:, j), -pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta!={zero,one}
        do j = 1, nrand
            call fmm_l2p_baseline(c, r, p, vscales, -three, src_l(:, j), pt5, &
                & dst_v(j))
            call fmm_l2p(c, r, p, vscales_rel, -three, src_l(:, j), pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_l2p

! Check L2P Bessel
subroutine check_l2p_bessel(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    real(dp) :: y(3, nx), r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), src_l((p+1)**2, nrand), dst_v(nrand), &
        & dst_v2(nrand), err, threshold
    logical :: ok
    integer :: i, j, iseed(4)
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = 1.25d0 * alpha * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_m randomly
    call dlarnv(3, iseed, nrand*((p+1)**2), src_l)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check L2P"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = dble(p+1) * 2d-15
    ! Check against the baseline
    do i = 1, nx
        c = y(:, i)
        ! Check alpha=beta=zero
        do j = 1, nrand
            call fmm_l2p_bessel(c, r, p, vscales, zero, src_l(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=one, beta=zero
        do j = 1, nrand
            call fmm_l2p_bessel_baseline(c, r, p, vscales, one, src_l(:, j), zero, &
                & dst_v(j))
            call fmm_l2p_bessel(c, r, p, vscales, one, src_l(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta=zero
        do j = 1, nrand
            call fmm_l2p_bessel_baseline(c, r, p, vscales, -three, src_l(:, j), &
                & zero, dst_v(j))
            call fmm_l2p_bessel(c, r, p, vscales, -three, src_l(:, j), zero, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta=one
        do j = 1, nrand
            call fmm_l2p_bessel(c, r, p, vscales, zero, src_l(:, j), one, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha=zero, beta!={zero,one}
        do j = 1, nrand
            call fmm_l2p_bessel_baseline(c, r, p, vscales, zero, src_l(:, j), -pt5, &
                & dst_v(j))
            call fmm_l2p_bessel(c, r, p, vscales, zero, src_l(:, j), -pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        ! Check alpha!={zero,one}, beta!={zero,one}
        do j = 1, nrand
            call fmm_l2p_bessel_baseline(c, r, p, vscales, -three, src_l(:, j), pt5, &
                & dst_v(j))
            call fmm_l2p_bessel(c, r, p, vscales, -three, src_l(:, j), pt5, &
                & dst_v2(j))
        end do
        err = dnrm2(nrand, dst_v-dst_v2, 1) / dnrm2(nrand, dst_v, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_l2p_bessel

! Check adjoint L2P
subroutine check_l2p_adj(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    real(dp) :: y(3, nx), r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), dst_l((p+1)**2), dst_l2((p+1)**2), err, &
        & threshold
    logical :: ok
    integer :: i
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = 1.25d0 * alpha * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check adjoint L2P"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = dble(p+2) * 2d-16
    ! Check against the baseline
    do i = 1, nx
        c = y(:, i)
        call fmm_l2p_adj(c, zero, r, p, vscales_rel, zero, dst_l2)
        err = dnrm2((p+1)**2, dst_l2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        call fmm_l2p_adj_baseline(c, one, r, p, vscales, zero, dst_l)
        call fmm_l2p_adj(c, one, r, p, vscales_rel, zero, dst_l2)
        err = dnrm2((p+1)**2, dst_l-dst_l2, 1) / dnrm2((p+1)**2, dst_l, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        call fmm_l2p_adj_baseline(c, -three, r, p, vscales, zero, dst_l)
        call fmm_l2p_adj(c, -three, r, p, vscales_rel, zero, dst_l2)
        err = dnrm2((p+1)**2, dst_l-dst_l2, 1) / dnrm2((p+1)**2, dst_l, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        call fmm_l2p_adj(c, zero, r, p, vscales_rel, one, dst_l2)
        err = dnrm2((p+1)**2, dst_l-dst_l2, 1) / dnrm2((p+1)**2, dst_l, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        call fmm_l2p_adj_baseline(c, zero, r, p, vscales, -pt5, dst_l)
        call fmm_l2p_adj(c, zero, r, p, vscales_rel, -pt5, dst_l2)
        err = dnrm2((p+1)**2, dst_l-dst_l2, 1) / dnrm2((p+1)**2, dst_l, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
        call fmm_l2p_adj_baseline(c, -three, r, p, vscales, pt5, dst_l)
        call fmm_l2p_adj(c, -three, r, p, vscales_rel, pt5, dst_l2)
        err = dnrm2((p+1)**2, dst_l-dst_l2, 1) / dnrm2((p+1)**2, dst_l, 1)
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_l2p_adj

! Check L2P Bessel gradient
subroutine check_l2p_bessel_grad(p)
    ! Inputs
    integer, intent(in) :: p
    ! Local variables
    real(dp) :: x(3), r, src_c(3), src_r, dst_r, dst_v, dst_v2, dst_v3
    real(dp) :: src_l((p+1)**2), dst_l((p+2)**2), vscales((p+2)**2), &
        & vscales_rel((p+2)**2), v4pi2lp1(p+2), vcnk((2*p+3)*(p+2)), &
        & tmp(p+2), work(p+2), src_si(p+2), si(p+1), di(p+1), basloc((p+1)**2), &
        & dbsloc(3, (p+1)**2), vplm((p+1)**2), vcos(p+1), vsin(p+1), &
        & dst_g1(3), dst_g2(3), dst_l_grad((p+2)**2, 3)
    complex(dp) :: work_complex(p+2)
    real(dp), external :: dnrm2
    integer :: i, l, m, indm
    type(ddx_params_type) :: params
    type(ddx_constants_type) :: constants
    ! Set params
    params % lmax = p
    constants % nbasis = (p+1)**2
    ! Compute special FMM constants
    call ylmscale(p+1, vscales, v4pi2lp1, vscales_rel)
    allocate(constants % vscales((p+2)**2))
    constants % vscales = vscales
    x = 7.5d-1 * (/one, zero, zero/)
    r = dnrm2(3, x, 1)
    src_r = one
    src_l = zero
    do i = 0, p
        !src_l(i*i+i+1) = 10d0 ** (-i)
        src_l(i*i+1:i*i+2*i+1) = 10d0 ** (-i)
    end do
    !src_l = one
    call modified_spherical_bessel_first_kind(p, r, si, di, &
        & work_complex)
    call modified_spherical_bessel_first_kind(p+1, src_r, src_si, work, &
        & work_complex)
    call dbasis(params, constants, x/r, basloc, dbsloc, vplm, vcos, vsin)
    dst_g1 = zero
    do l = 0, p
        do m = -l, l
            indm = l*l+l+1+m
            dst_g1 = dst_g1 + src_l(indm)/src_si(l+1)/r*(di(l+1)*x* &
                & basloc(indm)+si(l+1)*dbsloc(:, indm))
        end do
    end do
    write(*, *) dst_g1
    call fmm_l2p_bessel_grad(x, src_r, p, vscales, one, src_l, zero, dst_g2)
    write(*, *) dst_g2
    deallocate(constants % vscales)
end subroutine check_l2p_bessel_grad

! Check M2M
subroutine check_m2m(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    integer, parameter :: nrand=10
    real(dp) :: y(3, nx), r, dst_r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), vcnk((2*p+1)*(p+1)), &
        & m2l_ztranslate_coef(p+1, 1, 1), &
        & m2l_ztranslate_adj_coef(1, 1, p+1), &
        & src_m((p+1)**2, nrand), dst_m((p+1)**2, nrand), &
        & dst_m2((p+1)**2, nrand), err, threshold
    logical :: ok
    integer :: i, j, iseed(4), ndst_m
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = alpha * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    call fmm_constants(p, p, 0, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_m randomly
    ndst_m = nrand * ((p+1)**2)
    call dlarnv(3, iseed, ndst_m, src_m)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check M2M"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = (p+3) * 1d-15
    ! Check against baseline implementation
    do i = 1, nx
        c = y(:, i)
        dst_r = r + dnrm2(3, c, 1)
        ! Check alpha=beta=zero
        do j = 1, nrand
            call fmm_m2m_rotation(c, r, dst_r, p, vscales, vcnk, zero, &
                & src_m(:, j), zero, dst_m2(:, j))
        end do
        err = dnrm2(ndst_m, dst_m2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha=one, beta=zero
        dst_m = zero
        do j = 1, nrand
            call fmm_m2m_baseline(c, r, dst_r, p, vscales, src_m(:, j), &
                & dst_m(:, j))
            call fmm_m2m_rotation(c, r, dst_r, p, vscales, vcnk, one, &
                & src_m(:, j), zero, dst_m2(:, j))
        end do
        err = dnrm2(ndst_m, dst_m-dst_m2, 1) / dnrm2(ndst_m, dst_m, 1)
        ok = err .le. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        ! Check alpha!={zero,one}, beta=zero
        do j = 1, nrand
            call fmm_m2m_rotation(c, r, dst_r, p, vscales, vcnk, -three, &
                & src_m(:, j), zero, dst_m2(:, j))
        end do
        err = dnrm2(ndst_m, three*dst_m+dst_m2, 1) / dnrm2(ndst_m, dst_m, 1)
        ok = err .le. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha=zero, beta=one
        dst_m2 = dst_m
        do j = 1, nrand
            call fmm_m2m_rotation(c, r, dst_r, p, vscales, vcnk, zero, &
                & src_m(:, j), one, dst_m2(:, j))
        end do
        err = dnrm2(ndst_m, dst_m-dst_m2, 1) / dnrm2(ndst_m, dst_m, 1)
        ok = err .le. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha=zero, beta!={zero,one}
        dst_m2 = dst_m
        do j = 1, nrand
            call fmm_m2m_rotation(c, r, dst_r, p, vscales, vcnk, zero, &
                & src_m(:, j), -pt5, dst_m2(:, j))
        end do
        err = dnrm2((p+1)**2, pt5*dst_m+dst_m2, 1) / &
            & dnrm2((p+1)**2, dst_m, 1) / four
        ok = err .le. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha!={zero,one}, beta!={zero,one}
        dst_m2 = two * dst_m
        do j = 1, nrand
            call fmm_m2m_rotation(c, r, dst_r, p, vscales, vcnk, -three, &
                & src_m(:, j), pt5, dst_m2(:, j))
        end do
        err = dnrm2((p+1)**2, two*dst_m+dst_m2, 1) / &
            & dnrm2((p+1)**2, dst_m, 1) / four
        ok = err .le. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2m

subroutine check_m2m_bessel(p)
    ! Inputs
    integer, intent(in) :: p
    ! Local variables
    real(dp) :: x(3), src_c(3), src_r, dst_r, dst_v, dst_v2, dst_v3
    real(dp) :: src_m((p+1)**2), dst_m((p+1)**2), vscales((p+1)**2), &
        & vscales_rel((p+1)**2), v4pi2lp1(p+1), vcnk((2*p+1)*(p+1)), &
        tmp(p+1), work(p+1), kappa
    complex(dp) :: work_complex(max(2, p+1))
    real(dp), external :: dnrm2
    integer :: i
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    x = 7.5d0 * (/one, -1.1d0, one/)
    src_c = 1.01d0 * (/-one, one, one/)
    src_r = one
    dst_r = src_r + dnrm2(3, src_c, 1)
    src_m = zero
    do i = 3, p
        !src_m(i*i+i+1) = 10d0 ** (-i)
        !src_m(i*i+1:i*i+2*i+1) = 10d0 ** (-i)
    end do
    src_m = one
    kappa = 1d-2
    call fmm_m2p_bessel_baseline(kappa*(x-src_c), kappa*src_r, p, vscales, one, src_m, &
        & zero, dst_v)
    call fmm_m2m_bessel_rotation(src_c, src_r, dst_r, kappa, p, vscales, vcnk, one, &
        & src_m, zero, dst_m)
    call fmm_m2p_bessel_baseline(kappa*x, kappa*dst_r, p, vscales, one, dst_m, zero, &
        & dst_v2)
    call fmm_m2m_bessel_rotation(-src_c, dst_r, src_r, kappa, p, vscales, vcnk, one, &
        & dst_m, zero, src_m)
    call fmm_m2p_bessel_baseline(kappa*(x-src_c), kappa*src_r, p, vscales, one, src_m, &
        & zero, dst_v3)
    write(*, *) dst_v, dst_v2, dst_v3
end subroutine check_m2m_bessel

! Check adjoint M2M
subroutine check_m2m_adj(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    integer, parameter :: nrand=10
    real(dp) :: y(3, nx), r, dst_r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), vcnk((2*p+1)*(p+1)), &
        & m2l_ztranslate_coef(p+1, 1, 1), &
        & m2l_ztranslate_adj_coef(1, 1, p+1), &
        & src_m((p+1)**2, nrand), src_m2((p+1)**2, nrand), &
        & dst_m((p+1)**2, nrand), err, tmp(nrand, nrand)
    logical :: ok
    integer :: i, j, iseed(4), ndst_m
    real(dp), external :: dnrm2
    real(dp), parameter :: threshold=4d-15
    ! Copy templated points with a proper multiplier
    y = alpha * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    call fmm_constants(p, p, 0, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_m randomly
    ndst_m = nrand * ((p+1)**2)
    call dlarnv(3, iseed, ndst_m, src_m)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check adjoint M2M"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    ! Check against direct M2M
    do i = 1, nx
        c = y(:, i)
        dst_r = r + dnrm2(3, y(:, i), 1)
        do j = 1, nrand
            call fmm_m2m_rotation(c, r, dst_r, p, vscales, vcnk, one, &
                & src_m(:, j), zero, dst_m(:, j))
            call fmm_m2m_rotation_adj(-c, dst_r, r, p, vscales, vcnk, one, &
                & dst_m(:, j), zero, src_m2(:, j))
        end do
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, one, dst_m, (p+1)**2, &
            & dst_m, (p+1)**2, zero, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1)
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, -one, src_m2, (p+1)**2, &
            & src_m, (p+1)**2, one, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1) / err
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2m_adj

! Check adjoint M2M Bessel
subroutine check_m2m_bessel_adj(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    integer, parameter :: nrand=10
    real(dp) :: y(3, nx), r, dst_r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), vcnk((2*p+1)*(p+1)), &
        & m2l_ztranslate_coef(p+1, 1, 1), &
        & m2l_ztranslate_adj_coef(1, 1, p+1), &
        & src_m((p+1)**2, nrand), src_m2((p+1)**2, nrand), &
        & dst_m((p+1)**2, nrand), err, tmp(nrand, nrand), kappa
    logical :: ok
    integer :: i, j, iseed(4), ndst_m
    real(dp), external :: dnrm2
    real(dp), parameter :: threshold=4d-15
    ! Copy templated points with a proper multiplier
    y = alpha * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    call fmm_constants(p, p, 0, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    kappa = 1d+1
    ! Generate src_m randomly
    ndst_m = nrand * ((p+1)**2)
    call dlarnv(3, iseed, ndst_m, src_m)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check adjoint M2M Bessel"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    ! Check against direct M
    ! Ignore i=1 for now
    do i = 2, nx
        c = y(:, i)
        dst_r = r + dnrm2(3, y(:, i), 1)
        do j = 1, nrand
            call fmm_m2m_bessel_rotation(c, r, dst_r, kappa, p, vscales, vcnk, one, &
                & src_m(:, j), zero, dst_m(:, j))
            call fmm_m2m_bessel_rotation_adj(-c, dst_r, r, kappa, p, vscales, vcnk, one, &
                & dst_m(:, j), zero, src_m2(:, j))
        end do
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, one, dst_m, (p+1)**2, &
            & dst_m, (p+1)**2, zero, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1)
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, -one, src_m2, (p+1)**2, &
            & src_m, (p+1)**2, one, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1) / err
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2m_bessel_adj

! Check L2L
subroutine check_l2l(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    integer, parameter :: nrand=10
    real(dp) :: y(3, nx), r, src_r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), vfact(2*p+1), &
        & src_l((p+1)**2, nrand), dst_l((p+1)**2, nrand), &
        & dst_l2((p+1)**2, nrand), err
    logical :: ok
    integer :: i, j, iseed(4), ndst_l
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = alpha * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    vfact(1) = one
    do i = 2, 2*p+1
        vfact(i) = vfact(i-1) * sqrt(dble(i-1))
    end do
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_l randomly
    ndst_l = nrand * ((p+1)**2)
    call dlarnv(3, iseed, ndst_l, src_l)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check L2L"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    ! Check against baseline implementation
    do i = 1, nx
        c = y(:, i)
        src_r = r + dnrm2(3, c, 1)
        ! Check alpha=beta=zero
        do j = 1, nrand
            call fmm_l2l_rotation(c, src_r, r, p, vscales, vfact, zero, &
                & src_l(:, j), zero, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, dst_l2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha=one, beta=zero
        dst_l = zero
        do j = 1, nrand
            call fmm_l2l_baseline(c, src_r, r, p, vscales, src_l(:, j), &
                & dst_l(:, j))
            call fmm_l2l_rotation(c, src_r, r, p, vscales, vfact, one, &
                & src_l(:, j), zero, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, dst_l-dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 4d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha!={zero,one}, beta=zero
        do j = 1, nrand
            call fmm_l2l_rotation(c, src_r, r, p, vscales, vfact, -three, &
                & src_l(:, j), zero, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, three*dst_l+dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 4d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha=zero, beta=one
        dst_l2 = dst_l
        do j = 1, nrand
            call fmm_l2l_rotation(c, src_r, r, p, vscales, vfact, zero, &
                & src_l(:, j), one, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, dst_l-dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 4d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha=zero, beta!={zero,one}
        dst_l2 = dst_l
        do j = 1, nrand
            call fmm_l2l_rotation(c, src_r, r, p, vscales, vfact, zero, &
                & src_l(:, j), -pt5, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, pt5*dst_l+dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 4d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha!={zero,one}, beta!={zero,one}
        dst_l2 = two * dst_l
        do j = 1, nrand
            call fmm_l2l_rotation(c, src_r, r, p, vscales, vfact, -three, &
                & src_l(:, j), pt5, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, two*dst_l+dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 4d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_l2l

subroutine check_l2l_bessel(p)
    ! Inputs
    integer, intent(in) :: p
    ! Local variables
    real(dp) :: x(3), src_c(3), src_r, dst_r, dst_v, dst_v2, dst_v3
    real(dp) :: src_l((p+1)**2), dst_l((p+1)**2), vscales((p+1)**2), &
        & vscales_rel((p+1)**2), v4pi2lp1(p+1), vcnk((2*p+1)*(p+1)), &
        tmp(p+1), work(p+1), kappa
    complex(dp) :: work_complex(max(2, p+1))
    real(dp), external :: dnrm2
    integer :: i
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    x = 7.5d-1 * (/one, -1.1d0, one/)
    src_c = 1.01d0 * (/-one, one, one/)
    src_r = one
    dst_r = src_r + dnrm2(3, src_c, 1)
    src_l = zero
    do i = 0, p
        !src_l(i*i+i+1) = 10d0 ** (-i)
        !src_l(i*i+1:i*i+2*i+1) = 10d0 ** (-i)
    end do
    src_l = one
    kappa = 1d-1
    call fmm_l2p_bessel_baseline(kappa*(x-src_c), kappa*src_r, p, vscales, one, src_l, &
        & zero, dst_v)
    call fmm_l2l_bessel_rotation(src_c, src_r, dst_r, kappa, p, vscales, vcnk, one, &
        & src_l, zero, dst_l)
    call fmm_l2p_bessel_baseline(kappa*x, kappa*dst_r, p, vscales, one, dst_l, zero, &
        & dst_v2)
    call fmm_l2l_bessel_rotation(-src_c, dst_r, src_r, kappa, p, vscales, vcnk, one, &
        & dst_l, zero, src_l)
    call fmm_l2p_bessel_baseline(kappa*(x-src_c), kappa*src_r, p, vscales, one, src_l, &
        & zero, dst_v3)
    write(*, *) dst_v, dst_v2, dst_v3
end subroutine check_l2l_bessel

! Check adjoint l2l
subroutine check_l2l_adj(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    integer, parameter :: nrand=10
    real(dp) :: y(3, nx), r, dst_r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), vfact(2*p+1), &
        & src_l((p+1)**2, nrand), src_l2((p+1)**2, nrand), &
        & dst_l((p+1)**2, nrand), err, tmp(nrand, nrand)
    logical :: ok
    integer :: i, j, iseed(4), ndst_l
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = alpha * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    vfact(1) = one
    do i = 2, 2*p+1
        vfact(i) = vfact(i-1) * sqrt(dble(i-1))
    end do
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_l randomly
    ndst_l = nrand * ((p+1)**2)
    call dlarnv(3, iseed, ndst_l, src_l)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check adjoint L2L"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    ! Check against direct l2l
    do i = 1, nx
        c = y(:, i)
        dst_r = r + dnrm2(3, y(:, i), 1)
        do j = 1, nrand
            call fmm_l2l_rotation(c, r, dst_r, p, vscales, vfact, one, &
                & src_l(:, j), zero, dst_l(:, j))
            call fmm_l2l_rotation_adj(-c, dst_r, r, p, vscales, vfact, one, &
                & dst_l(:, j), zero, src_l2(:, j))
        end do
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, one, dst_l, (p+1)**2, &
            & dst_l, (p+1)**2, zero, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1)
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, -one, src_l2, (p+1)**2, &
            & src_l, (p+1)**2, one, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1) / err
        ok = err .lt. 1d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_l2l_adj

! Check adjoint l2l Bessel
subroutine check_l2l_bessel_adj(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    integer, parameter :: nrand=10
    real(dp) :: y(3, nx), r, dst_r, c(3), vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2), vfact(2*p+1), &
        & src_l((p+1)**2, nrand), src_l2((p+1)**2, nrand), &
        & dst_l((p+1)**2, nrand), err, tmp(nrand, nrand), kappa
    logical :: ok
    integer :: i, j, iseed(4), ndst_l
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = alpha * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    vfact(1) = one
    do i = 2, 2*p+1
        vfact(i) = vfact(i-1) * sqrt(dble(i-1))
    end do
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    kappa = 1d+1
    ! Generate src_l randomly
    ndst_l = nrand * ((p+1)**2)
    call dlarnv(3, iseed, ndst_l, src_l)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check adjoint L2L Bessel"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    ! Check against direct l2l
    ! Ignore i=1 for now
    do i = 2, nx
        c = y(:, i)
        dst_r = r + dnrm2(3, y(:, i), 1)
        do j = 1, nrand
            call fmm_l2l_bessel_rotation(c, r, dst_r, kappa, p, vscales, vfact, one, &
                & src_l(:, j), zero, dst_l(:, j))
            call fmm_l2l_bessel_rotation_adj(-c, dst_r, r, kappa, p, vscales, vfact, one, &
                & dst_l(:, j), zero, src_l2(:, j))
        end do
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, one, dst_l, (p+1)**2, &
            & dst_l, (p+1)**2, zero, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1)
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, -one, src_l2, (p+1)**2, &
            & src_l, (p+1)**2, one, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1) / err
        ok = err .lt. 1d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_l2l_bessel_adj

! Check M2L
subroutine check_m2l(pm, pl, alpha)
    ! Inputs
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: alpha
    ! Local variables
    integer, parameter :: nrand=10
    real(dp) :: y(3, nx), r, dst_r, c(3), vscales((pm+pl+1)**2), &
        & v4pi2lp1(pm+pl+1), &
        & vscales_rel((pm+pl+1)**2), vcnk((2*(pm+pl)+1)*(pm+pl+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), &
        & src_m((pm+1)**2, nrand), dst_l((pl+1)**2, nrand), &
        & dst_l2((pl+1)**2, nrand), err
    logical :: ok
    integer :: i, j, iseed(4), nsrc_m, ndst_l
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = 3d0 * alpha / sqrt(13d-2) * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(pm+pl, vscales, v4pi2lp1, vscales_rel)
    call fmm_constants(pm+pl, pm, pl, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_m randomly
    nsrc_m = nrand * ((pm+1)**2)
    ndst_l = nrand * ((pl+1)**2)
    call dlarnv(3, iseed, nsrc_m, src_m)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check M2L"
    print "(A,I0)", " pm=", pm
    print "(A,I0)", " pl=", pl
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    ! Check against the baseline, i=1 is ignored since it y(:,1)=zero
    do i = 2, nx
        c = y(:, i)
        !dst_r = dnrm2(3, c, 1) - two*r
        dst_r = r
        ! Check alpha=beta=zero
        do j = 1, nrand
            call fmm_m2l_rotation(c, r, dst_r, pm, pl, vscales, &
                & m2l_ztranslate_coef, zero, src_m(:, j), zero, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, dst_l2, 1)
        ok = err .eq. zero
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha=one, beta=zero
        dst_l = zero
        do j = 1, nrand
            call fmm_m2l_baseline(c, r, dst_r, pm, pl, vscales, src_m(:, j), &
                & dst_l(:, j))
            call fmm_m2l_rotation(c, r, dst_r, pm, pl, vscales, &
                & m2l_ztranslate_coef, one, src_m(:, j), zero, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, dst_l-dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 1d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        !if(.not. ok) stop 1
        if (.not. ok) then
            print *, src_m
            print *, dst_l2
            stop 1
        end if
        ! Check alpha!={zero,one}, beta=zero
        do j = 1, nrand
            call fmm_m2l_rotation(c, r, dst_r, pm, pl, vscales, &
                & m2l_ztranslate_coef, -three, src_m(:, j), zero, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, three*dst_l+dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 1d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha=zero, beta=one
        dst_l2 = dst_l
        do j = 1, nrand
            call fmm_m2l_rotation(c, r, dst_r, pm, pl, vscales, &
                & m2l_ztranslate_coef, zero, src_m(:, j), one, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, dst_l-dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 1d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha=zero, beta!={zero,one}
        dst_l2 = dst_l
        do j = 1, nrand
            call fmm_m2l_rotation(c, r, dst_r, pm, pl, vscales, &
                & m2l_ztranslate_coef, zero, src_m(:, j), -pt5, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, pt5*dst_l+dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 1d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
        ! Check alpha!={zero,one}, beta!={zero,one}
        dst_l2 = two * dst_l
        do j = 1, nrand
            call fmm_m2l_rotation(c, r, dst_r, pm, pl, vscales, &
                & m2l_ztranslate_coef, -three, src_m(:, j), pt5, dst_l2(:, j))
        end do
        err = dnrm2(ndst_l, two*dst_l+dst_l2, 1) / dnrm2(ndst_l, dst_l, 1)
        ok = err .le. 1d-14
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if(.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2l

! Check adjoint m2l
subroutine check_m2l_adj(pm, pl, alpha)
    ! Inputs
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: alpha
    ! Local variables
    integer, parameter :: nrand=10
    real(dp) :: y(3, nx), r, dst_r, c(3), vscales((pm+pl+1)**2), &
        & v4pi2lp1(pm+pl+1), &
        & vscales_rel((pm+pl+1)**2), vcnk((2*(pm+pl)+1)*(pm+pl+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), &
        & src_m((pm+1)**2, nrand), src_m2((pm+1)**2, nrand), &
        & dst_l((pl+1)**2, nrand), err, tmp(nrand, nrand), threshold
    logical :: ok
    integer :: i, j, iseed(4), nsrc_m, ndst_l
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = 3d0 * alpha / sqrt(13d-2) * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(pm+pl, vscales, v4pi2lp1, vscales_rel)
    call fmm_constants(pm+pl, pm, pl, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    ! Generate src_m randomly
    nsrc_m = nrand * ((pm+1)**2)
    ndst_l = nrand * ((pl+1)**2)
    call dlarnv(3, iseed, nsrc_m, src_m)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check adjoint M2L"
    print "(A,I0)", " pm=", pm
    print "(A,I0)", " pl=", pl
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = (pm+pl+3) * 4d-16
    ! Check against the baseline, i=1 is ignored since it y(:,1)=zero
    do i = 2, nx
        c = y(:, i)
        !dst_r = r + dnrm2(3, y(:, i), 1)
        dst_r = r
        do j = 1, nrand
            call fmm_m2l_rotation(c, r, dst_r, pm, pl, vscales, &
                & m2l_ztranslate_coef, one, src_m(:, j), zero, dst_l(:, j))
            call fmm_m2l_rotation_adj(-c, dst_r, r, pl, pm, vscales, &
                & m2l_ztranslate_adj_coef, one, dst_l(:, j), zero, &
                & src_m2(:, j))
        end do
        call dgemm('T', 'N', nrand, nrand, (pl+1)**2, one, dst_l, (pl+1)**2, &
            & dst_l, (pl+1)**2, zero, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1)
        call dgemm('T', 'N', nrand, nrand, (pm+1)**2, -one, src_m2, &
            & (pm+1)**2, src_m, (pm+1)**2, one, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1) / err
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2l_adj

subroutine check_m2l_bessel(p)
    ! Inputs
    integer, intent(in) :: p
    ! Local variables
    real(dp) :: x(3), src_c(3), src_r, dst_r, dst_v, dst_v2
    real(dp) :: src_m((p+1)**2), dst_l((p+1)**2), vscales((p+1)**2), &
        & vscales_rel((p+1)**2), v4pi2lp1(p+1), vcnk((2*p+1)*(p+1)), &
        tmp(p+1), work(p+1), kappa
    complex(dp) :: work_complex(max(2, p+1))
    real(dp), external :: dnrm2
    integer :: i
    ! Compute special FMM constants
    call ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    src_c = 3d0 * (/0d0, 1d0, 0d0/)
    src_r = one
    x = 5d-1 * (/one, -1.1d0, one/)
    !x = zero
    dst_r = one
    src_m = zero
    do i = 0, p
        !src_m(i*i+i+1) = 10d0 ** (-i)
        !src_m(i*i+1:i*i+2*i+1) = 10d0 ** (-i)
    end do
    src_m = one
    kappa = 1d0
    call fmm_m2p_bessel_baseline(kappa*(x-src_c), kappa*src_r, p, vscales, one, src_m, &
        & zero, dst_v)
    call fmm_m2l_bessel_rotation(src_c, src_r, dst_r, kappa, p, vscales, vcnk, one, &
        & src_m, zero, dst_l)
    call fmm_l2p_bessel_baseline(kappa*x, kappa*dst_r, p, vscales, one, dst_l, zero, &
        & dst_v2)
    write(*, *) dst_v, dst_v2
end subroutine check_m2l_bessel

! Check adjoint m2l
subroutine check_m2l_bessel_adj(p, alpha)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: alpha
    ! Local variables
    integer, parameter :: nrand=10
    real(dp) :: y(3, nx), r, dst_r, c(3), vscales((2*p+1)**2), &
        & v4pi2lp1(2*p+1), &
        & vscales_rel((2*p+1)**2), vcnk((4*p+1)*(2*p+1)), &
        & m2l_ztranslate_coef(p+1, p+1, p+1), &
        & m2l_ztranslate_adj_coef(p+1, p+1, p+1), &
        & src_m((p+1)**2, nrand), src_m2((p+1)**2, nrand), &
        & dst_l((p+1)**2, nrand), err, tmp(nrand, nrand), threshold, kappa
    logical :: ok
    integer :: i, j, iseed(4), nsrc_m, ndst_l
    real(dp), external :: dnrm2
    ! Copy templated points with a proper multiplier
    y = 3d0 * alpha / sqrt(13d-2) * x
    r = abs(alpha)
    ! Compute special FMM constants
    call ylmscale(2*p, vscales, v4pi2lp1, vscales_rel)
    call fmm_constants(2*p, p, p, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    ! Init random seed
    iseed = (/0, 0, 0, 1/)
    kappa = 1d+1
    ! Generate src_m randomly
    nsrc_m = nrand * ((p+1)**2)
    ndst_l = nrand * ((p+1)**2)
    call dlarnv(3, iseed, nsrc_m, src_m)
    ! Print header
    print "(/,A)", repeat("=", 40)
    print "(A)", "Check adjoint M2L Bessel"
    print "(A,I0)", " p=", p
    print "(A,ES24.16E3)", " alpha=", alpha
    print "(A)", repeat("=", 40)
    print "(A)", "  i | ok | err(i)"
    print "(A)", repeat("=", 40)
    threshold = (2*p+3) * 4d-16
    ! Check against the baseline, i=1 is ignored since it y(:,1)=zero
    do i = 2, nx
        c = y(:, i)
        !dst_r = r + dnrm2(3, y(:, i), 1)
        dst_r = r
        do j = 1, nrand
            call fmm_m2l_bessel_rotation(c, r, dst_r, kappa, p, vscales, &
                & m2l_ztranslate_coef, one, src_m(:, j), zero, dst_l(:, j))
            call fmm_m2l_bessel_rotation_adj(-c, dst_r, r, kappa, p, vscales, &
                & m2l_ztranslate_adj_coef, one, dst_l(:, j), zero, &
                & src_m2(:, j))
        end do
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, one, dst_l, (p+1)**2, &
            & dst_l, (p+1)**2, zero, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1)
        call dgemm('T', 'N', nrand, nrand, (p+1)**2, -one, src_m2, &
            & (p+1)**2, src_m, (p+1)**2, one, tmp, nrand)
        err = dnrm2(nrand*nrand, tmp, 1) / err
        ok = err .lt. threshold
        print "(I3.2,A,L3,A,D12.4)", i, " |", ok, " | ", err
        if (.not. ok) stop 1
    end do
    print "(A)", repeat("=", 40)
end subroutine check_m2l_bessel_adj

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

!subroutine check_tree_m2m(p, alpha)
!    ! Inputs
!    integer, intent(in) :: p
!    real(dp), intent(in) :: alpha
!    ! Local variables
!    integer, parameter :: nsph = 10
!    integer :: i, indi, j, k, order(nsph), istatus
!    real(dp) :: csph(3, nsph), rsph(nsph), src(3, nsph), csph2(3, nsph), &
!        & rsph2(nsph), sph_m((p+1)**2, nsph), node_m((p+1)**2, 2*nsph-1), &
!        & node_m2((p+1)**2, 2*nsph-1), vscales((p+1)**2), vfact(2*p+1), &
!        & rel_err, sph_m2((p+1)**2, nsph), full_norm, diff_norm, &
!        & node_m3((p+1)**2, 2*nsph-1)
!    type(ddx_type) :: ddx_data
!    logical :: ok
!    real(dp), external :: dnrm2
!    real(dp) :: full_mat((p+1)**2, 2*nsph-1, (p+1)**2, nsph)
!    ! Scale inputs
!    csph(:, 1) = alpha * (/1d0, 1d0, 1d0/)
!    csph(:, 2) = alpha * (/2d0, 2d0, 2d0/)
!    csph(:, 3) = alpha * (/1d0, 1d0, 3d0/)
!    csph(:, 4) = alpha * (/2d0, 2d0, 4d0/)
!    csph(:, 5) = alpha * (/1d0, 1d0, 5d0/)
!    csph(:, 6) = alpha * (/2d0, 2d0, 6d0/)
!    csph(:, 7) = alpha * (/1d0, 1d0, 7d0/)
!    csph(:, 8) = alpha * (/2d0, 2d0, 8d0/)
!    csph(:, 9) = alpha * (/1d0, 1d0, 9d0/)
!    csph(:, 10) = alpha * (/2d0, 2d0, 10d0/)
!    rsph = abs(alpha) * (/2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, &
!        & 1d0, 1.1d0/)
!    src(:, 1) = alpha * (/1.1d0, 0.9d0, 1d0/)
!    src(:, 2) = alpha * (/2.1d0, 1.9d0, 2.1d0/)
!    src(:, 3) = alpha * (/1.1d0, 1.1d0, 3.1d0/)
!    src(:, 4) = alpha * (/1.9d0, 1.9d0, 4.1d0/)
!    src(:, 5) = alpha * (/1.1d0, 0.9d0, 5d0/)
!    src(:, 6) = alpha * (/1.9d0, 2.1d0, 5.9d0/)
!    src(:, 7) = alpha * (/1.1d0, 0.9d0, 7d0/)
!    src(:, 8) = alpha * (/2.1d0, 2.1d0, 8d0/)
!    src(:, 9) = alpha * (/1.1d0, 0.9d0, 9d0/)
!    src(:, 10) = alpha * (/2.1d0, 1.9d0, 9.9d0/)
!    ! Reorder inputs
!    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
!    do i = 1, nsph
!        csph2(:, i) = csph(:, order(i))
!        rsph2(i) = rsph(order(i))
!    end do
!    ddx_data % nclusters = 2*nsph-1
!    ddx_data % pl = -1
!    ddx_data % pm = p
!    ddx_data % nsph = nsph
!    ! Allocate space for a tree
!    allocate(ddx_data % cluster(2, 2*nsph-1), ddx_data % children(2, 2*nsph-1), &
!        & ddx_data % parent(2*nsph-1), ddx_data % cnode(3, 2*nsph-1), &
!        & ddx_data % rnode(2*nsph-1), ddx_data % snode(nsph), &
!        & ddx_data % order(nsph), ddx_data % vscales((p+1)**2), &
!        & ddx_data % vfact(2*p+1), ddx_data % vscales_rel((p+1)**2), &
!        & ddx_data % v4pi2lp1(p+1), ddx_data % vcnk((2*p+1)*(p+1)), &
!        & ddx_data % m2l_ztranslate_coef(p+1, 0, 0), &
!        & ddx_data % m2l_ztranslate_adj_coef(0, 0, p+1), stat=istatus)
!    if(istatus .ne. 0) stop "ALLOC FAILED"
!    ! Build a recursive inertial binary tree
!    call tree_rib_build(nsph, csph2, rsph2, ddx_data % order, &
!        & ddx_data % cluster, ddx_data % children, ddx_data % parent, &
!        & ddx_data % cnode, ddx_data % rnode, ddx_data % snode)
!    ! Get constants, corresponding to given maximum degree of spherical
!    ! harmonics
!    call ylmscale(p, ddx_data % vscales, ddx_data % v4pi2lp1, &
!        & ddx_data % vscales_rel)
!    ddx_data % vfact(1) = one
!    do i = 2, 2*p+1
!        ddx_data % vfact(i) = ddx_data % vfact(i-1) * sqrt(dble(i-1))
!    end do
!    call fmm_constants(p, p, -1, ddx_data % vcnk, &
!        & ddx_data % m2l_ztranslate_coef, &
!        & ddx_data % m2l_ztranslate_adj_coef)
!    ! Init input harmonics
!    do i = 1, nsph
!        call fmm_p2m(src(:, i)-csph(:, i), one, rsph(i), p, &
!            & ddx_data % vscales, zero, sph_m(:, i))
!    end do
!    ! Get reference result of M2M operation
!    do i = 1, 2*nsph-1
!        node_m(:, i) = zero
!        do j = ddx_data % cluster(1, i), ddx_data % cluster(2, i)
!            k = ddx_data % order(j)
!            call fmm_m2m_rotation(csph2(:, k)-ddx_data % cnode(:, i), &
!                & rsph2(k), ddx_data % rnode(i), p, ddx_data % vscales, &
!                & ddx_data % vcnk, one, sph_m(:, k), one, node_m(:, i))
!        end do
!    end do
!    ! Check tree_m2m_rotation
!    if (iprint .gt. 0) then
!        write(*, "(/,A)") repeat("=", 40)
!        write(*, "(A)") "Check tree_m2m_rotation"
!        write(*, "(4x,A,I0)") "p = ", p
!        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!        write(*, "(4x,A,A)") "err = || [tree M2M] - [plain M2M] ||", &
!            & " / || [plain M2M] ||"
!        write(*, "(A)") repeat("=", 40)
!        write(*, "(A)") " ok | err"
!        write(*, "(A)") repeat("=", 40)
!    end if
!    node_m2 = one
!    do i = 1, nsph
!        node_m2(:, ddx_data % snode(i)) = sph_m(:, i)
!    end do
!    call tree_m2m_rotation(ddx_data, node_m2)
!    rel_err = dnrm2(((p+1)**2)*(2*nsph-1), node_m-node_m2, 1) / &
!        & dnrm2(((p+1)**2)*(2*nsph-1), node_m, 1)
!    ok = rel_err .le. 1d-15
!    write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!    write(*, "(A,/)") repeat("=", 40)
!    if (.not. ok) stop 1
!    ! Allocate space for transfer matrices
!    ddx_data % m2m_reflect_mat_size = (p+1)*(2*p+1)*(2*p+3)/3
!    ddx_data % m2m_ztranslate_mat_size = (p+1)*(p+2)*(p+3)/6
!    allocate(ddx_data % m2m_reflect_mat(ddx_data % m2m_reflect_mat_size, &
!        & ddx_data % nclusters-1), stat=istatus)
!    if (istatus .ne. 0) stop "Allocation failed"
!    allocate(ddx_data % m2m_ztranslate_mat(ddx_data % m2m_ztranslate_mat_size, &
!        & ddx_data % nclusters-1), stat=istatus)
!    if (istatus .ne. 0) stop "Allocation failed"
!    ! Compute transfer matrices
!    call tree_m2m_reflection_get_mat(ddx_data)
!    ! Check tree_m2m_reflection_use_mat
!    write(*, "(/,A)") repeat("=", 40)
!    write(*, "(A)") "Check tree_m2m_reflection_get_mat"
!    write(*, "(A)") "Check tree_m2m_reflection_use_mat"
!    write(*, "(4x,A,I0)") "p = ", p
!    write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!    write(*, "(4x,A,A)") "err = || [tree M2M] - [plain M2M] ||", &
!        & " / || [plain M2M] ||"
!    write(*, "(A)") repeat("=", 40)
!    write(*, "(A)") " ok | err"
!    write(*, "(A)") repeat("=", 40)
!    node_m2 = one
!    do i = 1, nsph
!        node_m2(:, ddx_data % snode(i)) = sph_m(:, i)
!    end do
!    call tree_m2m_reflection_use_mat(ddx_data, node_m2)
!    rel_err = dnrm2(((p+1)**2)*(2*nsph-1), node_m-node_m2, 1) / &
!        & dnrm2(((p+1)**2)*(2*nsph-1), node_m, 1)
!    ok = rel_err .le. 1d-15
!    write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!    write(*, "(A,/)") repeat("=", 40)
!    if (.not. ok) stop 1
!    ! Check tree_m2m_reflection_use_mat_adj
!    write(*, "(/,A)") repeat("=", 40)
!    write(*, "(A)") "Check tree_m2m_reflection_use_mat_adj"
!    write(*, "(4x,A,I0)") "p = ", p
!    write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!    write(*, "(4x,A,A)") "err = || [tree M2M] - [plain M2M] ||", &
!        & " / || [plain M2M] ||"
!    write(*, "(A)") repeat("=", 40)
!    write(*, "(A)") " ok | err"
!    write(*, "(A)") repeat("=", 40)
!    do i = 1, nsph
!        do j = 1, (p+1)**2
!            node_m2 = zero
!            node_m2(j, ddx_data % snode(i)) = one
!            call tree_m2m_reflection_use_mat(ddx_data, node_m2)
!            full_mat(:, :, j, i) = node_m2
!        end do
!    end do
!    full_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
!    do i = 1, 2*nsph-1
!        do j = 1, (p+1)**2
!            node_m2 = zero
!            node_m2(j, i) = one
!            call tree_m2m_reflection_use_mat_adj(ddx_data, node_m2)
!            do k = 1, nsph
!                full_mat(j, i, :, k) = full_mat(j, i, :, k) - &
!                    & node_m2(:, ddx_data % snode(k))
!            end do
!        end do
!    end do
!    diff_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
!    rel_err = diff_norm / full_norm
!    ok = rel_err .le. 1d-15
!    write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!    write(*, "(A,/)") repeat("=", 40)
!    if (.not. ok) stop 1
!    ! Check tree_m2m_rotation_adj
!    write(*, "(/,A)") repeat("=", 40)
!    write(*, "(A)") "Check tree_m2m_rotation_adj"
!    write(*, "(4x,A,I0)") "p = ", p
!    write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!    write(*, "(4x,A,A)") "err = || [tree M2M] - [plain M2M] ||", &
!        & " / || [plain M2M] ||"
!    write(*, "(A)") repeat("=", 40)
!    write(*, "(A)") " ok | err"
!    write(*, "(A)") repeat("=", 40)
!    node_m2 = node_m
!    call tree_m2m_reflection_use_mat_adj(ddx_data, node_m2)
!    node_m3 = node_m
!    call tree_m2m_rotation_adj(ddx_data, node_m3)
!    rel_err = dnrm2((2*nsph-1)*((p+1)**2), node_m2-node_m3, 1) / &
!        dnrm2((2*nsph-1)*((p+1)**2), node_m2, 1)
!    ok = rel_err .le. 1d-15
!    write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!    write(*, "(A,/)") repeat("=", 40)
!    if (.not. ok) stop 1
!    ! Deallocate tree
!    call ddfree(ddx_data)
!end subroutine check_tree_m2m
!
!subroutine check_tree_l2l(p, alpha)
!    ! Inputs
!    integer, intent(in) :: p, iprint
!    real(dp), intent(in) :: alpha, threshold
!    ! Local variables
!    integer, parameter :: nsph = 10
!    integer :: i, j, k, order(nsph), istatus
!    real(dp) :: csph(3, nsph), rsph(nsph), src(3, nsph), csph2(3, nsph), &
!        & rsph2(nsph), sph_l((p+1)**2, nsph), node_l((p+1)**2, 2*nsph-1), &
!        & node_l2((p+1)**2, 2*nsph-1), vscales((p+1)**2), vfact(2*p+1), &
!        & rel_err, sph_l2((p+1)**2, nsph), full_norm, diff_norm, &
!        & node_l3((p+1)**2, 2*nsph-1)
!    type(ddx_type) :: ddx_data
!    logical :: ok
!    real(dp), external :: dnrm2
!    real(dp) :: full_mat((p+1)**2, nsph, (p+1)**2, 2*nsph-1)
!    ! Scale inputs
!    csph(:, 1) = alpha * (/1d0, 1d0, 1d0/)
!    csph(:, 2) = alpha * (/2d0, 2d0, 2d0/)
!    csph(:, 3) = alpha * (/1d0, 1d0, 3d0/)
!    csph(:, 4) = alpha * (/2d0, 2d0, 4d0/)
!    csph(:, 5) = alpha * (/1d0, 1d0, 5d0/)
!    csph(:, 6) = alpha * (/2d0, 2d0, 6d0/)
!    csph(:, 7) = alpha * (/1d0, 1d0, 7d0/)
!    csph(:, 8) = alpha * (/2d0, 2d0, 8d0/)
!    csph(:, 9) = alpha * (/1d0, 1d0, 9d0/)
!    csph(:, 10) = alpha * (/2d0, 2d0, 10d0/)
!    rsph = abs(alpha) * (/2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, &
!        & 1d0, 1.1d0/)
!    src(:, 1) = alpha * (/1.1d0, 0.9d0, 1d0/)
!    src(:, 2) = alpha * (/2.1d0, 1.9d0, 2.1d0/)
!    src(:, 3) = alpha * (/1.1d0, 1.1d0, 3.1d0/)
!    src(:, 4) = alpha * (/1.9d0, 1.9d0, 4.1d0/)
!    src(:, 5) = alpha * (/1.1d0, 0.9d0, 5d0/)
!    src(:, 6) = alpha * (/1.9d0, 2.1d0, 5.9d0/)
!    src(:, 7) = alpha * (/1.1d0, 0.9d0, 7d0/)
!    src(:, 8) = alpha * (/2.1d0, 2.1d0, 8d0/)
!    src(:, 9) = alpha * (/1.1d0, 0.9d0, 9d0/)
!    src(:, 10) = alpha * (/2.1d0, 1.9d0, 9.9d0/)
!    ! Reorder inputs
!    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
!    do i = 1, nsph
!        csph2(:, i) = csph(:, order(i))
!        rsph2(i) = rsph(order(i))
!    end do
!    ddx_data % nclusters = 2*nsph-1
!    ddx_data % pl = p
!    ddx_data % pm = p
!    ddx_data % nsph = nsph
!    ! Allocate space for a tree
!    allocate(ddx_data % cluster(2, 2*nsph-1), ddx_data % children(2, 2*nsph-1), &
!        & ddx_data % parent(2*nsph-1), ddx_data % cnode(3, 2*nsph-1), &
!        & ddx_data % rnode(2*nsph-1), ddx_data % snode(nsph), &
!        & ddx_data % order(nsph), ddx_data % vscales((p+1)**2), &
!        & ddx_data % vfact(2*p+1), ddx_data % vscales_rel((p+1)**2), &
!        & ddx_data % v4pi2lp1(p+1), ddx_data % vcnk((2*p+1)*(p+1)), &
!        & ddx_data % m2l_ztranslate_coef(0, p+1, p+1), &
!        & ddx_data % m2l_ztranslate_adj_coef(p+1, p+1, 0), stat=istatus)
!    ! Build a recursive inertial binary tree
!    call tree_rib_build(nsph, csph2, rsph2, ddx_data % order, &
!        & ddx_data % cluster, ddx_data % children, ddx_data % parent, &
!        & ddx_data % cnode, ddx_data % rnode, ddx_data % snode)
!    ! Get constants, corresponding to given maximum degree of spherical
!    ! harmonics
!    call ylmscale(p, ddx_data % vscales, ddx_data % v4pi2lp1, &
!        & ddx_data % vscales_rel)
!    ddx_data % vfact(1) = one
!    do i = 2, 2*p+1
!        ddx_data % vfact(i) = ddx_data % vfact(i-1) * sqrt(dble(i-1))
!    end do
!    call fmm_constants(p, -1, p, ddx_data % vcnk, &
!        & ddx_data % m2l_ztranslate_coef, &
!        & ddx_data % m2l_ztranslate_adj_coef)
!    ! Init input harmonics
!    node_l = zero
!    do i = 1, nsph
!        call fmm_p2m(src(:, i)-csph(:, i), one, rsph(i), p, &
!            & ddx_data % vscales, zero, node_l(:, ddx_data % snode(i)))
!    end do
!    call tree_m2m_rotation(ddx_data, node_l)
!    ddx_data % pm = -1
!    ! Get reference result of L2L operation
!    sph_l = zero
!    do i = 1, 2*nsph-1
!        do j = ddx_data % cluster(1, i), ddx_data % cluster(2, i)
!            k = ddx_data % order(j)
!            call fmm_l2l_rotation(ddx_data % cnode(:, i)-csph2(:, k), &
!                & ddx_data % rnode(i), rsph2(k), p, ddx_data % vscales, &
!                & ddx_data % vfact, one, node_l(:, i), one, sph_l(:, k))
!        end do
!    end do
!    ! Check tree_l2l_rotation
!    if (iprint .gt. 0) then
!        write(*, "(/,A)") repeat("=", 40)
!        write(*, "(A)") "Check tree_l2l_rotation"
!        write(*, "(4x,A,I0)") "p = ", p
!        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!        write(*, "(4x,A,A)") "err = || [tree L2L] - [plain L2L] ||", &
!            & " / || [plain L2L] ||"
!        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
!        write(*, "(A)") repeat("=", 40)
!        write(*, "(A)") " ok | err"
!        write(*, "(A)") repeat("=", 40)
!    end if
!    node_l2 = node_l
!    call tree_l2l_rotation(ddx_data, node_l2)
!    do i = 1, nsph
!        sph_l2(:, i) = node_l2(:, ddx_data % snode(i))
!    end do
!    rel_err = dnrm2(((p+1)**2)*nsph, sph_l-sph_l2, 1) / &
!        & dnrm2(((p+1)**2)*nsph, sph_l, 1)
!    ok = rel_err .le. threshold
!    if (iprint .gt. 0) then
!        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!        write(*, "(A,/)") repeat("=", 40)
!    end if
!    if (.not. ok) stop 1
!    ! Allocate space for transfer matrices
!    ddx_data % l2l_reflect_mat_size = (p+1)*(2*p+1)*(2*p+3)/3
!    ddx_data % l2l_ztranslate_mat_size = (p+1)*(p+2)*(p+3)/6
!    allocate(ddx_data % l2l_reflect_mat(ddx_data % l2l_reflect_mat_size, &
!        & ddx_data % nclusters-1), stat=istatus)
!    if (istatus .ne. 0) stop "Allocation failed"
!    allocate(ddx_data % l2l_ztranslate_mat(ddx_data % l2l_ztranslate_mat_size, &
!        & ddx_data % nclusters-1), stat=istatus)
!    if (istatus .ne. 0) stop "Allocation failed"
!    ! Compute transfer matrices
!    call tree_l2l_reflection_get_mat(ddx_data)
!    ! Check tree_l2l_reflection_use_mat
!    if (iprint .gt. 0) then
!        write(*, "(/,A)") repeat("=", 40)
!        write(*, "(A)") "Check tree_l2l_reflection_get_mat"
!        write(*, "(A)") "Check tree_l2l_reflection_use_mat"
!        write(*, "(4x,A,I0)") "p = ", p
!        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!        write(*, "(4x,A,A)") "err = || [tree L2L] - [plain L2L] ||", &
!            & " / || [plain L2L] ||"
!        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
!        write(*, "(A)") repeat("=", 40)
!        write(*, "(A)") " ok | err"
!        write(*, "(A)") repeat("=", 40)
!    end if
!    node_l2 = node_l
!    call tree_l2l_reflection_use_mat(ddx_data, node_l2)
!    do i = 1, nsph
!        sph_l2(:, i) = node_l2(:, ddx_data % snode(i))
!    end do
!    rel_err = dnrm2(((p+1)**2)*nsph, sph_l-sph_l2, 1) / &
!        & dnrm2(((p+1)**2)*nsph, sph_l, 1)
!    ok = rel_err .le. threshold
!    if (iprint .gt. 0) then
!        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!        write(*, "(A,/)") repeat("=", 40)
!    end if
!    if (.not. ok) stop 1
!    ! Check tree_l2l_reflection_use_mat_adj
!    if (iprint .gt. 0) then
!        write(*, "(/,A)") repeat("=", 40)
!        write(*, "(A)") "Check tree_l2l_reflection_use_mat_adj"
!        write(*, "(4x,A,I0)") "p = ", p
!        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!        write(*, "(4x,A,A)") "err = || [tree L2L] - [plain L2L] ||", &
!            & " / || [plain L2L] ||"
!        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
!        write(*, "(A)") repeat("=", 40)
!        write(*, "(A)") " ok | err"
!        write(*, "(A)") repeat("=", 40)
!    end if
!    do i = 1, 2*nsph-1
!        do j = 1, (p+1)**2
!            node_l2 = zero
!            node_l2(j, i) = one
!            call tree_l2l_reflection_use_mat(ddx_data, node_l2)
!            do k = 1, nsph
!                full_mat(:, k, j, i) = node_l2(:, ddx_data % snode(k))
!            end do
!        end do
!    end do
!    full_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
!    do i = 1, nsph
!        do j = 1, (p+1)**2
!            node_l2 = one
!            do k = 1, nsph
!                node_l2(:, ddx_data % snode(k)) = zero
!            end do
!            node_l2(j, ddx_data % snode(i)) = one
!            call tree_l2l_reflection_use_mat_adj(ddx_data, node_l2)
!            full_mat(j, i, :, :) = full_mat(j, i, :, :) - node_l2
!        end do
!    end do
!    diff_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
!    rel_err = diff_norm / full_norm
!    ok = rel_err .le. threshold
!    if (iprint .gt. 0) then
!        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!        write(*, "(A,/)") repeat("=", 40)
!    end if
!    if (.not. ok) stop 1
!    ! Check tree_l2l_rotation_adj
!    if (iprint .gt. 0) then
!        write(*, "(/,A)") repeat("=", 40)
!        write(*, "(A)") "Check tree_l2l_rotation_adj"
!        write(*, "(4x,A,I0)") "p = ", p
!        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!        write(*, "(4x,A,A)") "err = || [tree L2L] - [plain L2L] ||", &
!            & " / || [plain L2L] ||"
!        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
!        write(*, "(A)") repeat("=", 40)
!        write(*, "(A)") " ok | err"
!        write(*, "(A)") repeat("=", 40)
!    end if
!    node_l2 = zero
!    do i = 1, nsph
!        node_l2(:, ddx_data % snode(i)) = node_l(:, ddx_data % snode(i))
!    end do
!    call tree_l2l_reflection_use_mat_adj(ddx_data, node_l2)
!    node_l3 = zero
!    do i = 1, nsph
!        node_l3(:, ddx_data % snode(i)) = node_l(:, ddx_data % snode(i))
!    end do
!    call tree_l2l_rotation_adj(ddx_data, node_l3)
!    rel_err = dnrm2((2*nsph-1)*((p+1)**2), node_l2-node_l3, 1) / &
!        dnrm2((2*nsph-1)*((p+1)**2), node_l2, 1)
!    ok = rel_err .le. threshold
!    if (iprint .gt. 0) then
!        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!        write(*, "(A,/)") repeat("=", 40)
!    end if
!    if (.not. ok) stop 1
!    ! Deallocate tree
!    call ddfree(ddx_data)
!end subroutine check_tree_l2l
!
!subroutine check_tree_m2l(pm, pl, alpha, threshold)
!    ! Inputs
!    integer, intent(in) :: pm, pl, iprint
!    real(dp), intent(in) :: alpha, threshold
!    ! Local variables
!    integer, parameter :: nsph = 10, lwork = 1000
!    integer :: i, j, k, l, order(nsph), istatus, iwork, jwork, work(3, lwork)
!    real(dp) :: csph(3, nsph), rsph(nsph), src(3, nsph), csph2(3, nsph), &
!        & rsph2(nsph), node_m((pm+1)**2, 2*nsph-1), &
!        & node_l((pl+1)**2, 2*nsph-1), vscales((pm+pl+1)**2), &
!        & vfact(2*(pm+pl)+1), rel_err, full_norm, diff_norm, &
!        & node_m2((pm+1)**2, 2*nsph-1), node_l2((pl+1)**2, 2*nsph-1)
!    type(ddx_type) :: ddx_data
!    logical :: ok
!    real(dp), external :: dnrm2
!    real(dp) :: full_mat((pl+1)**2, 2*nsph-1, (pm+1)**2, 2*nsph-1)
!    real(dp) :: far_p2p(nsph), far_p2p2(nsph)
!    ! Scale inputs
!    csph(:, 1) = alpha * (/1d0, 1d0, 1d0/)
!    csph(:, 2) = alpha * (/2d0, 2d0, 2d0/)
!    csph(:, 3) = alpha * (/1d0, 1d0, 3d0/)
!    csph(:, 4) = alpha * (/2d0, 2d0, 4d0/)
!    csph(:, 5) = alpha * (/1d0, 1d0, 5d0/)
!    csph(:, 6) = alpha * (/2d0, 2d0, 6d0/)
!    csph(:, 7) = alpha * (/1d0, 1d0, 7d0/)
!    csph(:, 8) = alpha * (/2d0, 2d0, 8d0/)
!    csph(:, 9) = alpha * (/1d0, 1d0, 9d0/)
!    csph(:, 10) = alpha * (/2d0, 2d0, 10d0/)
!    rsph = abs(alpha) * (/2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, &
!        & 1d0, 1.1d0/)
!    src(:, 1) = alpha * (/1.1d0, 0.9d0, 1d0/)
!    src(:, 2) = alpha * (/2.1d0, 1.9d0, 2.1d0/)
!    src(:, 3) = alpha * (/1.1d0, 1.1d0, 3.1d0/)
!    src(:, 4) = alpha * (/1.9d0, 1.9d0, 4.1d0/)
!    src(:, 5) = alpha * (/1.1d0, 0.9d0, 5d0/)
!    src(:, 6) = alpha * (/1.9d0, 2.1d0, 5.9d0/)
!    src(:, 7) = alpha * (/1.1d0, 0.9d0, 7d0/)
!    src(:, 8) = alpha * (/2.1d0, 2.1d0, 8d0/)
!    src(:, 9) = alpha * (/1.1d0, 0.9d0, 9d0/)
!    src(:, 10) = alpha * (/2.1d0, 1.9d0, 9.9d0/)
!    ! Reorder inputs
!    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
!    do i = 1, nsph
!        csph2(:, i) = csph(:, order(i))
!        rsph2(i) = rsph(order(i))
!    end do
!    ddx_data % nclusters = 2*nsph-1
!    ddx_data % pl = pl
!    ddx_data % pm = pm
!    ddx_data % nsph = nsph
!    ! Allocate space for a tree
!    allocate(ddx_data % cluster(2, 2*nsph-1), ddx_data % children(2, 2*nsph-1), &
!        & ddx_data % parent(2*nsph-1), ddx_data % cnode(3, 2*nsph-1), &
!        & ddx_data % rnode(2*nsph-1), ddx_data % snode(nsph), &
!        & ddx_data % order(nsph), ddx_data % vscales((pm+pl+1)**2), &
!        & ddx_data % vfact(2*(pm+pl)+1), ddx_data % vscales_rel((pm+pl+1)**2), &
!        & ddx_data % v4pi2lp1(pm+pl+1), &
!        & ddx_data % vcnk((2*(pm+pl)+1)*(pm+pl+1)), &
!        & ddx_data % m2l_ztranslate_coef(pm+1, pl+1, pl+1), &
!        & ddx_data % m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), &
!        & stat=istatus)
!    allocate(ddx_data % nfar(ddx_data % nclusters), &
!        & ddx_data % nnear(ddx_data % nclusters))
!    ! Build a recursive inertial binary tree
!    call tree_rib_build(nsph, csph2, rsph2, ddx_data % order, &
!        & ddx_data % cluster, ddx_data % children, ddx_data % parent, &
!        & ddx_data % cnode, ddx_data % rnode, ddx_data % snode)
!    ! Get list of neighbours for M2L operation
!    iwork = 0
!    call tree_get_farnear_work(ddx_data % nclusters, ddx_data % children, &
!        & ddx_data % cnode, ddx_data % rnode, lwork, iwork, jwork, work, &
!        & ddx_data % nnfar, ddx_data % nfar, ddx_data % nnnear, ddx_data % nnear)
!    if (iwork .le. jwork) then
!        write(*, "(A,A)") "Value of lwork, size of temporary buffer, ", &
!            & "is too low, please increase it"
!        stop 1
!    end if
!    allocate(ddx_data % far(ddx_data % nnfar), &
!        & ddx_data % sfar(ddx_data % nclusters+1), &
!        & ddx_data % near(ddx_data % nnnear), &
!        & ddx_data % snear(ddx_data % nclusters+1))
!    call tree_get_farnear(jwork, lwork, work, ddx_data % nclusters, &
!        & ddx_data % nnfar, ddx_data % nfar, ddx_data % sfar, ddx_data % far, &
!        & ddx_data % nnnear, ddx_data % nnear, ddx_data % snear, ddx_data % near)
!    ! Get far-field P2P for a reference result
!    do i = 1, nsph
!        far_p2p(i) = zero
!        do j = 1, nsph
!            ok = .true.
!            do k = ddx_data % snear(ddx_data % snode(i)), &
!                & ddx_data % snear(ddx_data % snode(i)+1)-1
!                if (ddx_data % near(k) .eq. ddx_data % snode(j)) ok = .false.
!            end do
!            if (ok) then
!                far_p2p(i) = far_p2p(i) + &
!                    & one/dnrm2(3, src(:, order(i))-src(:, order(j)), 1)
!            end if
!        end do
!    end do
!    ! Get constants, corresponding to given maximum degree of spherical
!    ! harmonics
!    call ylmscale(pm+pl, ddx_data % vscales, ddx_data % v4pi2lp1, &
!        & ddx_data % vscales_rel)
!    ddx_data % vfact(1) = one
!    do i = 2, 2*(pm+pl)+1
!        ddx_data % vfact(i) = ddx_data % vfact(i-1) * sqrt(dble(i-1))
!    end do
!    call fmm_constants(pm+pl, pm, pl, ddx_data % vcnk, &
!        & ddx_data % m2l_ztranslate_coef, &
!        & ddx_data % m2l_ztranslate_adj_coef)
!    ! Init input harmonics
!    do i = 1, nsph
!        call fmm_p2m(src(:, order(i))-csph2(:, i), one, rsph2(i), pm, &
!            & ddx_data % vscales, zero, node_m(:, ddx_data % snode(i)))
!    end do
!    ! Prepare M2M
!    call tree_m2m_rotation(ddx_data, node_m)
!    ! Check tree_m2l_rotation
!    if (iprint .gt. 0) then
!        write(*, "(/,A)") repeat("=", 40)
!        write(*, "(A)") "Check tree_m2l_rotation"
!        write(*, "(4x,A,I0)") "pm = ", pm
!        write(*, "(4x,A,I0)") "pl = ", pl
!        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!        write(*, "(4x,A,A)") "err = || P2M + [tree M2L] + L2P - P2P ||", &
!            & " / || P2P ||"
!        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
!        write(*, "(A)") repeat("=", 40)
!        write(*, "(A)") " ok | err"
!        write(*, "(A)") repeat("=", 40)
!    end if
!    call tree_m2l_rotation(ddx_data, node_m, node_l)
!    call tree_l2l_rotation(ddx_data, node_l)
!    do i = 1, nsph
!        call fmm_l2p(src(:, order(i))-csph2(:, i), rsph2(i), pl, &
!            & ddx_data % vscales_rel, one, node_l(:, ddx_data % snode(i)), &
!            & zero, far_p2p2(i))
!    end do
!    rel_err = dnrm2(nsph, far_p2p-far_p2p2, 1) / dnrm2(nsph, far_p2p, 1)
!    ok = rel_err .le. threshold
!    if (iprint .gt. 0) then
!        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!        write(*, "(A,/)") repeat("=", 40)
!    end if
!    if (.not. ok) stop 1
!    ! Reference value of M2L
!    call tree_m2l_rotation(ddx_data, node_m, node_l)
!    ! Allocate space for transfer matrices
!    ddx_data % m2l_reflect_mat_size = (max(pm,pl)+1)*(2*max(pm,pl)+1)* &
!        & (2*max(pm,pl)+3)/3
!    ddx_data % m2l_ztranslate_mat_size = (min(pm,pl)+1)*(min(pm,pl)+2)* &
!        & (3*max(pm,pl)+3-min(pm,pl))/6
!    allocate(ddx_data % m2l_reflect_mat(ddx_data % m2l_reflect_mat_size, &
!        & ddx_data % nnfar), stat=istatus)
!    if (istatus .ne. 0) stop "Allocation failed"
!    allocate(ddx_data % m2l_ztranslate_mat(ddx_data % m2l_ztranslate_mat_size, &
!        & ddx_data % nnfar), stat=istatus)
!    if (istatus .ne. 0) stop "Allocation failed"
!    ! Check tree_m2l_reflection_use_mat
!    if (iprint .gt. 0) then
!        write(*, "(/,A)") repeat("=", 40)
!        write(*, "(A)") "Check tree_m2l_reflection_get_mat"
!        write(*, "(A)") "Check tree_m2l_reflection_use_mat"
!        write(*, "(4x,A,I0)") "pm = ", pm
!        write(*, "(4x,A,I0)") "pl = ", pl
!        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!        write(*, "(4x,A,A)") "err = || [tree M2L] - [tree_m2l_rotation] ||", &
!            & " / || [tree_m2l_rotation] ||"
!        write(*, "(4x,A,ES23.16E3)") "threshold = ", 1d-15
!        write(*, "(A)") repeat("=", 40)
!        write(*, "(A)") " ok | err"
!        write(*, "(A)") repeat("=", 40)
!    end if
!    call tree_m2l_reflection_get_mat(ddx_data)
!    call tree_m2l_reflection_use_mat(ddx_data, node_m, node_l2)
!    rel_err = dnrm2((2*nsph-1)*((pl+1)**2), node_l2-node_l, 1) / &
!        & dnrm2((2*nsph-1)*((pl+1)**2), node_l, 1)
!    ok = rel_err .le. 1d-15
!    if (iprint .gt. 0) then
!        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!        write(*, "(A,/)") repeat("=", 40)
!    end if
!    if (.not. ok) stop 1
!    ! Check tree_m2l_reflection_use_mat_adj
!    if (iprint .gt. 0) then
!        write(*, "(/,A)") repeat("=", 40)
!        write(*, "(A)") "Check tree_m2l_reflection_use_mat_adj"
!        write(*, "(4x,A,I0)") "pm = ", pm
!        write(*, "(4x,A,I0)") "pl = ", pl
!        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!        write(*, "(4x,A,A)") "err = || [tree M2L] - [adj tree M2L]^T ||", &
!            & " / || [tree M2L] ||"
!        write(*, "(4x,A,ES23.16E3)") "threshold = ", 1d-15
!        write(*, "(A)") repeat("=", 40)
!        write(*, "(A)") " ok | err"
!        write(*, "(A)") repeat("=", 40)
!    end if
!    do i = 1, ddx_data % nclusters
!        do j = 1, (pm+1)**2
!            node_m2 = zero
!            node_m2(j, i) = one
!            call tree_m2l_reflection_use_mat(ddx_data, node_m2, node_l2)
!            full_mat(:, :, j, i) = node_l2
!        end do
!    end do
!    full_norm = dnrm2((((2*nsph-1)*((pm+1)*(pl+1)))**2), full_mat, 1)
!    do i = 1, ddx_data % nclusters
!        do j = 1, (pl+1)**2
!            node_l2 = zero
!            node_l2(j, i) = one
!            call tree_m2l_reflection_use_mat_adj(ddx_data, node_l2, node_m2)
!            full_mat(j, i, :, :) = full_mat(j, i, :, :) - node_m2
!        end do
!    end do
!    diff_norm = dnrm2((((2*nsph-1)*((pm+1)*(pl+1)))**2), full_mat, 1)
!    rel_err = diff_norm / full_norm
!    ok = rel_err .le. 1d-15
!    if (iprint .gt. 0) then
!        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!        write(*, "(A,/)") repeat("=", 40)
!    end if
!    if (.not. ok) stop 1
!    ! Check tree_m2l_rotation_adj
!    if (iprint .gt. 0) then
!        write(*, "(/,A)") repeat("=", 40)
!        write(*, "(A)") "Check tree_m2l_rotation_adj"
!        write(*, "(4x,A,I0)") "pm = ", pm
!        write(*, "(4x,A,I0)") "pl = ", pl
!        write(*, "(4x,A,ES24.16E3)") "alpha = ", alpha
!        write(*, "(4x,A,A)") "err = || [adj tree M2L] - [adj tree M2L] ||", &
!            & " / || [adj tree M2L] ||"
!        write(*, "(4x,A,ES23.16E3)") "threshold = ", 1d-15
!        write(*, "(A)") repeat("=", 40)
!        write(*, "(A)") " ok | err"
!        write(*, "(A)") repeat("=", 40)
!    end if
!    call tree_m2l_reflection_use_mat_adj(ddx_data, node_l, node_m)
!    call tree_m2l_rotation_adj(ddx_data, node_l, node_m2)
!    rel_err = dnrm2((2*nsph-1)*((pm+1)**2), node_m2-node_m, 1) / &
!        & dnrm2((2*nsph-1)*((pm+1)**2), node_m, 1)
!    ok = rel_err .le. 1d-15
!    if (iprint .gt. 0) then
!        write(*, "(L3,A,ES23.16E3)") ok, " | ", rel_err
!        write(*, "(A,/)") repeat("=", 40)
!    end if
!    if (.not. ok) stop 1
!    ! Deallocate tree
!    call ddfree(ddx_data)
!end subroutine check_tree_m2l
!
!subroutine check_tree_l2p(p, alpha, threshold)
!    ! Inputs
!    integer, intent(in) :: p, iprint
!    real(dp), intent(in) :: alpha, threshold
!    ! Local variables
!end subroutine check_tree_l2p
!
!subroutine check_tree_m2p
!end subroutine check_tree_m2p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Baseline implementations of certain routines

subroutine fmm_p2m_baseline(c, src_q, dst_r, p, vscales, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), t, rcoef
    integer :: n, ind
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Harmonics are available only if rho > 0
    if (rho .ne. zero) then
        rcoef = rho / dst_r
        t = src_q / dst_r
        ! Ignore input `m` in case of zero scaling factor
        if (beta .eq. zero) then
            do n = 0, p
                ind = n*n + n + 1
                dst_m(ind-n:ind+n) = t * vylm(ind-n:ind+n)
                t = t * rcoef
            end do
        ! Update `m` otherwise
        else
            do n = 0, p
                ind = n*n + n + 1
                dst_m(ind-n:ind+n) = beta*dst_m(ind-n:ind+n) + &
                    & t*vylm(ind-n:ind+n)
                t = t * rcoef
            end do
        end if
    ! Naive case of rho = 0
    else
        ! Ignore input `m` in case of zero scaling factor
        if (beta .eq. zero) then
            dst_m(1) = src_q / dst_r / sqrt4pi
            dst_m(2:) = zero
        ! Update `m` otherwise
        else
            dst_m(1) = beta*dst_m(1) + src_q/dst_r/sqrt4pi
            dst_m(2:) = beta * dst_m(2:)
        end if
    end if
end subroutine fmm_p2m_baseline

subroutine fmm_p2m_bessel_baseline(c, src_q, dst_r, p, vscales, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1), t
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), si(p+1), dst_sk(p+1), work(p+1)
    real(dp), parameter :: eight=8d0
    complex(dp) :: work_complex(max(2, p+1))
    integer :: n, ind
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Get Bessel function for the sphere
    call modified_spherical_bessel_first_kind(p, dst_r, dst_sk, work, &
        & work_complex)
    ! Harmonics are available only if rho > 0
    if (rho .ne. zero) then
        call modified_spherical_bessel_second_kind(p, rho, si, work, &
            & work_complex)
        ! Ignore input `m` in case of zero scaling factor
        if (beta .eq. zero) then
            do n = 0, p
                ind = n*n + n + 1
                t = eight * src_q * si(n+1) * dst_sk(n+1)
                dst_m(ind-n:ind+n) = t * vylm(ind-n:ind+n)
            end do
        ! Update `m` otherwise
        else
            do n = 0, p
                ind = n*n + n + 1
                t = eight * src_q * si(n+1) * dst_sk(n+1)
                dst_m(ind-n:ind+n) = beta*dst_m(ind-n:ind+n) + &
                    & t*vylm(ind-n:ind+n)
            end do
        end if
    ! Naive case of rho = 0
    else
        ! Ignore input `m` in case of zero scaling factor
        if (beta .eq. zero) then
            dst_m(1) = eight * src_q * dst_sk(1)
            dst_m(2:) = zero
        ! Update `m` otherwise
        else
            dst_m(1) = beta*dst_m(1) + eight*src_q*dst_sk(1)
            dst_m(2:) = beta * dst_m(2:)
        end if
    end if
end subroutine fmm_p2m_bessel_baseline

subroutine fmm_m2p_baseline(c, src_r, p, vscales, alpha, src_m, beta, dst_v)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Local variables
    real(dp) :: rho, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), rcoef, t, tmp
    integer :: n, ind
    real(dp), external :: dnrm2
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get radius and values of spherical harmonics
    call ylmbas_baseline(c, p, vscales, vylm, vplm, vcos, vsin)
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    rho = dnrm2(3, c, 1)
    if (rho .eq. zero) then
        return
    end if
    ! Compute the actual induced potential
    rcoef = src_r / rho
    t = alpha
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        tmp = dot_product(vylm(ind-n:ind+n), src_m(ind-n:ind+n))
        ! Here vscales(ind)**2 is 4*pi/sqrt(2n+1)
        dst_v = dst_v + t*tmp/vscales(ind)**2
    end do
end subroutine fmm_m2p_baseline

subroutine fmm_m2p_bessel_baseline(c, src_r, p, vscales, alpha, &
        & src_m, beta, dst_v)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Local variables
    real(dp) :: rho, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), rcoef, t, tmp
    real(dp) :: sk(p+1), dk(p+1), src_sk(p+1)
    complex(dp) :: work(max(2, p+1))
    integer :: n, ind
    real(dp), external :: dnrm2
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get radius and values of spherical harmonics
    call ylmbas_baseline(c, p, vscales, vylm, vplm, vcos, vsin)
    ! Get values of the second kind Bessel function
    rho = dnrm2(3, c, 1)
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    if (rho .eq. zero) then
        return
    end if
    call modified_spherical_bessel_second_kind(p, rho, sk, dk, work)
    call modified_spherical_bessel_second_kind(p, src_r, src_sk, dk, work)
    ! Compute the actual induced potential
    do n = 0, p
        ind = n*n + n + 1
        tmp = dot_product(vylm(ind-n:ind+n), src_m(ind-n:ind+n))
        dst_v = dst_v + alpha*sk(n+1)/src_sk(n+1)*tmp
    end do
end subroutine fmm_m2p_bessel_baseline

subroutine fmm_m2p_adj_baseline(c, src_q, dst_r, p, vscales, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), rcoef, t, tmp
    integer :: n, ind
    ! Scale output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    ! In case of zero src_q nothing else is required no matter what is the
    ! value of the induced potential
    if (src_q .eq. zero) then
        return
    end if
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    rho = dnrm2(3, c, 1)
    if (rho .eq. zero) then
        return
    end if
    ! Get radius and values of spherical harmonics
    call ylmbas_baseline(c, p, vscales, vylm, vplm, vcos, vsin)
    ! Compute actual induced potentials
    rcoef = dst_r / rho
    t = src_q
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        dst_m(ind-n:ind+n) = dst_m(ind-n:ind+n) + &
            & t/(vscales(ind)**2)*vylm(ind-n:ind+n)
    end do
end subroutine fmm_m2p_adj_baseline

subroutine fmm_p2l_baseline(c, src_q, dst_r, p, vscales, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)**2)
    ! Local variables
    real(dp) :: rho, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), t, rcoef
    integer :: n, ind
    ! Local expansion represents potential from the outside of the sphere of
    ! the given radius `dst_r`. In case `rho=zero` source particle is inside of
    ! the sphere no matter what radius `r` is, so we just ignore such a case.
    rho = dnrm2(3, c, 1)
    if (rho .eq. zero) then
        if (beta .eq. zero) then
            dst_l = zero
        else
            dst_l = beta * dst_l
        end if
        return
    end if
    ! Get radius and values of spherical harmonics
    call ylmbas_baseline(c, p, vscales, vylm, vplm, vcos, vsin)
    rcoef = dst_r / rho
    t = src_q / rho
    ! Ignore input `m` in case of zero scaling factor
    if (beta .eq. zero) then
        do n = 0, p
            ind = n*n + n + 1
            dst_l(ind-n:ind+n) = t * vylm(ind-n:ind+n)
            t = t * rcoef
        end do
    ! Update `m` otherwise
    else
        do n = 0, p
            ind = n*n + n + 1
            dst_l(ind-n:ind+n) = beta*dst_l(ind-n:ind+n) + t*vylm(ind-n:ind+n)
            t = t * rcoef
        end do
    end if
end subroutine fmm_p2l_baseline

subroutine fmm_l2p_baseline(c, src_r, p, vscales, alpha, src_l, beta, dst_v)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales((p+1)*(p+1)), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Local variables
    real(dp) :: tmp, tmp1, tmp2
    real(dp) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vplm((p+1)**2), vylm((p+1)**2), rcoef
    integer :: n, k, ind
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! In case `rho=zero` values of spherical harmonics make no sense and only
    ! spherical harmonic of degree 0 shall be taken into account
    rho = dnrm2(3, c, 1)
    if (rho .eq. zero) then
        dst_v = dst_v + alpha*src_l(1)/vscales(1)
    ! Compute the actual induced potential otherwise
    else
        ! Get radius and values of spherical harmonics
        call ylmbas_baseline(c, p, vscales, vylm, vplm, vcos, vsin)
        rcoef = rho / src_r
        t = alpha
        do n = 0, p
            ind = n*n + n + 1
            tmp = dot_product(vylm(ind-n:ind+n), src_l(ind-n:ind+n))
            dst_v = dst_v + t*tmp/(vscales(ind)**2)
            t = t * rcoef
        end do
    end if
end subroutine fmm_l2p_baseline

subroutine fmm_l2p_bessel_baseline(c, src_r, p, vscales, alpha, src_l, &
        & beta, dst_v)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales((p+1)*(p+1)), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Local variables
    real(dp) :: rho, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), t, tmp
    real(dp) :: si(p+1), di(p+1), src_si(p+1)
    complex(dp) :: work(max(2, p+1))
    integer :: n, ind
    real(dp), external :: dnrm2
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get radius and values of spherical harmonics
    call ylmbas_baseline(c, p, vscales, vylm, vplm, vcos, vsin)
    ! Get values of the second kind Bessel function
    rho = dnrm2(3, c, 1)
    call modified_spherical_bessel_first_kind(p, rho, si, di, work)
    call modified_spherical_bessel_first_kind(p, src_r, src_si, di, work)
    ! Compute the actual induced potential
    do n = 0, p
        ind = n*n + n + 1
        tmp = dot_product(vylm(ind-n:ind+n), src_l(ind-n:ind+n))
        dst_v = dst_v + alpha*si(n+1)/src_si(n+1)*tmp
    end do
end subroutine fmm_l2p_bessel_baseline

subroutine fmm_l2p_adj_baseline(c, src_q, dst_r, p, vscales, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)**2)
    ! Local variables
    real(dp) :: tmp, tmp1, tmp2
    real(dp) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vplm((p+1)**2), vylm((p+1)**2), rcoef
    integer :: n, k, ind
    ! Scale output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! In case of zero src_q nothing else is required no matter what is the
    ! value of the induced potential
    if (src_q .eq. zero) then
        return
    end if
    ! In case `rho=zero` values of spherical harmonics make no sense and only
    ! spherical harmonic of degree 0 shall be taken into account
    rho = dnrm2(3, c, 1)
    if (rho .eq. zero) then
        dst_l(1) = dst_l(1) + src_q/vscales(1)
    ! Compute the actual induced potential otherwise
    else
        ! Get radius and values of spherical harmonics
        call ylmbas_baseline(c, p, vscales, vylm, vplm, vcos, vsin)
        rcoef = rho / dst_r
        t = src_q
        do n = 0, p
            ind = n*n + n + 1
            !tmp = dot_product(vylm(ind-n:ind+n), src_l(ind-n:ind+n))
            !dst_v = dst_v + t*tmp/(vscales(ind)**2)
            dst_l(ind-n:ind+n) = dst_l(ind-n:ind+n) + &
                & t/(vscales(ind)**2)*vylm(ind-n:ind+n)
            t = t * rcoef
        end do
    end if
end subroutine fmm_l2p_adj_baseline

subroutine polleg_baseline(ctheta, stheta, p, vplm)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: ctheta, stheta
    ! Outputs
    real(dp), intent(out) :: vplm((p+1)**2)
    ! Temporary workspace
    real(dp) :: work(p+1)
    ! Local variables
    integer :: m, ind, l, ind2, vplm_ind
    real(dp) :: fact, pmm, pmm1, pmmo, pll, fm, fl
    ! Init aux factors
    fact = one
    pmm = one
    work(1) = ctheta
    fm = two * ctheta
    do m = 1, p
        work(m+1) = work(m) + fm
    end do
    ! This loop goes over non-negative upper index of P_l^m, namely m, and
    ! defines value for all l from m to p. Polynomials P_l^m are defined only
    ! for l >= |m|. Only positive values of m are filled. Here we
    ! define at first P_m^m, then P_{m+1}^m and then all remaining P_l^m for
    ! l from m+2 to p.
    do m = 0, p
        ! index of P_m^m
        ind = (m+1) * (m+1)
        ! Store P_m^m
        !print *, "l=", m, "z=", pmm
        vplm(ind) = pmm
        if (m .eq. p) then
            return
        end if
        fm = dble(m)
        ! index of P_{m+1}^m
        ind2 = ind + 2*m + 2
        ! P_{m+1}^m
        pmm1 = work(m+1) * pmm
        vplm(ind2) = pmm1
        ! Save value P_m^m for recursion
        pmmo = pmm
        ! Fill values of P_l^m
        do l = m+2, p
            ind2 = ind2 + 2*l
            ! pmm corresponds to P_{l-2}^m
            ! pmm1 corresponds to P_{l-1}^m
            fl = dble(l)
            ! Compute P_l^m
            pll = work(l)*pmm1 - (fl+fm-1)*pmm
            pll = pll / (fl-fm)
            ! Store P_l^m
            !vplm(l*l + l + m + 1) = pll
            vplm(ind2) = pll
            ! Save P_{l-1}^m and P_l^m for recursion
            pmm = pmm1
            pmm1 = pll
        end do
        ! Value of P_{m+1}^{m+1}
        pmm = -pmmo * fact * stheta
        fact = fact + two
    end do
end subroutine polleg_baseline

subroutine ylmbas_baseline(x, p, vscales, vylm, vplm, vcos, vsin)
    integer, intent(in) :: p
    real(dp), dimension(3), intent(in) :: x
    real(dp), dimension((p+1)**2), intent(in) :: vscales
    real(dp), dimension((p+1)**2), intent(out) :: vylm, vplm
    real(dp), dimension(p+1), intent(out) :: vcos, vsin
    integer :: l, m, ind
    real(dp) :: y(3), r, cthe, sthe, cphi, sphi, plm
    real(dp), external :: dnrm2
    r = dnrm2(3, x, 1)
    ! Check zero input, fake it with vector (0,0,1) as zero radius must be
    ! taken into account in some other routine, as spherical harmonics do not
    ! depend on radius
    if (r .eq. zero) then
        y = (/zero, zero, one/)
    ! Do normal calculations otherwise
    else
        y = x / r
    end if
    cthe = y(3)
    sthe = dnrm2(2, y, 1)
    if (sthe .ne. zero) then
        cphi = y(1) / sthe
        sphi = y(2) / sthe
        call trgev(cphi, sphi, p, vcos, vsin)
    else
        cphi = one
        sphi = zero
        vcos = one
        vsin = zero
    endif
    call polleg_baseline(cthe, sthe, p, vplm)
    do l = 0, p
        ind = l**2 + l + 1
        vylm(ind) = vscales(ind) * vplm(ind)
        do m = 1, l
            plm = vplm(ind+m) * vscales(ind+m)
            vylm(ind+m) = plm * vcos(m+1)
            vylm(ind-m) = plm * vsin(m+1)
        enddo
    enddo
end subroutine ylmbas_baseline

! M2M baseline translation (p^4 operations)
! Baseline in terms of operation count: p^4
subroutine fmm_m2m_baseline(c, src_r, dst_r, p, vscales, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1))
    real(dp), intent(in) :: src_m((p+1)*(p+1))
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    real(dp) :: r, r1, r2, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vplm((p+1)*(p+1)), fact(2*p+1), tmpk1, tmpk2, tmpk3, tmp1
    real(dp) :: tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, m, indj, indm, indn, indjn
    real(dp), external :: dnrm2
    stheta = dnrm2(2, c, 1)
    r = dnrm2(3, c, 1)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg_baseline(ctheta, stheta, p, vplm)
        r1 = src_r / dst_r
        r2 = r / dst_r
        pow_r1(1) = r1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Fill square roots of factorials
        fact(1) = 1
        do j = 2, 2*p+1
            fact(j) = sqrt(dble(j-1)) * fact(j-1)
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = vscales(indj) * fact(j-k+1) * fact(j+k+1)
                if (k .ne. 0) then
                    tmp1 = tmp1 * sqrt2
                end if
                do n = 0, j
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp2 = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                        & vscales(indjn) / vscales(indn)
                    do m = max(k+n-j, -n), min(k+j-n, n)
                        indm = indn + abs(m)
                        cphi = vcos(1+abs(m))
                        sphi = vsin(1+abs(m))
                        tmpk1 = tmp2 / fact(n-m+1) / fact(n+m+1) / &
                            & fact(j-n-k+m+1) / fact(j-n+k-m+1) * &
                            & vplm(indm) * vscales(indm)
                        if (mod(abs(k-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        tmpk2 = src_m(indjn+abs(k-m)) * cphi
                        if ((m .ge. 0) .and. (m .le. k)) then
                            sphi = -sphi
                        end if
                        tmpk3 = -src_m(indjn+abs(k-m)) * sphi
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt2
                            tmpk2 = tmpk2 + src_m(indjn-abs(k-m))*sphi
                            tmpk3 = tmpk3 + src_m(indjn-abs(k-m))*cphi
                        end if
                        if (m .gt. k) then
                            tmpk3 = -tmpk3
                        end if
                        dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            dst_m(indj-k) = dst_m(indj-k) + tmpk1*tmpk3
                        end if
                    end do
                end do
            end do
        end do
    else
        r1 = src_r / dst_r
        tmpk1 = r1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_m(k) = dst_m(k) + src_m(k)*tmpk1
            end do
            tmpk1 = tmpk1 * r1
        end do
    end if
end subroutine fmm_m2m_baseline

! M2M baseline translation (p^4 operations)
! Baseline in terms of operation count: p^4
subroutine fmm_m2m_bessel_baseline(c, src_r, dst_r, p, vscales, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    use complex_bessel
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1))
    real(dp), intent(in) :: src_m((p+1)*(p+1))
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    real(dp) :: r, r1, r2, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vplm((p+1)*(p+1)), fact(2*p+1), tmpk1, tmpk2, tmpk3, tmp1
    real(dp) :: tmp2, pow_r1(p+1), pow_r2(p+1)
    real(dp) :: r_si(p+1), src_sk(p+1), dst_sk(p+1), work(p+1)
    complex(dp) :: work_complex(p+1)
    integer :: j, k, n, m, indj, indm, indn, indjn
    real(dp), external :: dnrm2
    stheta = dnrm2(2, c, 1)
    r = dnrm2(3, c, 1)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg_baseline(ctheta, stheta, p, vplm)
        r1 = src_r / dst_r
        r2 = r! / dst_r
        pow_r1(1) = r1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        call modified_spherical_bessel_first_kind(p, r, r_si, work, &
            & work_complex)
        call modified_spherical_bessel_second_kind(p, src_r, src_sk, work, &
            & work_complex)
        call modified_spherical_bessel_second_kind(p, dst_r, dst_sk, work, &
            & work_complex)
        ! Fill square roots of factorials
        fact(1) = 1
        do j = 2, 2*p+1
            fact(j) = sqrt(dble(j-1)) * fact(j-1)
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = vscales(indj) * fact(j-k+1) * fact(j+k+1)
                if (k .ne. 0) then
                    tmp1 = tmp1 * sqrt2
                end if
                do n = 0, j
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    !tmp2 = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                    !    & vscales(indjn) / vscales(indn)
                    tmp2 = tmp1 * & !/ src_sk(j-n+1) * dst_sk(j+1) *
                        & r_si(n+1) / &
                        & vscales(indjn) / vscales(indn) !/ pow_r2(n+1)
                    !tmp2 = tmp1 * r_si(n+1) / &
                    !    & vscales(indjn) / vscales(indn)
                    do m = max(k+n-j, -n), min(k+j-n, n)
                        indm = indn + abs(m)
                        cphi = vcos(1+abs(m))
                        sphi = vsin(1+abs(m))
                        tmpk1 = tmp2 / fact(n-m+1) / fact(n+m+1) / &
                            & fact(j-n-k+m+1) / fact(j-n+k-m+1) * &
                            & vplm(indm) * vscales(indm)
                        if (mod(abs(k-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        tmpk2 = src_m(indjn+abs(k-m)) * cphi
                        if ((m .ge. 0) .and. (m .le. k)) then
                            sphi = -sphi
                        end if
                        tmpk3 = -src_m(indjn+abs(k-m)) * sphi
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt2
                            tmpk2 = tmpk2 + src_m(indjn-abs(k-m))*sphi
                            tmpk3 = tmpk3 + src_m(indjn-abs(k-m))*cphi
                        end if
                        if (m .gt. k) then
                            tmpk3 = -tmpk3
                        end if
                        dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            dst_m(indj-k) = dst_m(indj-k) + tmpk1*tmpk3
                        end if
                    end do
                end do
            end do
        end do
    else
        r1 = src_r / dst_r
        tmpk1 = r1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_m(k) = dst_m(k) + src_m(k)*tmpk1
            end do
            tmpk1 = tmpk1 * r1
        end do
    end if
end subroutine fmm_m2m_bessel_baseline

! Translate local expansion to another sphere
! Baseline in terms of operation count: p^4
subroutine fmm_l2l_baseline(c, src_r, dst_r, p, vscales, src_l, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm (of degree up to p)
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1))
    real(dp), intent(in) :: src_l((p+1)*(p+1))
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    real(dp) :: r, r1, r2, ctheta, stheta, cphi, sphi
    real(dp) :: vcos(p+1), vsin(p+1)
    real(dp) :: vplm((p+1)*(p+1)), fact(2*p+1), tmpk1, tmpk2, tmpk3
    real(dp) :: tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, m, indj, indmk, indn, indjn
    stheta = dnrm2(2, c, 1)
    r = dnrm2(3, c, 1)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg_baseline(ctheta, stheta, p, vplm)
        r1 = r / src_r
        r2 = dst_r / r
        pow_r1(1) = 1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Fill square roots of factorials
        fact(1) = 1
        do j = 2, 2*p+1
            fact(j) = sqrt(dble(j-1)) * fact(j-1)
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = pow_r2(j+1) / fact(j-k+1) / fact(j+k+1) * vscales(indj)
                if (k .ne. 0) then
                    tmp1 = tmp1 * sqrt2
                end if
                do n = j, p
                    indn = n*n + n + 1
                    indjn = (n-j)**2 + (n-j) + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indjn) / vscales(indn)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    do m = k+j-n, k+n-j
                        indmk = indjn + abs(m-k)
                        cphi = vcos(1+abs(m-k))
                        sphi = vsin(1+abs(m-k))
                        tmpk1 = tmp2 * fact(n-m+1) * fact(n+m+1) / &
                            & fact(n-j-m+k+1) / fact(n-j+m-k+1) * &
                            & vplm(indmk) * vscales(indmk)
                        if (mod(abs(k+abs(m-k)-abs(m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        tmpk2 = src_l(indn+abs(m)) * cphi
                        if ((m .ge. 0) .and. (m .le. k)) then
                            sphi = -sphi
                        end if
                        tmpk3 = -src_l(indn+abs(m)) * sphi
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt2
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt2
                            tmpk2 = tmpk2 + src_l(indn-abs(m))*sphi
                            tmpk3 = tmpk3 + src_l(indn-abs(m))*cphi
                        end if
                        if (m .lt. 0) then
                            tmpk3 = -tmpk3
                        end if
                        dst_l(indj+k) = dst_l(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            dst_l(indj-k) = dst_l(indj-k) + tmpk1*tmpk3
                        end if
                    end do
                end do
            end do
        end do
    else
        r1 = dst_r / src_r
        tmpk1 = 1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmpk1
            end do
            tmpk1 = tmpk1 * r1
        end do
    end if
end subroutine fmm_l2l_baseline

! Compute local expansion by given multipole expansion
! Baseline in terms of operation count: p^4
subroutine fmm_m2l_baseline(c, src_r, dst_r, pm, pl, vscales, src_m, dst_l)
! Parameters:
!   c: radius-vector from new (local) to old (multipole) centers of harmonics
!   src_r: radius of old (multipole) harmonics
!   dst_r: radius of new (local) harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm (of degree up to pl+pm)
!   src_m: expansion in old (multipole) harmonics
!   dst_l: expansion in new (local) harmonics
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: c(3), src_r, dst_r
    real(dp), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    real(dp), intent(in) :: src_m((pm+1)*(pm+1))
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    real(dp) :: r, r1, r2, ctheta, stheta, cphi, sphi
    real(dp) :: vcos(pm+pl+1), vsin(pm+pl+1)
    real(dp) :: vplm((pm+pl+1)*(pm+pl+1)), fact(2*(pm+pl)+1), tmpk1, tmpk2
    real(dp) :: tmpk3, tmp1, tmp2, pow_r1(pm+1), pow_r2(pl+1)
    integer :: j, k, n, m, indj, indmk, indn, indjn
    real(dp), external :: dnrm2
    stheta = dnrm2(2, c, 1)
    r = dnrm2(3, c, 1)
    ! r cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (r .eq. 0) then
        return
    end if
    ctheta = c(3) / r
    if (stheta .ne. 0) then
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / r
        call trgev(cphi, sphi, pm+pl, vcos, vsin)
    else
        cphi = 1
        sphi = 0
        vcos = 1
        vsin = 0
    end if
    call polleg_baseline(ctheta, stheta, pm+pl, vplm)
    r1 = src_r / r
    r2 = dst_r / r
    pow_r1(1) = r1
    pow_r2(1) = 1
    do j = 2, pm+1
        pow_r1(j) = pow_r1(j-1) * r1
    end do
    do j = 2, pl+1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    ! Fill square roots of factorials
    fact(1) = 1
    do j = 2, 2*(pm+pl)+1
        fact(j) = sqrt(dble(j-1)) * fact(j-1)
    end do
    do j = 0, pl
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * pow_r2(j+1) / fact(j-k+1) / fact(j+k+1)
            if (k .ne. 0) then
                tmp1 = tmp1 * sqrt2
            end if
            do n = 0, pm
                indn = n*n + n + 1
                indjn = (j+n)**2 + (j+n) + 1
                tmp2 = tmp1 * pow_r1(n+1) / vscales(indjn) / vscales(indn)
                if (mod(n, 2) .eq. 1) then
                    tmp2 = -tmp2
                end if
                do m = -n, n
                    indmk = indjn + abs(m-k)
                    cphi = vcos(1+abs(m-k))
                    sphi = vsin(1+abs(m-k))
                    tmpk1 = tmp2 / fact(n-m+1) / fact(n+m+1) * &
                        & fact(j+n-m+k+1) * fact(j+n+m-k+1) * vplm(indmk) * &
                        & vscales(indmk)
                    if (mod(abs(k+abs(m)-abs(k-m)), 4) .eq. 2) then
                        tmpk1 = -tmpk1
                    end if
                    tmpk2 = src_m(indn+abs(m)) * cphi
                    if ((m .ge. 0) .and. (m .le. k)) then
                        sphi = -sphi
                    end if
                    tmpk3 = -src_m(indn+abs(m)) * sphi
                    if (m .ne. k) then
                        tmpk1 = tmpk1 / sqrt2
                    end if
                    if (m .ne. 0) then
                        tmpk1 = tmpk1 / sqrt2
                        tmpk2 = tmpk2 + src_m(indn-abs(m))*sphi
                        tmpk3 = tmpk3 + src_m(indn-abs(m))*cphi
                    end if
                    if (m .lt. 0) then
                        tmpk3 = -tmpk3
                    end if
                    dst_l(indj+k) = dst_l(indj+k) + tmpk1*tmpk2
                    if (k .ne. 0) then
                        dst_l(indj-k) = dst_l(indj-k) + tmpk1*tmpk3
                    end if
                end do
            end do
        end do
    end do
end subroutine fmm_m2l_baseline

end program test_ddx_core

