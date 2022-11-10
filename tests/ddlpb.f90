!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/ddlpb.f90
!! Test of analytical derivatives against numerical
!!
!! @version 1.0.0
!! @author Abhinav Jha
!! @date 2021-09-01

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
! derivative_num_cosmo : Numerical derivatives for matrix A
! derivative_num_lpb   : Numerical derivatives for matrix B
! derivative_num_char  : Numerical derivatives for U_i^e(x_in)
! derivative_num_u     : Numerical derivatives for C1_C2
! derivative_num_g0    : Numerical derivatives for G0
! derivative_num_f0    : Numerical derivatives for F0
real(dp), allocatable :: derivative_num_cosmo(:, :), derivative_num_lpb(:,:),&
                         & derivative_num_char(:,:), &
                         & derivative_num_u(:,:), &
                         & derivative_num_g0(:,:), &
                         & derivative_num_f0(:,:)
! start_time   : Start time for the simumations
! finish_time  : Finish time for the simulations
! step         : Step size for the numerical derivatives
! relerr_cosmo : Relative error for matrix A
! relerr_lpb   : Relative error for matrix B
! relerr_char  : Relative error for U_i^e(x_in)
! relerr_der_u : Relative error for C1_C2
! relerr_g0    : Relative error for G0
! relerr_f0    : Relative error for F0
real(dp) :: start_time, finish_time, step, relerr_cosmo, relerr_lpb, relerr_char, relerr_der_u, &
            & relerr_g0, relerr_f0
! isph : Index for number of spheres
! i    : Index for derivative components (i = 1,2,3)
integer :: isph, i, icav, igrid, icav_g, ibasis
! vsin, vcos, vplm : Values used in basloc
! basloc : Y_lm
! dbasloc : Derivatives of Y_lm
real(dp), allocatable :: vsin(:), vcos(:), vplm(:), basloc(:), dbsloc(:,:), &
                        & phi_grid(:,:), gradphi_grid(:,:,:)
! unit_vector_evaluated_at_grid : Unit vector evaluated on grid points (\sum_lm(Y_lm(s_n)[X]_{ilm})
! unit_vector_ngrid_nsph        : Unit vector of size ngrid \times nsph
! unit_vector_nbasis_nsph       : Unit vector of size nbasis \times nsph
real(dp), allocatable:: unit_vector_evaluated_at_grid(:,:), unit_vector_ngrid_nsph(:,:), &
                        & unit_vector_nbasis_nsph(:,:), lpb_vector(:,:), &
                        & lpb_adjoint_vector(:,:)
! derivative_cosmo : Analytic derivative of matrix A
! derivative_lpb   : Analytic derivative of matrix B
! derivative_char  : Analytic derivative of U_i^e(x_in)
! derivative_u_r   : Analytic derivative of C1_C2 with Q constant and corresponding to ddCOSMO
! derivative_u_e   : Analytic derivative of C1_C2 with Q constant and corresponding to ddLPB
! derivative_g0    : Analytic derivative of G0 RHS
! derivative_f0    : Analytic derivative of F0 RHS
! diff_re          : Used for derivative in C1_C2 with Q not constant
real(dp), allocatable:: derivative_cosmo(:,:), derivative_lpb(:,:), derivative_char(:,:), &
                       & derivative_u_r(:,:), derivative_u_e(:,:), derivative_g0(:,:), &
                       & diff_re(:,:), derivative_f0(:,:)
! sum_cosmo_plus_h  : Computation of matrix A evaluated at x+h
! sum_cosmo_minus_h : Computation of matrix A evaluated at x-h
! sum_lpb_plus_h    : Computation of matrix B evaluated at x+h
! sum_lpb_minus_h   : Computation of matrix B evaluated at x-h
! sum_char_plus_h   : Computation of U_i^e(x_in) evaluated at x+h
! sum_char_minus_h  : Computation of U_i^e(x_in) evaluated at x-h
! sum_der_u_plus_h  : Computation of matrix C1_C2 evaluated at x+h
! sum_der_u_plus_h  : Computation of matrix C1_C2 evaluated at x-h
! sum_der_g0_plus_h : Computation of G0 evaluated at x+h
! sum_der_g0_minus_h: Computation of G0 evaluated at x-h
! sum_der_f0_plus_h : Computation of F0 evaluated at x+h
! sum_der_f0_minus_h: Computation of F0 evaluated at x-h
real(dp) :: sum_cosmo_plus_h, sum_lpb_plus_h, sum_cosmo_minus_h, sum_lpb_minus_h,&
            & sum_char_minus_h, &
            & sum_char_plus_h, sum_der_u_minus_h, sum_der_u_plus_h, bessel_argument, &
            & sum_der_g0_plus_h, sum_der_g0_minus_h, sum_der_f0_plus_h, &
            & sum_der_f0_minus_h, sum_lpb, sum_lpb_adjoint
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), psi(:, :), ef(:,:), &
                        & phi_n(:,:), hessian_cav(:,:,:), normal_hessian_cav(:,:)
real(dp), external :: dnrm2, ddot

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file 12"
call ddfromfile(fname, ddx_data, info)
if(info .ne. 0) stop "info != 0"

! lmax0 set to minimum of 6 or given lmax.
! nbasis0 set to minimum of 49 or given (lmax+1)^2.
! Maybe do it ddlpb_init(?)
lmax0 = MIN(6, ddx_data % lmax)
nbasis0 = MIN(49, ddx_data % nbasis)

! Allocation for various vectors
allocate(derivative_num_cosmo(3, ddx_data % nsph), &
    & derivative_num_lpb(3, ddx_data % nsph),&
    & derivative_num_char(3, ddx_data % nsph),&
    & derivative_num_u(3, ddx_data % nsph), &
    & derivative_num_g0(3, ddx_data % nsph), &
    & derivative_num_f0(3, ddx_data % nsph), &
    & vsin(ddx_data % lmax+1), vcos(ddx_data % lmax+1), &
    & vplm(ddx_data % nbasis), basloc(ddx_data % nbasis), &
    & dbsloc(3, ddx_data % nbasis), &
    & unit_vector_nbasis_nsph(ddx_data % nbasis, ddx_data % nsph), &
    & unit_vector_ngrid_nsph(ddx_data % ngrid, ddx_data % nsph), &
    & unit_vector_evaluated_at_grid(ddx_data % ngrid, ddx_data % nsph), &
    & derivative_cosmo(3, ddx_data % nsph), &
    & derivative_lpb(3, ddx_data % nsph), &
    & derivative_char(3, ddx_data % nsph), &
    & derivative_u_r(3, ddx_data % nsph), &
    & derivative_u_e(3, ddx_data % nsph), &
    & derivative_g0(3, ddx_data % nsph), &
    & derivative_f0(3, ddx_data % nsph), &
    & diff_re(ddx_data % nbasis, ddx_data % nsph), &
    & phi_cav(ddx_data % ncav), &
    & gradphi_cav(3, ddx_data % ncav), &
    & hessian_cav(3, 3, ddx_data % ncav), &
    & normal_hessian_cav(3, ddx_data % ncav), &
    & phi_grid(ddx_data % ngrid, ddx_data % nsph), &
    & gradphi_grid(3, ddx_data % ngrid, ddx_data % nsph), &
    & lpb_vector(ddx_data % nbasis, ddx_data % nsph), &
    & lpb_adjoint_vector(ddx_data % nbasis, ddx_data % nsph), &
    & phi_n(ddx_data % ngrid, ddx_data % nsph), &
    & ef(3, ddx_data % nsph), &
    & psi(ddx_data % nbasis, ddx_data % nsph))

! Allocation to unity
unit_vector_nbasis_nsph = one
unit_vector_ngrid_nsph = one
lpb_vector = one
lpb_adjoint_vector = one
sum_lpb = one
sum_lpb_adjoint = one

! Allocation to zero
diff_re = zero

derivative_cosmo = zero
derivative_lpb = zero
derivative_char = zero
derivative_u_e = zero
derivative_u_r = zero
derivative_g0 = zero
derivative_f0 = zero
phi_cav = zero
gradphi_cav = zero
psi = zero
phi_grid = zero
gradphi_grid = zero
ef = zero
phi_n = zero
normal_hessian_cav = zero

step = 0.00001

! Initialise SI, DI, SK, and DK in the ddLPB unit
call ddlpb_init(ddx_data)

call mkrhs(ddx_data % params, ddx_data % constants, ddx_data % workspace, &
    & phi_cav, gradphi_cav, hessian_cav, psi)
call wghpot(ddx_data, phi_cav, ddx_data % phi_grid, ddx_data % tmp_grid)
icav = 0
do isph = 1, ddx_data % nsph
  do igrid = 1, ddx_data % ngrid
    if(ddx_data % ui(igrid, isph) .gt. zero) then
      icav = icav + 1
      do i = 1, 3
        normal_hessian_cav(:, icav) = normal_hessian_cav(:,icav) +&
                                    & hessian_cav(:,i,icav)*ddx_data % cgrid(i,igrid)
      end do
    end if
  end do
end do

call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % vgrid, ddx_data % vgrid_nbasis, &
            & unit_vector_nbasis_nsph , ddx_data % nbasis, zero, &
            & unit_vector_evaluated_at_grid, &
            & ddx_data % ngrid)


! Compute unit_vector^T*B^k*unit_vector for x_1
icav_g = zero
do isph = 1, ddx_data % nsph
  ! Computation for matrix A
  call fdoka(ddx_data, isph, unit_vector_nbasis_nsph, &
             & unit_vector_evaluated_at_grid(:, isph), &
             & basloc, dbsloc, vplm, vcos, vsin, derivative_cosmo(:,isph))
  call fdokb(ddx_data, isph, unit_vector_nbasis_nsph, unit_vector_evaluated_at_grid, &
                  & basloc, dbsloc, vplm, vcos, vsin, derivative_cosmo(:, isph))
  ! Computation for matrix B
  call fdoka_b_xe(ddx_data, isph, unit_vector_nbasis_nsph, unit_vector_evaluated_at_grid(:, isph), &
                  & basloc, dbsloc, vplm, vcos, vsin, derivative_lpb(:,isph))
  call fdokb_b_xe(ddx_data, isph, unit_vector_nbasis_nsph, unit_vector_evaluated_at_grid, &
                  & basloc, dbsloc, vplm, vcos, vsin, derivative_lpb(:, isph))
  ! Computation for derivative of U_i^e(x_in)
  call fdoga(ddx_data, isph, unit_vector_ngrid_nsph, unit_vector_ngrid_nsph, &
             & derivative_char(:, isph))
  ! Computation for matrix C1_C2
  call fdouky(ddx_data, isph, &
                  & unit_vector_nbasis_nsph, &
                  & unit_vector_nbasis_nsph, &
                  & unit_vector_evaluated_at_grid, &
                  & unit_vector_evaluated_at_grid, &
                  & unit_vector_nbasis_nsph, &
                  & unit_vector_nbasis_nsph, &
                  & derivative_u_r(:, isph), &
                  & derivative_u_e(:, isph), &
                  & diff_re)
  call derivative_P(ddx_data, isph,&
                  & unit_vector_nbasis_nsph,&
                  & unit_vector_nbasis_nsph,&
                  & unit_vector_evaluated_at_grid,&
                  & unit_vector_evaluated_at_grid,&
                  & diff_re, &
                  & derivative_u_r(:, isph), &
                  & derivative_u_e(:, isph))
  ! Computation for G0
  call fdoga(ddx_data, isph, unit_vector_evaluated_at_grid, ddx_data % phi_grid, &
            & derivative_g0(:, isph))
  ! Computation for F0
  call fdouky_f0(ddx_data, isph , unit_vector_evaluated_at_grid, &
                 & unit_vector_nbasis_nsph, gradphi_cav, derivative_f0(:, isph))
  call fdoco(ddx_data, isph, unit_vector_evaluated_at_grid, gradphi_cav, &
             & normal_hessian_cav, icav_g, derivative_f0(:, isph))
end do
! NOTE: fdoga returns a positive summation
derivative_g0 = -derivative_g0
icav = 0
do isph = 1, ddx_data % nsph
  do igrid = 1, ddx_data % ngrid
    if(ddx_data % ui(igrid, isph) .ne. zero) then
      icav = icav + 1
      ddx_data % zeta(icav) = -ddx_data % wgrid(igrid) * &
                        & ddx_data % ui(igrid, isph) * ddot(ddx_data % nbasis, &
                        & ddx_data % vgrid(1, igrid), 1, &
                        & unit_vector_nbasis_nsph(1, isph), 1)
      derivative_g0(:, isph) = derivative_g0(:, isph) + &
                        & ddx_data % zeta(icav)*gradphi_cav(:, icav)
    end if
  end do
end do
call efld(ddx_data % ncav, ddx_data % zeta, ddx_data % ccav, &
                & ddx_data % nsph, ddx_data % csph, ef)
do isph = 1, ddx_data % nsph
  derivative_g0(:, isph) = derivative_g0(:, isph) - ef(:, isph)*ddx_data % charge(isph)
end do

!call debug_bessel(ddx_data, bessel_argument)


do isph = 1, ddx_data % nsph
    do i = 1, 3
        ! Set the initial sums to zero
        sum_cosmo_minus_h = zero
        sum_cosmo_plus_h = zero
        sum_lpb_minus_h = zero
        sum_lpb_plus_h = zero
        sum_char_minus_h = zero
        sum_char_plus_h = zero
        sum_der_u_minus_h = zero
        sum_der_u_plus_h = zero
        sum_der_g0_minus_h = zero
        sum_der_g0_plus_h = zero
        sum_der_f0_minus_h = zero
        sum_der_f0_plus_h = zero
        ! Set the centers to x + h
        ddx_data % csph(i, isph) = ddx_data % csph(i, isph) + step
        ! Call solve
        call solve(ddx_data, sum_cosmo_plus_h, sum_lpb_plus_h, sum_char_plus_h, &
                   & sum_der_u_plus_h, sum_der_g0_plus_h, sum_der_f0_plus_h)
        ! Set the center to x - h
        ddx_data % csph(i, isph) = ddx_data % csph(i, isph) - two*step
        ! Call solve
        call solve(ddx_data, sum_cosmo_minus_h, sum_lpb_minus_h, sum_char_minus_h, &
                   & sum_der_u_minus_h, sum_der_g0_minus_h, sum_der_f0_minus_h)
        ! Set the center to x
        ddx_data % csph(i, isph) = ddx_data % csph(i, isph) + step
        ! Numerical derivative = (f(x+h)-f(x-h))/(2*h)
        derivative_num_cosmo(i, isph) = (sum_cosmo_plus_h  - sum_cosmo_minus_h) / two / step
        derivative_num_lpb(i, isph)   = (sum_lpb_plus_h    - sum_lpb_minus_h) / two / step
        derivative_num_char(i, isph)  = (sum_char_plus_h   - sum_char_minus_h) / two / step
        derivative_num_u(i, isph)     = (sum_der_u_plus_h  - sum_der_u_minus_h) / two / step
        derivative_num_g0(i, isph)    = (sum_der_g0_plus_h - sum_der_g0_minus_h) / two / step
        derivative_num_f0(i, isph)    = (sum_der_f0_plus_h - sum_der_f0_minus_h) / two / step
    end do
end do
! Relative Errors
relerr_cosmo = dnrm2(3*ddx_data % nsph, derivative_num_cosmo - derivative_cosmo, 1) / &
    & dnrm2(3*ddx_data % nsph, derivative_cosmo, 1)
relerr_lpb = dnrm2(3*ddx_data % nsph, derivative_num_lpb - derivative_lpb, 1) / &
    & dnrm2(3*ddx_data % nsph, derivative_lpb, 1)
relerr_char = dnrm2(3*ddx_data % nsph, derivative_num_char - derivative_char, 1) / &
    & dnrm2(3*ddx_data % nsph, derivative_char, 1)
relerr_der_u = dnrm2(3*ddx_data % nsph, derivative_num_u - (derivative_u_e + derivative_u_r), 1) / &
    & dnrm2(3*ddx_data % nsph, (derivative_u_e + derivative_u_r), 1)
relerr_g0 = dnrm2(3*ddx_data % nsph, derivative_num_g0 - derivative_g0, 1) / &
    & dnrm2(3*ddx_data % nsph, derivative_g0, 1)
relerr_f0 = dnrm2(3*ddx_data % nsph, derivative_num_f0 - derivative_f0, 1) / &
    & dnrm2(3*ddx_data % nsph, derivative_f0, 1)

call matABx(ddx_data, ddx_data % n, unit_vector_nbasis_nsph, lpb_vector)
call bstarx(ddx_data, unit_vector_nbasis_nsph, lpb_adjoint_vector)

do ibasis = 1, ddx_data % nbasis
  do isph = 1, ddx_data % nsph
    sum_lpb = sum_lpb + lpb_vector(ibasis, isph)
    sum_lpb_adjoint = sum_lpb_adjoint + lpb_adjoint_vector(ibasis, isph)
  end do
end do


write(6,'(2A60)') 'Matrix A Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') derivative_cosmo(1,i), derivative_cosmo(2,i), derivative_cosmo(3,i), &
  & derivative_num_cosmo(1,i), derivative_num_cosmo(2,i), derivative_num_cosmo(3,i)
end do

write(6,'(2A60)') ' Matrix B Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') derivative_lpb(1,i), derivative_lpb(2,i), derivative_lpb(3,i), &
  & derivative_num_lpb(1,i), derivative_num_lpb(2,i), derivative_num_lpb(3,i)
end do

write(6,'(2A60)') ' U_i^e(x_in) Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') derivative_char(1,i), derivative_char(2,i), derivative_char(3,i), &
  & derivative_num_char(1,i), derivative_num_char(2,i), derivative_num_char(3,i)
end do

write(6,'(2A60)') ' Matrix C1_C2 Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') derivative_u_r(1,i) + derivative_u_e(1,i), &
                     & derivative_u_r(2,i) + derivative_u_e(2,i), &
                     & derivative_u_r(3,i) + derivative_u_e(3,i), &
                     & derivative_num_u(1,i), derivative_num_u(2,i), &
                     & derivative_num_u(3,i)
end do

write(6,'(2A60)') ' G0 Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') derivative_g0(1,i), derivative_g0(2,i), derivative_g0(3,i), &
  & derivative_num_g0(1,i), derivative_num_g0(2,i), derivative_num_g0(3,i)
end do

write(6,'(2A60)') ' F0 Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') derivative_f0(1,i), derivative_f0(2,i), derivative_f0(3,i), &
  & derivative_num_f0(1,i), derivative_num_f0(2,i), derivative_num_f0(3,i)
end do

! Deallocation
deallocate(derivative_num_cosmo, derivative_num_lpb, derivative_num_char, &
           & derivative_num_g0, &
           & derivative_num_f0, &
           & vcos, vsin, vplm, basloc, dbsloc, &
           & unit_vector_nbasis_nsph, unit_vector_evaluated_at_grid, &
           & derivative_u_e, derivative_u_r, &
           & derivative_cosmo, derivative_lpb, diff_re, derivative_num_u, &
           & derivative_g0, derivative_f0, ef, phi_n, hessian_cav, normal_hessian_cav)
call ddfree(ddx_data)

write(*, *) "Rel.error of A     :", relerr_cosmo
write(*, *) "Rel.error of B     :", relerr_lpb
write(*, *) "Rel.error of U_i^e :", relerr_char
write(*, *) "Rel.error of C1_C2 :", relerr_der_u
write(*, *) "Rel.error of G0    :", relerr_g0
write(*, *) "Rel.error of F0    :", relerr_f0
!if (relerr .gt. 1d-6) stop 1
contains

subroutine solve(ddx_data, sum_cosmo, sum_lpb, sum_char, sum_der_u, sum_g0, sum_f0)
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(out) :: sum_cosmo, sum_lpb, sum_char, sum_der_u, sum_g0, sum_f0
    ! Local variables
    ! ddx_data2  : New ddx_data with new coordinates
    type(ddx_type) :: ddx_data2
    ! unit_vector_n           : Unit vector of size n\times 1
    ! vector_cosmo            : y = Ax
    ! vector_lpb              : y = Bx
    ! unit_vector_nbasis_nsph : Unit vector of size nbasis \times nsph
    ! vector_e_c1_c2          : y = C1Xr + C2Xe
    ! vector_r_c1_c2          : y = C1Xr + C2Xe
    real(dp), allocatable :: unit_vector_n(:), vector_cosmo(:), vector_lpb(:),&
                             & unit_vector_nbasis_nsph(:,:),&
                             & vector_e_c1_c2(:,:), vector_r_c1_c2(:,:), vector_g0(:,:), &
                             & vector_f0(:,:), &
                             & unit_vector_ncav(:)
    real(dp) , dimension(ddx_data % ngrid, ddx_data % nsph):: phi_grid
    real(dp), allocatable :: phi_cav2(:), gradphi_cav2(:, :), psi2(:, :), &
                            & phi_grid2(:,:), tmp_grid2(:,:), hessian_cav2(:,:,:)
    ! i      : Index for n
    ! isph   : Index for number of sphers
    ! igrid  : Index for grid points
    ! ibasis : Index for number of basis
    integer :: i, isph, igrid, ibasis, jsph
    real(dp) :: v, vij(3), rijn

    ! Initialise new ddx_data with new centers coordinates
    call ddinit(ddx_data % nsph, ddx_data % charge, ddx_data % csph(1, :), &
        & ddx_data % csph(2, :), ddx_data % csph(3, :), ddx_data % rsph, &
        & ddx_data % model, ddx_data % lmax, ddx_data % ngrid, 0, &
        & ddx_data % fmm, ddx_data % pm, ddx_data % pl, &
        & ddx_data % fmm_precompute, ddx_data % iprint, ddx_data % se, &
        & ddx_data % eta, ddx_data % eps, ddx_data % kappa, &
        & ddx_data % matvecmem, &
        & ddx_data % tol, ddx_data % maxiter, &
        & ddx_data % ndiis, ddx_data % nproc, ddx_data2, info)
    ! Allocation
    allocate(unit_vector_n(ddx_data2 % n), vector_cosmo(ddx_data2 % n), vector_lpb(ddx_data2 % n), &
             & unit_vector_nbasis_nsph(ddx_data2 % nbasis, ddx_data2 % nsph), &
             & vector_r_c1_c2(ddx_data2 % nbasis, ddx_data2 % nsph), &
             & vector_e_c1_c2(ddx_data2 % nbasis, ddx_data2 % nsph), &
             & vector_g0(ddx_data2 % nbasis, ddx_data2 % nsph), &
             & vector_f0(ddx_data2 % nbasis, ddx_data2 % nsph), &
             & phi_cav2(ddx_data2 % ncav), &
             & phi_grid2(ddx_data2 % ngrid, ddx_data2 % nsph), &
             & gradphi_cav2(3, ddx_data2 % ncav), &
             & hessian_cav2(3, 3, ddx_data2 % ncav), &
             & psi2(ddx_data2 % nbasis, ddx_data2 % nsph), &
             & unit_vector_ncav(ddx_data2 % ncav), &
             & tmp_grid2(ddx_data2 % ngrid, ddx_data2 % nsph))
    ! Intialisation
    unit_vector_n = one
    unit_vector_nbasis_nsph = one
    unit_vector_ncav = one
    vector_cosmo = one
    vector_lpb = one
    vector_e_c1_c2 = zero
    vector_r_c1_c2 = zero
    vector_g0 = zero
    vector_f0 = zero
    gradphi_cav2 = zero
    tmp_grid2 = zero
    ! Call for matrix A
    call lx(ddx_data2, unit_vector_n, vector_cosmo)
    ! Call for matrix B
    call matABx(ddx_data2 , ddx_data2 % n, unit_vector_n, vector_lpb)
    ! Call for matrix C1_C2
    call C1_C2(ddx_data2, vector_r_c1_c2, vector_e_c1_c2, unit_vector_nbasis_nsph,&
               & unit_vector_nbasis_nsph)

    ! Call for G0
    call mkrhs(ddx_data2 % params, ddx_data2 % constants, ddx_data2 % workspace, &
        & phi_cav2, gradphi_cav2, hessian_cav2, psi2)
    call wghpot(ddx_data2, phi_cav2, ddx_data2 % phi_grid, ddx_data2 % tmp_grid)
    ! Call for F0
    call wghpot_debug(ddx_data2, gradphi_cav2, tmp_grid2)

    do isph = 1, ddx_data2 % nsph
      !! intrhs is a subroutine in ddx_operators
      !! @param[in]  isph : Sphere number, used for output
      !! @param[in]  g    : Intermediate right side g
      !! @param[out] g0   : Integrated right side Eq.(77) in QSM19.SISC
      call intrhs(ddx_data2 % iprint, ddx_data2 % ngrid, &
                ddx_data2 % lmax, ddx_data2 % vwgrid, ddx_data2 % vgrid_nbasis, &
                & isph, ddx_data2 % tmp_grid(:,isph), vector_g0(:, isph))
      call intrhs(ddx_data2 % iprint, ddx_data2 % ngrid, &
                ddx_data2 % lmax, ddx_data2 % vwgrid, ddx_data2 % vgrid_nbasis, &
                & isph, tmp_grid2(:,isph), vector_f0(:, isph))
    end do
    ! Sum for U_i^e(x_in)
    do isph = 1,ddx_data2 % nsph
      do igrid = 1,ddx_data2 % ngrid
        if (ddx_data2 % ui(igrid, isph) .gt. zero) then
          sum_char = sum_char + ddx_data2 % wgrid(igrid) * &
                        & ddx_data2 % ui(igrid, isph)
        end if
      end do
    end do
    ! Sum for matrix C1_C2 with Q constant
    do ibasis = 1, ddx_data2 % nbasis
      do isph = 1, ddx_data2 % nsph
        sum_der_u = sum_der_u + vector_e_c1_c2(ibasis, isph) + vector_r_c1_c2(ibasis, isph)
        sum_g0 = sum_g0 + vector_g0(ibasis, isph)
        sum_f0 = sum_f0 + vector_f0(ibasis, isph)
      end do
    end do
    ! Sum for matrix A and matrix B
    do i = 1, ddx_data2 % n
      sum_cosmo = sum_cosmo + vector_cosmo(i)
      sum_lpb = sum_lpb + vector_lpb(i)
    end do
    ! Deallocation
    deallocate(unit_vector_n, vector_cosmo, vector_lpb, unit_vector_nbasis_nsph, vector_r_c1_c2,&
               & vector_e_c1_c2, vector_g0, phi_grid2, psi2, gradphi_cav2, phi_cav2, &
               & unit_vector_ncav, vector_f0, tmp_grid2, hessian_cav2)
end subroutine solve

end program main
