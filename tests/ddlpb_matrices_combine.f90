!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/ddlpb.f90
!! Test of analytical derivatives against numerical
!!
!! @version 1.0.0
!! @author Abhinav Jha
!! @date 2021-10-21

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
! derivative_num_force  : Numerical derivatives for Force
! derivative_num_char   : Numerical derivatives for U_i^e(x_in)
real(dp), allocatable :: derivative_num_force(:, :), &
                         & derivative_num_char(:,:)
! start_time    : Start time for the simumations
! finish_time   : Finish time for the simulations
! step          : Step size for the numerical derivatives
! relerr_force  : Relative error for Force
! relerr_char   : Relative error for U_i^e(x_in)
real(dp) :: start_time, finish_time, step, relerr_force, &
            & relerr_char
! isph : Index for number of spheres
! i    : Index for derivative components (i = 1,2,3)
integer :: isph, i, icav, igrid, icav_g1, icav_g2, ibasis
! vsin, vcos, vplm : Values used in basloc
! basloc : Y_lm
! dbasloc : Derivatives of Y_lm
real(dp), allocatable :: vsin(:), vcos(:), vplm(:), basloc(:), dbsloc(:,:), &
                        & phi_grid(:,:), gradphi_grid(:,:,:)
! unit_vector_evaluated_at_grid : Unit vector evaluated on grid points (\sum_lm(Y_lm(s_n)[X]_{ilm})
! unit_vector_ngrid_nsph        : Unit vector of size ngrid \times nsph
! unit_vector_nbasis_nsph       : Unit vector of size nbasis \times nsph
real(dp), allocatable:: unit_vector_evaluated_at_grid(:,:), unit_vector_ngrid_nsph(:,:), &
                        & unit_vector_nbasis_nsph(:,:), Xadj_r(:,:), Xadj_e(:,:), &
                        & Xadj_r_sgrid(:,:), Xadj_e_sgrid(:,:)
! force : Analytic derivative of Force
! derivative_char  : Analytic derivative of U_i^e(x_in)
! diff_re          : Used for derivative in C1_C2 with Q not constant
real(dp), allocatable:: force(:,:), derivative_char(:,:), &
                       & diff_re(:,:)
! sum_esolv_plus_h      : Computation of matrix force evaluated at x+h
! sum_esolv_minus_h     : Computation of matrix force evaluated at x-h
! sum_char_plus_h       : Computation of U_i^e(x_in) evaluated at x+h
! sum_char_minus_h      : Computation of U_i^e(x_in) evaluated at x-h
real(dp) :: sum_esolv_plus_h, sum_esolv_minus_h,&
            & sum_char_minus_h, &
            & sum_char_plus_h
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
allocate(derivative_num_force(3, ddx_data % nsph), &
    & derivative_num_char(3, ddx_data % nsph),&
    & vsin(ddx_data % lmax+1), vcos(ddx_data % lmax+1), &
    & vplm(ddx_data % nbasis), basloc(ddx_data % nbasis), &
    & dbsloc(3, ddx_data % nbasis), &
    & Xadj_r(ddx_data % nbasis, ddx_data % nsph), &
    & Xadj_e(ddx_data % nbasis, ddx_data % nsph), &
    & Xadj_r_sgrid(ddx_data % ngrid, ddx_data % nsph), &
    & Xadj_e_sgrid(ddx_data % ngrid, ddx_data % nsph), &
    & unit_vector_nbasis_nsph(ddx_data % nbasis, ddx_data % nsph), &
    & unit_vector_ngrid_nsph(ddx_data % ngrid, ddx_data % nsph), &
    & unit_vector_evaluated_at_grid(ddx_data % ngrid, ddx_data % nsph), &
    & force(3, ddx_data % nsph), &
    & derivative_char(3, ddx_data % nsph), &
    & diff_re(ddx_data % nbasis, ddx_data % nsph), &
    & phi_cav(ddx_data % ncav), &
    & gradphi_cav(3, ddx_data % ncav), &
    & hessian_cav(3, 3, ddx_data % ncav), &
    & normal_hessian_cav(3, ddx_data % ncav), &
    & phi_grid(ddx_data % ngrid, ddx_data % nsph), &
    & gradphi_grid(3, ddx_data % ngrid, ddx_data % nsph), &
    & phi_n(ddx_data % ngrid, ddx_data % nsph), &
    & ef(3, ddx_data % nsph), &
    & psi(ddx_data % nbasis, ddx_data % nsph))

! Allocation to unity
unit_vector_nbasis_nsph = one
unit_vector_ngrid_nsph = one

! Allocation to zero
diff_re = zero

force = zero
derivative_char = zero
phi_cav = zero
gradphi_cav = zero
psi = zero
phi_grid = zero
gradphi_grid = zero
ef = zero
phi_n = zero
normal_hessian_cav = zero
Xadj_r = zero; Xadj_e = zero
Xadj_r_sgrid = zero; Xadj_e_sgrid = zero

step = 0.00001

! Initialise SI, DI, SK, and DK in the ddLPB unit
call ddlpb_init(ddx_data)

call mkrhs(ddx_data % params, ddx_data % constants, ddx_data % workspace, &
    & phi_cav, gradphi_cav, hessian_cav, psi)
call wghpot(ddx_data, phi_cav, ddx_data % phi_grid, ddx_data % tmp_grid)
call ddx_lpb_adjoint(ddx_data, psi, Xadj_r, Xadj_e)

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

call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % vgrid, ddx_data % vgrid_nbasis, &
            & Xadj_r , ddx_data % nbasis, zero, &
            & Xadj_r_sgrid, &
            & ddx_data % ngrid)

call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % vgrid, ddx_data % vgrid_nbasis, &
            & Xadj_e , ddx_data % nbasis, zero, &
            & Xadj_e_sgrid, &
            & ddx_data % ngrid)


! Compute unit_vector^T*B^k*unit_vector for x_1
do isph = 1, ddx_data % nsph
  ! Computation for matrix A
  call fdoka(ddx_data, isph, unit_vector_nbasis_nsph, &
             & Xadj_r_sgrid(:, isph), &
             & basloc, dbsloc, vplm, vcos, vsin, force(:,isph))
  call fdokb(ddx_data, isph, unit_vector_nbasis_nsph, Xadj_r_sgrid, &
                  & basloc, dbsloc, vplm, vcos, vsin, force(:, isph))
  ! Computation for matrix B
  call fdoka_b_xe(ddx_data, isph, unit_vector_nbasis_nsph, &
                  & Xadj_e_sgrid(:, isph), &
                  & basloc, dbsloc, vplm, vcos, vsin, force(:,isph))
  call fdokb_b_xe(ddx_data, isph, unit_vector_nbasis_nsph, &
                  & Xadj_e_sgrid, &
                  & basloc, dbsloc, vplm, vcos, vsin, force(:, isph))
  ! Computation for derivative of U_i^e(x_in)
  call fdoga(ddx_data, isph, unit_vector_ngrid_nsph, unit_vector_ngrid_nsph, &
             & derivative_char(:, isph))
  ! Computation for matrix C1_C2
  call fdouky(ddx_data, isph, &
                  & unit_vector_nbasis_nsph, &
                  & unit_vector_nbasis_nsph, &
                  & Xadj_r_sgrid, &
                  & Xadj_e_sgrid, &
                  & Xadj_r, &
                  & Xadj_e, &
                  & force(:, isph), &
                  & diff_re)
  call derivative_P(ddx_data, isph,&
                  & unit_vector_nbasis_nsph,&
                  & unit_vector_nbasis_nsph,&
                  & Xadj_r_sgrid,&
                  & Xadj_e_sgrid,&
                  & diff_re, &
                  & force(:, isph))
  ! Computation for G0
  call fdoga(ddx_data, isph, Xadj_r_sgrid, ddx_data % phi_grid, &
            & force(:, isph))
end do
! NOTE: fdoga returns a positive summation
force = -force
icav = 0
do isph = 1, ddx_data % nsph
  do igrid = 1, ddx_data % ngrid
    if(ddx_data % ui(igrid, isph) .ne. zero) then
      icav = icav + 1
      ddx_data % zeta(icav) = -ddx_data % wgrid(igrid) * &
                        & ddx_data % ui(igrid, isph) * ddot(ddx_data % nbasis, &
                        & ddx_data % vgrid(1, igrid), 1, &
                        & Xadj_r(1, isph), 1)
      force(:, isph) = force(:, isph) + &
                        & ddx_data % zeta(icav)*gradphi_cav(:, icav)
    end if
  end do
end do
call efld(ddx_data % ncav, ddx_data % zeta, ddx_data % ccav, &
                & ddx_data % nsph, ddx_data % csph, ef)
do isph = 1, ddx_data % nsph
  force(:, isph) = force(:, isph) &
                          & - ef(:, isph)*ddx_data % charge(isph)
end do

icav_g1 = zero
icav_g2 = zero
do isph = 1, ddx_data % nsph
   ! Computation for F0
  call fdouky_f0(ddx_data, isph ,Xadj_r, &
                & Xadj_r_sgrid, &
                & gradphi_cav, force(:, isph))
  call fdoco(ddx_data, isph, Xadj_r_sgrid, gradphi_cav, &
                & normal_hessian_cav, icav_g1, force(:, isph))
   ! Computation for F0
  call fdouky_f0(ddx_data, isph , Xadj_e, &
                & Xadj_e_sgrid, &
                & gradphi_cav, force(:, isph))
  call fdoco(ddx_data, isph, Xadj_e_sgrid, gradphi_cav, &
                & normal_hessian_cav, icav_g2, force(:, isph))
end do

call ddlpb_free(ddx_data)
do isph = 1, ddx_data % nsph
    do i = 1, 3
        ! Set the initial sums to zero
        sum_esolv_minus_h = zero
        sum_esolv_plus_h = zero
        sum_char_minus_h = zero
        sum_char_plus_h = zero
        ! Set the centers to x + h
        ddx_data % csph(i, isph) = ddx_data % csph(i, isph) + step
        ! Call solve
        call solve(ddx_data, sum_esolv_plus_h, sum_char_plus_h, &
                 & Xadj_r, Xadj_e)
        ! Set the center to x - h
        ddx_data % csph(i, isph) = ddx_data % csph(i, isph) - two*step
        ! Call solve
        call solve(ddx_data, sum_esolv_minus_h, sum_char_minus_h, &
                 & Xadj_r, Xadj_e)
        ! Set the center to x
        ddx_data % csph(i, isph) = ddx_data % csph(i, isph) + step
        ! Numerical derivative = (f(x+h)-f(x-h))/(2*h)
        derivative_num_force(i, isph)= (sum_esolv_plus_h  - sum_esolv_minus_h)/ two / step
        derivative_num_char(i, isph)  = (sum_char_plus_h   - sum_char_minus_h) / two / step
    end do
end do
! Relative Errors
relerr_force = dnrm2(3*ddx_data % nsph, derivative_num_force - force, 1) / &
    & dnrm2(3*ddx_data % nsph, force, 1)
relerr_char = dnrm2(3*ddx_data % nsph, derivative_num_char - derivative_char, 1) / &
    & dnrm2(3*ddx_data % nsph, derivative_char, 1)

write(6,'(2A60)') 'Force Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') force(1,i), force(2,i), &
  & force(3,i), &
  & derivative_num_force(1,i), derivative_num_force(2,i), &
  & derivative_num_force(3,i)
end do

write(6,'(2A60)') ' U_i^e(x_in) Analytical forces', 'Numerical forces'
do i = 1, ddx_data % nsph
  write(6,'(6E20.10)') derivative_char(1,i), derivative_char(2,i), derivative_char(3,i), &
  & derivative_num_char(1,i), derivative_num_char(2,i), derivative_num_char(3,i)
end do

! Deallocation
deallocate(derivative_num_force, derivative_num_char, &
           & vcos, vsin, vplm, basloc, dbsloc, &
           & Xadj_r, Xadj_e, Xadj_r_sgrid, Xadj_e_sgrid, &
           & unit_vector_nbasis_nsph, unit_vector_evaluated_at_grid, &
           & force, diff_re, &
           & ef, phi_n, hessian_cav, normal_hessian_cav)
call ddfree(ddx_data)

write(*, *) "Rel.error of Force :", relerr_force
write(*, *) "Rel.error of U_i^e :", relerr_char
!if (relerr .gt. 1d-6) stop 1
contains

subroutine solve(ddx_data, sum_esolv, sum_char, Xadj_r, Xadj_e)
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), dimension(ddx_data % nbasis, ddx_data % nsph) :: Xadj_r, Xadj_e
    real(dp), intent(out) :: sum_esolv, sum_char
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
                             & vector_f0(:,:)
    real(dp) , dimension(ddx_data % ngrid, ddx_data % nsph):: phi_grid
    real(dp), allocatable :: phi_cav2(:), gradphi_cav2(:, :), psi2(:, :), &
                            & phi_grid2(:,:), tmp_grid2(:,:), hessian_cav2(:,:,:)
    ! i      : Index for n
    ! isph   : Index for number of sphers
    ! igrid  : Index for grid points
    ! ibasis : Index for number of basis
    integer :: i, isph, igrid, ibasis, jsph, kep
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
             & tmp_grid2(ddx_data2 % ngrid, ddx_data2 % nsph))
    ! Intialisation
    unit_vector_n = one
    unit_vector_nbasis_nsph = one
    vector_cosmo = one
    vector_lpb = one
    vector_e_c1_c2 = zero
    vector_r_c1_c2 = zero
    vector_g0 = zero
    vector_f0 = zero
    gradphi_cav2 = zero
    tmp_grid2 = zero
    call ddlpb_init(ddx_data2)
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
        sum_esolv = sum_esolv + &
                   & Xadj_r(ibasis, isph)*(vector_g0(ibasis, isph) + &
                   & vector_f0(ibasis, isph)) + &
                   & Xadj_e(ibasis, isph)*vector_f0(ibasis, isph) - &
                   & Xadj_r(ibasis, isph)*vector_r_c1_c2(ibasis, isph) - &
                   & Xadj_e(ibasis, isph)*vector_e_c1_c2(ibasis, isph)
      end do
    end do

    kep = zero
    ! Sum for matrix A and matrix B
    do isph = 1, ddx_data2 % nsph
      do ibasis = 1, ddx_data2 % nbasis
        kep = kep + 1
        sum_esolv = sum_esolv - Xadj_r(ibasis, isph)*&
                 & vector_cosmo(kep) &
                 & - Xadj_e(ibasis, isph)*&
                 & vector_lpb(kep)
      end do
    end do
    call ddlpb_free(ddx_data2)
    ! Deallocation
    deallocate(unit_vector_n, vector_cosmo, vector_lpb, unit_vector_nbasis_nsph,&
               & vector_r_c1_c2,&
               & vector_e_c1_c2, vector_g0, phi_grid2, psi2, gradphi_cav2, phi_cav2, &
               & vector_f0, tmp_grid2, hessian_cav2)
end subroutine solve

end program main
