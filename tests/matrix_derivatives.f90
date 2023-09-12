!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/ddlpb.f90
!! Test of analytical derivatives against numerical for each matrices
!!
!! @version 1.0.0
!! @author Abhinav Jha
!! @date 2022-02-22

program main
use ddx_core
use ddx_operators
use ddx_solvers
use ddx
use ddx_lpb
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
type(ddx_error_type) :: ddx_error
character(len=255) :: dummy_file_name = ''
! derivative_num_A  : Numerical derivatives for matrix A
! derivative_num_B  : Numerical derivatives for matrix B
! derivative_num_Ui : Numerical derivatives for U_i^e(x_in) (Char. function)
! derivative_num_C  : Numerical derivatives for C1_C2
real(dp), allocatable :: derivative_num_A(:, :), derivative_num_B(:,:),&
                         & derivative_num_Ui(:,:), &
                         & derivative_num_C(:,:)
! step         : Step size for the numerical derivatives
! relerr_A     : Relative error for matrix A
! relerr_B     : Relative error for matrix B
! relerr_Ui    : Relative error for U_i^e(x_in)
! relerr_C1_C2 : Relative error for C1_C2
real(dp) :: step, relerr_A = zero, relerr_B = zero, relerr_Ui = zero,&
            & relerr_C1_C2 = zero, tol = zero
! isph   : Index for number of spheres
! i      : Index for derivative components (i = 1,2,3)
! ibasis : Index for number of basis
integer :: isph, i
! vsin, vcos, vplm : Values used in basloc
! basloc : Y_lm
! dbasloc : Derivatives of Y_lm
real(dp), allocatable :: vsin(:), vcos(:), vplm(:), basloc(:), dbsloc(:,:)
! random_vector_two_evaluated_at_grid : Vector evaluated on grid points
!                                 (\sum_lm(Y_lm(s_n)[X]_{ilm})
! vector_ngrid_nsph        : Vector of size ngrid \times nsph
! vector_nbasis_nsph       : Vector of size nbasis \times nsph
real(dp), allocatable:: random_vector_two_evaluated_at_grid(:,:),&
                        & vector_ngrid_nsph(:,:), &
                        & random_vector_nbasis_nsph_one(:,:), &
                        & random_vector_nbasis_nsph_two(:,:)
! derivative_A     : Analytic derivative of matrix A
! derivative_B     : Analytic derivative of matrix B
! derivative_Ui    : Analytic derivative of U_i^e(x_in)
! derivative_C1_C2 : Analytic derivative of C1_C2
! diff_re          : Used for derivative in C1_C2 with Q not constant
real(dp), allocatable:: derivative_A(:,:), derivative_B(:,:), &
                       & derivative_Ui(:,:), &
                       & derivative_C1_C2(:,:), &
                       & diff_re(:,:)
real(dp), allocatable :: charges(:)
! sum_A_plus_h  : Computation of matrix A evaluated at x+h
! sum_A_minus_h : Computation of matrix A evaluated at x-h
! sum_B_plus_h    : Computation of matrix B evaluated at x+h
! sum_B_minus_h   : Computation of matrix B evaluated at x-h
! sum_Ui_plus_h   : Computation of U_i^e(x_in) evaluated at x+h
! sum_Ui_minus_h  : Computation of U_i^e(x_in) evaluated at x-h
! sum_C_plus_h  : Computation of matrix C1_C2 evaluated at x+h
! sum_C_plus_h  : Computation of matrix C1_C2 evaluated at x-h
real(dp) :: sum_A_plus_h, sum_B_plus_h, sum_A_minus_h, sum_B_minus_h,&
            & sum_Ui_minus_h, &
            & sum_Ui_plus_h, sum_C_minus_h, sum_C_plus_h, &
            & lmax0, nbasis0
real(dp), external :: dnrm2, ddot

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file 12"
call ddfromfile(fname, ddx_data, tol, charges, ddx_error)
call check_error(ddx_error)

! lmax0 set to minimum of 6 or given lmax.
! nbasis0 set to minimum of 49 or given (lmax+1)^2.
! Maybe do it ddlpb_init(?)
lmax0 = MIN(6, ddx_data % params % lmax)
nbasis0 = MIN(49, ddx_data % constants % nbasis)

! Allocation for various vectors
allocate(derivative_num_A(3, ddx_data % params % nsph), &
    & derivative_num_B(3, ddx_data % params % nsph),&
    & derivative_num_Ui(3, ddx_data % params % nsph),&
    & derivative_num_C(3, ddx_data % params % nsph), &
    & vsin(ddx_data % params % lmax+1), vcos(ddx_data % params % lmax+1), &
    & vplm(ddx_data % constants % nbasis), &
    & basloc(ddx_data % constants % nbasis), &
    & dbsloc(3, ddx_data % constants % nbasis), &
    & random_vector_nbasis_nsph_one(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & random_vector_nbasis_nsph_two(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & vector_ngrid_nsph(ddx_data % params % ngrid, ddx_data % params % nsph), &
    & random_vector_two_evaluated_at_grid(ddx_data % params % ngrid, &
    & ddx_data % params % nsph), &
    & derivative_A(3, ddx_data % params % nsph), &
    & derivative_B(3, ddx_data % params % nsph), &
    & derivative_Ui(3, ddx_data % params % nsph), &
    & derivative_C1_C2(3, ddx_data % params % nsph), &
    & diff_re(ddx_data % constants % nbasis, ddx_data % params % nsph))

! Random vectors
random_vector_nbasis_nsph_one = one
random_vector_nbasis_nsph_two = one
random_vector_nbasis_nsph_one(1,1) = 2.0
random_vector_nbasis_nsph_two(2,1) = 3.0

! For Ui the unit vector works well because the contraction is from
! one side
vector_ngrid_nsph = one

! Allocation to zero
diff_re = zero

derivative_A = zero
derivative_B = zero
derivative_Ui = zero
derivative_C1_C2 = zero

step = 1d-7

call dgemm('T', 'N', ddx_data % params % ngrid, ddx_data % params % nsph, &
            & ddx_data % constants % nbasis, one, ddx_data % constants % vgrid, ddx_data % constants % vgrid_nbasis, &
            & random_vector_nbasis_nsph_two, &
            & ddx_data % constants % nbasis, zero, &
            & random_vector_two_evaluated_at_grid, &
            & ddx_data % params % ngrid)


! Compute unit_vector^T*B^k*unit_vector for x_1
do isph = 1, ddx_data % params % nsph
  ! Computation of derivative for matrix A
  call contract_grad_L(ddx_data % params, ddx_data % constants, isph, &
             & random_vector_nbasis_nsph_one, &
             & random_vector_two_evaluated_at_grid, &
             & basloc, dbsloc, vplm, vcos, vsin, derivative_A(:,isph))
  ! Computation of derivative for matrix B
  call contract_grad_B(ddx_data % params, ddx_data % constants, &
                 & isph, random_vector_nbasis_nsph_one, random_vector_two_evaluated_at_grid, &
                 & derivative_B(:,isph))
  ! Computation for derivative of U_i^e(x_in)
  call contract_grad_U(ddx_data % params, ddx_data % constants, &
            & isph, vector_ngrid_nsph, vector_ngrid_nsph, &
            & derivative_Ui(:, isph))
end do

! Computation for matrix C1_C2
call contract_grad_C(ddx_data % params, ddx_data % constants, &
             & ddx_data % workspace, &
             & random_vector_nbasis_nsph_one, &
             & random_vector_nbasis_nsph_one, &
             & random_vector_two_evaluated_at_grid, &
             & random_vector_two_evaluated_at_grid, &
             & random_vector_nbasis_nsph_two, &
             & random_vector_nbasis_nsph_two, &
             & derivative_C1_C2, &
             & diff_re, ddx_error)

do isph = 1, ddx_data % params % nsph
    do i = 1, 3
        ! Set the initial sums to zero
        sum_A_minus_h = zero
        sum_A_plus_h = zero
        sum_B_minus_h = zero
        sum_B_plus_h = zero
        sum_Ui_minus_h = zero
        sum_Ui_plus_h = zero
        sum_C_minus_h = zero
        sum_C_plus_h = zero
        ! Set the centers to x + h
        ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) &
                                           & + step
        ! Call solve
        call test_solve(ddx_data, sum_A_plus_h, sum_B_plus_h, sum_Ui_plus_h, &
                   & sum_C_plus_h)
        ! Set the center to x - h
        ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) &
                                           & - two*step
        ! Call solve
        call test_solve(ddx_data, sum_A_minus_h, sum_B_minus_h, sum_Ui_minus_h, &
                   & sum_C_minus_h)
        ! Set the center to x
        ddx_data % params % csph(i, isph) = ddx_data % params % csph(i, isph) &
                                          & + step
        ! Numerical derivative = (f(x+h)-f(x-h))/(2*h)
        derivative_num_A(i, isph) = (sum_A_plus_h - sum_A_minus_h) / two / step
        derivative_num_B(i, isph) = (sum_B_plus_h - sum_B_minus_h) / two / step
        derivative_num_Ui(i, isph)= (sum_Ui_plus_h - sum_Ui_minus_h) / two / step
        derivative_num_C(i, isph) = (sum_C_plus_h  - sum_C_minus_h) / two / step
    end do
end do
! Relative Errors
relerr_A = dnrm2(3*ddx_data % params % nsph,derivative_num_A-derivative_A,1) / &
    & dnrm2(3*ddx_data % params % nsph, derivative_A, 1)
relerr_B = dnrm2(3*ddx_data % params % nsph,derivative_num_B-derivative_B,1) / &
    & dnrm2(3*ddx_data % params % nsph, derivative_B, 1)
relerr_Ui = dnrm2(3*ddx_data % params % nsph, derivative_num_Ui - derivative_Ui, 1) / &
    & dnrm2(3*ddx_data % params % nsph, derivative_Ui, 1)
relerr_C1_C2 = dnrm2(3*ddx_data % params % nsph, &
    & derivative_num_C - (derivative_C1_C2), 1) / &
    & dnrm2(3*ddx_data % params % nsph, (derivative_C1_C2), 1)


write(6,'(2A60)') 'Matrix A Analytical forces', 'Numerical forces'
do i = 1, ddx_data % params % nsph
  write(6,'(6E20.10)') derivative_A(1,i), derivative_A(2,i), derivative_A(3,i), &
  & derivative_num_A(1,i), derivative_num_A(2,i), derivative_num_A(3,i)
end do

write(6,'(2A60)') ' Matrix B Analytical forces', 'Numerical forces'
do i = 1, ddx_data % params % nsph
  write(6,'(6E20.10)') derivative_B(1,i), derivative_B(2,i), derivative_B(3,i), &
  & derivative_num_B(1,i), derivative_num_B(2,i), derivative_num_B(3,i)
end do

write(6,'(2A60)') ' U_i^e(x_in) Analytical forces', 'Numerical forces'
do i = 1, ddx_data % params % nsph
  write(6,'(6E20.10)') derivative_Ui(1,i), derivative_Ui(2,i), derivative_Ui(3,i), &
  & derivative_num_Ui(1,i), derivative_num_Ui(2,i), derivative_num_Ui(3,i)
end do

write(6,'(2A60)') ' Matrix C1_C2 Analytical forces', 'Numerical forces'
do i = 1, ddx_data % params % nsph
  write(6,'(6E20.10)') derivative_C1_C2(1,i), &
                     & derivative_C1_C2(2,i), &
                     & derivative_C1_C2(3,i), &
                     & derivative_num_C(1,i), derivative_num_C(2,i), &
                     & derivative_num_C(3,i)
end do


! Deallocation
deallocate(derivative_num_A, derivative_num_B, derivative_num_Ui, &
           & vcos, vsin, vplm, basloc, dbsloc, &
           & random_vector_nbasis_nsph_two, &
           & random_vector_nbasis_nsph_one, &
           & random_vector_two_evaluated_at_grid, &
           & derivative_C1_C2, &
           & derivative_A, derivative_B, diff_re, derivative_num_C, charges)
call deallocate_model(ddx_data, ddx_error)

write(*, *) "Rel. error of A     :", relerr_A
write(*, *) "Rel. error of B     :", relerr_B
write(*, *) "Rel. error of U_i^e :", relerr_Ui
write(*, *) "Rel. error of C1_C2 :", relerr_C1_C2

if (relerr_A .gt. 1d-6) then
    write(*,*) 'Error in computing derivatives of A, Rel.Error : ', relerr_A
    stop 1
endif
if (relerr_B .gt. 1d-6) then
    write(*,*) 'Error in computing derivatives of B, Rel.Error : ', relerr_B
    stop 1
endif
if (relerr_Ui .gt. 1d-6) then
    write(*,*) 'Error in computing derivatives of Ui, Rel.Error : ', relerr_Ui
    stop 1
endif
if (relerr_C1_C2 .gt. 1d-6) then
    write(*,*) 'Error in computing derivatives of C1C2, Rel.Error : ',&
                & relerr_C1_C2
    stop 1
endif

contains

subroutine test_solve(ddx_data, sum_der_A, sum_der_B, sum_der_Ui, sum_der_C1_C2)
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(inout) :: sum_der_A, sum_der_B, sum_der_Ui, sum_der_C1_C2
    ! Local variables
    ! ddx_data2  : New ddx_data with new coordinates
    type(ddx_type) :: ddx_data2
    ! vector_n           : Unit vector of size n\times 1
    ! vector_cosmo            : y = Ax
    ! vector_lpb              : y = Bx
    ! vector_nbasis_nsph : Unit vector of size nbasis \times nsph
    ! vector_e_c1_c2          : y = C1Xr + C2Xe
    ! vector_r_c1_c2          : y = C1Xr + C2Xe
    real(dp), allocatable :: random_vector_n_one(:), &
                             & random_vector_n_two(:), &
                             & random_vector_nbasis_nsph_one(:,:), &
                             & random_vector_nbasis_nsph_two(:,:), &
                             & random_vector_C_one(:,:, :), &
                             & vector_cosmo(:), vector_lpb(:),&
                             & vector_c1_c2(:,:, :),&
                             & zero_vector(:,:)
    ! i      : Index for n
    ! isph   : Index for number of sphers
    ! igrid  : Index for grid points
    ! ibasis : Index for number of basis
    integer :: i, isph, igrid, ibasis
    type(ddx_error_type) :: ddx_error

    ! Initialise new ddx_data with new centers coordinates
    call allocate_model(ddx_data % params % nsph, &
        & ddx_data % params % csph(1, :), &
        & ddx_data % params % csph(2, :), &
        & ddx_data % params % csph(3, :), &
        & ddx_data % params % rsph, &
        & ddx_data % params % model, &
        & ddx_data % params % lmax, &
        & ddx_data % params % ngrid, 0, &
        & ddx_data % params % fmm, &
        & ddx_data % params % pm, &
        & ddx_data % params % pl, &
        & ddx_data % params % se, &
        & ddx_data % params % eta, &
        & ddx_data % params % eps, &
        & ddx_data % params % kappa, 0, &
        & ddx_data % params % maxiter, &
        & ddx_data % params % jacobi_ndiis, &
        & ddx_data % params % nproc, dummy_file_name, ddx_data2, ddx_error)
    ! Allocation
    allocate(random_vector_n_one(ddx_data2 % constants % n), &
             & random_vector_n_two(ddx_data2 % constants % n), &
             & random_vector_nbasis_nsph_one(ddx_data2 % constants % nbasis,&
             & ddx_data2 % params % nsph),&
             & random_vector_nbasis_nsph_two(ddx_data2 % constants % nbasis,&
             & ddx_data2 % params % nsph),&
             & random_vector_C_one(ddx_data2 % constants % nbasis,&
             & ddx_data2 % params % nsph, 2),&
             & vector_cosmo(ddx_data2 % constants % n), &
             & vector_lpb(ddx_data2 % constants % n), &
             & vector_c1_c2(ddx_data2 % constants % nbasis, &
             &  ddx_data2 % params % nsph, 2), &
             & zero_vector(ddx_data2 % constants%  nbasis, &
             & ddx_data2 % params % nsph))
    ! Intialisation
    random_vector_n_one = one
    random_vector_n_two = one
    random_vector_n_one(1) = two
    random_vector_n_two(2) = three

    random_vector_nbasis_nsph_one = one
    random_vector_nbasis_nsph_two = one
    random_vector_nbasis_nsph_one(1,1) = two
    random_vector_nbasis_nsph_two(2,1) = three

    random_vector_C_one(:,:,1) = random_vector_nbasis_nsph_one
    random_vector_C_one(:,:,2) = random_vector_nbasis_nsph_one


    vector_cosmo = one
    vector_lpb = one
    vector_c1_c2 = zero
    zero_vector = zero
    ! Call for matrix A
    call lx(ddx_data2 % params, ddx_data2 % constants, &
          & ddx_data2 % workspace, random_vector_n_one, vector_cosmo, ddx_error)
    ! Call for matrix B
    call bx(ddx_data2 % params, ddx_data2 % constants, &
              & ddx_data2 % workspace, &
              & random_vector_n_one, vector_lpb, ddx_error)
    call cx(ddx_data2 % params, ddx_data2 % constants, &
                 & ddx_data2 % workspace, &
                 & random_vector_C_one, &
                 & vector_c1_c2, ddx_error)
    ! Sum for U_i^e(x_in)
    do isph = 1,ddx_data2 %  params % nsph
      do igrid = 1,ddx_data2 %  params % ngrid
        if (ddx_data2 %  constants % ui(igrid, isph) .gt. zero) then
          sum_der_Ui = sum_der_Ui + ddx_data2 %  constants % wgrid(igrid) * &
                        & ddx_data2 %  constants % ui(igrid, isph)
        end if
      end do
    end do
    ! Sum for matrix C1_C2 with Q constant
    do ibasis = 1, ddx_data2 % constants % nbasis
      do isph = 1, ddx_data2 %  params % nsph
        sum_der_C1_C2 = sum_der_C1_C2 &
                    & + random_vector_nbasis_nsph_two(ibasis, isph)* &
                    & vector_c1_c2(ibasis, isph, 1)&
                    & + random_vector_nbasis_nsph_two(ibasis, isph)* &
                    & vector_c1_c2(ibasis, isph, 2)
      end do
    end do
    ! Sum for matrix A and matrix B
    do i = 1, ddx_data2 % constants % n
        sum_der_A = sum_der_A + random_vector_n_two(i)*vector_cosmo(i)
        sum_der_B = sum_der_B + random_vector_n_two(i)*vector_lpb(i)
    end do
    ! Deallocation
    deallocate(random_vector_n_one, vector_cosmo, &
               & random_vector_n_two, &
               & random_vector_nbasis_nsph_one, &
               & random_vector_nbasis_nsph_two, &
               & random_vector_C_one, &
               & vector_lpb, vector_c1_c2)
end subroutine test_solve

end program main
