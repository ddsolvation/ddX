!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/ddlpb_esolv_iterations.f90
!! Test for solvation energy and number of outer iterations for different
!! values of input parameter
!! NOTE: This test is not a definitive test. If the values are different
!!       maybe the new ones are correct as they might be improved values.
!!       But, the default values should be used as a benchmark
!!
!! @version 1.0.0
!! @author Abhinav Jha
!! @date 2022-03-28

program main
use ddx_core
use ddx_operators
use ddx_solvers
use ddx
use ddx_lpb
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
type(ddx_error_type) :: error

real(dp) :: esolv, default_value, tol
integer :: i, istatus, default_lmax_val, n_iter
real(dp), allocatable :: default_epsilon(:), default_eta(:), &
                       & default_kappa(:), default_lmax(:)
integer, allocatable :: default_iter_epsilon(:), default_iter_eta(:), &
                       & default_iter_kappa(:), default_iter_lmax(:)
real(dp), allocatable :: charges(:)

real(dp), external :: dnrm2
character(len=255) :: dummy_file_name = ''

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, ddx_data, tol, charges, error)
call check_error(error)

! Allocation for variable vectors
! default_"variable_name" : These are the precomputed values
! computed_"variable_name" : These are the computed values
allocate(default_epsilon(4), default_eta(4), &
       & default_kappa(4), default_lmax(4), &
       & default_iter_epsilon(4), default_iter_eta(4), &
       & default_iter_kappa(4), default_iter_lmax(4), stat=istatus)

if (istatus.ne.0) write(6,*) 'Allocation failed'

!Default values precomputed
!epsilon_solv : 2, 20, 200, 2000

default_epsilon = (/ -5.3280230267698165E-004, -9.7406452041931936E-004, &
                   & -1.0243211234919017E-003, -1.0294180585685288E-003 /)

! obsolete case: lmax0 = min(6, lmax), ngrid ~ 200
! default_epsilon = (/ -5.3518110117345332E-004, -9.7393853923451006E-004, &
!                    & -1.0237954253809903E-003, -1.0288501533655906E-003 /)

!eta : 0.0001, 0.001, 0.01, 0.1
default_eta = (/ -1.0151699327713233E-003, -1.0152224349631500E-003, &
               & -1.0153111171958898E-003, -1.0155598899696948E-003 /)

! obsolete case: lmax0 = min(6, lmax), ngrid ~ 200
! default_eta = (/ -1.0144763156304426E-003, -1.0144763156304426E-003, &
!                & -1.0144761325115853E-003, -1.0151060052969220E-003 /)

!kappa : 0.5, 0.25, 0.16667, 0.125
default_kappa = (/ -1.0233753433615445E-003, -1.0197139938641997E-003, &
                 & -1.0175938367002952E-003, -1.0162805489436332E-003 /)

! obsolete case: lmax0 = min(6, lmax), ngrid ~ 200
! default_kappa = (/ -1.0228739915625511E-003, -1.0192321310712657E-003,&
!                  & -1.0171251065945882E-003, -1.0158211256043079E-003 /)

!lmax : 2, 4, 8, 16
default_lmax = (/ -1.0000499299268406E-003, -1.0131496803324205E-003, &
                & -1.0163417969913872E-003, -1.0180392518042718E-003 /)

! obsolete case: lmax0 = min(6, lmax), ngrid ~ 200
! default_lmax = (/ -9.9977268971430206E-004, -1.0128441259120659E-003, &
!               & -1.0157913843224611E-003, -1.0177420746553952E-003 /)

!Default values precomputed
!epsilon_solv : 2, 20, 200, 2000
default_iter_epsilon = (/ 7, 8, 7, 7 /)

!eta : 0.0001, 0.001, 0.01, 0.1
default_iter_eta = (/ 8, 8, 8, 8 /)
!kappa : 0.5, 0.25, 0.16667, 0.125
default_iter_kappa = (/ 7, 7, 8, 8 /)
!lmax : 2, 4, 8, 16
default_iter_lmax = (/ 7, 7, 8, 8 /)

! Initial values
esolv = zero
n_iter = zero
! Computation for different eps_solv
write(*,*) 'Varying values of epsilon_solv'
do i = 1, 4
  default_value = 0.2*(10**i)
  write(*,*) 'epsilon_solv : ', default_value
  esolv = zero
  call solve(ddx_data, esolv, n_iter, default_value, &
           & ddx_data % params % eta, ddx_data % params % kappa, &
           & ddx_data % params % lmax, tol, charges)
  call check_values(default_epsilon(i), esolv)
  call check_iter_values(default_iter_epsilon(i), n_iter)
end do


! Computation for different eta
write(*,*) 'Varying values of eta'
do i = 1, 4
  default_value = 0.00001*(10**i)
  write(*,*) 'eta : ', default_value
  esolv = zero
  call solve(ddx_data, esolv, n_iter, ddx_data % params % eps, &
           & default_value, ddx_data % params % kappa, &
           & ddx_data % params % lmax, tol, charges)
  call check_values(default_eta(i), esolv)
  call check_iter_values(default_iter_eta(i), n_iter)
end do

! Computation for different kappa
write(*,*) 'Varying values of kappa'
do i = 1, 4
  default_value = 1.0/(2.0*i)
  write(*,*) 'kappa : ', default_value
  esolv = zero
  call solve(ddx_data, esolv, n_iter, ddx_data % params % eps, &
           & ddx_data % params % eta, default_value, &
           & ddx_data % params % lmax, tol, charges)
  call check_values(default_kappa(i), esolv)
  call check_iter_values(default_iter_kappa(i), n_iter)
end do

! Computation for different lmax
write(*,*) 'Varying values of lmax'
do i = 1, 4
  default_lmax_val = 2**i
  write(*,*) 'lmax : ', default_lmax_val
  esolv = zero
  call solve(ddx_data, esolv, n_iter, ddx_data % params % eps, &
           & ddx_data % params % eta, &
           & ddx_data % params % kappa, default_lmax_val, tol, charges)
  call check_values(default_lmax(i), esolv)
  call check_iter_values(default_iter_lmax(i), n_iter)
end do

deallocate(default_epsilon, default_eta, &
       & default_kappa, default_lmax, &
       & default_iter_epsilon, default_iter_eta, &
       & default_iter_kappa, default_iter_lmax, charges, stat = istatus)

if (istatus.ne.0) write(6,*) 'Deallocation failed'

call ddfree(ddx_data, error)

contains

subroutine solve(ddx_data, esolv_in, n_iter, epsilon_solv, eta, kappa, lmax, tol, charges)
    implicit none
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(in) :: charges(ddx_data % params % nsph)
    real(dp), intent(inout) :: esolv_in
    integer, intent(inout)  :: n_iter
    real(dp), intent(in) :: epsilon_solv
    real(dp), intent(in) :: eta
    real(dp), intent(in) :: kappa
    integer, intent(in)  :: lmax
    real(dp), intent(in) :: tol

    type(ddx_state_type) :: state
    type(ddx_type) :: ddx_data2
    type(ddx_error_type) :: error2
    real(dp), allocatable :: phi_cav2(:)
    real(dp), allocatable :: gradphi_cav2(:,:)
    real(dp), allocatable :: hessianphi_cav2(:,:,:)
    real(dp), allocatable :: psi2(:,:)
    real(dp), allocatable :: force2(:,:)

    call ddinit(ddx_data % params % nsph, &
        & ddx_data % params % csph(1, :), &
        & ddx_data % params % csph(2, :), ddx_data % params % csph(3, :), ddx_data % params % rsph, &
        & ddx_data % params % model, lmax, ddx_data % params % ngrid, 0, &
        & ddx_data % params % fmm, ddx_data % params % pm, ddx_data % params % pl, &
        & ddx_data % params % se, &
        & eta, epsilon_solv, kappa, 0,&
        & ddx_data % params % maxiter, &
        & ddx_data % params % jacobi_ndiis, &
        & ddx_data % params % nproc, dummy_file_name, ddx_data2, error2)
    call check_error(error2)

    ! the state depends on lmax, so it is allocated here
    call allocate_state(ddx_data2 % params, ddx_data2 % constants, state, error2)
    call check_error(error2)

    allocate(phi_cav2(ddx_data2 % constants % ncav), gradphi_cav2(3, ddx_data2 % constants % ncav), &
            & hessianphi_cav2(3, 3, ddx_data2 % constants % ncav), &
            & psi2(ddx_data2 % constants % nbasis, ddx_data2 % params % nsph), &
            & force2(3, ddx_data2 % params % nsph))

    gradphi_cav2 = zero; phi_cav2 = zero
    hessianphi_cav2 = zero; psi2 = zero; force2 = zero

    call mkrhs(ddx_data2 % params, ddx_data2 % constants, ddx_data2 % workspace, &
        &  1, phi_cav2, 1, gradphi_cav2, 1, hessianphi_cav2, psi2, charges)
    gradphi_cav2 = - gradphi_cav2
    call ddsolve(ddx_data2, state, phi_cav2, gradphi_cav2, hessianphi_cav2, &
        & psi2, tol, esolv_in, force2, error2)
    call check_error(error2)

    n_iter = state % x_lpb_niter
    deallocate(phi_cav2, gradphi_cav2, hessianphi_cav2, psi2, force2)

    call deallocate_state(state, error2)
    call ddfree(ddx_data2, error2)
end subroutine solve

! This subroutine checks if the default and computed values are same
subroutine check_values(default_value, computed_value)
    real(dp), intent(in) :: default_value
    real(dp), intent(in) :: computed_value
    if(abs(default_value - computed_value) .gt. 1d-9) then
      write(*,*) 'Issue in value, default : ', default_value, ', computed : ',&
                & computed_value
      stop 1
    endif
end subroutine check_values

! This subroutine checks if the default and computed values are same
subroutine check_iter_values(default_value, computed_value)
    integer, intent(in) :: default_value
    integer, intent(in) :: computed_value
    if(default_value .ne. computed_value) then
      write(*,*) 'Issue in iterative value, default : ', default_value, ', computed : ',&
                & computed_value
      stop 1
    endif
end subroutine check_iter_values


end program main


