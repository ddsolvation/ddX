!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! Test of adjoint for each matrices
!!
!! @version 1.0.0
!! @author Abhinav Jha
!! @date 2022-03-29

program main
use ddx_core
use ddx_operators
use ddx_lpb
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
integer :: info
! isph   : Index for number of spheres
! i      : Index for derivative components (i = 1,2,3)
! ibasis : Index for number of basis
integer :: isph, i, ibasis, iprint


real(dp) :: check_A_one, check_A_two, check_A_three, check_A_four, tol
real(dp) :: check_A_star_one, check_A_star_two, check_A_star_three,&
          & check_A_star_four
real(dp) :: check_B_one, check_B_two, check_B_three, check_B_four
real(dp) :: check_B_star_one, check_B_star_two, check_B_star_three,&
          & check_B_star_four
real(dp) :: check_C_one, check_C_two
real(dp) :: check_C_star_one, check_C_star_two
real(dp), allocatable:: random_vector_nbasis_nsph_one(:,:), &
                        & random_vector_nbasis_nsph_two(:,:), &
                        & random_vector_nbasis_nsph_three(:,:), &
                        & random_vector_nbasis_nsph_four(:,:), &
                        & zero_vector(:,:)


real(dp), allocatable :: random_vector_n_one(:), &
                       & random_vector_n_two(:), &
                       & random_vector_n_three(:), &
                       & random_vector_n_four(:), &
                       & random_vector_C_one(:,:,:),&
                       & random_vector_C_two(:,:,:),&
                       & vector_A_one(:), &
                       & vector_A_two(:), &
                       & vector_A_three(:), &
                       & vector_A_four(:), &
                       & vector_B_one(:), &
                       & vector_B_two(:), &
                       & vector_B_three(:), &
                       & vector_B_four(:), &
                       & vector_C1_one(:,:), &
                       & vector_C1_two(:,:), &
                       & vector_C2_one(:,:), &
                       & vector_C2_two(:,:), &
                       & vector_C_one(:,:,:), &
                       & vector_C_two(:,:,:), &
                       & vector_A_star_one(:), &
                       & vector_A_star_two(:), &
                       & vector_A_star_three(:), &
                       & vector_A_star_four(:), &
                       & vector_B_star_one(:), &
                       & vector_B_star_two(:), &
                       & vector_B_star_three(:), &
                       & vector_B_star_four(:), &
                       & vector_C1_star_one(:,:), &
                       & vector_C1_star_two(:,:), &
                       & vector_C2_star_one(:,:), &
                       & vector_C2_star_two(:,:), &
                       & vector_C_star_one(:,:,:), &
                       & vector_C_star_two(:,:,:)
real(dp) :: lmax0, nbasis0

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file 12"
call ddfromfile(fname, ddx_data, tol, iprint, info)
if(info .ne. 0) stop "info != 0"

! lmax0 set to minimum of 6 or given lmax.
! nbasis0 set to minimum of 49 or given (lmax+1)^2.
! Maybe do it ddlpb_init(?)
lmax0 = MIN(6, ddx_data % params % lmax)
nbasis0 = MIN(49, ddx_data % constants % nbasis)

! Allocation for various vectors
allocate(random_vector_nbasis_nsph_one(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & random_vector_nbasis_nsph_two(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & random_vector_nbasis_nsph_three(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & random_vector_nbasis_nsph_four(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & random_vector_C_one(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph, 2), &
    & random_vector_C_two(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph, 2), &
    & random_vector_n_one(ddx_data % constants % n), &
    & random_vector_n_two(ddx_data % constants % n), &
    & random_vector_n_three(ddx_data % constants % n), &
    & random_vector_n_four(ddx_data % constants % n), &
    & vector_A_one(ddx_data % constants % n), &
    & vector_A_two(ddx_data % constants % n), &
    & vector_A_three(ddx_data % constants % n), &
    & vector_A_four(ddx_data % constants % n), &
    & vector_B_one(ddx_data % constants % n), &
    & vector_B_two(ddx_data % constants % n), &
    & vector_B_three(ddx_data % constants % n), &
    & vector_B_four(ddx_data % constants % n), &
    & vector_C1_one(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & vector_C1_two(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & vector_C2_one(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & vector_C2_two(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & vector_C_one(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph, 2), &
    & vector_C_two(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph, 2), &
    & vector_A_star_one(ddx_data % constants % n), &
    & vector_A_star_two(ddx_data % constants % n), &
    & vector_A_star_three(ddx_data % constants % n), &
    & vector_A_star_four(ddx_data % constants % n), &
    & vector_B_star_one(ddx_data % constants % n), &
    & vector_B_star_two(ddx_data % constants % n), &
    & vector_B_star_three(ddx_data % constants % n), &
    & vector_B_star_four(ddx_data % constants % n), &
    & vector_C1_star_one(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & vector_C1_star_two(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & vector_C2_star_one(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & vector_C2_star_two(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph), &
    & vector_C_star_one(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph, 2), &
    & vector_C_star_two(ddx_data % constants % nbasis, &
    & ddx_data % params % nsph, 2), &
    & zero_vector(ddx_data % constants %  nbasis, ddx_data % params % nsph))

! Initialise
vector_A_one = one
vector_A_two = one
vector_A_three = one
vector_A_four = one
vector_B_one = one
vector_B_two = one
vector_B_three = one
vector_B_four = one
vector_C1_one = one
vector_C1_two = one
vector_C2_one = one
vector_C2_two = one
vector_C_one = one
vector_C_two = one
vector_A_star_one = one
vector_A_star_two = one
vector_A_star_three = one
vector_A_star_four = one
vector_B_star_one = one
vector_B_star_two = one
vector_B_star_three = one
vector_B_star_four = one
vector_C1_star_one = one
vector_C1_star_two = one
vector_C2_star_one = one
vector_C2_star_two = one
vector_C_star_one = one
vector_C_star_two = one
zero_vector = zero

check_A_one = zero
check_A_two = zero
check_A_three = zero
check_A_four = zero
check_B_one = zero
check_B_two = zero
check_B_three = zero
check_B_four = zero
check_C_one = zero
check_C_two = zero

check_A_star_one = zero
check_A_star_two = zero
check_A_star_three = zero
check_A_star_four = zero
check_B_star_one = zero
check_B_star_two = zero
check_B_star_three = zero
check_B_star_four = zero
check_C_star_one = zero
check_C_star_two = zero

! Random vectors
random_vector_nbasis_nsph_one = one
random_vector_nbasis_nsph_two = one
random_vector_nbasis_nsph_three = one
random_vector_nbasis_nsph_four = one

random_vector_n_one = one
random_vector_n_two = one
random_vector_n_three = one
random_vector_n_four = one

! random_vector_n_one   = [2,1,1,1,...]
! random_vector_n_two   = [1,3,1,1,...]
! random_vector_n_three = [1,1,0,1,...]
! random_vector_n_four  = [0,0,0,1,...]
random_vector_n_one(1) = two
random_vector_n_two(2) = three
random_vector_n_three(3) = zero
random_vector_n_four(1) = zero
random_vector_n_four(2) = zero
random_vector_n_four(3) = zero

random_vector_nbasis_nsph_one(1,1) = 2.0
random_vector_nbasis_nsph_two(2,1) = 3.0
random_vector_nbasis_nsph_three(3,1) = 0.0
random_vector_nbasis_nsph_four(1,1) = 0.0
random_vector_nbasis_nsph_four(2,1) = 0.0
random_vector_nbasis_nsph_four(3,1) = 0.0

random_vector_C_one(:,:,1) = random_vector_nbasis_nsph_one(:,:)
random_vector_C_one(:,:,2) = random_vector_nbasis_nsph_two(:,:)
random_vector_C_two(:,:,1) = random_vector_nbasis_nsph_three(:,:)
random_vector_C_two(:,:,2) = random_vector_nbasis_nsph_four(:,:)

! Call for matrix A
call lx(ddx_data % params, ddx_data % constants, &
          & ddx_data % workspace, random_vector_n_one, vector_A_one)
call lx(ddx_data % params, ddx_data % constants, &
          & ddx_data % workspace, random_vector_n_two, vector_A_two)
call lx(ddx_data % params, ddx_data % constants, &
          & ddx_data % workspace, random_vector_n_three, vector_A_three)
call lx(ddx_data % params, ddx_data % constants, &
          & ddx_data % workspace, random_vector_n_four, vector_A_four)


! Call for matrix B
call bx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_one, vector_B_one)
call bx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_two, vector_B_two)
call bx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_three, vector_B_three)
call bx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_four, vector_B_four)

! Call for C1 and C2
call cx(ddx_data % params, ddx_data % constants, &
                 & ddx_data % workspace, &
                 & random_vector_C_one, &
                 & vector_C_one)

call cx(ddx_data % params, ddx_data % constants, &
                 & ddx_data % workspace, &
                 & random_vector_C_two, &
                 & vector_C_two)
! Call for matrix Astar
call lstarx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_four, vector_A_star_one)
call lstarx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_three, vector_A_star_two)
call lstarx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_two, vector_A_star_three)
call lstarx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_one, vector_A_star_four)

! Call for matrix Bstar
call bstarx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_four, vector_B_star_one)
call bstarx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_three, vector_B_star_two)
call bstarx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_two, vector_B_star_three)
call bstarx(ddx_data % params, ddx_data % constants, &
      & ddx_data % workspace, &
      & random_vector_n_one, vector_B_star_four)

! Call for C1 and C2 star
! |C1* C1*||X3|
! |C2* C2*||X4|
call cstarx(ddx_data % params, ddx_data % constants, &
                 & ddx_data % workspace, &
                 & random_vector_C_two, &
                 & vector_C_star_one)

call cstarx(ddx_data % params, ddx_data % constants, &
                 & ddx_data % workspace, &
                 & random_vector_C_one, &
                 & vector_C_star_two)

!Compute the contraction
do i = 1, ddx_data % constants % n
  check_A_one = check_A_one + random_vector_n_four(i)*vector_A_one(i)
  check_A_two = check_A_two + random_vector_n_three(i)*vector_A_two(i)
  check_A_three = check_A_three + random_vector_n_two(i)*vector_A_three(i)
  check_A_four = check_A_four + random_vector_n_one(i)*vector_A_four(i)
  check_A_star_one = check_A_star_one + &
                   & random_vector_n_one(i)*vector_A_star_one(i)
  check_A_star_two = check_A_star_two + &
                  & random_vector_n_two(i)*vector_A_star_two(i)
  check_A_star_three = check_A_star_three + &
                  & random_vector_n_three(i)*vector_A_star_three(i)
  check_A_star_four = check_A_star_four + &
                  & random_vector_n_four(i)*vector_A_star_four(i)

  check_B_one = check_B_one + random_vector_n_four(i)*vector_B_one(i)
  check_B_two = check_B_two + random_vector_n_three(i)*vector_B_two(i)
  check_B_three = check_B_three + random_vector_n_two(i)*vector_B_three(i)
  check_B_four = check_B_four + random_vector_n_one(i)*vector_B_four(i)
  check_B_star_one = check_B_star_one + &
                   & random_vector_n_one(i)*vector_B_star_one(i)
  check_B_star_two = check_B_star_two + &
                  & random_vector_n_two(i)*vector_B_star_two(i)
  check_B_star_three = check_B_star_three + &
                  & random_vector_n_three(i)*vector_B_star_three(i)
  check_B_star_four = check_B_star_four + &
                  & random_vector_n_four(i)*vector_B_star_four(i)
end do

do ibasis = 1, ddx_data % constants % nbasis
  do isph = 1, ddx_data % params % nsph
    check_C_one = check_C_one + random_vector_C_two(ibasis, isph, 1) &
                & * vector_C_one(ibasis, isph, 1) + &
                & random_vector_C_two(ibasis, isph, 2) &
                & * vector_C_one(ibasis, isph, 2)
    check_C_two = check_C_two + random_vector_C_one(ibasis, isph, 1) &
                & * vector_C_two(ibasis, isph, 1) + &
                & random_vector_C_one(ibasis, isph, 2) &
                & * vector_C_two(ibasis, isph, 2)
    check_C_star_one = check_C_star_one + &
                & random_vector_C_one(ibasis, isph, 1) &
                & * vector_C_star_one(ibasis, isph, 1) + &
                & random_vector_C_one(ibasis, isph, 2) &
                & * vector_C_star_one(ibasis, isph, 2)
    check_C_star_two = check_C_star_two + &
                & random_vector_C_two(ibasis, isph, 1) &
                & * vector_C_star_two(ibasis, isph, 1) + &
                & random_vector_C_two(ibasis, isph, 2) &
                & * vector_C_star_two(ibasis, isph, 2)
  end do
end do

! Deallocation
deallocate(random_vector_nbasis_nsph_two, &
           & random_vector_nbasis_nsph_one, &
           & random_vector_nbasis_nsph_three, &
           & random_vector_nbasis_nsph_four, &
           & random_vector_n_one, &
           & random_vector_n_two, &
           & random_vector_n_three, &
           & random_vector_n_four, &
           & random_vector_C_one, &
           & random_vector_C_two, &
           & vector_A_one, &
           & vector_A_two, &
           & vector_A_three, &
           & vector_A_four, &
           & vector_B_one, &
           & vector_B_two, &
           & vector_B_three, &
           & vector_B_four, &
           & vector_A_star_one, &
           & vector_A_star_two, &
           & vector_A_star_three, &
           & vector_A_star_four, &
           & vector_B_star_one, &
           & vector_B_star_two, &
           & vector_B_star_three, &
           & vector_B_star_four, &
           & vector_C1_one, &
           & vector_C1_two, &
           & vector_C2_one, &
           & vector_C2_two, &
           & vector_C_one, &
           & vector_C_two, &
           & vector_C1_star_one, &
           & vector_C1_star_two, &
           & vector_C2_star_one, &
           & vector_C2_star_two, &
           & vector_C_star_one, &
           & vector_C_star_two, &
           & zero_vector)
call ddfree(ddx_data)

write(*, *) "y4(A)x1  :", check_A_one, ", x1(A*)y4  :", check_A_star_one
write(*, *) "y3(A)x2  :", check_A_two, ", x2(A*)y3  :", check_A_star_two
write(*, *) "y2(A)x3  :", check_A_three, ", x3(A*)y2  :", check_A_star_three
write(*, *) "y1(A)x4  :", check_A_four, ", x4(A*)y1  :", check_A_star_four
write(*, *) "y4(B)x1  :", check_B_one, ", x1(B*)y4  :", check_B_star_one
write(*, *) "y3(B)x2  :", check_B_two, ", x2(B*)y3  :", check_B_star_two
write(*, *) "y2(B)x3  :", check_B_three, ", x3(B*)y2  :", check_B_star_three
write(*, *) "y1(B)x4  :", check_B_four, ", x4(B*)y1  :", check_B_star_four
write(*, *) "y1(C)x1  :", check_C_one, ", x1(C*)y1  :", check_C_star_one
write(*, *) "y2(C)x2  :", check_C_two, ", x2(C*)y2  :", check_C_star_two

if(abs(check_A_one - check_A_star_one) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix A not equal! y(A)x : ', check_A_one, &
            & ', x(A*)y : ', check_A_star_one
  stop 1
endif


if(abs(check_A_two - check_A_star_two) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix A not equal! y(A)x : ', check_A_two, &
            & ', x(A*)y : ', check_A_star_two
  stop 1
endif


if(abs(check_A_three - check_A_star_three) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix A not equal! y(A)x : ', check_A_three, &
            & ', x(A*)y : ', check_A_star_three
  stop 1
endif


if(abs(check_A_four - check_A_star_four) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix A not equal! y(A)x : ', check_A_four, &
            & ', x(A*)y : ', check_A_star_four
  stop 1
endif


if(abs(check_B_one - check_B_star_one) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix B not equal! y(B)x : ', check_B_one, &
            & ', x(B*)y : ', check_B_star_one
  stop 1
endif


if(abs(check_B_two - check_B_star_two) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix B not equal! y(A)x : ', check_B_two, &
            & ', x(B*)y : ', check_B_star_two
  stop 1
endif


if(abs(check_B_three - check_B_star_three) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix B not equal! y(C)x : ', check_B_three, &
            & ', x(B*)y : ', check_B_star_three
  stop 1
endif


if((check_B_four - check_B_star_four) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix B not equal! y(B)x : ', check_B_four, &
            & ', x(B*)y : ', check_B_star_four
  stop 1
endif

if((check_C_one - check_C_star_one) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix C not equal! y(C)x : ', check_C_one, &
            & ', x(C*)y : ', check_C_star_one
  stop 1
endif

if(abs(check_C_two - check_C_star_two) .gt. 1e-12) then
  write(*,*) 'Contraction for matrix C not equal! y(C)x : ', check_C_two, &
            & ', x(C*)y : ', check_C_star_two
  stop 1
endif



end program main
