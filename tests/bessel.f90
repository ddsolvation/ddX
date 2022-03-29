!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/test_bessel.f90
!! Test of Bessel functions of first kind and second kind with varying arguments
!!
!! @version 1.0.0
!! @author Abhinav Jha
!! @date 2022-03-28

program main
use ddx_harmonics
implicit none

!! lmax_final :Maximal number of spherical harmonics, 2^lmax
integer :: lmax_final = 6, lmax = 0, iter_lmax = 0, l = 0, istatus = 1
!! argument_counter : Final value of argument, 10^argument_counter
integer :: iter_arg = 1, argument_counter = 1
!! Arguments for Bessel function of first kind and second kind
real(dp) :: argument_SI = 0.0 , argument_SI_plus_h = 0.0, argument_SI_minus_h = 0.0
real(dp) :: argument_SK = 0.0 , argument_SK_plus_h = 0.0, argument_SK_minus_h = 0.0
!! Step size : 10^(-6)
!! relerr_SI : Relative error for Bessel function of first kind
!! relerr_SK : Relative error for Bessel function of second kind
real(dp) :: step_size = 1d-7, relerr_SI = 0.0, relerr_SK = 0.0

real(dp), allocatable :: SI(:), DI(:), SK(:), DK(:)
real(dp), allocatable :: SI_plus_h(:), DI_plus_h(:), SK_plus_h(:), DK_plus_h(:)
real(dp), allocatable :: SI_minus_h(:), DI_minus_h(:), SK_minus_h(:), DK_minus_h(:)
real(dp), allocatable :: num_der_SI(:), num_der_SK(:)
!! Temporary workspace
complex(dp), allocatable :: work(:)

real(dp), external :: dnrm2

!Loop over for different values of lmax
do iter_lmax = 1, lmax_final
  !Resetting the values to zero
  argument_SI = 0.0
  argument_SI_plus_h = 0.0
  argument_SI_minus_h = 0.0

  argument_SK = 0.0
  argument_SK_plus_h = 0.0
  argument_SK_minus_h = 0.0

  lmax = 0
  !Loop over for different values of argument
  do iter_arg = -4, argument_counter
    !For Bessel function of first kind we only look at values less than 1,
    !because it acts on the inside of the sphere
    argument_SI = 10.0**(iter_arg)
    argument_SI_plus_h = argument_SI + step_size
    argument_SI_minus_h = argument_SI - step_size
    !For Bessel function of second kind we only look at values greater than 1,
    !because it acts on the outside of the sphere
    argument_SK = 10.0**(-iter_arg)
    argument_SK_plus_h = argument_SK + step_size
    argument_SK_minus_h = argument_SK - step_size

    !lmax = 2^n, n: 1,...,6
    lmax = 2**iter_lmax
    !Allocation of vectors
    allocate(SI(0:lmax), DI(0:lmax), SK(0:lmax), DK(0:lmax), &
       & SI_plus_h(0:lmax), DI_plus_h(0:lmax), &
       & SK_plus_h(0:lmax), DK_plus_h(0:lmax), &
       & SI_minus_h(0:lmax), DI_minus_h(0:lmax), &
       & SK_minus_h(0:lmax), DK_minus_h(0:lmax), &
       & num_der_SI(0:lmax), num_der_SK(0:lmax), &
       & work(max(2, lmax + 1)), stat=istatus)
    if ( istatus.ne.0 ) then
        write(*,*) 'Allocation failed !'
        stop
    endif

    write(6,*) 'lmax        : ', lmax
    write(6,*) 'argument_SI : ', argument_SI
    write(6,*) 'argument_SK : ', argument_SK
    write(6,*) 'h           : ', step_size
    !Calling at argument
    call modified_spherical_bessel_first_kind(lmax, &
                            & argument_SI, SI, DI, work)

    call modified_spherical_bessel_second_kind(lmax, &
                            & argument_SK, SK, DK, work)

    !Calling at argument + h
    call modified_spherical_bessel_first_kind(lmax, &
                            & argument_SI_plus_h, SI_plus_h, DI_plus_h, work)

    call modified_spherical_bessel_second_kind(lmax, &
                            & argument_SK_plus_h, SK_plus_h, DK_plus_h, work)

    !Calling at argument - h
    call modified_spherical_bessel_first_kind(lmax, &
                            & argument_SI_minus_h, SI_minus_h, DI_minus_h, work)

    call modified_spherical_bessel_second_kind(lmax, &
                            & argument_SK_minus_h, SK_minus_h, DK_minus_h, work)

    !Compute central difference numerical derivatives
    do l = 0, lmax
      num_der_SI(l) = (SI_plus_h(l)-SI_minus_h(l))/step_size/2

      num_der_SK(l) = (SK_plus_h(l)-SK_minus_h(l))/step_size/2
    end do

    !Output the derivatives
    write(6,*) 'l ', 'Analytical Derivatives (SI) ', 'Numerical Derivatives'
    do l = 0, lmax
      write(6,*) l, DI(l), num_der_SI(l)
    end do

     write(6,*) 'l ', 'Analytical Derivatives (SK) ', 'Numerical Derivatives'
    do l = 0, lmax
      write(6,*) l, DK(l), num_der_SK(l)
    end do

    !Relative error
    relerr_SI = dnrm2(lmax, num_der_SI - DI, 1) / dnrm2(lmax, DI, 1)
    relerr_SK = dnrm2(lmax, num_der_SK - DK, 1) / dnrm2(lmax, DK, 1)
    if (relerr_SI .gt. 1d-7) then
        write(*,*) 'Error in Bessel function of first kind'
        stop 1
    endif
    if (relerr_SK .gt. 1d-7) then
        write(*,*) 'Error in Bessel function of second kind'
        stop 1
    endif
    deallocate(SI, DI, SK, DK, &
       & SI_plus_h, DI_plus_h, &
       & SK_plus_h, DK_plus_h, &
       & SI_minus_h, DI_minus_h, &
       & SK_minus_h, DK_minus_h, &
       & num_der_SI, num_der_SK, &
       & work)
  end do
end do

end program main
