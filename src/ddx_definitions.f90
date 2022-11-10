!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_definitions.f90
!! Compile-time constants are defined here.
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-06-15

!> Compile-time constants and definitions
module ddx_definitions

! Disable implicit types
implicit none

!> Kind of double precision
integer, parameter :: dp = kind(1.0d0)

!! Compile-time constants
real(dp), parameter :: zero = 0d0, one = 1d0, two = 2d0, three = 3d0
real(dp), parameter :: four = 4d0, pt5 = 5d-1
real(dp), parameter :: sqrt2 = sqrt(two)
real(dp), parameter :: sqrt3 = sqrt(three)
real(dp), parameter :: pi4 = atan(one)
real(dp), parameter :: pi = four * pi4
real(dp), parameter :: fourpi = four * pi
real(dp), parameter :: twopi = two * pi
real(dp), parameter :: sqrt4pi = four * sqrt(pi4)
real(dp), parameter :: machine_eps = epsilon(zero)
real(dp), parameter :: toang = 0.52917721092d0
real(dp), parameter :: tokcal = 627.509469d0
real(dp), parameter :: tokj = 2625.509469d0
real(dp), parameter :: tobohr = one / toang
!> Number of supported Lebedev grids
integer, parameter :: nllg = 32
!> Number of grid points of each Lebedev grid
integer, parameter :: ng0(nllg) = (/ 6, 14, 26, 38, 50, 74, 86, 110, 146, &
    & 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, &
    & 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 /)
!> Names of ddX models
character(len=255), parameter :: model_str(3) = (/ "COSMO", "PCM  ", "LPB  " /)

end module ddx_definitions

