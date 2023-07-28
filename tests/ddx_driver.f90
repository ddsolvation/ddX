!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/ddx_driver.f90
!! Main ddx driver.
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

program main
use ddx_core
use ddx_operators
use ddx_solvers
use ddx
implicit none

character(len=255) :: finname, foutname, tmpstr
type(ddx_type) :: ddx_data
type(ddx_error_type) :: error
type(ddx_state_type) :: state
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), &
    & hessianphi_cav(:, :, :), psi(:, :), force(:, :), charges(:)
real(dp) :: tol, threshold, esolv, esolv2, fnorm, fdiff, ftmp(3)
integer :: i, isph, istatus
real(dp), external :: dnrm2

! Read input file name
call getarg(1, finname)
! Read output file name
call getarg(2, foutname)
! Read accuracy threshold
call getarg(3, tmpstr)
read(tmpstr, *) threshold
! Init input from a file
call ddfromfile(finname, ddx_data, tol, charges, error)
call check_error(error)
call allocate_state(ddx_data % params, ddx_data % constants, state, error)
call check_error(error)

! Allocate resources
allocate(phi_cav(ddx_data % constants % ncav), gradphi_cav(3, ddx_data % constants % ncav), &
    & hessianphi_cav(3, 3, ddx_data % constants % ncav), &
    & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), force(3, ddx_data % params % nsph), &
    & stat=istatus)
if(istatus .ne. 0) call test_error(-1, "Allocation failed")
! Prepare host-code-related entities
call mkrhs(ddx_data % params, ddx_data % constants, ddx_data % workspace, 1, &
    & phi_cav, 1, gradphi_cav, 1, hessianphi_cav, psi, charges)
! Use the solver
call ddsolve_legacy(ddx_data, state, phi_cav, -gradphi_cav, hessianphi_cav, psi, tol, esolv, force, error)
call check_error(error)
! compute the second contribution to the forces
call grad_phi_for_charges(ddx_data % params, ddx_data % constants, &
    & ddx_data % workspace, state, charges, force, error)
call check_error(error)
! Open output file for reading
open(unit=100, file=foutname, form='formatted', access='sequential')
! Skip 
read(100, *) ! Using provided ...
read(100, *) ! mkrhs time ...
read(100, *) ! ddsolve time ...
read(100, "(A)") tmpstr
read(tmpstr(16:), *) esolv2
if(abs(esolv2-esolv) .gt. threshold*abs(esolv)) then
    close(100)
    call test_error(-1, "Energy is different")
end if
! Compare forces
if (ddx_data % params % force .eq. 1) then
    read(100, *) ! Full forces ...
    fnorm = dnrm2(3*ddx_data % params % nsph, force, 1)
    do isph = 1, ddx_data % params % nsph
        read(100, *) i, ftmp(:)
        force(:, isph) = force(:, isph) - ftmp
    end do
    fdiff = dnrm2(3*ddx_data % params % nsph, force, 1)
    if(fdiff .gt. threshold*fnorm) then
        close(100)
        call test_error(-1, "Forces are different")
    end if
end if
! Close output file and deallocate resources
close(100)
deallocate(phi_cav, gradphi_cav, psi, force, charges, stat=istatus)
if(istatus .ne. 0) call test_error(-1, "Deallocation failed")
call deallocate_state(state, error)
call deallocate_model(ddx_data, error)

contains
!> Print error message and exit with provided error code
subroutine test_error(code, message)
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    write(0, "(A,A)") "ERROR: ", message
    write(0, "(A,I2)") "CODE:  ", code
    stop -1
end subroutine test_error
end program main


