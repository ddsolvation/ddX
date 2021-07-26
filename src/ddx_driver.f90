!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_driver.f90
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

character(len=255) :: fname
type(ddx_type) :: ddx_data
integer :: iprint, info
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), psi(:, :), force(:, :)
real(dp) :: tol, esolv, start_time, finish_time
integer :: i, j, isph

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, ddx_data, tol, iprint, info)
if(info .ne. 0) stop "info != 0"
allocate(phi_cav(ddx_data % constants % ncav), gradphi_cav(3, ddx_data % constants % ncav), &
    & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), force(3, ddx_data % params % nsph))
call cpu_time(start_time)
call mkrhs(ddx_data, phi_cav, gradphi_cav, psi)
call cpu_time(finish_time)
write(*, "(A,ES11.4E2,A)") "mkrhs time:", finish_time-start_time, " seconds"
call cpu_time(start_time)
call ddsolve(ddx_data, phi_cav, gradphi_cav, psi, tol, esolv, force, info)
call cpu_time(finish_time)
write(*, "(A,ES11.4E2,A)") "ddsolve time:", finish_time-start_time, " seconds"
write(*, "(A,ES25.16E3)") "ddsolve esolv:", esolv
if (ddx_data % params % force .eq. 1) then
write(*, *) "Full forces"
    do isph = 1, ddx_data % params % nsph
        write(6,'(1x,i5,3ES25.16E3)') isph, force(:,isph)
    end do
end if
deallocate(phi_cav, gradphi_cav, psi, force)
call ddfree(ddx_data)

end program main

