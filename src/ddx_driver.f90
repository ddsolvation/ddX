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
use ddx_lpb
use omp_lib
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
integer :: iprint, info
integer :: phi_flag=1, grad_flag=0, hessian_flag=0
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), &
    & hessianphi_cav(:, :, :), psi(:, :), force(:, :)
real(dp) :: tol, esolv, start_time, finish_time
integer :: i, j, isph

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, ddx_data, tol, iprint, info)
if(info .ne. 0) stop "info != 0"

allocate(phi_cav(ddx_data % constants % ncav), &
    & psi(ddx_data % constants % nbasis, ddx_data % params % nsph))

if(ddx_data % params % force .eq. 1) then
    allocate(gradphi_cav(3, ddx_data % constants % ncav), &
        & hessianphi_cav(3, 3, ddx_data % constants % ncav), &
        & force(3, ddx_data % params % nsph))
    phi_flag = 1
    hessian_flag = 1
end if

start_time = omp_get_wtime()
call mkrhs(ddx_data, phi_flag, phi_cav, grad_flag, gradphi_cav, hessian_flag, &
    & hessianphi_cav, psi)
finish_time = omp_get_wtime()
write(*, "(A,ES11.4E2,A)") "mkrhs time:", finish_time-start_time, " seconds"
start_time = omp_get_wtime()
call ddsolve(ddx_data, phi_cav, gradphi_cav, hessianphi_cav, psi, tol, esolv, &
    & force, info)
finish_time = omp_get_wtime()
! Print info depending on iprint flag
if (iprint .gt. 0) then
    ! Print info on the primal ddPCM system
    if (ddx_data % params % model .eq. 2) then
        ! Print each iteration if needed
        if (iprint .gt. 1) then
            do i = 1, ddx_data % phieps_niter
                print "(A,I4,A,ES20.14)", "iter=", i, &
                    & " relative difference: ", ddx_data % phieps_rel_diff(i)
            end do
        end if
        ! Print number of iterations and time
        print "(A,ES11.4E2,A)", "ddpcm step time: ", ddx_data % phieps_time, &
            & " seconds"
        print "(A,I4)", "ddpcm step iterations: ", &
            & ddx_data % phieps_niter
    end if
    ! Print info on the primal ddCOSMO system
    ! Print each iteration if needed
    if (iprint .gt. 1) then
        do i = 1, ddx_data % xs_niter
            print "(A,I4,A,ES20.14)", "iter=", i, &
                & " relative difference: ", ddx_data % xs_rel_diff(i)
        end do
    end if
    ! Print number of iterations and time
    print "(A,ES11.4E2,A)", "ddcosmo step time: ", ddx_data % xs_time, &
        & " seconds"
    print "(A,I4)", "ddcosmo step iterations: ", ddx_data % xs_niter
    ! Print info on the adjoint solver
    if (ddx_data % params % force .eq. 1) then
        ! Print info on the adjoint ddCOSMO system
        ! Print each iteration if needed
        if (iprint .gt. 1) then
            do i = 1, ddx_data % s_niter
                print "(A,I4,A,ES20.14)", "iter=", i, &
                    & " relative difference: ", ddx_data % s_rel_diff(i)
            end do
        end if
        ! Print number of iterations and time
        print "(A,ES11.4E2,A)", "adjoint ddcosmo step time: ", &
            & ddx_data % s_time, " seconds"
        print "(A,I4)", "adjoint ddcosmo step iterations: ", &
            & ddx_data % s_niter
        ! Print info on the adjoint ddPCM system
        if (ddx_data % params % model .eq. 2) then
            ! Print each iteration if needed
            if (iprint .gt. 1) then
                do i = 1, ddx_data % y_niter
                    print "(A,I4,A,ES20.14)", "iter=", i, &
                        & " relative difference: ", ddx_data % y_rel_diff(i)
                end do
            end if
            ! Print number of iterations and time
            print "(A,ES11.4E2,A)", "adjoint ddpcm step time: ", &
                & ddx_data % y_time, " seconds"
            print "(A,I4)", "adjoint ddpcm step iterations: ", &
                & ddx_data % y_niter
        end if
    end if
end if
write(*, "(A,ES11.4E2,A)") "ddx_driver time:", finish_time-start_time, " seconds"
write(*, "(A,ES25.16E3)") "Solvation energy:", esolv
write(*, "(A,ES25.16E3)") "Solvation energy (kJ/mol):", esolv*2625.5002d0
if (ddx_data % params % force .eq. 1) then
write(*, *) "Full forces"
    do isph = 1, ddx_data % params % nsph
        write(6,'(1x,i5,3ES25.16E3)') isph, force(:,isph)
    end do
end if
deallocate(phi_cav, psi)
if (ddx_data % params % force .eq. 1) then
    deallocate(gradphi_cav, hessianphi_cav, force)
end if
call ddfree(ddx_data)

end program main

