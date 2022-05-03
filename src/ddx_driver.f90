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

!> Standalone application of ddX
program main
! Get the high-level ddX-module 
use ddx
! Enable OpenMP
use omp_lib
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
type(ddx_state_type) :: state
integer :: iprint, info
integer :: phi_flag=1, grad_flag=1, hessian_flag=1
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), &
    & hessianphi_cav(:, :, :), psi(:, :), force(:, :)
real(dp) :: tol, esolv, start_time, finish_time
integer :: i, j, isph

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, ddx_data, tol, iprint, info)
if(info .ne. 0) stop "info != 0"

call ddx_init_state(ddx_data % params, ddx_data % constants, state)

call print_header(iprint,ddx_data%params)
! determine needed arrays
if (ddx_data % params % model .eq. 3) then
    phi_flag = 1
    grad_flag = 1
    hessian_flag = 0
    if (ddx_data % params % force .eq. 1) hessian_flag = 1
else
    phi_flag = 1
    grad_flag = 0
    hessian_flag = 0
    if (ddx_data % params % force .eq. 1) grad_flag = 1
end if

! allocate for everything
! it is possible to pass to mkrhs unallocated arrays as optional arguments
! but only since fortran 2008
allocate(phi_cav(ddx_data % constants % ncav), &
    & gradphi_cav(3, ddx_data % constants % ncav), &
    & hessianphi_cav(3, 3, ddx_data % constants % ncav), &
    & force(3, ddx_data % params % nsph), &
    & psi(ddx_data % constants % nbasis, ddx_data % params % nsph))

start_time = omp_get_wtime()
call mkrhs(ddx_data % params, ddx_data % constants, ddx_data % workspace, &
    & phi_flag, phi_cav, grad_flag, gradphi_cav, hessian_flag, hessianphi_cav, &
    & psi)
finish_time = omp_get_wtime()
write(*, "(A,ES11.4E2,A)") " mkrhs time:", finish_time-start_time, " seconds"

start_time = omp_get_wtime()
call ddsolve(ddx_data, state,  phi_cav, gradphi_cav, hessianphi_cav, psi, &
    & tol, esolv, force, info)
finish_time = omp_get_wtime()

! Print info depending on iprint flag
if (iprint .gt. 0) then
    ! Print info on the primal ddPCM system
    if (ddx_data % params % model .eq. 2) then
        ! Print each iteration if needed
        if (iprint .gt. 1 .and. ddx_data % params % itersolver.eq.1) then
            do i = 1, state % phieps_niter
                print " (A,I4,A,ES20.14)", "iter=", i, &
                    & " relative difference: ", state % phieps_rel_diff(i)
            end do
        end if
        ! Print number of iterations and time
        print " (A,ES11.4E2,A)", " ddpcm step time: ", state % phieps_time, &
            & " seconds"
        print "(A,I4)", " ddpcm step iterations: ", &
            & state % phieps_niter
    end if
    ! Print info on the primal ddLPB system
    if (ddx_data % params % model .eq. 3) then
        ! Print each iteration if needed
        if (iprint .gt. 1 .and. ddx_data % params % itersolver.eq.1) then
            do i = 1, state % phieps_niter
                print " (A,I4,A,ES20.14)", "iter=", i, &
                    & " relative difference: ", state % x_lpb_rel_diff(i)
            end do
        end if
        ! Print number of iterations and time
        print " (A,ES11.4E2,A)", " ddlpb step time: ", state % x_lpb_time, &
            & " seconds"
        print "(A,I4)", " ddlpb step iterations: ", &
            & state % x_lpb_niter
    end if
    ! Print info on the primal ddCOSMO system
    ! Print each iteration if needed
    if (iprint .gt. 1 .and. ddx_data % params % itersolver.eq.1) then
        do i = 1, state % xs_niter
            print "(A,I4,A,ES20.14)", " iter=", i, &
                & " relative difference: ", state % xs_rel_diff(i)
        end do
    end if
    ! Print number of iterations and time
    print "(A,ES11.4E2,A)", " ddcosmo step time: ", state % xs_time, &
        & " seconds"
    print "(A,I4)", " ddcosmo step iterations: ", state % xs_niter
    ! Print info on the adjoint solver
    if (ddx_data % params % force .eq. 1 .and. ddx_data % params % itersolver.eq.1) then
        ! Print info on the adjoint ddCOSMO system
        ! Print each iteration if needed
        if (iprint .gt. 1) then
            do i = 1, state % s_niter
                print "(A,I4,A,ES20.14)", " iter=", i, &
                    & " relative difference: ", state % s_rel_diff(i)
            end do
        end if
        ! Print number of iterations and time
        print "(A,ES11.4E2,A)", " adjoint ddcosmo step time: ", &
            & state % s_time, " seconds"
        print "(A,I4)", " adjoint ddcosmo step iterations: ", &
            & state % s_niter
        ! Print info on the adjoint ddPCM system
        if (ddx_data % params % model .eq. 2) then
            ! Print each iteration if needed
            if (iprint .gt. 1 .and. ddx_data % params % itersolver.eq.1) then
                do i = 1, state % y_niter
                    print "(A,I4,A,ES20.14)", " iter=", i, &
                        & " relative difference: ", state % y_rel_diff(i)
                end do
            end if
            ! Print number of iterations and time
            print "(A,ES11.4E2,A)", " adjoint ddpcm step time: ", &
                & state % y_time, " seconds"
            print "(A,I4)", " adjoint ddpcm step iterations: ", &
                & state % y_niter
        end if
        ! Print info on the adjoint ddLPB system
        if (ddx_data % params % model .eq. 3) then
            ! Print each iteration if needed
            if (iprint .gt. 1 .and. ddx_data % params % itersolver.eq.1) then
                do i = 1, state % y_niter
                    print "(A,I4,A,ES20.14)", " iter=", i, &
                        & " relative difference: ", state % x_adj_lpb_rel_diff(i)
                end do
            end if
            ! Print number of iterations and time
            print "(A,ES11.4E2,A)", " adjoint ddlpb step time: ", &
                & state % x_adj_lpb_time, " seconds"
            print "(A,I4)", " adjoint ddlpb step iterations: ", &
                & state % x_adj_lpb_niter
        end if
    end if
end if
write(*, "(A,ES11.4E2,A)") " ddx_driver time:", finish_time-start_time, " seconds"
write(*, "(A,ES25.16E3)") " Solvation energy:", esolv
write(*, "(A,ES25.16E3)") " Solvation energy (kJ/mol):", esolv*tokj
if (ddx_data % params % force .eq. 1) then
    write(*, *) " Full forces (kJ/mol/A)"
    do isph = 1, ddx_data % params % nsph
        write(6,'(1x,i5,3ES25.16E3)') isph, force(:,isph)*tokj/toang
    end do
end if

! deallocation
deallocate(psi, phi_cav, gradphi_cav, hessianphi_cav, force)

call ddx_free_state(state)
call ddfree(ddx_data)

end program main

