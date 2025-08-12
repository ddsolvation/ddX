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

!> Standalone application of ddX: this program can compute the solvation
!! energy and forces for a solute made of point charges. All the
!! relevant steps are here outlined.
program main
use ddx
use ddx_multipolar_solutes
use omp_lib
implicit none

character(len=255) :: fname
character(len=2047) :: banner
type(ddx_type) :: ddx_data
type(ddx_state_type) :: state
type(ddx_error_type) :: ddx_error
type(ddx_electrostatics_type) :: electrostatics
real(dp), allocatable :: psi(:, :), force(:, :), charges(:), &
    & multipoles(:, :)
real(dp) :: tol, esolv, start_time, finish_time
integer :: i, isph, info
100 format(1X,A,ES11.4E2,A)
200 format(1X,A,I4,A,ES21.14)
300 format(1X,A,I4)
400 format(1X,A,ES25.16E3)
500 format(1X,I5,3ES25.16E3)

! Get the input file from the command line arguments
call get_command_argument(1, fname)
write(6, *) "Using provided file ", trim(fname), " as a config file"

! STEP 1: Initialization of the model.
! Read the input file and call "allocate_model" to initialize the model.
! The model is a container for all the parameters, precomputed constants
! and preallocated workspaces.
start_time = omp_get_wtime()
call ddfromfile(fname, ddx_data, tol, charges, ddx_error)
finish_time = omp_get_wtime()
write(*, 100) "Initialization time:", finish_time - start_time, " seconds"
call check_error(ddx_error)

! STEP 2: Initialization of the state.
! The state is a high level object which is related to solving the
! solvation problem for a given solute. A state must be used with a
! model (ddx_data). Different states can be used at the same time with
! a given model, for instance when solving for different solutes,
! or for different states of the solute.
call allocate_state(ddx_data % params, ddx_data % constants, state, ddx_error)
call check_error(ddx_error)

! Print the ddX banner
call get_banner(banner)
write(6, *) trim(banner)

! STEP 3: Building the required RHSs and electrostatic properties.
!
! dd models require two different kind of RHS. The primal RHS which is
! based on electrostatic properties, and the adjoint RHS which is a
! representation of the density in spherical harmonics (with extra steps).
! Furthermore, when computing the forces, electrostatic properties of
! order +1 are required to compute the force contribution from the RHS.
!
! ddLPB is slightly more complicated as its primal RHS requires also
! the electric field, and if forces are to be computed, requires the
! electric field gradient.
!
! For convenience all the required electrostatic properties are computed
! together in STEP 3a. The adjoint RHS is computed in STEP 3b.

! STEP 3a: electrostatic properties for primal RHS and for the forces.

start_time = omp_get_wtime()

! The ddx_multipolar_solute provides routines to compute the
! electrostatic properties up to field gradient for multipolar
! distributions of arbitrary order, provided that they are given in
! real spherical harmonics. Here we have charges, so we need to convert
! them to monopoles. For charges the conversion is simply a scaling by
! 1/sqrt(four*pi). For more complicated cases, up to the 16-poles, check
! https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
allocate(multipoles(1, ddx_data % params % nsph), stat=info)
if (info .ne. 0) then
    write(6, *) "Allocation failed in ddx_driver"
    stop 1
end if
multipoles(1, :) = charges/sqrt4pi

! compute the required electrostatics properties for a multipolar solute
call multipole_electrostatics(ddx_data % params, ddx_data % constants, &
    & ddx_data % workspace, multipoles, 0, electrostatics, ddx_error)

finish_time = omp_get_wtime()
write(*, 100) "Electrostatic properties time:", finish_time-start_time, &
    & " seconds"

! STEP 3b: adjoint RHS.
start_time = omp_get_wtime()
allocate(psi(ddx_data % constants % nbasis, ddx_data % params % nsph), &
    & stat=info)
if (info .ne. 0) then
    write(6, *) "Allocation failed in ddx_driver"
    stop 1
end if
call multipole_psi(ddx_data % params, multipoles, 0, psi)
finish_time = omp_get_wtime()
write(*, 100) "Psi time:", finish_time-start_time, " seconds"

call time_push()
if (ddx_data % params % force .eq. 1) then
    allocate(force(3, ddx_data % params % nsph), stat=info)
    if (info .ne. 0) then
        write(6, *) "Allocation failed in ddx_driver"
        stop 1
    end if
    call ddrun(ddx_data, state, electrostatics, psi, tol, esolv, ddx_error, &
        & force)
else
    call ddrun(ddx_data, state, electrostatics, psi, tol, esolv, ddx_error)
end if
call check_error(ddx_error)
call time_pull("ddrun")

if (ddx_data % params % force .eq. 1) write(*, 100) &
    & "solvation force terms time:", state % force_time, " seconds"

! STEP 8: if required compute the solute specific contributions to the
! forces.
if (ddx_data % params % force .eq. 1) then
    start_time = omp_get_wtime()
    call time_push()
    call multipole_force_terms(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, 0, multipoles, force, ddx_error)
    call time_pull("solute forces")
    call check_error(ddx_error)
    finish_time = omp_get_wtime()
    write(*, 100) "multipolar force terms time:", &
        & finish_time - start_time, " seconds"
end if

! From now on, we print the results.

! Print info on the primal ddPCM system
if (ddx_data % params % model .eq. 2) then
    ! Print each iteration if needed
    do i = 1, state % phieps_niter
        write(6, 200) "iter=", i, " relative difference: ", &
            & state % phieps_rel_diff(i)
    end do
    ! Print number of iterations and time
    write(6, 100) "ddpcm step time: ", state % phieps_time, &
        & " seconds"
    write(6, 300) "ddpcm step iterations: ", &
        & state % phieps_niter
end if

! Print info on the primal ddLPB system
if (ddx_data % params % model .eq. 3) then
    ! Print each iteration if needed
    do i = 1, state % x_lpb_niter
        write(6, 200) "iter=", i, " relative difference: ", &
            & state % x_lpb_rel_diff(i)
    end do
    ! Print number of iterations and time
    write(6, 100) "ddlpb step time: ", state % x_lpb_time, &
        & " seconds, of which:"
    write(6, 100) "    ddcosmo: ", state % xs_time, " seconds"
    write(6, 100) "    hsp: ", state % hsp_time, " seconds"
    write(6, 100) "    coupling terms: ", state % x_lpb_time &
        & - state % xs_time - state % hsp_time, " seconds"
    write(6, 300) "ddlpb step iterations: ", state % x_lpb_niter
end if

! Print info on the primal ddCOSMO system
! Print each iteration if needed
if ((ddx_data % params % model.eq.1) .or. &
        & (ddx_data % params % model.eq.2)) then
    do i = 1, state % xs_niter
        write(6, 200) "iter=", i, " relative difference: ", &
            & state % xs_rel_diff(i)
    end do
    ! Print number of iterations and time
    write(6, 100) "ddcosmo step time: ", state % xs_time, " seconds"
    write(6, 300) "ddcosmo step iterations: ", state % xs_niter
end if

! Print info on the adjoint solver
if (ddx_data % params % force .eq. 1) then
    ! Print info on the adjoint ddCOSMO system
    if ((ddx_data % params % model.eq.1) .or. &
            & (ddx_data % params % model.eq.2)) then
        do i = 1, state % s_niter
            write(6, 200) "iter=", i, " relative difference: ", &
                & state % s_rel_diff(i)
        end do
        ! Print number of iterations and time
        write(6, 100) "adjoint ddcosmo step time: ", state % s_time, " seconds"
        write(6, 300) "adjoint ddcosmo step iterations: ", state % s_niter
    end if

    ! Print info on the adjoint ddPCM system
    if (ddx_data % params % model .eq. 2) then
        ! Print each iteration if needed
        do i = 1, state % y_niter
            write(6, 200) "iter=", i, " relative difference: ", &
                & state % y_rel_diff(i)
        end do
        write(6, 100) "adjoint ddpcm step time: ", state % y_time, " seconds"
        write(6, 300) "adjoint ddpcm step iterations: ", state % y_niter
    end if

    ! Print info on the adjoint ddLPB system
    if (ddx_data % params % model .eq. 3) then
        do i = 1, state % x_adj_lpb_niter
            write(6, 200) "iter=", i, " relative difference: ", &
                & state % x_adj_lpb_rel_diff(i)
        end do
        ! Print number of iterations and time
        write(6, 100) "adjoint ddlpb step time: ", &
            & state % x_adj_lpb_time, " seconds, of which:"
        write(6, 100) "    adjoint ddcosmo: ", state % s_time, " seconds"
        write(6, 100) "    adjoint hsp: ", state % hsp_adj_time, " seconds"
        write(6, 100) "    adjoint coupling terms: ", state % x_adj_lpb_time &
            & - state % s_time - state % hsp_adj_time, " seconds"
        write(6, 200) "adjoint ddlpb step iterations: ", &
            & state % x_adj_lpb_niter
    end if
end if

write(6, 400) "Solvation energy (Hartree):", esolv
write(6, 400) "Solvation energy (kcal/mol):", esolv*tokcal
if (ddx_data % params % force .eq. 1) then
    write(6, *) "Full forces (kcal/mol/A)"
    do isph = 1, ddx_data % params % nsph
        write(6, 500) isph, force(:, isph)*tokcal/toang
    end do
end if

! Clean the workspace.

deallocate(psi, multipoles, charges, stat=info)
if (info .ne. 0) then
    write(6, *) "Deallocation failed in ddx_driver"
    stop 1
end if
if (allocated(force)) then
    deallocate(force, stat=info)
    if (info .ne. 0) then
        write(6, *) "Deallocation failed in ddx_driver"
        stop 1
    end if
end if

call deallocate_electrostatics(electrostatics, ddx_error)
call deallocate_state(state, ddx_error)
call deallocate_model(ddx_data, ddx_error)

end program main
