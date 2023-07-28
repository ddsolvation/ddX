!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx.f90
!! Main driver routine of the ddX with all the per-model solvers
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

!> High-level module of the ddX software
module ddx
! Get ddcosmo-module
use ddx_cosmo
! Get ddpcm-module
use ddx_pcm
! Get ddlpb-module
use ddx_lpb
implicit none

contains

!> Main solver routine
!!
!! Solves the solvation problem, computes the energy, and if required
!! computes the forces.
!!
!! @param[in] ddx_data: ddX object with all input information
!! @param[inout] state: ddx state (contains RHSs and solutions)
!! @param[in] electrostatics: electrostatic property container
!! @param[in] psi: RHS of the adjoint problem
!! @param[in] tol: tolerance for the linear system solvers
!! @param[out] esolv: solvation energy
!! @param[inout] error: ddX error
!! @param[out] force: Analytical forces (optional argument, only if
!!             required)
!! @param[in] read_guess: optional argument, if true read the guess
!!            from the state object
!!
subroutine ddsolve(ddx_data, state, electrostatics, psi, tol, esolv, &
        & error, force, read_guess)
    type(ddx_type), intent(inout) :: ddx_data
    type(ddx_state_type), intent(inout) :: state
    type(ddx_electrostatics_type), intent(in) :: electrostatics
    real(dp), intent(in) :: psi(ddx_data % constants % nbasis, &
        & ddx_data % params % nsph)
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: esolv
    type(ddx_error_type), intent(inout) :: error
    real(dp), intent(out), optional :: force(3, ddx_data % params % nsph)
    logical, intent(in), optional :: read_guess
    ! local variables
    logical :: do_guess

    ! decide if the guess has to be read or must be done
    if (present(read_guess)) then
        do_guess = .not.read_guess
    else
        do_guess = .true.
    end if

    ! if the forces are to be computed, but the array is not passed, raise
    ! an error
    if ((.not.present(force)) .and. (ddx_data % params % force .eq. 1)) then
        call update_error(error, &
            & "ddsolve: forces are to be computed, but the optional force" // &
            & " array has not been passed.")
        return
    end if

    call setup(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, electrostatics, psi, error)
    if (error % flag .ne. 0) then
        call update_error(error, "ddsolve: setup returned an error, exiting")
        return
    end if

    ! solve the primal linear system
    if (do_guess) then
        call fill_guess(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, state, tol, error)
        if (error % flag .ne. 0) then
            call update_error(error, &
                & "ddsolve: fill_guess returned an error, exiting")
            return
        end if
    end if
    call solve(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, tol, error)
    if (error % flag .ne. 0) then
        call update_error(error, "ddsolve: solve returned an error, exiting")
        return
    end if

    ! compute the energy
    call energy(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, esolv, error)
    if (error % flag .ne. 0) then
        call update_error(error, "ddsolve: energy returned an error, exiting")
        return
    end if

    ! solve the primal linear system
    if (ddx_data % params % force .eq. 1) then
        if (do_guess) then
            call fill_guess_adjoint(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, tol, error)
            if (error % flag .ne. 0) then
                call update_error(error, &
                    & "ddsolve: fill_guess_adjoint returned an error, exiting")
                return
            end if
        end if
        call solve_adjoint(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, state, tol, error)
        if (error % flag .ne. 0) then
            call update_error(error, &
                & "ddsolve: solve_adjoint returned an error, exiting")
            return
        end if
    end if

    ! compute the forces
    if (ddx_data % params % force .eq. 1) then
        force = zero
        call solvation_force_terms(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, state, electrostatics, force, error)
        if (error % flag .ne. 0) then
            call update_error(error, &
                & "ddsolve: solvation_force_terms returned an error, exiting")
            return
        end if
    end if

end subroutine ddsolve

!> Main solver routine
!!
!! Solves the problem within COSMO model using a domain decomposition approach.
!!
!! @param[in] ddx_data: ddX object with all input information
!! @param[in] phi_cav: Potential at cavity points
!! @param[in] e_cav: Gradient of the eletric potential at cavity points
!! @param[in] hessianphi_cav: Hessian of the eletric potential at cavity points
!! @param[inout] state: ddx state (contains RHSs and solutions)
!! @param[in] psi: RHS of the adjoint problem
!! @param[in] tol
!! @param[out] esolv: Solvation energy
!! @param[out] force: Analytical forces
!! @param[inout] error: ddX error
!!
subroutine ddsolve_legacy(ddx_data, state, phi_cav, e_cav, hessianphi_cav, &
        & psi, tol, esolv, force, error)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(ddx_data % constants % ncav), &
        & e_cav(3, ddx_data % constants % ncav), &
        & hessianphi_cav(3, ddx_data % constants % ncav), &
        & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), tol
    ! Outputs
    real(dp), intent(out) :: esolv, force(3, ddx_data % params % nsph)
    type(ddx_error_type), intent(inout) :: error
    ! Find proper model
    select case(ddx_data % params % model)
        ! COSMO model
        case (1)
            call ddcosmo(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, phi_cav, psi, e_cav, &
                & tol, esolv, force, error)
        ! PCM model
        case (2)
            call ddpcm(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, phi_cav, psi, e_cav, &
                & tol, esolv, force, error)
        ! LPB model
        case (3)
            call ddlpb(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, phi_cav, e_cav, &
                & psi, tol, esolv, hessianphi_cav, force, error)
        ! Error case
        case default
            call update_error(error, "unsupported solvation " // &
                & " model in the dd solver.")
            return
    end select
end subroutine ddsolve_legacy

!> Setup the state for the different models
!!
!> @ingroup Fortran_interface
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] electrostatics: electrostatic property container
!! @param[in] psi: Representation of the solute potential in spherical
!!     harmonics, size (nbasis, nsph)
!! @param[inout] error: ddX error
!!
subroutine setup(params, constants, workspace, state, electrostatics, &
        & psi, error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_electrostatics_type), intent(in) :: electrostatics
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: error

    if (params % model .eq. 1) then
        call cosmo_setup(params, constants, workspace, state, &
            & electrostatics % phi_cav, psi, error)
    else if (params % model .eq. 2) then
        call pcm_setup(params, constants, workspace, state, &
            & electrostatics % phi_cav, psi, error)
    else if (params % model .eq. 3) then
        call lpb_setup(params, constants, workspace, state, &
            & electrostatics % phi_cav, electrostatics % e_cav, &
            & psi, error)
    else
        call update_error(error, "Unknow model in setup.")
        return
    end if

end subroutine setup

!> Do a guess for the primal linear system for the different models
!!
!> @ingroup Fortran_interface
!! @param[in] params: User specified parameters
!! @param[inout] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: tolerance
!! @param[inout] error: ddX error
!!
subroutine fill_guess(params, constants, workspace, state, tol, error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: error

    if (params % model .eq. 1) then
        call cosmo_guess(params, constants, workspace, state, error)
    else if (params % model .eq. 2) then
        call pcm_guess(params, constants, workspace, state, error)
    else if (params % model .eq. 3) then
        call lpb_guess(params, constants, workspace, state, tol, error)
    else
        call update_error(error, "Unknow model in fill_guess.")
        return
    end if

end subroutine fill_guess

!> Do a guess for the adjoint linear system for the different models
!!
!> @ingroup Fortran_interface
!! @param[in] params: User specified parameters
!! @param[inout] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: tolerance
!! @param[inout] error: ddX error
!!
subroutine fill_guess_adjoint(params, constants, workspace, state, tol, error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: error

    if (params % model .eq. 1) then
        call cosmo_guess_adjoint(params, constants, workspace, state, error)
    else if (params % model .eq. 2) then
        call pcm_guess_adjoint(params, constants, workspace, state, error)
    else if (params % model .eq. 3) then
        call lpb_guess_adjoint(params, constants, workspace, state, tol, error)
    else
        call update_error(error, "Unknow model in fill_guess_adjoint.")
        return
    end if

end subroutine fill_guess_adjoint

!> Solve the primal linear system for the different models
!!
!> @ingroup Fortran_interface
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!! @param[inout] error: ddX error
!!
subroutine solve(params, constants, workspace, state, tol, error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: error

    if (params % model .eq. 1) then
        call cosmo_solve(params, constants, workspace, state, tol, error)
    else if (params % model .eq. 2) then
        call pcm_solve(params, constants, workspace, state, tol, error)
    else if (params % model .eq. 3) then
        call lpb_solve(params, constants, workspace, state, tol, error)
    else
        call update_error(error, "Unknow model in solve.")
        return
    end if

end subroutine solve

!> Solve the adjoint linear system for the different models
!!
!> @ingroup Fortran_interface
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!! @param[inout] error: ddX error
!!
subroutine solve_adjoint(params, constants, workspace, state, tol, error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: error

    if (params % model .eq. 1) then
        call cosmo_solve_adjoint(params, constants, workspace, state, tol, error)
    else if (params % model .eq. 2) then
        call pcm_solve_adjoint(params, constants, workspace, state, tol, error)
    else if (params % model .eq. 3) then
        call lpb_solve_adjoint(params, constants, workspace, state, tol, error)
    else
        call update_error(error, "Unknow model in solve_adjoint.")
        return
    end if

end subroutine solve_adjoint

!> Compute the energy for the different models
!!
!> @ingroup Fortran_interface
!! @param[in] params: General options
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[in] state: ddx state (contains solutions and RHSs)
!! @param[out] solvation_energy: resulting energy
!! @param[inout] error: ddX error
!!
subroutine energy(params, constants, workspace, state, solvation_energy, error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(in) :: workspace
    type(ddx_state_type), intent(in) :: state
    type(ddx_error_type), intent(inout) :: error
    real(dp), intent(out) :: solvation_energy

    if (params % model .eq. 1) then
        call cosmo_energy(constants, state, solvation_energy, error)
    else if (params % model .eq. 2) then
        call pcm_energy(constants, state, solvation_energy, error)
    else if (params % model .eq. 3) then
        call lpb_energy(constants, state, solvation_energy, error)
    else
        call update_error(error, "Unknow model in energy.")
        return
    end if

end subroutine energy

!> Compute the solvation terms of the forces (solute aspecific) for the
!! different models. This must be summed to the solute specific term to get
!! the full forces
!!
!> @ingroup Fortran_interface
!! @param[in] params: General options
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: Solutions and relevant quantities
!! @param[in] electrostatics: Electrostatic properties container.
!! @param[out] force: Geometrical contribution to the forces
!! @param[inout] error: ddX error
!!
subroutine solvation_force_terms(params, constants, workspace, &
        & state, electrostatics, force, error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_electrostatics_type), intent(in) :: electrostatics
    real(dp), intent(out) :: force(3, params % nsph)
    type(ddx_error_type), intent(inout) :: error

    if (params % model .eq. 1) then
        call cosmo_solvation_force_terms(params, constants, workspace, &
            & state, electrostatics % e_cav, force, error)
    else if (params % model .eq. 2) then
        call pcm_solvation_force_terms(params, constants, workspace, &
            & state, electrostatics % e_cav, force, error)
    else if (params % model .eq. 3) then
        call lpb_solvation_force_terms(params, constants, workspace, &
            & state, electrostatics % g_cav, force, error)
    else
        call update_error(error, "Unknow model in solvation_force_terms.")
        return
    end if

end subroutine solvation_force_terms

end module ddx
