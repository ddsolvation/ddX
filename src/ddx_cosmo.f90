!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_cosmo.f90
!! COSMO solver
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

!> High-level subroutines for ddcosmo
module ddx_cosmo
! Get ddx-operators
use ddx_operators
use ddx_multipolar_solutes
implicit none

!> @defgroup Fortran_interface_ddcosmo Fortran interface: ddcosmo
!! Exposed ddcosmo modules in the Fortran API

contains

!> Compute the ddCOSMO energy
!!
!> @ingroup Fortran_interface_ddcosmo
!! @param[in] constants: Precomputed constants
!! @param[in] state: ddx state (contains solutions and RHSs)
!! @param[out] esolv: resulting energy
!! @param[inout] ddx_error: ddX error
!!
subroutine cosmo_energy(constants, state, esolv, ddx_error)
    implicit none
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(in) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(out) :: esolv
    real(dp), external :: ddot

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    esolv = pt5*ddot(constants % n, state % xs, 1, state % psi, 1)
end subroutine cosmo_energy

!> Given the potential at the cavity points, assemble the RHS for ddCOSMO
!> or for ddPCM.
!!
!> @ingroup Fortran_interface_ddcosmo
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[inout] state: ddx state
!! @param[in] phi_cav: electrostatic potential at the cavity points
!! @param[in] psi: representation of the solute density
!! @param[inout] ddx_error: ddX error
!!
subroutine cosmo_setup(params, constants, workspace, state, phi_cav, psi, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: phi_cav(constants % ncav)
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    call cav_to_spherical(params, constants, workspace, phi_cav, &
        & state % phi)
    state % phi = - state % phi
    state % phi_cav = phi_cav
    state % psi = psi
end subroutine cosmo_setup

!> Do a guess for the primal ddCOSMO linear system
!!
!> @ingroup Fortran_interface_ddcosmo
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[inout] ddx_error: ddX error
!!
subroutine cosmo_guess(params, constants, workspace, state, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error

    ! apply the diagonal preconditioner as a guess
    call ldm1x(params, constants, workspace, state % phi, state % xs, ddx_error)

end subroutine cosmo_guess

!> Do a guess for the adjoint ddCOSMO linear system
!!
!> @ingroup Fortran_interface_ddcosmo
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[inout] ddx_error: ddX error
!!
subroutine cosmo_guess_adjoint(params, constants, workspace, state, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error

    ! apply the diagonal preconditioner as a guess
    call ldm1x(params, constants, workspace, state % psi, state % s, ddx_error)

end subroutine cosmo_guess_adjoint

!> Solve the primal ddCOSMO linear system
!!
!> @ingroup Fortran_interface_ddcosmo
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: Tolerance for the linear system solver
!! @param[inout] ddx_error: ddX error
!!
subroutine cosmo_solve(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error
    ! local variables
    real(dp) :: start_time, finish_time

    state % xs_niter =  params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % phi, &
        & state % xs, state % xs_niter, state % xs_rel_diff, lx, ldm1x, &
        & hnorm, ddx_error)
    finish_time = omp_get_wtime()
    state % xs_time = finish_time - start_time

end subroutine cosmo_solve

!> Solve the adjoint ddCOSMO linear system
!!
!> @ingroup Fortran_interface_ddcosmo
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: Tolerance for the linear system solver
!! @param[inout] ddx_error: ddX error
!!
subroutine cosmo_solve_adjoint(params, constants, workspace, state, tol, &
        & ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error
    ! local variables
    real(dp) :: start_time, finish_time

    state % s_niter = params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % psi, &
        & state % s, state % s_niter, state % s_rel_diff, lstarx, ldm1x, &
        & hnorm, ddx_error)
    finish_time = omp_get_wtime()
    state % s_time = finish_time - start_time

    state % q = state % s

    call cosmo_derivative_setup(params, constants, workspace, state)

end subroutine cosmo_solve_adjoint

!> Compute the solvation term of the forces (solute aspecific). This must
!> be summed to the solute specific term to get the full forces.
!!
!> @ingroup Fortran_interface_ddcosmo
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] e_cav: electric field, size (3, ncav)
!! @param[inout] force: force term
!! @param[inout] ddx_error: ddX error
!!
subroutine cosmo_solvation_force_terms(params, constants, workspace, &
        & state, e_cav, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: e_cav(3, constants % ncav)
    real(dp), intent(inout) :: force(3, params % nsph)
    ! local variables
    real(dp) :: start_time, finish_time
    integer :: isph

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    start_time = omp_get_wtime()

    force = zero
    do isph = 1, params % nsph
        call contract_grad_l(params, constants, isph, state % xs, &
            & state % sgrid, workspace % tmp_vylm(:, 1), &
            & workspace % tmp_vdylm(:, :, 1), workspace % tmp_vplm(:, 1), &
            & workspace % tmp_vcos(:, 1), workspace % tmp_vsin(:, 1), &
            & force(:, isph))
        call contract_grad_u(params, constants, isph, state % sgrid, &
            & state % phi_grid, force(:, isph))
    end do

    force = -pt5*force

    call zeta_grad(params, constants, state, e_cav, force)

    finish_time = omp_get_wtime()
    state % force_time = finish_time - start_time

end subroutine cosmo_solvation_force_terms

!> This routines precomputes the intermediates to be used in the evaluation
!! of ddCOSMO analytical derivatives.
!!
!! @param[in] params: ddx parameters
!! @param[in] constant: ddx constants
!! @param[inout] workspace: ddx workspaces
!! @param[inout] state: ddx state
!!
subroutine cosmo_derivative_setup(params, constants, workspace, state)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state

    real(dp), external :: ddot
    integer :: icav, isph, igrid

    ! dummy operation on unused interface arguments
    if (allocated(workspace % tmp_pot)) continue

    ! Get values of S on the grid
    call ddeval_grid_work(constants % nbasis, params % ngrid, params % nsph, &
        & constants % vgrid, constants % vgrid_nbasis, one, state % s, zero, &
        & state % sgrid)
    ! Get the values of phi on the grid
    call ddcav_to_grid_work(params % ngrid, params % nsph, constants % ncav, &
        & constants % icav_ia, constants % icav_ja, state % phi_cav, &
        & state % phi_grid)

    ! assemble the intermediate zeta: S weighted by U evaluated on the
    ! exposed grid points.
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if (constants % ui(igrid, isph) .ne. zero) then
                icav = icav + 1
                state % zeta(icav) = constants % wgrid(igrid) * &
                    & constants % ui(igrid, isph) * ddot(constants % nbasis, &
                    & constants % vgrid(1, igrid), 1, state % s(1, isph), 1)
            end if
        end do
    end do

end subroutine cosmo_derivative_setup

end module ddx_cosmo

