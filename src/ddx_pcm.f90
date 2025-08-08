!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_pcm.f90
!! PCM solver
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

!> High-level subroutines for ddpcm
module ddx_pcm
! Get ddx-operators
use ddx_operators
implicit none

!> @defgroup Fortran_interface_ddpcm Fortran interface: ddpcm
!! Exposed ddpcm modules in the Fortran API

contains

!> Compute the ddPCM energy
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] constants: Precomputed constants
!! @param[in] state: ddx state (contains solutions and RHSs)
!! @param[out] esolv: resulting energy
!! @param[inout] ddx_error: ddX error
!!
subroutine pcm_energy(constants, state, esolv, ddx_error)
    implicit none
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(in) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(out) :: esolv
    real(dp), external :: ddot
    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue
    esolv = pt5*ddot(constants % n, state % xs, 1, state % psi, 1)
end subroutine pcm_energy

!> Given the potential at the cavity points, assemble the RHS for ddCOSMO
!> or for ddPCM.
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[inout] state: ddx state
!! @param[in] phi_cav: electrostatic potential at the cavity points
!! @param[in] psi: representation of the solute density
!! @param[inout] ddx_error: ddX error
!!
subroutine pcm_setup(params, constants, workspace, state, phi_cav, psi, ddx_error)
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
    state % rhs_done = .true.
    state % adjoint_rhs_done = .true.
end subroutine pcm_setup

!> Do a guess for the primal ddPCM linear system
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!!
subroutine pcm_guess(params, constants, workspace, state, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error

    state % xs = zero
    call prec_repsx(params, constants, workspace, state % phi, &
        & state % phieps, ddx_error)

end subroutine pcm_guess

!> Do a guess for the adjoint ddPCM linear system
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!!
subroutine pcm_guess_adjoint(params, constants, workspace, state, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error

    state % y = zero
    call ldm1x(params, constants, workspace, state % psi, state % s, ddx_error)

end subroutine pcm_guess_adjoint

!> Solve the ddPCM primal linear system
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!! @param[inout] ddx_error     : ddX error
!!
subroutine pcm_solve(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: tol
    ! local variables
    real(dp) :: start_time, finish_time

    call rinfx(params, constants, workspace, state % phi, state % phiinf, &
        & ddx_error)

    state % phieps_niter = params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % phiinf, &
        & state % phieps, state % phieps_niter, state % phieps_rel_diff, &
        & repsx, prec_repsx, hnorm, ddx_error)
    finish_time = omp_get_wtime()
    state % phieps_time = finish_time - start_time

    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddpcm_solve: solver for ddPCM " // &
            & "system did not converge, exiting")
        return
    end if

    state % xs_niter =  params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % phieps, &
        & state % xs, state % xs_niter, state % xs_rel_diff, lx, ldm1x, &
        & hnorm, ddx_error)
    finish_time = omp_get_wtime()
    state % xs_time = finish_time - start_time

    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddpcm_solve: solver for ddCOSMO " // &
            & "system did not converge, exiting")
        return
    end if

end subroutine pcm_solve

!> Solve the ddPCM adjpint linear system
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!! @param[inout] ddx_error: ddX error
!!
subroutine pcm_solve_adjoint(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error
    ! local variables
    real(dp) :: start_time, finish_time

    state % s_niter =  params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % psi, &
        & state % s, state % s_niter, state % s_rel_diff, lstarx, ldm1x, &
        & hnorm, ddx_error)
    finish_time = omp_get_wtime()
    state % s_time = finish_time - start_time

    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddpcm_solve_adjoint: solver for ddCOSMO " // &
            & "system did not converge, exiting")
        return
    end if

    state % y_niter = params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % s, state % y, &
        & state % y_niter, state % y_rel_diff, repsstarx, prec_repsstarx, &
        & hnorm, ddx_error)
    finish_time = omp_get_wtime()
    state % y_time = finish_time - start_time

    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddpcm_solve_adjoint: solver for ddPCM " // &
            & "system did not converge, exiting")
        return
    end if

    ! compute the real adjoint solution and store it in Q
    state % q = state % s - fourpi/(params % eps - one)*state % y

    call pcm_derivative_setup(params, constants, workspace, state)

end subroutine pcm_solve_adjoint

!> Compute the solvation contribution to the ddPCM forces
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions and relevant quantities
!! @param[in] e_cav: electric field, size (3, ncav)
!! @param[out] force       : Geometrical contribution to the forces
!! @param[inout] ddx_error: ddX error
!!
subroutine pcm_solvation_force_terms(params, constants, workspace, &
        & state, e_cav, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: e_cav(3, constants % ncav)
    real(dp), intent(out) :: force(3, params % nsph)

    real(dp) :: start_time, finish_time
    integer :: isph

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    start_time = omp_get_wtime()

    call gradr(params, constants, workspace, state % g, state % ygrid, force)

    do isph = 1, params % nsph
        call contract_grad_L(params, constants, isph, state % xs, &
            & state % sgrid, workspace % tmp_vylm(:, 1), &
            & workspace % tmp_vdylm(:, :, 1), workspace % tmp_vplm(:, 1), &
            & workspace % tmp_vcos(:, 1), workspace % tmp_vsin(:, 1), &
            & force(:, isph))
        call contract_grad_U(params, constants, isph, state % qgrid, &
            & state % phi_grid, force(:, isph))
    end do
    force = - pt5 * force

    call zeta_grad(params, constants, state, e_cav, force)

    finish_time = omp_get_wtime()
    state % force_time = finish_time - start_time

end subroutine pcm_solvation_force_terms

!> This routines precomputes the intermediates to be used in the evaluation
!! of ddCOSMO analytical derivatives.
!!
!! @param[in] params: ddx parameters
!! @param[in] constant: ddx constants
!! @param[inout] workspace: ddx workspaces
!! @param[inout] state: ddx state
!!
subroutine pcm_derivative_setup(params, constants, workspace, state)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state

    real(dp), external :: ddot
    integer :: icav, isph, igrid

    ! dummy operation on unused interface arguments
    if (allocated(workspace % tmp_pot)) continue

    ! Get grid values of S and Y
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & constants % nbasis, one, constants % vgrid, constants % vgrid_nbasis, &
        & state % s, constants % nbasis, zero, state % sgrid, params % ngrid)
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & constants % nbasis, one, constants % vgrid, constants % vgrid_nbasis, &
        & state % y, constants % nbasis, zero, state % ygrid, params % ngrid)

    ! Get the values of phi on the grid
    call ddcav_to_grid_work(params % ngrid, params % nsph, constants % ncav, &
        & constants % icav_ia, constants % icav_ja, state % phi_cav, &
        & state % phi_grid)

    state % g = - (state % phieps - state % phi)
    state % qgrid = state % sgrid - fourpi/(params % eps - one)*state % ygrid

    ! assemble the intermediate zeta: S weighted by U evaluated on the
    ! exposed grid points.
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if (constants % ui(igrid, isph) .ne. zero) then
                icav = icav + 1
                state % zeta(icav) = constants % wgrid(igrid) * &
                    & constants % ui(igrid, isph) * ddot(constants % nbasis, &
                    & constants % vgrid(1, igrid), 1, state % q(1, isph), 1)
            end if
        end do
    end do

end subroutine pcm_derivative_setup

end module ddx_pcm
