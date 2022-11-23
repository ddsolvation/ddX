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

!> ddPCM solver
!!
!! Solves the problem within PCM model using a domain decomposition approach.
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] phi_cav: Potential at cavity points, size (ncav)
!! @param[in] psi: Representation of the solute potential in spherical
!!     harmonics, size (nbasis, nsph)
!! @param[in] tol: Tolerance for the linear system solver
!! @param[out] esolv: Solvation energy
!! @param[out] force: Solvation contribution to the forces
!!
subroutine ddpcm(params, constants, workspace, state, phi_cav, &
        & psi, tol, esolv, force)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(constants % ncav), &
        & psi(constants % nbasis, params % nsph), tol
    real(dp), intent(out) :: esolv, force(3, params % nsph)
    real(dp), external :: ddot

    call ddpcm_setup(params, constants, workspace, state, phi_cav, psi)
    call ddpcm_guess(params, constants, workspace, state)
    call ddpcm_solve(params, constants, workspace, state, tol)

    ! Compute the solvation energy
    esolv = pt5*ddot(constants % n, state % xs, 1, psi, 1)

    ! Get forces if needed
    if (params % force .eq. 1) then
        ! solve the adjoint
        call ddpcm_guess_adjoint(params, constants, workspace, state)
        call ddpcm_solve_adjoint(params, constants, workspace, state, tol)

        ! evaluate the solvent unspecific contribution analytical derivatives
        force = zero
        call ddpcm_solvation_force_terms(params, constants, workspace, &
            & state, force)
    end if

end subroutine ddpcm

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
!!
subroutine ddpcm_setup(params, constants, workspace, state, phi_cav, psi)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(constants % ncav)
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)
    call cav_to_spherical(params, constants, workspace, phi_cav, &
        & state % phi)
    state % phi = - state % phi
    state % phi_cav = phi_cav
    state % psi = psi
end subroutine ddpcm_setup

!> Do a guess for the primal ddPCM linear system
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!!
subroutine ddpcm_guess(params, constants, workspace, state)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state

    state % xs = zero
    call prec_repsx(params, constants, workspace, state % phi, state % phieps)

end subroutine ddpcm_guess

!> Do a guess for the adjoint ddPCM linear system
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!!
subroutine ddpcm_guess_adjoint(params, constants, workspace, state)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state

    state % y = zero
    call ldm1x(params, constants, workspace, state % psi, state % s)

end subroutine ddpcm_guess_adjoint

!> Solve the ddPCM primal linear system
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!!
subroutine ddpcm_solve(params, constants, workspace, state, tol)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    ! local variables
    real(dp) :: start_time, finish_time

    call rinfx(params, constants, workspace, state % phi, state % phiinf)

    state % phieps_niter = params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % phiinf, &
        & state % phieps, state % phieps_niter, state % phieps_rel_diff, &
        & repsx, prec_repsx, hnorm)
    finish_time = omp_get_wtime()
    state % phieps_time = finish_time - start_time

    if (workspace % error_flag .ne. 0) then
        workspace % error_message = "ddpcm_energy: solver for ddPCM " // &
            & "system did not converge"
        return
    end if

    state % xs_niter =  params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % phieps, &
        & state % xs, state % xs_niter, state % xs_rel_diff, lx, ldm1x, hnorm)
    finish_time = omp_get_wtime()
    state % xs_time = finish_time - start_time

    if (workspace % error_flag .ne. 0) then
        workspace % error_message = "ddpcm_energy: solver for ddCOSMO " // &
            & "system did not converge"
        return
    end if

end subroutine ddpcm_solve

!> Solve the ddPCM adjpint linear system
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!!
subroutine ddpcm_solve_adjoint(params, constants, workspace, state, tol)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    ! local variables
    real(dp) :: start_time, finish_time

    state % s_niter =  params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % psi, &
        & state % s, state % s_niter, state % s_rel_diff, lstarx, ldm1x, &
        & hnorm)
    finish_time = omp_get_wtime()
    state % s_time = finish_time - start_time

    if (workspace % error_flag .ne. 0) then
        workspace % error_message = "ddpcm_energy: solver for adjoint " // &
            & "ddCOSMO system did not converge"
        return
    end if

    state % y_niter = params % maxiter
    start_time = omp_get_wtime()
    call jacobi_diis(params, constants, workspace, tol, state % s, state % y, &
        & state % y_niter, state % y_rel_diff, repsstarx, prec_repsstarx, &
        & hnorm)
    finish_time = omp_get_wtime()
    state % y_time = finish_time - start_time

    if (workspace % error_flag .ne. 0) then
        workspace % error_message = "ddpcm_energy: solver for adjoint " // &
            & "ddPCM system did not converge"
        return
    end if

end subroutine ddpcm_solve_adjoint

!> Compute the solvation contribution to the ddPCM forces
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions and relevant quantities
!! @param[out] force       : Geometrical contribution to the forces
!!
subroutine ddpcm_solvation_force_terms(params, constants, workspace, &
        & state, force)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(out) :: force(3, params % nsph)
    ! local
    real(dp), external :: ddot
    integer :: icav, isph, igrid

    call ddcav_to_grid_work(params % ngrid, params % nsph, constants % ncav, &
        & constants % icav_ia, constants % icav_ja, state % phi_cav, &
        & state % phi_grid)
    ! Get grid values of S and Y
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & constants % nbasis, one, constants % vgrid, constants % vgrid_nbasis, &
        & state % s, constants % nbasis, zero, state % sgrid, params % ngrid)
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & constants % nbasis, one, constants % vgrid, constants % vgrid_nbasis, &
        & state % y, constants % nbasis, zero, state % ygrid, params % ngrid)

    state % g = - (state % phieps - state % phi)
    state % q = state % s - fourpi/(params % eps - one)*state % y
    state % qgrid = state % sgrid - fourpi/(params % eps - one)*state % ygrid

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

    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .ne. zero) then
                icav = icav + 1
                state % zeta(icav) = constants % wgrid(igrid) * &
                    & constants % ui(igrid, isph) * ddot(constants % nbasis, &
                    & constants % vgrid(1, igrid), 1, &
                    & state % q(1, isph), 1)
            end if
        end do
    end do

end subroutine ddpcm_solvation_force_terms

end module ddx_pcm
