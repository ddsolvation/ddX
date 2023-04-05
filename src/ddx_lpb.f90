!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_core.f90
!! Core routines and parameters of entire ddX software
!!
!! @version 1.0.0
!! @author Abhinav Jha and Michele Nottoli
!! @date 2021-02-25

!> High-level subroutines for ddlpb
module ddx_lpb
! Get ddx-operators
use ddx_operators
implicit none

!> @defgroup Fortran_interface_ddlpb Fortran interface: ddlpb

contains

!> ddLPB solver
!!
!! Solves the LPB problem using a domain decomposition approach.
!!
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] phi_cav: Potential at cavity points, size (ncav)
!! @params[in] gradphi_cav: Electric field at cavity points, size (3, ncav)
!! @param[in] psi: Representation of the solute potential in spherical
!!     harmonics, size (nbasis, nsph)
!! @param[in] tol: Tolerance for the linear system solver
!! @param[out] esolv: Solvation energy
!! @param[out] force: Solvation contribution to the forces
!!
subroutine ddlpb(params, constants, workspace, state, phi_cav, gradphi_cav, &
        & hessianphi_cav, psi, tol, esolv, force)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(constants % ncav), &
        & gradphi_cav(3, constants % ncav), &
        & hessianphi_cav(3, 3, constants % ncav), &
        & psi( constants % nbasis,  params % nsph), tol
    real(dp), intent(out) :: esolv, force(3, params % nsph)

    call ddlpb_setup(params, constants, workspace, state, phi_cav, &
        & gradphi_cav, psi)
    if (workspace % error_flag .eq. 1) return
    call ddlpb_guess(params, constants, workspace, state, tol)
    if (workspace % error_flag .eq. 1) return
    call ddlpb_solve(params, constants, workspace, state, tol)
    if (workspace % error_flag .eq. 1) return

    ! Compute the solvation energy
    call ddlpb_energy(constants, state, esolv)

    ! Get forces if needed
    if(params % force .eq. 1) then
        call ddlpb_guess_adjoint(params, constants, workspace, state, tol)
        if (workspace % error_flag .eq. 1) return
        call ddlpb_solve_adjoint(params, constants, workspace, state, tol)
        if (workspace % error_flag .eq. 1) return
        ! TODO: (if easy) remove hessianphi_cav
        call ddlpb_solvation_force_terms(params, constants, workspace, &
            & state, hessianphi_cav, force)
        if (workspace % error_flag .eq. 1) return
    endif

end subroutine ddlpb

!> Given the potential and the electric field at the cavity points,
!> assemble the RHS for ddLPB
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[inout] state: ddx state
!! @param[in] phi_cav: electrostatic potential at the cavity points
!! @param[in] gradphi_cav: electrostatic field at the cavity points
!! @param[in] psi: representation of the solute density
!!
subroutine ddlpb_setup(params, constants, workspace, state, phi_cav, &
        & gradphi_cav, psi)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(constants % ncav)
    real(dp), intent(in) :: gradphi_cav(3, constants % ncav)
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)

    state % psi = psi
    state % rhs_adj_lpb(:, :, 1) = psi
    ! state % rhs_adj_lpb(:, :, 1) = psi/fourpi
    state % rhs_adj_lpb(:, :, 2) = 0.0d0

    !! Setting initial values to zero
    state % g_lpb = zero
    state % f_lpb = zero
    state % phi_grid = zero

    ! Unwrap sparsely stored potential at cavity points phi_cav into phi_grid
    ! and multiply it by characteristic function at cavity points ui
    call ddcav_to_grid_work(params % ngrid, params % nsph, &
        & constants % ncav, constants % icav_ia, &
        & constants % icav_ja, phi_cav, state % phi_grid)
    workspace % tmp_cav = phi_cav * constants % ui_cav
    call ddcav_to_grid_work(params % ngrid, params % nsph, &
        & constants % ncav, constants % icav_ia, &
        & constants % icav_ja, workspace % tmp_cav, &
        & workspace % tmp_grid)
    state % g_lpb = - workspace % tmp_grid

    ! store gradphi_cav for later use in the forces
    state % gradphi_cav = gradphi_cav

    ! wghpot_f : Intermediate computation of F_0 Eq.(75) from QSM19.SISC
    call wghpot_f(params, constants, workspace, gradphi_cav, state % f_lpb)

    ! Setting of the local variables
    state % rhs_lpb = zero

    ! integrate RHS
    call ddintegrate(params % nsph, constants % nbasis, &
        & params % ngrid, constants % vwgrid, &
        & constants % vgrid_nbasis, state % g_lpb, state % rhs_lpb(:,:,1))
    call ddintegrate(params % nsph, constants % nbasis, &
        & params % ngrid, constants % vwgrid, &
        & constants % vgrid_nbasis, state % f_lpb, state % rhs_lpb(:,:,2))
    state % rhs_lpb(:,:,1) = state % rhs_lpb(:,:,1) + state % rhs_lpb(:,:,2)

end subroutine ddlpb_setup

!> Compute the ddLPB energy
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] constants: Precomputed constants
!! @param[in] state: ddx state (contains solutions and RHSs)
!! @param[out] esolv: resulting energy
!!
subroutine ddlpb_energy(constants, state, esolv)
    implicit none
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(in) :: state
    real(dp), intent(out) :: esolv
    real(dp), external :: ddot
    ! TODO: sort once and for all the fourpi issue
    esolv = pt5*ddot(constants % n, state % x_lpb(:,:,1), 1, state % psi, 1)
end subroutine ddlpb_energy

!> Do a guess for the primal ddLPB linear system
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params: User specified parameters
!! @param[inout] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: tolerance
!!
subroutine ddlpb_guess(params, constants, workspace, state, tol)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol

    ! set the inner tolerance
    constants % inner_tol =  sqrt(tol)

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call prec_tx(params, constants, workspace, state % rhs_lpb, state % x_lpb)
    state % x_lpb = state % x_lpb / fourpi

end subroutine ddlpb_guess

!> Do a guess for the adjoint ddLPB linear system
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: tolerance
!!
subroutine ddlpb_guess_adjoint(params, constants, workspace, state, tol)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol

    ! set the inner tolerance
    constants % inner_tol =  sqrt(tol)

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call prec_tstarx(params, constants, workspace, state % rhs_adj_lpb, &
        & state % x_adj_lpb)

    state % x_adj_lpb = state % x_adj_lpb / fourpi
end subroutine ddlpb_guess_adjoint

!> Solve the ddLPB primal linear system
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!!
subroutine ddlpb_solve(params, constants, workspace, state, tol)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    real(dp) :: start_time

    state % x_lpb_niter = params % maxiter
    state % x_lpb = state % x_lpb * fourpi

    ! solve LS using Jacobi/DIIS
    start_time = omp_get_wtime()
    call jacobi_diis_external(params, constants, workspace, &
        & 2*constants % n, tol, state % rhs_lpb, state % x_lpb, &
        & state % x_lpb_niter, state % x_lpb_rel_diff, cx, prec_tx, rmsnorm)
    if (workspace % error_flag .ne. 0) then
        workspace % error_message = 'Jacobi solver failed to converge " // &
            & "in ddlpb_solve'
        return
    end if
    state % x_lpb_time = omp_get_wtime() - start_time
    ! the integral operators of ddLPB are defined according to a
    ! different convention, so we need to scale the density by 1/fourpi
    state % x_lpb = state % x_lpb / fourpi
end subroutine ddlpb_solve


!> Solve the adjoint ddLPB linear system
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!!
subroutine ddlpb_solve_adjoint(params, constants, workspace, state, tol)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol

    real(dp) :: start_time

    state % x_adj_lpb_niter = params % maxiter
    constants % inner_tol = sqrt(tol)

    state % x_adj_lpb = state % x_adj_lpb * fourpi

    ! solve adjoint LS using Jacobi/DIIS
    start_time = omp_get_wtime()
    call jacobi_diis_external(params, constants, workspace, &
        & 2*constants % n, tol, state % rhs_adj_lpb, state % x_adj_lpb, &
        & state % x_adj_lpb_niter, state % x_adj_lpb_rel_diff, &
        & cstarx, prec_tstarx, rmsnorm)
    if (workspace % error_flag .ne. 0) then
        workspace % error_message = 'Jacobi solver failed to ' // &
            & 'converge in ddlpb_solve_adjoint'
        return
    end if
    state % x_adj_lpb_time = omp_get_wtime() - start_time
    ! the integral operators of ddLPB are defined according to a
    ! different convention, so we need to scale the adjoint density
    ! by 1/fourpi
    state % x_adj_lpb = state % x_adj_lpb / fourpi
end subroutine ddlpb_solve_adjoint

!!
!! Wrapper routine for the computation of ddPCM forces. It makes the
!! interface easier to implement. If a fine control is needed, the
!! worker routine should be directly called.
!!
!! @param[in] params         : General options
!! @param[in] constants      : Precomputed constants
!! @param[inout] workspace   : Preallocated workspaces
!! @param[inout] state       : Solutions and relevant quantities
!! @param[in] hessianphi_cav : Electric field gradient at the grid points
!! @param[out] force         : Geometrical contribution to the forces
!!
subroutine ddlpb_solvation_force_terms(params, constants, workspace, &
        & state, hessianphi_cav, force)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: hessianphi_cav(3, 3, constants % ncav)
    real(dp), intent(out) :: force(3, params % nsph)

    ! local
    real(dp), dimension(constants % nbasis) :: basloc, vplm
    real(dp), dimension(3, constants % nbasis) :: dbasloc
    real(dp), dimension(params % lmax + 1) :: vsin, vcos

    ! large local are allocatable
    real(dp), allocatable :: ef(:,:), xadj_r_sgrid(:,:), xadj_e_sgrid(:,:), &
        & normal_hessian_cav(:,:), diff_re(:,:), scaled_xr(:,:)
    integer :: isph, icav, icav_gr, icav_ge, igrid, istat
    integer :: i
    real(dp), external :: ddot, dnrm2
    real(dp) :: start_time, finish_time

    start_time = omp_get_wtime()
    allocate(ef(3, params % nsph), &
        & xadj_r_sgrid(params % ngrid, params % nsph), &
        & xadj_e_sgrid(params % ngrid, params % nsph), &
        & normal_hessian_cav(3, constants % ncav), &
        & diff_re(constants % nbasis, params % nsph), &
        & scaled_xr(constants % nbasis, params % nsph), stat=istat)
    if (istat.ne.0) then
        workspace % error_message = 'Allocation failed in ddlpb_force_worker'
        workspace % error_flag = 1
        return
    end if

    state % x_lpb = state % x_lpb * fourpi

    diff_re = zero
    vsin = zero
    vcos = zero
    vplm = zero
    basloc = zero
    dbasloc = zero
    ef = zero
    force = zero

    ! Compute the derivative of the normal derivative of psi_0
    normal_hessian_cav = zero
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .gt. zero) then
                icav = icav + 1
                do i = 1, 3
                    normal_hessian_cav(:, icav) = normal_hessian_cav(:,icav) +&
                        & hessianphi_cav(:,i,icav)*constants % cgrid(i,igrid)
                end do
            end if
        end do
    end do

    ! Call dgemm to integrate the adjoint solution on the grid points
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & one, constants % vgrid, constants % vgrid_nbasis, &
        & state % x_adj_lpb(:, :, 1), constants % nbasis, zero, &
        & Xadj_r_sgrid, params % ngrid)
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & one, constants % vgrid, constants % vgrid_nbasis, &
        & state % x_adj_lpb(:,:,2), constants % nbasis, zero, &
        & Xadj_e_sgrid, params % ngrid)

    ! Scale by the factor of 1/(4Pi/(2l+1))
    scaled_Xr = state % x_lpb(:,:,1)
    call convert_ddcosmo(params, constants, -1, scaled_Xr)

    !$omp parallel do default(none) shared(params,constants,workspace, &
    !$omp scaled_xr,xadj_r_sgrid,state,force,xadj_e_sgrid) &
    !$omp private(isph,basloc,dbasloc,vplm,vcos,vsin) &
    !$omp schedule(static,1)
    do isph = 1, params % nsph
        ! Compute A^k*Xadj_r, using Subroutine from ddCOSMO
        call contract_grad_L(params, constants, isph, scaled_Xr, Xadj_r_sgrid, &
            & basloc, dbasloc, vplm, vcos, vsin, force(:,isph))
        ! Compute B^k*Xadj_e
        call contract_grad_B(params, constants, isph, &
            & state % x_lpb(:,:,2), Xadj_e_sgrid, force(:, isph))
        ! Computation of G0
        call contract_grad_U(params, constants, isph, Xadj_r_sgrid, &
            & state % phi_grid, force(:, isph))
    end do
    ! Compute C1 and C2 contributions
    diff_re = zero
    call contract_grad_C(params, constants, workspace, state % x_lpb(:,:,1), &
        & state % x_lpb(:,:,2), Xadj_r_sgrid, Xadj_e_sgrid, &
        & state % x_adj_lpb(:,:,1), state % x_adj_lpb(:,:,2), force, &
        & diff_re)
    ! Computation of G0 continued

    ! NOTE: contract_grad_U returns a positive summation
    if (workspace % error_flag .eq. 1) return
    force = -force
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .ne. zero) then
                icav = icav + 1
                state % zeta(icav) = constants % wgrid(igrid) * &
                    & constants % ui(igrid, isph) * ddot(constants % nbasis, &
                    & constants % vgrid(1, igrid), 1, &
                    & state % x_adj_lpb(1, isph, 1), 1)
            end if
        end do
    end do

    icav_gr = zero
    icav_ge = zero
    ! Computation of F0
    call contract_grad_f(params, constants, workspace, &
        & state % x_adj_lpb(:,:,1) + state % x_adj_lpb(:,:,2), &
        & Xadj_r_sgrid + xadj_e_sgrid, state % gradphi_cav, &
        & normal_hessian_cav, icav_gr, force, state)
    if (workspace % error_flag .eq. 1) return

    force = pt5*force

    deallocate(ef, xadj_r_sgrid, xadj_e_sgrid, normal_hessian_cav, &
        & diff_re, scaled_xr, stat=istat)
    if (istat.ne.0) then
        workspace % error_message = 'Deallocation failed in ddlpb_force_worker'
        workspace % error_flag = 1
        return
    end if
    finish_time = omp_get_wtime()
    state % force_time = finish_time - start_time

    ! restore x_lpb
    state % x_lpb = state % x_lpb / fourpi

end subroutine ddlpb_solvation_force_terms

end module ddx_lpb
