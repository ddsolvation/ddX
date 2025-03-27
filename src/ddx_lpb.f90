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

!> Given the potential and the electric field at the cavity points,
!> assemble the RHS for ddLPB
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[inout] state: ddx state
!! @param[in] phi_cav: electrostatic potential at the cavity points
!! @param[in] e_cav: electrostatic field at the cavity points
!! @param[in] psi: representation of the solute density
!! @param[inout] ddx_error: ddX error
!!
subroutine lpb_setup(params, constants, workspace, state, phi_cav, &
        & e_cav, psi, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: phi_cav(constants % ncav)
    real(dp), intent(in) :: e_cav(3, constants % ncav)
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)

    state % psi = psi
    state % rhs_adj_lpb(:, :, 1) = psi
    ! in ddLPB the polarization density X is expanded according to a
    ! slightly different definition which lacks the 4pi/(2l + 1) factors
    ! which are present in ddCOSMO and ddPCM. This affect the adjoint
    ! density as well. To have a consistent adjoint density with ddCOSMO
    ! and ddPCM, we have to scale the RHS.
    call convert_ddcosmo(params, constants, -1, state % rhs_adj_lpb(:,:,1))
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

    ! store the gradient of the potential (- electric field)
    state % gradphi_cav = - e_cav

    ! wghpot_f : Intermediate computation of F_0 Eq.(75) from QSM19.SISC
    call wghpot_f(params, constants, workspace, state % gradphi_cav, &
        & state % f_lpb, ddx_error)

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

end subroutine lpb_setup

!> Compute the ddLPB energy
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] constants: Precomputed constants
!! @param[in] state: ddx state (contains solutions and RHSs)
!! @param[out] esolv: resulting energy
!! @param[inout] ddx_error: ddX error
!!
subroutine lpb_energy(constants, state, esolv, ddx_error)
    implicit none
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(in) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(out) :: esolv
    real(dp), external :: ddot

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    esolv = pt5*ddot(constants % n, state % x_lpb(:,:,1), 1, state % psi, 1)
end subroutine lpb_energy

!> Do a guess for the primal ddLPB linear system
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params: User specified parameters
!! @param[inout] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: tolerance
!! @param[inout] ddx_error: ddX error
!!
subroutine lpb_guess(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error

    ! set the inner tolerance
    constants % inner_tol =  sqrt(tol)

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call prec_tx(params, constants, workspace, state % rhs_lpb, &
        & state % x_lpb, ddx_error)

end subroutine lpb_guess

!> Do a guess for the adjoint ddLPB linear system
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: tolerance
!! @param[inout] ddx_error: ddX error
!!
subroutine lpb_guess_adjoint(params, constants, workspace, state, tol, &
        & ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error

    ! set the inner tolerance
    constants % inner_tol =  sqrt(tol)

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call prec_tstarx(params, constants, workspace, state % rhs_adj_lpb, &
        & state % x_adj_lpb, ddx_error)

end subroutine lpb_guess_adjoint

!> Solve the ddLPB primal linear system
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!!
!! @param[inout] ddx_error: ddX error
subroutine lpb_solve(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp) :: start_time

    state % x_lpb_niter = params % maxiter
    workspace % xs_time = zero
    workspace % hsp_time = zero

    ! solve LS using Jacobi/DIIS
    start_time = omp_get_wtime()
    call jacobi_diis_external(params, constants, workspace, &
        & 2*constants % n, tol, state % rhs_lpb, state % x_lpb, &
        & state % x_lpb_niter, state % x_lpb_rel_diff, cx, prec_tx, rmsnorm, &
        & ddx_error)
    if (ddx_error % flag.ne. 0) then
        call update_error(ddx_error, 'Jacobi solver failed to converge " // &
            & "in ddlpb_solve')
        return
    end if
    state % x_lpb_time = omp_get_wtime() - start_time
    ! in ddLPB the polarization density X is expanded according to a
    ! slightly different definition which lacks the 4pi/(2l + 1) factors
    ! which are present in ddCOSMO and ddPCM. Before exiting, we scale
    ! the density so that it is consistent with the ddCOSMO and ddPCM
    ! ones.
    call convert_ddcosmo(params, constants, -1, state % x_lpb(:,:,1))

    ! put the timings in the right places
    state % xs_time = workspace % xs_time
    state % hsp_time = workspace % hsp_time
end subroutine lpb_solve

!> Solve the adjoint ddLPB linear system
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!! @param[inout] ddx_error: ddX error
!!
subroutine lpb_solve_adjoint(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error

    real(dp) :: start_time

    state % x_adj_lpb_niter = params % maxiter
    workspace % s_time = zero
    workspace % hsp_adj_time = zero

    ! solve adjoint LS using Jacobi/DIIS
    start_time = omp_get_wtime()
    call jacobi_diis_external(params, constants, workspace, &
        & 2*constants % n, tol, state % rhs_adj_lpb, state % x_adj_lpb, &
        & state % x_adj_lpb_niter, state % x_adj_lpb_rel_diff, &
        & cstarx, prec_tstarx, rmsnorm, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, 'Jacobi solver failed to ' // &
            & 'converge in ddlpb_solve_adjoint')
        return
    end if
    state % x_adj_lpb_time = omp_get_wtime() - start_time

    state % q = state % x_adj_lpb(:, :, 1)

    ! put the timings in the right places
    state % s_time = workspace % s_time
    state % hsp_adj_time = workspace % hsp_adj_time

    call lpb_derivative_setup(params, constants, workspace, state, ddx_error)

end subroutine lpb_solve_adjoint

!> Compute the solvation terms of the forces (solute aspecific). This
!! must be summed to the solute specific term to get the full forces.
!!
!> @ingroup Fortran_interface_ddlpb
!! @param[in] params         : General options
!! @param[in] constants      : Precomputed constants
!! @param[inout] workspace   : Preallocated workspaces
!! @param[inout] state       : Solutions and relevant quantities
!! @param[in] hessianphi_cav : Electric field gradient at the grid points
!! @param[out] force         : Geometrical contribution to the forces
!! @param[inout] ddx_error: ddX error
!!
subroutine lpb_solvation_force_terms(params, constants, workspace, &
        & state, hessianphi_cav, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: hessianphi_cav(3, 3, constants % ncav)
    real(dp), intent(out) :: force(3, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error

    ! local
    real(dp), dimension(constants % nbasis) :: basloc, vplm
    real(dp), dimension(3, constants % nbasis) :: dbasloc
    real(dp), dimension(params % lmax + 1) :: vsin, vcos

    ! large local are allocatable
    real(dp), allocatable :: normal_hessian_cav(:,:), &
        & diff_re(:,:), scaled_xr(:,:)
    integer :: isph, icav, igrid, istat
    integer :: i
    real(dp), external :: ddot, dnrm2
    real(dp) :: start_time, finish_time

    start_time = omp_get_wtime()
    allocate(normal_hessian_cav(3, constants % ncav), &
        & diff_re(constants % nbasis, params % nsph), &
        & scaled_xr(constants % nbasis, params % nsph), stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, 'Allocation failed in ddlpb_force_worker')
        return
    end if

    ! we have to insert again the factor that we removed in ddlpb_solve
    call convert_ddcosmo(params, constants, 1, state % x_lpb(:, :, 1))

    diff_re = zero
    vsin = zero
    vcos = zero
    vplm = zero
    basloc = zero
    dbasloc = zero
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

    ! Scale by the factor of 1/(4Pi/(2l+1))
    scaled_Xr = state % x_lpb(:,:,1)
    call convert_ddcosmo(params, constants, -1, scaled_Xr)

    !$omp parallel do default(none) shared(params,constants,workspace, &
    !$omp scaled_xr,state,force) private(isph,basloc,dbasloc,vplm,vcos,vsin) &
    !$omp schedule(static,1)
    do isph = 1, params % nsph
        ! Compute A^k*Xadj_r, using Subroutine from ddCOSMO
        call contract_grad_L(params, constants, isph, scaled_Xr, &
            & state % x_adj_r_grid, basloc, dbasloc, vplm, vcos, vsin, &
            & force(:,isph))
        ! Compute B^k*Xadj_e
        call contract_grad_B(params, constants, isph, state % x_lpb(:,:,2), &
            & state % x_adj_e_grid, force(:, isph))
        ! Computation of G0
        call contract_grad_U(params, constants, isph, state % x_adj_r_grid, &
            & state % phi_grid, force(:, isph))
    end do
    ! Compute C1 and C2 contributions
    diff_re = zero
    call contract_grad_C(params, constants, workspace, state % x_lpb(:,:,1), &
        & state % x_lpb(:,:,2), state % x_adj_r_grid, state % x_adj_e_grid, &
        & state % x_adj_lpb(:,:,1), state % x_adj_lpb(:,:,2), force, &
        & diff_re, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "contract_grad_C returned an error, exiting")
        return
    end if
    ! Computation of G0 continued

    ! NOTE: contract_grad_U returns a positive summation
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

    ! Computation of F0
    call contract_grad_f(params, constants, workspace, &
        & state % x_adj_lpb(:,:,1) + state % x_adj_lpb(:,:,2), &
        & state % x_adj_re_grid, state % gradphi_cav, &
        & normal_hessian_cav, force, state, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "contract_grad_f returned an error, exiting")
        return
    end if

    force = pt5*force

    deallocate(normal_hessian_cav, diff_re, scaled_xr, stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, 'Deallocation failed in ddlpb_force_worker')
        return
    end if

    ! and now we remove again the factor, so that the density is
    ! restored.
    call convert_ddcosmo(params, constants, -1, state % x_lpb(:, :, 1))

    ! Finally add the zeta - electric field contribution. Since the
    ! subroutine takes as input the electric field, we temporarily swap
    ! the sign of the gradient of the potential (this trick avoids
    ! dangerous memory allocation)
    state % gradphi_cav = - state % gradphi_cav
    call zeta_grad(params, constants, state, state % gradphi_cav, force)
    state % gradphi_cav = - state % gradphi_cav

    finish_time = omp_get_wtime()
    state % force_time = finish_time - start_time

end subroutine lpb_solvation_force_terms

!> This routines precomputes two intermediates for its later usage in
!! the computation of analytical derivatives (forces or other).
!!
!! @param[in] params: ddx parameters
!! @param[in] constant: ddx constants
!! @param[inout] workspace: ddx workspaces
!! @param[inout] state: ddx state
!! @param[inout] ddx_error: ddX error
!!
subroutine lpb_derivative_setup(params, constants, workspace, state, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error

    real(dp), allocatable :: coefy_d(:, :, :)
    integer :: istat, jsph, igrid0, icav, isph, igrid, l0, m0, ind0, &
        & inode, indl
    real(dp) :: rijn, sum_int, fac
    real(dp) :: work(constants % lmax0 + 1), coef(constants % nbasis0), &
        & vij(3), vtij(3), sij(3)
    complex(dp) :: work_complex(constants % lmax0 + 1)

    ! expand the adjoint solutions at the Lebedev points and compute
    ! their sum
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & constants % nbasis, one, constants % vgrid, &
        & constants % vgrid_nbasis, state % x_adj_lpb(:, :, 1), &
        & constants % nbasis, zero, state % x_adj_r_grid, params % ngrid)
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & constants % nbasis, one, constants % vgrid, &
        & constants % vgrid_nbasis, state % x_adj_lpb(:, :, 2), &
        & constants % nbasis, zero, state % x_adj_e_grid, params % ngrid)
    state % x_adj_re_grid = state % x_adj_r_grid + state % x_adj_e_grid

    if (params % fmm .eq. 0) then
        allocate(coefy_d(constants % ncav, params % ngrid, params % nsph), &
            & stat=istat)
        if (istat.ne.0) then
            call update_error(ddx_error, "allocation ddx_error in " // &
                & "ddlpb_derivative_setup")
            return
        end if
        coefy_d = zero
        do jsph = 1, params % nsph
            do igrid0 = 1, params % ngrid
                icav = zero
                do isph = 1, params % nsph
                    do igrid = 1, params % ngrid
                        if(constants % ui(igrid, isph) .gt. zero) then
                            icav = icav + 1
                            vij = params % csph(:,isph) + params % rsph(isph) &
                                & *constants % cgrid(:,igrid) &
                                & - params % csph(:,jsph)
                            rijn = sqrt(dot_product(vij,vij))
                            if (rijn.ne.zero) then
                                sij = vij/rijn
                            else
                                sij = one
                            end if
                            do l0 = 0, constants % lmax0
                                do m0 = -l0, l0
                                    ind0 = l0**2 + l0 + m0 + 1
                                    coef(ind0) = &
                                        & constants % vgrid(ind0, igrid0) * &
                                        & constants % C_ik(l0, jsph)
                                end do
                            end do
                            vtij = vij*params % kappa
                            call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                                & constants % vscales, &
                                & constants % sk_ri(:, jsph), one, &
                                & coef, zero, coefy_d(icav, igrid0, jsph), &
                                & work_complex, work)
                        end if
                    end do
                end do
            end do
        end do
        do jsph = 1, params % nsph
            do igrid0 = 1, params % ngrid
                icav = zero
                sum_int = zero
                do isph = 1, params % nsph
                    do igrid = 1, params % ngrid
                        if(constants % ui(igrid, isph) .gt. zero) then
                            icav = icav + 1
                            sum_int = sum_int + coefy_d(icav, igrid0, jsph) &
                                & *state % x_adj_re_grid(igrid, isph) &
                                & *constants % wgrid(igrid) &
                                & *constants % ui(igrid, isph)
                        end if
                    end do
                end do
                state % phi_n(igrid0, jsph) = &
                    & - (params % epsp/params % eps)*sum_int
            end do
        end do
    else
        ! Adjoint integration from spherical harmonics to grid points is not needed
        ! here as ygrid already contains grid values, we just need to scale it by
        ! weights of grid points
        do isph = 1, params % nsph
            workspace % tmp_grid(:, isph) = state % x_adj_re_grid(:, isph) * &
                & constants % wgrid(:) * constants % ui(:, isph)
        end do
        ! Adjoint FMM
        call tree_m2p_bessel_adj(params, constants, constants % lmax0, one, &
            & workspace % tmp_grid, zero, params % lmax, workspace % tmp_sph)
        call tree_l2p_bessel_adj(params, constants, one, &
            & workspace % tmp_grid, zero, workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l, workspace % tmp_node_m)
        call tree_m2m_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_m)
        ! Properly load adjoint multipole harmonics into tmp_sph
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph(1:constants % nbasis0, isph) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph) + &
                    & workspace % tmp_node_m(1:constants % nbasis0, inode)
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph(1:indl, isph) = &
                    & workspace % tmp_sph(1:indl, isph) + &
                    & workspace % tmp_node_m(:, inode)
            end do
        end if
        ! Scale by C_ik
        do isph = 1, params % nsph
            do l0 = 0, constants % lmax0
                ind0 = l0*l0 + l0 + 1
                workspace % tmp_sph(ind0-l0:ind0+l0, isph) = &
                    & workspace % tmp_sph(ind0-l0:ind0+l0, isph) * &
                    & constants % C_ik(l0, isph)
            end do
        end do
        ! Multiply by vgrid
        call dgemm('T', 'N', params % ngrid, params % nsph, &
            & constants % nbasis0, -params % epsp/params % eps, &
            & constants % vgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, &
            & state % phi_n, params % ngrid)
    end if

    ! assemble the pseudo-dipoles
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if (constants % ui(igrid, isph) .gt. zero) then
                icav = icav + 1
                fac = - constants % ui(igrid, isph)*constants % wgrid(igrid) &
                    & *state % phi_n(igrid, isph)
                state % zeta_dip(1, icav) = fac*constants % cgrid(1, igrid)
                state % zeta_dip(2, icav) = fac*constants % cgrid(2, igrid)
                state % zeta_dip(3, icav) = fac*constants % cgrid(3, igrid)
            end if
        end do
    end do

    if (allocated(coefy_d)) then
        deallocate(coefy_d, stat=istat)
        if (istat.ne.0) then
            call update_error(ddx_error, &
                & "deallocation ddx_error in ddlpb_derivative_setup")
            return
        end if
    end if

end subroutine lpb_derivative_setup

end module ddx_lpb
