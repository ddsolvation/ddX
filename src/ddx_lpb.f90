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
!!
!! Logical variables for iterations
!!
!! Local variables and their definitions that will be used throughout this file
!! isph    : Index for the sphere i
!! jsph    : Index for the sphere j
!! ksph    : Index for the sphere k
!! igrid   : Index for the grid points
!! icav    : Index for the external grid points
!! ibasis  : Index for the basis
!! ibasis0 : Index for the fixed basis (nbasis0)
!! l       : Index for lmax
!! m       : Index for m:-l,...,l
!! ind     : l^2 + l + m + 1
!! l0      : Index for lmax0
!! m0      : Index for m0:-l0,...,l0
!! ind0    : l0^2 + l0 + m0 + 1
!! ineigh  : Index over Row space of neighbors
!! rijn    : r_j*r_1^j(x_i^n) = |x_i^n-x_j|
!! tij     : r_1^j(x_i^n)
!! xij     : chi_j(x_i^n)
!! oij     : omega^\eta_ijn
!! vij     : x_i^n-x_j
!! sij     : e^j(x_i^n)
!! SI_rijn : Besssel function of first kind for rijn
!! DI_rijn : Derivative of Bessel function of first kind for rijn
!! tlow    : Lower bound for switch region
!! thigh   : Upper bound for switch region
!! basloc  : Y_lm(s_n)
!! dbasloc : Derivative of Y_lm(s_n)
!! i       : Index for dimension x,y,z

contains

!> ddLPB solver
!!
!! Wrapper routine for the solution of the direct ddLPB linear
!! system. It makes the interface easier to implement. If a fine
!! control is needed, the worker routine should be directly called.
!!
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] phi_cav      : Electric potential at the grid points
!! @param[in] gradphi_cav  : Electric field at the grid points
!! @param[in] tol          : Tolerance for the iterative solvers
!!
subroutine ddlpb_solve(params, constants, workspace, state, phi_cav, &
        & gradphi_cav, tol)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(constants % ncav)
    real(dp), intent(in) :: gradphi_cav(3, constants % ncav)
    real(dp), intent(in) :: tol

    state % x_lpb_niter = params % maxiter
    call ddlpb_solve_worker(params, constants, workspace, &
        & phi_cav, gradphi_cav, state % g_lpb, state % f_lpb, &
        & state % phi_grid, state % x_lpb, state % x_lpb_niter, &
        & state % x_lpb_time, state % x_lpb_rel_diff, tol)
end subroutine ddlpb_solve

!!
!! Wrapper routine for the solution of the adjoint ddPCM linear
!! system. It makes the interface easier to implement. If a fine
!! control is needed, the worker routine should be directly called.
!!
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] psi          : Representation of the solute's density
!! @param[in] tol          : Tolerance for the iterative solvers
!!
subroutine ddlpb_adjoint(params, constants, workspace, state, psi, tol)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)
    real(dp), intent(in) :: tol

    state % x_adj_lpb_niter = params % maxiter
    call ddlpb_adjoint_worker(params, constants, &
      & workspace, psi, tol, state % x_adj_lpb, state % x_adj_lpb_niter, &
      & state % x_adj_lpb_time, state % x_adj_lpb_rel_diff)
end subroutine ddlpb_adjoint

!!
!! Wrapper routine for the computation of ddPCM forces. It makes the
!! interface easier to implement. If a fine control is needed, the
!! worker routine should be directly called.
!!
!! @param[in] params         : General options
!! @param[in] constants      : Precomputed constants
!! @param[inout] workspace   : Preallocated workspaces
!! @param[inout] state       : Solutions and relevant quantities
!! @param[in] phi_cav        : Electric potential at the grid points
!! @param[in] gradphi_cav    : Electric field at the grid points
!! @param[in] hessianphi_cav : Electric field gradient at the grid points
!! @param[in] psi            : Representation of the solute's density
!! @param[out] force         : Geometrical contribution to the forces
!!
subroutine ddlpb_force(params, constants, workspace, state, phi_cav, &
        & gradphi_cav, hessianphi_cav, psi, force)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(constants % ncav)
    real(dp), intent(in) :: gradphi_cav(3, constants % ncav)
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)
    real(dp), intent(in) :: hessianphi_cav(3, 3, constants % ncav)
    real(dp), intent(out) :: force(3, params % nsph)

    call ddlpb_force_worker(params, constants, workspace, &
        & hessianphi_cav, state % phi_grid, gradphi_cav, &
        & state % x_lpb, state % x_adj_lpb, state % zeta, force)
end subroutine ddlpb_force

subroutine ddlpb_guess(params, constants, state)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(inout) :: state
    ! TODO
end subroutine ddlpb_guess

!!
!! A standalone ddLPB calculation happens here
!! @param[in] ddx_data : dd Data
!! @param[in] phi      : Boundary conditions
!! @param[in] psi      : Electrostatic potential vector.
!! @param[in] gradphi  : Gradient of phi
!! @param[in] hessianphi  : Hessian of phi
!! @param[out] esolv   : Electrostatic solvation energy
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
    real(dp), external :: ddot

    ! TODO: find a consistent way to do the guess

    call ddlpb_solve(params, constants, workspace, state, phi_cav, &
        & gradphi_cav, tol)
    if (workspace % error_flag .eq. 1) return

    ! Compute the solvation energy
    ! note that psi is divided by fourpi
    esolv = pt5*ddot(constants % n, state % x_lpb(:,:,1), 1, psi, 1) &
        & /fourpi

    ! Get forces if needed
    if(params % force .eq. 1) then
        call ddlpb_adjoint(params, constants, workspace, state, psi, tol)
        if (workspace % error_flag .eq. 1) return
        call ddlpb_force(params, constants, workspace, state, phi_cav, &
            & gradphi_cav, hessianphi_cav, psi, force)
        if (workspace % error_flag .eq. 1) return
    endif

end subroutine ddlpb

subroutine ddlpb_solve_worker(params, constants, workspace, phi_cav, &
        & gradphi_cav, g_lpb, f_lpb, phi_grid, x_lpb, x_lpb_niter, &
        & x_lpb_time, x_lpb_rel_diff, tol)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: phi_cav(constants % ncav)
    real(dp), intent(in) :: gradphi_cav(3, constants % ncav)
    real(dp), intent(inout) :: x_lpb(constants % nbasis, params % nsph, 2)
    integer, intent(inout) :: x_lpb_niter
    real(dp), intent(out) :: x_lpb_time
    real(dp), intent(out) :: x_lpb_rel_diff(x_lpb_niter)
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: g_lpb(params % ngrid, params % nsph)
    real(dp), intent(out) :: f_lpb(params % ngrid, params % nsph)
    real(dp), intent(out) :: phi_grid(params % ngrid, params % nsph)

    integer :: istat
    real(dp) :: start_time

    real(dp), allocatable :: rhs(:,:,:)

    allocate(rhs(constants % nbasis, params % nsph, 2), stat=istat)
    if (istat.ne.0) then
        workspace % error_message = 'Allocation failed in ddlpb_solve_worker'
        workspace % error_flag = 1
        return
    end if

    ! set the initial convergence for the microiterations, it will be
    ! adjusted depending on the external one
    constants % inner_tol =  sqrt(tol)

    !! Setting initial values to zero
    g_lpb = zero
    f_lpb = zero
    phi_grid = zero

    ! Unwrap sparsely stored potential at cavity points phi_cav into phi_grid
    ! and multiply it by characteristic function at cavity points ui
    call ddcav_to_grid_work(params % ngrid, params % nsph, &
        & constants % ncav, constants % icav_ia, &
        & constants % icav_ja, phi_cav, phi_grid)
    workspace % tmp_cav = phi_cav * constants % ui_cav
    call ddcav_to_grid_work(params % ngrid, params % nsph, &
        & constants % ncav, constants % icav_ia, &
        & constants % icav_ja, workspace % tmp_cav, &
        & workspace % tmp_grid)
    g_lpb = - workspace % tmp_grid

    ! wghpot_f : Intermediate computation of F_0 Eq.(75) from QSM19.SISC
    call wghpot_f(params, constants, workspace, gradphi_cav, f_lpb)

    ! Setting of the local variables
    rhs = zero

    ! integrate RHS
    call ddintegrate(params % nsph, constants % nbasis, &
        & params % ngrid, constants % vwgrid, &
        & constants % vgrid_nbasis, g_lpb, rhs(:,:,1))
    call ddintegrate(params % nsph, constants % nbasis, &
        & params % ngrid, constants % vwgrid, &
        & constants % vgrid_nbasis, f_lpb, rhs(:,:,2))
    rhs(:,:,1) = rhs(:,:,1) + rhs(:,:,2)

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call prec_tx(params, constants, workspace, rhs, x_lpb)

    ! solve LS using Jacobi/DIIS
    start_time = omp_get_wtime()
    call jacobi_diis_external(params, constants, workspace, &
        & 2*constants % n, tol, rhs, x_lpb, x_lpb_niter, x_lpb_rel_diff, &
        & cx, prec_tx, rmsnorm)
    if (workspace % error_flag .ne. 0) then
        workspace % error_message = 'Jacobi solver failed to converge in ddlpb_solve_worker'
        return
    end if
    x_lpb_time = omp_get_wtime() - start_time

    deallocate(rhs, stat=istat)
    if (istat .ne. 0) then
        workspace % error_message = 'Deallocation failed in ddlpb_solve_worker'
        workspace % error_flag = 1
        return
    end if
end subroutine ddlpb_solve_worker

subroutine ddlpb_adjoint_worker(params, constants, workspace, psi, tol, &
        & x_adj_lpb, x_adj_lpb_niter, x_adj_lpb_time, x_adj_lpb_rel_diff)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: x_adj_lpb(constants % nbasis, params % nsph, 2)
    integer, intent(inout) :: x_adj_lpb_niter
    real(dp), intent(out) :: x_adj_lpb_time
    real(dp), intent(out) :: x_adj_lpb_rel_diff(x_adj_lpb_niter)

    real(dp), allocatable :: rhs(:,:,:)
    real(dp) :: start_time
    integer :: istat

    allocate(rhs(constants % nbasis, params % nsph, 2), stat=istat)
    if (istat.ne.0) then
        workspace % error_message = 'Allocation failed in ddlpb_adjoint_worker'
        workspace % error_flag = 1
        return
    end if

    !! Use a tighter tolerance for the microiterations to ensure convergence
    constants % inner_tol =  sqrt(tol)

    ! Psi shall be divided by a factor 4pi for the LPB case
    ! It is intended to take into account this constant in the LPB

    ! set up the RHS
    rhs(:,:,1) = psi/fourpi
    rhs(:,:,2) = zero

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call prec_tstarx(params, constants, workspace, rhs, x_adj_lpb)

    ! solve adjoint LS using Jacobi/DIIS
    start_time = omp_get_wtime()
    call jacobi_diis_external(params, constants, workspace, &
        & 2*constants % n, tol, rhs, x_adj_lpb, x_adj_lpb_niter, &
        & x_adj_lpb_rel_diff, cstarx, prec_tstarx, &
        & rmsnorm)
    if (workspace % error_flag .ne. 0) then
        workspace % error_message = 'Jacobi solver failed to converge in ddlpb_adjoint_worker'
        return
    end if
    x_adj_lpb_time = omp_get_wtime() - start_time
    deallocate(rhs, stat=istat)
    if (istat.ne.0) then
        workspace % error_message = 'Deallocation failed in ddlpb_adjoint_worker'
        workspace % error_flag = 1
        return
    end if
end subroutine ddlpb_adjoint_worker

subroutine ddlpb_force_worker(params, constants, workspace, hessian, &
        & phi_grid, gradphi, x, x_adj, zeta, force)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace

    real(dp), dimension(3, 3, constants % ncav), intent(in) :: hessian
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: phi_grid
    real(dp), dimension(3, constants % ncav), intent(in) :: gradphi
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x, x_adj
    real(dp), dimension(3, params % nsph), intent(out) :: force
    real(dp), intent(out) :: zeta(constants % ncav)

    ! local
    real(dp), dimension(constants % nbasis) :: basloc, vplm
    real(dp), dimension(3, constants % nbasis) :: dbasloc
    real(dp), dimension(params % lmax + 1) :: vsin, vcos

    ! large local are allocatable
    real(dp), allocatable :: ef(:,:), xadj_r_sgrid(:,:), xadj_e_sgrid(:,:), &
        & normal_hessian_cav(:,:), diff_re(:,:), scaled_xr(:,:)
    integer :: isph, icav, icav_gr, icav_ge, igrid, istat
    integer :: i, inode, jnear, inear, jnode, jsph
    real(dp), external :: ddot, dnrm2
    real(dp) :: tcontract_gradi_Lik, tcontract_gradi_Lji, &
        & tcontract_gradi_Bik, tcontract_gradi_Bji, tcontract_grad_U, &
        & tcontract_grad_C_worker1, tcontract_grad_C_worker2
    real(dp) :: d(3), dnorm, tmp1, tmp2

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
                        & hessian(:,i,icav)*constants % cgrid(i,igrid)
                end do
            end if
        end do
    end do

    ! Call dgemm to integrate the adjoint solution on the grid points
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & one, constants % vgrid, constants % vgrid_nbasis, x_adj(:,:,1), &
        & constants % nbasis, zero, Xadj_r_sgrid, params % ngrid)
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & one, constants % vgrid, constants % vgrid_nbasis, x_adj(:,:,2), &
        & constants % nbasis, zero, Xadj_e_sgrid, params % ngrid)

    ! Scale by the factor of 1/(4Pi/(2l+1))
    scaled_Xr = x(:,:,1)
    call convert_ddcosmo(params, constants, -1, scaled_Xr)

    !$omp parallel do default(none) shared(params,constants,workspace, &
    !$omp scaled_xr,xadj_r_sgrid,x,force,xadj_e_sgrid,phi_grid) &
    !$omp private(isph,basloc,dbasloc,vplm,vcos,vsin) &
    !$omp schedule(static,1)
    do isph = 1, params % nsph
        ! Compute A^k*Xadj_r, using Subroutine from ddCOSMO
        call contract_grad_L(params, constants, isph, scaled_Xr, Xadj_r_sgrid, &
            & basloc, dbasloc, vplm, vcos, vsin, force(:,isph))
        ! Compute B^k*Xadj_e
        call contract_grad_B(params, constants, workspace, isph, x(:,:,2), &
            & Xadj_e_sgrid, basloc, dbasloc, vplm, vcos, vsin, force(:, isph))
        ! Computation of G0
        call contract_grad_U(params, constants, isph, Xadj_r_sgrid, phi_grid, &
            & force(:, isph))
    end do
    ! Compute C1 and C2 contributions
    diff_re = zero
    call contract_grad_C(params, constants, workspace, x(:,:,1), x(:,:,2), &
        & Xadj_r_sgrid, Xadj_e_sgrid, x_adj(:,:,1), x_adj(:,:,2), force, &
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
                zeta(icav) = -constants % wgrid(igrid) * &
                    & constants % ui(igrid, isph) * ddot(constants % nbasis, &
                    & constants % vgrid(1, igrid), 1, &
                    & x_adj(1, isph, 1), 1)
                force(:, isph) = force(:, isph) + &
                    & zeta(icav)*gradphi(:, icav)
            end if
        end do
    end do
    !! Last term where we compute gradients of potential at centers of atoms
    !! spawned by intermediate zeta.
    if(params % fmm .eq. 1) then
        !! This step can be substituted by a proper dgemm if zeta
        !! intermediate is converted from cavity points to all grid points
        !! with zero values at internal grid points
        ! P2M step
        icav = 0
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_m(:, inode) = zero
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
                icav = icav + 1
                call fmm_p2m(constants % cgrid(:, igrid), zeta(icav), &
                    & one, params % pm, constants % vscales, one, &
                    & workspace % tmp_node_m(:, inode))
            end do
            workspace % tmp_node_m(:, inode) = workspace % tmp_node_m(:, inode) / &
                & params % rsph(isph)
        end do
        ! M2M, M2L and L2L translations
        call tree_m2m_rotation(params, constants, workspace % tmp_node_m)
        call tree_m2l_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call tree_l2l_rotation(params, constants, workspace % tmp_node_l)
        ! Now compute near-field FMM gradients
        ! Cycle over all spheres
        icav = 0
        workspace % tmp_efld = zero
        do isph = 1, params % nsph
            ! Cycle over all external grid points
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
                icav = icav + 1
                ! Cycle over all near-field admissible pairs of spheres,
                ! including pair (isph, isph) which is a self-interaction
                inode = constants % snode(isph)
                do jnear = constants % snear(inode), constants % snear(inode+1)-1
                    ! Near-field interactions are possible only between leaf
                    ! nodes, which must contain only a single input sphere
                    jnode = constants % near(jnear)
                    jsph = constants % order(constants % cluster(1, jnode))
                    d = params % csph(:, isph) + &
                        & constants % cgrid(:, igrid)*params % rsph(isph) - &
                        & params % csph(:, jsph)
                    dnorm = dnrm2(3, d, 1)
                    workspace % tmp_efld(:, jsph) = workspace % tmp_efld(:, jsph) + &
                        & zeta(icav)*d/(dnorm**3)
                end do
            end do
        end do
        ! Take into account far-field FMM gradients only if pl > 0
        if (params % pl .gt. 0) then
            tmp1 = one / sqrt(3d0) / constants % vscales(1)
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                tmp2 = tmp1 / params % rsph(isph)
                workspace % tmp_efld(3, isph) = workspace % tmp_efld(3, isph) + &
                    & tmp2*workspace % tmp_node_l(3, inode)
                workspace % tmp_efld(1, isph) = workspace % tmp_efld(1, isph) + &
                    & tmp2*workspace % tmp_node_l(4, inode)
                workspace % tmp_efld(2, isph) = workspace % tmp_efld(2, isph) + &
                    & tmp2*workspace % tmp_node_l(2, inode)
            end do
        end if
        do isph = 1, params % nsph
            force(:, isph) = force(:, isph) + &
                & workspace % tmp_efld(:, isph)*params % charge(isph)
        end do
    ! Naive quadratically scaling implementation
    else
        ! This routines actually computes -grad, not grad
        call efld(constants % ncav, zeta, constants % ccav, params % nsph, &
            & params % csph, workspace % tmp_efld)
        do isph = 1, params % nsph
            force(:, isph) = force(:, isph) - &
                & workspace % tmp_efld(:, isph)*params % charge(isph)
        end do
    end if

    icav_gr = zero
    icav_ge = zero
    ! Computation of F0
    call contract_grad_f(params, constants, workspace, x_adj(:,:,1), &
        & Xadj_r_sgrid, gradphi, normal_hessian_cav, icav_gr, force)
    if (workspace % error_flag .eq. 1) return
    call contract_grad_f(params, constants, workspace, x_adj(:,:,2), &
        & Xadj_e_sgrid, gradphi, normal_hessian_cav, icav_ge, force)
    if (workspace % error_flag .eq. 1) return

    force = - pt5*force

    deallocate(ef, xadj_r_sgrid, xadj_e_sgrid, normal_hessian_cav, &
        & diff_re, scaled_xr, stat=istat)
    if (istat.ne.0) then
        workspace % error_message = 'Deallocation failed in ddlpb_force_worker'
        workspace % error_flag = 1
        return
    end if
end subroutine ddlpb_force_worker

end module ddx_lpb
