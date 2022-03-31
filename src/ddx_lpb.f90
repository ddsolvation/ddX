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
!! ddLPB calculation happens here
!! @param[in] ddx_data    : dd Data
!! @param[in] phi_cav     : Boundary conditions
!! @param[in] gradphi_cav : Gradient of phi
!! @param[in] hessianphi  : Hessian of phi
!! @param[in] psi         : Electrostatic potential vector
!! @param[in] tol         : Tolerance
!! @param[out] esolv      : Electrostatic solvation energy
!! @param[out] force      : Forces
!! @param[out] info       : Information
!!
subroutine ddlpb(ddx_data, phi_cav, gradphi_cav, hessianphi_cav, psi, tol, esolv, &
    & force, info)
    ! main ddLPB
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), dimension(ddx_data % constants % ncav), intent(in) :: phi_cav
    real(dp), dimension(3, ddx_data % constants % ncav), intent(in) :: gradphi_cav
    real(dp), dimension(3,3, ddx_data % constants % ncav), intent(in) :: hessianphi_cav
    real(dp), dimension(ddx_data % constants % nbasis, &
        & ddx_data % params % nsph), intent(in) :: psi
    real(dp), intent(in) :: tol
    ! Outputs
    real(dp), intent(out) :: esolv
    real(dp), dimension(3, ddx_data % params % nsph), intent(out) :: force
    ! internal
    integer, intent(out) :: info
    integer :: isph, igrid, istatus
    !
    ! x(:,:,1): X_r Reaction potential solution (Laplace equation)
    ! x(:,:,2): X_e Extended potential solution (HSP equation)
    !

    !! Setting initial values to zero
    ddx_data % g_lpb = zero
    ddx_data % f_lpb = zero
    ddx_data % phi_grid = zero

    ! Unwrap sparsely stored potential at cavity points phi_cav into phi_grid
    ! and multiply it by characteristic function at cavity points ui
    call ddcav_to_grid_work(ddx_data % params % ngrid, ddx_data % params % nsph, &
        & ddx_data % constants % ncav, ddx_data % constants % icav_ia, &
        & ddx_data % constants % icav_ja, phi_cav, ddx_data % phi_grid)
    ddx_data % workspace % tmp_cav = phi_cav * ddx_data % constants % ui_cav
    call ddcav_to_grid_work(ddx_data % params % ngrid, ddx_data % params % nsph, &
        & ddx_data % constants % ncav, ddx_data % constants % icav_ia, &
        & ddx_data % constants % icav_ja, ddx_data % workspace % tmp_cav, &
        & ddx_data % workspace % tmp_grid)
    ddx_data % g_lpb = - ddx_data % workspace % tmp_grid

    ! wghpot_f : Intermediate computation of F_0 Eq.(75) from QSM19.SISC
    call wghpot_f(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, gradphi_cav, ddx_data % f_lpb)

    ! use a tighter tolerance for microiterations
    ddx_data % constants % inner_tol =  tol/100.0d0

    ddx_data % xs_niter = ddx_data % params % maxiter
    ! Call the subroutine to solve for Esolv
    call ddlpb_energy(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, ddx_data % g_lpb, ddx_data % f_lpb, &
        & ddx_data % x_lpb, ddx_data % xs_niter, ddx_data % xs_rel_diff, ddx_data % xs_time, &
        & tol, esolv, info)

    ! Start the Force computation
    if(ddx_data % params % force .eq. 1) then
      ddx_data % s_niter = ddx_data % params % maxiter
      ! Call the subroutine adjoint to solve the adjoint solution
      call ddlpb_adjoint(ddx_data % params, ddx_data % constants, &
          & ddx_data % workspace, psi, tol, ddx_data % x_adj_lpb, &
          ddx_data % s_niter, ddx_data % s_rel_diff, ddx_data % s_time, info)

      !Call the subroutine to evaluate derivatives
      call ddlpb_force(ddx_data % params, ddx_data % constants, &
          & ddx_data % workspace, hessianphi_cav, ddx_data % phi_grid, gradphi_cav, &
          & ddx_data % x_lpb, ddx_data % x_adj_lpb, ddx_data % zeta, force, info)

    endif

!   deallocate(x, x_adj, g, f, phi_grid, stat = istatus)
!   if (istatus.ne.0) write(6,*) 'ddlpb deallocation failed'

end subroutine ddlpb

!> Solve primal ddLPB system to find solvation energy
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  workspace  : Input workspace
!! @param[in]  g          : Intermediate matrix for computation of g0
!! @param[in]  f          : Intermediate matrix for computation of f0
!! @param[out] x          : Solution
!! @param[out] n_iter     : Number of iterations
!! @param[out] xs_rel_diff: Relative difference
!! @param[out] xs_time    : Time
!! @param[in]  tol        : Tolerance for solver
!! @param[out] esolv      : Solvation energy
!! @param[out] info       : Information
subroutine ddlpb_energy(params, constants, workspace, g, f, &
    & x, n_iter, xs_rel_diff, xs_time, tol , esolv, info)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: g, f
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: x
    real(dp), intent(in) :: tol
    integer, intent(inout) :: n_iter
    !! Outputs
    real(dp), intent(out) :: xs_rel_diff(n_iter), xs_time
    real(dp), intent(out) :: esolv
    integer, intent(out) :: info
    integer  :: isph, istat
    real(dp) :: start_time, finish_time
    real(dp), allocatable :: rhs(:,:,:)

    allocate(rhs(constants % nbasis, params % nsph, 2), stat=istat)
    if (istat.ne.0) stop 1

    ! Setting of the local variables
    rhs = zero

    ! integrate RHS
    call ddintegrate(params % nsph, constants % nbasis, &
        & params % ngrid, constants % vwgrid, &
        & constants % vgrid_nbasis, g, rhs(:,:,1))
    call ddintegrate(params % nsph, constants % nbasis, &
        & params % ngrid, constants % vwgrid, &
        & constants % vgrid_nbasis, f, rhs(:,:,2))
    rhs(:,:,1) = rhs(:,:,1) + rhs(:,:,2)

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call lpb_direct_prec(params, constants, workspace, rhs, x)

    ! solve LS using Jacobi/DIIS
    call cpu_time(start_time)
    call jacobi_diis_external(params, constants, workspace, 2*constants % n, &
        & tol, rhs, x, n_iter, xs_rel_diff, lpb_direct_matvec, &
        & lpb_direct_prec, rmsnorm, info)
    call cpu_time(finish_time)
    xs_time = finish_time - start_time

    esolv = zero
    do isph = 1, params % nsph
      esolv = esolv + pt5*params % charge(isph)*x(1,isph,1)*(one/sqrt4pi)
    end do
    deallocate(rhs)
    if (istat.ne.0) stop 1
end subroutine ddlpb_energy

!> Solve adjoint ddLPB system
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  workspace  : Input workspace
!! @param[in]  psi        : psi_r
!! @param[in]  tol        : Tolerance
!! @param[out] x_adj      : Adjoint solution
!! @param[out] info       : Information

subroutine ddlpb_adjoint(params, constants, workspace, psi, tol, x_adj, &
                       & n_iter, s_rel_diff, s_time, info)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: psi
    real(dp), intent(in) :: tol
    !!Output
    integer, intent(inout) :: n_iter
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: x_adj
    real(dp), intent(out) :: s_time
    real(dp), intent(out) :: s_rel_diff(n_iter)
    integer, intent(out) :: info

    real(dp), allocatable :: rhs(:,:,:)
    integer :: istat
    real(dp) :: start_time, finish_time

    allocate(rhs(constants % nbasis, params % nsph, 2), stat=istat)
    if (istat.ne.0) stop 1

    ! set up the RHS
    rhs(:,:,1) = psi
    rhs(:,:,2) = zero

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call lpb_adjoint_prec(params, constants, workspace, rhs, x_adj)

    call cpu_time(start_time)
    ! solve adjoint LS using Jacobi/DIIS
    call jacobi_diis_external(params, constants, workspace, 2*constants % n, &
        & tol, rhs, x_adj, n_iter, s_rel_diff, lpb_adjoint_matvec, &
        & lpb_adjoint_prec, rmsnorm, info)
    call cpu_time(finish_time)
    s_time = finish_time - start_time

    deallocate(rhs)
    if (istat.ne.0) stop 1

end subroutine ddlpb_adjoint

!> Compute ddLPB forces
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  workspace  : Input workspace
!! @param[in]  hessian    : Hessian of Psi
!! @param[in]  phi_grid   : Phi evaluated at the grid point
!! @param[in]  gradphi    : Gradient of phi
!! @param[in]  x          : Solution vector
!! @param[in]  x_adj      : Adjoint solution vector
!! @param[in]  zeta       : Intermediate calculation
!! @param[out] force      : Force
!! @param[out] info       : Information
subroutine ddlpb_force(params, constants, workspace, hessian, phi_grid, gradphi, &
        & x, x_adj, zeta, force, info)
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace

    real(dp), dimension(3, 3, constants % ncav), intent(in) :: hessian
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: phi_grid
    real(dp), dimension(3, constants % ncav), intent(in)    :: gradphi
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x, x_adj
    real(dp), dimension(3, params % nsph), intent(out) :: force
    real(dp), intent(out) :: zeta(constants % ncav)
    integer, intent(out) :: info

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
    real(dp) :: tcontract_gradi_Lik, tcontract_gradi_Lji, tcontract_gradi_Bik, tcontract_gradi_Bji, tcontract_grad_U, &
              & tcontract_grad_C_worker1, tcontract_grad_C_worker2
    real(dp) :: d(3), dnorm, tmp1, tmp2

    allocate(ef(3, params % nsph), &
        & xadj_r_sgrid(params % ngrid, params % nsph), &
        & xadj_e_sgrid(params % ngrid, params % nsph), &
        & normal_hessian_cav(3, constants % ncav), &
        & diff_re(constants % nbasis, params % nsph), &
        & scaled_xr(constants % nbasis, params % nsph), stat=istat)
    if (istat.ne.0) stop 1

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

    tcontract_gradi_Lik = zero
    tcontract_gradi_Lji = zero
    tcontract_gradi_Bik = zero
    tcontract_gradi_Bji = zero
    tcontract_grad_C_worker2 = zero
    tcontract_grad_C_worker1 = zero
    tcontract_grad_U = zero
    do isph = 1, params % nsph
        ! Compute A^k*Xadj_r, using Subroutine from ddCOSMO
        call contract_grad_L(params, constants, isph, scaled_Xr, Xadj_r_sgrid, &
            & basloc, dbasloc, vplm, vcos, vsin, force(:,isph))
        ! Compute B^k*Xadj_e
        call contract_grad_B(params, constants, workspace, isph, x(:,:,2), Xadj_e_sgrid, &
            & basloc, dbasloc, vplm, vcos, vsin, force(:, isph))
        ! Computation of G0
        call contract_grad_U(params, constants, isph, Xadj_r_sgrid, phi_grid, &
            & force(:, isph))
    end do
    ! Compute C1 and C2 contributions
    diff_re = zero
    call contract_grad_C(params, constants, workspace, x(:,:,1), x(:,:,2), Xadj_r_sgrid, &
        & Xadj_e_sgrid, x_adj(:,:,1), x_adj(:,:,2), force, diff_re)
    ! Computation of G0 continued
    ! NOTE: contract_grad_U returns a positive summation
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

    tcontract_grad_C_worker2 = zero
    tcontract_grad_U = zero
    icav_gr = zero
    icav_ge = zero
    ! Computation of F0
    call contract_grad_f(params, constants, workspace, x_adj(:,:,1), Xadj_r_sgrid, &
                  & gradphi, normal_hessian_cav, icav_gr, force)
    call contract_grad_f(params, constants, workspace, x_adj(:,:,2), Xadj_e_sgrid, &
                  & gradphi, normal_hessian_cav, icav_ge, force)

    force = pt5*force

    deallocate(ef, xadj_r_sgrid, xadj_e_sgrid, normal_hessian_cav, &
        & diff_re, scaled_xr, stat=istat)
    if (istat.ne.0) stop 1
end subroutine ddlpb_force

!> Free ddLPB constants
!!
subroutine ddlpb_free(ddx_data)
  implicit none
  type(ddx_type), intent(inout) :: ddx_data
  integer :: istatus
  if(allocated(ddx_data % constants % SI_ri)) then
    deallocate(ddx_data % constants % SI_ri, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % DI_ri)) then
    deallocate(ddx_data % constants % DI_ri, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % SK_ri)) then
    deallocate(ddx_data % constants % SK_ri, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % DK_ri)) then
    deallocate(ddx_data % constants % DK_ri, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if


  if(allocated(ddx_data % constants % diff_ep_adj)) then
    deallocate(ddx_data % constants % diff_ep_adj, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if


  if(allocated(ddx_data % constants % termimat)) then
    deallocate(ddx_data % constants % termimat, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % coefY)) then
    deallocate(ddx_data % constants % coefY, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % C_ik)) then
    deallocate(ddx_data % constants % C_ik, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % coefvec)) then
    deallocate(ddx_data % constants % coefvec, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % Pchi)) then
    deallocate(ddx_data % constants % Pchi, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

end subroutine ddlpb_free

end module ddx_lpb
