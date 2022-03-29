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

real(dp) :: t0, t1, tt0, tt1

contains

!!
!! ddLPB calculation happens here
!! @param[in] ddx_data : dd Data
!! @param[in] phi      : Boundary conditions
!! @param[in] psi      : Electrostatic potential vector.
!! @param[in] gradphi  : Gradient of phi
!! @param[in] hessianphi  : Hessian of phi
!! @param[out] esolv   : Electrostatic solvation energy
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
    real(dp), allocatable :: x(:,:,:), x_adj(:,:,:)
    ! phi_grid: Phi evaluated at grid points
    real(dp), allocatable :: g(:,:), f(:,:), phi_grid(:, :)

    allocate(x(ddx_data % constants % nbasis, ddx_data % params % nsph, 2), &
        & x_adj(ddx_data % constants % nbasis, ddx_data % params % nsph, 2), &
        & g(ddx_data % params % ngrid, ddx_data % params % nsph), &
        & f(ddx_data % params % ngrid, ddx_data % params % nsph), &
        & phi_grid(ddx_data % params % ngrid, ddx_data % params % nsph), &
        & stat = istatus)
    if (istatus.ne.0) write(6,*) 'ddlpb allocation failed'

    !! Setting initial values to zero
    g = zero
    f = zero
    phi_grid = zero

    t0 = omp_get_wtime()
    ! Unwrap sparsely stored potential at cavity points phi_cav into phi_grid
    ! and multiply it by characteristic function at cavity points ui
    call ddcav_to_grid_work(ddx_data % params % ngrid, ddx_data % params % nsph, &
        & ddx_data % constants % ncav, ddx_data % constants % icav_ia, &
        & ddx_data % constants % icav_ja, phi_cav, phi_grid)
    ddx_data % workspace % tmp_cav = phi_cav * ddx_data % constants % ui_cav
    call ddcav_to_grid_work(ddx_data % params % ngrid, ddx_data % params % nsph, &
        & ddx_data % constants % ncav, ddx_data % constants % icav_ia, &
        & ddx_data % constants % icav_ja, ddx_data % workspace % tmp_cav, &
        & ddx_data % workspace % tmp_grid)
    g = - ddx_data % workspace % tmp_grid
    t1 = omp_get_wtime()
    write(6,*) '@wghpot', t1 - t0

    ! wghpot_f : Intermediate computation of F_0 Eq.(75) from QSM19.SISC
    t0 = omp_get_wtime()
    call wghpot_f(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, gradphi_cav,f)
    t1 = omp_get_wtime()
    write(6,*) '@wghpot_f', t1 - t0

    ! use a tighter tolerance for microiterations
    ddx_data % constants % inner_tol =  tol/100.0d0

    ! Call the subroutine to solve for Esolv
    t0 = omp_get_wtime()
    call ddx_lpb_solve(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, g, f, x, tol, esolv)
    t1 = omp_get_wtime()
    write(6,*) '@direct_ls', t1 - t0

    ! Start the Force computation
    if(ddx_data % params % force .eq. 1) then
      write(*,*) 'Computation of Forces for ddLPB'
      ! Call the subroutine adjoint to solve the adjoint solution
      t0 = omp_get_wtime()
      call ddx_lpb_adjoint(ddx_data % params, ddx_data % constants, &
          & ddx_data % workspace, psi, tol, x_adj)
      t1 = omp_get_wtime()

      !Call the subroutine to evaluate derivatives
      t0 = omp_get_wtime()
      call ddx_lpb_force(ddx_data % params, ddx_data % constants, &
          & ddx_data % workspace, hessianphi_cav, phi_grid, gradphi_cav, &
          & x, x_adj, ddx_data % zeta, force)
      t1 = omp_get_wtime()
      write(6,*) '@forces', t1 - t0

    endif

    deallocate(x, x_adj, g, f, phi_grid, stat = istatus)
    if (istatus.ne.0) write(6,*) 'ddlpb deallocation failed'

end subroutine ddlpb

!
! Computation for Solvation energy
! @param[in]  ddx_data   : Input data file
! @param[in]  g          : Intermediate matrix for computation of g0
! @param[in]  f          : Intermediate matrix for computation of f0
! @param[out] x          : Solution
! @param[out] esolv      : Solvation energy
subroutine ddx_lpb_solve(params, constants, workspace, g, f, &
    & x, tol, esolv)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: g, f
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: x
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: esolv
    integer  :: n_iter, isph, istat, info
    real(dp), allocatable :: rhs(:,:,:)

    allocate(rhs(constants % nbasis, params % nsph, 2), stat=istat)
    if (istat.ne.0) stop 1

    ! Setting of the local variables
    rhs = zero

    ! integrate RHS
    tt0 = omp_get_wtime()
    call intrhs(params % nsph, constants % nbasis, &
        & params % ngrid, constants % vwgrid, &
        & constants % vgrid_nbasis, g, rhs(:,:,1))
    call intrhs(params % nsph, constants % nbasis, &
        & params % ngrid, constants % vwgrid, &
        & constants % vgrid_nbasis, f, rhs(:,:,2))
    tt1 = omp_get_wtime()
    rhs(:,:,1) = rhs(:,:,1) + rhs(:,:,2)

    ! call prtsph('direct rhs', constants % nbasis, params % lmax, &
    !    & 2*params % nsph, 0, rhs)

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call lpb_direct_prec(params, constants, workspace, rhs, x)

    ! solve LS using Jacobi/DIIS
    n_iter = params % maxiter
    call jacobi_diis_external(params, constants, workspace, 2*constants % n, &
        & tol, rhs, x, n_iter, lpb_direct_matvec, lpb_direct_prec, rmsnorm, info)

    ! call prtsph('direct sol', constants % nbasis, params % lmax, &
    !    & 2*params % nsph, 0, x)

    esolv = zero
    do isph = 1, params % nsph
      esolv = esolv + pt5*params % charge(isph)*x(1,isph,1)*(one/sqrt4pi)
    end do
    deallocate(rhs)
    if (istat.ne.0) stop 1
end subroutine ddx_lpb_solve
!
! Computation of Adjoint
! @param[in] ddx_data: Input data file
! @param[in] psi     : psi_r
subroutine ddx_lpb_adjoint(params, constants, workspace, psi, tol, x_adj)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: psi
    real(dp), intent(in) :: tol
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: x_adj
    real(dp), allocatable :: rhs(:,:,:)
    integer :: n_iter, info, istat
    real(dp), dimension(params % maxiter) :: x_rel_diff

    allocate(rhs(constants % nbasis, params % nsph, 2), stat=istat)
    if (istat.ne.0) stop 1

    ! set up the RHS
    rhs(:,:,1) = psi
    rhs(:,:,2) = zero

    ! call prtsph('adjoint rhs', constants % nbasis, params % lmax, &
    !    & 2*params % nsph, 0, rhs)

    ! guess
    workspace % ddcosmo_guess = zero
    workspace % hsp_guess = zero
    call lpb_adjoint_prec(params, constants, workspace, rhs, x_adj)

    ! solve adjoint LS using Jacobi/DIIS
    n_iter = params % maxiter
    call jacobi_diis_external(params, constants, workspace, 2*constants % n, &
        & tol, rhs, x_adj, n_iter, lpb_adjoint_matvec, lpb_adjoint_prec, rmsnorm, info)

    ! call prtsph('adjoint sol', constants % nbasis, params % lmax, &
    !    & 2*params % nsph, 0, x_adj)

    deallocate(rhs)
    if (istat.ne.0) stop 1

end subroutine ddx_lpb_adjoint

!
! Computation for Solvation energy
! @param[in]  ddx_data   : Input data file
! @param[in]  hessian    : Hessian of Psi
! @param[in]  phi_grid   : Phi evaluated at the grid point
! @param[in]  gradphi  : Gradient of phi
! @param[in]  Xr         : Solution corresponding to COSMO
! @param[in]  Xe         : Solution corresponding to HSP
! @param[in]  Xadj_r     : Solution corresponding to adjoint of the COSMO
! @param[in]  Xadj_e     : Solution corresponding to adjoint of the HSP
! @param[out] force      : Force
subroutine ddx_lpb_force(params, constants, workspace, hessian, phi_grid, gradphi, &
        & x, x_adj, zeta, force)
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
    real(dp) :: tfdoka, tfdokb, tfdoka_xe, tfdokb_xe, tfdoga, tp, tfdouky
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
    tt0 = omp_get_wtime()
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
    write(6,*) '@forces@normal_psi', omp_get_wtime() - tt0

    ! Call dgemm to integrate the adjoint solution on the grid points
    tt0 = omp_get_wtime()
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & one, constants % vgrid, constants % vgrid_nbasis, x_adj(:,:,1), &
        & constants % nbasis, zero, Xadj_r_sgrid, params % ngrid)
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & one, constants % vgrid, constants % vgrid_nbasis, x_adj(:,:,2), &
        & constants % nbasis, zero, Xadj_e_sgrid, params % ngrid)
    write(6,*) '@forces@xi', omp_get_wtime() - tt0

    ! Scale by the factor of 1/(4Pi/(2l+1))
    scaled_Xr = x(:,:,1)
    call convert_ddcosmo(params, constants, -1, scaled_Xr)

    tfdoka = zero
    tfdokb = zero
    tfdoka_xe = zero
    tfdokb_xe = zero
    tfdouky = zero
    tp = zero
    tfdoga = zero
    do isph = 1, params % nsph
        ! Compute A^k*Xadj_r, using Subroutine from ddCOSMO
        tt1 = omp_get_wtime()
        call fdoka(params, constants, isph, scaled_Xr, Xadj_r_sgrid(:, isph), &
            & basloc, dbasloc, vplm, vcos, vsin, force(:,isph))
        tfdoka = tfdoka + omp_get_wtime() - tt1
        tt1 = omp_get_wtime()
        call fdokb(params, constants, isph, scaled_Xr, Xadj_r_sgrid, basloc, &
            & dbasloc, vplm, vcos, vsin, force(:, isph))
        tfdokb = tfdokb + omp_get_wtime() - tt1
        ! Compute B^k*Xadj_e
        tt1 = omp_get_wtime()
        call fdoka_b_xe(params, constants, workspace, isph, x(:,:,2), &
            & Xadj_e_sgrid(:, isph), basloc, dbasloc, vplm, vcos, vsin, &
            & force(:,isph))
        tfdoka_xe = tfdoka_xe + omp_get_wtime() - tt1
        tt1 = omp_get_wtime()
        call fdokb_b_xe(params, constants, workspace, isph, x(:,:,2), Xadj_e_sgrid, &
            & basloc, dbasloc, vplm, vcos, vsin, force(:, isph))
        tfdokb_xe = tfdokb_xe + omp_get_wtime() - tt1
        ! Computation of G0
        tt1 = omp_get_wtime()
        call fdoga(params, constants, isph, Xadj_r_sgrid, phi_grid, &
            & force(:, isph))
        tfdoga = tfdoga + omp_get_wtime() - tt1
    end do
    ! Compute C1 and C2 contributions
    tt1 = omp_get_wtime()
    diff_re = zero
    call fdouky(params, constants, workspace, x(:,:,1), x(:,:,2), Xadj_r_sgrid, &
        & Xadj_e_sgrid, x_adj(:,:,1), x_adj(:,:,2), force, diff_re)
    tfdouky = tfdouky + omp_get_wtime() - tt1
    tt1 = omp_get_wtime()
    call derivative_P(params, constants, workspace, x(:,:,1), x(:,:,2), Xadj_r_sgrid, &
        & Xadj_e_sgrid, diff_re, force)
    tp = tp + omp_get_wtime() - tt1
    write(6,*) '@forces@fdoka', tfdoka
    write(6,*) '@forces@fdokb', tfdokb
    write(6,*) '@forces@fdoka_xe', tfdoka_xe
    write(6,*) '@forces@fdoka_xe', tfdokb_xe
    write(6,*) '@forces@fdouky', tfdouky
    write(6,*) '@forces@p', tp
    write(6,*) '@forces@fdoga', tfdoga
    ! Computation of G0 continued
    ! NOTE: fdoga returns a positive summation
    force = -force
    tt0 = omp_get_wtime()
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
    write(6,*) '@forces@rhs', omp_get_wtime() - tt0

    tfdouky = zero
    tfdoga = zero
    icav_gr = zero
    icav_ge = zero
    ! Computation of F0
    tt1 = omp_get_wtime()
    call fdouky_f0(params, constants, workspace, x_adj(:,:,1), Xadj_r_sgrid, &
        & gradphi, force)
    tfdouky = tfdouky + omp_get_wtime() - tt1
    ! Computation of F0
    tt1 = omp_get_wtime()
    call fdouky_f0(params, constants, workspace, x_adj(:,:,2), Xadj_e_sgrid, &
        & gradphi, force)
    tfdouky = tfdouky + omp_get_wtime() - tt1
    tt1 = omp_get_wtime()
    call fdoco(params, constants, workspace, Xadj_r_sgrid, gradphi, &
        & normal_hessian_cav, icav_gr, force)
    tfdoga = tfdoga + omp_get_wtime() - tt1
    tt1 = omp_get_wtime()
    call fdoco(params, constants, workspace, Xadj_e_sgrid, gradphi, &
        & normal_hessian_cav, icav_ge, force)
    tfdoga = tfdoga + omp_get_wtime() - tt1
    write(6,*) '@forces@fdouky', tfdouky
    write(6,*) '@forces@fdoco', tfdoga

    force = pt5*force

    deallocate(ef, xadj_r_sgrid, xadj_e_sgrid, normal_hessian_cav, &
        & diff_re, scaled_xr, stat=istat)
    if (istat.ne.0) stop 1
end subroutine ddx_lpb_force

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
