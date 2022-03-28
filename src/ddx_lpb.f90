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

module ddx_lpb
use ddx_core
use ddx_operators
use ddx_solvers
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


!!
!! Find intermediate F0 in the RHS of the ddLPB model given in Eq.(82)
!! @param[in]  gradphi : Gradient of psi_0
!! @param[out] f       : Intermediate calculation of F0
!!
subroutine wghpot_f(params, constants, workspace, gradphi, f)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: gradphi(3, constants % ncav)
    !! Temporary buffers
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Outputs
    real(dp), intent(out) :: f(params % ngrid, params % nsph)
    !! Local variables
    integer :: isph, ig, ic, ind, ind0, jg, l, m, jsph
    real(dp) :: nderphi, sumSijn, rijn, coef_Ylm, sumSijn_pre, termi, &
        & termk, term, tmp1, tmp2
    real(dp), dimension(3) :: sijn, vij, vtij
    real(dp) :: rho, ctheta, stheta, cphi, sphi, start_time, finish_time
    real(dp), allocatable :: SK_rijn(:), DK_rijn(:), c0(:,:), c1(:,:)
    integer :: l0, m0, icav, istatus, indl, inode
    complex(dp) :: work_complex(constants % lmax0 + 1)
    real(dp) :: work(constants % lmax0 + 1)


    allocate(SK_rijn(0:constants % lmax0), DK_rijn(0:constants % lmax0), &
        & c0(constants % nbasis0, params % nsph), c1(constants % nbasis0, params % nsph))

    ic = 0
    f = zero
    c0 = zero
    c1 = zero

    ! Compute c0 Eq.(98) QSM19.SISC
    do isph = 1, params % nsph
        do ig = 1, params % ngrid
            if (constants % ui(ig, isph) .ne. zero) then
                ic = ic + 1
                nderphi = dot_product(gradphi(:, ic), constants % cgrid(:,ig))
                c0(:, isph) = c0(:, isph) + &
                    & constants % wgrid(ig)*constants % ui(ig,isph)*&
                    & nderphi*constants % vgrid(1:constants % nbasis0,ig)
                do l0 = 0, constants % lmax0
                    ind0 = l0*l0 + 1
                    ind = ind0 + 2*l0
                    c1(ind0:ind, isph) = constants % C_ik(l0, isph) * &
                        & c0(ind0:ind, isph)
                end do
            end if
        end do
    end do

    ! Computation of F0 using above terms
    ! icav: External grid poitns
    if (params % fmm .eq. 0) then
        icav = 0
        do isph = 1, params % nsph
            do ig = 1, params % ngrid
                if (constants % ui(ig,isph).gt.zero) then
                    icav = icav + 1
                    sumSijn = zero
                    ! Loop to compute Sijn
                    do jsph = 1, params % nsph
                        sumSijn_pre = sumSijn
                        vij  = params % csph(:,isph) + &
                            & params % rsph(isph)*constants % cgrid(:,ig) - &
                            & params % csph(:,jsph)
                        vtij = vij*params % kappa
                        call fmm_m2p_bessel_work(vtij, &
                            & constants % lmax0, constants % vscales, &
                            & constants % SK_ri(:, jsph), one, &
                            & c1(:, jsph), one, sumSijn, work_complex, work)
                    end do
                    ! Here Intermediate value of F_0 is computed Eq. (99)
                    ! Mutilplication with Y_lm and weights will happen afterwards
                    !write(6,*) sumSijn, epsp, eps, ddx_data % constants % ui(ig,isph)
                    f(ig,isph) = -(params % epsp/params % eps)*constants % ui(ig,isph) * sumSijn
                end if
            end do
        end do
    else
        ! Load input harmonics into tree data
        workspace % tmp_node_m = zero
        workspace % tmp_node_l = zero
        workspace % tmp_sph = zero
        workspace % tmp_sph(1:constants % nbasis0, :) = c1(:, :)
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(1:constants % nbasis0, inode) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph)
                workspace % tmp_node_m(constants % nbasis0+1:, inode) = zero
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(:, inode) = workspace % tmp_sph(1:indl, isph)
            end do
        end if
        ! Do FMM operations
        call tree_m2m_bessel_rotation(params, constants, workspace % tmp_node_m)
        call tree_m2l_bessel_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation(params, constants, workspace % tmp_node_l)
        call tree_l2p_bessel(params, constants, one, workspace % tmp_node_l, zero, &
            & workspace % tmp_grid)
        call tree_m2p_bessel(params, constants, constants % lmax0, one, &
            & params % lmax, workspace % tmp_sph, one, &
            & workspace % tmp_grid)
        f = -(params % epsp/params % eps) * constants % ui * workspace % tmp_grid
    end if

    deallocate(SK_rijn, DK_rijn, c0, c1)

end subroutine wghpot_f

subroutine lx_nodiag_incore(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    integer :: isph, jsph, ij
    y = zero
    !$omp parallel do default(none) shared(params,constants,x,y) &
    !$omp private(isph,ij,jsph)
    do isph = 1, params % nsph
        do ij = constants % inl(isph), constants % inl(isph + 1) - 1
            jsph = constants % nl(ij)
            call dgemv('n', constants % nbasis, constants % nbasis, one, &
                & constants % l(:,:,ij), constants % nbasis, x(:,jsph), 1, &
                & one, y(:,isph), 1)
        end do
    end do
end subroutine lx_nodiag_incore

subroutine lstarx_nodiag_incore(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    integer :: isph, jsph, ij, indmat
    y = zero
    !$omp parallel do default(none) shared(params,constants,x,y) &
    !$omp private(isph,ij,jsph,indmat)
    do isph = 1, params % nsph
        do ij = constants % inl(isph), constants % inl(isph + 1) - 1
            jsph = constants % nl(ij)
            indmat = constants % itrnl(ij)
            call dgemv('t', constants % nbasis, constants % nbasis, one, &
                & constants % l(:,:,indmat), constants % nbasis, x(:,jsph), 1, &
                & one, y(:,isph), 1)
        end do
    end do
end subroutine lstarx_nodiag_incore

subroutine bx_incore(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    integer :: isph, jsph, ij
    y = x
    !$omp parallel do default(none) shared(params,constants,x,y) &
    !$omp private(isph,ij,jsph)
    do isph = 1, params % nsph
        do ij = constants % inl(isph), constants % inl(isph + 1) - 1
            jsph = constants % nl(ij)
            call dgemv('n', constants % nbasis, constants % nbasis, one, &
                & constants % b(:,:,ij), constants % nbasis, x(:,jsph), 1, &
                & one, y(:,isph), 1)
        end do
    end do
end subroutine bx_incore

subroutine bx_nodiag_incore(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    integer :: isph, jsph, ij
    y = zero
    !$omp parallel do default(none) shared(params,constants,x,y) &
    !$omp private(isph,ij,jsph)
    do isph = 1, params % nsph
        do ij = constants % inl(isph), constants % inl(isph + 1) - 1
            jsph = constants % nl(ij)
            call dgemv('n', constants % nbasis, constants % nbasis, one, &
                & constants % b(:,:,ij), constants % nbasis, x(:,jsph), 1, &
                & one, y(:,isph), 1)
        end do
    end do
end subroutine bx_nodiag_incore

subroutine bstarx_incore(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    integer :: isph, jsph, ij, indmat
    y = x
    !$omp parallel do default(none) shared(params,constants,x,y) &
    !$omp private(isph,ij,jsph,indmat)
    do isph = 1, params % nsph
        do ij = constants % inl(isph), constants % inl(isph + 1) - 1
            jsph = constants % nl(ij)
            indmat = constants % itrnl(ij)
            call dgemv('t', constants % nbasis, constants % nbasis, one, &
                & constants % b(:,:,indmat), constants % nbasis, x(:,jsph), 1, &
                & one, y(:,isph), 1)
        end do
    end do
end subroutine bstarx_incore

subroutine bstarx_nodiag_incore(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    integer :: isph, jsph, ij, indmat
    y = zero
    !$omp parallel do default(none) shared(params,constants,x,y) &
    !$omp private(isph,ij,jsph,indmat)
    do isph = 1, params % nsph
        do ij = constants % inl(isph), constants % inl(isph + 1) - 1
            jsph = constants % nl(ij)
            indmat = constants % itrnl(ij)
            call dgemv('t', constants % nbasis, constants % nbasis, one, &
                & constants % b(:,:,indmat), constants % nbasis, x(:,jsph), 1, &
                & one, y(:,isph), 1)
        end do
    end do
end subroutine bstarx_nodiag_incore

!
! Subroutine used for the GMRES solver
! @param[in]      n : Size of the matrix
! @param[in]      x : Input vector
! @param[in, out] y : y=A*x
!
subroutine bx(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    integer :: isph, iproc

    y = zero
    !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
    !$omp private(isph,iproc) schedule(static,1)
    do isph = 1, params % nsph
      iproc = omp_get_thread_num() + 1
      call calcv2_lpb(params, constants, isph, workspace % tmp_pot(:, iproc), x, &
          & workspace % tmp_vylm, workspace % tmp_vplm, workspace % tmp_vcos, &
          & workspace % tmp_vsin, workspace % tmp_bessel)
      call intrhs(1, constants % nbasis, params % ngrid, constants % vwgrid, &
          & constants % vgrid_nbasis, workspace % tmp_pot(:, iproc), y(:,isph))
      y(:,isph) = - y(:,isph) + x(:,isph)
    end do
end subroutine bx

subroutine bx_nodiag(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    integer :: isph, iproc

    y = zero
    !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
    !$omp private(isph,iproc) schedule(static,1)
    do isph = 1, params % nsph
      iproc = omp_get_thread_num() + 1
      call calcv2_lpb(params, constants, isph, workspace % tmp_pot(:, iproc), x, &
          & workspace % tmp_vylm, workspace % tmp_vplm, workspace % tmp_vcos, &
          & workspace % tmp_vsin, workspace % tmp_bessel)
      call intrhs(1, constants % nbasis, params % ngrid, constants % vwgrid, &
          & constants % vgrid_nbasis, workspace % tmp_pot(:, iproc), y(:,isph))
      y(:,isph) = - y(:,isph)
    end do
end subroutine bx_nodiag

subroutine bx_prec(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    y = x
end subroutine bx_prec


!!
!! Scale the ddCOSMO solution vector
!! @param[in]      direction : Direction of the scaling
!! @param[in, out] vector    : ddCOSMO solution vector
!!
subroutine convert_ddcosmo(params, constants, direction, vector)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: direction
    !! Outputs
    real(dp), intent(inout) :: vector(constants % nbasis, params % nsph)
    !! Local variables
    integer :: isph, l, m, ind
    real(dp) :: fac
    !! The code
    do isph = 1, params % nsph
        do l = 0, params % lmax
            ind = l*l + l + 1
            fac = fourpi / (two*dble(l) + one) 
            if (direction .eq. -1) fac = one / fac
            do m = -l, l
                vector(ind+m, isph) = fac*vector(ind+m, isph)
            end do
        end do
    end do
end subroutine convert_ddcosmo

  !
  ! Intermediate computation of BX_e
  ! @param[in]      isph   : Number of the sphere
  ! @param[in, out] pot    : Array of size ngrid
  ! @param[in]      x      : Input vector (Usually X_e)
  ! @param[in, out] basloc : Used to compute spherical harmonic
  ! @param[in, out] vplm   : Used to compute spherical harmonic
  ! @param[in, out] vcos   : Used to compute spherical harmonic
  ! @param[in, out] vsin   : Used to compute spherical harmonic
  !
  subroutine calcv2_lpb (params, constants, isph, pot, x, basloc, vplm, vcos, &
          & vsin, bessel_work)
  type(ddx_params_type), intent(in)  :: params
  type(ddx_constants_type), intent(in)  :: constants
  integer, intent(in) :: isph
  real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
  real(dp), dimension(params % ngrid), intent(inout) :: pot
  real(dp), dimension(constants % nbasis), intent(inout) :: basloc
  real(dp), dimension(constants % nbasis), intent(inout) :: vplm
  real(dp), dimension(params % lmax+1), intent(inout) :: vcos
  real(dp), dimension(params % lmax+1), intent(inout) :: vsin
  complex(dp), intent(out) :: bessel_work(max(2, params % lmax+1))
  complex(dp) :: work_complex(params % lmax+1)
  real(dp) :: work(params % lmax+1)
  real(dp), dimension(constants % nbasis) :: fac_cosmo, fac_hsp
  real(dp), dimension(params % ngrid) :: pot2
  integer :: its, ij, jsph
  real(dp) :: rho, ctheta, stheta, cphi, sphi
  real(dp) :: vij(3), sij(3), vtij(3)
  real(dp) :: vvij, tij, xij, oij

  pot = zero
  pot2 = zero
  do its = 1, params % ngrid
    if (constants % ui(its,isph).lt.one) then
      do ij = constants % inl(isph), constants % inl(isph+1)-1
        jsph = constants % nl(ij)

        ! compute geometrical variables
        vij  = params % csph(:,isph) + params % rsph(isph)*constants % cgrid(:,its) - params % csph(:,jsph)
        vvij = sqrt(dot_product(vij,vij))
        tij  = vvij/params % rsph(jsph)

        if ( tij.lt.( one + (params % se+one)/two*params % eta ) ) then
          sij = vij/vvij
          xij = fsw(tij, params % se, params % eta)
          if (constants % fi(its,isph).gt.one) then
            oij = xij/constants % fi(its, isph)
          else
            oij = xij
          end if
          vtij = vij*params % kappa
          call fmm_l2p_bessel_work(vtij, &
              & params % lmax, constants % vscales, &
              & constants % SI_ri(:, jsph), oij, x(:, jsph), one, &
              & pot(its), work_complex, work)
        end if
      end do
    end if
  end do
  endsubroutine calcv2_lpb


  !
  ! Intermediate calculation in calcv2_lpb subroutine
  ! @param[in]  rijn    : Radius of sphers x_ijn
  ! @param[in]  ri      : Radius of sphers x_i
  ! @param[in]  isph    : Index of sphere
  ! @param[in]  basloc  : Spherical Harmonic
  ! @param[out] fac_hsp : Return bessel function ratio multiplied by 
  !                       the spherical harmonic Y_l'm'. Array of size nylm
  !
  subroutine inthsp(params, constants, rijn, ri, isph, basloc, fac_hsp, work)
  implicit none
  type(ddx_params_type), intent(in)  :: params
  type(ddx_constants_type), intent(in)  :: constants
  integer, intent(in) :: isph
  real(dp), intent(in) :: rijn, ri
  real(dp), dimension(constants % nbasis), intent(in) :: basloc
  real(dp), dimension(constants % nbasis), intent(inout) :: fac_hsp
  complex(dp), intent(out) :: work(max(2, params % lmax+1))
  real(dp), dimension(0:params % lmax) :: SI_rijn, DI_rijn
  integer :: l, m, ind

  SI_rijn = 0
  DI_rijn = 0
  fac_hsp = 0

  ! Computation of modified spherical Bessel function values
  call modified_spherical_bessel_first_kind(params % lmax, &
      & rijn*params % kappa, SI_rijn, DI_rijn, work)

  do l = 0, params % lmax
    do  m = -l, l
      ind = l*l + l + 1 + m
      fac_hsp(ind) = SI_rijn(l)/constants % SI_ri(l,isph)*basloc(ind)
    end do
  end do
  endsubroutine inthsp

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
    write(6,*) '@direct@intrhs', tt1 - tt0

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

! This routine is not needed when using a Jacobi/DIIS solver for the
! LPB linear systems. It is needed if using a GMRES solver.
subroutine lpb_adjoint_matvec_full(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y
    real(dp), dimension(constants % nbasis, params % nsph, 2) :: scratch

    call lpb_adjoint_matvec(params, constants, workspace, x, y)
    call lstarx(params, constants, workspace, x(:,:,1), scratch(:,:,1))
    call convert_ddcosmo(params, constants, -1, scratch(:,:,1))
    call bstarx(params, constants, workspace, x(:,:,2), scratch(:,:,2))
    y = y + scratch

end subroutine lpb_adjoint_matvec_full

! This routine is not needed when using a Jacobi/DIIS solver for the
! LPB linear systems. It is needed if using a GMRES solver.
subroutine lpb_direct_matvec_full(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y
    real(dp), dimension(constants % nbasis, params % nsph, 2) :: scratch

    call lpb_direct_matvec(params, constants, workspace, x, y)
    scratch(:,:,2) = x(:,:,1)
    call convert_ddcosmo(params, constants, -1, scratch(:,:,2))
    call lx(params, constants, workspace, scratch(:,:,2), scratch(:,:,1))
    call bx(params, constants, workspace, x(:,:,2), scratch(:,:,2))
    y = y + scratch

end subroutine lpb_direct_matvec_full

subroutine lpb_adjoint_matvec(params, constants, workspace, x, y)
    implicit none
    ! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y
    ! local
    real(dp), dimension(0:constants % lmax0) :: SK_rijn, DK_rijn
    real(dp), dimension(constants % nbasis) :: basloc, vplm
    real(dp), dimension(params % lmax + 1) :: vcos, vsin
    complex(dp) :: bessel_work(max(2, params % lmax+1))
    complex(dp) :: work_complex(constants % lmax0+1)
    real(dp) :: work(constants % lmax0+1)
    integer :: isph, igrid, jsph, ind, l0, m0, ind0, indl, inode, l, m
    real(dp), dimension(3) :: vij, vtij
    real(dp) :: val, epsilon_ratio
    real(dp), allocatable :: scratch(:,:), scratch0(:,:)

    allocate(scratch(constants % nbasis, params % nsph), &
        & scratch0(constants % nbasis0, params % nsph))

    epsilon_ratio = params % epsp/params % eps

    tt0 = omp_get_wtime()
    ! TODO: maybe use ddeval_grid for code consistency
    scratch = - x(:,:,1) - x(:,:,2)
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & one, constants % vwgrid, constants % vgrid_nbasis, scratch, &
        & constants % nbasis, zero, workspace % tmp_grid, params % ngrid)
    tt1 = omp_get_wtime()
    write(6,*) '@adjoint@matvec1', tt1 - tt0

    tt0 = omp_get_wtime()
    scratch0 = zero
    if (params % fmm .eq. 0) then
        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, isph).gt.zero) then
                    val = workspace % tmp_grid(igrid,isph) &
                        & *constants % ui(igrid,isph)
                    ! quadratically scaling loop
                    do jsph = 1, params % nsph
                        vij  = params % csph(:,isph) + &
                            & params % rsph(isph)*constants % cgrid(:,igrid) - &
                            & params % csph(:,jsph)
                        vtij = vij * params % kappa
                        call fmm_m2p_bessel_adj_work(vtij, val, &
                            & constants % SK_ri(:, jsph), constants % lmax0, &
                            & constants % vscales, one, scratch0(:, jsph), &
                            & work_complex, work)
                    end do
                end if
            end do
        end do
    else
        ! Multiply by characteristic function
        workspace % tmp_grid = workspace % tmp_grid * constants % ui
        workspace % tmp_sph = zero
        ! Do FMM operations adjointly
        call tree_m2p_bessel_adj(params, constants, constants % lmax0, &
            & one, workspace % tmp_grid, zero, params % lmax, workspace % tmp_sph)
        call tree_l2p_bessel_adj(params, constants, one, &
            & workspace % tmp_grid, zero, workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l, workspace % tmp_node_m)
        call tree_m2m_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_m)
        ! Adjointly move tree multipole harmonics into output
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                scratch0(:, isph) = workspace % tmp_sph(1:constants % nbasis0, isph) + &
                    & workspace % tmp_node_m(1:constants % nbasis0, inode)
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                scratch0(1:indl, isph) = &
                    & workspace % tmp_sph(1:indl, isph) + &
                    & workspace % tmp_node_m(:, inode)
                scratch0(indl+1:, isph) = zero
            end do
        end if
    end if
    ! Scale by C_ik
    do isph = 1, params % nsph
        do l0 = 0, constants % lmax0
            ind0 = l0*l0 + l0 + 1
            scratch0(ind0-l0:ind0+l0, isph) = scratch0(ind0-l0:ind0+l0, isph) * &
                & constants % C_ik(l0, isph)
        end do
    end do
    tt1 = omp_get_wtime()
    write(6,*) '@adjoint@matvec2', tt1 - tt0

    tt0 = omp_get_wtime()
    do jsph = 1, params % nsph
        call dgemv('n', constants % nbasis, constants % nbasis0, one, &
            & constants % pchi(1,1,jsph), constants % nbasis, &
            & scratch0(1,jsph), 1, zero, scratch(1,jsph), 1)
    end do
    tt1 = omp_get_wtime()
    write(6,*) '@adjoint@matvec3', tt1 - tt0

    tt0 = omp_get_wtime()
    do isph = 1, params % nsph
        do l = 0, params % lmax
            do m = -l, l
                ind = l**2 + l + m + 1
                y(ind,isph,1) = - (epsilon_ratio*l*scratch(ind,isph)) &
                    & /params % rsph(isph)
                y(ind,isph,2) = constants % termimat(l,isph)*scratch(ind,isph)
          end do
        end do
    end do
    tt1 = omp_get_wtime()
    write(6,*) '@adjoint@matvec4', tt1 - tt0

end subroutine lpb_adjoint_matvec

! Perform |Yr| = |C1 C2|*|Xr|
!         |Ye|   |C1 C2| |Xe|
! @param[in] ddx_data : dd Data
! @param[in] x        : Input array X (Xr, Xe)
! @param[out] y       : Matrix-vector product result Y
subroutine lpb_direct_matvec(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y

    integer :: isph, jsph, igrid, ind, l, m, ind0
    real(dp), dimension(3) :: sijn ,vij, vtij
    real(dp) :: term, rho, ctheta, stheta, cphi, sphi, rijn, val
    real(dp), dimension(constants % nbasis) :: basloc, vplm
    real(dp), dimension(params % lmax + 1) :: vcos, vsin
    real(dp), dimension(constants % lmax0 + 1) :: SK_rijn, DK_rijn
    complex(dp) :: work_complex(constants % lmax0 + 1)
    real(dp) :: work(constants % lmax0 + 1)
    integer :: indl, inode

    real(dp), allocatable :: diff_re(:,:), diff0(:,:)
    allocate(diff_re(constants % nbasis, params % nsph), &
        & diff0(constants % nbasis0, params % nsph))

    ! diff_re = params % epsp/eps*l1/ri*Xr - i'(ri)/i(ri)*Xe,
    tt0 = omp_get_wtime()
    diff_re = zero
    !$omp parallel do default(none) shared(params,diff_re, &
    !$omp constants,x) private(jsph,l,m,ind)
    do jsph = 1, params % nsph
      do l = 0, params % lmax
        do m = -l,l
          ind = l**2 + l + m + 1
          diff_re(ind,jsph) = (params % epsp/params % eps)* &
              & (l/params % rsph(jsph))*x(ind,jsph,1) &
              & - constants % termimat(l,jsph)*x(ind,jsph,2)
        end do
      end do
    end do
    tt1 = omp_get_wtime()
    write(6,*) '@direct@matvec1', tt1 - tt0

    ! diff0 = Pchi * diff_er, linear scaling
    !$omp parallel do default(none) shared(constants,params, &
    !$omp diff_re,diff0) private(jsph)
    do jsph = 1, params % nsph
      call dgemv('t', constants % nbasis, constants % nbasis0, one, &
          & constants % pchi(1,1,jsph), constants % nbasis, &
          & diff_re(1,jsph), 1, zero, diff0(1,jsph), 1)
    end do
    tt1 = omp_get_wtime()
    write(6,*) '@direct@matvec2', tt1 - tt0

    ! Multiply diff0 by C_ik inplace
    do isph = 1, params % nsph
        do l = 0, constants % lmax0
            ind0 = l*l+l+1
            diff0(ind0-l:ind0+l, isph) = diff0(ind0-l:ind0+l, isph) * &
                & constants % C_ik(l, isph)
        end do
    end do
    ! avoiding N^2 storage, this code does not use the cached coefY
    tt0 = omp_get_wtime()
    y(:,:,1) = zero
    if (params % fmm .eq. 0) then
        !$omp parallel do default(none) shared(params,constants, &
        !$omp diff0,y) private(isph,igrid,val,vij,rijn,sijn,SK_rijn, &
        !$omp DK_rijn,work,rho,ctheta,stheta,cphi,sphi,basloc,vplm, &
        !$omp vcos,vsin,term,ind0,ind,vtij,work_complex)
        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                if (constants % ui(igrid,isph).gt.zero) then
                    val = zero

                    ! quadratically scaling loop
                    do jsph = 1, params % nsph
                        vij  = params % csph(:,isph) + &
                            & params % rsph(isph)*constants % cgrid(:,igrid) - &
                            & params % csph(:,jsph)
                        vtij = vij * params % kappa
                        call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & one, diff0(:, jsph), one, val, work_complex, work)
                    end do
                    do ind = 1, constants % nbasis
                        y(ind,isph,1) = y(ind,isph,1) + val*&
                          & constants % ui(igrid,isph)*&
                          & constants % vwgrid(ind,igrid)
                    end do
                end if
            end do
        end do
    else
        ! Load input harmonics into tree data
        workspace % tmp_node_m = zero
        workspace % tmp_node_l = zero
        workspace % tmp_sph = zero
        do isph = 1, params % nsph
            do l = 0, constants % lmax0
                ind0 = l*l+l+1
                workspace % tmp_sph(ind0-l:ind0+l, isph) = &
                    & diff0(ind0-l:ind0+l, isph)
            end do
        end do
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(1:constants % nbasis0, inode) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph)
                workspace % tmp_node_m(constants % nbasis0+1:, inode) = zero
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(:, inode) = workspace % tmp_sph(1:indl, isph)
            end do
        end if
        ! Do FMM operations
        call tree_m2m_bessel_rotation(params, constants, workspace % tmp_node_m)
        call tree_m2l_bessel_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation(params, constants, workspace % tmp_node_l)
        workspace % tmp_grid = zero
        call tree_l2p_bessel(params, constants, one, workspace % tmp_node_l, zero, &
            & workspace % tmp_grid)
        call tree_m2p_bessel(params, constants, constants % lmax0, one, &
            & params % lmax, workspace % tmp_sph, one, &
            & workspace % tmp_grid)

        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                do ind = 1, constants % nbasis
                    y(ind,isph,1) = y(ind,isph,1) + workspace % tmp_grid(igrid, isph)*&
                        & constants % vwgrid(ind, igrid)*&
                        & constants % ui(igrid,isph)
                end do
            end do
        end do
    end if
    tt1 = omp_get_wtime()
    write(6,*) '@direct@matvec3', tt1 - tt0

    y(:,:,2) = y(:,:,1)
    deallocate(diff_re, diff0)

end subroutine lpb_direct_matvec

subroutine lpb_adjoint_prec(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph, 2)
    real(dp), intent(inout) :: y(constants % nbasis, params % nsph, 2)
    integer :: n_iter, info
    real(dp) :: r_norm
    real(dp), dimension(params % maxiter) :: x_rel_diff

    tt0 = omp_get_wtime()
    y(:,:,1) = x(:,:,1)
    call convert_ddcosmo(params, constants, 1, y(:,:,1))
    n_iter = params % maxiter
    if (params % incore) then
        call jacobi_diis(params, constants, workspace, constants % inner_tol, y(:,:,1), &
            & workspace % ddcosmo_guess, n_iter, x_rel_diff, lstarx_nodiag_incore, ldm1x, hnorm, info)
    else
        call jacobi_diis(params, constants, workspace, constants % inner_tol, y(:,:,1), &
            & workspace % ddcosmo_guess, n_iter, x_rel_diff, lstarx_nodiag, ldm1x, hnorm, info)
    end if
    if (info.ne.0) then
        write(*,*) 'lpb_adjoint_prec: [1] ddCOSMO failed to converge'
        stop 1
    end if
    y(:,:,1) = workspace % ddcosmo_guess
    tt1 = omp_get_wtime()
    write(6,*) '@adjoint@cosmo', tt1 - tt0, n_iter

    tt0 = omp_get_wtime()
    n_iter = params % maxiter
    if (params % incore) then
        ! call gmresr(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
        !     & n_iter, r_norm, bstarx_incore, info)
        call jacobi_diis(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
            & n_iter, x_rel_diff, bstarx_nodiag_incore, bx_prec, hnorm, info)
    else
        ! call gmresr(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
        !    & n_iter, r_norm, bstarx, info)
        call jacobi_diis(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
            & n_iter, x_rel_diff, bstarx_nodiag, bx_prec, hnorm, info)
    end if
    if (info.ne.0) then
        write(*,*) 'lpb_adjoint_prec: [1] HSP failed to converge'
        stop 1
    end if
    y(:,:,2) = workspace % hsp_guess
    tt1 = omp_get_wtime()
    write(6,*) '@adjoint@hsp', tt1 - tt0, n_iter

end subroutine lpb_adjoint_prec

! apply preconditioner
! |Yr| = |A^-1 0 |*|Xr|
! |Ye|   |0 B^-1 | |Xe|
! @param[in] ddx_data : dd Data
! @param[in] x        : Input array
! @param[out] y       : Linear system solution at current iteration
subroutine lpb_direct_prec(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph, 2)
    real(dp), intent(inout) :: y(constants % nbasis, params % nsph, 2)
    integer :: n_iter, info
    real(dp) :: r_norm
    real(dp), dimension(params % maxiter) :: x_rel_diff

    ! perform A^-1 * Yr
    tt0 = omp_get_wtime()
    n_iter = params % maxiter
    if (params % incore) then
        call jacobi_diis(params, constants, workspace, constants % inner_tol, x(:,:,1), &
            & workspace % ddcosmo_guess, n_iter, x_rel_diff, lx_nodiag_incore, ldm1x, hnorm, info)
    else
        call jacobi_diis(params, constants, workspace, constants % inner_tol, x(:,:,1), &
            & workspace % ddcosmo_guess, n_iter, x_rel_diff, lx_nodiag, ldm1x, hnorm, info)
    end if
    if (info.ne.0) then
        write(*,*) 'lpb_direct_prec: [1] ddCOSMO failed to converge'
        stop 1
    end if

    ! Scale by the factor of (2l+1)/4Pi
    y(:,:,1) = workspace % ddcosmo_guess
    call convert_ddcosmo(params, constants, 1, y(:,:,1))
    tt1 = omp_get_wtime()
    write(6,*) '@direct@cosmo', tt1 - tt0, n_iter

    ! perform B^-1 * Ye
    tt0 = omp_get_wtime()
    n_iter = params % maxiter
    if (params % incore) then
        ! call gmresr(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
        !    & n_iter, r_norm, bx_incore, info)
        call jacobi_diis(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
            & n_iter, x_rel_diff, bx_nodiag_incore, bx_prec, hnorm, info)
    else
        ! call gmresr(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
        !  & n_iter, r_norm, bx, info)
        call jacobi_diis(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
            & n_iter, x_rel_diff, bx_nodiag, bx_prec, hnorm, info)
    end if
    y(:,:,2) = workspace % hsp_guess
    tt1 = omp_get_wtime()
    write(6,*) '@direct@hsp', tt1 - tt0, n_iter

    if (info.ne.0) then
        write(*,*) 'lpb_direct_prec: [1] HSP failed to converge'
        stop 1
    end if
end subroutine lpb_direct_prec


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

subroutine bstarx(params, constants, workspace, x, y)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    integer :: isph, igrid, iproc
    y = zero
    !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
    !$omp private(isph) schedule(static,1)
    do isph = 1, params % nsph
        call dgemv('t', constants % nbasis, params % ngrid, one, constants % vgrid, &
            & constants % vgrid_nbasis, x(:, isph), 1, zero, &
            & workspace % tmp_grid(:, isph), 1)
    end do
    !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
    !$omp private(isph,iproc) schedule(static,1)
    do isph = 1, params % nsph
        iproc = omp_get_thread_num() + 1
        call adjrhs_lpb(params, constants, workspace, isph, workspace % tmp_grid, &
            & y(:, isph), workspace % tmp_vylm(:, iproc), workspace % tmp_vplm(:, iproc), &
            & workspace % tmp_vcos(:, iproc), workspace % tmp_vsin(:, iproc), &
            & workspace % tmp_bessel(:, iproc))
        y(:,isph)  = - y(:,isph) + x(:,isph)
    end do
end subroutine bstarx

subroutine bstarx_nodiag(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    integer :: isph, igrid, iproc
    y = zero
    !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
    !$omp private(isph) schedule(static,1)
    do isph = 1, params % nsph
        call dgemv('t', constants % nbasis, params % ngrid, one, constants % vgrid, &
            & constants % vgrid_nbasis, x(:, isph), 1, zero, &
            & workspace % tmp_grid(:, isph), 1)
    end do
    !$omp parallel do default(none) shared(params,constants,workspace,y) &
    !$omp private(isph,iproc) schedule(static,1)
    do isph = 1, params % nsph
        iproc = omp_get_thread_num() + 1
        call adjrhs_lpb(params, constants, workspace, isph, workspace % tmp_grid, &
            & y(:, isph), workspace % tmp_vylm(:, iproc), workspace % tmp_vplm(:, iproc), &
            & workspace % tmp_vcos(:, iproc), workspace % tmp_vsin(:, iproc), &
            & workspace % tmp_bessel(:, iproc))
        y(:,isph) = - y(:,isph)
    end do
end subroutine bstarx_nodiag

!
! Taken from ddx_core routine adjrhs
! Called from bstarx
! Compute the Adjoint matix B*x
!
subroutine adjrhs_lpb(params, constants, workspace, isph, xi, vlm, basloc, &
    & vplm, vcos, vsin, tmp_bessel)
    implicit none
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: xi
    real(dp), dimension((params % lmax+1)**2), intent(out) :: vlm
    real(dp), dimension((params % lmax+1)**2), intent(out) :: basloc, vplm
    real(dp), dimension(params % lmax+1), intent(out) :: vcos, vsin
    complex(dp), dimension(max(2, params % lmax+1)), intent(out) :: tmp_bessel
    type(ddx_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: isph
    !! local
    integer :: ij, jsph, ig, l, ind, m
    real(dp) :: vji(3), vvji, tji, sji(3), xji, oji, fac
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    real(dp), dimension(constants % nbasis) :: fac_hsp

    !loop over neighbors of i-sphere
    do ij = constants % inl(isph), constants % inl(isph+1)-1
      !j-sphere is neighbor
      jsph = constants % nl(ij)
      !loop over integration points
      do ig = 1, params % ngrid
        !compute t_n^ji = | r_j + \rho_j s_n - r_i | / \rho_i
        vji  = params % csph(:,jsph) + params % rsph(jsph)* &
              & constants % cgrid(:,ig) - params % csph(:,isph)
        vvji = sqrt(dot_product(vji,vji))
        tji  = vvji/params % rsph(isph)
        !point is INSIDE i-sphere (+ transition layer)
        if ( tji.lt.( one + (params % se+one)/two*params % eta ) ) then
          !compute s_n^ji
          sji = vji/vvji
          call ylmbas(sji, rho, ctheta, stheta, cphi, &
                        & sphi, params % lmax, &
                        & constants % vscales, basloc, &
                        & vplm, vcos, vsin)
          call inthsp_adj(params, constants, vvji, params % rsph(isph), isph, &
              & basloc, fac_hsp, tmp_bessel)
          !compute \chi( t_n^ji )
          xji = fsw( tji, params % se, params % eta )
          !compute W_n^ji
          if ( constants % fi(ig,jsph).gt.one ) then
            oji = xji/ constants % fi(ig,jsph)
          else
            oji = xji
          endif
          !compute w_n * xi(n,j) * W_n^ji
          fac = constants % wgrid(ig) * xi(ig,jsph) * oji
          !loop over l
          do l = 0, params % lmax
            ind  = l*l + l + 1
            !loop over m
              do m = -l,l
                vlm(ind+m) = vlm(ind+m) + fac*fac_hsp(ind+m)
              enddo
          enddo
        endif
      enddo
    enddo
end subroutine adjrhs_lpb

!
! Intermediate calculation in adjrhs_lpb subroutine
! @param[in]  rijn    : Radius of sphers x_ijn
! @param[in]  ri      : Radius of sphers x_i
! @param[in]  isph    : Index of sphere
! @param[in]  basloc  : Spherical Harmonic
! @param[out] fac_hsp : Return bessel function ratio multiplied by 
!                       the spherical harmonic Y_l'm'. Array of size nylm
!
subroutine inthsp_adj(params, constants, rjin, rj, jsph, basloc, &
          & fac_hsp, work)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: jsph
    real(dp), intent(in) :: rjin, rj
    real(dp), dimension(constants % nbasis), intent(in) :: basloc
    real(dp), dimension(constants % nbasis), intent(inout) :: fac_hsp
    complex(dp), intent(out) :: work(max(2, params % lmax+1))
    real(dp), dimension(0:params % lmax) :: SI_rjin, DI_rjin
    integer :: l, m, ind

    SI_rjin = 0
    DI_rjin = 0
    fac_hsp = 0

    ! Computation of modified spherical Bessel function values
    call modified_spherical_bessel_first_kind(params % lmax, &
        & rjin*params % kappa, SI_rjin, DI_rjin, work)

    do l = 0, params % lmax
      do  m = -l, l
        ind = l*l + l + 1 + m
        fac_hsp(ind) = SI_rjin(l)/constants % SI_ri(l,jsph)*basloc(ind)
      end do
    end do
endsubroutine inthsp_adj

! Subroutine to compute K^A counterpart for the HSP equation. Similar to fdoka.
! @param[in]  ddx_data  : Data type
! @param[in]  isph      : Index of sphere
! @param[in]  Xe        : Solution vector Xe
! @param[in]  Xadj_e    : Adjoint solution on evaluated on grid points Xadj_e_sgrid
! @param[in]  basloc    : Spherical harmonics Y_lm
! @param[in]  dbasloc   : Derivative of spherical harmonics \nabla^i(Y_lm)
! @param[in]  vplm      : Argument to call ylmbas
! @param[in]  vcos      : Argument to call ylmbas
! @param[in]  vsin      : Argument to call ylmbas
! @param[out] force_e   : Force of adjoint part
subroutine fdoka_b_xe(params, constants, workspace, isph, Xe, Xadj_e, basloc, &
    & dbasloc, vplm, vcos, vsin, force_e)
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: isph
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xe
    real(dp), dimension(params % ngrid), intent(in) :: Xadj_e
    real(dp), dimension(constants % nbasis), intent(inout) :: basloc, vplm
    real(dp), dimension(3, constants % nbasis), intent(inout) :: dbasloc
    real(dp), dimension(params % lmax+1), intent(inout) :: vcos, vsin
    real(dp), dimension(3), intent(inout) :: force_e

    ! Local Variables
    integer :: igrid, ineigh, jsph, l, ind, m
    real(dp), dimension(0:params % lmax) :: SI_rijn
    real(dp), dimension(0:params % lmax) :: DI_rijn
    ! beta   : Eq.(53) Stamm.etal.18
    ! tlow   : Lower bound for switch region
    ! thigh  : Upper bound for switch region
    ! f1     : First factor in alpha computation
    ! f2     : Second factor in alpha computation
    real(dp)  :: rijn, tij, beta, tlow, thigh, xij, oij, f1, f2, f3
    ! alpha : Eq.(52) Stamm.etal.18
    ! va    : Eq.(54) Stamm.etal.18
    real(dp)  :: vij(3), sij(3), alpha(3), va(3), rj, vtij(3)
    real(dp), external :: dnrm2
    real(dp) :: work(params % lmax+1)
    complex(dp) :: work_complex(params % lmax+1)

    SI_rijn = 0
    DI_rijn = 0

    tlow  = one - pt5*(one - params % se)*params % eta
    thigh = one + pt5*(one + params % se)*params % eta

    ! Loop over grid points
    do igrid = 1, params % ngrid
        va = zero
        do ineigh = constants % inl(isph), constants % inl(isph+1) - 1
            jsph = constants % nl(ineigh)
            vij  = params % csph(:,isph) + &
                & params % rsph(isph)*constants % cgrid(:,igrid) - &
                & params % csph(:,jsph)
            rijn = dnrm2(3, vij, 1)
            tij  = rijn/params % rsph(jsph)
            rj = params % rsph(jsph)
            if (tij.ge.thigh) cycle
            ! Computation of modified spherical Bessel function values
            sij  = vij/rijn
            vtij = vij*params % kappa
            call fmm_l2p_bessel_grad(vtij, params % rsph(jsph)*params % kappa, &
                & params % lmax, constants % vscales, params % kappa, &
                & Xe(:, jsph), zero, alpha)
            call fmm_l2p_bessel_work(vtij, params % lmax, constants % vscales, &
                & constants % SI_ri(:, jsph), one, Xe(:, jsph), zero, beta, &
                & work_complex, work)
            xij = fsw(tij, params % se, params % eta)
            if (constants % fi(igrid,isph).gt.one) then
                oij = xij/constants % fi(igrid,isph)
                f2  = -oij/constants % fi(igrid,isph)
            else
                oij = xij
                f2  = zero
            end if
            f1 = oij
            va(:) = va(:) + f1*alpha(:) + beta*f2*constants % zi(:,igrid,isph)
            if (tij .gt. tlow) then
                f3 = beta*dfsw(tij,params % se,params % eta)/params % rsph(jsph)
                if (constants % fi(igrid,isph).gt.one) &
                    & f3 = f3/constants % fi(igrid,isph)
                va(:) = va(:) + f3*sij(:)
            end if
        end do
        force_e = force_e - constants % wgrid(igrid)*Xadj_e(igrid)*va(:)
    end do
end subroutine fdoka_b_xe

!
! Subroutine to compute K^A+K^C counterpart for the HSP equation. Similar to fdokb.
! @param[in]  ddx_data  : Data type
! @param[in]  isph      : Index of sphere
! @param[in]  Xe        : Solution vector Xe
! @param[in]  Xadj_e    : Adjoint solution on evaluated on grid points Xadj_e_sgrid
! @param[in]  basloc    : Spherical harmonics Y_lm
! @param[in]  dbasloc   : Derivative of spherical harmonics \nabla^i(Y_lm)
! @param[in]  vplm      : Argument to call ylmbas
! @param[in]  vcos      : Argument to call ylmbas
! @param[in]  vsin      : Argument to call ylmbas
! @param[out] force_e   : Force of adjoint part
subroutine fdokb_b_xe(params, constants, workspace, isph, Xe, Xadj_e, basloc, dbasloc, &
    & vplm, vcos, vsin, force_e)
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    integer, intent(in)    :: isph
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xe
    real(dp), dimension(params % ngrid, params % nsph),  intent(in) :: Xadj_e
    real(dp), dimension(constants % nbasis), intent(inout) :: basloc,  vplm
    real(dp), dimension(3, constants % nbasis), intent(inout) :: dbasloc
    real(dp), dimension(params % lmax+1), intent(inout) :: vcos, vsin
    real(dp), dimension(3), intent(inout) :: force_e

    ! Local Variables
    ! jk : Row pointer over kth row
    integer :: igrid, jsph, ksph, ineigh, l, m, ind, jk
    real(dp), dimension(0:params % lmax) :: SI_rjin, SI_rjkn
    real(dp), dimension(0:params % lmax) :: DI_rjin, DI_rjkn

    logical :: proc
    ! fac     : \delta_fj_n*\omega^\eta_ji
    ! f1      : First factor in alpha computation
    ! f2      : Second factor in alpha computation
    ! beta_ji : Eq.(57) Stamm.etal.18
    ! dj      : Before Eq.(10) Stamm.etal.18
    real(dp)  :: rjin, tji, xji, oji, fac, f1, f2, beta_ji, dj, tlow, thigh
    real(dp)  :: b, beta_jk, g1, g2, rjkn, tjk, xjk
    ! alpha : Eq.(56) Stamm.etal.18
    ! vb    : Eq.(60) Stamm.etal.18
    ! vc    : Eq.(59) Stamm.etal.18
    real(dp)  :: vji(3), sji(3), vjk(3), sjk(3), alpha(3), vb(3), vc(3), &
        & vtji(3), vtjk(3)
    ! rho    : Argument for ylmbas
    ! ctheta : Argument for ylmbas
    ! stheta : Argument for ylmbas
    ! cphi   : Argument for ylmbas
    ! sphi   : Argument for ylmbas
    real(dp) :: rho, ctheta, stheta, cphi, sphi, ri, arg_bessel

    real(dp), external :: dnrm2
    real(dp) :: work(params % lmax+1)
    complex(dp) :: work_complex(params % lmax+1)

    SI_rjin = 0
    DI_rjin = 0
    SI_rjkn = 0
    DI_rjkn = 0

    tlow  = one - pt5*(one - params % se)*params % eta
    thigh = one + pt5*(one + params % se)*params % eta

    do igrid = 1, params % ngrid
        vb = zero
        vc = zero
        do ineigh = constants % inl(isph), constants % inl(isph+1) - 1
            jsph = constants % nl(ineigh)
            vji  = params % csph(:,jsph) + &
                    & params % rsph(jsph)*constants % cgrid(:,igrid) - &
                    & params % csph(:,isph)
            rjin = dnrm2(3, vji, 1)
            ri = params % rsph(isph)
            tji  = rjin/ri

            if (tji.gt.thigh) cycle

            sji  = vji/rjin
            vtji = vji*params % kappa
            call fmm_l2p_bessel_grad(vtji, params % rsph(isph)*params % kappa, &
                & params % lmax, constants % vscales, params % kappa, &
                & Xe(:, isph), zero, alpha)
            xji = fsw(tji,params % se,params % eta)
            if (constants % fi(igrid,jsph).gt.one) then
                oji = xji/constants % fi(igrid,jsph)
            else
                oji = xji
            end if
            f1 = oji
            vb = vb + f1*alpha*Xadj_e(igrid,jsph)
            if (tji .gt. tlow) then
                call fmm_l2p_bessel_work(vtji, params % lmax, &
                    & constants % vscales, constants % SI_ri(:, isph), one, &
                    & Xe(:, isph), zero, beta_ji, work_complex, work)
            if (constants % fi(igrid,jsph) .gt. one) then
                dj  = one/constants % fi(igrid,jsph)
                fac = dj*xji
                proc = .false.
                b    = zero
                do jk = constants % inl(jsph), constants % inl(jsph+1) - 1
                    ksph = constants % nl(jk)
                    vjk  = params % csph(:,jsph) + &
                        & params % rsph(jsph)*constants % cgrid(:,igrid) - &
                        & params % csph(:,ksph)
                    rjkn = dnrm2(3, vjk, 1)
                    tjk  = rjkn/params % rsph(ksph)
                    if (ksph.ne.isph) then
                        if (tjk .le. thigh) then
                            proc = .true.
                            vtjk = vjk*params % kappa
                            call fmm_l2p_bessel_work(vtjk, params % lmax, &
                                & constants % vscales, &
                                & constants % SI_ri(:, ksph), one, Xe(:, ksph), &
                                & zero, beta_jk, work_complex, work)
                            xjk = fsw(tjk, params % se, params % eta)
                            b = b + beta_jk*xjk
                        end if
                    end if
                end do
                if (proc) then
                    g1 = dj*dj*dfsw(tji,params % se,params % eta) &
                        & /params % rsph(isph)
                    g2 = g1*Xadj_e(igrid,jsph)*b
                    vc = vc + g2*sji
                end if
              else
                  dj  = one
                  fac = zero
              end if
              f2 = (one-fac)*dj*dfsw(tji,params % se,params % eta) &
                  & /params % rsph(isph)
              vb = vb + f2*Xadj_e(igrid,jsph)*beta_ji*sji
            end if
        end do
        force_e = force_e + constants % wgrid(igrid)*(vb - vc)
    end do
end subroutine fdokb_b_xe

real(dp) function compute_beta(params, constants, workspace, SI_rijn, rijn, jsph, Xe, basloc)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    !! Temporary buffers
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(0:params % lmax), intent(in) :: SI_rijn
    real(dp), dimension(constants % nbasis), intent(in) :: basloc
    real(dp), dimension(constants % nbasis), intent(in)   :: Xe
    integer, intent(in) :: jsph
    real(dp), intent(in) :: rijn

    integer :: l, m, ind
    real(dp)  :: ss, fac
    ss = zero

    ! loop over l
    do l = 0, params % lmax
      do m = -l, l
        ind = l*l + l + m + 1
        fac = SI_rijn(l)/constants % SI_ri(l,jsph)
        ss = ss + fac*basloc(ind)*Xe(ind)
      end do
    end do

    compute_beta = ss
end function compute_beta

!
! Subroutine to compute the derivative of U_i(x_in) and \bf{k}^j_l0(x_in)Y^j_l0m0(x_in)
! fdouky : Force Derivative of U_i^e(x_in), k_l0, and Y_l0m0
! @param[in] ddx_data     : Data Type
! @param[in] ksph         : Derivative with respect to x_k
! @param[in] Xr           : Solution of the Laplace problem
! @param[in] Xe           : Solution of the HSP problem
! @param[in] Xadj_r_sgrid : Solution of the Adjoint Laplace problem evaluated at the
!                           grid
! @param[in] Xadj_e_sgrid : Solution of the Adjoint HSP problem evaluated at the grid
! @param[in] Xadj_r       : Adjoint solution of the Laplace problem
! @param[in] Xadj_e       : Adjoint solution of the HSP problem
! @param[inout] force     : Force
! @param[out] diff_re     : epsilon_1/epsilon_2 * l'/r_j[Xr]_jl'm'
!                         - (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
subroutine fdouky(params, constants, workspace, Xr, Xe, Xadj_r_sgrid, &
    & Xadj_e_sgrid, Xadj_r, Xadj_e, force, diff_re)
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xr, Xe
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: &
        & Xadj_r_sgrid, Xadj_e_sgrid
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: &
        & Xadj_r, Xadj_e
    real(dp), dimension(3, params % nsph), intent(inout) :: force
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: &
        & diff_re
    real(dp), external :: dnrm2
    ! Local variable
    integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, icav, indl, inode, &
        & ksph, knode, jnode, knear, jsph_node, istat
    ! val_dim3 : Intermediate value array of dimension 3
    real(dp), dimension(3) :: sij, vij, val_dim3, vtij
    ! val   : Intermediate variable to compute diff_ep
    ! f1    : Intermediate variable for derivative of coefY_der
    ! f2    : Intermediate variable for derivative of coefY_der
    real(dp) :: val, f1, f2
    ! large local are allocatable
    ! phi_in : sum_{j=1}^N diff0_j * coefY_j
    ! diff_ep_dim3 : 3 dimensional couterpart of diff_ep
    ! sum_dim3 : Storage of sum
    ! diff0       : dot_product([PU_j]_l0m0^l'm', l'/r_j[Xr]_jl'm' -
    !                        (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm')
    real(dp), allocatable :: diff_ep_dim3(:,:), phi_in(:,:), sum_dim3(:,:,:), &
        & diff0(:,:), diff1(:,:), diff1_grad(:,:,:), l2l_grad(:,:,:)
    real(dp) :: termi, termk, rijn
    ! basloc : Y_lm(s_n)
    ! vplm   : Argument to call ylmbas
    real(dp),  dimension(constants % nbasis):: basloc, vplm
    ! dbasloc : Derivative of Y_lm(s_n)
    real(dp),  dimension(3, constants % nbasis):: dbasloc
    ! vcos   : Argument to call ylmbas
    ! vsin   : Argument to call ylmbas
    real(dp),  dimension(params % lmax+1):: vcos, vsin
    real(dp), dimension(0:params % lmax) :: SK_rijn, DK_rijn
    complex(dp) :: work_complex(constants % lmax0+2)
    real(dp) :: work(constants % lmax0+2)

    allocate(diff_ep_dim3(3, constants % ncav), &
        & phi_in(params % ngrid, params % nsph), &
        & sum_dim3(3, constants % nbasis, params % nsph), &
        & diff0(constants % nbasis0, params % nsph), &
        & diff1(constants % nbasis0, params % nsph), &
        & diff1_grad((constants % lmax0+2)**2, 3, params % nsph), &
        & l2l_grad((params % pl+2)**2, 3, params % nsph), stat=istat)
    if (istat.ne.0) stop 1

    ! Setting initial values to zero
    SK_rijn = zero
    DK_rijn = zero

    diff_re = zero
    ! Compute l'/r_j[Xr]_jl'm' -(i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
    do jsph = 1, params % nsph
      do l = 0, params % lmax
        do m = -l,l
          ind = l**2 + l + m + 1
          diff_re(ind,jsph) = (params % epsp/params % eps)*(l/params % rsph(jsph)) * &
                & Xr(ind,jsph) - constants % termimat(l,jsph)*Xe(ind,jsph)
        end do
      end do
    end do

    ! diff0 = Pchi * diff_re, linear scaling
    diff0 = zero
    do jsph = 1, params % nsph
      do l0 = 0, constants % lmax0
        do ind0 = l0*l0+1, l0*l0+2*l0+1
          diff0(ind0, jsph) = dot_product(diff_re(:,jsph), &
              & constants % Pchi(:,ind0, jsph))
          diff1(ind0, jsph) = diff0(ind0, jsph) * constants % C_ik(l0, jsph)
        end do
      end do
      ! Prepare diff1_grad
      call fmm_m2m_bessel_grad(constants % lmax0, constants % SK_ri(:, jsph), &
          & constants % vscales, constants % vcnk, diff1(:, jsph), &
          & diff1_grad(:, :, jsph))
    end do

    if (params % fmm .eq. 0) then
        ! phi_in = diff0 * coefY
        ! Here, summation over j takes place
        phi_in = zero
        icav = 0
        do isph = 1, params % nsph
          do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .gt. zero) then
              ! Extrenal grid point
              icav = icav + 1
              val = zero
              do jsph = 1, params % nsph
                !do ind0 = 1, constants % nbasis0
                !!====== This place requirs coefY, that is not precomputed anymore
                !  val = val + diff0(ind0,jsph)*constants % coefY(icav,ind0,jsph)
                !end do
                vij = params % csph(:, isph) + &
                    & params % rsph(isph)*constants % cgrid(:, igrid) - &
                    & params % csph(:, jsph)
                vtij = vij * params % kappa
                call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                    & constants % vscales, constants % SK_ri(:, jsph), one, &
                    & diff1(:, jsph), one, val, work_complex, work)
              end do
              phi_in(igrid, isph) = val
            end if
          end do
        end do
    else
        ! phi_in
        ! Load input harmonics into tree data
        workspace % tmp_sph = zero
        workspace % tmp_sph(1:constants % nbasis0, :) = diff1(:, :)
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(1:constants % nbasis0, inode) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph)
                workspace % tmp_node_m(constants % nbasis0+1:, inode) = zero
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(:, inode) = workspace % tmp_sph(1:indl, isph)
            end do
        end if
        ! Do FMM operations
        call tree_m2m_bessel_rotation(params, constants, workspace % tmp_node_m)
        call tree_m2l_bessel_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation(params, constants, workspace % tmp_node_l)
        call tree_l2p_bessel(params, constants, one, workspace % tmp_node_l, zero, &
            & phi_in)
        call tree_m2p_bessel(params, constants, constants % lmax0, one, &
            & params % lmax, workspace % tmp_sph, one, &
            & phi_in)
        ! Make phi_in zero at internal grid points
        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, isph) .eq. zero) then
                    phi_in(igrid, isph) = zero
                end if
            end do
        end do
        ! Get gradients of the L2L
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_sph_l(:, isph) = workspace % tmp_node_l(:, inode)
            call fmm_l2l_bessel_grad(params % pl, &
                & constants % SI_ri(:, isph), constants % vscales, &
                & constants % vcnk, workspace % tmp_node_l(:, inode), &
                & l2l_grad(:, :, isph))
        end do
        workspace % tmp_sph = Xadj_r + Xadj_e
        call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
            & one, constants % vwgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, &
            & workspace % tmp_grid, params % ngrid)
        workspace % tmp_grid = workspace % tmp_grid * constants % ui
        ! Adjoint FMM with output tmp_sph2(:, :) which stores coefficients of
        ! harmonics of degree up to lmax+1
        call tree_m2p_bessel_nodiag_adj(params, constants, constants % lmax0+1, one, &
            & workspace % tmp_grid, zero, params % lmax+1, workspace % tmp_sph2)
        call tree_l2p_bessel_adj(params, constants, one, workspace % tmp_grid, zero, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_node_m)
        call tree_m2m_bessel_rotation_adj(params, constants, workspace % tmp_node_m)
        ! Properly load adjoint multipole harmonics into tmp_sph2 that holds
        ! harmonics of a degree up to lmax+1
        if(constants % lmax0+1 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph2(1:(constants % lmax0+2)**2, isph) = &
                    & workspace % tmp_sph2(1:(constants % lmax0+2)**2, isph) + &
                    & workspace % tmp_node_m(1:(constants % lmax0+2)**2, inode)
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph2(1:indl, isph) = &
                    & workspace % tmp_sph2(1:indl, isph) + &
                    & workspace % tmp_node_m(:, inode)
            end do
        end if
    end if

    do ksph = 1, params % nsph
        ! Computation of derivative of U_i^e(x_in)
        call fdoga(params, constants, ksph, Xadj_r_sgrid, phi_in, force(:, ksph))
        call fdoga(params, constants, ksph, Xadj_e_sgrid, phi_in, force(:, ksph))

        ! Aleksandr: my loop for the diff_ep_dim3
        diff_ep_dim3 = zero
        ! At first isph=ksph, jsph!=ksph
        icav = constants % icav_ia(ksph) - 1
        if (params % fmm .eq. 0) then
            do igrid = 1, params % ngrid
              if (constants % ui(igrid, ksph) .eq. zero) cycle
              icav = icav + 1
              do jsph = 1, params % nsph
                if (jsph .eq. ksph) cycle
                  vij  = params % csph(:,ksph) + &
                      & params % rsph(ksph)*constants % cgrid(:,igrid) - &
                      & params % csph(:,jsph)
                  vtij = vij * params % kappa
                  !call fmm_m2p_bessel_grad(vtij, &
                  !    & params % rsph(jsph)*params % kappa, &
                  !    & constants % lmax0, &
                  !    & constants % vscales, params % kappa, diff1(:, jsph), one, &
                  !    & diff_ep_dim3(:, icav))
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, diff1_grad(:, 1, jsph), one, &
                      & diff_ep_dim3(1, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, diff1_grad(:, 2, jsph), one, &
                      & diff_ep_dim3(2, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, diff1_grad(:, 3, jsph), one, &
                      & diff_ep_dim3(3, icav), work_complex, work)
              end do
            end do
        else
            knode = constants % snode(ksph)
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, ksph) .eq. zero) cycle
                icav = icav + 1
                ! Far-field
                call dgemv('T', (params % pl+2)**2, 3, -params % kappa, &
                    & l2l_grad(1, 1, ksph), &
                    & (params % pl+2)**2, constants % vgrid(1, igrid), 1, &
                    & one, diff_ep_dim3(1, icav), 1)
                ! Near-field
                do knear = constants % snear(knode), constants % snear(knode+1)-1
                    jnode = constants % near(knear)
                    do jsph_node = constants % cluster(1, jnode), &
                        & constants % cluster(2, jnode)
                        jsph = constants % order(jsph_node)
                        if (jsph .eq. ksph) cycle
                        vij = params % csph(:, ksph) - params % csph(:, jsph) + &
                            & params % rsph(ksph)*constants % cgrid(:, igrid)
                        vtij = vij * params % kappa
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, diff1_grad(:, 1, jsph), one, &
                            & diff_ep_dim3(1, icav), work_complex, work)
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, diff1_grad(:, 2, jsph), one, &
                            & diff_ep_dim3(2, icav), work_complex, work)
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, diff1_grad(:, 3, jsph), one, &
                            & diff_ep_dim3(3, icav), work_complex, work)
                    end do
                end do
            end do
        end if
        sum_dim3 = zero
        icav = constants % icav_ia(ksph) - 1
        do igrid =1, params % ngrid
            if(constants % ui(igrid, ksph) .gt. zero) then
                icav = icav + 1
                do ind = 1, constants % nbasis
                    sum_dim3(:,ind,ksph) = sum_dim3(:,ind,ksph) + &
                        & constants % coefvec(igrid, ind, ksph)*diff_ep_dim3(:,icav)
                end do
            end if
        end do
        do ind = 1, constants % nbasis
            force(:, ksph) = force(:, ksph) + &
                & sum_dim3(:, ind, ksph)*(Xadj_r(ind, ksph) + &
                & Xadj_e(ind, ksph))
        end do

        ! Now jsph=ksph and isph!=ksph
        if (params % fmm .eq. 0) then
            diff_ep_dim3 = zero
            sum_dim3 = zero
            do isph = 1, params % nsph
              if (isph .eq. ksph) cycle
              icav = constants % icav_ia(isph) - 1
              do igrid = 1, params % ngrid
                  if (constants % ui(igrid, isph) .eq. zero) cycle
                  icav = icav + 1
                  vij  = params % csph(:,isph) + &
                      & params % rsph(isph)*constants % cgrid(:,igrid) - &
                      & params % csph(:,ksph)
                  vtij = vij * params % kappa
                  !call fmm_m2p_bessel_grad(vij * params % kappa, &
                  !    & params % rsph(ksph)*params % kappa, &
                  !    & constants % lmax0, &
                  !    & constants % vscales, -params % kappa, diff1(:, ksph), one, &
                  !    & diff_ep_dim3(:, icav))
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, diff1_grad(:, 1, ksph), one, &
                      & diff_ep_dim3(1, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, diff1_grad(:, 2, ksph), one, &
                      & diff_ep_dim3(2, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, diff1_grad(:, 3, ksph), one, &
                      & diff_ep_dim3(3, icav), work_complex, work)
              end do
            end do
            icav = zero
            do isph = 1, params % nsph
              do igrid =1, params % ngrid
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  do ind = 1, constants % nbasis
                    sum_dim3(:,ind,isph) = sum_dim3(:,ind,isph) + &
                        & constants % coefvec(igrid, ind, isph)*diff_ep_dim3(:,icav)
                  end do
                end if
              end do
            end do
            ! Computation of derivative of \bf(k)_j^l0(x_in)\times Y^j_l0m0(x_in)
            do isph = 1, params % nsph
              do ind = 1, constants % nbasis
                force(:, ksph) = force(:, ksph) &
                    & + sum_dim3(:, ind, isph)*(Xadj_r(ind, isph) &
                    & + Xadj_e(ind, isph))
              end do
            end do
        else
            call dgemv('T', (constants % lmax0+2)**2, 3, params % kappa, &
                & diff1_grad(1, 1, ksph), (constants % lmax0+2)**2, &
                & workspace % tmp_sph2(1, ksph), 1, one, force(1, ksph), 1)
        end if
    end do
    deallocate(diff_ep_dim3, phi_in, sum_dim3, diff1_grad, l2l_grad, &
        & diff0, diff1, stat=istat)
    if (istat.ne.0) stop 1
end subroutine fdouky

!
! Subroutine to compute the derivative of U_i(x_in) and \bf{k}^j_l0(x_in)Y^j_l0m0(x_in)
! in F0
! fdouky_f0 : Force Derivative of U_i^e(x_in), k_l0, and Y_l0m0 in F0
! @param[in] ddx_data     : Data Type
! @param[in] ksph         : Derivative with respect to x_k
! @param[in] sol_sgrid    : Solution of the Adjoint problem evaluated at the grid
! @param[in] sol_adj      : Adjoint solution
! @param[in] gradpsi      : Gradient of Psi_0
! @param[inout] force     : Force
subroutine fdouky_f0(params, constants, workspace, sol_adj, sol_sgrid, &
    & gradpsi, force)
    ! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: &
        & sol_adj
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: sol_sgrid
    real(dp), dimension(3, constants % ncav), intent(in) :: gradpsi
    real(dp), dimension(3, params % nsph), intent(inout) :: force

    ! local
    real(dp), external :: dnrm2
    integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, icav, ksph, &
        & knode, jnode, jsph_near, knear, jsph_node, indl, inode, istat
    ! val_dim3 : Intermediate value array of dimension 3
    real(dp), dimension(3) :: sij, vij, val_dim3, vtij
    ! val     : Intermediate variable to compute diff_ep
    ! f1      : Intermediate variable for derivative of coefY_der
    ! f2      : Intermediate variable for derivative of coefY_der
    ! nderpsi : Derivative of psi on grid points
    real(dp) :: val, f1, f2, nderpsi, sum_int

    ! local allocatable
    ! phi_in : sum_{j=1}^N diff0_j * coefY_j
    ! diff_ep_dim3 : 3 dimensional couterpart of diff_ep
    ! sum_dim3 : Storage of sum
    ! coefY_der : Derivative of k_l0 and Y_l0m0
    ! Debug purpose
    ! These variables can be taken from the subroutine update_rhs
    ! diff0       : dot_product([PU_j]_l0m0^l'm', l'/r_j[Xr]_jl'm' -
    !                        (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm')
    ! sum_Sjin : \sum_j [S]_{jin} Eq.~(97) [QSM20.SISC]
    ! c0 : \sum_{n=1}^N_g w_n U_j^{x_nj}\partial_n psi_0(x_nj)Y_{l0m0}(s_n)
    real(dp), allocatable :: phi_in(:,:), diff_ep_dim3(:,:), &
        & sum_dim3(:,:,:), coefy_der(:,:,:,:), diff0(:,:), sum_sjin(:,:), &
        & c0_d(:,:), c0_d1(:,:), c0_d1_grad(:,:,:), l2l_grad(:,:,:)

    real(dp) :: termi, termk, rijn
    ! vplm   : Argument to call ylmbas
    real(dp),  dimension(constants % nbasis):: basloc, vplm
    real(dp),  dimension(3, constants % nbasis):: dbasloc
    ! vcos   : Argument to call ylmbas
    ! vsin   : Argument to call ylmbas
    real(dp),  dimension(params % lmax+1):: vcos, vsin
    real(dp), dimension(0:params % lmax) :: SK_rijn, DK_rijn
    complex(dp) :: work_complex(constants % lmax0 + 2)
    real(dp) :: work(constants % lmax0 + 2)

    allocate(phi_in(params % ngrid, params % nsph), &
        & diff_ep_dim3(3, constants % ncav), &
        & sum_dim3(3, constants % nbasis, params % nsph), &
        & coefY_der(3, constants % ncav, constants % nbasis0, params % nsph), &
        & c0_d(constants % nbasis0, params % nsph), &
        & c0_d1(constants % nbasis0, params % nsph), &
        & c0_d1_grad((constants % lmax0+2)**2, 3, params % nsph), &
        & sum_Sjin(params % ngrid, params % nsph), &
        & l2l_grad((params % pl+2)**2, 3, params % nsph), &
        & diff0(constants % nbasis0, params % nsph), stat=istat)
    if (istat.ne.0) stop 1


    ! Setting initial values to zero
    SK_rijn = zero
    DK_rijn = zero
    coefY_der = zero
    c0_d = zero
    c0_d1 = zero
    c0_d1_grad = zero

    icav = zero
    do isph = 1, params % nsph
      do igrid= 1, params % ngrid
        if ( constants % ui(igrid,isph) .gt. zero ) then
          icav = icav + 1
          nderpsi = dot_product( gradpsi(:,icav),constants % cgrid(:,igrid) )
          c0_d(:, isph) = c0_d(:,isph) + &
                       & constants % wgrid(igrid)* &
                       & constants % ui(igrid,isph)*&
                       & nderpsi* &
                       & constants % vgrid(1:constants % nbasis0,igrid)
          do l0 = 0, constants % lmax0
              ind0 = l0*l0 + l0 + 1
              c0_d1(ind0-l0:ind0+l0, isph) = c0_d(ind0-l0:ind0+l0, isph) * &
                  & constants % C_ik(l0, isph)
          end do
          ! Prepare c0_d1_grad
          call fmm_m2m_bessel_grad(constants % lmax0, constants % SK_ri(:, isph), &
              & constants % vscales, constants % vcnk, c0_d1(:, isph), &
              & c0_d1_grad(:, :, isph))
        end if
      end do
    end do

    if (params % fmm .eq. 0) then
        ! Compute [S]_{jin}
        icav = 0
        do isph = 1, params % nsph
          do igrid = 1, params % ngrid
            if (constants % ui(igrid,isph).gt.zero) then
              icav = icav + 1
              sum_int = zero
              ! Loop to compute Sijn
              do jsph = 1, params % nsph
                vij = params % csph(:, isph) + &
                    & params % rsph(isph)*constants % cgrid(:, igrid) - &
                    & params % csph(:, jsph)
                vtij = vij * params % kappa
                call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                    & constants % vscales, constants % SK_ri(:, jsph), one, &
                    & c0_d1(:, jsph), one, sum_int, work_complex, work)
              end do ! End of loop jsph
              sum_Sjin(igrid,isph) = -(params % epsp/params % eps)*sum_int
            end if
          end do ! End of loop igrid
        end do ! End of loop isph
    else
        ! Load input harmonics into tree data
        workspace % tmp_sph = zero
        workspace % tmp_sph(1:constants % nbasis0, :) = c0_d1(:, :)
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(1:constants % nbasis0, inode) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph)
                workspace % tmp_node_m(constants % nbasis0+1:, inode) = zero
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(:, inode) = workspace % tmp_sph(1:indl, isph)
            end do
        end if
        ! Do FMM operations
        call tree_m2m_bessel_rotation(params, constants, workspace % tmp_node_m)
        call tree_m2l_bessel_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation(params, constants, workspace % tmp_node_l)
        call tree_l2p_bessel(params, constants, -params % epsp/params % eps, workspace % tmp_node_l, zero, &
            & sum_Sjin)
        call tree_m2p_bessel(params, constants, constants % lmax0, -params % epsp/params % eps, &
            & params % lmax, workspace % tmp_sph, one, &
            & sum_Sjin)
        ! Make phi_in zero at internal grid points
        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, isph) .eq. zero) then
                    sum_Sjin(igrid, isph) = zero
                end if
            end do
        end do
        ! Get gradients of the L2L
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_sph_l(:, isph) = workspace % tmp_node_l(:, inode)
            call fmm_l2l_bessel_grad(params % pl, &
                & constants % SI_ri(:, isph), constants % vscales, &
                & constants % vcnk, workspace % tmp_node_l(:, inode), &
                & l2l_grad(:, :, isph))
        end do
        workspace % tmp_sph = sol_adj
        call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
            & one, constants % vwgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, &
            & workspace % tmp_grid, params % ngrid)
        workspace % tmp_grid = workspace % tmp_grid * constants % ui
        ! Adjoint FMM with output tmp_sph2(:, :) which stores coefficients of
        ! harmonics of degree up to lmax+1
        call tree_m2p_bessel_nodiag_adj(params, constants, constants % lmax0+1, one, &
            & workspace % tmp_grid, zero, params % lmax+1, workspace % tmp_sph2)
        call tree_l2p_bessel_adj(params, constants, one, workspace % tmp_grid, zero, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_node_m)
        call tree_m2m_bessel_rotation_adj(params, constants, workspace % tmp_node_m)
        ! Properly load adjoint multipole harmonics into tmp_sph2 that holds
        ! harmonics of a degree up to lmax+1
        if(constants % lmax0+1 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph2(1:(constants % lmax0+2)**2, isph) = &
                    & workspace % tmp_sph2(1:(constants % lmax0+2)**2, isph) + &
                    & workspace % tmp_node_m(1:(constants % lmax0+2)**2, inode)
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph2(1:indl, isph) = &
                    & workspace % tmp_sph2(1:indl, isph) + &
                    & workspace % tmp_node_m(:, inode)
            end do
        end if
    end if
    do ksph = 1, params % nsph
        ! Computation of derivative of U_i^e(x_in)
        call fdoga(params, constants, ksph, sol_sgrid, sum_Sjin, force(:, ksph))
        ! Aleksandr: my loop for the diff_ep_dim3
        diff_ep_dim3 = zero
        ! At first isph=ksph, jsph!=ksph
        icav = constants % icav_ia(ksph) - 1
        if (params % fmm .eq. 0) then
            do igrid = 1, params % ngrid
              if (constants % ui(igrid, ksph) .eq. zero) cycle
              icav = icav + 1
              do jsph = 1, params % nsph
                if (jsph .eq. ksph) cycle
                  vij  = params % csph(:,ksph) + &
                      & params % rsph(ksph)*constants % cgrid(:,igrid) - &
                      & params % csph(:,jsph)
                  vtij = vij * params % kappa
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, c0_d1_grad(:, 1, jsph), one, &
                      & diff_ep_dim3(1, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, c0_d1_grad(:, 2, jsph), one, &
                      & diff_ep_dim3(2, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, c0_d1_grad(:, 3, jsph), one, &
                      & diff_ep_dim3(3, icav), work_complex, work)
              end do
            end do
        else
            knode = constants % snode(ksph)
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, ksph) .eq. zero) cycle
                icav = icav + 1
                ! Far-field
                call dgemv('T', (params % pl+2)**2, 3, -params % kappa, &
                    & l2l_grad(1, 1, ksph), &
                    & (params % pl+2)**2, constants % vgrid(1, igrid), 1, &
                    & one, diff_ep_dim3(1, icav), 1)
                ! Near-field
                do knear = constants % snear(knode), constants % snear(knode+1)-1
                    jnode = constants % near(knear)
                    do jsph_node = constants % cluster(1, jnode), &
                        & constants % cluster(2, jnode)
                        jsph = constants % order(jsph_node)
                        if (jsph .eq. ksph) cycle
                        vij = params % csph(:, ksph) - params % csph(:, jsph) + &
                            & params % rsph(ksph)*constants % cgrid(:, igrid)
                        vtij = vij * params % kappa
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, c0_d1_grad(:, 1, jsph), one, &
                            & diff_ep_dim3(1, icav), work_complex, work)
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, c0_d1_grad(:, 2, jsph), one, &
                            & diff_ep_dim3(2, icav), work_complex, work)
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, c0_d1_grad(:, 3, jsph), one, &
                            & diff_ep_dim3(3, icav), work_complex, work)
                    end do
                end do
            end do
        end if
        sum_dim3 = zero
        icav = constants % icav_ia(ksph) - 1
        do igrid =1, params % ngrid
            if(constants % ui(igrid, ksph) .gt. zero) then
                icav = icav + 1
                do ind = 1, constants % nbasis
                    sum_dim3(:,ind,ksph) = sum_dim3(:,ind,ksph) + &
                        & -(params % epsp/params % eps)* &
                        & constants % coefvec(igrid, ind, ksph)*diff_ep_dim3(:,icav)
                end do
            end if
        end do
        do ind = 1, constants % nbasis
            force(:, ksph) = force(:, ksph) + &
                & sum_dim3(:, ind, ksph)*sol_adj(ind, ksph)
        end do

        ! Now jsph=ksph and isph!=ksph
        if (params % fmm .eq. 0) then
            diff_ep_dim3 = zero
            sum_dim3 = zero
            do isph = 1, params % nsph
              if (isph .eq. ksph) cycle
              icav = constants % icav_ia(isph) - 1
              do igrid = 1, params % ngrid
                  if (constants % ui(igrid, isph) .eq. zero) cycle
                  icav = icav + 1
                  vij  = params % csph(:,isph) + &
                      & params % rsph(isph)*constants % cgrid(:,igrid) - &
                      & params % csph(:,ksph)
                  vtij = vij * params % kappa
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, c0_d1_grad(:, 1, ksph), one, &
                      & diff_ep_dim3(1, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, c0_d1_grad(:, 2, ksph), one, &
                      & diff_ep_dim3(2, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, c0_d1_grad(:, 3, ksph), one, &
                      & diff_ep_dim3(3, icav), work_complex, work)
              end do
            end do

            icav = zero
            do isph = 1, params % nsph
              do igrid =1, params % ngrid
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  do ind = 1, constants % nbasis
                    sum_dim3(:,ind,isph) = sum_dim3(:,ind,isph) + &
                                          & -(params % epsp/params % eps)* &
                                          & constants % coefvec(igrid, ind, isph)*diff_ep_dim3(:,icav)
                  end do
                end if
              end do
            end do

            ! Computation of derivative of \bf(k)_j^l0(x_in)\times Y^j_l0m0(x_in)
            do isph = 1, params % nsph
              do ind = 1, constants % nbasis
                force(:, ksph) = force(:, ksph) + sum_dim3(:, ind, isph)*sol_adj(ind, isph)
              end do
            end do
        else
            call dgemv('T', (constants % lmax0+2)**2, 3, -params % epsp/params % eps*params % kappa, &
                & c0_d1_grad(1, 1, ksph), (constants % lmax0+2)**2, &
                & workspace % tmp_sph2(1, ksph), 1, one, force(1, ksph), 1)
        end if
    end do

    deallocate(phi_in, diff_ep_dim3, sum_dim3, coefY_der, c0_d, c0_d1, &
          & c0_d1_grad, sum_Sjin, l2l_grad, diff0, stat=istat)
    if (istat.ne.0) stop 1

end subroutine fdouky_f0


! Subroutine to calculate the third derivative term in C1_C2 matrix,
! namely the derivative of PU_i
! @param[in]  ddx_data     : Input data file
! @param[in]  Xr           : Solution of the Laplace problem
! @param[in]  Xe           : Solution of the HSP problem
! @param[in]  Xadj_r_sgrid : Adjoint Laplace solution evaluated at grid point
! @param[in]  Xadj_e_sgrid : Adjoint HSP solution evaluated at grid point
! @param[in]  diff_re      : l'/r_j[Xr]_jl'm' -(i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
! @param[out] force        : Force
subroutine derivative_P(params, constants, workspace, Xr, Xe, Xadj_r_sgrid, &
    & Xadj_e_sgrid, diff_re, force)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xr, &
        & Xe, diff_re
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: &
        & Xadj_r_sgrid, Xadj_e_sgrid
    real(dp), dimension(3, params % nsph), intent(inout) :: force
    ! Local variable
    ! igrid0: Index for grid point n0
    integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, igrid0, icav, &
        & indl, inode, istat
    ! term  : SK_rijn/SK_rj
    ! termi : DI_ri/SI_ri
    ! termk : DK_ri/SK_ri
    ! sum_int : Intermediate sum
    ! sum_r   : Intermediate sum for Laplace
    ! sum_e   : Intermediate sum for HSP
    real(dp) :: term, termi, termk, sum_int, sum_r, sum_e
    real(dp) :: rijn
    real(dp)  :: vij(3), sij(3), vtij(3)

    ! local allocatable
    ! phi_n_r : Phi corresponding to Laplace problem
    ! phi_n_e : Phi corresponding to HSP problem
    ! coefY_d : sum_{l0m0} C_ik*term*Y_l0m0^j(x_in)*Y_l0m0(s_n)
    ! diff_re_sgrid : diff_re evaluated at grid point
    real(dp), allocatable :: phi_n_r(:,:), phi_n_e(:,:), coefY_d(:,:,:), &
        & diff_re_sgrid(:,:)

    ! basloc : Y_lm(s_n)
    ! vplm   : Argument to call ylmbas
    real(dp),  dimension(constants % nbasis):: basloc, vplm
    ! dbasloc : Derivative of Y_lm(s_n)
    real(dp),  dimension(3, constants % nbasis):: dbasloc
    ! vcos   : Argument to call ylmbas
    ! vsin   : Argument to call ylmbas
    real(dp),  dimension(params % lmax+1):: vcos, vsin
    ! SK_rijn : Besssel function of first kind for rijn
    ! DK_rijn : Derivative of Besssel function of first kind for rijn
    real(dp), dimension(0:params % lmax) :: SK_rijn, DK_rijn
    real(dp) :: coef(constants % nbasis0), work(constants % lmax0+1)
    complex(dp) :: work_complex(constants % lmax0+1)

    ! TODO: coefY_d is N^2 storage, it prevents the application of ddLPB forces
    ! to large systems. Probably it is not needed and should be removed one day.
    allocate(phi_n_r(params % ngrid, params % nsph), &
        & phi_n_e(params % ngrid, params % nsph), &
        & coefY_d(constants % ncav, params % ngrid, params % nsph), &
        & diff_re_sgrid(params % ngrid, params % nsph), stat=istat)
    if (istat.ne.0) stop 1

    sum_int = zero
    sum_r = zero
    sum_e = zero
    phi_n_r = zero
    phi_n_e = zero
    coefY_d = zero
    diff_re_sgrid = zero
    basloc = zero
    vplm = zero
    dbasloc = zero
    vcos = zero
    vsin = zero
    SK_rijn = zero
    DK_rijn = zero

    if (params % fmm .eq. 0) then
        ! Compute  summation over l0, m0
        ! Loop over the sphers j
        do jsph = 1, params % nsph
          ! Loop over the grid points n0
          do igrid0 = 1, params % ngrid
            icav = zero
            ! Loop over spheres i
            do isph = 1, params % nsph
              ! Loop over grid points n
              do igrid = 1, params % ngrid
                ! Check for U_i^{eta}(x_in)
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  vij  = params % csph(:,isph) + &
                         & params % rsph(isph)*constants % cgrid(:,igrid) - &
                         & params % csph(:,jsph)
                  rijn = sqrt(dot_product(vij,vij))
                  sij = vij/rijn

                  do l0 = 0, constants % lmax0
                    do m0 = -l0,l0
                      ind0 = l0**2 + l0 + m0 + 1
                      coef(ind0) = constants % vgrid(ind0, igrid0) * &
                          & constants % C_ik(l0, jsph)
                    end do
                  end do
                  vtij = vij*params % kappa
                  call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                      & constants % vscales, constants % SK_ri(:, jsph), one, &
                      & coef, zero, coefY_d(icav, igrid0, jsph), work_complex, work)
                end if
              end do
            end do
          end do
        end do

        ! Compute phi_in
        ! Loop over spheres j
        do jsph = 1, params % nsph
          ! Loop over grid points n0
          do igrid0 = 1, params % ngrid
            icav = zero
            sum_r = zero
            sum_e = zero
            ! Loop over sphers i
            do isph = 1, params % nsph
              ! Loop over grid points n
              do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  sum_r = sum_r + coefY_d(icav, igrid0, jsph)*Xadj_r_sgrid(igrid, isph) &
                          & * constants % wgrid(igrid)*constants % ui(igrid, isph)
                  sum_e = sum_e + coefY_d(icav, igrid0, jsph)*Xadj_e_sgrid(igrid, isph) &
                          & * constants % wgrid(igrid)*constants % ui(igrid, isph)
                end if
              end do
            end do
            phi_n_r(igrid0, jsph) = sum_r
            phi_n_e(igrid0, jsph) = sum_e
          end do
        end do
    else
        ! Compute phi_n_r at first
        ! Adjoint integration from spherical harmonics to grid points is not needed
        ! here as ygrid already contains grid values, we just need to scale it by
        ! weights of grid points
        do isph = 1, params % nsph
            workspace % tmp_grid(:, isph) = Xadj_r_sgrid(:, isph) * &
                & constants % wgrid(:) * constants % ui(:, isph)
        end do
        ! Adjoint FMM
        call tree_m2p_bessel_adj(params, constants, constants % lmax0, one, &
            & workspace % tmp_grid, zero, params % lmax, workspace % tmp_sph)
        call tree_l2p_bessel_adj(params, constants, one, workspace % tmp_grid, &
            & zero, workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l, workspace % tmp_node_m)
        call tree_m2m_rotation_adj(params, constants, workspace % tmp_node_m)
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
        call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis0, &
            & one, constants % vgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, phi_n_r, params % ngrid)
        ! Compute phi_n_e now
        ! Adjoint integration from spherical harmonics to grid points is not needed
        ! here as ygrid already contains grid values, we just need to scale it by
        ! weights of grid points
        do isph = 1, params % nsph
            workspace % tmp_grid(:, isph) = Xadj_e_sgrid(:, isph) * &
                & constants % wgrid(:) * constants % ui(:, isph)
        end do
        ! Adjoint FMM
        call tree_m2p_bessel_adj(params, constants, constants % lmax0, one, &
            & workspace % tmp_grid, zero, params % lmax, workspace % tmp_sph)
        call tree_l2p_bessel_adj(params, constants, one, workspace % tmp_grid, &
            & zero, workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l, workspace % tmp_node_m)
        call tree_m2m_rotation_adj(params, constants, workspace % tmp_node_m)
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
        call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis0, &
            & one, constants % vgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, phi_n_e, params % ngrid)
    end if
    call dgemm('T', 'N', params % ngrid, params % nsph, &
              & constants % nbasis, one, constants % vgrid, constants % vgrid_nbasis, &
              & diff_re , constants % nbasis, zero, diff_re_sgrid, &
              & params % ngrid)
    do isph = 1, params % nsph
        call fdoga(params, constants, isph, diff_re_sgrid, phi_n_r, force(:, isph))
        call fdoga(params, constants, isph, diff_re_sgrid, phi_n_e, force(:, isph))
    end do

    deallocate( phi_n_r, phi_n_e, coefY_d, diff_re_sgrid, stat=istat)
    if (istat.ne.0) stop 1
  end subroutine derivative_P

! Subroutine to calculate the derivative of C_{0l0m0}^j
! @param[in]  ddx_data           : Input data file
! @param[in]  ksph               : Derivative with respect to x_k
! @param[in]  sol_sgrid          : Solution of the Adjoint problem evaluated at the grid
! @param[in]  gradpsi            : Gradient of Psi_0
! @param[in]  normal_hessian_cav : Normal of the Hessian evaluated at cavity points
! @param[in]  icav_g             : Index of outside cavity point
! @param[out] force              : Force corresponding to HSP problem
subroutine fdoco(params, constants, workspace, sol_sgrid, gradpsi, &
    & normal_hessian_cav, icav_g, force)
    ! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: sol_sgrid
    real(dp), dimension(3, constants % ncav), intent(in) :: gradpsi
    real(dp), dimension(3, constants % ncav), intent(in) :: normal_hessian_cav
    integer, intent(inout) :: icav_g
    real(dp), dimension(3, params % nsph), intent(inout) :: force

    ! local
    integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, igrid0, icav, &
        & indl, inode, istat
    ! term  : SK_rijn/SK_rj
    ! termi : DI_ri/SI_ri
    ! termk : DK_ri/SK_ri
    ! sum_int : Intermediate sum
    ! hessian_contribution :
    real(dp) :: term, termi, termk, sum_int, hessian_contribution(3), nderpsi
    real(dp) :: rijn
    real(dp)  :: vij(3), sij(3), vtij(3)

    ! local allocatable
    ! phi_n : Phi corresponding to Laplace problem
    ! coefY_d : sum_{l0m0} C_ik*term*Y_l0m0^j(x_in)*Y_l0m0(s_n)
    ! gradpsi_grid : gradpsi evaluated at grid point
    real(dp), allocatable :: phi_n(:,:), phi_n2(:,:), coefY_d(:,:,:), &
        & gradpsi_grid(:,:)

    ! vplm   : Argument to call ylmbas
    real(dp),  dimension(constants % nbasis):: basloc, vplm
    real(dp),  dimension(3, constants % nbasis):: dbasloc
    ! vcos   : Argument to call ylmbas
    ! vsin   : Argument to call ylmbas
    real(dp),  dimension(params % lmax+1):: vcos, vsin
    real(dp), dimension(0:params % lmax) :: SK_rijn, DK_rijn
    real(dp) :: coef(constants % nbasis0), work(constants % lmax0+1)
    complex(dp) :: work_complex(constants % lmax0+1)
    real(dp), external :: dnrm2

    allocate(phi_n(params % ngrid, params % nsph), &
        & phi_n2(params % ngrid, params % nsph), &
        & coefY_d(constants % ncav, params % ngrid, params % nsph), &
        & gradpsi_grid(params % ngrid, params % nsph), stat=istat)
    if (istat.ne.0) stop 1

    ! Intial allocation of vectors
    sum_int = zero
    phi_n = zero
    coefY_d = zero
    gradpsi_grid = zero
    basloc = zero
    vplm = zero
    dbasloc = zero
    vcos = zero
    vsin = zero
    SK_rijn = zero
    DK_rijn = zero

    if (params % fmm .eq. 0) then
        ! Compute  summation over l0, m0
        ! Loop over the sphers j
        do jsph = 1, params % nsph
          ! Loop over the grid points n0
          do igrid0 = 1, params % ngrid
            icav = zero
            ! Loop over spheres i
            do isph = 1, params % nsph
              ! Loop over grid points n
              do igrid = 1, params % ngrid
                ! Check for U_i^{eta}(x_in)
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  vij  = params % csph(:,isph) + &
                         & params % rsph(isph)*constants % cgrid(:,igrid) - &
                         & params % csph(:,jsph)
                  rijn = sqrt(dot_product(vij,vij))
                  sij = vij/rijn
                  ! Loop over l0
                  do l0 = 0, constants % lmax0
      !              term = SK_rijn(l0)/constants % SK_ri(l0,jsph)
                    ! Loop over m0
                    do m0 = -l0,l0
                      ind0 = l0**2 + l0 + m0 + 1
                      coef(ind0) = constants % vgrid(ind0, igrid0) * &
                          & constants % C_ik(l0, jsph)
                    end do ! End of loop m0
                  end do! End of loop l0
                  vtij = vij*params % kappa
                  call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                      & constants % vscales, constants % SK_ri(:, jsph), one, &
                      & coef, zero, coefY_d(icav, igrid0, jsph), work_complex, work)
                end if
              end do ! End of loop igrid
            end do! End of loop isph
          end do ! End of loop igrid0
        end do ! End of loop jsph
        ! Compute phi_in
        ! Loop over spheres j
        do jsph = 1, params % nsph
          ! Loop over grid points n0
          do igrid0 = 1, params % ngrid
            icav = zero
            sum_int = zero
            ! Loop over sphers i
            do isph = 1, params % nsph
              ! Loop over grid points n
              do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  sum_int = sum_int + coefY_d(icav, igrid0, jsph)*sol_sgrid(igrid, isph) &
                          & * constants % wgrid(igrid)*constants % ui(igrid, isph)
                end if
              end do
            end do
            phi_n(igrid0, jsph) = -(params % epsp/params % eps)*sum_int
          end do! End of loop j
        end do ! End of loop igrid
    else
        ! Adjoint integration from spherical harmonics to grid points is not needed
        ! here as ygrid already contains grid values, we just need to scale it by
        ! weights of grid points
        do isph = 1, params % nsph
            workspace % tmp_grid(:, isph) = sol_sgrid(:, isph) * &
                & constants % wgrid(:) * constants % ui(:, isph)
        end do
        ! Adjoint FMM
        call tree_m2p_bessel_adj(params, constants, constants % lmax0, one, &
            & workspace % tmp_grid, zero, params % lmax, workspace % tmp_sph)
        call tree_l2p_bessel_adj(params, constants, one, workspace % tmp_grid, &
            & zero, workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l, workspace % tmp_node_m)
        call tree_m2m_rotation_adj(params, constants, workspace % tmp_node_m)
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
        call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis0, &
            & -params % epsp/params % eps, constants % vgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, phi_n, params % ngrid)
    end if
    icav = zero
    do isph = 1, params % nsph
      do igrid = 1, params % ngrid
        if(constants % ui(igrid, isph) .gt. zero) then
          icav = icav + 1
          nderpsi = dot_product( gradpsi(:,icav),constants % cgrid(:,igrid) )
          gradpsi_grid(igrid, isph) = nderpsi
        end if
      end do ! End of loop igrid
    end do ! End of loop i
    do isph = 1, params % nsph
        call fdoga(params, constants, isph, gradpsi_grid, phi_n, force(:, isph))
        ! Compute the Hessian contributions
        do igrid = 1, params % ngrid
          if(constants % ui(igrid, isph) .gt. zero) then
            icav_g = icav_g + 1
            force(:, isph) = force(:, isph) + constants % wgrid(igrid)*constants % ui(igrid, isph)*&
                             & phi_n(igrid, isph)*normal_hessian_cav(:, icav_g)
          end if
        end do
        call fdops(params, constants, workspace, isph, phi_n, force(:, isph))
    end do

    deallocate(phi_n, phi_n2, coefY_d, gradpsi_grid, stat=istat)
    if (istat.ne.0) stop 1
end subroutine fdoco

  !
  ! fdops : Force derivative of potential at spheres
  ! @param[in]  ddx_data : Input data file
  ! @param[in]  phi_n    : phi_n^j
  ! @param[in]  ksph     : Derivative with respect to x_k
  ! @param[out] force    : Force
  subroutine fdops(params, constants, workspace, ksph, phi_n, force)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    !! Temporary buffers
    type(ddx_workspace_type), intent(inout) :: workspace
  integer, intent(in) :: ksph
  real(dp),  dimension(params % ngrid, params % nsph), intent(in)    :: phi_n
  real(dp),  dimension(3), intent(inout) :: force
  !
  ! Local variables
  integer :: jsph, i,j, igrid
  ! sum_int            : Intermediate sum
  ! normal_hessian_cav : Normal derivative of Hessian of psi
  ! vij_vijT           : vij*vij^T
  real(dp)  :: sum_int(3), vij(3), rijn, normal_hessian_cav(3), vij_vijT(3,3)
  ! identity_matrix : Identity matrix of size 3x3
  ! hessianv        : Hessian of psi evaluated by centers
  real(dp) :: identity_matrix(3,3), hessianv(3,3)
  real(dp), external :: dnrm2
  ! Create Identity matrix
  identity_matrix = zero
  do i = 1, 3
    identity_matrix(i,i) = one
  end do

  sum_int = zero
  ! Loop over spheres
  do jsph = 1, params % nsph
    ! Loop over grid points
    do igrid = 1, params % ngrid
      if(constants % ui(igrid, jsph) .ne. zero) then
        vij_vijT = zero
        hessianv = zero
        normal_hessian_cav = zero
        vij = zero
        vij = params % csph(:, jsph) + &
           & params % rsph(jsph)*constants % cgrid(:, igrid) - &
           & params % csph(:, ksph)
        do i = 1,3
          do j = 1,3
            vij_vijT(i,j) = vij(i)*vij(j)
          end do
        end do
        rijn = dnrm2(3, vij, 1)
        hessianv = 3*vij_vijT/(rijn**5)- &
                 & identity_matrix/(rijn**3)
        do i = 1, 3
          normal_hessian_cav = normal_hessian_cav + &
                             & hessianv(:,i)*constants % cgrid(i,igrid)
        end do
        sum_int = sum_int + &
                 & constants % ui(igrid, jsph)* &
                 & constants % wgrid(igrid)* &
                 & phi_n(igrid, jsph)* &
                 & normal_hessian_cav
      end if
    end do
  end do
  force = force - params % charge(ksph)*sum_int
  end subroutine fdops

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
