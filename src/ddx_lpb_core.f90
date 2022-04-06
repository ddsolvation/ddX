!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!!
!! @version 1.0.0
!! @author Abhinav Jha and Michele Nottoli
!! @date 2022-03-31

!> Core routines and parameters specific to ddLPB
module ddx_lpb_core
! Get the core-routines
use ddx_core
!
contains

!> Find intermediate F0 in the RHS of the ddLPB
!! model given in Eq.(82)
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  workspace  : Input workspace
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

!> Intermediate computation of Bx
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]      isph   : Number of the sphere
!! @param[in, out] pot    : Array of size ngrid
!! @param[in]      x      : Input vector (Usually X_e)
!! @param[in, out] basloc : Used to compute spherical harmonic
!! @param[in, out] vplm   : Used to compute spherical harmonic
!! @param[in, out] vcos   : Used to compute spherical harmonic
!! @param[in, out] vsin   : Used to compute spherical harmonic
!!
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

!> Scale the ddCOSMO solution vector
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  direction  : Direction of the scaling
!! @param[in, out] vector : ddCOSMO solution vector
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


!> Intermediate calculation in calcv2_lpb subroutine
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  rijn       : Radius of sphers x_ijn
!! @param[in]  ri         : Radius of sphers x_i
!! @param[in]  isph       : Index of sphere
!! @param[in]  basloc     : Spherical Harmonic
!! @param[out] fac_hsp    : Return bessel function ratio multiplied by
!!                       the spherical harmonic Y_l'm'. Array of size nylm
!!
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


!> Intermediate for computing Bstarx
!! Compute the Adjoint matix B*x
!! Taken from ddx_core routine adjrhs
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  workspace  : Input workspace
!! @param[in]      isph   : Number of the sphere
!! @param[in, out] basloc : Used to compute spherical harmonic
!! @param[in, out] vplm   : Used to compute spherical harmonic
!! @param[in, out] vcos   : Used to compute spherical harmonic
!! @param[in, out] vsin   : Used to compute spherical harmonic
!!
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


!> Intermediate calculation in adjrhs_lpb subroutine
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  rjin       : Radius of sphers x_jin
!! @param[in]  rj         : Radius of sphers x_j
!! @param[in]  jsph       : Index of sphere
!! @param[in]  basloc     : Spherical Harmonic
!! @param[out] fac_hsp    : Return bessel function ratio multiplied by
!!                       the spherical harmonic Y_l'm'. Array of size nylm
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
end module ddx_lpb_core
