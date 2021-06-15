!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_constants.f90
!! This file contains the type ddx_constants_type to hold all required run-time
!! constants and routines to actually compute these constants. These constants
!! include scaling factors for spherical harmonics, certain combinatorial
!! numbers, hierarchical trees for the FMM and so on. Everything that can be
!! precomputed to reduce ddX timings shall be stored in the ddx_constants_type.
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-06-15

!> Run-time constants
module ddx_constants
! Get ddx_params_type and all compile-time definitions
use ddx_parameters
! Get harmonics-related functions
use ddx_harmonics

! Disable implicit types
implicit none

type ddx_constants_type
    !> Number modeling spherical harmonics per sphere.
    integer :: nbasis
    !> Total number of modeling degrees of freedom.
    !!
    !! This is equal to `nsph*nbasis`.
    integer :: n
    !> Maximal degree of used real normalized spherical harmonics.
    !!
    !! For example, if FMM is
    !! used, its M2L operation requires computing spherical harmonics of all
    !! degrees up to `pm+pl`. If `force=1` then this parameter might be
    !! increased by 1 because certain implementations of analytical gradients
    !! might rely on spherical harmonics of +1 degree.
    integer :: dmax
    !> Total number of used real spherical harmonics and a size of `vscales`.
    integer :: nscales
    !> Scales of real normalized spherical harmonics of degree up to dmax.
    !!
    !! This array has a dimension (nscales).
    real(dp), allocatable :: vscales(:)
    !> Array of values 4pi/(2l+1), dimension is (dmax+1). Referenced only if
    !!      fmm=1, but allocated and computed in any case.
    real(dp), allocatable :: v4pi2lp1(:)
    !> Relative scales of real normalized spherical harmonics.
    !!
    !! Each values is multiplied by a corresponding 4pi/(2l+1). Dimension of
    !! this array is (nscales). Referenced only if fmm=1, but allocated and
    !! computed in any case.
    real(dp), allocatable :: vscales_rel(:)
    !> Number of precomputed square roots of factorials.
    !!
    !! Just like with `dmax` parameter, number of used factorials is either
    !! `2*lmax+1` or `2*(pm+pl)+1` depending on whether FMM is used or not.
    integer :: nfact
    !> Array of square roots of factorials of a dimension (nfact).
    real(dp), allocatable :: vfact(:)
    !> Array of square roots of combinatorial numbers C_n^k.
    !!
    !! Dimension of this array is ((2*dmax+1)*(dmax+1)). Allocated, computed
    !! and referenced only if fmm=1.
    real(dp), allocatable :: vcnk(:)
    !> Array of common M2L coefficients for any OZ translation.
    !!
    !! This array has a dimension (pm+1, pl+1, pl+1). Allocated, computed and
    !! referenced only if fmm=1.
    real(dp), allocatable :: m2l_ztranslate_coef(:, :, :)
    !> Array of common M2L coefficients for any adjoint OZ translation.
    !!
    !! Dimension of this array is (pl+1, pl+1, pm+1). It is allocated, computed
    !! and referenced only if fmm=1.
    real(dp), allocatable :: m2l_ztranslate_adj_coef(:, :, :)
    !> Coordinates of Lebedev quadrature points of a dimension (3, ngrid).
    real(dp), allocatable :: cgrid(:, :)
    !> Weights of Lebedev quadrature points of a dimension (ngrid).
    real(dp), allocatable :: wgrid(:)
    !> Maximal degree of spherical harmonics evaluated at Lebedev grid points.
    !!
    !! Although we use spherical harmonics of degree up to `dmax`, only
    !! spherical harmonics of degree up to `lmax` and `pl` are evaluated
    !! at Lebedev grid points. In the case `force=1` this degrees might be
    !! increased by 1 depending on implementation of gradients.
    integer :: vgrid_dmax
    !> Number of spherical harmonics evaluated at Lebedev grid points.
    integer :: vgrid_nbasis
    !> Values of spherical harmonics at Lebedev grid points.
    !!
    !! Dimensions of this array are (vgrid_nbasis, ngrid)
    real(dp), allocatable :: vgrid(:, :)
    !> Weighted values of spherical harmonics at Lebedev grid points.
    !!
    !! vwgrid(:, igrid) = vgrid(:, igrid) * wgrid(igrid)
    !! Dimension of this array is (vgrid_nbasis, ngrid).
    real(dp), allocatable :: vwgrid(:, :)
    !> Values of L2P at grid points. Dimension is (vgrid_nbasis, ngrid).
    real(dp), allocatable :: l2grid(:, :)
    !> Upper limit on a number of neighbours per sphere. This value is just an
    !!      upper bound that is not guaranteed to be the actual maximum.
    integer :: nngmax
    !> List of intersecting spheres in a CSR format. Dimension is (nsph+1).
    integer, allocatable :: inl(:)
    !> List of intersecting spheres in a CSR format. Dimension is
    !!      (nsph*nngmax).
    integer, allocatable :: nl(:)
    !> Values of a characteristic function f at all grid points of all spheres.
    !!      Dimension is (ngrid, npsh).
    real(dp), allocatable :: fi(:, :)
    !> Values of a characteristic function U at all grid points of all spheres.
    !!      Dimension is (ngrid, nsph).
    real(dp), allocatable :: ui(:, :)
    !> Derivative of the characteristic function U at all grid points of all
    !!      spheres. Dimension is (3, ngrid, nsph).
    real(dp), allocatable :: zi(:, :, :)
    !> Number of external Lebedev grid points on a molecular surface.
    integer :: ncav
    !> Number of external Lebedev grid points on each sphere.
    integer, allocatable :: ncav_sph(:)
    !> Coordinates of external Lebedev grid points. Dimension is (3, ncav).
    real(dp), allocatable :: ccav(:, :)
    !> Row indexes in CSR format of all cavity points. Dimension is (nsph+1).
    integer, allocatable :: icav_ia(:)
    !> Column indexes in CSR format of all cavity points. Dimension is (ncav).
    integer, allocatable :: icav_ja(:)
    !> Preconditioner for an operator R_eps. Allocated and computed only for
    !!      the PCM model (model=3). Dimension is (nbasis, nbasis, nsph).
    real(dp), allocatable :: rx_prc(:, :, :)
    !! Cluster tree information that is allocated and computed only if fmm=1.
    !> Reordering of spheres for better locality. This array has a dimension
    !!      (nsph) and is allocated/used only if fmm=1.
    integer, allocatable :: order(:)
    !> Number of clusters. Defined only if fmm=1.
    integer :: nclusters
    !> The first and the last spheres of each node. Dimension of this array is
    !!      (2, nclusters) and it is allocated/used only if fmm=1.
    integer, allocatable :: cluster(:, :)
    !> Children of each cluster. Dimension is (2, nclusters). Allocated and
    !!      used only if fmm=1.
    integer, allocatable :: children(:, :)
    !> Parent of each cluster. Dimension is (nclusters). Allocated and used
    !!      only if fmm=1.
    integer, allocatable :: parent(:)
    !> Center of bounding sphere of each cluster. Dimension is (3, nclusters).
    !!      This array is allocated and used only if fmm=1.
    real(dp), allocatable :: cnode(:, :)
    !> Radius of bounding sphere of each cluster. Dimension is (nclusters).
    !!      This array is allocated and used only if fmm=1.
    real(dp), allocatable :: rnode(:)
    !> Which leaf node contains only given input sphere. Dimension is (nsph).
    !!      This array is allocated and used only if fmm=1.
    integer, allocatable :: snode(:)
    !> Total number of far admissible pairs. Defined only if fmm=1.
    integer :: nnfar
    !> Total number of near admissible pairs. Defined only if fmm=1.
    integer :: nnnear
    !> Number of admissible far pairs for each node. Dimension is (nclusters).
    !!      This array is allocated and used only if fmm=1.
    integer, allocatable :: nfar(:)
    !> Number of admissible near pairs for each node. Dimension is (nclusters).
    !!      This array is allocated and used only if fmm=1.
    integer, allocatable :: nnear(:)
    !> Arrays of admissible far pairs. Dimension is (nnfar). This array is
    !!      allocated and used only if fmm=1.
    integer, allocatable :: far(:)
    !> Arrays of admissible near pairs. Dimension is (nnnear). This array is
    !!      allocated and used only if fmm=1.
    integer, allocatable :: near(:)
    !> Index of the first element of array of all admissible far pairs stored
    !!      in the array `far`. Dimension is (nclusters+1). This array is
    !!      allocated and used only if fmm=1.
    integer, allocatable :: sfar(:)
    !> Index of the first element of array of all admissible near pairs stored
    !!      in the array `near`. Dimension is (nclusters+1). This array is
    !!      allocated and used only if fmm=1.
    integer, allocatable :: snear(:)
    !> Number of near-field M2P interactions with cavity points. Defined only
    !!      if fmm=1.
    integer :: nnear_m2p
    !> Maximal degree of near-field M2P spherical harmonics. Defined only if
    !!      fmm=1.
    integer :: m2p_lmax
    !> Number of spherical harmonics used for near-field M2P. Defined only if
    !!      fmm=1.
    integer :: m2p_nbasis
    !> Number of spherical harmonics of degree up to lmax+1 used for
    !!      computation of forces (gradients). Allocated and used only if
    !!      fmm=1.
    integer :: grad_nbasis
    !> Flag if there were an error
    integer :: error_flag
    !> Last error message
    character(len=255) :: error_message
end type ddx_constants_type

contains

!> Compute all necessary constants
!! @param[in] params: Object containing all inputs.
!! @param[out] constants: Object containing all constants.
!! @param[out] info: flag of succesfull exit
!!      = 0: Succesfull exit
!!      = -1: params is in error state
!!      = 1: Allocation of memory failed.
subroutine constants_init(params, constants, info)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    !! Outputs
    type(ddx_constants_type), intent(out) :: constants
    integer, intent(out) :: info
    !! Local variables
    integer :: i, alloc_size, l, indl, igrid
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    real(dp), allocatable :: vplm(:), vcos(:), vsin(:)
    !! The code
    ! Check if params are OK
    if (params % error_flag .eq. 1) then
        constants % error_flag = 1
        constants % error_message = "params is in error state"
        info = -1
        return
    end if
    ! Maximal number of modeling spherical harmonics
    constants % nbasis = (params % lmax+1) ** 2
    ! Maximal number of modeling degrees of freedom
    constants % n = params % nsph * constants % nbasis
    ! Calculate dmax, vgrid_dmax, m2p_lmax, m2p_nbasis and grad_nbasis
    if (params % fmm .eq. 0) then
        constants % dmax = params % lmax
        constants % vgrid_dmax = params % lmax
        ! Other constants are not referenced if fmm=0
        constants % m2p_lmax = -1
        constants % m2p_nbasis = -1
        constants % grad_nbasis = -1
    else
        ! If forces are required then we need the M2P of a degree lmax+1 for
        ! the near-field analytical gradients
        if (params % force .eq. 1) then
            constants % dmax = max(params % pm+params % pl, &
                & params % lmax+1)
            constants % m2p_lmax = params % lmax + 1
            constants % grad_nbasis = (params % lmax+2) ** 2
        else
            constants % dmax = max(params % pm+params % pl, &
                & params % lmax)
            constants % m2p_lmax = params % lmax
            constants % grad_nbasis = -1
        end if
        constants % vgrid_dmax = max(params % pl, params % lmax)
        constants % m2p_nbasis = (constants % m2p_lmax+1) ** 2
    end if
    ! Compute sizes of vgrid, vfact and vscales
    constants % vgrid_nbasis = (constants % vgrid_dmax+1) ** 2
    constants % nfact = 2*constants % dmax+1
    constants % nscales = (constants % dmax+1) ** 2
    ! Allocate space for scaling factors of spherical harmonics
    allocate(constants % vscales(constants % nscales), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`vscales` allocation failed"
        info = 1
        return
    end if
    allocate(constants % v4pi2lp1(constants % dmax+1), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`v4pi2lp1` allocation failed"
        info = 1
        return
    end if
    allocate(constants % vscales_rel(constants % nscales), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`vscales_rel` allocation failed"
        info = 1
        return
    end if
    ! Compute scaling factors of spherical harmonics
    call ylmscale(constants % dmax, constants % vscales, &
        & constants % v4pi2lp1, constants % vscales_rel)
    ! Allocate square roots of factorials
    allocate(constants % vfact(constants % nfact), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`vfact` allocation failed"
        info = 1
        return
    end if
    ! Compute square roots of factorials
    constants % vfact(1) = 1
    do i = 2, constants % nfact
        constants % vfact(i) = constants % vfact(i-1) * sqrt(dble(i-1))
    end do
    ! Allocate square roots of combinatorial numbers C_n^k and M2L OZ
    ! translation coefficients
    if (params % fmm .eq. 1) then
        alloc_size = 2*constants % dmax + 1
        alloc_size = alloc_size * (constants % dmax+1)
        ! Allocate C_n^k
        allocate(constants % vcnk(alloc_size), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "`vcnk` allocation failed"
            info = 1
            return
        end if
        ! Allocate M2L OZ coefficients
        allocate(constants % m2l_ztranslate_coef(&
            & (params % pm+1), (params % pl+1), (params % pl+1)), &
            & stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "`m2l_ztranslate_coef` " // &
                & "allocation failed"
            info = 1
            return
        end if
        ! Allocate adjoint M2L OZ coefficients
        allocate(constants % m2l_ztranslate_adj_coef(&
            & (params % pl+1), (params % pl+1), (params % pm+1)), &
            & stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "`m2l_ztranslate_adj_coef` " // &
                & "allocation failed"
            info = 1
            return
        end if
        ! Compute combinatorial numbers C_n^k and M2L OZ translate coefficients
        call fmm_constants(constants % dmax, params % pm, &
            & params % pl, constants % vcnk, &
            & constants % m2l_ztranslate_coef, &
            & constants % m2l_ztranslate_adj_coef)
    end if
    ! Allocate space for Lebedev grid coordinates and weights
    allocate(constants % cgrid(3, params % ngrid), &
        & constants % wgrid(params % ngrid), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`cgrid` and `wgrid` " // &
            & "allocations failed"
        info = 1
        return
    end if
    ! Get weights and coordinates of Lebedev grid points
    call llgrid(params % ngrid, constants % wgrid, constants % cgrid)
    ! Allocate space for values of non-weighted and weighted spherical
    ! harmonics and L2P at Lebedev grid points
    allocate(constants % vgrid(constants % vgrid_nbasis, params % ngrid), &
        & constants % vwgrid(constants % vgrid_nbasis, params % ngrid), &
        & constants % l2grid(constants % vgrid_nbasis, params % ngrid), &
        & stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`vgrid`, `wgrid` and " // &
            & "`l2grid` allocations failed"
        info = 1
        return
    end if
    allocate(vplm(constants % vgrid_nbasis), vcos(constants % vgrid_dmax+1), &
        & vsin(constants % vgrid_dmax+1), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`vplm`, `vcos` and " // &
            & "`vsin` allocations failed"
        info = 1
        return
    end if
    ! Compute non-weighted and weighted spherical harmonics and the single
    ! layer operator at Lebedev grid points
    do igrid = 1, params % ngrid
        call ylmbas(constants % cgrid(:, igrid), rho, ctheta, stheta, cphi, &
            & sphi, constants % vgrid_dmax, constants % vscales, &
            & constants % vgrid(:, igrid), vplm, vcos, vsin)
        constants % vwgrid(:, igrid) = constants % wgrid(igrid) * &
            & constants % vgrid(:, igrid)
        do l = 0, params % pl
            indl = l*l + l + 1
            constants % l2grid(indl-l:indl+l, igrid) = &
                & constants % vgrid(indl-l:indl+l, igrid) / &
                & constants % vscales(indl)**2
        end do
    end do
    deallocate(vplm, vcos, vsin, stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`vplm`, `vcos` and " // &
            & "`vsin` deallocations failed"
        info = 1
        return
    end if
    ! Clean error state of constants to proceed with geometry
    constants % error_flag = 0
    constants % error_message = ""
    ! Generate geometry-related constants
    call constants_geometry_init(params, constants, info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "error in constants_geometry_init:" // &
            & constants % error_message
        ! The only possible error is related to allocation-deallocation, so
        ! output info is 1
        info = 1
        return
    end if
end subroutine constants_init

!> Initialize geometry-related constants like list of neighbouring spheres
!!
!! This routine does not check error state of params or constants as it is
!! intended for a low-level use.
!!
!! @param[in] params: Object containing all inputs.
!! @param[inout] constants: Object containing all constants.
!! @param[out] info: flag of succesfull exit
!!      = 0: Succesfull exit
!!      = -1: params is in error state
!!      = 1: Allocation of memory failed.
subroutine constants_geometry_init(params, constants, info)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    !! Outputs
    type(ddx_constants_type), intent(out) :: constants
    integer, intent(out) :: info
    !! Local variables
    real(dp) :: swthr, v(3), maxv, ssqv, vv, r, t
    integer :: nngmax, i, lnl, isph, jsph, inear, igrid
    integer, allocatable :: tmp_nl(:)
    !! The code
    ! Upper bound of switch region. Defines intersection criterion for spheres
    swthr = one + (params % se+one)*params % eta/two
    ! Build list of neighbours in CSR format
    nngmax = 1
    allocate(constants % inl(params % nsph+1), &
        & constants % nl(params % nsph*nngmax), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`inl` and `nl` allocations failed"
        info = 1
        return
    end if
    i = 1
    lnl = 0
    do isph = 1, params % nsph
        constants % inl(isph) = lnl + 1
        do jsph = 1, params % nsph
            if (isph .ne. jsph) then
                v = params % csph(:, isph)-params % csph(:, jsph)
                maxv = max(abs(v(1)), abs(v(2)), abs(v(3)))
                ssqv = (v(1)/maxv)**2 + (v(2)/maxv)**2 + (v(3)/maxv)**2
                vv = maxv * sqrt(ssqv)
                ! Take regularization parameter into account with respect to
                ! shift se. It is described properly by the upper bound of a
                ! switch region `swthr`.
                r = params % rsph(isph) + swthr*params % rsph(jsph)
                if (vv .le. r) then
                    constants % nl(i) = jsph
                    i  = i + 1
                    lnl = lnl + 1
                    ! Extend ddx_data % nl if needed
                    if (i .gt. params % nsph*nngmax) then
                        allocate(tmp_nl(params % nsph*nngmax), stat=info)
                        if (info .ne. 0) then
                            constants % error_flag = 1
                            constants % error_message = "`tmp_nl` " // &
                                & "allocation failed"
                            info = 1
                            return
                        end if
                        tmp_nl(1:params % nsph*nngmax) = &
                            & constants % nl(1:params % nsph*nngmax)
                        deallocate(constants % nl, stat=info)
                        if (info .ne. 0) then
                            constants % error_flag = 1
                            constants % error_message = "`nl` " // &
                                & "deallocation failed"
                            info = 1
                            return
                        end if
                        nngmax = nngmax + 10
                        allocate(constants % nl(params % nsph*nngmax), &
                            & stat=info)
                        if (info .ne. 0) then
                            constants % error_flag = 1
                            constants % error_message = "`nl` " // &
                                & "allocation failed"
                            info = 1
                            return
                        end if
                        constants % nl(1:params % nsph*(nngmax-10)) = &
                            & tmp_nl(1:params % nsph*(nngmax-10))
                        deallocate(tmp_nl, stat=info)
                        if (info .ne. 0) then
                            constants % error_flag = 1
                            constants % error_message = "`tmp_nl` " // &
                                & "deallocation failed"
                            info = 1
                            return
                        end if
                    end if
                end if
            end if
        end do
    end do
    constants % inl(params % nsph+1) = lnl+1
    constants % nngmax = nngmax
    ! Allocate space for characteristic functions fi and ui
    allocate(constants % fi(params % ngrid, params % nsph), &
        & constants % ui(params % ngrid, params % nsph), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`fi` and `ui` allocations failed"
        info = 1
        return
    end if
    constants % fi = zero
    constants % ui = zero
    ! Allocate space for force-related arrays
    if (params % force .eq. 1) then
        allocate(constants % zi(3, params % ngrid, params % nsph), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "`zi` allocation failed"
            info = 1
            return
        endif
        constants % zi = zero
    end if
    ! Build arrays fi, ui, zi
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            ! Loop over neighbours of i-th sphere
            do inear = constants % inl(isph), constants % inl(isph+1)-1
                ! Neighbour's index
                jsph = constants % nl(inear)
                ! Compute t_n^ij
                v = params % csph(:, isph) - params % csph(:, jsph) + &
                    & params % rsph(isph)*constants % cgrid(:, igrid)
                maxv = max(abs(v(1)), abs(v(2)), abs(v(3)))
                ssqv = (v(1)/maxv)**2 + (v(2)/maxv)**2 + (v(3)/maxv)**2
                vv = maxv * sqrt(ssqv)
                t = vv / params % rsph(jsph)
                ! Accumulate characteristic function \chi(t_n^ij)
                constants % fi(igrid, isph) = constants % fi(igrid, isph) + &
                    & fsw(t, params % se, params % eta)
                ! Check if gradients are required
                if (params % force .eq. 1) then
                    ! Check if t_n^ij belongs to switch region
                    if ((t .lt. swthr) .and. (t .gt. swthr-params % eta)) then
                        ! Accumulate gradient of characteristic function \chi
                        vv = dfsw(t, params % se, params % eta) / &
                            & params % rsph(jsph) / vv
                        constants % zi(:, igrid, isph) = &
                            & constants % zi(:, igrid, isph) + vv*v
                    end if
                end if
            enddo
            ! Compute characteristic function of a molecular surface ui
            if (constants % fi(igrid, isph) .le. one) then
                constants % ui(igrid, isph) = one - constants % fi(igrid, isph)
            end if
        enddo
    enddo
    ! Build cavity array. At first get total count for each sphere
    allocate(constants % ncav_sph(params % nsph), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`ncav_sph` allocation failed"
        info = 1
        return
    endif
    do isph = 1, params % nsph
        constants % ncav_sph(isph) = 0
        ! Loop over integration points
        do i = 1, params % ngrid
            ! Positive contribution from integration point
            if (constants % ui(i, isph) .gt. zero) then
                constants % ncav_sph(isph) = constants % ncav_sph(isph) + 1
            end if
        end do
    end do
    constants % ncav = sum(constants % ncav_sph)
    ! Allocate cavity array and CSR format for indexes of cavities
    allocate(constants % ccav(3, constants % ncav), &
        & constants % icav_ia(params % nsph+1), &
        & constants % icav_ja(constants % ncav), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`ccav`, `icav_ia` and " // &
            & "`icav_ja` allocations failed"
        info = 1
        return
    endif
    ! Get actual cavity coordinates and indexes in CSR format
    constants % icav_ia(1) = 1
    i = 1
    do isph = 1, params % nsph
        constants % icav_ia(isph+1) = constants % icav_ia(isph) + &
            & constants % ncav_sph(isph)
        ! Loop over integration points
        do igrid = 1, params % ngrid
            ! Ignore zero contribution
            if (constants % ui(igrid, isph) .eq. zero) cycle
            ! Store coordinates
            constants % ccav(:, i) = params % csph(:, isph) + &
                & params % rsph(isph)*constants % cgrid(:, igrid)
            ! Store index
            constants % icav_ja(i) = igrid
            ! Advance cavity array index
            i = i + 1
        end do
    end do
    ! RX_PRC
    ! ORDER
    ! NCLUSTERS
    ! CLUSTER
    ! CHILDREN
    ! PARENT
    ! CNODE
    ! SNODE
    ! NNFAR
    ! NNNEAR
    ! NFAR
    ! NNEAR
    ! FAR
    ! NEAR
    ! SFAR
    ! SNEAR
    ! NNEAR_M2P
end subroutine constants_geometry_init

!> Update geometry-related constants like list of neighbouring spheres
!!
!! This procedure simply deletes current geometry constants and initializes new
!! ones.
!!
!! @param[in] params: Object containing all inputs.
!! @param[inout] constants: Object containing all constants.
!! @param[out] info: flag of succesfull exit
!!      = 0: Succesfull exit
!!      = -1: params is in error state
!!      = 1: Allocation of memory failed.
subroutine constants_geometry_update(params, constants, info)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    !! Outputs
    type(ddx_constants_type), intent(out) :: constants
    integer, intent(out) :: info
    !! Local variables
    !! The code
    ! Enforce error for now (not yet implemented)
    stop 1
end subroutine constants_geometry_update

!> Switching function
!!
!! This is an implementation of \f$ \chi(t) \f$ with a shift \f$ se \f$:
!! \f[
!!      \chi(t) = \left\{ \begin{array}{ll} 0, & \text{if} \quad x \geq
!!      1 \\ p_\eta(x), & \text{if} \quad 1-\eta < x < 1 \\ 1, & \text{if}
!!      \quad x \leq 1-\eta \end{array} \right.
!! \f]
!! where \f$ x = t - \frac{1+se}{2} \eta \f$ is a shifted coordinate and
!! \f[
!!      p_\eta(x) = \frac{1}{\eta^5} (1-t)^3 (6t^2 + (15\eta-12)t + (10\eta^2
!!      -15\eta+6))
!! \f]
!! is a smoothing polynomial of the 5th degree.
!! In the case shift \f$ se=1 \f$ the switch function \f$ \chi(t) \f$ is 1 for
!! \f$ t \in [0,1] \f$, is 0 for \f$ t \in [1+\eta, \infty) \f$ and varies in
!! \f$ [1, 1+\eta] \f$ which is an external shift.
!! In the case shift \f$ se=-1 \f$ the switch function \f$ \chi(t) \f$ is 1 for
!! \f$ t \in [0,1-\eta] \f$, is 0 for \f$ t \in [1, \infty) \f$ and varies in
!! \f$ [1-\eta, 1] \f$ which is an internal shift.
!! In the case shift \f$ se=0 \f$ the switch function \f$ \chi(t) \f$ is 1 for
!! \f$ t \in [0,1-\eta/2] \f$, is 0 for \f$ t \in [1+\eta/2, \infty) \f$ and
!! varies in \f$ [1-\eta/2, 1+\eta/2] \f$ which is a centered shift.
!!
!! @param[in] t: real non-negative input value
!! @param[in] se: shift
!! @param[in] eta: regularization parameter \f$ \eta \f$
real(dp) function fsw(t, se, eta)
    ! Inputs
    real(dp), intent(in) :: t, se, eta
    ! Local variables
    real(dp) :: a, b, c, x
    ! Apply shift:
    !   se =  0   =>   t - eta/2  [ CENTERED ]
    !   se =  1   =>   t - eta    [ EXTERIOR ]
    !   se = -1   =>   t          [ INTERIOR ]
    x = t - (se*pt5 + pt5)*eta
    ! a must be in range (0,eta) for the switch region
    a = one - x
    ! Define switch function chi
    ! a <= 0 corresponds to the outside region
    if (a .le. zero) then
        fsw = zero
    ! a >= eta corresponds to the inside region
    else if (a .ge. eta) then
        fsw = one
    ! Switch region wher x is in range (1-eta,1)
    else
        ! Normalize a to region (0,1)
        a = a / eta
        b = 6d0*a - 15
        c = b*a + 10
        fsw = a * a * a * c
    end if
end function fsw

!> Derivative of a switching function
!!
!! This is an implementation of \f$ \chi'(t) \f$ with a shift \f$ se \f$:
!! \f[
!!      \chi'(t) = \left\{ \begin{array}{ll} 0, & \text{if} \quad x \geq
!!      1 \\ p'_\eta(x), & \text{if} \quad 1-\eta < x < 1 \\ 0, & \text{if}
!!      \quad x \leq 1-\eta \end{array} \right.
!! \f]
!! where \f$ x = t - \frac{1+se}{2} \eta \f$ is a shifted coordinate and
!! \f[
!!      p'_\eta(x) = -\frac{30}{\eta^5} (1-t)^2 (t-1+\eta)^2
!! \f]
!! is a derivative of the smoothing polynomial.
!! In the case shift \f$ se=1 \f$ the derivative \f$ \chi'(t) \f$ is 0 for
!! \f$ t \in [0,1] \cup [1+\eta, \infty) \f$ and varies in
!! \f$ [1, 1+\eta] \f$ which is an external shift.
!! In the case shift \f$ se=-1 \f$ the derivative \f$ \chi'(t) \f$ is 0 for
!! \f$ t \in [0,1-\eta] \cup [1, \infty) \f$ and varies in
!! \f$ [1-\eta, 1] \f$ which is an internal shift.
!! In the case shift \f$ se=0 \f$ the derivative \f$ \chi'(t) \f$ is 0 for
!! \f$ t \in [0,1-\eta/2] \cup [1+\eta/2, \infty) \f$ and
!! varies in \f$ [1-\eta/2, 1+\eta/2] \f$ which is a centered shift.
!!
!! @param[in] t: real non-negative input value
!! @param[in] se: shift
!! @param[in] eta: regularization parameter \f$ \eta \f$
real(dp) function dfsw(t, se, eta)
    ! Inputs
    real(dp), intent(in) :: t, se, eta
    ! Local variables
    real(dp) :: flow, x
    real(dp), parameter :: f30=30.0d0
    ! Apply shift:
    !   s =  0   =>   t - eta/2  [ CENTERED ]
    !   s =  1   =>   t - eta    [ EXTERIOR ]
    !   s = -1   =>   t          [ INTERIOR ]
    x = t - (se + 1.d0)*eta / 2.d0
    ! Lower bound of switch region
    flow = one - eta
    ! Define derivative of switch function chi
    if (x .ge. one) then
        dfsw = zero
    else if (x .le. flow) then
        dfsw = zero
    else
        dfsw = -f30 * (( (one-x)*(x-one+eta) )**2) / (eta**5)
    endif
end function dfsw

!> Compute scaling factors of real normalized spherical harmonics
!!
!! Output values of scaling factors of \f$ Y_\ell^m \f$ harmonics are filled
!! only for non-negative \f$ m \f$ since scaling factor of \f$ Y_\ell^{-m} \f$
!! is the same as scaling factor of \f$ Y_\ell^m \f$.
!!
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[out] vscales: Array of scaling factors. Dimension is `(p+1)**2`
!! @param[out] vscales: Array of values 4pi/(2l+1). Dimension is `p+1`
!! @param[out] vscales_rel: Array of relative scaling factors.
!!      Dimension is `(p+1)**2`.
subroutine ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    ! Input
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2)
    ! Local variables
    real(dp) :: tmp, twolp1
    integer :: l, ind, m
    twolp1 = one
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        tmp = fourpi / twolp1
        v4pi2lp1(l+1) = tmp
        tmp = sqrt(tmp)
        vscales_rel(ind) = tmp
        vscales(ind) = one / tmp
        twolp1 = twolp1 + two
        tmp = vscales(ind) * sqrt2
        ! m != 0
        do m = 1, l
            tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
            vscales(ind+m) = tmp
            vscales(ind-m) = tmp
            vscales_rel(ind+m) = tmp * v4pi2lp1(l+1)
            vscales_rel(ind-m) = vscales_rel(ind+m)
        end do
    end do
end subroutine ylmscale

!> Compute FMM-related constants
!!
!! @param[in] dmax: Maximal degree of spherical harmonics to be evaluated.
!!      `dmax` >= 0
!! @param[in] pm: Maximal degree of the multipole expansion. `pm` >= 0.
!! @param[in] pl: Maximal degree of the local expansion. `pl` >= 0.
!! @param[out] vcnk: Array of squre roots of combinatorial factors C_n^k.
!!      Dimension is `(2*dmax+1)*(dmax+1)`.
!! @param[out] m2l_ztranslate_coef: Constants for M2L translation over OZ axis.
!!      Dimension is `(pm+1, pl+1, pl+1)`.
!! @param[out] m2l_ztranslate_coef: Constants for adjoint M2L translation over
!!      OZ axis. Dimension is `(pl+1, pl+1, pm+1)`.
subroutine fmm_constants(dmax, pm, pl, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    ! Inputs
    integer, intent(in) :: dmax, pm, pl
    ! Outputs
    real(dp), intent(out) :: vcnk((2*dmax+1)*(dmax+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1)
    ! Local variables
    integer :: i, indi, j, k, n, indjn
    real(dp) :: tmp1
    ! Compute combinatorial numbers C_n^k for n=0..dmax
    ! C_0^0 = 1
    vcnk(1) = one
    do i = 2, 2*dmax+1
        ! Offset to the C_{i-2}^{i-2}, next item to be stored is C_{i-1}^0
        indi = (i-1) * i / 2
        ! C_{i-1}^0 = 1
        vcnk(indi+1) = one
        ! C_{i-1}^{i-1} = 1
        vcnk(indi+i) = one
        ! C_{i-1}^{j-1} = C_{i-2}^{j-1} + C_{i-2}^{j-2}
        ! Offset to C_{i-3}^{i-3} is indi-i+1
        do j = 2, i-1
            vcnk(indi+j) = vcnk(indi-i+j+1) + vcnk(indi-i+j)
        end do
    end do
    ! Get square roots of C_n^k. sqrt(one) is one, so no need to update C_n^0
    ! and C_n^n
    do i = 3, 2*dmax+1
        indi = (i-1) * i / 2
        do j = 2, i-1
            vcnk(indi+j) = sqrt(vcnk(indi+j))
        end do
    end do
    ! Fill in m2l_ztranslate_coef and m2l_ztranslate_adj_coef
    do j = 0, pl
        do k = 0, j
            tmp1 = one
            do n = k, pm
                indjn = (j+n)*(j+n+1)/2 + 1
                m2l_ztranslate_coef(n-k+1, k+1, j-k+1) = &
                    & tmp1 * vcnk(indjn+j-k) * vcnk(indjn+j+k)
                m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1) = &
                    & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)
                tmp1 = -tmp1
            end do
        end do
    end do
end subroutine fmm_constants

end module ddx_constants

