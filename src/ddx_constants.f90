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
use omp_lib, only : omp_get_wtime

implicit none

!> @defgroup Fortran_interface_core Fortran interface: core routines

!> Container for precomputed constants
type ddx_constants_type
    !> Number modeling spherical harmonics per sphere.
    integer :: nbasis
    !> Total number of modeling degrees of freedom.
    !! This is equal to `nsph*nbasis`.
    integer :: n
    !> Maximal degree of used real normalized spherical harmonics.
    !! For example, if FMM is
    !! used, its M2L operation requires computing spherical harmonics of all
    !! degrees up to `pm+pl`. If `force=1` then this parameter might be
    !! increased by 1 because certain implementations of analytical gradients
    !! might rely on spherical harmonics of +1 degree.
    integer :: dmax
    !> Total number of used real spherical harmonics and a size of `vscales`.
    integer :: nscales
    !> Scales of real normalized spherical harmonics of degree up to dmax.
    !! This array has a dimension (nscales).
    real(dp), allocatable :: vscales(:)
    !> Array of values 4pi/(2l+1), dimension is (dmax+1). Referenced only if
    !! fmm=1, but allocated and computed in any case.
    real(dp), allocatable :: v4pi2lp1(:)
    !> Relative scales of real normalized spherical harmonics.
    !! Each values is multiplied by a corresponding 4pi/(2l+1). Dimension of
    !! this array is (nscales). Referenced only if fmm=1, but allocated and
    !! computed in any case.
    real(dp), allocatable :: vscales_rel(:)
    !> Number of precomputed square roots of factorials.
    !! Just like with `dmax` parameter, number of used factorials is either
    !! `2*lmax+1` or `2*(pm+pl)+1` depending on whether FMM is used or not.
    integer :: nfact
    !> Array of square roots of factorials of a dimension (nfact).
    real(dp), allocatable :: vfact(:)
    !> Array of square roots of combinatorial numbers C_n^k.
    !! Dimension of this array is ((2*dmax+1)*(dmax+1)). Allocated, computed
    !! and referenced only if fmm=1 or if force=1
    real(dp), allocatable :: vcnk(:)
    !> Array of common M2L coefficients for any OZ translation.
    !! This array has a dimension (pm+1, pl+1, pl+1). Allocated, computed and
    !! referenced only if fmm=1.
    real(dp), allocatable :: m2l_ztranslate_coef(:, :, :)
    !> Array of common M2L coefficients for any adjoint OZ translation.
    !! Dimension of this array is (pl+1, pl+1, pm+1). It is allocated, computed
    !! and referenced only if fmm=1.
    real(dp), allocatable :: m2l_ztranslate_adj_coef(:, :, :)
    !> Coordinates of Lebedev quadrature points of a dimension (3, ngrid).
    real(dp), allocatable :: cgrid(:, :)
    !> Weights of Lebedev quadrature points of a dimension (ngrid).
    real(dp), allocatable :: wgrid(:)
    !> Maximal degree of spherical harmonics evaluated at Lebedev grid points.
    !! Although we use spherical harmonics of degree up to `dmax`, only
    !! spherical harmonics of degree up to `lmax` and `pl` are evaluated
    !! at Lebedev grid points. In the case `force=1` this degrees might be
    !! increased by 1 depending on implementation of gradients.
    integer :: vgrid_dmax
    !> Number of spherical harmonics evaluated at Lebedev grid points.
    integer :: vgrid_nbasis
    !> Values of spherical harmonics at Lebedev grid points.
    !! Dimensions of this array are (vgrid_nbasis, ngrid)
    real(dp), allocatable :: vgrid(:, :)
    !> Weighted values of spherical harmonics at Lebedev grid points.
    !! vwgrid(:, igrid) = vgrid(:, igrid) * wgrid(igrid)
    !! Dimension of this array is (vgrid_nbasis, ngrid).
    real(dp), allocatable :: vwgrid(:, :)
    !> Values of L2P/M2P at grid points. Dimension is (vgrid_nbasis, ngrid).
    real(dp), allocatable :: vgrid2(:, :)
    !> LPB value max of lmax and 6
    integer :: lmax0
    !> LPB value max of nbasis and 49
    integer :: nbasis0
    !> LPB matrix, Eq. (87) from [QSM19.SISC]. Dimension is
    !! (nbasis, nbasis0, nsph)
    real(dp), allocatable :: Pchi(:, :, :)
    !> LPB value (i'_l0(r_j)/i_l0(r_j)-k'_l0(r_j)/k_l0(r_j))^{-1}. Dimension
    !! is ???
    real(dp), allocatable :: C_ik(:, :)
    !> LPB Bessel function of the first kind. Dimension is (dmax+1, nsph).
    real(dp), allocatable :: SI_ri(:, :)
    !> LPB Derivative of Bessel function of the first kind. Dimension is
    !! (dmax+1, nsph).
    real(dp), allocatable :: DI_ri(:, :)
    !> LPB Bessel function of the second kind. Dimension is (lmax+2, nsph).
    real(dp), allocatable :: SK_ri(:, :)
    !> LPB Derivative Bessel function of the second kind. Dimension is
    !! (lmax+2, nsph).
    real(dp), allocatable :: DK_ri(:, :)
    !> LPB value i'_l(r_j)/i_l(r_j). Dimension is (lmax, nsph).
    real(dp), allocatable :: termimat(:, :)
    !> LPB B matrix for doing incore BX product
    real(dp), allocatable :: b(:,:,:)
    !> ddCOSMO L matrix fo doing incore LX product
    real(dp), allocatable :: l(:,:,:)
    !> Upper limit on a number of neighbours per sphere. This value is just an
    !! upper bound that is not guaranteed to be the actual maximum.
    integer :: nngmax
    !> List of intersecting spheres in a CSR format. Dimension is (nsph+1).
    integer, allocatable :: inl(:)
    !> List of intersecting spheres in a CSR format. Dimension is
    !! (nsph*nngmax).
    integer, allocatable :: nl(:)
    !> transpose list of intersecting spheres
    integer, allocatable :: itrnl(:)
    !> Values of a characteristic function f at all grid points of all spheres.
    !! Dimension is (ngrid, npsh).
    real(dp), allocatable :: fi(:, :)
    !> Values of a characteristic function U at all grid points of all spheres.
    !! Dimension is (ngrid, nsph).
    real(dp), allocatable :: ui(:, :)
    !> Values of a characteristic function U at cavity points. Dimension is
    !! (ncav).
    real(dp), allocatable :: ui_cav(:)
    !> Derivative of the characteristic function U at all grid points of all
    !! spheres. Dimension is (3, ngrid, nsph).
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
    !! the PCM model (model=3). Dimension is (nbasis, nbasis, nsph).
    real(dp), allocatable :: rx_prc(:, :, :)
    !! Cluster tree information that is allocated and computed only if fmm=1.
    !> Reordering of spheres for better locality. This array has a dimension
    !! (nsph) and is allocated/used only if fmm=1.
    integer, allocatable :: order(:)
    !> Number of clusters. Defined only if fmm=1.
    integer :: nclusters
    !> The first and the last spheres of each node. Dimension of this array is
    !! (2, nclusters) and it is allocated/used only if fmm=1.
    integer, allocatable :: cluster(:, :)
    !> Children of each cluster. Dimension is (2, nclusters). Allocated and
    !! used only if fmm=1.
    integer, allocatable :: children(:, :)
    !> Parent of each cluster. Dimension is (nclusters). Allocated and used
    !! only if fmm=1.
    integer, allocatable :: parent(:)
    !> Center of bounding sphere of each cluster. Dimension is (3, nclusters).
    !! This array is allocated and used only if fmm=1.
    real(dp), allocatable :: cnode(:, :)
    !> Radius of bounding sphere of each cluster. Dimension is (nclusters).
    !! This array is allocated and used only if fmm=1.
    real(dp), allocatable :: rnode(:)
    !> Which leaf node contains only given input sphere. Dimension is (nsph).
    !! This array is allocated and used only if fmm=1.
    integer, allocatable :: snode(:)
    !> Bessel function for bounding sphere of each cluster. Dimension is
    !! (pm, nclusters). This array is allocated and used only if fmm=1.
    real(dp), allocatable :: SK_rnode(:, :)
    !> Bessel function for bounding sphere of each cluster. Dimension is
    !! (pm, nclusters). This array is allocated and used only if fmm=1.
    real(dp), allocatable :: SI_rnode(:, :)
    !> Total number of far admissible pairs. Defined only if fmm=1.
    integer :: nnfar
    !> Total number of near admissible pairs. Defined only if fmm=1.
    integer :: nnnear
    !> Number of admissible far pairs for each node. Dimension is (nclusters).
    !! This array is allocated and used only if fmm=1.
    integer, allocatable :: nfar(:)
    !> Number of admissible near pairs for each node. Dimension is (nclusters).
    !! This array is allocated and used only if fmm=1.
    integer, allocatable :: nnear(:)
    !> Arrays of admissible far pairs. Dimension is (nnfar). This array is
    !! allocated and used only if fmm=1.
    integer, allocatable :: far(:)
    !> Arrays of admissible near pairs. Dimension is (nnnear). This array is
    !! allocated and used only if fmm=1.
    integer, allocatable :: near(:)
    !> Index of the first element of array of all admissible far pairs stored
    !! in the array `far`. Dimension is (nclusters+1). This array is
    !! allocated and used only if fmm=1.
    integer, allocatable :: sfar(:)
    !> Index of the first element of array of all admissible near pairs stored
    !! in the array `near`. Dimension is (nclusters+1). This array is
    !! allocated and used only if fmm=1.
    integer, allocatable :: snear(:)
    !> Number of near-field M2P interactions with cavity points. Defined only
    !! if fmm=1.
    integer :: nnear_m2p
    !> Maximal degree of near-field M2P spherical harmonics. Defined only if
    !! fmm=1.
    integer :: m2p_lmax
    !> Number of spherical harmonics used for near-field M2P. Defined only if
    !! fmm=1.
    integer :: m2p_nbasis
    !> Number of spherical harmonics of degree up to lmax+1 used for
    !! computation of forces (gradients). Allocated and used only if
    !! fmm=1.
    integer :: grad_nbasis
    !> Inner tolerance for microiterations done when using ddLPB
    real(dp) :: inner_tol
    !> Whether the diagonal of the matrices has to be used in the mvp for
    !! ddCOSMO, ddPCM or inner ddLPB iterations
    logical  :: dodiag
    !> Flag if there were an error
    integer :: error_flag = 2
    !> Last error message
    character(len=255) :: error_message
end type ddx_constants_type

contains

!> Compute all necessary constants
!> @ingroup Fortran_interface_core
!!
!! @param[in] params: Object containing all inputs.
!! @param[out] constants: Object containing all constants.
!!
subroutine constants_init(params, constants)
    use complex_bessel
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    !! Outputs
    type(ddx_constants_type), intent(out) :: constants
    !! Local variables
    integer :: i, alloc_size, l, indl, igrid, isph, ind, icav, l0, m0, ind0, &
        & jsph, ibasis, ibasis0, NZ, ierr, info
    real(dp) :: rho, ctheta, stheta, cphi, sphi, termi, termk, term, rijn, &
        & sijn(3), vij(3), val, s1, s2
    real(dp), allocatable :: vplm(:), vcos(:), vsin(:), vylm(:), SK_rijn(:), &
        & DK_rijn(:)
    complex(dp), allocatable :: bessel_work(:)
    complex(dp) :: z
    !! The code
    ! Clean error state of constants to proceed with geometry
    constants % error_flag = 0
    constants % error_message = ""
    ! Check if params are OK
    if (params % error_flag .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "constants_init: `params` is in " // &
            & "error state"
        return
    end if
    ! activate inner iterations diagonal in mvp for debugging purposes only.
    ! could be useful for different linear solvers.
    constants % dodiag = .false.
    ! Maximal number of modeling spherical harmonics
    constants % nbasis = (params % lmax+1) ** 2
    ! Maximal number of modeling degrees of freedom
    constants % n = params % nsph * constants % nbasis
    ! Calculate dmax, vgrid_dmax, m2p_lmax, m2p_nbasis and grad_nbasis
    if (params % fmm .eq. 0) then
        if (params % force .eq. 1) then
            constants % dmax = params % lmax + 1
            constants % vgrid_dmax = params % lmax + 1
        else
            constants % dmax = params % lmax
            constants % vgrid_dmax = params % lmax
        end if
        ! Other constants are not referenced if fmm=0
        constants % m2p_lmax = -1
        constants % m2p_nbasis = -1
        constants % grad_nbasis = -1
    else
        ! If forces are required then we need the M2P of a degree lmax+1 for
        ! the near-field analytical gradients
        if (params % force .eq. 1) then
            constants % dmax = max(params % pm+params % pl+1, &
                & params % lmax+1)
            constants % m2p_lmax = params % lmax + 1
            constants % grad_nbasis = (params % lmax+2) ** 2
            constants % vgrid_dmax = max(params % pl, params % lmax) + 1
        else
            constants % dmax = max(params % pm+params % pl, &
                & params % lmax)
            constants % m2p_lmax = params % lmax
            constants % grad_nbasis = -1
            constants % vgrid_dmax = max(params % pl, params % lmax)
        end if
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
        constants % error_message = "constants_init: `vscales` allocation " &
            & // "failed"
        return
    end if
    allocate(constants % v4pi2lp1(constants % dmax+1), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "constants_init: `v4pi2lp1` allocation " &
            & // "failed"
        return
    end if
    allocate(constants % vscales_rel(constants % nscales), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "constants_init: `vscales_rel` " // &
            & "allocation failed"
        return
    end if
    ! Compute scaling factors of spherical harmonics
    call ylmscale(constants % dmax, constants % vscales, &
        & constants % v4pi2lp1, constants % vscales_rel)
    ! Allocate square roots of factorials
    allocate(constants % vfact(constants % nfact), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "constants_init: `vfact` allocation failed"
        return
    end if
    ! Compute square roots of factorials
    constants % vfact(1) = 1
    do i = 2, constants % nfact
        constants % vfact(i) = constants % vfact(i-1) * sqrt(dble(i-1))
    end do
    ! Allocate square roots of combinatorial numbers C_n^k and M2L OZ
    ! translation coefficients
    if (params % fmm .eq. 1 .or. &
        & (params % model .eq. 3 .and. params % force .eq. 1)) then
        alloc_size = 2*constants % dmax + 1
        alloc_size = alloc_size * (constants % dmax+1)
        ! Allocate C_n^k
        allocate(constants % vcnk(alloc_size), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_init: `vcnk` allocation " &
                & // "failed"
            return
        end if
        ! Allocate M2L OZ coefficients
        allocate(constants % m2l_ztranslate_coef(&
            & (params % pm+1), (params % pl+1), (params % pl+1)), &
            & stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_init: " // &
                & "`m2l_ztranslate_coef` allocation failed"
            return
        end if
        ! Allocate adjoint M2L OZ coefficients
        allocate(constants % m2l_ztranslate_adj_coef(&
            & (params % pl+1), (params % pl+1), (params % pm+1)), &
            & stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_init: " // &
                & "`m2l_ztranslate_adj_coef` allocation failed"
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
        constants % error_message = "constants_init: `cgrid` and `wgrid` " // &
            & "allocations failed"
        return
    end if
    ! Get weights and coordinates of Lebedev grid points
    call llgrid(params % ngrid, constants % wgrid, constants % cgrid)
    ! Allocate space for values of non-weighted and weighted spherical
    ! harmonics and L2P at Lebedev grid points
    allocate(constants % vgrid(constants % vgrid_nbasis, params % ngrid), &
        & constants % vwgrid(constants % vgrid_nbasis, params % ngrid), &
        & stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "constants_init: `vgrid`, `wgrid` and" // &
            & " allocations failed"
        return
    end if
    allocate(vplm(constants % vgrid_nbasis), vcos(constants % vgrid_dmax+1), &
        & vsin(constants % vgrid_dmax+1), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "constants_init: `vplm`, `vcos` and " // &
            & "`vsin` allocations failed"
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
    end do
    if (params % fmm .eq. 1) then
        allocate(constants % vgrid2( &
            & constants % vgrid_nbasis, params % ngrid), &
            & stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_init: `vgrid2` " // &
                & "allocation failed"
            return
        end if
        do igrid = 1, params % ngrid
            do l = 0, constants % vgrid_dmax
                indl = l*l + l + 1
                constants % vgrid2(indl-l:indl+l, igrid) = &
                    & constants % vgrid(indl-l:indl+l, igrid) / &
                    & constants % vscales(indl)**2
            end do
        end do
    end if
    ! Generate geometry-related constants (required by the LPB code)
    call constants_geometry_init(params, constants)
    if (constants % error_flag .ne. 0) return
    ! Precompute LPB-related constants
    if (params % model .eq. 3) then
        constants % lmax0 = params % lmax
        constants % nbasis0 = constants % nbasis
        allocate(vylm(constants % vgrid_nbasis), &
            & SK_rijn(0:constants % lmax0), DK_rijn(0:constants % lmax0), &
            & bessel_work(constants % dmax+2), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_init: `vylm`, `SK_rijn` and " // &
                & "`DK_rijn` allocations failed"
            return
        end if
        allocate(constants % SI_ri(0:constants % dmax+1, params % nsph))
        allocate(constants % DI_ri(0:constants % dmax+1, params % nsph))
        allocate(constants % SK_ri(0:params % lmax+1, params % nsph))
        allocate(constants % DK_ri(0:params % lmax+1, params % nsph))
        allocate(constants % Pchi(constants % nbasis, constants % nbasis0, &
            & params % nsph))
        allocate(constants % C_ik(0:params % lmax, params % nsph))
        allocate(constants % termimat(0:params % lmax, params % nsph))
        SK_rijn = zero
        DK_rijn = zero
        do isph = 1, params % nsph
            ! We compute Bessel functions of degrees 0..lmax+1 because the
            ! largest degree is required for forces
            call modified_spherical_bessel_first_kind(constants % dmax+1, &
                & params % rsph(isph)*params % kappa,&
                & constants % SI_ri(:, isph), constants % DI_ri(:, isph), &
                & bessel_work)
            call modified_spherical_bessel_second_kind(params % lmax+1, &
                & params % rsph(isph)*params % kappa, &
                & constants % SK_ri(:, isph), constants % DK_ri(:, isph), &
                & bessel_work)
            ! Compute matrix PU_i^e(x_in)
            ! Previous implementation in update_rhs. Made it in ddinit, so as to use
            ! it in Forces as well.
            call mkpmat(params, constants, isph, constants % Pchi(:, :, isph))
            ! Compute i'_l(r_i)/i_l(r_i)
            do l = 0, params % lmax
                constants % termimat(l, isph) = constants % DI_ri(l, isph) / &
                    & constants % SI_ri(l, isph) * params % kappa
            end do
            ! Compute (i'_l0/i_l0 - k'_l0/k_l0)^(-1) is computed in Eq.(97)
            do l0 = 0, constants % lmax0
                termi = constants % DI_ri(l0, isph) / &
                    & constants % SI_ri(l0, isph) * params % kappa
                termk = constants % DK_ri(l0, isph)/ &
                    & constants % SK_ri(l0, isph) * params % kappa
                constants % C_ik(l0, isph) = one / (termi-termk)
            end do
        end do
        if (params % fmm .eq. 1) then
            allocate(constants % SI_rnode(params % pm+1, constants % nclusters))
            allocate(constants % SK_rnode(params % pm+1, constants % nclusters))
            do i = 1, constants % nclusters
                z = constants % rnode(i) * params % kappa
                s1 = sqrt(two/(pi*real(z)))
                s2 = sqrt(pi/(two*real(z)))
                call cbesk(z, pt5, 1, 1, bessel_work(1), NZ, ierr)
                constants % SK_rnode(1, i) = s1 * real(bessel_work(1))
                call cbesi(z, pt5, 1, 1, bessel_work(1), NZ, ierr)
                constants % SI_rnode(1, i) = s2 * real(bessel_work(1))
                if (params % pm .gt. 0) then
                    call cbesk(z, 1.5d0, 1, params % pm, bessel_work(2:), NZ, ierr)
                    constants % SK_rnode(2:, i) = s1 * real(bessel_work(2:params % pm+1))
                    call cbesi(z, 1.5d0, 1, params % pm, bessel_work(2:), NZ, ierr)
                    constants % SI_rnode(2:, i) = s2 * real(bessel_work(2:params % pm+1))
                end if
            end do
        end if
        deallocate(vylm, SK_rijn, DK_rijn, bessel_work, stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_init: `vylm`, `SK_rijn` and " // &
                & "`DK_rijn` deallocations failed"
            return
        end if
        ! if doing incore build nonzero blocks of B
        if (params % matvecmem .eq. 1) then
            call build_b(constants, params)
        end if
    end if
    deallocate(vplm, vcos, vsin, stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "constants_init: `vplm`, `vcos` and " // &
            & "`vsin` deallocations failed"
        return
    end if
    ! if doing incore build nonzero blocks of L
    if (params % matvecmem .eq. 1) then
        call build_itrnl(constants, params)
        call build_l(constants, params)
    end if
end subroutine constants_init

!> Build the transposed neighbor list
subroutine build_itrnl(constants, params)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    integer :: isph, ij, jsph, ji, istat

    allocate(constants % itrnl(constants % inl(params % nsph + 1)), stat=istat)
    if (istat.ne.0) then
        constants % error_message = 'Allocation failed in build_itrnl'
        constants % error_flag = 1
        return
    end if

    do isph = 1, params % nsph
        do ij = constants % inl(isph), constants % inl(isph + 1) - 1
            jsph = constants % nl(ij)
            do ji = constants % inl(jsph), constants % inl(jsph + 1) - 1
                if (constants % nl(ji) .eq. isph) exit
            end do
            constants % itrnl(ij) = ji
        end do
    end do
end subroutine build_itrnl

!> Allocate and build the ddCOSMO sparse matrix, only if incore is set
subroutine build_l(constants, params)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    integer :: isph, ij, jsph, igrid, l, m, ind, info
    real(dp), dimension(3) :: vij, sij
    real(dp) :: vvij, tij, xij, oij, rho, ctheta, stheta, cphi, sphi, &
        & fac, tt, thigh
    real(dp), dimension(constants % nbasis) :: vylm, vplm
    real(dp), dimension(params % lmax + 1) :: vcos, vsin
    real(dp), dimension(constants % nbasis, params % ngrid) :: scratch
    real(dp) :: t

    allocate(constants % l(constants % nbasis, constants % nbasis, &
        & constants % inl(params % nsph + 1)), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = 'Allocation failed in build_l'
        return
    end if

    thigh = one + pt5*(params % se + one)*params % eta

    t = omp_get_wtime()
    !$omp parallel do default(none) shared(params,constants,thigh) &
    !$omp private(isph,ij,jsph,scratch,igrid,vij,vvij,tij,sij,xij,oij, &
    !$omp rho,ctheta,stheta,cphi,sphi,vylm,vplm,vcos,vsin,l,fac,ind,m,tt)
    do isph = 1, params % nsph
        do ij = constants % inl(isph), constants % inl(isph + 1) - 1
            jsph = constants % nl(ij)
            scratch = zero
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, isph).eq.one) cycle
                vij = params % csph(:, isph) &
                    & + params % rsph(isph)*constants % cgrid(:,igrid) &
                    & - params % csph(:, jsph)
                vvij = sqrt(dot_product(vij, vij))
                tij = vvij/params % rsph(jsph)
                if (tij.lt.thigh .and. tij.gt.zero) then
                    sij = vij/vvij
                    xij = fsw(tij, params % se, params % eta)
                    if (constants % fi(igrid, isph).gt.one) then
                        oij = xij/constants % fi(igrid, isph)
                    else
                        oij = xij
                    end if
                    call ylmbas(sij, rho, ctheta, stheta, cphi, sphi, &
                        & params % lmax, constants % vscales, vylm, vplm, &
                        & vcos, vsin)
                    tt = oij
                    do l = 0, params % lmax
                        ind = l*l + l + 1
                        fac = - tt/(constants % vscales(ind)**2)
                        do m = -l, l
                            scratch(ind + m, igrid) = fac*vylm(ind + m)
                        end do
                        tt = tt*tij
                    end do
                end if
            end do
            call dgemm('n', 't', constants % nbasis, constants % nbasis, params % ngrid, &
                & one, constants % vwgrid, constants % vgrid_nbasis, scratch, &
                & constants % nbasis, zero, constants % l(:,:,ij), constants % nbasis)
        end do
    end do
end subroutine build_l

!> Allocate and build the HSP sparse matrix, only if incore is set
subroutine build_b(constants, params)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    integer :: isph, ij, jsph, igrid, l, m, ind
    real(dp), dimension(3) :: vij, sij
    real(dp) :: vvij, tij, xij, oij, rho, ctheta, stheta, cphi, sphi, &
        & fac, vvtij, thigh
    real(dp), dimension(constants % nbasis) :: vylm, vplm
    real(dp), dimension(params % lmax + 1) :: vcos, vsin
    complex(dp), dimension(max(2, params % lmax + 1)) :: bessel_work
    real(dp), dimension(0:params % lmax) :: SI_rijn, DI_rijn
    real(dp), dimension(constants % nbasis, params % ngrid) :: scratch
    real(dp) :: t
    integer :: info

    allocate(constants % b(constants % nbasis, constants % nbasis, &
        & constants % inl(params % nsph + 1)), stat=info)

    thigh = one + pt5*(params % se + one)*params % eta

    t = omp_get_wtime()
    !$omp parallel do default(none) shared(params,constants,thigh) &
    !$omp private(isph,ij,jsph,scratch,igrid,vij,vvij,tij,sij,xij,oij, &
    !$omp rho,ctheta,stheta,cphi,sphi,vylm,vplm,vcos,vsin,si_rijn,di_rijn, &
    !$omp vvtij,l,fac,ind,m,bessel_work)
    do isph = 1, params % nsph
        do ij = constants % inl(isph), constants % inl(isph + 1) - 1
            jsph = constants % nl(ij)
            scratch = zero
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, isph).eq.one) cycle
                vij = params % csph(:, isph) &
                    & + params % rsph(isph)*constants % cgrid(:,igrid) &
                    & - params % csph(:, jsph)
                vvij = sqrt(dot_product(vij, vij))
                tij = vvij/params % rsph(jsph)
                if (tij.lt.thigh .and. tij.gt.zero) then
                    sij = vij/vvij
                    xij = fsw(tij, params % se, params % eta)
                    if (constants % fi(igrid, isph).gt.one) then
                        oij = xij/constants % fi(igrid, isph)
                    else
                        oij = xij
                    end if
                    call ylmbas(sij, rho, ctheta, stheta, cphi, sphi, &
                        & params % lmax, constants % vscales, vylm, vplm, &
                        & vcos, vsin)
                    SI_rijn = 0
                    DI_rijn = 0
                    vvtij = vvij*params % kappa
                    call modified_spherical_bessel_first_kind(params % lmax, &
                        & vvtij, SI_rijn, DI_rijn, bessel_work)
                    do l = 0, params % lmax
                        fac = - oij*SI_rijn(l)/constants % SI_ri(l, jsph)
                        ind = l*l + l + 1
                        do m = -l, l
                            scratch(ind + m, igrid) = fac*vylm(ind + m)
                        end do
                    end do
                end if
            end do
            call dgemm('n', 't', constants % nbasis, constants % nbasis, params % ngrid, &
                & one, constants % vwgrid, constants % vgrid_nbasis, scratch, &
                & constants % nbasis, zero, constants % b(:,:,ij), constants % nbasis)
        end do
    end do
end subroutine build_b

!> Computation of P_chi
!!
!! @param[in]  isph : Sphere number
!! @param[out] pmat : Matrix of size nbasis X (lmax0+1)^2, Fixed lmax0
!!
subroutine mkpmat(params, constants, isph, pmat)
    type(ddx_params_type), intent(in)  :: params
    type(ddx_constants_type), intent(in)  :: constants
    integer,  intent(in) :: isph
    real(dp), dimension(constants % nbasis, constants % nbasis0), intent(inout) :: pmat
    integer :: l, m, ind, l0, m0, ind0, its
    real(dp)  :: f, f0
    pmat(:,:) = zero
    do its = 1, params % ngrid
        if (constants % ui(its,isph).ne.0) then
            do l = 0, params % lmax
                ind = l*l + l + 1
                do m = -l,l
                    f = constants % wgrid(its) * constants % vgrid(ind+m,its) * constants % ui(its,isph)
                    do l0 = 0, constants % lmax0
                        ind0 = l0*l0 + l0 + 1
                        do m0 = -l0, l0
                            f0 = constants % vgrid(ind0+m0,its)
                            pmat(ind+m,ind0+m0) = pmat(ind+m,ind0+m0) + f * f0
                        end do
                    end do
                end do
            end do
        end if
    end do
end subroutine mkpmat

!> Initialize geometry-related constants like list of neighbouring spheres
!!
!! This routine does not check error state of params or constants as it is
!! intended for a low-level use.
!!
!! @param[in] params: Object containing all inputs.
!! @param[inout] constants: Object containing all constants.
subroutine constants_geometry_init(params, constants)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    !! Outputs
    type(ddx_constants_type), intent(inout) :: constants
    !! Local variables
    real(dp) :: swthr, v(3), maxv, ssqv, vv, r, t
    integer :: i, isph, jsph, inear, igrid, iwork, jwork, lwork, &
        & old_lwork, icav, info
    integer, allocatable :: tmp_nl(:), work(:, :), tmp_work(:, :)
    real(dp) :: start_time
    !! The code
    ! Prepare FMM structures if needed
    start_time = omp_get_wtime()
    if (params % fmm .eq. 1) then
        ! Allocate space for a cluster tree
        allocate(constants % order(params % nsph), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `order` " &
                & // "allocation failed"
            return
        end if
        constants % nclusters = 2*params % nsph - 1
        allocate(constants % cluster(2, constants % nclusters), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `cluster` " &
                & // "allocation failed"
            return
        end if
        allocate(constants % children(2, constants % nclusters), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: " // &
                & "`children` allocation failed"
            return
        endif
        allocate(constants % parent(constants % nclusters), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `parent` " &
                & // "allocation failed"
            return
        endif
        allocate(constants % cnode(3, constants % nclusters), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `cnode` " &
                & // "allocation failed"
            return
        endif
        allocate(constants % rnode(constants % nclusters), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `rnode` " &
                & // "allocation failed"
            return
        endif
        allocate(constants % snode(params % nsph), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `snode` " &
                & // "allocation failed"
            return
        endif
        allocate(constants % nfar(constants % nclusters), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `nfar` " &
                & // "allocation failed"
            return
        endif
        allocate(constants % nnear(constants % nclusters), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `nnear` " &
                & // "allocation failed"
            return
        endif
        ! Get the tree
        call tree_rib_build(params % nsph, params % csph, params % rsph, &
            & constants % order, constants % cluster, constants % children, &
            & constants % parent, constants % cnode, constants % rnode, &
            & constants % snode, constants % error_message, &
            & constants % error_flag)
        if (params % error_flag .ne. 0) return
        ! Get number of far and near admissible pairs
        iwork = 0
        jwork = 1
        lwork = 1
        allocate(work(3, lwork), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `work` " &
                & // "allocation failed"
            return
        end if
        do while (iwork .le. jwork)
            allocate(tmp_work(3, lwork), stat=info)
            if (info .ne. 0) then
                constants % error_flag = 1
                constants % error_message = "constants_geometry_init: " // &
                    & "`tmp_work` allocation failed"
                return
            end if
            tmp_work = work
            deallocate(work, stat=info)
            if (info .ne. 0) then
                constants % error_flag = 1
                constants % error_message = "constants_geometry_init: " // &
                    & "`work` deallocation failed"
                return
            end if
            old_lwork = lwork
            lwork = old_lwork + 1000*params % nsph
            allocate(work(3, lwork), stat=info)
            if (info .ne. 0) then
                constants % error_flag = 1
                constants % error_message = "constants_geometry_init: " // &
                    & "`work` allocation failed"
                return
            end if
            work(:, 1:old_lwork) = tmp_work
            deallocate(tmp_work, stat=info)
            if (info .ne. 0) then
                constants % error_flag = 1
                constants % error_message = "constants_geometry_init: " // &
                    & "`tmp_work` deallocation failed"
                return
            end if
            call tree_get_farnear_work(constants % nclusters, &
                & constants % children, constants % cnode, &
                & constants % rnode, lwork, iwork, jwork, work, &
                & constants % nnfar, constants % nfar, constants % nnnear, &
                & constants % nnear)
        end do
        allocate(constants % sfar(constants % nclusters+1), &
            & constants % snear(constants % nclusters+1), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `sfar` " // &
                & "and `snear` allocations failed"
            return
        end if
        allocate(constants % far(constants % nnfar), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `far` " // &
                & "allocation failed"
            return
        end if
        allocate(constants % near(constants % nnnear), stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `near` " // &
                & "allocation failed"
            return
        end if
        call tree_get_farnear(jwork, lwork, work, constants % nclusters, &
            & constants % nnfar, constants % nfar, constants % sfar, &
            & constants % far, constants % nnnear, constants % nnear, &
            & constants % snear, constants % near)
        deallocate(work, stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `work` " // &
                & "deallocation failed"
            return
        end if
    end if
    ! Upper bound of switch region. Defines intersection criterion for spheres
    swthr = one + (params % se+one)*params % eta/two
    ! Assemble neighbor list
    if (params % fmm .eq. 1) then
        call neighbor_list_init_fmm(params, constants)
    else
        call neighbor_list_init(params, constants)
    end if
    ! Allocate space for characteristic functions fi and ui
    allocate(constants % fi(params % ngrid, params % nsph), &
        & constants % ui(params % ngrid, params % nsph), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`fi` and `ui` allocations failed"
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
            return
        endif
        constants % zi = zero
    end if
    ! Build arrays fi, ui, zi
    !$omp parallel do default(none) shared(params,constants,swthr) &
    !$omp private(isph,igrid,jsph,v,maxv,ssqv,vv,t) schedule(dynamic)
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
        return
    endif
    !$omp parallel do default(none) shared(params,constants) &
    !$omp private(isph,i) schedule(dynamic)
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
        return
    endif
    ! Allocate space for characteristic functions ui at cavity points
    allocate(constants % ui_cav(constants % ncav), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`ui_cav` allocations failed"
        return
    end if
    ! Get actual cavity coordinates and indexes in CSR format and fill in
    ! ui_cav aray
    constants % icav_ia(1) = 1
    icav = 1
    do isph = 1, params % nsph
        constants % icav_ia(isph+1) = constants % icav_ia(isph) + &
            & constants % ncav_sph(isph)
        ! Loop over integration points
        do igrid = 1, params % ngrid
            ! Ignore zero contribution
            if (constants % ui(igrid, isph) .eq. zero) cycle
            ! Store coordinates
            constants % ccav(:, icav) = params % csph(:, isph) + &
                & params % rsph(isph)* &
                & constants % cgrid(:, igrid)
            ! Store index
            constants % icav_ja(icav) = igrid
            ! Store characteristic function
            constants % ui_cav(icav) = constants % ui(igrid, isph)
            ! Advance cavity array index
            icav = icav + 1
        end do
    end do
    ! Compute diagonal preconditioner for PCM equations
    if (params % model .eq. 2) then
        allocate(constants % rx_prc( &
            & constants % nbasis, constants % nbasis, params % nsph), &
            & stat=info)
        if (info .ne. 0) then
            constants % error_flag = 1
            constants % error_message = "constants_geometry_init: `rx_prc` " &
                & // "allocation failed"
            return
        endif
        call mkprec(params % lmax, constants % nbasis, params % nsph, &
            & params % ngrid, params % eps, constants % ui, &
            & constants % wgrid, constants % vgrid, constants % vgrid_nbasis, &
            & constants % rx_prc, info, constants % error_message)
        if (info .ne. 0) then
            constants % error_flag = 1
            return
        end if
    end if
end subroutine constants_geometry_init

!> Build the neighbor list using a N^2 code
subroutine neighbor_list_init(params, constants)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    real(dp) :: swthr, v(3), maxv, ssqv, vv, r
    integer :: nngmax, i, lnl, isph, jsph, info
    integer, allocatable :: tmp_nl(:)
    ! Upper bound of switch region. Defines intersection criterion for spheres
    swthr = (params % se+one)*params % eta/two
    ! Build list of neighbours in CSR format
    nngmax = 1
    allocate(constants % inl(params % nsph+1), &
        & constants % nl(params % nsph*nngmax), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`inl` and `nl` allocations failed"
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
                ! switch region `swthr`. Note that we take the largest
                ! switching region as we need a symmetric neighbor list to
                ! use it also for the adjoint products.
                r = params % rsph(isph) + params % rsph(jsph) &
                    & + swthr * max(params % rsph(isph), params % rsph(jsph))
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
                            return
                        end if
                        tmp_nl(1:params % nsph*nngmax) = &
                            & constants % nl(1:params % nsph*nngmax)
                        deallocate(constants % nl, stat=info)
                        if (info .ne. 0) then
                            constants % error_flag = 1
                            constants % error_message = "`nl` " // &
                                & "deallocation failed"
                            return
                        end if
                        nngmax = nngmax + 10
                        allocate(constants % nl(params % nsph*nngmax), &
                            & stat=info)
                        if (info .ne. 0) then
                            constants % error_flag = 1
                            constants % error_message = "`nl` " // &
                                & "allocation failed"
                            return
                        end if
                        constants % nl(1:params % nsph*(nngmax-10)) = &
                            & tmp_nl(1:params % nsph*(nngmax-10))
                        deallocate(tmp_nl, stat=info)
                        if (info .ne. 0) then
                            constants % error_flag = 1
                            constants % error_message = "`tmp_nl` " // &
                                & "deallocation failed"
                            return
                        end if
                    end if
                end if
            end if
        end do
    end do
    constants % inl(params % nsph+1) = lnl+1
    constants % nngmax = nngmax
end subroutine neighbor_list_init

!> Build the neighbor list using a linear scaling code (only if FMMs
!! are enabled)
subroutine neighbor_list_init_fmm(params, constants)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    real(dp) :: swthr, v(3), maxv, ssqv, vv, r
    integer :: nngmax, i, lnl, isph, jsph, inode, jnode, j, k, info
    integer, allocatable :: tmp_nl(:)
    ! Upper bound of switch region. Defines intersection criterion for spheres
    swthr = (params % se+one)*params % eta/two
    ! Build list of neighbours in CSR format
    nngmax = 1
    allocate(constants % inl(params % nsph+1), &
        & constants % nl(params % nsph*nngmax), stat=info)
    if (info .ne. 0) then
        constants % error_flag = 1
        constants % error_message = "`inl` and `nl` allocations failed"
        return
    end if
    i = 1
    lnl = 0
    ! loop over leaf clusters
    do isph = 1, params % nsph
        constants % inl(isph) = lnl + 1
        inode = constants % snode(isph)
        ! loop over near clusters
        do j = constants % snear(inode), constants % snear(inode+1) - 1
            jnode = constants % near(j)
            do k = constants % cluster(1, jnode), constants % cluster(2, jnode)
                jsph = constants % order(k)
                if (isph .ne. jsph) then
                    v = params % csph(:, isph)-params % csph(:, jsph)
                    maxv = max(abs(v(1)), abs(v(2)), abs(v(3)))
                    ssqv = (v(1)/maxv)**2 + (v(2)/maxv)**2 + (v(3)/maxv)**2
                    vv = maxv * sqrt(ssqv)
                    ! Take regularization parameter into account with respect to
                    ! shift se. It is described properly by the upper bound of a
                    ! switch region `swthr`. Note that we take the largest
                    ! switching region as we need a symmetric neighbor list to
                    ! use it also for the adjoint products.
                    r = params % rsph(isph) + params % rsph(jsph) &
                        & + swthr * max(params % rsph(isph), params % rsph(jsph))
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
                                return
                            end if
                            tmp_nl(1:params % nsph*nngmax) = &
                                & constants % nl(1:params % nsph*nngmax)
                            deallocate(constants % nl, stat=info)
                            if (info .ne. 0) then
                                constants % error_flag = 1
                                constants % error_message = "`nl` " // &
                                    & "deallocation failed"
                                return
                            end if
                            nngmax = nngmax + 10
                            allocate(constants % nl(params % nsph*nngmax), &
                                & stat=info)
                            if (info .ne. 0) then
                                constants % error_flag = 1
                                constants % error_message = "`nl` " // &
                                    & "allocation failed"
                                return
                            end if
                            constants % nl(1:params % nsph*(nngmax-10)) = &
                                & tmp_nl(1:params % nsph*(nngmax-10))
                            deallocate(tmp_nl, stat=info)
                            if (info .ne. 0) then
                                constants % error_flag = 1
                                constants % error_message = "`tmp_nl` " // &
                                    & "deallocation failed"
                                return
                            end if
                        end if
                    end if
                end if
            end do
        end do
    end do
    constants % inl(params % nsph+1) = lnl+1
    constants % nngmax = nngmax
end subroutine neighbor_list_init_fmm

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
subroutine constants_geometry_update(params, constants)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    !! Outputs
    type(ddx_constants_type), intent(out) :: constants
    !! Local variables
    !! The code
    ! Enforce error for now (not yet implemented)
    constants % error_message = 'constants geometry update is not yet implemented'
    constants % error_flag = 1
end subroutine constants_geometry_update

!> Switching function
!!
!! This is an implementation of \f$ \chi(t) \f$ with a shift \f$ se \f$:
!! \f[
!!      \chi_\eta(t) = \left\{ \begin{array}{ll} 0, & \text{if} \quad x \geq
!!      1 \\ p_\eta(x), & \text{if} \quad 1-\eta < x < 1 \\ 1, & \text{if}
!!      \quad x \leq 1-\eta \end{array} \right.
!! \f]
!! where \f$ x = t - \frac{1+se}{2} \eta \f$ is a shifted coordinate and
!! \f[
!!      p_\eta(x) = \frac{1}{\eta^5} (1-t)^3 (6t^2 + (15\eta-12)t + (10\eta^2
!!      -15\eta+6))
!! \f]
!! is a smoothing polynomial of the 5th degree.
!! In the case shift \f$ se=1 \f$ the switch function \f$ \chi_\eta(t) \f$ is 1 for
!! \f$ t \in [0,1] \f$, is 0 for \f$ t \in [1+\eta, \infty) \f$ and varies in
!! \f$ [1, 1+\eta] \f$ which is an external shift.
!! In the case shift \f$ se=-1 \f$ the switch function \f$ \chi_\eta(t) \f$ is 1 for
!! \f$ t \in [0,1-\eta] \f$, is 0 for \f$ t \in [1, \infty) \f$ and varies in
!! \f$ [1-\eta, 1] \f$ which is an internal shift.
!! In the case shift \f$ se=0 \f$ the switch function \f$ \chi_\eta(t) \f$ is 1 for
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
!!      \chi_\eta'(t) = \left\{ \begin{array}{ll} 0, & \text{if} \quad x \geq
!!      1 \\ p'_\eta(x), & \text{if} \quad 1-\eta < x < 1 \\ 0, & \text{if}
!!      \quad x \leq 1-\eta \end{array} \right.
!! \f]
!! where \f$ x = t - \frac{1+se}{2} \eta \f$ is a shifted coordinate and
!! \f[
!!      p'_\eta(x) = -\frac{30}{\eta^5} (1-t)^2 (t-1+\eta)^2
!! \f]
!! is a derivative of the smoothing polynomial.
!! In the case shift \f$ se=1 \f$ the derivative \f$ \chi_\eta'(t) \f$ is 0 for
!! \f$ t \in [0,1] \cup [1+\eta, \infty) \f$ and varies in
!! \f$ [1, 1+\eta] \f$ which is an external shift.
!! In the case shift \f$ se=-1 \f$ the derivative \f$ \chi_\eta'(t) \f$ is 0 for
!! \f$ t \in [0,1-\eta] \cup [1, \infty) \f$ and varies in
!! \f$ [1-\eta, 1] \f$ which is an internal shift.
!! In the case shift \f$ se=0 \f$ the derivative \f$ \chi_\eta'(t) \f$ is 0 for
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

!> Compute preconditioner
!!
!! assemble the diagonal blocks of the reps matrix
!! then invert them to build the preconditioner
subroutine mkprec(lmax, nbasis, nsph, ngrid, eps, ui, wgrid, vgrid, &
        & vgrid_nbasis, rx_prc, info, error_message)
    !! Inputs
    integer, intent(in) :: lmax, nbasis, nsph, ngrid, vgrid_nbasis
    real(dp), intent(in) :: eps, ui(ngrid, nsph), wgrid(ngrid), &
        & vgrid(vgrid_nbasis, ngrid)
    !! Output
    real(dp), intent(out) :: rx_prc(nbasis, nbasis, nsph)
    integer, intent(out) :: info
    character(len=255), intent(out) :: error_message
    !! Local variables
    integer :: isph, lm, ind, l1, m1, ind1, igrid
    real(dp)  :: f, f1
    integer, allocatable :: ipiv(:)
    real(dp),  allocatable :: work(:)
    external :: dgetrf, dgetri
    ! Allocation of temporaries
    allocate(ipiv(nbasis), work(nbasis**2), stat=info)
    if (info .ne. 0) then
        error_message = "mkprec: `ipiv` and `work` allocation failed"
        info = 1
        return
    endif
    ! Init
    rx_prc = zero
    ! Off-diagonal part
    do isph = 1, nsph
        do igrid = 1, ngrid
            f = twopi * ui(igrid, isph) * wgrid(igrid)
            do l1 = 0, lmax
                ind1 = l1*l1 + l1 + 1
                do m1 = -l1, l1
                    f1 = f * vgrid(ind1+m1, igrid) / (two*dble(l1) + one)
                    do lm = 1, nbasis
                        rx_prc(lm, ind1+m1, isph) = &
                            & rx_prc(lm, ind1+m1, isph) + f1*vgrid(lm, igrid)
                    end do
                end do
            end do
        end do
    end do
    ! add diagonal
    f = twopi * (eps+one) / (eps-one)
    do isph = 1, nsph
        do lm = 1, nbasis
            rx_prc(lm, lm, isph) = rx_prc(lm, lm, isph) + f
        end do
    end do
    ! invert the blocks
    do isph = 1, nsph
        call dgetrf(nbasis, nbasis, rx_prc(:, :, isph), nbasis, ipiv, info)
        if (info .ne. 0) then 
            error_message = "mkprec: dgetrf failed"
            info = 1
            return
        end if
        call dgetri(nbasis, rx_prc(:, :, isph), nbasis, ipiv, work, &
            & nbasis**2, info)
        if (info .ne. 0) then 
            error_message = "mkprec: dgetri failed"
            info = 1
            return
        end if
    end do
    ! Cleanup temporaries
    deallocate(work, ipiv, stat=info)
    if (info .ne. 0) then
        error_message = "mkprec: `ipiv` and `work` deallocation failed"
        info = 1
        return
    end if
  ! Cleanup error message if there were no error
    error_message = ""
end subroutine mkprec

!> Build a recursive inertial binary tree
!!
!! Uses inertial bisection in a recursive manner until each leaf node has only
!! one input sphere inside. Number of tree nodes is always 2*nsph-1.
!!
!!
!! @param[in] nsph: Number of input spheres
!! @param[in] csph: Centers of input spheres
!! @param[in] rsph: Radii of input spheres
!! @param[out] order: Ordering of input spheres to make spheres of one cluster
!!      contiguous.
!! @param[out] cluster: Indexes in `order` array of the first and the last
!!      spheres of each cluster
!! @param[out] children: Indexes of the first and the last children nodes. If
!!      both indexes are equal this means there is only 1 child node and if
!!      value of the first child node is 0, there are no children nodes.
!! @param[out] parent: Parent of each cluster. Value 0 means there is no parent
!!      node which corresponds to the root node (with index 1).
!! @param[out] cnode: Center of a bounding sphere of each node
!! @param[out] rnode: Radius of a bounding sphere of each node
!! @param[out] snode: Array of leaf nodes containing input spheres
subroutine tree_rib_build(nsph, csph, rsph, order, cluster, children, parent, &
        & cnode, rnode, snode, error_message, error_flag)
    ! Inputs
    integer, intent(in) :: nsph
    real(dp), intent(in) :: csph(3, nsph), rsph(nsph)
    ! Outputs
    integer, intent(out) :: order(nsph), cluster(2, 2*nsph-1), &
        & children(2, 2*nsph-1), parent(2*nsph-1), snode(nsph)
    real(dp), intent(out) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    integer, intent(out) :: error_flag
    character(len=255), intent(out) :: error_message
    ! Local variables
    integer :: nclusters, i, j, n, s, e, div
    real(dp) :: r, r1, r2, c(3), c1(3), c2(3), d, maxc, ssqc
    !! At first construct the tree
    nclusters = 2*nsph - 1
    ! Init the root node
    cluster(1, 1) = 1
    cluster(2, 1) = nsph
    parent(1) = 0
    ! Index of the first unassigned node
    j = 2
    ! Init spheres ordering
    do i = 1, nsph
        order(i) = i
    end do
    ! Divide nodes until a single spheres
    do i = 1, nclusters
        ! The first sphere in i-th node
        s = cluster(1, i)
        ! The last sphere in i-th node
        e = cluster(2, i)
        ! Number of spheres in i-th node
        n = e - s + 1
        ! Divide only if there are 2 or more spheres
        if (n .gt. 1) then
            ! Use inertial bisection to reorder spheres and cut into 2 halfs
            call tree_rib_node_bisect(nsph, csph, n, order(s:e), div, &
                & error_message, error_flag)
            if (error_flag .ne. 0) return
            ! Assign the first half to the j-th node
            cluster(1, j) = s
            cluster(2, j) = s + div - 1
            ! Assign the second half to the (j+1)-th node
            cluster(1, j+1) = s + div
            cluster(2, j+1) = e
            ! Update list of children of i-th node
            children(1, i) = j
            children(2, i) = j + 1
            ! Set parents of new j-th and (j+1)-th node
            parent(j) = i
            parent(j+1) = i
            ! Shift index of the first unassigned node
            j = j + 2
        ! Set information for a leaf node
        else
            ! No children nodes
            children(:, i) = 0
            ! i-th node contains sphere ind(s)
            snode(order(s)) = i
        end if
    end do
    !! Compute bounding spheres of each node of the tree
    ! Bottom-to-top pass over all nodes of the tree
    do i = nclusters, 1, -1
        ! In case of a leaf node just use corresponding input sphere as a
        ! bounding sphere of the node
        if (children(1, i) .eq. 0) then
            ! Get correct index of the corresponding input sphere
            j = order(cluster(1, i))
            ! Get corresponding center and radius
            cnode(:, i) = csph(:, j)
            rnode(i) = rsph(j)
        ! In case of a non-leaf node get minimal sphere that contains bounding
        ! spheres of children nodes
        else
            ! The first child
            j = children(1, i)
            c1 = cnode(:, j)
            r1 = rnode(j)
            ! The second child
            j = children(2, i)
            c2 = cnode(:, j)
            r2 = rnode(j)
            ! Distance between centers of bounding spheres of children nodes
            c = c1 - c2
            ! Compute distance by scale and sum of scaled squares
            maxc = max(abs(c(1)), abs(c(2)), abs(c(3)))
            if (maxc .eq. zero) then
                d = zero
            else
                d = (c(1)/maxc)**2 + (c(2)/maxc)**2 + (c(3)/maxc)**2
                d = maxc * sqrt(d)
            end if
            ! If sphere #2 is completely inside sphere #1 use the first sphere
            if (r1 .ge. (r2+d)) then
                c = c1
                r = r1
            ! If sphere #1 is completely inside sphere #2 use the second sphere
            else if(r2 .ge. (r1+d)) then
                c = c2
                r = r2
            ! Otherwise use special formula to find a minimal sphere
            else
                r = (r1+r2+d) / 2
                c = c2 + c/d*(r-r2)
            end if
            ! Assign computed bounding sphere
            cnode(:, i) = c
            rnode(i) = r
        end if
    end do
end subroutine tree_rib_build

!> Divide given cluster of spheres into two subclusters by inertial bisection
!!
!! @param[in] nsph: Number of all input spheres
!! @param[in] csph: Centers of all input spheres
!! @param[in] n: Number of spheres in a given cluster
!! @param[inout] order: Indexes of spheres in a given cluster. On exit, indexes
!!      `order(1:div)` correspond to the first subcluster and indexes
!!      `order(div+1:n)` correspond to the second subcluster.
!! @param[out] div: Break point of `order` array between two clusters.
subroutine tree_rib_node_bisect(nsph, csph, n, order, div, error_message, &
        & error_flag)
    ! Inputs
    integer, intent(in) :: nsph, n
    real(dp), intent(in) :: csph(3, nsph)
    ! Outputs
    integer, intent(inout) :: order(n)
    integer, intent(out) :: div
    integer, intent(out) :: error_flag
    character(len=255), intent(out) :: error_message
    ! Local variables
    real(dp) :: c(3),  s(3)
    real(dp), allocatable :: tmp_csph(:, :), work(:)
    external :: dgesvd
    integer :: i, l, r, lwork, info, istat
    integer, allocatable :: tmp_order(:)

    allocate(tmp_csph(3, n), tmp_order(n), stat=istat)
    if (istat.ne.0) then
        error_message = 'Allocation failed in tree_node_bisect'
        error_flag = 1
        return
    end if

    ! Get average coordinate
    c = zero
    do i = 1, n
        c = c + csph(:, order(i))
    end do
    c = c / n
    ! Get coordinates minus average in a contiguous array
    do i = 1, n
        tmp_csph(:, i) = csph(:, order(i)) - c
    end do
    !! Find right singular vectors
    ! Get proper size of temporary workspace
    lwork = -1
    call dgesvd('N', 'O', 3, n, tmp_csph, 3, s, tmp_csph, 3, tmp_csph, 3, &
        & s, lwork, info)
    lwork = int(s(1))
    allocate(work(lwork), stat=istat)
    if (istat.ne.0) then
        error_message = 'Allocation failed in tree_node_bisect'
        error_flag = 1
        return
    end if
    ! Get right singular vectors
    call dgesvd('N', 'O', 3, n, tmp_csph, 3, s, tmp_csph, 3, tmp_csph, 3, &
        & work, lwork, info)
    if (info.ne.0) then
        error_message = 'DGESVD failed in tree_node_bisect'
        error_flag = 1
        return
    end if
    deallocate(work, stat=istat)
    if (istat.ne.0) then
        error_message = 'Deallocation failed in tree_node_bisect'
        error_flag = 1
        return
    end if
    !! Sort spheres by sign of the above scalar product, which is equal to
    !! the leading right singular vector scaled by the leading singular value.
    !! However, we only care about sign, so we take into account only the
    !! leading right singular vector.
    ! First empty index from the left
    l = 1
    ! First empty index from the right
    r = n
    ! Cycle over values of the singular vector
    do i = 1, n
        ! Positive scalar products are moved to the beginning of `order`
        if (tmp_csph(1, i) .ge. zero) then
            tmp_order(l) = order(i)
            l = l + 1
        ! Negative scalar products are moved to the end of `order`
        else
            tmp_order(r) = order(i)
            r = r - 1
        end if
    end do
    ! Set divider and update order
    div = r
    order = tmp_order
    deallocate(tmp_csph, tmp_order, stat=istat)
    if (istat.ne.0) then
        error_message = 'Deallocation failed in tree_node_bisect'
        error_flag = 1
        return
    end if
end subroutine tree_rib_node_bisect

!> Find near and far admissible pairs of tree nodes and store it in work array
!!
!! @param[in] n: number of nodes
!! @param[in] children: first and last children of each cluster. Value 0 means
!!      no children (leaf node).
!! @param[in] cnode: center of bounding sphere of each cluster (node) of tree
!! @param[in] rnode: radius of bounding sphere of each cluster (node) of tree
!! @param[in] lwork: size of work array in dimension 2
!! @param[inout] iwork: index of current pair of nodes that needs to be checked
!!      for admissibility. Must be 0 for the first call of this subroutine. If
!!      on exit iwork is less or equal to jwork, that means lwork was too
!!      small, please reallocate work array and copy all the values into new
!!      array and then run procedure again.
!! @param[inout] jwork: amount of stored possible admissible pairs of nodes.
!!      Please read iwork comments.
!! @param[inout] work: all the far and near pairs will be stored here
!! @param[out] nnfar: total amount of far admissible pairs. valid only if iwork
!!      is greater than jwork on exit.
!! @param[out] nfar: amount of far admissible pairs for each node. valid only
!!      if iwork is greater than jwork on exit.
!! @param[out] nnnear: total amount of near admissible pairs. valid only if
!!      iwork is greater than jwork on exit
!! @param[out] nnear: amount of near admissible pairs for each node. valid only
!!      if iwork is greater than jwork on exit
subroutine tree_get_farnear_work(n, children, cnode, rnode, lwork, iwork, &
        & jwork, work, nnfar, nfar, nnnear, nnear)
    ! Inputs
    integer, intent(in) :: n, children(2, n), lwork
    real(dp), intent(in) :: cnode(3, n), rnode(n)
    ! Outputs
    integer, intent(inout) :: iwork, jwork, work(3, lwork)
    integer, intent(out) :: nnfar, nfar(n), nnnear, nnear(n)
    ! Local variables
    integer :: j(2), npairs, k1, k2
    real(dp) :: c(3), r, d, dmax, dssq
    ! iwork is current temporary item in work array to process
    if (iwork .eq. 0) then
        work(1, 1) = 1
        work(2, 1) = 1
        iwork = 1
        jwork = 1
    end if
    ! jwork is total amount of temporary items in work array
    do while (iwork .le. jwork)
        j = work(1:2, iwork)
        c = cnode(:, j(1)) - cnode(:, j(2))
        r = rnode(j(1)) + rnode(j(2)) + max(rnode(j(1)), rnode(j(2)))
        !r = rnode(j(1)) + rnode(j(2))
        dmax = max(abs(c(1)), abs(c(2)), abs(c(3)))
        dssq = (c(1)/dmax)**2 + (c(2)/dmax)**2 + (c(3)/dmax)**2
        d = dmax * sqrt(dssq)
        !d = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
        ! If node has no children then assume itself for purpose of finding
        ! far-field and near-filed interactions with children nodes of another
        ! node
        npairs = max(1, children(2, j(1))-children(1, j(1))+1) * &
            & max(1, children(2, j(2))-children(1, j(2))+1)
        if (d .ge. r) then
            ! Mark as far admissible pair
            work(3, iwork) = 1
        else if (npairs .eq. 1) then
            ! Mark as near admissible pair if both nodes are leaves
            work(3, iwork) = 2
        else if (jwork+npairs .gt. lwork) then
            ! Exit procedure, since work array was too small
            return
        else
            ! Mark as non-admissible pair and check all pairs of children nodes
            ! or pairs of one node (if it is a leaf node) with children of
            ! another node
            work(3, iwork) = 0
            if (children(1, j(1)) .eq. 0) then
                k1 = j(1)
                do k2 = children(1, j(2)), children(2, j(2))
                    jwork = jwork + 1
                    work(1, jwork) = k1
                    work(2, jwork) = k2
                end do
            else if(children(1, j(2)) .eq. 0) then
                k2 = j(2)
                do k1 = children(1, j(1)), children(2, j(1))
                    jwork = jwork + 1
                    work(1, jwork) = k1
                    work(2, jwork) = k2
                end do
            else
                do k1 = children(1, j(1)), children(2, j(1))
                    do k2 = children(1, j(2)), children(2, j(2))
                        jwork = jwork + 1
                        work(1, jwork) = k1
                        work(2, jwork) = k2
                    end do
                end do
            end if
        end if
        iwork = iwork + 1
    end do
    nfar = 0
    nnear = 0
    do iwork = 1, jwork
        if (work(3, iwork) .eq. 1) then
            nfar(work(1, iwork)) = nfar(work(1, iwork)) + 1
        else if (work(3, iwork) .eq. 2) then
            nnear(work(1, iwork)) = nnear(work(1, iwork)) + 1
        end if
    end do
    iwork = jwork + 1
    nnfar = sum(nfar)
    nnnear = sum(nnear)
end subroutine tree_get_farnear_work

!> Get near and far admissible pairs from work array of tree_get_farnear_work
!! Works only for binary tree
subroutine tree_get_farnear(jwork, lwork, work, n, nnfar, nfar, sfar, far, &
        & nnnear, nnear, snear, near)
!   Parameters:
!   jwork: Total number of checked pairs in work array
!   lwork: Total length of work array
!   work: Work array itself
!   n: Number of nodes
!   nnfar: Total number of all far-field interactions
!   nfar: Number of far-field interactions of each node
!   sfar: Index in far array of first far-field node for each node
!   far: Indexes of far-field nodes
!   nnnear: Total number of all near-field interactions
!   nnear: Number of near-field interactions of each node
!   snear: Index in near array of first near-field node for each node
!   near: Indexes of near-field nodes
    integer, intent(in) :: jwork, lwork, work(3, lwork), n, nnfar, nnnear
    integer, intent(in) :: nfar(n), nnear(n)
    integer, intent(out) :: sfar(n+1), far(nnfar), snear(n+1), near(nnnear)
    integer :: i, j
    integer, allocatable :: cfar(:), cnear(:)

    allocate(cfar(n+1), cnear(n+1))
    sfar(1) = 1
    snear(1) = 1
    do i = 2, n+1
        sfar(i) = sfar(i-1) + nfar(i-1)
        snear(i) = snear(i-1) + nnear(i-1)
    end do
    cfar = sfar
    cnear = snear
    do i = 1, jwork
        if (work(3, i) .eq. 1) then
            ! Far
            j = work(1, i)
            if ((j .gt. n) .or. (j .le. 0)) then
                write(*,*) "ALARM", j
            end if
            far(cfar(j)) = work(2, i)
            cfar(j) = cfar(j) + 1
        else if (work(3, i) .eq. 2) then
            ! Near
            j = work(1, i)
            if ((j .gt. n) .or. (j .le. 0)) then
                write(*,*) "ALARM", j
            end if
            near(cnear(j)) = work(2, i)
            cnear(j) = cnear(j) + 1
        end if
    end do

    deallocate(cfar, cnear)

end subroutine tree_get_farnear

!> Deallocate the constants
!> @ingroup Fortran_interface_core
!!
!! @param[out] constants: Precomputed constants
!!
subroutine constants_free(constants)
    implicit none
    type(ddx_constants_type), intent(out) :: constants
    integer :: istat

    istat = 0

    if (allocated(constants % vscales)) then
        deallocate(constants % vscales, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`vscales` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % v4pi2lp1)) then
        deallocate(constants % v4pi2lp1, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`v4pi2lp1` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % vscales_rel)) then
        deallocate(constants % vscales_rel, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`vscales_rel` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % vfact)) then
        deallocate(constants % vfact, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`vfact` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % vcnk)) then
        deallocate(constants % vcnk, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`vcnk` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % m2l_ztranslate_coef)) then
        deallocate(constants % m2l_ztranslate_coef, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`m2l_ztranslate_coef` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % m2l_ztranslate_adj_coef)) then
        deallocate(constants % m2l_ztranslate_adj_coef, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`m2l_ztranslate_adj_coef` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % cgrid)) then
        deallocate(constants % cgrid, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`cgrid` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % wgrid)) then
        deallocate(constants % wgrid, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`wgrid` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % vgrid)) then
        deallocate(constants % vgrid, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`vgrid` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % vwgrid)) then
        deallocate(constants % vwgrid, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`vwgrid` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % vgrid2)) then
        deallocate(constants % vgrid2, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`vgrid2` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % pchi)) then
        deallocate(constants % pchi, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`pchi` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % c_ik)) then
        deallocate(constants % c_ik, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`c_ik` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % si_ri)) then
        deallocate(constants % si_ri, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`si_ri` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % di_ri)) then
        deallocate(constants % di_ri, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`di_ri` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % sk_ri)) then
        deallocate(constants % sk_ri, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`sk_ri` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % dk_ri)) then
        deallocate(constants % dk_ri, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`dk_ri` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % termimat)) then
        deallocate(constants % termimat, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`termimat` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % b)) then
        deallocate(constants % b, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`b` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % l)) then
        deallocate(constants % l, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`l` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % inl)) then
        deallocate(constants % inl, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`inl` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % nl)) then
        deallocate(constants % nl, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`nl` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % itrnl)) then
        deallocate(constants % itrnl, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`itrnl` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % fi)) then
        deallocate(constants % fi, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`fi` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % ui)) then
        deallocate(constants % ui, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`ui` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % ui_cav)) then
        deallocate(constants % ui_cav, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`ui_cav` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % zi)) then
        deallocate(constants % zi, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`zi` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % ncav_sph)) then
        deallocate(constants % ncav_sph, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`ncav_sph` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % ccav)) then
        deallocate(constants % ccav, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`ccav` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % icav_ia)) then
        deallocate(constants % icav_ia, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`icav_ia` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % icav_ja)) then
        deallocate(constants % icav_ja, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`icav_ja` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % rx_prc)) then
        deallocate(constants % rx_prc, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`rx_prc` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % order)) then
        deallocate(constants % order, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`order` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % cluster)) then
        deallocate(constants % cluster, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`cluster` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % children)) then
        deallocate(constants % children, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`children` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % parent)) then
        deallocate(constants % parent, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`parent` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % cnode)) then
        deallocate(constants % cnode, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`cnode` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % rnode)) then
        deallocate(constants % rnode, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`rnode` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % snode)) then
        deallocate(constants % snode, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`snode` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % sk_rnode)) then
        deallocate(constants % sk_rnode, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`sk_rnode` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % si_rnode)) then
        deallocate(constants % si_rnode, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`si_rnode` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % nfar)) then
        deallocate(constants % nfar, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`nfar` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % nnear)) then
        deallocate(constants % nnear, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`nnear` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % far)) then
        deallocate(constants % far, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`far` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % near)) then
        deallocate(constants % near, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`near` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % sfar)) then
        deallocate(constants % sfar, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`sfar` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
    if (allocated(constants % snear)) then
        deallocate(constants % snear, stat=istat)
        if (istat .ne. 0) then
            constants % error_message = "`snear` deallocation failed!"
            constants % error_flag = 1
            return
        end if
    end if
end subroutine constants_free

end module ddx_constants

