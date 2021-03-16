!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_core.f90
!! Core routines and parameters of entire ddX software
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

!> Core routines and parameters of ddX software
module ddx_core
! Enable OpenMP
use omp_lib
implicit none

!> Kind of double precision
integer, parameter :: dp = kind(1d0)

!! Compile-time constants
real(dp), parameter :: zero = 0d0, one = 1d0, two = 2d0, three = 3d0
real(dp), parameter :: four = 4d0, pt5 = 5d-1
real(dp), parameter :: sqrt2 = sqrt(two)
real(dp), parameter :: sqrt3 = sqrt(three)
real(dp), parameter :: pi4 = atan(one)
real(dp), parameter :: pi = four * pi4
real(dp), parameter :: fourpi = four * pi
real(dp), parameter :: twopi = two * pi
real(dp), parameter :: sqrt4pi = four * sqrt(pi4)
real(dp), parameter :: machine_eps = epsilon(zero)
real(dp), parameter :: toang = 0.52917721092d0
real(dp), parameter :: tokcal = 627.509469d0
real(dp), parameter :: tobohr = one / toang
!> Number of supported Lebedev grids
integer, parameter :: nllg = 32
!> Number of grid points of each Lebedev grid
integer, parameter :: ng0(nllg) = (/ 6, 14, 26, 38, 50, 74, 86, 110, 146, &
    & 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, &
    & 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 /)

!> Main ddX type that stores all required information
type ddx_type
    !!!!!!!!!!!! Parameters
    !> Model to use 1 for cosmo, 2 for pcm, 3 for lpb.
    integer :: model
    !> Verbosity level. 0 for no output, 1 for ... Larger value means more
    !!      output.
    integer :: iprint
    !> Number of OpenMP threads to be used. Currently, only nproc=1 is
    !!      supported as ddX is sequential right now.
    integer :: nproc
    !> Number of atoms in the molecule.
    integer :: nsph
    !> Charges of atoms of a dimension (nsph).
    real(dp), allocatable :: charge(:)
    !> Centers of atoms of a dimension (3, nsph).
    real(dp), allocatable :: csph(:, :)
    !> Array of radii of atoms of a dimension (nsph).
    real(dp), allocatable :: rsph(:)
    !> Dielectric permittivity inside cavity.
    real(dp) :: epsin
    !> Dielectric permittivity outside cavity.
    real(dp) :: epsout
    !> Relative dielectric permittivity.
    real(dp) :: eps
    !> Debye H\"{u}ckel parameter. Referenced only in LPB model (model=3)
    real(dp) :: kappa
    !> Shift of a characteristic function. -1 for interior, 0 for centered and
    !!      1 for outer regularization.
    real(dp) :: se
    !> Regularization parameter.
    real(dp) :: eta
    !> Iterative solver to be used. 1 for Jacobi/DIIS.
    !!
    !! Other solvers might be added later.
    integer :: itersolver
    !> Relative threshold for an iterative solvers.
    real(dp) :: tol
    !> Maximum number of iterations for an iterative solver.
    integer :: maxiter
    !> Number of extrapolation points for Jacobi/DIIS solver. Referenced only
    !!      if Jacobi solver is used.
    integer :: ndiis
    !> Maximal degree of modeling spherical harmonics.
    integer :: lmax
    !> Number modeling spherical harmonics per sphere.
    integer :: nbasis
    !> Total number of modeling degrees of freedom.
    !!
    !! This is equal to `nsph*nbasis`.
    integer :: n
    !> Number of Lebedev grid points on each sphere.
    integer :: ngrid
    !> Whether to compute (1) or not (0) analytical forces.
    integer :: force
    !> Enable (1) or disable (0) use of FMM techniques.
    integer :: fmm
    !> Maximum degree of spherical harmonics for M (multipole) expansion.
    !!
    !! If this value is -1 then no far-field FMM interactions are performed.
    integer :: pm
    !> Maximum degree of spherical harmonics for L (local) expansion.
    !!
    !! If this value is -1 then no far-field FMM interactions are performed.
    integer :: pl
    !> Flag if FMM transfer matrices have to be precomputed
    integer :: fmm_precompute
    !!!!!!!!!!!!!! Constants, initialized by ddinit
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
    !! Reflection matrices for M2M, M2L and L2L operations. Allocated and used
    !! only if fmm=1 and fmm_precompute=1
    !> Size of each M2M reflection matrix. Defined only if fmm=1 and
    !!      fmm_precompute=1.
    integer :: m2m_reflect_mat_size
    !> Array of M2M reflection matrices. Dimension is (m2m_reflect_mat_size,
    !!      nclusters). Allocated and used only if fmm=1 and fmm_precompue=1.
    real(dp), allocatable :: m2m_reflect_mat(:, :)
    !> Size of each M2M OZ translation matrix. Defined only if fmm=1 and
    !!      fmm_precompute=1.
    integer :: m2m_ztranslate_mat_size
    !> Array of M2M OZ translation matrices. Dimension is
    !!      (m2m_ztranslate_mat_size, nclusters). Allocated and used if fmm=1
    !!      and fmm_precompute=1.
    real(dp), allocatable :: m2m_ztranslate_mat(:, :)
    !> Size of each L2L reflection matrix. Defined only if fmm=1 and
    !!      fmm_precompute=1.
    integer :: l2l_reflect_mat_size
    !> Array of L2L reflection matrices. Dimension is (l2l_reflect_mat_size,
    !!      nclusters). Allocated and used if fmm=1 and fmm_precompute=1.
    real(dp), allocatable :: l2l_reflect_mat(:, :)
    !> Size of each L2L OZ translation matrix. Defined only if fmm=1 and
    !!      fmm_precompute=1.
    integer :: l2l_ztranslate_mat_size
    !> Array of L2L OZ translation matrices. Dimension is
    !!      (l2l_ztranslate_mat_size, nclusters). Allocated and used only if
    !!      fmm=1 and fmm_precompute=1.
    real(dp), allocatable :: l2l_ztranslate_mat(:, :)
    !> Size of each M2L reflection matrix. Defined only if fmm=1 and
    !!      fmm_precompute=1.
    integer :: m2l_reflect_mat_size
    !> Array of M2L reflection matrices. Dimension is (m2l_reflect_mat_size,
    !!      nclusters). Allocated and used if fmm=1 and fmm_precompute=1.
    real(dp), allocatable :: m2l_reflect_mat(:, :)
    !> Size of each M2L OZ translation matrix. Defined only if fmm=1 and
    !!      fmm_precompute=1.
    integer :: m2l_ztranslate_mat_size
    !> Array of M2L OZ translation matrices. Dimension is
    !!      (m2l_ztranslate_mat_size, nclusters). Allocated and used only if
    !!      fmm=1 and fmm_precompute=1.
    real(dp), allocatable :: m2l_ztranslate_mat(:, :)
    !> Number of near-field M2P interactions with cavity points. Defined only
    !!      if fmm=1.
    integer :: nnear_m2p
    !> Maximal degree of near-field M2P spherical harmonics. Defined only if
    !!      fmm=1.
    integer :: m2p_lmax
    !> Number of spherical harmonics used for near-field M2P. Defined only if
    !!      fmm=1.
    integer :: m2p_nbasis
    !> Near-field M2P matrices of a dimension (m2p_nbasis, nnear_m2p).
    !!      Allocated and used only if fmm=1 and fmm_precompute=1.
    real(dp), allocatable :: m2p_mat(:, :)
    !> Potential at all grid points. Dimension is (ngrid, nsph). Allocated
    !!      and used by COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: phi_grid(:, :)
    !> Variable \f$ \Phi \f$ of a dimension (nbasis, nsph). Allocated and used
    !!      by COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: phi(:, :)
    !> Variable \f$ \Phi_\infty \f$ of a dimension (nbasis, nsph). Allocated
    !!      and used only by PCM (model=2) model.
    real(dp), allocatable :: phiinf(:, :)
    !> Variable \f$ \Phi_\varepsilon \f$ of a dimension (nbasis, nsph).
    !!      Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: phieps(:, :)
    !> Solution of the ddCOSMO system of a dimension (nbasis, nsph). Allocated
    !!      and used by COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: xs(:, :)
    !> Solution of the adjoint ddCOSMO system of a dimension (nbasis, nsph).
    !!      Allocated and used by COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: s(:, :)
    !> Values of s at grid points. Dimension is (ngrid, nsph). Allocated and
    !!      used by COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: sgrid(:, :)
    !> Solution of the adjoint ddPCM system of a dimension (nbasis, nsph).
    !!      Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: y(:, :)
    !> Values of y at grid points. Dimension is (ngrid, nsph). Allocated and
    !!      used only by PCM (model=2) model.
    real(dp), allocatable :: ygrid(:, :)
    !> Shortcut of \f$ \Phi_\varepsilon - \Phi \f$
    real(dp), allocatable :: g(:, :)
    !> s-4pi/(eps-1)y of dimension (nbasis, nsph)
    real(dp), allocatable :: q(:, :)
    !> Values of q at grid points. Dimension is (ngrid, nsph)
    real(dp), allocatable :: qgrid(:, :)
    !> Zeta intermediate for forces. Dimension is (ncav)
    real(dp), allocatable :: zeta(:)
    !> Temporary workspace for associated legendre polynomials. Dimension is
    !!      (vgrid_nbasis, nproc).
    real(dp), allocatable :: tmp_vplm(:, :)
    !> Temporary workspace for an array of cosinuses of a dimension
    !!      (vgrid_dmax+1, nproc).
    real(dp), allocatable :: tmp_vcos(:, :)
    !> Temporary workspace for an array of sinuses of a dimension
    !!      (vgrid_dmax+1, nproc).
    real(dp), allocatable :: tmp_vsin(:, :)
    !> Temporary workspace for multipole coefficients of a degree up to lmax
    !!      of each sphere. Dimension is (nbasis, nsph).
    real(dp), allocatable :: tmp_sph(:, :)
    !> Number of spherical harmonics of degree up to lmax+1 used for
    !!      computation of forces (gradients). Allocated and used only if
    !!      fmm=1.
    integer :: grad_nbasis
    !> Temporary workspace for multipole coefficients of a degree up to lmax+1
    !!      of each sphere. Dimension is (grad_nbasis, nsph). Allocated and
    !!      used only if fmm=1.
    real(dp), allocatable :: tmp_sph2(:, :)
    !> Temporary workspace for a gradient of M2M of harmonics of a degree up to
    !!      lmax+1 of each sphere. Dimension is ((grad_nbasis, 3, nsph).
    real(dp), allocatable :: tmp_sph_grad(:, :, :)
    !> Temporary workspace for local coefficients of a degree up to pl
    !!      of each sphere. Dimension is ((pl+1)**2, nsph).
    real(dp), allocatable :: tmp_sph_l(:, :)
    !> Temporary workspace for a gradient of L2L of harmonics of a degree up to
    !!      pl of each sphere. Dimension is ((pl+1)**2, 3, nsph).
    real(dp), allocatable :: tmp_sph_l_grad(:, :, :)
    !> Temporary workspace for multipole coefficients of each node. Dimension
    !!      is ((pm+1)**2, nsph)
    real(dp), allocatable :: tmp_node_m(:, :)
    !> Temporary workspace for local coefficients of each node. Dimension is
    !!      ((pl+1)**2, nsph)
    real(dp), allocatable :: tmp_node_l(:, :)
    !> Temporary workspace for grid values of each sphere. Dimension is
    !!      (ngrid, nsph).
    real(dp), allocatable :: tmp_grid(:, :)
end type ddx_type

contains

!> Initialize ddX input with a full set of parameters
!!
!! @param[in] nsph: Number of atoms. n > 0.
!! @param[in] charge: Charges of atoms. Dimension is `(n)`
!! @param[in] x: \f$ x \f$ coordinates of atoms. Dimension is `(n)`
!! @param[in] y: \f$ y \f$ coordinates of atoms. Dimension is `(n)`
!! @param[in] z: \f$ z \f$ coordinates of atoms. Dimension is `(n)`
!! @param[in] rvdw: Van-der-Waals radii of atoms. Dimension is `(n)`
!! @param[in] model: Choose model: 1 for COSMO, 2 for PCM and 3 for LPB
!! @param[in] lmax: Maximal degree of modeling spherical harmonics. `lmax` >= 0
!! @param[inout] ngrid: Approximate number of Lebedev grid points on input and
!!      actual number of grid points on exit. `ngrid` >= 0
!! @param[in] force: 1 if forces are required and 0 otherwise
!! @param[in] fmm: 1 to use FMM acceleration and 0 otherwise
!! @param[in] pm: Maximal degree of multipole spherical harmonics. Ignored in
!!      the case `fmm=0`. Value -1 means no far-field FMM interactions are
!!      computed. `pm` >= -1.
!! @param[in] pl: Maximal degree of local spherical harmonics. Ignored in
!!      the case `fmm=0`. Value -1 means no far-field FMM interactions are
!!      computed. `pl` >= -1.
!! @param[in] fmm_precompute: 1 to precompute FMM matrices, 0 to compute them
!!      on demand.
!! @param[in] iprint: Level of printing debug info.
!! @param[in] se: Shift of characteristic function. -1 for interior, 0 for
!!      centered and 1 for outer regularization
!! @param[in] eta: Regularization parameter
!! @param[in] eps: Relative dielectric permittivity
!! @param[in] kappa: Debye-H\"{u}ckel parameter
!! @param[in] itersolver: Iterative solver to be used. 1 for Jacobi iterative
!!      solver. Other solvers might be added later.
!! @param[in] tol: Relative error threshold for an iterative solver. tol must
!!      in range [1d-14, 1].
!! @param[in] maxiter: Maximum number of iterations for an iterative solver.
!!      maxiter > 0.
!! @param[in] ndiis: Number of extrapolation points for Jacobi/DIIS solver.
!!      ndiis >= 0.
!! @param[inout] nproc: Number of OpenMP threads to be used where applicable.
!!      nproc >= 0. If input nproc=0 then on exit this value contains actual
!!      number of threads that ddX will use. If OpenMP support in ddX is
!!      disabled, only possible input values are 0 or 1 and both inputs lead
!!      to the same output nproc=1 since the library is not parallel.
!! @param[out] ddx_data: Object containing all inputs
!! @param[out] info: flag of succesfull exit
!!      = 0: Succesfull exit
!!      < 0: If info=-i then i-th argument had an illegal value
!!      > 0: Allocation of a buffer for the output ddx_data failed
subroutine ddinit(nsph, charge, x, y, z, rvdw, model, lmax, ngrid, force, &
        & fmm, pm, pl, fmm_precompute, iprint, se, eta, eps, kappa, &
        & itersolver, tol, maxiter, ndiis, nproc, ddx_data, info)
    ! Inputs
    integer, intent(in) :: nsph, model, lmax, force, fmm, pm, pl, &
        & fmm_precompute, iprint, itersolver, maxiter, ndiis
    real(dp), intent(in):: charge(nsph), x(nsph), y(nsph), z(nsph), &
        & rvdw(nsph), se, eta, eps, kappa, tol
    ! Output
    type(ddx_type), target, intent(out) :: ddx_data
    integer, intent(out) :: info
    ! Inouts
    integer, intent(inout) :: ngrid, nproc
    ! Local variables
    integer :: istatus, i, indi, j, ii, inear, jnear, igrid, jgrid, isph, &
        & jsph, lnl, l, indl, ithread, k, n, indjn
    real(dp) :: v(3), maxv, ssqv, vv, t, swthr, fac, r, start_time, &
        & finish_time, tmp1
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    real(dp), allocatable :: sphcoo(:, :)
    double precision, external :: dnrm2
    integer :: iwork, jwork, lwork, old_lwork, nngmax
    integer, allocatable :: work(:, :), tmp_work(:, :), tmp_nl(:)
    logical :: use_omp
    ! Pointers for OMP, that can not use fields of a derived type in a depend
    ! clause. This is done to compute all the ddx_data fields in a parallel way
    real(dp), pointer :: p_vscales(:), p_v4pi2lp1(:), p_vscales_rel(:), &
        & p_vfact(:), p_vcnk(:), p_m2l_ztranslate_coef(:, :, :), &
        & p_m2l_ztranslate_adj_coef(:, :, :), p_cgrid(:, :), p_wgrid(:), &
        & p_vgrid(:, :), p_vwgrid(:, :), p_l2grid(:, :), p_tmp_vplm(:, :), &
        & p_tmp_vcos(:, :), p_tmp_vsin(:, :), p_fi(:, :), p_ui(:, :), &
        & p_zi(:, :, :)
    ! Reset info
    info = 0
    !!! Check input parameters
    ! Number of atoms
    if (nsph .le. 0) then
        write(*, "(A)") "ddinit: wrong value of a parameter `nsph`"
        info = -1
        return
    end if
    ddx_data % nsph = nsph
    ! Allocation of charge, csph and rsph
    allocate(ddx_data % charge(nsph), ddx_data % csph(3, nsph), &
        & ddx_data % rsph(nsph), stat=istatus)
    if (istatus .ne. 0) then
        write(*, "(A)") "ddinit: allocation (charge, csph, rsph) failed!"
        info = 1
        return
    end if
    ! No checks for NaNs and Infs in input coordinates
    ddx_data % charge = charge
    ddx_data % csph(1, :) = x
    ddx_data % csph(2, :) = y
    ddx_data % csph(3, :) = z
    ddx_data % rsph = rvdw
    ! Model, 1=COSMO, 2=PCM, 3=LPB
    if ((model .lt. 1) .or. (model .gt. 3)) then
        write(*, "(A)") "ddinit: wrong value of a parameter `model`"
        info = -7
        return
    end if
    ddx_data % model = model
    ! Degree of modeling spherical harmonics
    if (lmax .lt. 0) then
        write(*, "(A)") "ddinit: wrong value of a parameter `lmax`"
        info = -8
        return
    end if
    ddx_data % lmax = lmax
    ! Total number of modeling spherical harmonics
    ddx_data % nbasis = (lmax+1)**2
    ! Total number of modeling degrees of freedom
    ddx_data % n = nsph * ddx_data % nbasis
    ! Number of Lebedev grid points
    if (ngrid .lt. 0) then
        write(*, "(A)") "ddinit: wrong value of a parameter `ngrid`"
        info = -9
        return
    end if
    ! Actual value of ngrid will be calculated later
    ! Check if forces are needed
    if ((force .lt. 0) .or. (force .gt. 1)) then
        write(*, "(A)") "ddinit: wrong value of a parameter `force`"
        info = -10
        return
    end if
    ddx_data % force = force
    ! Check if FMM-acceleration is needed
    if ((fmm .lt. 0) .or. (fmm .gt. 1)) then
        write(*, "(A)") "ddinit: wrong value of a parameter `fmm`"
        info = -11
        return
    end if
    ddx_data % fmm = fmm
    ! Printing flag
    if (iprint .lt. 0) then
        write(*, "(A)") "ddinit: wrong value of a parameter `iprint`"
        info = -15
        return
    end if
    ddx_data % iprint = iprint
    ! Shift of a regularization
    if ((se .lt. -one) .or. (se .gt. one)) then
        write(*, "(A)") "ddinit: wrong value of a parameter `se`"
        info = -16
        return
    end if
    ddx_data % se = se
    ! Regularization parameter
    if ((eta .lt. zero) .or. (eta .gt. one)) then
        write(*, "(A)") "ddinit: wrong value of a parameter `eta`"
        info = -17
        return
    end if
    ddx_data % eta = eta
    ! Relative dielectric permittivity
    if (eps .lt. zero) then
        write(*, "(A)") "ddinit: wrong value of a parameter `eps`"
        info = -18
        return
    end if
    ddx_data % eps = eps
    ! Debye-H\"{u}ckel parameter (only used in ddLPB)
    if ((model .eq. 3) .and. (kappa .lt. zero)) then
        write(*, "(A)") "ddinit: wrong value of a parameter `kappa`"
        info = -19
        return
    end if
    ddx_data % kappa = kappa
    ! Iterative solver, 1=Jacobi
    if (itersolver .ne. 1) then
        write(*, "(A)") "ddinit: wrong value of a parameter `tol`"
        info = -20
        return
    end if
    ddx_data % itersolver = itersolver
    ! Relative threshold for an iterative solver
    if ((tol .lt. 1d-14) .or. (tol .gt. one)) then
        write(*, "(A)") "ddinit: wrong value of a parameter `tol`"
        info = -21
        return
    end if
    ddx_data % tol = tol
    ! Maximum number of iterations
    if (maxiter .le. 0) then
        write(*, "(A)") "ddinit: wrong value of a parameter `maxiter`"
        info = -22
        return
    end if
    ddx_data % maxiter = maxiter
    ! Number of DIIS extrapolation points (ndiis=25 works)
    if (ndiis .lt. 0) then
        write(*, "(A)") "ddinit: wrong value of a parameter `ndiis`"
        info = -23
        return
    end if
    ddx_data % ndiis = ndiis
    ! Number of OpenMP threads to be used
    ! As of now only sequential version with nproc=1 is supported, providing
    ! nproc=0 means it is up to ddX to decide on parallelism which will set
    ! it to nproc=1
    if (nproc .ne. 1 .and. nproc .ne. 0) then
        write(*, "(A)") "ddinit: wrong value of parameter `nproc`"
        info = -24
        return
    end if
    ! Only 1 thread is used (due to disabled OpenMP)
    nproc = 1
    ddx_data % nproc = nproc
    ! Set FMM parameters if FMM is needed
    if (ddx_data % fmm .eq. 1) then
        ! Maximal degree of multipole spherical harmonics. Value -1 means no 
        ! far-field interactions are to be computed, only near-field
        ! interactions are taken into account.
        if (pm .lt. -1) then
            write(*, "(A)") "ddinit: wrong value of a parameter `pm`"
            info = -12
            return
        end if
        ! Maximal degree of local spherical harmonics. Value -1 means no 
        ! far-field interactions are to be computed, only near-field
        ! interactions are taken into account.
        if (pl .lt. -1) then
            write(*, "(A)") "ddinit: wrong value of a parameter `pl`"
            info = -13
            return
        end if
        ! If far-field interactions are to be ignored
        if ((pl .eq. -1) .or. (pm .eq. -1)) then
            ddx_data % pm = -1
            ddx_data % pl = -1
        ! If far-field interactions are to be taken into account
        else
            ddx_data % pm = pm
            ddx_data % pl = pl
        end if
        ! Whether FMM translation matrices are to be computed. This might
        ! require lots of RAM.
        if ((fmm_precompute .lt. 0) .or. (fmm_precompute .gt. 1)) then
            write(*, "(A)") &
                & "ddinit: wrong value of a parameter `fmm_precompute`"
            info = -14
            return
        end if
        ddx_data % fmm_precompute = fmm_precompute
    else
        ! These values are ignored if fmm flag is 0
        ddx_data % pm = -2
        ddx_data % pl = -2
        ddx_data % fmm_precompute = -1
    end if
    !!! Generate constants and auxiliary arrays
    ! Compute sizes of auxiliary arrays for `fmm=0`
    if (ddx_data % fmm .eq. 0) then
        ddx_data % dmax = ddx_data % lmax
        ddx_data % vgrid_dmax = ddx_data % lmax
    ! Compute sizes of arrays if `fmm=1`
    else
        ! If forces are required then we need M2P of degree lmax+1 for
        ! near-field analytical gradients
        if (force .eq. 1) then
            ddx_data % dmax = max(ddx_data % pm+ddx_data % pl, &
                & ddx_data % lmax+1)
            ddx_data % m2p_lmax = ddx_data % lmax + 1
        else
            ddx_data % dmax = max(ddx_data % pm+ddx_data % pl, ddx_data % lmax)
            ddx_data % m2p_lmax = ddx_data % lmax
        end if
        ddx_data % vgrid_dmax = max(ddx_data % pl, ddx_data % lmax)
        ddx_data % m2p_nbasis = (ddx_data % m2p_lmax+1) ** 2
        ddx_data % grad_nbasis = (ddx_data % lmax+2) ** 2
    end if
    ! Compute sizes of vgrid, vfact and vscales that are used for both 
    ddx_data % vgrid_nbasis = (ddx_data % vgrid_dmax+1) ** 2
    ddx_data % nfact = 2*ddx_data % dmax+1
    ddx_data % nscales = (ddx_data % dmax+1) ** 2
    ! Get nearest number of Lebedev grid points
    igrid = 0
    inear = 100000
    do i = 1, nllg
        jnear = abs(ng0(i)-ngrid)
        if (jnear .lt. inear) then
            inear = jnear
            igrid = i
        end if
    end do
    ! Update inout parameter `ngrid` also
    ngrid = ng0(igrid)
    ddx_data % ngrid = ngrid
    !! Now all constants are initialized, continue with memory allocations
    ! Allocate space for scaling factors of spherical harmonics
    allocate(ddx_data % vscales(ddx_data % nscales), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [2] allocation failed!"
        info = 2
        return
    end if
    p_vscales => ddx_data % vscales
    allocate(ddx_data % v4pi2lp1(ddx_data % dmax+1), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [2] allocation failed!"
        info = 2
        return
    end if
    p_v4pi2lp1 => ddx_data % v4pi2lp1
    allocate(ddx_data % vscales_rel(ddx_data % nscales), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [2] allocation failed!"
        info = 2
        return
    end if
    p_vscales_rel => ddx_data % vscales_rel
    ! Allocate square roots of factorials
    allocate(ddx_data % vfact(ddx_data % nfact), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [3] allocation failed!"
        info = 3
        return
    end if
    p_vfact => ddx_data % vfact
    ! Allocate space for Lebedev grid coordinates and weights
    allocate(ddx_data % cgrid(3, ddx_data % ngrid), &
        & ddx_data % wgrid(ddx_data % ngrid), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [4] allocation failed!"
        info = 4
        return
    end if
    p_cgrid => ddx_data % cgrid
    p_wgrid => ddx_data % wgrid
    ! Allocate temporary arrays for spherical coordinates
    allocate(sphcoo(5, ddx_data % nproc), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [5] allocation failed!"
        info = 5
        return
    end if
    ! Allocate temporary arrays for associated Legendre polynomials and
    ! arrays of cosinuses and sinuses
    allocate(ddx_data % tmp_vplm(ddx_data % vgrid_nbasis, ddx_data % nproc), &
        & ddx_data % tmp_vcos(ddx_data % vgrid_dmax+1, ddx_data % nproc), &
        & ddx_data % tmp_vsin(ddx_data % vgrid_dmax+1, ddx_data % nproc), &
        & stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [5] allocation failed!"
        info = 5
        return
    end if
    p_tmp_vplm => ddx_data % tmp_vplm
    p_tmp_vcos => ddx_data % tmp_vcos
    p_tmp_vsin => ddx_data % tmp_vsin
    ! Allocate space for values of non-weighted and weighted spherical
    ! harmonics at Lebedev grid points
    allocate(ddx_data % vgrid(ddx_data % vgrid_nbasis, ddx_data % ngrid), &
        & ddx_data % vwgrid(ddx_data % vgrid_nbasis, ddx_data % ngrid), &
        & ddx_data % l2grid(ddx_data % vgrid_nbasis, ddx_data % ngrid), &
        & stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [6] allocation failed!"
        info = 6
        return
    end if
    p_vgrid => ddx_data % vgrid
    p_vwgrid => ddx_data % vwgrid
    p_l2grid => ddx_data % l2grid
    ! Allocate space for characteristic functions fi and ui
    allocate(ddx_data % fi(ddx_data % ngrid, ddx_data % nsph), &
        & ddx_data % ui(ddx_data % ngrid, ddx_data % nsph), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) 'ddinit : [9] allocation failed !'
        info = 9
        return
    end if
    ddx_data % fi = zero
    ddx_data % ui = zero
    p_fi => ddx_data % fi
    p_ui => ddx_data % ui
    ! Allocate space for force-related arrays
    if (force .eq. 1) then
        allocate(ddx_data % zi(3, ddx_data % ngrid, ddx_data % nsph), &
            & stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) 'ddinit : [10] allocation failed !'
            info = 10
            return
        endif
        ddx_data % zi = zero
        p_zi => ddx_data % zi
    end if
    ! Compute scaling factors of spherical harmonics
    call ylmscale(ddx_data % dmax, p_vscales, p_v4pi2lp1, p_vscales_rel)
    ! Compute square roots of factorials
    p_vfact(1) = 1
    do i = 2, ddx_data % nfact
        p_vfact(i) = p_vfact(i-1) * sqrt(dble(i-1))
    end do
    ! Get weights and coordinates of Lebedev grid points
    call llgrid(ddx_data % ngrid, p_wgrid, p_cgrid)
    ! Compute non-weighted and weighted spherical harmonics and the single
    ! layer operator at Lebedev grid points
    ithread = 1
    do igrid = 1, ddx_data % ngrid
        call ylmbas2(p_cgrid(:, igrid), sphcoo(:, ithread), &
            & ddx_data % vgrid_dmax, p_vscales, &
            & p_vgrid(:, igrid), p_tmp_vplm(:, ithread), &
            & p_tmp_vcos(:, ithread), p_tmp_vsin(:, ithread))
        p_vwgrid(:, igrid) = p_wgrid(igrid) * p_vgrid(:, igrid)
        do l = 0, ddx_data % pl
            indl = l*l + l + 1
            p_l2grid(indl-l:indl+l, igrid) = &
                & p_vgrid(indl-l:indl+l, igrid) / p_vscales(indl)**2
        end do
    end do
    ! Upper bound of switch region. Defines intersection criterion for spheres
    swthr = one + (ddx_data % se+one)*ddx_data % eta/two
    ! Build list of neighbours in CSR format
    nngmax = 1
    allocate(ddx_data % inl(nsph+1), ddx_data % nl(nsph*nngmax), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) 'ddinit : [8] allocation failed !'
        info = 8
        return
    end if
    i = 1
    lnl = 0
    do isph = 1, ddx_data % nsph
        ddx_data % inl(isph) = lnl + 1
        do jsph = 1, ddx_data % nsph
            if (isph .ne. jsph) then
                v = ddx_data % csph(:, isph)-ddx_data % csph(:, jsph)
                maxv = max(abs(v(1)), abs(v(2)), abs(v(3)))
                ssqv = (v(1)/maxv)**2 + (v(2)/maxv)**2 + (v(3)/maxv)**2
                vv = maxv * sqrt(ssqv)
                ! Take regularization parameter into account with respect to
                ! shift se. It is described properly by the upper bound of a
                ! switch region `swthr`.
                r = rvdw(isph) + swthr*rvdw(jsph)
                if (vv .le. r) then
                    ddx_data % nl(i) = jsph
                    i  = i + 1
                    lnl = lnl + 1
                    ! Extend ddx_data % nl if needed
                    if (i .gt. nsph*nngmax) then
                        allocate(tmp_nl(nsph*nngmax), stat=istatus)
                        if (istatus .ne. 0) then
                            !write(*, *) 'ddinit : [8] allocation failed!'
                            info = 8
                            return
                        end if
                        tmp_nl(1:nsph*nngmax) = ddx_data % nl(1:nsph*nngmax)
                        deallocate(ddx_data % nl, stat=istatus)
                        if (istatus .ne. 0) then
                            ! write(*, *) 'ddinit : [8] deallocation failed!'
                            info = 8
                            return
                        end if
                        nngmax = nngmax + 10
                        allocate(ddx_data % nl(nsph*nngmax), stat=istatus)
                        if (istatus .ne. 0) then
                            !write(*, *) 'ddinit : [8] allocation failed!'
                            info = 8
                            return
                        end if
                        ddx_data % nl(1:nsph*(nngmax-10)) = &
                            & tmp_nl(1:nsph*(nngmax-10))
                        deallocate(tmp_nl, stat=istatus)
                        if (istatus .ne. 0) then
                            ! write(*, *) 'ddinit : [8] deallocation failed!'
                            info = 8
                            return
                        end if
                    end if
                end if
            end if
        end do
    end do
    ddx_data % inl(nsph+1) = lnl+1
    ddx_data % nngmax = nngmax
    ! Build arrays fi, ui, zi
    do isph = 1, ddx_data % nsph
        do igrid = 1, ddx_data % ngrid
            ! Loop over neighbours of i-th sphere
            do inear = ddx_data % inl(isph), ddx_data % inl(isph+1)-1
                ! Neighbour's index
                jsph = ddx_data % nl(inear)
                ! Compute t_n^ij
                v = ddx_data % csph(:, isph) - ddx_data % csph(:, jsph) + &
                    & rvdw(isph)*ddx_data % cgrid(:, igrid)
                maxv = max(abs(v(1)), abs(v(2)), abs(v(3)))
                ssqv = (v(1)/maxv)**2 + (v(2)/maxv)**2 + (v(3)/maxv)**2
                vv = maxv * sqrt(ssqv)
                t = vv / ddx_data % rsph(jsph)
                ! Accumulate characteristic function \chi(t_n^ij)
                ddx_data % fi(igrid, isph) = ddx_data % fi(igrid, isph) + &
                    & fsw(t, ddx_data % se, ddx_data % eta)
                ! Check if gradients are required
                if (ddx_data % force .eq. 1) then
                    ! Check if t_n^ij belongs to switch region
                    if ((t .lt. swthr) .and. (t .gt. swthr-ddx_data % eta)) then
                        ! Accumulate gradient of characteristic function \chi
                        fac = dfsw(t, ddx_data % se, ddx_data % eta) / rvdw(jsph)
                        ddx_data % zi(:, igrid, isph) = &
                            & ddx_data % zi(:, igrid, isph) + fac/vv*v
                    end if
                end if
            enddo
            ! Compute characteristic function of a molecular surface ui
            if (ddx_data % fi(igrid, isph) .le. one) then
                ddx_data % ui(igrid, isph) = one - ddx_data % fi(igrid, isph)
            end if
        enddo
    enddo
    ! Debug printing
    if (iprint .ge. 4) then
        call ptcart('fi', ngrid, nsph, 0, ddx_data % fi)
        call ptcart('ui', ngrid, nsph, 0, ddx_data % ui)
    end if
    ! Build cavity array. At first get total count for each sphere
    allocate(ddx_data % ncav_sph(nsph), stat=istatus)
    if (istatus .ne. 0) then
        !write(*,*)'ddinit : [11] allocation failed!'
        info = 11
        return
    endif
    do isph = 1, nsph
        ddx_data % ncav_sph(isph) = 0
        ! Loop over integration points
        do i = 1, ngrid
            ! Positive contribution from integration point
            if (ddx_data % ui(i, isph) .gt. zero) then
                ddx_data % ncav_sph(isph) = ddx_data % ncav_sph(isph) + 1
            end if
        end do
    end do
    ddx_data % ncav = sum(ddx_data % ncav_sph)
    ! Allocate cavity array and CSR format for indexes of cavities
    allocate(ddx_data % ccav(3, ddx_data % ncav), ddx_data % icav_ia(nsph+1), &
        & ddx_data % icav_ja(ddx_data % ncav), stat=istatus)
    if (istatus .ne. 0) then
        !write(*,*)'ddinit : [11] allocation failed!'
        info = 11
        return
    endif
    ! Get actual cavity coordinates and indexes in CSR format
    ddx_data % icav_ia(1) = 1
    i = 1
    do isph = 1, nsph
        ddx_data % icav_ia(isph+1) = ddx_data % icav_ia(isph) + &
            & ddx_data % ncav_sph(isph)
        ! Loop over integration points
        do igrid = 1, ngrid
            ! Ignore zero contribution
            if (ddx_data % ui(igrid, isph) .eq. zero) cycle
            ! Store coordinates
            ddx_data % ccav(:, i) = ddx_data % csph(:, isph) + &
                & rvdw(isph)*ddx_data % cgrid(:, igrid)
            ! Store index
            ddx_data % icav_ja(i) = igrid
            ! Advance cavity array index
            i = i + 1
        end do
    end do
    ! Create preconditioner for PCM
    if (model .eq. 2) then
        allocate(ddx_data % rx_prc(ddx_data % nbasis, ddx_data % nbasis, nsph), &
            & stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [12] allocation failed!'
            info = 12
            return
        endif
        call mkprec(ddx_data)
    end if
    ! Again some debug output
    1100  format(t3,i8,3f14.6)
    if (iprint .ge. 4) then
        write(6, *) '   external cavity points:'
        do ii = 1, ddx_data % ncav
            write(6,1100) ii, ddx_data % ccav(:, ii)
        end do
        write(6, *)
    end if
    !! Prepare FMM structures if needed
    if (fmm .eq. 1) then
        ! Allocate space for a cluster tree
        allocate(ddx_data % order(nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [13] allocation failed!'
            info = 13
            return
        endif
        ddx_data % nclusters = 2*nsph - 1
        allocate(ddx_data % cluster(2, ddx_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [14] allocation failed!'
            info = 14
            return
        endif
        allocate(ddx_data % children(2, ddx_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [15] allocation failed!'
            info = 15
            return
        endif
        allocate(ddx_data % parent(ddx_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [16] allocation failed!'
            info = 16
            return
        endif
        allocate(ddx_data % cnode(3, ddx_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [17] allocation failed!'
            info = 17
            return
        endif
        allocate(ddx_data % rnode(ddx_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [18] allocation failed!'
            info = 18
            return
        endif
        allocate(ddx_data % snode(ddx_data % nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [19] allocation failed!'
            info = 19
            return
        endif
        allocate(ddx_data % nfar(ddx_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [20] allocation failed!'
            info = 20
            return
        endif
        allocate(ddx_data % nnear(ddx_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [21] allocation failed!'
            info = 21
            return
        endif
        ! Get the tree
        call tree_rib_build(nsph, ddx_data % csph, ddx_data % rsph, &
            & ddx_data % order, ddx_data % cluster, ddx_data % children, &
            & ddx_data % parent, ddx_data % cnode, ddx_data % rnode, &
            & ddx_data % snode)
        ! Get number of far and near admissible pairs
        iwork = 0
        jwork = 1
        lwork = 1
        allocate(work(3, lwork), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Could not allocate space for workspace"
            !write(*, *) "Attempted size is", lwork*24, " bytes"
            !write(*,*)'ddinit : [22] allocation failed!'
            info = 22
            return
        end if
        do while (iwork .le. jwork)
            allocate(tmp_work(3, lwork), stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Could not allocate space for workspace"
                !write(*, *) "Attempted size is", lwork*24, " bytes"
                !write(*,*)'ddinit : [23] allocation failed!'
                info = 23
                return
            end if
            tmp_work = work
            deallocate(work, stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Error in deallocation of workspace"
                info = 24
                return
            end if
            old_lwork = lwork
            lwork = old_lwork + 1000*nsph
            allocate(work(3, lwork), stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Could not allocate space for workspace"
                !write(*, *) "Attempted size is", lwork*24, " bytes"
                info = 25
                return
            end if
            work(:, 1:old_lwork) = tmp_work
            deallocate(tmp_work, stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Error in deallocation of workspace"
                info = 26
                return
            end if
            call tree_get_farnear_work(ddx_data % nclusters, &
                & ddx_data % children, ddx_data % cnode, &
                & ddx_data % rnode, lwork, iwork, jwork, work, &
                & ddx_data % nnfar, ddx_data % nfar, ddx_data % nnnear, &
                & ddx_data % nnear)
        end do
        allocate(ddx_data % sfar(ddx_data % nclusters+1), &
            & ddx_data % snear(ddx_data % nclusters+1), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation"
            info = 27
            return
        end if
        allocate(ddx_data % far(ddx_data % nnfar), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation"
            info = 28
            return
        end if
        allocate(ddx_data % near(ddx_data % nnnear), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation"
            info = 29
            return
        end if
        call tree_get_farnear(jwork, lwork, work, ddx_data % nclusters, &
            & ddx_data % nnfar, ddx_data % nfar, ddx_data % sfar, ddx_data % far, &
            & ddx_data % nnnear, ddx_data % nnear, ddx_data % snear, &
            & ddx_data % near)
        deallocate(work, stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in deallocation of workspace"
            info = 30
            return
        end if
        ! Allocate square roots of combinatorial numbers C_n^k
        allocate(ddx_data % vcnk((2*ddx_data % dmax+1)*(ddx_data % dmax+1)), &
            & stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "ddinit: [3] allocation failed!"
            info = 3
            return
        end if
        p_vcnk => ddx_data % vcnk
        ! Allocate space for M2L OZ translation coefficients
        allocate(ddx_data % m2l_ztranslate_coef(&
            & (ddx_data % pm+1), (ddx_data % pl+1), (ddx_data % pl+1)), &
            & stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "ddinit: [3] allocation failed!"
            info = 3
            return
        end if
        p_m2l_ztranslate_coef => ddx_data % m2l_ztranslate_coef
        ! Allocate space for adjoint M2L OZ translation coefficients
        allocate(ddx_data % m2l_ztranslate_adj_coef(&
            & (ddx_data % pl+1), (ddx_data % pl+1), (ddx_data % pm+1)), &
            & stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "ddinit: [3] allocation failed!"
            info = 3
            return
        end if
        p_m2l_ztranslate_adj_coef => ddx_data % m2l_ztranslate_adj_coef
        ! Compute combinatorial numbers C_n^k and M2L OZ translate coefficients
        call fmm_constants(ddx_data % dmax, ddx_data % pm, ddx_data % pl, &
            & p_vcnk, p_m2l_ztranslate_coef, p_m2l_ztranslate_adj_coef)
        allocate(ddx_data % tmp_sph(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of a temporary"
            info = 38
            return
        end if
        allocate(ddx_data % tmp_sph2(ddx_data % grad_nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of a temporary"
            info = 38
            return
        end if
        allocate(ddx_data % tmp_sph_grad(ddx_data % grad_nbasis, 3, nsph), &
            & stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of a temporary"
            info = 38
            return
        end if
        allocate(ddx_data % tmp_sph_l((ddx_data % pl+1)**2, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of a temporary"
            info = 38
            return
        end if
        allocate(ddx_data % tmp_sph_l_grad((ddx_data % pl+1)**2, 3, nsph), &
            & stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of a temporary"
            info = 38
            return
        end if
        allocate(ddx_data % tmp_node_m((ddx_data % pm+1)**2, &
            & ddx_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of a temporary"
            info = 30
            return
        end if
        allocate(ddx_data % tmp_node_l((ddx_data % pl+1)**2, &
            & ddx_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of a temporary"
            info = 30
            return
        end if
        if (fmm_precompute .eq. 1) then
            ddx_data % m2m_ztranslate_mat_size = (pm+1)*(pm+2)*(pm+3)/6
            ddx_data % m2m_reflect_mat_size = (pm+1)*(2*pm+1)*(2*pm+3)/3
            ddx_data % l2l_ztranslate_mat_size = (pl+1)*(pl+2)*(pl+3)/6
            ddx_data % l2l_reflect_mat_size = (pl+1)*(2*pl+1)*(2*pl+3)/3
            ddx_data % m2l_ztranslate_mat_size = (min(pm,pl)+1) * &
                & (min(pm,pl)+2) * (3*max(pm,pl)+3-min(pm,pl)) / 6
            ddx_data % m2l_reflect_mat_size = max( &
                & ddx_data % m2m_reflect_mat_size, &
                & ddx_data % l2l_reflect_mat_size)
            allocate(ddx_data % m2m_ztranslate_mat( &
                & ddx_data % m2m_ztranslate_mat_size, ddx_data % nclusters), &
                & stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Error in allocation of M2M matrices"
                info = 31
                return
            end if
            allocate(ddx_data % m2m_reflect_mat( &
                & ddx_data % m2m_reflect_mat_size, ddx_data % nclusters), &
                & stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Error in allocation of M2M matrices"
                info = 32
                return
            end if
            allocate(ddx_data % l2l_ztranslate_mat( &
                & ddx_data % l2l_ztranslate_mat_size, ddx_data % nclusters), &
                & stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Error in allocation of L2L matrices"
                info = 33
                return
            end if
            allocate(ddx_data % l2l_reflect_mat( &
                & ddx_data % l2l_reflect_mat_size, ddx_data % nclusters), &
                & stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Error in allocation of L2L matrices"
                info = 34
                return
            end if
            allocate(ddx_data % m2l_ztranslate_mat( &
                & ddx_data % m2l_ztranslate_mat_size, ddx_data % nnfar), &
                & stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Error in allocation of M2L matrices"
                info = 35
                return
            end if
            allocate(ddx_data % m2l_reflect_mat( &
                & ddx_data % m2l_reflect_mat_size, ddx_data % nnfar), &
                & stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Error in allocation of M2L matrices"
                info = 36
                return
            end if
            ! Precompute M2M, L2L and M2L matrices
            call tree_m2m_reflection_get_mat(ddx_data)
            call tree_l2l_reflection_get_mat(ddx_data)
            call tree_m2l_reflection_get_mat(ddx_data)
            ! Get near-field M2P data
            ddx_data % nnear_m2p = 0
            do isph = 1, nsph
                ddx_data % nnear_m2p = ddx_data % nnear_m2p + &
                    & ddx_data % ncav_sph(isph)*ddx_data % nnear( &
                    & ddx_data % snode(isph))
            end do
            allocate(ddx_data % m2p_mat(ddx_data % m2p_nbasis, &
                & ddx_data % nnear_m2p), stat=istatus)
            if (istatus .ne. 0) then
                !write(*, *) "Error in allocation of M2P matrices"
                info = 37
                return
            end if
            ! Precompute M2P matrices
            call tree_m2p_get_mat(ddx_data)
        end if
    end if
    !! Per-model allocations
    ! COSMO model
    if (model .eq. 1) then
        allocate(ddx_data % phi_grid(ddx_data % ngrid, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % phi(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % xs(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % s(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % sgrid(ddx_data % ngrid, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % zeta(ddx_data % ncav), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
    ! PCM model
    else if (model .eq. 2) then
        allocate(ddx_data % phi_grid(ddx_data % ngrid, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % phi(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % phiinf(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % phieps(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % xs(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % s(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % sgrid(ddx_data % ngrid, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % y(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % ygrid(ddx_data % ngrid, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % g(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % q(ddx_data % nbasis, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % qgrid(ddx_data % ngrid, nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
        allocate(ddx_data % zeta(ddx_data % ncav), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) "Error in allocation of M2P matrices"
            info = 38
            return
        end if
    ! LPB model
    else if (model .eq. 3) then
    end if
    allocate(ddx_data % tmp_grid(ddx_data % ngrid, nsph), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "Error in allocation of a temporary"
        info = 38
        return
    end if
    !! Debug outputs
    ! Some format data that I will throw away as soon as this code works
    1000 format(t3,'neighbours of sphere ',i6)
    1010 format(t5,12i6)
    ! Debug printing
    if (iprint .ge. 4) then
        write(*, "(A)") '    inl:'
        write(*, "(10i6)") ddx_data % inl(1:nsph+1)
        write(*,*)
        do isph = 1, ddx_data % nsph
            write(6,1000) isph
            write(6,1010) ddx_data % nl(ddx_data % inl(isph): &
                & ddx_data % inl(isph+1)-1)
        end do
        write(6,*)
    end if
end subroutine ddinit

!> Read configuration from a file
!!
!! @param[in] fname: Filename containing all the required info
!! @param[out] ddx_data: Object containing all inputs
!! @param[out] info: flag of succesfull exit
!!      = 0: Succesfull exit
!!      < 0: If info=-i then i-th argument had an illegal value
!!      > 0: Allocation of a buffer for the output ddx_data failed
subroutine ddfromfile(fname, ddx_data, info)
    ! Input
    character(len=*), intent(in) :: fname
    ! Outputs
    type(ddx_type), intent(out) :: ddx_data
    integer, intent(out) :: info
    ! Local variables
    integer :: iprint, nproc, model, lmax, ngrid, force, fmm, pm, pl, &
        & fmm_precompute, nsph, i, itersolver, maxiter, ndiis, istatus
    real(dp) :: eps, se, eta, kappa, tol
    real(dp), allocatable :: charge(:), x(:), y(:), z(:), rvdw(:)
    !! Read all the parameters from the file
    ! Open a configuration file
    open(unit=100, file=fname, form='formatted', access='sequential')
    ! Printing flag
    read(100, *) iprint
    if(iprint .lt. 0) then
        write(*, "(3A)") "Error on the 1st line of a config file ", fname, &
            & ": `iprint` must be a non-negative integer value."
        stop 1
    end if
    ! Number of OpenMP threads to be used
    read(100, *) nproc
    if(nproc .lt. 0) then
        write(*, "(3A)") "Error on the 2nd line of a config file ", fname, &
            & ": `nproc` must be a positive integer value."
        stop 1
    end if
    ! Model to be used: 1 for COSMO, 2 for PCM and 3 for LPB
    read(100, *) model
    if((model .lt. 1) .or. (model .gt. 3)) then
        write(*, "(3A)") "Error on the 3rd line of a config file ", fname, &
            & ": `model` must be an integer of a value 1, 2 or 3."
        stop 1
    end if
    ! Max degree of modeling spherical harmonics
    read(100, *) lmax
    if(lmax .lt. 0) then
        write(*, "(3A)") "Error on the 4th line of a config file ", fname, &
            & ": `lmax` must be a non-negative integer value."
        stop 1
    end if
    ! Approximate number of Lebedev points
    read(100, *) ngrid
    if(ngrid .lt. 0) then
        write(*, "(3A)") "Error on the 5th line of a config file ", fname, &
            & ": `ngrid` must be a non-negative integer value."
        stop 1
    end if
    ! Dielectric permittivity constant of the solvent
    read(100, *) eps
    if(eps .lt. zero) then
        write(*, "(3A)") "Error on the 6th line of a config file ", fname, &
            & ": `eps` must be a non-negative floating point value."
        stop 1
    end if
    ! Shift of the regularized characteristic function
    read(100, *) se
    if((se .lt. -one) .or. (se .gt. one)) then
        write(*, "(3A)") "Error on the 7th line of a config file ", fname, &
            & ": `se` must be a floating point value in a range [-1, 1]."
        stop 1
    end if
    ! Regularization parameter
    read(100, *) eta
    if((eta .lt. zero) .or. (eta .gt. one)) then
        write(*, "(3A)") "Error on the 8th line of a config file ", fname, &
            & ": `eta` must be a floating point value in a range [0, 1]."
        stop 1
    end if
    ! Debye H\"{u}ckel parameter
    read(100, *) kappa
    if(kappa .lt. zero) then
        write(*, "(3A)") "Error on the 9th line of a config file ", fname, &
            & ": `kappa` must be a non-negative floating point value."
        stop 1
    end if
    ! Iterative solver of choice. Only one is supported as of now.
    read(100, *) itersolver
    if(itersolver .ne. 1) then
        write(*, "(3A)") "Error on the 10th line of a config file ", fname, &
            & ": `itersolver` must be an integer value of a value 1."
        stop 1
    end if
    ! Relative convergence threshold for the iterative solver
    read(100, *) tol
    if((tol .lt. 1d-14) .or. (tol .gt. one)) then
        write(*, "(3A)") "Error on the 11th line of a config file ", fname, &
            & ": `tol` must be a floating point value in a range [1d-14, 1]."
        stop 1
    end if
    ! Maximum number of iterations for the iterative solver
    read(100, *) maxiter
    if((maxiter .le. 0)) then
        write(*, "(3A)") "Error on the 12th line of a config file ", fname, &
            & ": `maxiter` must be a positive integer value."
        stop 1
    end if
    ! Number of extrapolation points for Jacobi/DIIS solver
    read(100, *) ndiis
    if((ndiis .lt. 0)) then
        write(*, "(3A)") "Error on the 13th line of a config file ", fname, &
            & ": `ndiis` must be a non-negative integer value."
        stop 1
    end if
    ! Whether to compute (1) or not (0) forces as analytical gradients
    read(100, *) force
    if((force .lt. 0) .or. (force .gt. 1)) then
        write(*, "(3A)") "Error on the 14th line of a config file ", fname, &
            & ": `force` must be an integer value of a value 0 or 1."
        stop 1
    end if
    ! Whether to use (1) or not (0) the FMM to accelerate computations
    read(100, *) fmm
    if((fmm .lt. 0) .or. (fmm .gt. 1)) then
        write(*, "(3A)") "Error on the 15th line of a config file ", fname, &
            & ": `fmm` must be an integer value of a value 0 or 1."
        stop 1
    end if
    ! Max degree of multipole spherical harmonics for the FMM
    read(100, *) pm
    if(pm .lt. 0) then
        write(*, "(3A)") "Error on the 16th line of a config file ", fname, &
            & ": `pm` must be a non-negative integer value."
        stop 1
    end if
    ! Max degree of local spherical harmonics for the FMM
    read(100, *) pl
    if(pl .lt. 0) then
        write(*, "(3A)") "Error on the 17th line of a config file ", fname, &
            & ": `pl` must be a non-negative integer value."
        stop 1
    end if
    ! Whether to precompute (1) or obtain on demand (0) the FMM translations
    read(100, *) fmm_precompute
    if((fmm_precompute .lt. 0) .or. (fmm_precompute .gt. 1)) then
        write(*, "(3A)") "Error on the 18th line of a config file ", fname, &
            & ": `fmm_precompute` must be an integer value of a value 0 or 1."
        stop 1
    end if
    ! Number of input spheres
    read(100, *) nsph
    if(nsph .le. 0) then
        write(*, "(3A)") "Error on the 19th line of a config file ", fname, &
            & ": `nsph` must be a positive integer value."
        stop 1
    end if
    ! Coordinates, radii and charges
    allocate(charge(nsph), x(nsph), y(nsph), z(nsph), rvdw(nsph), stat=istatus)
    if(istatus .ne. 0) then
        write(*, "(2A)") "Could not allocate space for coordinates, radii ", &
            & "and charges of atoms"
        stop 1
    end if
    do i = 1, nsph
        read(100, *) charge(i), x(i), y(i), z(i), rvdw(i)
    end do
    ! Finish reading
    close(100)
    !! Convert Angstrom input into Bohr units
    x = x * tobohr
    y = y * tobohr
    z = z * tobohr
    rvdw = rvdw * tobohr
    !! Initialize ddx_data object
    call ddinit(nsph, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, &
        & pm, pl, fmm_precompute, iprint, se, eta, eps, kappa, itersolver, &
        & tol, maxiter, ndiis, nproc, ddx_data, info)
    !! Clean local temporary data
    deallocate(charge, x, y, z, rvdw, stat=istatus)
    if(istatus .ne. 0) then
        write(*, "(2A)") "Could not deallocate space for coordinates, ", &
            & "radii and charges of atoms"
        stop 1
    end if
end subroutine ddfromfile

!> Deallocate object with corresponding data
!!
!! @param[inout] ddx_data: object to deallocate
subroutine ddfree(ddx_data)
    ! Input/output
    type(ddx_type), intent(inout) :: ddx_data
    ! Local variables
    integer :: istatus
    if (allocated(ddx_data % charge)) then
        deallocate(ddx_data % charge, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [1] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % csph)) then
        deallocate(ddx_data % csph, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [1] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % rsph)) then
        deallocate(ddx_data % rsph, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [2] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % vscales)) then
        deallocate(ddx_data % vscales, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [3] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % v4pi2lp1)) then
        deallocate(ddx_data % v4pi2lp1, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [3] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % vscales_rel)) then
        deallocate(ddx_data % vscales_rel, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [3] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % vfact)) then
        deallocate(ddx_data % vfact, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [4] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % vcnk)) then
        deallocate(ddx_data % vcnk, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [4] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % m2l_ztranslate_coef)) then
        deallocate(ddx_data % m2l_ztranslate_coef, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [4] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % m2l_ztranslate_adj_coef)) then
        deallocate(ddx_data % m2l_ztranslate_adj_coef, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [4] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % cgrid)) then
        deallocate(ddx_data % cgrid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [5] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % wgrid)) then
        deallocate(ddx_data % wgrid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [6] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % vgrid)) then
        deallocate(ddx_data % vgrid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [7] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % vwgrid)) then
        deallocate(ddx_data % vwgrid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [8] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % l2grid)) then
        deallocate(ddx_data % l2grid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [8+] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % inl)) then
        deallocate(ddx_data % inl, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [9] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % nl)) then
        deallocate(ddx_data % nl, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [10] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % fi)) then
        deallocate(ddx_data % fi, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [11] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % ui)) then
        deallocate(ddx_data % ui, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [12] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % zi)) then
        deallocate(ddx_data % zi, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [13] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % ncav_sph)) then
        deallocate(ddx_data % ncav_sph, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [14] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % icav_ia)) then
        deallocate(ddx_data % icav_ia, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [14] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % icav_ja)) then
        deallocate(ddx_data % icav_ja, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [14] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % ccav)) then
        deallocate(ddx_data % ccav, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [14] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % rx_prc)) then
        deallocate(ddx_data % rx_prc, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [15] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(ddx_data % order)) then
        deallocate(ddx_data % order, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [16] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % cluster)) then
        deallocate(ddx_data % cluster, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [17] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % children)) then
        deallocate(ddx_data % children, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [18] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % parent)) then
        deallocate(ddx_data % parent, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [19] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % cnode)) then
        deallocate(ddx_data % cnode, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [20] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % rnode)) then
        deallocate(ddx_data % rnode, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [21] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % snode)) then
        deallocate(ddx_data % snode, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [22] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % nfar)) then
        deallocate(ddx_data % nfar, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [23] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % nnear)) then
        deallocate(ddx_data % nnear, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [24] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % sfar)) then
        deallocate(ddx_data % sfar, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [25] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % snear)) then
        deallocate(ddx_data % snear, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [26] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % far)) then
        deallocate(ddx_data % far, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [27] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % near)) then
        deallocate(ddx_data % near, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [28] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % tmp_sph)) then
        deallocate(ddx_data % tmp_sph, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [28] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % tmp_node_m)) then
        deallocate(ddx_data % tmp_node_m, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [28] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % tmp_node_l)) then
        deallocate(ddx_data % tmp_node_l, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [28] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % tmp_grid)) then
        deallocate(ddx_data % tmp_grid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [28] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % m2m_ztranslate_mat)) then
        deallocate(ddx_data % m2m_ztranslate_mat, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [29] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % m2m_reflect_mat)) then
        deallocate(ddx_data % m2m_reflect_mat, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [30] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % l2l_ztranslate_mat)) then
        deallocate(ddx_data % l2l_ztranslate_mat, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [31] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % l2l_reflect_mat)) then
        deallocate(ddx_data % l2l_reflect_mat, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [32] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % m2l_ztranslate_mat)) then
        deallocate(ddx_data % m2l_ztranslate_mat, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [33] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % m2l_reflect_mat)) then
        deallocate(ddx_data % m2l_reflect_mat, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [34] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % m2p_mat)) then
        deallocate(ddx_data % m2p_mat, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [35] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % phi)) then
        deallocate(ddx_data % phi, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [36] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % phiinf)) then
        deallocate(ddx_data % phiinf, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [36] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % phieps)) then
        deallocate(ddx_data % phieps, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [36] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % xs)) then
        deallocate(ddx_data % xs, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [36] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % s)) then
        deallocate(ddx_data % s, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [36] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % y)) then
        deallocate(ddx_data % y, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [36] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % ygrid)) then
        deallocate(ddx_data % ygrid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [36] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % g)) then
        deallocate(ddx_data % g, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [36] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(ddx_data % tmp_grid)) then
        deallocate(ddx_data % tmp_grid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [36] deallocation failed!"
            stop 1
        endif
    end if
end subroutine ddfree

!> Print array of spherical harmonics
!!
!! Prints (nbasis, ncol) array
!!
!! @param[in] label: Label to print
!! @param[in] nbasis: Number of rows of input x. nbasis >= 0
!! @param[in] lmax: Maximal degree of corresponding spherical harmonics.
!!      (lmax+1)**2 = nbasis
!! @param[in] ncol: Number of columns to print
!! @param[in] icol: This number is only for printing purposes
!! @param[in] x: Actual data to print
subroutine prtsph(label, nbasis, lmax, ncol, icol, x)
    ! Inputs
    character (len=*), intent(in) :: label
    integer, intent(in) :: nbasis, lmax, ncol, icol
    real(dp), intent(in) :: x(nbasis, ncol)
    ! Local variables
    integer :: l, m, ind, noff, nprt, ic, j
    ! Print header:
    if (ncol .eq. 1) then
        write (6,'(3x,a,1x,"(column ",i4")")') label, icol
    else
        write (6,'(3x,a)') label
    endif
    ! Print entries:
    if (ncol .eq. 1) then
        do l = 0, lmax
            ind = l*l + l + 1
            do m = -l, l
                write(6,1000) l, m, x(ind+m, 1)
            end do
        end do
    else
        noff = mod(ncol, 5)
        nprt = max(ncol-noff, 0)
        do ic = 1, nprt, 5
            write(6,1010) (j, j = ic, ic+4)
            do l = 0, lmax
                ind = l*l + l + 1
                do m = -l, l
                    write(6,1020) l, m, x(ind+m, ic:ic+4)
                end do
            end do
        end do
        write (6,1010) (j, j = nprt+1, nprt+noff)
        do l = 0, lmax
            ind = l*l + l + 1
            do m = -l, l
                write(6,1020) l, m, x(ind+m, nprt+1:nprt+noff)
            end do
        end do
    end if
    1000 format(1x,i3,i4,f14.8)
    1010 format(8x,5i14)
    1020 format(1x,i3,i4,5f14.8)
end subroutine prtsph

!> Print array of quadrature points
!!
!! Prints (ngrid, ncol) array
!!
!! @param[in] label: Label to print
!! @param[in] ngrid: Number of rows of input x. ngrid >= 0
!! @param[in] ncol: Number of columns to print
!! @param[in] icol: This number is only for printing purposes
!! @param[in] x: Actual data to print
subroutine ptcart(label, ngrid, ncol, icol, x)
    ! Inputs
    character (len=*), intent(in) :: label
    integer, intent(in) :: ngrid, ncol, icol
    real(dp), intent(in) :: x(ngrid, ncol)
    ! Local variables
    integer :: ig, noff, nprt, ic, j
    ! Print header :
    if (ncol .eq. 1) then
        write (6,'(3x,a,1x,"(column ",i4")")') label, icol
    else
        write (6,'(3x,a)') label
    endif
    ! Print entries :
    if (ncol .eq. 1) then
        do ig = 1, ngrid
            write(6,1000) ig, x(ig, 1)
        enddo
    else
        noff = mod(ncol, 5)
        nprt = max(ncol-noff, 0)
        do ic = 1, nprt, 5
            write(6,1010) (j, j = ic, ic+4)
            do ig = 1, ngrid
                write(6,1020) ig, x(ig, ic:ic+4)
            end do
        end do
        write (6,1010) (j, j = nprt+1, nprt+noff)
        do ig = 1, ngrid
            write(6,1020) ig, x(ig, nprt+1:nprt+noff)
        end do
    end if
    !
    1000 format(1x,i5,f14.8)
    1010 format(6x,5i14)
    1020 format(1x,i5,5f14.8)
    !
end subroutine ptcart

!> Print dd Solution vector
!!
!! @param[in] label : Label to print
!! @param[in] vector: Vector to print

subroutine print_ddvector(ddx_data, label, vector)
    implicit none
    type(ddx_type), intent(in)  :: ddx_data
    character(len=*) :: label
    real(dp) :: vector(ddx_data % nbasis, ddx_data % nsph)
    integer :: isph, lm
  
    write(6,*) label
    do isph = 1, ddx_data % nsph
      do lm = 1, ddx_data % nbasis
        write(6,'(F15.8)') vector(lm,isph)
      end do
    end do
    return
end subroutine print_ddvector


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

subroutine fmm_constants(dmax, pm, pl, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    integer, intent(in) :: dmax, pm, pl
    real(dp), intent(out) :: vcnk((2*dmax+1)*(dmax+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1)
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

!> Compute arrays of \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$
!!
!! All values are computed recurrently from input \f$ \cos(\phi) \f$ and \f$
!! \sin(\phi) \f$ without accessing arccos or arcsin functions.
!!
!! @param[in] cphi: \f$ \cos(\phi) \f$. -1 <= `cphi` <= 1
!! @param[in] sphi: \f$ \sin(\phi) \f$. -1 <= `sphi` <= 1
!! @param[in] p: Maximal value of \f$ m \f$, for which to compute \f$ \cos(m
!!      \phi) \f$ and \f$ \sin(m\phi) \f$. `p` >= 0
!! @param[out] vcos: Array of \f$ \cos(m\phi) \f$ for \f$ m=0..p \f$. Dimension
!!      is `(p+1)`
!! @param[out] vsin: Array of \f$ \sin(m\phi) \f$ for \f$ m=0..p \f$. Dimension
!!      is `(p+1)`
subroutine trgev(cphi, sphi, p, vcos, vsin)
    ! Inputs
    real(dp), intent(in) :: cphi, sphi
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: vcos(p+1), vsin(p+1)
    ! Local variables
    integer :: m
    real(dp) :: c4phi, s4phi
    !! Treat values of p from 0 to 3 differently
    select case(p)
        ! Do nothing if p < 0
        case (:-1)
            return
        ! p = 0
        case (0)
            vcos(1) = one
            vsin(1) = zero
            return
        ! p = 1
        case (1)
            vcos(1) = one
            vsin(1) = zero
            vcos(2) = cphi
            vsin(2) = sphi
            return
        ! p = 2
        case (2)
            vcos(1) = one
            vsin(1) = zero
            vcos(2) = cphi
            vsin(2) = sphi
            vcos(3) = cphi**2 - sphi**2
            vsin(3) = 2 * cphi * sphi
            return
        ! p = 3
        case (3)
            vcos(1) = one
            vsin(1) = zero
            vcos(2) = cphi
            vsin(2) = sphi
            vcos(3) = cphi**2 - sphi**2
            vsin(3) = 2 * cphi * sphi
            vcos(4) = vcos(3)*cphi - vsin(3)*sphi
            vsin(4) = vcos(3)*sphi + vsin(3)*cphi
            return
        ! p >= 4
        case default
            vcos(1) = one
            vsin(1) = zero
            vcos(2) = cphi
            vsin(2) = sphi
            vcos(3) = cphi**2 - sphi**2
            vsin(3) = 2 * cphi * sphi
            vcos(4) = vcos(3)*cphi - vsin(3)*sphi
            vsin(4) = vcos(3)*sphi + vsin(3)*cphi
            vcos(5) = vcos(3)**2 - vsin(3)**2
            vsin(5) = 2 * vcos(3) * vsin(3)
            c4phi = vcos(5)
            s4phi = vsin(5)
    end select
    ! Define cos(m*phi) and sin(m*phi) recurrently 4 values at a time
    do m = 6, p-2, 4
        vcos(m:m+3) = vcos(m-4:m-1)*c4phi - vsin(m-4:m-1)*s4phi
        vsin(m:m+3) = vcos(m-4:m-1)*s4phi + vsin(m-4:m-1)*c4phi
    end do
    ! Work with leftover
    select case(m-p)
        case (-1)
            vcos(p-1) = vcos(p-2)*cphi - vsin(p-2)*sphi
            vsin(p-1) = vcos(p-2)*sphi + vsin(p-2)*cphi
            vcos(p) = vcos(p-2)*vcos(3) - vsin(p-2)*vsin(3)
            vsin(p) = vcos(p-2)*vsin(3) + vsin(p-2)*vcos(3)
            vcos(p+1) = vcos(p-2)*vcos(4) - vsin(p-2)*vsin(4)
            vsin(p+1) = vcos(p-2)*vsin(4) + vsin(p-2)*vcos(4)
        case (0)
            vcos(p) = vcos(p-1)*cphi - vsin(p-1)*sphi
            vsin(p) = vcos(p-1)*sphi + vsin(p-1)*cphi
            vcos(p+1) = vcos(p-1)*vcos(3) - vsin(p-1)*vsin(3)
            vsin(p+1) = vcos(p-1)*vsin(3) + vsin(p-1)*vcos(3)
        case (1)
            vcos(p+1) = vcos(p)*cphi - vsin(p)*sphi
            vsin(p+1) = vcos(p)*sphi + vsin(p)*cphi
    end select
end subroutine trgev

!> Compute associated Legendre polynomials
!!
!! Only polynomials \f$ P_\ell^m (\cos \theta) \f$ with non-negative parameter
!! \f$ m \f$ are computed. Implemented via following recurrent formulas:
!! \f{align}{
!!      &P_0^0(\cos \theta) = 1\\
!!      &P_{m+1}^{m+1}(\cos \theta) = -(2m+1) \sin \theta P_m^m(\cos \theta) \\
!!      &P_{m+1}^m(\cos \theta) = \cos \theta (2m+1) P_m^m(\cos \theta) \\
!!      &P_\ell^m(\cos \theta) = \frac{1}{\ell-m} \left( \cos \theta (2\ell-1)
!!      P_{\ell-1}^m(\cos \theta) - (\ell+m-1)P_{\ell-2}^m(\cos \theta)
!!      \right), \quad \forall \ell \geq m+2.
!! \f}
!!
!! @param[in] ctheta: \f$ \cos(\theta) \f$. -1 <= `ctheta` <= 1
!! @param[in] stheta: \f$ \sin(\theta) \f$. 0 <= `stheta` <= 1
!! @param[in] p: Maximal degree of polynomials to compute. `p` >= 0
!! @param[out] vplm: Values of associated Legendre polynomials. Dimension is
!!      `(p+1)**2`
subroutine polleg(ctheta, stheta, p, vplm)
    ! Inputs
    real(dp), intent(in) :: ctheta, stheta
    integer, intent(in) :: p
    ! Outputs
    real(dp), intent(out) :: vplm((p+1)**2)
    ! Temporary workspace
    real(dp) :: work(p+1)
    ! Call corresponding work routine
    call polleg_work(ctheta, stheta, p, vplm, work)
end subroutine polleg

!> Compute associated Legendre polynomials
!!
!! Only polynomials \f$ P_\ell^m (\cos \theta) \f$ with non-negative parameter
!! \f$ m \f$ are computed. Implemented via following recurrent formulas:
!! \f{align}{
!!      &P_0^0(\cos \theta) = 1\\
!!      &P_{m+1}^{m+1}(\cos \theta) = -(2m+1) \sin \theta P_m^m(\cos \theta) \\
!!      &P_{m+1}^m(\cos \theta) = \cos \theta (2m+1) P_m^m(\cos \theta) \\
!!      &P_\ell^m(\cos \theta) = \frac{1}{\ell-m} \left( \cos \theta (2\ell-1)
!!      P_{\ell-1}^m(\cos \theta) - (\ell+m-1)P_{\ell-2}^m(\cos \theta)
!!      \right), \quad \forall \ell \geq m+2.
!! \f}
!!
!! @param[in] ctheta: \f$ \cos(\theta) \f$. -1 <= `ctheta` <= 1
!! @param[in] stheta: \f$ \sin(\theta) \f$. 0 <= `stheta` <= 1
!! @param[in] p: Maximal degree of polynomials to compute. `p` >= 0
!! @param[out] vplm: Values of associated Legendre polynomials. Dimension is
!!      `(p+1)**2`
!! @param[out] work: Temporary workspace of a size (p+1)
subroutine polleg_work_old(ctheta, stheta, p, vplm, work)
    ! Inputs
    real(dp), intent(in) :: ctheta, stheta
    integer, intent(in) :: p
    ! Outputs
    real(dp), intent(out) :: vplm((p+1)**2)
    ! Temporary workspace
    real(dp), intent(out) :: work(p+1)
    ! Local variables
    integer :: m, ind, l, ind2, vplm_ind
    real(dp) :: fact, pmm, pmm1, pmmo, pll, fm, fl
    ! Init aux factors
    fact = one
    pmm = one
    work(1) = ctheta
    fm = two * ctheta
    do m = 1, p
        work(m+1) = work(m) + fm
    end do
    ! This loop goes over non-negative upper index of P_l^m, namely m, and
    ! defines value for all l from m to p. Polynomials P_l^m are defined only
    ! for l >= |m|. Only positive values of m are filled. Here we
    ! define at first P_m^m, then P_{m+1}^m and then all remaining P_l^m for
    ! l from m+2 to p.
    do m = 0, p
        ! index of P_m^m
        ind = (m+1) * (m+1)
        ! Store P_m^m
        vplm(ind) = pmm
        if (m .eq. p) then
            return
        end if
        fm = dble(m)
        ! index of P_{m+1}^m
        ind2 = ind + 2*m + 2
        ! P_{m+1}^m
        pmm1 = work(m+1) * pmm
        vplm(ind2) = pmm1
        ! Save value P_m^m for recursion
        pmmo = pmm
        ! Fill values of P_l^m
        do l = m+2, p
            ind2 = ind2 + 2*l
            ! pmm corresponds to P_{l-2}^m
            ! pmm1 corresponds to P_{l-1}^m
            fl = dble(l)
            ! Compute P_l^m
            pll = work(l)*pmm1 - (fl+fm-1)*pmm
            pll = pll / (fl-fm)
            ! Store P_l^m
            !vplm(l*l + l + m + 1) = pll
            vplm(ind2) = pll
            ! Save P_{l-1}^m and P_l^m for recursion
            pmm = pmm1
            pmm1 = pll
        end do
        ! Value of P_{m+1}^{m+1}
        pmm = -pmmo * fact * stheta
        fact = fact + two
    end do
end subroutine polleg_work_old

subroutine polleg_work2(ctheta, stheta, p, vplm, work)
    ! Inputs
    real(dp), intent(in) :: ctheta, stheta
    integer, intent(in) :: p
    ! Outputs
    real(dp), intent(out) :: vplm((p+1)**2)
    ! Temporary workspace
    real(dp), intent(out) :: work(p+1)
    ! Local variables
    integer :: indlm1, l, indl, m, ind
    real(dp) :: tmp, plm2, plm1, plm, ctantheta, fml
    ! Easy cases
    select case (p)
        case (0)
            vplm(1) = one
            return
        case (1)
            vplm(1) = one
            vplm(3) = ctheta
            vplm(4) = -stheta
            return
    end select
    ! Now p >= 2
    vplm(1) = one
    vplm(3) = ctheta
    vplm(4) = -stheta
    ctantheta = ctheta / stheta
    ! Prepare scalars (2*m+2)*ctantheta
    !select case (p)
    !    case (2)
    !        work(1) = two * ctantheta
    work(1) = two * ctantheta
    select case (p)
        case (2)
        case (3)
            work(2) = two * work(1)
        case (4)
            work(2) = two * work(1)
            work(3) = three * work(1)
        case (5)
            work(2) = two * work(1)
            work(3) = three * work(1)
            work(4) = four * work(1)
        case (6:)
            work(2) = two * work(1)
            work(3) = three * work(1)
            work(4) = four * work(1)
            do m = 4, p-5, 4
                work(m+1:m+4) = work(m-3:m) + work(4)
            end do
            select case (p-m)
                case (2)
                    work(p-1) = work(p-2) + work(1)
                case (3)
                    work(p-2:p-1) = work(p-4:p-3) + work(2)
                case (4)
                    work(p-3:p-1) = work(p-6:p-4) + work(3)
            end select
    end select
    ! Index of P_{l-1}^{l-1} for the next loop
    indlm1 = 4
    do l = 2, p
        ! index of P_l^l is (l+1)**2
        indl = indlm1 + 2*l + 1
        tmp = dble(2*l-1) * vplm(indlm1)
        ! P_l^l
        plm2 = - tmp * stheta
        ! P_l^{l-1}
        plm1 = tmp * ctheta
        vplm(indl) = plm2
        vplm(indl-1) = plm1
        !ind = indl-2
        fml = - two * dble(l)
        do m = l-2, 0, -1
            fml = fml - two*dble(m+1)
            ! P_l^m
            tmp = work(m+1)*plm1 + plm2
            plm = tmp / fml
            plm2 = plm1
            plm1 = plm
            vplm(indl+m-l) = plm
        end do
        indlm1 = indl
    end do
end subroutine polleg_work2

subroutine polleg_work(ctheta, stheta, p, vplm, work)
    ! Inputs
    real(dp), intent(in) :: ctheta, stheta
    integer, intent(in) :: p
    ! Outputs
    real(dp), intent(out) :: vplm((p+1)**2)
    ! Temporary workspace
    real(dp), intent(out) :: work(p+1)
    ! Local variables
    integer :: m, l, indlm1, indlm2, indl
    real(dp) :: tmp, tmp1, tmp2, tmp3, tmp4, tmp5, fl
    ! Easy cases
    select case (p)
        case (0)
            vplm(1) = one
            return
        case (1)
            vplm(1) = one
            vplm(3) = ctheta
            vplm(4) = -stheta
            return
        case (2)
            vplm(1) = one
            vplm(3) = ctheta
            vplm(4) = -stheta
            tmp1 = three * stheta
            tmp2 = ctheta * ctheta
            vplm(7) = 1.5d0*tmp2 - pt5
            vplm(8) = -tmp1 * ctheta
            vplm(9) = tmp1 * stheta
            return
        case (3)
            vplm(1) = one
            vplm(3) = ctheta
            vplm(4) = -stheta
            tmp1 = three * stheta
            tmp2 = ctheta * ctheta
            vplm(7) = 1.5d0*tmp2 - pt5
            vplm(8) = -tmp1 * ctheta
            vplm(9) = tmp1 * stheta
            tmp3 = 2.5d0 * tmp2 - pt5
            vplm(13) = tmp3*ctheta - ctheta
            vplm(14) = -tmp1 * tmp3
            tmp4 = -5d0 * stheta
            vplm(15) = tmp4 * vplm(8)
            vplm(16) = tmp4 * vplm(9)
            return
    end select
    ! Now p >= 4
    ! At first take into account l=0..3
    vplm(1) = one
    vplm(3) = ctheta
    vplm(4) = -stheta
    tmp1 = three * stheta
    tmp2 = ctheta * ctheta
    vplm(7) = 1.5d0*tmp2 - pt5
    vplm(8) = -tmp1 * ctheta
    vplm(9) = tmp1 * stheta
    tmp3 = 2.5d0 * tmp2 - pt5
    vplm(13) = tmp3*ctheta - ctheta
    vplm(14) = -tmp1 * tmp3
    tmp4 = -5d0 * stheta
    vplm(15) = tmp4 * vplm(8)
    vplm(16) = tmp4 * vplm(9)
    ! Set indexes of P_3^0 and P_2^0
    indl = 13
    indlm1 = 7
    ! Precompute constants
    tmp = -stheta * pt5
    do m = 1, p-2
        work(m) = tmp / dble(m)
    end do
    ! Loop over l >= 4
    do l = 4, p
        ! Index of P_{l-2}^0
        indlm2 = indlm1
        ! Index of P_{l-1}^0
        indlm1 = indl
        ! Index of P_l^0
        indl = indl + 2*l
        fl = dble(l)
        tmp1 = two*fl - one
        tmp2 = fl-1
        tmp3 = tmp1 * ctheta
        ! Set P_l^0
        vplm(indl) = tmp3*vplm(indlm1) - tmp2*vplm(indlm2)
        vplm(indl) = vplm(indl) / fl
        tmp4 = two * tmp2
        tmp5 = fl*fl - fl
        ! Set P_l^m for m=1..l-2
        do m = 1, l-2
            tmp4 = tmp4 + two
            tmp5 = tmp5 + tmp4
            vplm(indl+m) = work(m) * (vplm(indlm1+m+1) + &
                & tmp5*vplm(indlm1+m-1))
        end do
        ! m = l-1
        vplm(indl+l-1) = tmp3 * vplm(indlm1+l-1)
        ! m = l
        vplm(indl+l) = -stheta * tmp1 * vplm(indlm1+l-1)
    end do
end subroutine polleg_work

!> Convert input cartesian coordinate into spherical coordinate
!!
!! Output coordinate \f$ (\rho, \theta, \phi) \f$ is presented by \f$ (\rho,
!! \cos \theta, \sin \theta, \cos \phi, \sin\phi) \f$.
!!
!! @param[in] x: Cartesian coordinate
!! @param[out] rho: \f$ \rho \f$
!! @param[out] ctheta: \f$ \cos \theta \f$
!! @param[out] stheta: \f$ \sin \theta \f$
!! @param[out] cphi: \f$ \cos \phi \f$
!! @param[out] sphi: \f$ \sin \phi \f$
subroutine carttosph(x, rho, ctheta, stheta, cphi, sphi)
    ! Input
    real(dp), intent(in) :: x(3)
    ! Output
    real(dp), intent(out) :: rho, ctheta, stheta, cphi, sphi
    ! Local variables
    real(dp) :: max12, ssq12
    ! Check x(1:2) = 0
    if ((x(1) .eq. zero) .and. (x(2) .eq. zero)) then
        rho = abs(x(3))
        ctheta = sign(one, x(3))
        stheta = zero
        cphi = one
        sphi = zero
        return
    end if
    ! In other situations use sum-of-scaled-squares technique
    ! Get norm of x(1:2) and cphi with sphi outputs
    if (abs(x(2)) .gt. abs(x(1))) then
        max12 = abs(x(2))
        ssq12 = one + (x(1)/x(2))**2
    else
        max12 = abs(x(1))
        ssq12 = one + (x(2)/x(1))**2
    end if
    stheta = max12 * sqrt(ssq12)
    cphi = x(1) / stheta
    sphi = x(2) / stheta
    ! Then compute rho, ctheta and stheta outputs
    if (abs(x(3)) .gt. max12) then
        rho = one + ssq12*(max12/x(3))**2
        rho = abs(x(3)) * sqrt(rho)
        stheta = stheta / rho
        ctheta = x(3) / rho
    else
        rho = ssq12 + (x(3)/max12)**2
        rho = max12 * sqrt(rho)
        stheta = stheta / rho
        ctheta = x(3) / rho
    end if
end subroutine carttosph

!> Compute all spherical harmonics up to a given degree at a given point
!!
!! Spherical harmonics are computed for a point \f$ x / \|x\| \f$. Cartesian
!! coordinate of input `x` is translated into a spherical coordinate \f$ (\rho,
!! \theta, \phi) \f$ that is represented by \f$ \rho, \cos \theta, \sin \theta,
!! \cos \phi \f$ and \f$ \sin \phi \f$. If \f$ \rho=0 \f$ nothing is computed,
!! only zero \f$ \rho \f$ is returned without doing anything else. If \f$
!! \rho>0 \f$ values \f$ \cos \theta \f$ and \f$ \sin \theta \f$ are computed.
!! If \f$ \sin \theta \ne 0 \f$ then \f$ \cos \phi \f$ and \f$ \sin \phi \f$
!! are computed.
!! Auxiliary values of associated Legendre polynomials \f$ P_\ell^m(\theta) \f$
!! are computed along with \f$ \cos (m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!! @param[in] x: Target point
!! @param[out] rho: Euclidian length of `x`
!! @param[out] ctheta: \f$ -1 \leq \cos \theta \leq 1\f$
!! @param[out] stheta: \f$ 0 \leq \sin \theta \leq 1\f$
!! @param[out] cphi: \f$ -1 \leq \cos \phi \leq 1\f$
!! @param[out] sphi: \f$ -1 \leq \sin \phi \leq 1\f$
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[in] vscales: Scaling factors of real normalized spherical harmonics.
!!      Dimension is `(p+1)**2`
!! @param[out] vylm: Values of spherical harmonics \f$ Y_\ell^m(x) \f$.
!!      Dimension is `(p+1)**2`
!! @param[out] vplm: Values of associated Legendre polynomials \f$ P_\ell^m(
!!      \theta) \f$. Dimension is `(p+1)**2`
!! @param[out] vcos: Array of alues of \f$ \cos(m\phi) \f$ of a dimension
!!      `(p+1)`
!! @param[out] vsin: array of values of \f$ \sin(m\phi) \f$ of a dimension
!!      `(p+1)`
subroutine ylmbas(x, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Inputs
    real(dp), intent(in) :: x(3)
    integer, intent(in) :: p
    real(dp), intent(in) :: vscales((p+1)**2)
    ! Outputs
    real(dp), intent(out) :: rho, ctheta, stheta, cphi, sphi
    real(dp), intent(out) :: vylm((p+1)**2), vplm((p+1)**2)
    real(dp), intent(out) :: vcos(p+1), vsin(p+1)
    ! Local variables
    integer :: l, m, ind
    real(dp) :: max12, ssq12, tmp
    ! Get rho cos(theta), sin(theta), cos(phi) and sin(phi) from the cartesian
    ! coordinates of x. To support full range of inputs we do it via a scale
    ! and a sum of squares technique.
    ! At first we compute x(1)**2 + x(2)**2
    if (x(1) .eq. zero) then
        max12 = abs(x(2))
        ssq12 = one
    else if (abs(x(2)) .gt. abs(x(1))) then
        max12 = abs(x(2))
        ssq12 = one + (x(1)/x(2))**2
    else
        max12 = abs(x(1))
        ssq12 = one + (x(2)/x(1))**2
    end if
    ! Then we compute rho
    if (x(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(x(3)) .gt. max12) then
        rho = one + ssq12 *(max12/x(3))**2
        rho = abs(x(3)) * sqrt(rho)
    else
        rho = ssq12 + (x(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case x=0 just exit without setting any other variable
    if (rho .eq. zero) then
        return
    end if
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = x(1) / stheta
        sphi = x(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = x(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                vylm(1) = vscales(1)
                return
            case (1)
                vylm(1) = vscales(1)
                vylm(2) = - vscales(4) * stheta * sphi
                vylm(3) = vscales(3) * ctheta
                vylm(4) = - vscales(4) * stheta * cphi
                return
        end select
        ! Evaluate associated Legendre polynomials, size of temporary workspace
        ! here needs only (p+1) elements so we use vcos as a temporary
        call polleg_work(ctheta, stheta, p, vplm, vcos)
        ! Evaluate cos(m*phi) and sin(m*phi) arrays
        call trgev(cphi, sphi, p, vcos, vsin)
        ! Construct spherical harmonics
        ! l = 0
        vylm(1) = vscales(1)
        ! l = 1
        vylm(2) = - vscales(4) * stheta * sphi
        vylm(3) = vscales(3) * ctheta
        vylm(4) = - vscales(4) * stheta * cphi
        ind = 3
        do l = 2, p
            ! Offset of a Y_l^0 harmonic in vplm and vylm arrays
            ind = ind + 2*l
            !l**2 + l + 1
            ! m = 0 implicitly uses `vcos(1) = 1`
            vylm(ind) = vscales(ind) * vplm(ind)
            do m = 1, l
                ! only P_l^m for non-negative m is used/defined
                tmp = vplm(ind+m) * vscales(ind+m)
                ! m > 0
                vylm(ind+m) = tmp * vcos(m+1)
                ! m < 0
                vylm(ind-m) = tmp * vsin(m+1)
            end do
        end do
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! Set spherical coordinates
        cphi = one
        sphi = zero
        ctheta = sign(one, x(3))
        stheta = zero
        ! Set output arrays vcos and vsin
        vcos = one
        vsin = zero
        ! Evaluate spherical harmonics. P_l^m = 0 for m > 0. In the case m = 0
        ! it depends if l is odd or even. Additionally, vcos = one and vsin =
        ! zero for all elements
        vylm = zero
        vplm = zero
        do l = 0, p, 2
            ind = l**2 + l + 1
            ! only case m = 0
            vplm(ind) = one
            vylm(ind) = vscales(ind)
        end do
        do l = 1, p, 2
            ind = l**2 + l + 1
            ! only case m = 0
            vplm(ind) = ctheta
            vylm(ind) = ctheta * vscales(ind)
        end do
    end if
end subroutine ylmbas

!> Compute all spherical harmonics up to a given degree at a given point
!!
!! Attempt to improve previous version.
!! Spherical harmonics are computed for a point \f$ x / \|x\| \f$. Cartesian
!! coordinate of input `x` is translated into a spherical coordinate \f$ (\rho,
!! \theta, \phi) \f$ that is represented by \f$ \rho, \cos \theta, \sin \theta,
!! \cos \phi \f$ and \f$ \sin \phi \f$. If \f$ \rho=0 \f$ nothing is computed,
!! only zero \f$ \rho \f$ is returned without doing anything else. If \f$
!! \rho>0 \f$ values \f$ \cos \theta \f$ and \f$ \sin \theta \f$ are computed.
!! If \f$ \sin \theta \ne 0 \f$ then \f$ \cos \phi \f$ and \f$ \sin \phi \f$
!! are computed.
!! Auxiliary values of associated Legendre polynomials \f$ P_\ell^m(\theta) \f$
!! are computed along with \f$ \cos (m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!! @param[in] x: Target point
!! @param[out] rho: Euclidian length of `x`
!! @param[out] ctheta: \f$ -1 \leq \cos \theta \leq 1\f$
!! @param[out] stheta: \f$ 0 \leq \sin \theta \leq 1\f$
!! @param[out] cphi: \f$ -1 \leq \cos \phi \leq 1\f$
!! @param[out] sphi: \f$ -1 \leq \sin \phi \leq 1\f$
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[in] vscales: Scaling factors of real normalized spherical harmonics.
!!      Dimension is `(p+1)**2`
!! @param[out] vylm: Values of spherical harmonics \f$ Y_\ell^m(x) \f$.
!!      Dimension is `(p+1)**2`
!! @param[out] vplm: Values of associated Legendre polynomials \f$ P_\ell^m(
!!      \theta) \f$. Dimension is `(p+1)**2`
!! @param[out] vcos: Array of alues of \f$ \cos(m\phi) \f$ of a dimension
!!      `(p+1)`
!! @param[out] vsin: array of values of \f$ \sin(m\phi) \f$ of a dimension
!!      `(p+1)`
subroutine ylmbas2(x, sphcoo, p, vscales, vylm, vplm, vcos, vsin)
    ! Inputs
    real(dp), intent(in) :: x(3)
    integer, intent(in) :: p
    real(dp), intent(in) :: vscales((p+1)**2)
    ! Outputs
    real(dp), intent(out) :: sphcoo(5)
    real(dp), intent(out) :: vylm((p+1)**2), vplm((p+1)**2)
    real(dp), intent(out) :: vcos(p+1), vsin(p+1)
    ! Local variables
    integer :: l, m, ind
    real(dp) :: max12, ssq12, tmp, rho, ctheta, stheta, cphi, sphi
    ! Get rho cos(theta), sin(theta), cos(phi) and sin(phi) from the cartesian
    ! coordinates of x. To support full range of inputs we do it via a scale
    ! and a sum of squares technique.
    ! At first we compute x(1)**2 + x(2)**2
    if (x(1) .eq. zero) then
        max12 = abs(x(2))
        ssq12 = one
    else if (abs(x(2)) .gt. abs(x(1))) then
        max12 = abs(x(2))
        ssq12 = one + (x(1)/x(2))**2
    else
        max12 = abs(x(1))
        ssq12 = one + (x(2)/x(1))**2
    end if
    ! Then we compute rho
    if (x(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(x(3)) .gt. max12) then
        rho = one + ssq12*(max12/x(3))**2
        rho = abs(x(3)) * sqrt(rho)
    else
        rho = ssq12 + (x(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case x=0 just exit without setting any other variable
    if (rho .eq. zero) then
        sphcoo = zero
        return
    end if
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Evaluate cos(m*phi) and sin(m*phi) arrays
        cphi = x(1) / stheta
        sphi = x(2) / stheta
        call trgev(cphi, sphi, p, vcos, vsin)
        ! Normalize ctheta and stheta
        ctheta = x(3) / rho
        stheta = stheta / rho
        ! Evaluate associated Legendre polynomials
        call polleg(ctheta, stheta, p, vplm)
        ! Construct spherical harmonics
        do l = 0, p
            ! Offset of a Y_l^0 harmonic in vplm and vylm arrays
            ind = l**2 + l + 1
            ! m = 0 implicitly uses `vcos(1) = 1`
            vylm(ind) = vscales(ind) * vplm(ind)
            do m = 1, l
                ! only P_l^m for non-negative m is used/defined
                tmp = vplm(ind+m) * vscales(ind+m)
                ! m > 0
                vylm(ind+m) = tmp * vcos(m+1)
                ! m < 0
                vylm(ind-m) = tmp * vsin(m+1)
            end do
        end do
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! Set spherical coordinates
        cphi = one
        sphi = zero
        ctheta = sign(one, x(3))
        stheta = zero
        ! Set output arrays vcos and vsin
        vcos = one
        vsin = zero
        ! Evaluate spherical harmonics. P_l^m = 0 for m > 0. In the case m = 0
        ! it depends if l is odd or even. Additionally, vcos = one and vsin =
        ! zero for all elements
        vylm = zero
        vplm = zero
        do l = 0, p, 2
            ind = l**2 + l + 1
            ! only case m = 0
            vplm(ind) = one
            vylm(ind) = vscales(ind)
        end do
        do l = 1, p, 2
            ind = l**2 + l + 1
            ! only case m = 0
            vplm(ind) = ctheta
            vylm(ind) = ctheta * vscales(ind)
        end do
    end if
    ! Set output spherical coordinates
    sphcoo(1) = rho
    sphcoo(2) = ctheta
    sphcoo(3) = stheta
    sphcoo(4) = cphi
    sphcoo(5) = sphi
end subroutine ylmbas2

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

!> Integrate against spherical harmonics
!!
!! Integrate by Lebedev spherical quadrature. This function can be simply
!! substituted by a matrix-vector product.
!!
!! TODO: use dgemv. Loop of this cycle can be easily substituted by a dgemm.
!!
!! @param[in] ngrid: Number of Lebedev grid points. `ngrid` > 0
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[in] vwgrid: Values of spherical harmonics at Lebedev grid points,
!!      multiplied by weights of grid points. Dimension is ((p+1)**2, ngrid)
!! @param[in] isph: Index of given sphere. Used for debug purpose.
!! @param[in] x: Input values at grid points of the sphere. Dimension is
!!      (ngrid)
!! @param[out] xlm: Output spherical harmonics. Dimension is ((p+1)**2)
subroutine intrhs(iprint, ngrid, p, vwgrid, ldvwgrid, isph, x, xlm)
    ! Inputs
    integer, intent(in) :: iprint, ngrid, p, ldvwgrid, isph
    real(dp), intent(in) :: vwgrid(ldvwgrid, ngrid)
    real(dp), intent(in) :: x(ngrid)
    ! Output
    real(dp), intent(out) :: xlm((p+1)**2)
    ! Local variables
    integer :: igrid
    ! Initialize
    xlm = zero
    ! Accumulate over integration points
    do igrid = 1, ngrid
        xlm = xlm + vwgrid(:(p+1)**2,igrid)*x(igrid)
    end do
    ! Printing (these functions are not implemented yet)
    if (iprint .ge. 5) then
        call ptcart('pot', ngrid, 1, isph, x)
        call prtsph('vlm', (p+1)**2, p, 1, isph, xlm)
    end if
end subroutine intrhs

!> Compute first derivatives of spherical harmonics
!!
!! @param[in] x:
!! @param[out] basloc:
!! @param[out] dbsloc:
!! @param[out] vplm:
!! @param[out] vcos:
!! @param[out] vsin:
!!
!!
!! TODO: rewrite code and fill description. Computing sqrt(one-cthe*cthe)
!! reduces effective range of input double precision values. cthe*cthe for
!! cthe=1d+155 is NaN.
subroutine dbasis(ddx_data, x, basloc, dbsloc, vplm, vcos, vsin)
    type(ddx_type) :: ddx_data
    real(dp), dimension(3),        intent(in)    :: x
    real(dp), dimension((ddx_data % lmax+1)**2),   intent(inout) :: basloc, vplm
    real(dp), dimension(3,(ddx_data % lmax +1)**2), intent(inout) :: dbsloc
    real(dp), dimension(ddx_data % lmax+1),   intent(inout) :: vcos, vsin
    integer :: l, m, ind
    real(dp)  :: cthe, sthe, cphi, sphi, plm, fln, pp1, pm1, pp, VC, VS
    real(dp)  :: et(3), ep(3)
    !     get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
    !     coordinates of x.
    cthe = x(3)
    sthe = sqrt(one - cthe*cthe)
    !     not ( NORTH or SOUTH pole )
    if ( sthe.ne.zero ) then
        cphi = x(1)/sthe
        sphi = x(2)/sthe
        !     NORTH or SOUTH pole
    else
        cphi = one
        sphi = zero
    end if
    !     evaluate the derivatives of theta and phi:
    et(1) = cthe*cphi
    et(2) = cthe*sphi
    et(3) = -sthe
    !     not ( NORTH or SOUTH pole )
    if( sthe.ne.zero ) then
        ep(1) = -sphi/sthe
        ep(2) = cphi/sthe
        ep(3) = zero
        !     NORTH or SOUTH pole
    else
        ep(1) = zero
        ep(2) = one
        ep(3) = zero
    end if
    VC=zero
    VS=cthe
    !     evaluate the generalized legendre polynomials. Temporary workspace
    !       is of size (p+1) here, so we use vcos for that purpose
    call polleg_work( cthe, sthe, ddx_data % lmax, vplm, vcos )
    !
    !     evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is 
    !     pointless if z = 1, as the only non vanishing terms will be the 
    !     ones with m=0.
    !
    !     not ( NORTH or SOUTH pole )
    if ( sthe.ne.zero ) then
        call trgev( cphi, sphi, ddx_data % lmax, vcos, vsin )
        !     NORTH or SOUTH pole
    else
        vcos = one
        vsin = zero
    end if
    !     now build the spherical harmonics. we will distinguish m=0,
    !     m>0 and m<0:
    !
    basloc = zero
    dbsloc = zero
    do l = 0, ddx_data % lmax
        ind = l*l + l + 1
        ! m = 0
        fln = ddx_data % vscales(ind)   
        basloc(ind) = fln*vplm(ind)
        if (l.gt.0) then
            dbsloc(:,ind) = fln*vplm(ind+1)*et(:)
        else
            dbsloc(:,ind) = zero
        end if
        !dir$ simd
        do m = 1, l
            fln = ddx_data % vscales(ind+m)
            plm = fln*vplm(ind+m)   
            pp1 = zero
            if (m.lt.l) pp1 = -pt5*vplm(ind+m+1)
            pm1 = pt5*(dble(l+m)*dble(l-m+1)*vplm(ind+m-1))
            pp  = pp1 + pm1   
            !
            !         m > 0
            !         -----
            !
            basloc(ind+m) = plm*vcos(m+1)
            !          
            !         not ( NORTH or SOUTH pole )
            if ( sthe.ne.zero ) then
                ! 
                dbsloc(:,ind+m) = -fln*pp*vcos(m+1)*et(:) - dble(m)*plm*vsin(m+1)*ep(:)
                !
                !            
                !         NORTH or SOUTH pole
            else
                !                  
                dbsloc(:,ind+m) = -fln*pp*vcos(m+1)*et(:) - fln*pp*ep(:)*VC
                !
                !
            endif
            !
            !         m < 0
            !         -----
            !
            basloc(ind-m) = plm*vsin(m+1)
            !          
            !         not ( NORTH or SOUTH pole )
            if ( sthe.ne.zero ) then
                ! 
                dbsloc(:,ind-m) = -fln*pp*vsin(m+1)*et(:) + dble(m)*plm*vcos(m+1)*ep(:)
                !
                !         NORTH or SOUTH pole
            else
                !                  
                dbsloc(:,ind-m) = -fln*pp*vsin(m+1)*et(:) - fln*pp*ep(:)*VS
                !            
            endif
            !  
        enddo
    enddo
end subroutine dbasis

! Purpose : compute
!
!                               l
!             sum   4pi/(2l+1) t  * Y_l^m( s ) * sigma_l^m
!             l,m                                           
!
!           which is need to compute action of COSMO matrix L.
!------------------------------------------------------------------------------------------------
!
!!> TODO
real(dp) function intmlp( ddx_data, t, sigma, basloc )
!  
      implicit none
      type(ddx_type) :: ddx_data
      real(dp), intent(in) :: t
      real(dp), dimension((ddx_data % lmax+1)**2), intent(in) :: sigma, basloc
!
      integer :: l, ind
      real(dp)  :: tt, ss, fac
!
!------------------------------------------------------------------------------------------------
!
!     initialize t^l
      tt = one
!
!     initialize
      ss = zero
!
!     loop over l
      do l = 0, ddx_data % lmax
!      
        ind = l*l + l + 1
!
!       update factor 4pi / (2l+1) * t^l
        fac = tt / ddx_data % vscales(ind)**2
!
!       contract over l,m and accumulate
        ss = ss + fac * dot_product( basloc(ind-l:ind+l), &
                                     sigma( ind-l:ind+l)   )
!
!       update t^l
        tt = tt*t
!        
      enddo
!      
!     redirect
      intmlp = ss
!
!
end function intmlp

! Purpose : weigh potential at cavity points by characteristic function "ui"
!------------------------------------------------------------------------------------------------
!> TODO
subroutine wghpot( ddx_data, phi, phi_grid, g)
!
      implicit none
!
    type(ddx_type) :: ddx_data
      real(dp), dimension(ddx_data % ncav),       intent(in)  :: phi
      real(dp), dimension(ddx_data % ngrid, ddx_data % nsph), intent(out) :: g
      real(dp), dimension(ddx_data % ngrid, ddx_data % nsph), intent(out) :: phi_grid
!
    integer isph, ig, ic
!
!------------------------------------------------------------------------------------------------
!
!   initialize
    ic = 0 ; g(:,:)=0.d0
    phi_grid = zero
!      
!   loop over spheres
    do isph = 1, ddx_data % nsph
!
!   loop over points
      do ig = 1, ddx_data % ngrid
!
!       nonzero contribution from point
        if ( ddx_data % ui(ig,isph).ne.zero ) then
!
!         advance cavity point counter
          ic = ic + 1
          phi_grid(ig, isph) = phi(ic)
!            
!         weigh by (negative) characteristic function
          g(ig,isph) = -ddx_data % ui(ig,isph) * phi(ic)
        endif
!          
       enddo
   enddo
end subroutine wghpot

! Purpose : compute H-norm
!------------------------------------------------------------------------------------------------
!> TODO
subroutine hsnorm( ddx_data, u, unorm )
!          
      implicit none
      type(ddx_type) :: ddx_data
      real(dp), dimension((ddx_data % lmax+1)**2), intent(in)    :: u
      real(dp),                    intent(inout) :: unorm
!
      integer :: l, m, ind
      real(dp)  :: fac
!
!------------------------------------------------------------------------------------------------
!
!     initialize
      unorm = zero
!      
!     loop over l
      do l = 0, ddx_data % lmax
!      
!       first index associated to l
        ind = l*l + l + 1
!
!       scaling factor
        fac = one/(one + dble(l))
!
!       loop over m
        do m = -l, l
!
!         accumulate
          unorm = unorm + fac*u(ind+m)*u(ind+m)
!          
        enddo
      enddo
!
!     the much neglected square root
      unorm = sqrt(unorm)
!
      return
!
!
end subroutine hsnorm

! compute the h^-1/2 norm of the increment on each sphere, then take the
! rms value.
!-------------------------------------------------------------------------------
!
!> TODO
real(dp) function hnorm(ddx_data, x)
    type(ddx_type), intent(in) :: ddx_data
      real(dp),  dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: x
!
      integer                                     :: isph, istatus
      real(dp)                                      :: vrms, vmax
      real(dp), allocatable                         :: u(:)
!
!-------------------------------------------------------------------------------
!
!     allocate workspace
      allocate( u(ddx_data % nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'hnorm: allocation failed !'
        stop
      endif
!
!     loop over spheres
      do isph = 1, ddx_data % nsph
!
!       compute norm contribution
        call hsnorm(ddx_data, x(:,isph), u(isph))
      enddo
!
!     compute rms of norms
      call rmsvec( ddx_data % nsph, u, vrms, vmax )
!
!     return value
      hnorm = vrms
!
!     deallocate workspace
      deallocate( u , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'hnorm: deallocation failed !'
        stop
      endif
!
!
end function hnorm

!------------------------------------------------------------------------------
! Purpose : compute root-mean-square and max norm
!------------------------------------------------------------------------------
subroutine rmsvec( n, v, vrms, vmax )
!
      implicit none
      integer,               intent(in)    :: n
      real(dp),  dimension(n), intent(in)    :: v
      real(dp),                intent(inout) :: vrms, vmax
!
      integer :: i
      real(dp), parameter :: zero=0.0d0
!      
!------------------------------------------------------------------------------
!      
!     initialize
      vrms = zero
      vmax = zero
!
!     loop over entries
      do i = 1,n
!
!       max norm
        vmax = max(vmax,abs(v(i)))
!
!       rms norm
        vrms = vrms + v(i)*v(i)
!        
      enddo
!
!     the much neglected square root
      vrms = sqrt(vrms/dble(n))
!      
      return
!      
!      
endsubroutine rmsvec

!-----------------------------------------------------------------------------------
! Purpose : compute
!
!   v_l^m = v_l^m +
!
!               4 pi           l
!     sum  sum  ---- ( t_n^ji )  Y_l^m( s_n^ji ) W_n^ji [ \xi_j ]_n
!      j    n   2l+1
!
! which is related to the action of the adjont COSMO matrix L^* in the following
! way. Set
!
!   [ \xi_j ]_n = sum  Y_l^m( s_n ) [ s_j ]_l^m
!                 l,m
!
! then
!
!   v_l^m = -   sum    (L^*)_ij s_j
!             j \ne i 
!
! The auxiliary quantity [ \xi_j ]_l^m needs to be computed explicitly.
!-----------------------------------------------------------------------------------
!
!> TODO
subroutine adjrhs( ddx_data, isph, xi, vlm, basloc, vplm, vcos, vsin )
!
      implicit none
      type(ddx_type), intent(in) :: ddx_data
      integer,                       intent(in)    :: isph
      real(dp), dimension(ddx_data % ngrid, ddx_data % nsph), intent(in)    :: xi
      real(dp), dimension((ddx_data % lmax+1)**2),     intent(inout) :: vlm
      real(dp), dimension((ddx_data % lmax+1)**2),     intent(inout) :: basloc, vplm
      real(dp), dimension(ddx_data % lmax+1),     intent(inout) :: vcos, vsin
!
      integer :: ij, jsph, ig, l, ind, m
      real(dp)  :: vji(3), vvji, tji, sji(3), xji, oji, fac, ffac, t
      real(dp) :: rho, ctheta, stheta, cphi, sphi
      real(dp) :: work((ddx_data % lmax+1)*(ddx_data % lmax+1)+3*ddx_data % lmax)
!      
!-----------------------------------------------------------------------------------
!
!     loop over neighbors of i-sphere
      do ij = ddx_data % inl(isph), ddx_data % inl(isph+1)-1
!
!       j-sphere is neighbor
        jsph = ddx_data % nl(ij)
!
!       loop over integration points
        do ig = 1, ddx_data % ngrid
!        
!         compute t_n^ji = | r_j + \rho_j s_n - r_i | / \rho_i
          vji  = ddx_data % csph(:,jsph) + ddx_data % rsph(jsph)* &
              & ddx_data % cgrid(:,ig) - ddx_data % csph(:,isph)
          vvji = sqrt(dot_product(vji,vji))
          tji  = vvji/ddx_data % rsph(isph)
!
!         point is INSIDE i-sphere (+ transition layer)
!         ---------------------------------------------
          if ( tji.lt.( one + (ddx_data % se+one)/two*ddx_data % eta ) ) then
!                  
!           compute s_n^ji
            !sji = vji/vvji
!
!           compute \chi( t_n^ji )
            xji = fsw( tji, ddx_data % se, ddx_data % eta )
!
!           compute W_n^ji
            if ( ddx_data % fi(ig,jsph).gt.one ) then
!                    
              oji = xji/ ddx_data % fi(ig,jsph)
!              
            else
!                    
              oji = xji
!              
            endif
!            
!           compute Y_l^m( s_n^ji )
            !call ylmbas(sji, rho, ctheta, stheta, cphi, sphi, &
            !    & ddx_data % lmax, ddx_data % vscales, basloc, vplm, &
            !    & vcos, vsin )
!            
!           initialize ( t_n^ji )^l
            !t = one
!            
!           compute w_n * xi(n,j) * W_n^ji
            fac = ddx_data % wgrid(ig) * xi(ig,jsph) * oji

            call fmm_l2p_adj_work(vji, fac, ddx_data % rsph(isph), &
                & ddx_data % lmax, ddx_data % vscales_rel, one, vlm, work)
!            
!           loop over l
            !do l = 0, ddx_data % lmax
!            
            !  ind  = l*l + l + 1
!
!             compute 4pi / (2l+1) * ( t_n^ji )^l * w_n * xi(n,j) * W_n^ji
            !  ffac = fac*t/ ddx_data % vscales(ind)**2
!
!             loop over m
            !  do m = -l,l
!              
            !    vlm(ind+m) = vlm(ind+m) + ffac*basloc(ind+m)
!                
            !  enddo
!
!             update ( t_n^ji )^l
            !  t = t*tji
!              
            !enddo
!            
          endif
        enddo
      enddo
!
!
end subroutine adjrhs

!------------------------------------------------------------------------
! Purpose : compute
!
!   \Phi( n ) =
!     
!                       4 pi           l
!     sum  W_n^ij  sum  ---- ( t_n^ij )  Y_l^m( s_n^ij ) [ \sigma_j ]_l^m
!      j           l,m  2l+1
!
! which is related to the action of the COSMO matrix L in the following
! way :
!
!   -   sum    L_ij \sigma_j = sum  w_n Y_l^m( s_n ) \Phi( n ) 
!     j \ne i                   n
!
! This second step is performed by routine "intrhs".
!------------------------------------------------------------------------
!
!> TODO
subroutine calcv( ddx_data, first, isph, pot, sigma, basloc, vplm, vcos, vsin )
!
    type(ddx_type) :: ddx_data
      logical,                        intent(in)    :: first
      integer,                        intent(in)    :: isph
      real(dp), dimension((ddx_data % lmax+1)**2, ddx_data % nsph), intent(in)    :: sigma
      real(dp), dimension(ddx_data % ngrid),       intent(inout) :: pot
      real(dp), dimension((ddx_data % lmax+1)**2),      intent(inout) :: basloc
      real(dp), dimension((ddx_data % lmax+1)**2),      intent(inout) :: vplm
      real(dp), dimension(ddx_data % lmax+1),      intent(inout) :: vcos
      real(dp), dimension(ddx_data % lmax+1),      intent(inout) :: vsin
!
      integer :: its, ij, jsph
      real(dp)  :: vij(3), sij(3)
      real(dp)  :: vvij, tij, xij, oij, stslm, stslm2, stslm3, &
          & thigh, rho, ctheta, stheta, cphi, sphi
      real(dp) :: work((ddx_data % lmax+1)*(ddx_data % lmax+1)+3*ddx_data % lmax)
!
!------------------------------------------------------------------------
      thigh = one + pt5*(ddx_data % se + one)*ddx_data % eta
!
!     initialize
      pot(:) = zero
!
!     if 1st iteration of Jacobi method, then done!
      if ( first )  return
!
!     loop over grid points
      do its = 1, ddx_data % ngrid
!
!       contribution from integration point present
        if ( ddx_data % ui(its,isph).lt.one ) then
!
!         loop over neighbors of i-sphere
          do ij = ddx_data % inl(isph), ddx_data % inl(isph+1)-1
!
!           neighbor is j-sphere
            jsph = ddx_data % nl(ij)
!            
!           compute t_n^ij = | r_i + \rho_i s_n - r_j | / \rho_j
            vij  = ddx_data % csph(:,isph) + ddx_data % rsph(isph)* &
                & ddx_data % cgrid(:,its) - ddx_data % csph(:,jsph)
            vvij = sqrt( dot_product( vij, vij ) )
            tij  = vvij / ddx_data % rsph(jsph) 
!
!           point is INSIDE j-sphere
!           ------------------------
            if ( tij.lt.thigh .and. tij.gt.zero) then
!!
!!             compute s_n^ij = ( r_i + \rho_i s_n - r_j ) / | ... |
!              sij = vij / vvij
!!            
!!             compute \chi( t_n^ij )
              xij = fsw( tij, ddx_data % se, ddx_data % eta )
!!
!!             compute W_n^ij
              if ( ddx_data % fi(its,isph).gt.one ) then
!!
                oij = xij / ddx_data % fi(its,isph)
!!
              else
!!
                oij = xij
!!
              endif
!!
!!             compute Y_l^m( s_n^ij )
!              !call ylmbas( sij, basloc, vplm, vcos, vsin )
!              call ylmbas(sij, rho, ctheta, stheta, cphi, sphi, ddx_data % lmax, &
!                  & ddx_data % vscales, basloc, vplm, vcos, vsin)
!!                    
!!             accumulate over j, l, m
!              pot(its) = pot(its) + oij * intmlp( ddx_data, tij, sigma(:,jsph), basloc )
!!              
                call fmm_l2p_work(vij, ddx_data % rsph(jsph), ddx_data % lmax, &
                    & ddx_data % vscales_rel, oij, sigma(:, jsph), one, pot(its), &
                    & work)
            endif
          end do
        end if
      end do
!      
!      
!      
end subroutine calcv
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
! Purpose : compute
!
! \xi(n,i) = 
!
!  sum w_n U_n^i Y_l^m(s_n) [S_i]_l^m
!  l,m
!
!------------------------------------------------------------------------
subroutine ddmkxi( ddx_data, s, xi)
!
    type(ddx_type) :: ddx_data
       real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in)    :: s
       real(dp), dimension(ddx_data % ncav),      intent(inout) :: xi
!
       integer :: its, isph, ii
!
       ii = 0
       do isph = 1, ddx_data % nsph
         do its = 1, ddx_data % ngrid
           if (ddx_data % ui(its,isph) .gt. zero) then
             ii     = ii + 1
             xi(ii) = ddx_data % ui(its,isph)*dot_product(ddx_data % vwgrid(:,its),s(:,isph))
           end if
         end do
       end do
!
       return
end subroutine ddmkxi
!
!------------------------------------------------------------------------
! Purpose : compute
!
! \zeta(n,i) = 
!
!  1/2 f(\eps) sum w_n U_n^i Y_l^m(s_n) [S_i]_l^m
!              l,m
!
!------------------------------------------------------------------------
subroutine ddmkzeta( ddx_data, s, zeta)
!
    type(ddx_type) :: ddx_data
       real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in)    :: s
       real(dp), dimension(ddx_data % ncav),      intent(inout) :: zeta
!
       integer :: its, isph, ii
!
       ii = 0
       do isph = 1, ddx_data % nsph
         do its = 1, ddx_data % ngrid
           if (ddx_data % ui(its,isph) .gt. zero) then
             ii     = ii + 1
             zeta(ii) = ddx_data % ui(its,isph)* &
                 & dot_product(ddx_data % vwgrid(:,its),s(:,isph))
           end if
         end do
       end do
!
       zeta = pt5*((ddx_data % eps-one)/ddx_data % eps)*zeta
       return
end subroutine ddmkzeta

!> Compute preconditioner
!!
!! assemble the diagonal blocks of the reps matrix
!! then invert them to build the preconditioner
subroutine mkprec(ddx_data)
    ! Inouts
    type(ddx_type), intent(inout) :: ddx_data
    integer :: isph, lm, ind, l1, m1, ind1, its, istatus
    real(dp)  :: f, f1
    integer, allocatable :: ipiv(:)
    real(dp),  allocatable :: work(:)
    external :: dgetrf, dgetri
    ! Allocation of temporaries
    allocate(ipiv(ddx_data % nbasis),work(ddx_data % nbasis**2),stat=istatus)
    if (istatus.ne.0) then
        write(*,*) 'mkprec : allocation failed !'
        stop
    endif
    ! Init
    ddx_data % rx_prc = zero
    ! Off-diagonal part
    do isph = 1, ddx_data % nsph
        do its = 1, ddx_data % ngrid
            f = twopi* ddx_data % ui(its,isph) * ddx_data % wgrid(its)
            do l1 = 0, ddx_data % lmax
                ind1 = l1*l1 + l1 + 1
                do m1 = -l1, l1
                    f1 = f*ddx_data % vgrid(ind1 + m1,its)/(two*dble(l1) + one)
                    do lm = 1, ddx_data % nbasis
                        ddx_data % rx_prc(lm,ind1 + m1,isph) = ddx_data % rx_prc(lm,ind1 + m1,isph) + &
                            & f1*ddx_data % vgrid(lm,its)
                    end do
                end do
            end do
        end do
    end do
    ! add diagonal
    f = twopi*(ddx_data % eps + one)/(ddx_data % eps - one)
    do isph = 1, ddx_data % nsph
        do lm = 1, ddx_data % nbasis
            ddx_data % rx_prc(lm,lm,isph) = ddx_data % rx_prc(lm,lm,isph) + f
        end do
    end do
    ! invert the blocks
    do isph = 1, ddx_data % nsph
        call dgetrf(ddx_data % nbasis, ddx_data % nbasis, &
            & ddx_data % rx_prc(:,:,isph), ddx_data % nbasis, ipiv, istatus)
        if (istatus.ne.0) then 
            write(6,*) 'lu failed in mkprc'
            stop
        end if
        call dgetri(ddx_data % nbasis, ddx_data % rx_prc(:,:,isph), &
            & ddx_data % nbasis, ipiv, work, ddx_data % nbasis**2, istatus)
        if (istatus.ne.0) then 
            write(6,*) 'inversion failed in mkprc'
            stop
        end if
    end do
    ! Cleanup temporaries
    deallocate(work, ipiv, stat=istatus)
    if (istatus.ne.0) then
        write(*,*) 'mkprec : deallocation failed !'
        stop
    end if
endsubroutine mkprec

!> Accumulate a multipole expansion induced by a particle of a given charge
!!
!! Computes the following sums:
!! \f[
!!      \forall \ell=0, \ldots, p, \quad \forall m=-\ell, \ldots, \ell : \quad
!!      M_\ell^m = \beta M_\ell^m + \frac{q \|c\|^\ell}{r^{\ell+1}}
!!      Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of output harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \beta \f$ is a scaling factor, \f$ q \f$ and \f$
!! c \f$ are a charge and coordinates of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{\ell+1} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^{-\ell-1} \f$.
!!
!! @param[in] c: Radius-vector from the particle to the center of harmonics
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output multipole spherical harmonics
!! @param[in] p: Maximal degree of output multipole spherical harmonics
!! @param[in] vscales: Normalization constants for spherical harmonics
!! @param[in] beta: Scaling factor for `dst_m`
!! @param[inout] dst_m: Multipole coefficients
!!
!! @sa fmm_m2p
subroutine fmm_p2m(c, src_q, dst_r, p, vscales, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Workspace
    real(dp) :: work(2*(p+1)*(p+2))
    ! Call corresponding work routine
    call fmm_p2m_work(c, src_q, dst_r, p, vscales, beta, dst_m, work)
end subroutine fmm_p2m

!> Accumulate a multipole expansion induced by a particle of a given charge
!!
!! Computes the following sums:
!! \f[
!!      \forall \ell=0, \ldots, p, \quad \forall m=-\ell, \ldots, \ell : \quad
!!      M_\ell^m = \beta M_\ell^m + \frac{q \|c\|^\ell}{r^{\ell+1}}
!!      Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of output harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \beta \f$ is a scaling factor, \f$ q \f$ and \f$
!! c \f$ are a charge and coordinates of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{\ell+1} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^{-\ell-1} \f$.
!!
!! @param[in] c: Radius-vector from the particle to the center of harmonics
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output multipole spherical harmonics
!! @param[in] p: Maximal degree of output multipole spherical harmonics
!! @param[in] vscales: Normalization constants for spherical harmonics
!! @param[in] beta: Scaling factor for `dst_m`
!! @param[inout] dst_m: Multipole coefficients
!! @param[out] work: Temporary workspace of a size (2*(p+1)*(p+2))
!!
!! @sa fmm_m2p
subroutine fmm_p2m_work(c, src_q, dst_r, p, vscales, beta, dst_m, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Workspace
    real(dp), intent(out), target :: work(2*(p+1)*(p+2))
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, t, rcoef
    integer :: n, ind, vylm1, vylm2, vplm1, vplm2, vcos1, vcos2, vsin1, vsin2
    ! In case src_q is zero just scale output properly
    if (src_q .eq. zero) then
        ! Zero init output if beta is also zero
        if (beta .eq. zero) then
            dst_m = zero
        ! Scale output by beta otherwise
        else
            dst_m = beta * dst_m
        end if
        ! Exit subroutine
        return
    end if
    ! Now src_q is non-zero
    ! Mark first and last elements of all work subarrays
    vylm1 = 1
    vylm2 = (p+1)**2
    vplm1 = vylm2 + 1
    vplm2 = 2 * vylm2
    vcos1 = vplm2 + 1
    vcos2 = vcos1 + p
    vsin1 = vcos2 + 1
    vsin2 = vsin1 + p
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, &
        & work(vylm1:vylm2), work(vplm1:vplm2), work(vcos1:vcos2), &
        & work(vsin1:vsin2))
    ! Harmonics are available only if rho > 0
    if (rho .ne. zero) then
        rcoef = rho / dst_r
        t = src_q / dst_r
        ! Ignore input `dst_m` in case of a zero scaling factor beta
        if (beta .eq. zero) then
            do n = 0, p
                ind = n*n + n + 1
                ! Array vylm is in the beginning of work array
                dst_m(ind-n:ind+n) = t * work(ind-n:ind+n)
                t = t * rcoef
            end do
        ! Update `dst_m` otherwise
        else
            do n = 0, p
                ind = n*n + n + 1
                ! Array vylm is in the beginning of work array
                dst_m(ind-n:ind+n) = beta*dst_m(ind-n:ind+n) + &
                    & t*work(ind-n:ind+n)
                t = t * rcoef
            end do
        end if
    ! Naive case of rho = 0
    else
        ! Ignore input `dst_m` in case of a zero scaling factor beta
        if (beta .eq. zero) then
            dst_m(1) = src_q / dst_r / sqrt4pi
            dst_m(2:) = zero
        ! Update `m` otherwise
        else
            dst_m(1) = beta*dst_m(1) + src_q/dst_r/sqrt4pi
            dst_m(2:) = beta * dst_m(2:)
        end if
    end if
end subroutine fmm_p2m_work

!> Accumulate potential, induced by multipole spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `v`
!! @param[inout] v: Value of induced potential
subroutine fmm_m2p(c, r, p, vscales_rel, alpha, src_m, beta, v)
    ! Inputs
    real(dp), intent(in) :: c(3), r, vscales_rel((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: v
    ! Temporary workspace
    real(dp) :: work((p+1)*(p+1)+3*p)
    ! Call corresponding work routine
    call fmm_m2p_work(c, r, p, vscales_rel, alpha, src_m, beta, v, work)
end subroutine fmm_m2p

!> Accumulate potential, induced by multipole spherical harmonics
!!
!! This function relies on a user-provided temporary workspace
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `v`
!! @param[inout] v: Value of induced potential
!! @param[out] work: Temporary workspace of size (p+1)*(p+1)+3*p
subroutine fmm_m2p_work(c, src_r, p, vscales_rel, alpha, src_m, &
        & beta, dst_v, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, vscales_rel((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Workspace
    real(dp), intent(out), target :: work((p+1)*(p+1)+3*p)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, tmp4, tmp5, fl, max12, ssq12
    integer :: vplm1, vplm2, vcos1, vcos2, vsin1, vsin2, vminv1, vminv2, &
        & l, m, indl, indlm1, indlm2
    real(dp), pointer :: vplm(:), vcos(:), vsin(:), vminv(:)
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    if (rho .eq. zero) then
        return
    end if
    ! Compute the actual induced potential
    rcoef = src_r / rho
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                ! l = 0
                dst_v = dst_v + alpha*rcoef*vscales_rel(1)*src_m(1)
                return
            case (1)
                ! l = 0
                tmp = src_m(1) * vscales_rel(1)
                ! l = 1
                tmp2 = ctheta * src_m(3) * vscales_rel(3)
                tmp3 = vscales_rel(4) * stheta
                tmp2 = tmp2 - tmp3*sphi*src_m(2)
                tmp2 = tmp2 - tmp3*cphi*src_m(4)
                dst_v = dst_v + alpha*rcoef*(tmp+rcoef*tmp2)
                return
        end select
        ! Mark first and last elements of all work subarrays
        vplm1 = 1
        vplm2 = (p+1)**2
        vplm => work(vplm1:vplm2)
        vcos1 = vplm2 + 1
        vcos2 = vcos1 + p
        vcos => work(vcos1:vcos2)
        vsin1 = vcos2 + 1
        vsin2 = vsin1 + p
        vsin => work(vsin1:vsin2)
        vminv1 = vsin2 + 1
        vminv2 = vminv1 + p-3
        vminv => work(vminv1:vminv2)
        ! Evaluate cos(m*phi) and sin(m*phi) arrays
        call trgev(cphi, sphi, p, vcos, vsin)
        ! Prepare -stheta*pt5/m
        tmp2 = -stheta * pt5
        do m = 1, p-2
            vminv(m) = tmp2 / dble(m)
        end do
        ! Construct spherical harmonics
        ! l = 0
        vplm(1) = one
        tmp = src_m(1) * vscales_rel(1)
        ! l = 1
        vplm(3) = ctheta
        vplm(4) = -stheta
        tmp2 = ctheta * src_m(3) * vscales_rel(3)
        tmp2 = tmp2 - stheta*vscales_rel(4)*(sphi*src_m(2)+cphi*src_m(4))
        tmp = tmp + rcoef*tmp2
        t = rcoef
        indl = 3
        indlm1 = 1
        do l = 2, p
            t = t * rcoef
            ! Offset of a P_{l-2}^0 harmonic
            indlm2 = indlm1
            ! Offset of a P_{l-1}^0 harmonic
            indlm1 = indl
            ! Offset of a P_l^0 harmonic
            indl = indl + 2*l
            ! Some temp constants
            fl = dble(l)
            tmp1 = two*fl - one
            tmp2 = fl - one
            tmp3 = tmp1 * ctheta
            tmp4 = two * tmp2
            tmp5 = fl*fl - fl
            ! Y_l^0
            vplm(indl) = tmp3*vplm(indlm1) - tmp2*vplm(indlm2)
            vplm(indl) = vplm(indl) / fl
            tmp2 = vscales_rel(indl) * vplm(indl) * src_m(indl)
            ! Y_l^m for m=1..l-2
            do m = 1, l-2
                tmp4 = tmp4 + two
                tmp5 = tmp5 + tmp4
                ! P_l^m
                vplm(indl+m) = vminv(m) * (vplm(indlm1+m+1) + &
                    & tmp5*vplm(indlm1+m-1))
                tmp2 = tmp2 + vplm(indl+m)*(vcos(m+1)*src_m(indl+m)+ &
                    & vsin(m+1)*src_m(indl-m))*vscales_rel(indl+m)
            end do
            ! m = l-1
            vplm(indl+l-1) = tmp3 * vplm(indlm1+l-1)
            tmp2 = tmp2 + vplm(indl+l-1)*(vcos(l)*src_m(indl+l-1)+ &
                & vsin(l)*src_m(indl-l+1))*vscales_rel(indl+l-1)
            ! m = l
            vplm(indl+l) = -stheta * tmp1 * vplm(indlm1+l-1)
            tmp2 = tmp2 + vplm(indl+l)*(vcos(l+1)*src_m(indl+l)+ &
                & vsin(l+1)*src_m(indl-l))*vscales_rel(indl+l)
            tmp = tmp + t*tmp2
        end do
        dst_v = dst_v + alpha*rcoef*tmp
    ! Case of x(1:2) = 0 and x(3) != 0
    else
    ! TODO: Add proper check for this routine, I am leaving here stop
    ! instruction to indicate that it is not tested
        stop "fmm_m2p_work untested branch"
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^m = ctheta^m in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef.
        rcoef = sign(rcoef, c(3))
        ! Init t and proceed with accumulation of a potential
        t = alpha
        indl = 1
        do l = 0, p
            ! Index of Y_l^0
            indl = indl + 2*l
            ! Update t
            t = t * rcoef
            ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
            dst_v = dst_v + t*src_m(indl)*vscales_rel(indl)
        end do
    end if
end subroutine fmm_m2p_work

!> Adjoint M2P operation
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output multipole spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole coefficients. Dimension is `(p+1)**2`
subroutine fmm_m2p_adj(c, src_q, dst_r, p, vscales_rel, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales_rel((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Workspace
    real(dp) :: work((p+1)*(p+1)+3*p)
    ! Call corresponding work routine
    call fmm_m2p_adj_work(c, src_q, dst_r, p, vscales_rel, beta, dst_m, work)
end subroutine fmm_m2p_adj

!> Adjoint M2P operation
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output multipole spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[out] work: Temporary workspace of a size ((p+1)*(p+1)+3*p)
subroutine fmm_m2p_adj_work(c, src_q, dst_r, p, vscales_rel, beta, dst_m, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales_rel((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Workspace
    real(dp), intent(out), target :: work((p+1)*(p+1)+3*p)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, tmp4, tmp5, fl, max12, ssq12
    integer :: vplm1, vplm2, vcos1, vcos2, vsin1, vsin2, vminv1, vminv2, &
        & l, m, indl, indlm1, indlm2
    real(dp), pointer :: vplm(:), vcos(:), vsin(:), vminv(:)
    ! In case src_q is zero just scale output properly
    if (src_q .eq. zero) then
        ! Zero init output if beta is also zero
        if (beta .eq. zero) then
            dst_m = zero
        ! Scale output by beta otherwise
        else
            dst_m = beta * dst_m
        end if
        ! Exit subroutine
        return
    end if
    ! Now src_q is non-zero
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    if (rho .eq. zero) then
        return
    end if
    ! Compute actual induced potentials
    rcoef = dst_r / rho
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        t = src_q * rcoef
        select case(p)
            case (0)
                if (beta .eq. zero) then
                    dst_m(1) = t * vscales_rel(1)
                else
                    dst_m(1) = beta*dst_m(1) + t*vscales_rel(1)
                end if
                return
            case (1)
                if (beta .eq. zero) then
                    ! l = 0
                    dst_m(1) = t * vscales_rel(1)
                    ! l = 1
                    t = t * rcoef
                    tmp = -vscales_rel(4) * stheta
                    tmp2 = t * tmp
                    dst_m(2) = tmp2 * sphi
                    dst_m(3) = t * vscales_rel(3) * ctheta
                    dst_m(4) = tmp2 * cphi
                else
                    ! l = 0
                    dst_m(1) = beta*dst_m(1) + t*vscales_rel(1)
                    ! l = 1
                    t = t * rcoef
                    tmp = -vscales_rel(4) * stheta
                    tmp2 = t * tmp
                    dst_m(2) = beta*dst_m(2) + tmp2*sphi
                    dst_m(3) = beta*dst_m(3) + t*vscales_rel(3)*ctheta
                    dst_m(4) = beta*dst_m(4) + tmp2*cphi
                end if
                return
        end select
        ! Mark first and last elements of all work subarrays
        vplm1 = 1
        vplm2 = (p+1)**2
        vplm => work(vplm1:vplm2)
        vcos1 = vplm2 + 1
        vcos2 = vcos1 + p
        vcos => work(vcos1:vcos2)
        vsin1 = vcos2 + 1
        vsin2 = vsin1 + p
        vsin => work(vsin1:vsin2)
        vminv1 = vsin2 + 1
        vminv2 = vminv1 + p-3
        vminv => work(vminv1:vminv2)
        ! Evaluate cos(m*phi) and sin(m*phi) arrays
        call trgev(cphi, sphi, p, vcos, vsin)
        ! Prepare -stheta*pt5/m
        tmp2 = -stheta * pt5
        do m = 1, p-2
            vminv(m) = tmp2 / dble(m)
        end do
        ! Overwrite output
        if (beta .eq. zero) then
            ! Construct spherical harmonics
            ! l = 0
            vplm(1) = one
            dst_m(1) = t * vscales_rel(1)
            ! l = 1
            t = t * rcoef
            vplm(3) = ctheta
            vplm(4) = -stheta
            tmp = t * vscales_rel(4) * vplm(4)
            dst_m(2) = tmp * sphi
            dst_m(3) = t * vscales_rel(3) * ctheta
            dst_m(4) = tmp * cphi
            indl = 3
            indlm1 = 1
            do l = 2, p
                t = t * rcoef
                ! Offset of a P_{l-2}^0 harmonic
                indlm2 = indlm1
                ! Offset of a P_{l-1}^0 harmonic
                indlm1 = indl
                ! Offset of a P_l^0 harmonic
                indl = indl + 2*l
                ! Some temp constants
                fl = dble(l)
                tmp1 = two*fl - one
                tmp2 = fl - one
                tmp3 = tmp1 * ctheta
                tmp4 = two * tmp2
                tmp5 = fl*fl - fl
                ! P_l^0
                vplm(indl) = tmp3*vplm(indlm1) - tmp2*vplm(indlm2)
                vplm(indl) = vplm(indl) / fl
                dst_m(indl) = t * vscales_rel(indl) * vplm(indl)
                ! P_l^m for m=1..l-2
                do m = 1, l-2
                    tmp4 = tmp4 + two
                    tmp5 = tmp5 + tmp4
                    ! P_l^m
                    vplm(indl+m) = vminv(m) * (vplm(indlm1+m+1) + &
                        & tmp5*vplm(indlm1+m-1))
                    tmp = t * vplm(indl+m) * vscales_rel(indl+m)
                    dst_m(indl+m) = tmp * vcos(m+1)
                    dst_m(indl-m) = tmp * vsin(m+1)
                end do
                ! m = l-1
                vplm(indl+l-1) = tmp3 * vplm(indlm1+l-1)
                tmp = t * vplm(indl+l-1) * vscales_rel(indl+l-1)
                dst_m(indl+l-1) = tmp * vcos(l)
                dst_m(indl-l+1) = tmp * vsin(l)
                ! m = l
                vplm(indl+l) = -stheta * tmp1 * vplm(indlm1+l-1)
                tmp = t * vplm(indl+l) * vscales_rel(indl+l)
                dst_m(indl+l) = tmp * vcos(l+1)
                dst_m(indl-l) = tmp * vsin(l+1)
            end do
        else
            ! Construct spherical harmonics
            ! l = 0
            vplm(1) = one
            dst_m(1) = beta*dst_m(1) + t*vscales_rel(1)
            ! l = 1
            t = t * rcoef
            vplm(3) = ctheta
            vplm(4) = -stheta
            tmp = t * vscales_rel(4) * vplm(4)
            dst_m(2) = beta*dst_m(2) + tmp*sphi
            dst_m(3) = beta*dst_m(3) + t*vscales_rel(3)*ctheta
            dst_m(4) = beta*dst_m(4) + tmp*cphi
            indl = 3
            indlm1 = 1
            do l = 2, p
                t = t * rcoef
                ! Offset of a P_{l-2}^0 harmonic
                indlm2 = indlm1
                ! Offset of a P_{l-1}^0 harmonic
                indlm1 = indl
                ! Offset of a P_l^0 harmonic
                indl = indl + 2*l
                ! Some temp constants
                fl = dble(l)
                tmp1 = two*fl - one
                tmp2 = fl - one
                tmp3 = tmp1 * ctheta
                tmp4 = two * tmp2
                tmp5 = fl*fl - fl
                ! Y_l^0
                vplm(indl) = tmp3*vplm(indlm1) - tmp2*vplm(indlm2)
                vplm(indl) = vplm(indl) / fl
                dst_m(indl) = beta*dst_m(indl) + t*vscales_rel(indl)*vplm(indl)
                ! Y_l^m for m=1..l-2
                do m = 1, l-2
                    tmp4 = tmp4 + two
                    tmp5 = tmp5 + tmp4
                    ! P_l^m
                    vplm(indl+m) = vminv(m) * (vplm(indlm1+m+1) + &
                        & tmp5*vplm(indlm1+m-1))
                    tmp = t * vplm(indl+m) * vscales_rel(indl+m)
                    dst_m(indl+m) = beta*dst_m(indl+m) + tmp*vcos(m+1)
                    dst_m(indl-m) = beta*dst_m(indl-m) + tmp*vsin(m+1)
                end do
                ! m = l-1
                vplm(indl+l-1) = tmp3 * vplm(indlm1+l-1)
                tmp = t * vplm(indl+l-1) * vscales_rel(indl+l-1)
                dst_m(indl+l-1) = beta*dst_m(indl+l-1) + tmp*vcos(l)
                dst_m(indl-l+1) = beta*dst_m(indl-l+1) + tmp*vsin(l)
                ! m = l
                vplm(indl+l) = -stheta * tmp1 * vplm(indlm1+l-1)
                tmp = t * vplm(indl+l) * vscales_rel(indl+l)
                dst_m(indl+l) = beta*dst_m(indl+l) + tmp*vcos(l+1)
                dst_m(indl-l) = beta*dst_m(indl-l) + tmp*vsin(l+1)
            end do
        end if
    ! Case of x(1:2) = 0 and x(3) != 0
    else
    ! TODO: Add proper check for this routine, I am leaving here stop
    ! instruction to indicate that it is not tested
        stop "fmm_m2p_adj_work untested branch"
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^m = ctheta^m in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef.
        rcoef = sign(rcoef, c(3))
        ! Init t and proceed with accumulation of a potential
        t = src_q
        indl = 1
        if (beta .eq. zero) then
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Update t
                t = t * rcoef
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_m(indl) = t * vscales_rel(indl)
                dst_m(indl-l:indl-1) = zero
                dst_m(indl+1:indl+l) = zero
            end do
        else
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Update t
                t = t * rcoef
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_m(indl) = beta*dst_m(indl) + t*vscales_rel(indl)
                dst_m(indl-l:indl-1) = beta * dst_m(indl-l:indl-1)
                dst_m(indl+1:indl+l) = beta * dst_m(indl+1:indl+l)
            end do
        end if
    end if
end subroutine fmm_m2p_adj_work

!> Accumulate potentials, induced by each multipole spherical harmonic
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[inout] mat: Values of potentials induced by each spherical harmonic
subroutine fmm_m2p_mat(c, r, p, vscales, mat)
    ! Inputs
    real(dp), intent(in) :: c(3), r, vscales((p+1)**2)
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: mat((p+1)**2)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), rcoef, t, tmp
    integer :: n, ind
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    if (rho .eq. zero) then
        return
    end if
    ! Compute actual induced potentials
    rcoef = r / rho
    t = one
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        mat(ind-n:ind+n) = t / vscales(ind)**2 * vylm(ind-n:ind+n)
    end do
end subroutine fmm_m2p_mat

!> Accumulate a local expansion induced by a particle of a given charge
!!
!! Computes the following sums:
!! \f[
!!      \forall \ell=0, \ldots, p, \quad \forall m=-\ell, \ldots, \ell : \quad
!!      M_\ell^m = \beta M_\ell^m + \frac{q \|c\|^\ell}{r^{\ell+1}}
!!      Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of output harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \beta \f$ is a scaling factor, \f$ q \f$ and \f$
!! c \f$ are a charge and coordinates of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{\ell+1} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^{-\ell-1} \f$.
!!
!! @param[in] c: Radius-vector from the particle to the center of harmonics
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output local spherical harmonics
!! @param[in] p: Maximal degree of output local spherical harmonics
!! @param[in] vscales: Normalization constants for spherical harmonics
!! @param[in] beta: Scaling factor for `dst_l`
!! @param[inout] dst_l: Local coefficients
!!
!! @sa fmm_l2p
subroutine fmm_p2l(c, src_q, dst_r, p, vscales, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)**2)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), t, rcoef
    integer :: n, ind
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Local expansion represents potential from the outside of the sphere of
    ! the given radius `dst_r`. In case `rho=zero` source particle is inside of
    ! the sphere no matter what radius `r` is, so we just ignore such a case.
    if (rho .eq. zero) then
        if (beta .eq. zero) then
            dst_l = zero
        else
            dst_l = beta * dst_l
        end if
        return
    end if
    rcoef = dst_r / rho
    t = src_q / rho
    ! Ignore input `m` in case of zero scaling factor
    if (beta .eq. zero) then
        do n = 0, p
            ind = n*n + n + 1
            dst_l(ind-n:ind+n) = t * vylm(ind-n:ind+n)
            t = t * rcoef
        end do
    ! Update `m` otherwise
    else
        do n = 0, p
            ind = n*n + n + 1
            dst_l(ind-n:ind+n) = beta*dst_l(ind-n:ind+n) + t*vylm(ind-n:ind+n)
            t = t * rcoef
        end do
    end if
end subroutine fmm_p2l

!> Accumulate potential, induced by local spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Local coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_v`
!! @param[inout] dst_v: Value of induced potential
subroutine fmm_l2p(c, src_r, p, vscales, alpha, src_l, beta, dst_v)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, vscales((p+1)*(p+1)), alpha, &
        & src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Workspace
    real(dp) :: work((p+1)*(p+1)+3*p)
    ! Call corresponding work routine
    call fmm_l2p_work(c, src_r, p, vscales, alpha, src_l, beta, dst_v, work)
end subroutine fmm_l2p

!> Accumulate potential, induced by local spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Local coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_v`
!! @param[inout] dst_v: Value of induced potential
!! @param[out] work: Temporary workspace of a size (p+1)**2+3*p
subroutine fmm_l2p_work(c, src_r, p, vscales_rel, alpha, src_l, beta, dst_v, &
        & work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, vscales_rel((p+1)*(p+1)), alpha, &
        & src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Workspace
    real(dp), intent(out), target :: work((p+1)*(p+1)+3*p)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, tmp4, tmp5, fl, max12, ssq12
    integer :: vplm1, vplm2, vcos1, vcos2, vsin1, vsin2, vminv1, vminv2, &
        & l, m, indl, indlm1, indlm2
    real(dp), pointer :: vplm(:), vcos(:), vsin(:), vminv(:)
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential depends only on the
    ! leading 1-st spherical harmonic
    if (rho .eq. zero) then
        dst_v = dst_v + alpha*src_l(1)*vscales_rel(1)
        return
    end if
    ! Compute the actual induced potential
    rcoef = rho / src_r
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                ! l = 0
                dst_v = dst_v + alpha*vscales_rel(1)*src_l(1)
                return
            case (1)
                ! l = 0
                tmp = src_l(1) * vscales_rel(1)
                ! l = 1
                tmp2 = ctheta * src_l(3) * vscales_rel(3)
                tmp3 = vscales_rel(4) * stheta
                tmp2 = tmp2 - tmp3*sphi*src_l(2)
                tmp2 = tmp2 - tmp3*cphi*src_l(4)
                dst_v = dst_v + alpha*(tmp+rcoef*tmp2)
                return
        end select
        ! Mark first and last elements of all work subarrays
        vplm1 = 1
        vplm2 = (p+1)**2
        vplm => work(vplm1:vplm2)
        vcos1 = vplm2 + 1
        vcos2 = vcos1 + p
        vcos => work(vcos1:vcos2)
        vsin1 = vcos2 + 1
        vsin2 = vsin1 + p
        vsin => work(vsin1:vsin2)
        vminv1 = vsin2 + 1
        vminv2 = vminv1 + p-3
        vminv => work(vminv1:vminv2)
        ! Evaluate cos(m*phi) and sin(m*phi) arrays
        call trgev(cphi, sphi, p, vcos, vsin)
        ! Prepare -stheta*pt5/m
        tmp2 = -stheta * pt5
        do m = 1, p-2
            vminv(m) = tmp2 / dble(m)
        end do
        ! Construct spherical harmonics
        ! l = 0
        vplm(1) = one
        tmp = src_l(1) * vscales_rel(1)
        ! l = 1
        vplm(3) = ctheta
        vplm(4) = -stheta
        tmp2 = ctheta * src_l(3) * vscales_rel(3)
        tmp2 = tmp2 - stheta*vscales_rel(4)*(sphi*src_l(2)+cphi*src_l(4))
        tmp = tmp + rcoef*tmp2
        t = rcoef
        indl = 3
        indlm1 = 1
        do l = 2, p
            t = t * rcoef
            ! Offset of a P_{l-2}^0 harmonic
            indlm2 = indlm1
            ! Offset of a P_{l-1}^0 harmonic
            indlm1 = indl
            ! Offset of a P_l^0 harmonic
            indl = indl + 2*l
            ! Some temp constants
            fl = dble(l)
            tmp1 = two*fl - one
            tmp2 = fl - one
            tmp3 = tmp1 * ctheta
            tmp4 = two * tmp2
            tmp5 = fl*fl - fl
            ! Y_l^0
            vplm(indl) = tmp3*vplm(indlm1) - tmp2*vplm(indlm2)
            vplm(indl) = vplm(indl) / fl
            tmp2 = vscales_rel(indl) * vplm(indl) * src_l(indl)
            ! Y_l^m for m=1..l-2
            do m = 1, l-2
                tmp4 = tmp4 + two
                tmp5 = tmp5 + tmp4
                ! P_l^m
                vplm(indl+m) = vminv(m) * (vplm(indlm1+m+1) + &
                    & tmp5*vplm(indlm1+m-1))
                tmp2 = tmp2 + vplm(indl+m)*(vcos(m+1)*src_l(indl+m)+ &
                    & vsin(m+1)*src_l(indl-m))*vscales_rel(indl+m)
            end do
            ! m = l-1
            vplm(indl+l-1) = tmp3 * vplm(indlm1+l-1)
            tmp2 = tmp2 + vplm(indl+l-1)*(vcos(l)*src_l(indl+l-1)+ &
                & vsin(l)*src_l(indl-l+1))*vscales_rel(indl+l-1)
            ! m = l
            vplm(indl+l) = -stheta * tmp1 * vplm(indlm1+l-1)
            tmp2 = tmp2 + vplm(indl+l)*(vcos(l+1)*src_l(indl+l)+ &
                & vsin(l+1)*src_l(indl-l))*vscales_rel(indl+l)
            tmp = tmp + t*tmp2
        end do
        dst_v = dst_v + alpha*tmp
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^m = ctheta^m in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef.
        rcoef = sign(rcoef, c(3))
        ! Init t and proceed with accumulation of a potential
        t = alpha
        indl = 1
        do l = 0, p
            ! Index of Y_l^0
            indl = indl + 2*l
            ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
            dst_v = dst_v + t*src_l(indl)*vscales_rel(indl)
            ! Update t
            t = t * rcoef
        end do
    end if
end subroutine fmm_l2p_work

!> Adjoint of L2P
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output local spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Local coefficients. Dimension is `(p+1)**2`
subroutine fmm_l2p_adj(c, src_q, dst_r, p, vscales_rel, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales_rel((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)**2)
    ! Workspace
    real(dp) :: work((p+1)*(p+1)+3*p)
    ! Call corresponding work routine
    call fmm_l2p_adj_work(c, src_q, dst_r, p, vscales_rel, beta, dst_l, work)
end subroutine fmm_l2p_adj

!> Adjoint L2P
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output local spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Local coefficients. Dimension is `(p+1)**2`
!! @param[out] work: Temporary workspace of a size (p+1)**2+3*p
subroutine fmm_l2p_adj_work(c, src_q, dst_r, p, vscales_rel, beta, dst_l, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales_rel((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)**2)
    ! Workspace
    real(dp), intent(out), target :: work((p+1)*(p+1)+3*p)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, tmp4, tmp5, fl, max12, ssq12
    integer :: vplm1, vplm2, vcos1, vcos2, vsin1, vsin2, vminv1, vminv2, &
        & l, m, indl, indlm1, indlm2
    real(dp), pointer :: vplm(:), vcos(:), vsin(:), vminv(:)
    ! In case src_q is zero just scale output properly
    if (src_q .eq. zero) then
        ! Zero init output if beta is also zero
        if (beta .eq. zero) then
            dst_l = zero
        ! Scale output by beta otherwise
        else
            dst_l = beta * dst_l
        end if
        ! Exit subroutine
        return
    end if
    ! Now src_q is non-zero
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential depends only on the
    ! leading 1-st spherical harmonic, so in adjoint way we update only the leading
    if (rho .eq. zero) then
        if (beta .eq. zero) then
            dst_l(1) = src_q * vscales_rel(1)
            dst_l(2:) = zero
        else
            dst_l(1) = beta*dst_l(1) + src_q*vscales_rel(1)
            dst_l(2:) = beta * dst_l(2:)
        end if
        return
    end if
    ! Compute the actual induced potential
    rcoef = rho / dst_r
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                ! l = 0
                if (beta .eq. zero) then
                    dst_l(1) = src_q * vscales_rel(1)
                else
                    dst_l(1) = beta*dst_l(1) + src_q*vscales_rel(1)
                end if
                return
            case (1)
                ! l = 0 and l = 1
                tmp = src_q * rcoef
                tmp2 = -tmp * vscales_rel(4) * stheta
                if (beta .eq. zero) then
                    dst_l(1) = src_q * vscales_rel(1)
                    dst_l(2) = tmp2 * sphi
                    dst_l(3) = tmp * vscales_rel(3) * ctheta
                    dst_l(4) = tmp2 * cphi
                else
                    dst_l(1) = beta*dst_l(1) + src_q*vscales_rel(1)
                    dst_l(2) = beta*dst_l(2) + tmp2*sphi
                    dst_l(3) = beta*dst_l(3) + tmp*vscales_rel(3)*ctheta
                    dst_l(4) = beta*dst_l(4) + tmp2*cphi
                end if
                return
        end select
        ! Mark first and last elements of all work subarrays
        vplm1 = 1
        vplm2 = (p+1)**2
        vplm => work(vplm1:vplm2)
        vcos1 = vplm2 + 1
        vcos2 = vcos1 + p
        vcos => work(vcos1:vcos2)
        vsin1 = vcos2 + 1
        vsin2 = vsin1 + p
        vsin => work(vsin1:vsin2)
        vminv1 = vsin2 + 1
        vminv2 = vminv1 + p-3
        vminv => work(vminv1:vminv2)
        ! Evaluate cos(m*phi) and sin(m*phi) arrays
        call trgev(cphi, sphi, p, vcos, vsin)
        ! Prepare -stheta*pt5/m
        tmp2 = -stheta * pt5
        do m = 1, p-2
            vminv(m) = tmp2 / dble(m)
        end do
        ! Overwrite output
        if (beta .eq. zero) then
            ! Construct spherical harmonics
            ! l = 0
            vplm(1) = one
            dst_l(1) = src_q * vscales_rel(1)
            ! l = 1
            vplm(3) = ctheta
            vplm(4) = -stheta
            t = src_q * rcoef
            tmp = -t * vscales_rel(4) * stheta
            dst_l(2) = tmp * sphi
            dst_l(3) = t * vscales_rel(3) * ctheta
            dst_l(4) = tmp * cphi
            indl = 3
            indlm1 = 1
            do l = 2, p
                t = t * rcoef
                ! Offset of a P_{l-2}^0 harmonic
                indlm2 = indlm1
                ! Offset of a P_{l-1}^0 harmonic
                indlm1 = indl
                ! Offset of a P_l^0 harmonic
                indl = indl + 2*l
                ! Some temp constants
                fl = dble(l)
                tmp1 = two*fl - one
                tmp2 = fl - one
                tmp3 = tmp1 * ctheta
                tmp4 = two * tmp2
                tmp5 = fl*fl - fl
                ! P_l^0
                vplm(indl) = tmp3*vplm(indlm1) - tmp2*vplm(indlm2)
                vplm(indl) = vplm(indl) / fl
                dst_l(indl) = t * vscales_rel(indl) * vplm(indl)
                ! P_l^m for m=1..l-2
                do m = 1, l-2
                    tmp4 = tmp4 + two
                    tmp5 = tmp5 + tmp4
                    ! P_l^m
                    vplm(indl+m) = vminv(m) * (vplm(indlm1+m+1) + &
                        & tmp5*vplm(indlm1+m-1))
                    tmp = t * vplm(indl+m) * vscales_rel(indl+m)
                    dst_l(indl+m) = tmp * vcos(m+1)
                    dst_l(indl-m) = tmp * vsin(m+1)
                end do
                ! m = l-1
                vplm(indl+l-1) = tmp3 * vplm(indlm1+l-1)
                tmp = t * vplm(indl+l-1) * vscales_rel(indl+l-1)
                dst_l(indl+l-1) = tmp * vcos(l)
                dst_l(indl-l+1) = tmp * vsin(l)
                ! m = l
                vplm(indl+l) = -stheta * tmp1 * vplm(indlm1+l-1)
                tmp = t * vplm(indl+l) * vscales_rel(indl+l)
                dst_l(indl+l) = tmp * vcos(l+1)
                dst_l(indl-l) = tmp * vsin(l+1)
            end do
        else
            ! Construct spherical harmonics
            ! l = 0
            vplm(1) = one
            dst_l(1) = beta*dst_l(1) + src_q*vscales_rel(1)
            ! l = 1
            vplm(3) = ctheta
            vplm(4) = -stheta
            t = src_q * rcoef
            tmp = -t * vscales_rel(4) * stheta
            dst_l(2) = beta*dst_l(2) + tmp*sphi
            dst_l(3) = beta*dst_l(3) + t*vscales_rel(3)*ctheta
            dst_l(4) = beta*dst_l(4) + tmp*cphi
            indl = 3
            indlm1 = 1
            do l = 2, p
                t = t * rcoef
                ! Offset of a P_{l-2}^0 harmonic
                indlm2 = indlm1
                ! Offset of a P_{l-1}^0 harmonic
                indlm1 = indl
                ! Offset of a P_l^0 harmonic
                indl = indl + 2*l
                ! Some temp constants
                fl = dble(l)
                tmp1 = two*fl - one
                tmp2 = fl - one
                tmp3 = tmp1 * ctheta
                tmp4 = two * tmp2
                tmp5 = fl*fl - fl
                ! P_l^0
                vplm(indl) = tmp3*vplm(indlm1) - tmp2*vplm(indlm2)
                vplm(indl) = vplm(indl) / fl
                dst_l(indl) = beta*dst_l(indl) + t*vscales_rel(indl)*vplm(indl)
                ! P_l^m for m=1..l-2
                do m = 1, l-2
                    tmp4 = tmp4 + two
                    tmp5 = tmp5 + tmp4
                    ! P_l^m
                    vplm(indl+m) = vminv(m) * (vplm(indlm1+m+1) + &
                        & tmp5*vplm(indlm1+m-1))
                    tmp = t * vplm(indl+m) * vscales_rel(indl+m)
                    dst_l(indl+m) = beta*dst_l(indl+m) + tmp*vcos(m+1)
                    dst_l(indl-m) = beta*dst_l(indl-m) + tmp*vsin(m+1)
                end do
                ! m = l-1
                vplm(indl+l-1) = tmp3 * vplm(indlm1+l-1)
                tmp = t * vplm(indl+l-1) * vscales_rel(indl+l-1)
                dst_l(indl+l-1) = beta*dst_l(indl+l-1) + tmp*vcos(l)
                dst_l(indl-l+1) = beta*dst_l(indl-l+1) + tmp*vsin(l)
                ! m = l
                vplm(indl+l) = -stheta * tmp1 * vplm(indlm1+l-1)
                tmp = t * vplm(indl+l) * vscales_rel(indl+l)
                dst_l(indl+l) = beta*dst_l(indl+l) + tmp*vcos(l+1)
                dst_l(indl-l) = beta*dst_l(indl-l) + tmp*vsin(l+1)
            end do
        end if
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^m = ctheta^m in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef.
        rcoef = sign(rcoef, c(3))
        ! Init t and proceed with accumulation of a potential
        t = src_q
        indl = 1
        ! Overwrite output
        if (beta .eq. zero) then
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_l(indl) = t * vscales_rel(indl)
                dst_l(indl-l:indl-1) = zero
                dst_l(indl+1:indl+l) = zero
                ! Update t
                t = t * rcoef
            end do
        ! Update output
        else
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_l(indl) = beta*dst_l(indl) + t*vscales_rel(indl)
                dst_l(indl-l:indl-1) = beta * dst_l(indl-l:indl-1)
                dst_l(indl+1:indl+l) = beta * dst_l(indl+1:indl+l)
                ! Update t
                t = t * rcoef
            end do
        end if
    end if
end subroutine fmm_l2p_adj_work

!> Transform coefficients of spherical harmonics to a new cartesion system
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$ R \f$
!! is a '(p+1)**2'-by-'(p+1)**2' matrix that transforms coefficients of the
!! source spherical harmonics to the destination spherical harmonics while
!! preserving values of any function on a unit sphere defined as a linear
!! combination of the source spherical harmonics.
!!
!! Input 3-by-3 matrix `r1` must be an orthogonal matrix \f$ R_1 \f$ of a
!! transform of new cartesion coordinates \f$ (\widetilde{y}, \widetilde{z},
!! \widetilde{x}) \f$ into initial cartesian coordinates \f$ (y, z, x) \f$.
!! This is due to the following equalities:
!! \f{align}{
!!      Y_1^{-1} (\theta, \phi) &= \sqrt{\frac{3}{4\pi}} \sin \theta \sin \phi
!!      = \sqrt{\frac{3}{4\pi}} y, \\ Y_1^0 (\theta, \phi) &=
!!      \sqrt{\frac{3}{4\pi}} \cos \theta = \sqrt{\frac{3}{4\pi}} z, \\
!!      Y_1^1 (\theta, \phi) &= \sqrt{\frac{3}{4\pi}} \sin \theta \cos \phi =
!!      \sqrt{\frac{3}{4\pi}} x.
!! \f}
!! So, to find a column-vector \f$ \widetilde{c} \f$ of coefficients of
!! spherical harmonics \f$ Y_1^{-1}, Y_1^0 \f$ and \f$ Y_1^1 \f$ in a new
!! system of coordinates \f$ (\widetilde{y}, \widetilde{z}, \widetilde{x}) \f$
!! the following system needs to be solved:
!! \f[
!!      \widetilde{c}^\top \cdot \begin{bmatrix} Y_1^{-1}
!!      (\widetilde{\theta}, \widetilde{\phi}) \\ Y_1^0 (\widetilde{\theta},
!!      \widetilde{\phi}) \\ Y_1^1 (\widetilde{\theta}, \widetilde{\phi})
!!      \end{bmatrix} = c ^\top \cdot \begin{bmatrix} Y_1^{-1} (\theta, \phi)
!!      \\ Y_1^0 (\theta, \phi) \\ Y_1^1 (\theta, \phi) \end{bmatrix}.
!! \f]
!! The solution has the following obvious form:
!! \f[
!!      \widetilde{c} = R_1^\top c.
!! \f]
!!
!! Translation of spherical harmonics of all other degrees is computed
!! recursively as is described in the following source:
!!      @cite ir-realharms-1996
!!      @cite ir-realharms-1998
!!
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] r1: Transformation from new to old cartesian coordinates
!! @param[in] alpha: Scalar multipler for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of transformed spherical harmonics
!!
!! @sa fmm_sph_transform_get_mat, fmm_sph_transform_use_mat
subroutine fmm_sph_transform(p, r1, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1(-1:1, -1:1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(2*p+1)*(2*p+3))
    ! Call corresponding work routine
    call fmm_sph_transform_work(p, r1, alpha, src, beta, dst, work)
end subroutine fmm_sph_transform

!> Transform spherical harmonics to a new cartesion system of coordinates
!!
!! This function implements @ref fmm_sph_transform with predefined values of
!! parameters \p alpha=one and \p beta=zero.
!! 
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] r1: Transformation from new to old cartesian coordinates
!! @param[in] alpha: Scalar multipler for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of transformed spherical harmonics
!! @param[out] work: Temporary workspace of a size (2*(2*p+1)*(2*p+3))
!!
!! @sa fmm_sph_transform_get_mat, fmm_sph_transform_use_mat
subroutine fmm_sph_transform_work(p, r1, alpha, src, beta, dst, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1(-1:1, -1:1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(out) :: dst((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(2*p+1)*(2*p+3))
    ! Local variables
    real(dp) :: u, v, w
    integer :: l, m, n, ind
    ! Pointers for a workspace
    real(dp), pointer :: r_prev(:, :), r(:, :), scal_uvw_m(:), scal_u_n(:), &
        & scal_v_n(:), scal_w_n(:), r_swap(:, :)
    ! In case alpha is zero just scale output
    if (alpha .eq. zero) then
        ! Set output to zero if beta is also zero
        if (beta .eq. zero) then
            dst = zero
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-zero
    ! In case beta is zero output is just overwritten without being read
    if (beta .eq. zero) then
        ! Compute rotations/reflections
        ! l = 0
        dst(1) = alpha * src(1)
        if (p .eq. 0) then
            return
        end if
        ! l = 1
        dst(2) = alpha * dot_product(r1(-1:1,-1), src(2:4))
        dst(3) = alpha * dot_product(r1(-1:1,0), src(2:4))
        dst(4) = alpha * dot_product(r1(-1:1,1), src(2:4))
        if (p .eq. 1) then
            return
        end if
        ! Set pointers
        n = 2*p + 1
        l = n ** 2
        r_prev(-p:p, -p:p) => work(1:l)
        m = 2*l
        r(-p:p, -p:p) => work(l+1:m)
        l = m + n
        scal_uvw_m(-p:p) => work(m+1:l)
        m = l + n
        scal_u_n(-p:p) => work(l+1:m)
        l = m + n
        scal_v_n(-p:p) => work(m+1:l)
        m = l + n
        scal_w_n(-p:p) => work(l+1:m)
        ! l = 2
        r(2, 2) = (r1(1, 1)*r1(1, 1) - r1(1, -1)*r1(1, -1) - &
            & r1(-1, 1)*r1(-1, 1) + r1(-1, -1)*r1(-1, -1)) / two
        r(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
        r(2, 0) = sqrt3 / two * (r1(1, 0)*r1(1, 0) - r1(-1, 0)*r1(-1, 0))
        r(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
        r(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
        r(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
        r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
        r(1, 0) = sqrt3 * r1(1, 0) * r1(0, 0)
        r(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
        r(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
        r(0, 2) = sqrt3 / two * (r1(0, 1)*r1(0, 1) - r1(0, -1)*r1(0, -1))
        r(0, 1) = sqrt3 * r1(0, 1) * r1(0, 0)
        r(0, 0) = (three*r1(0, 0)*r1(0, 0)-one) / two
        r(0, -1) = sqrt3 * r1(0, -1) * r1(0, 0)
        r(0, -2) = sqrt3 * r1(0, 1) * r1(0, -1)
        r(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
        r(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
        r(-1, 0) = sqrt3 * r1(-1, 0) * r1(0, 0)
        r(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(0, -1)*r1(-1, 0)
        r(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
        r(-2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
        r(-2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
        r(-2, 0) = sqrt3 * r1(1, 0) * r1(-1, 0)
        r(-2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
        r(-2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
        do m = -2, 2
            dst(7+m) = alpha * dot_product(r(-2:2, m), src(5:9))
        end do
        ! l > 2
        do l = 3, p
            ! Swap previous and current rotation matrices
            r_swap => r_prev
            r_prev => r
            r => r_swap
            ! Prepare scalar factors
            scal_uvw_m(0) = dble(l)
            do m = 1, l-1
                u = sqrt(dble(l*l-m*m))
                scal_uvw_m(m) = u
                scal_uvw_m(-m) = u
            end do
            u = two * dble(l)
            u = sqrt(dble(u*(u-one)))
            scal_uvw_m(l) = u
            scal_uvw_m(-l) = u
            scal_u_n(0) = dble(l)
            scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
            scal_w_n(0) = zero
            do n = 1, l-2
                u = sqrt(dble(l*l-n*n))
                scal_u_n(n) = u
                scal_u_n(-n) = u
                v = dble(l+n)
                v = sqrt(v*(v-one)) / two
                scal_v_n(n) = v
                scal_v_n(-n) = v
                w = dble(l-n)
                w = -sqrt(w*(w-one)) / two
                scal_w_n(n) = w
                scal_w_n(-n) = w
            end do
            u = sqrt(dble(2*l-1))
            scal_u_n(l-1) = u
            scal_u_n(1-l) = u
            scal_u_n(l) = zero
            scal_u_n(-l) = zero
            v = sqrt(dble((2*l-1)*(2*l-2))) / two
            scal_v_n(l-1) = v
            scal_v_n(1-l) = v
            v = sqrt(dble(2*l*(2*l-1))) / two
            scal_v_n(l) = v
            scal_v_n(-l) = v
            scal_w_n(l-1) = zero
            scal_w_n(l) = zero
            scal_w_n(-l) = zero
            scal_w_n(1-l) = zero
            ind = l*l + l + 1
            ! m = l, n = l
            v = r1(1, 1)*r_prev(l-1, l-1) - r1(1, -1)*r_prev(l-1, 1-l) - &
                & r1(-1, 1)*r_prev(1-l, l-1) + r1(-1, -1)*r_prev(1-l, 1-l)
            r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind+l) = src(ind+l) * r(l, l)
            ! m = l, n = -l
            v = r1(1, 1)*r_prev(1-l, l-1) - r1(1, -1)*r_prev(1-l, 1-l) + &
                & r1(-1, 1)*r_prev(l-1, l-1) - r1(-1, -1)*r_prev(l-1, 1-l)
            r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
            ! m = l, n = l-1
            u = r1(0, 1)*r_prev(l-1, l-1) - r1(0, -1)*r_prev(l-1, 1-l)
            v = r1(1, 1)*r_prev(l-2, l-1) - r1(1, -1)*r_prev(l-2, 1-l) - &
                & r1(-1, 1)*r_prev(2-l, l-1) + r1(-1, -1)*r_prev(2-l, 1-l)
            r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
            ! m = l, n = 1-l
            u = r1(0, 1)*r_prev(1-l, l-1) - r1(0, -1)*r_prev(1-l, 1-l)
            v = r1(1, 1)*r_prev(2-l, l-1) - r1(1, -1)*r_prev(2-l, 1-l) + &
                & r1(-1, 1)*r_prev(l-2, l-1) - r1(-1, -1)*r_prev(l-2, 1-l)
            r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
            ! m = l, n = 1
            u = r1(0, 1)*r_prev(1, l-1) - r1(0, -1)*r_prev(1, 1-l)
            v = r1(1, 1)*r_prev(0, l-1) - r1(1, -1)*r_prev(0, 1-l)
            w = r1(1, 1)*r_prev(2, l-1) - r1(1, -1)*r_prev(2, 1-l) + &
                & r1(-1, 1)*r_prev(-2, l-1) - r1(-1, -1)*r_prev(-2, 1-l)
            r(1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, l) = r(1, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
            ! m = l, n = -1
            u = r1(0, 1)*r_prev(-1, l-1) - r1(0, -1)*r_prev(-1, 1-l)
            v = r1(-1, 1)*r_prev(0, l-1) - r1(-1, -1)*r_prev(0, 1-l)
            w = r1(1, 1)*r_prev(-2, l-1) - r1(1, -1)*r_prev(-2, 1-l) - &
                & r1(-1, 1)*r_prev(2, l-1) + r1(-1, -1)*r_prev(2, 1-l)
            r(-1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, l) = r(-1, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
            ! m = l, n = 0
            u = r1(0, 1)*r_prev(0, l-1) - r1(0, -1)*r_prev(0, 1-l)
            v = r1(1, 1)*r_prev(1, l-1) - r1(1, -1)*r_prev(1, 1-l) + &
                & r1(-1, 1)*r_prev(-1, l-1) - r1(-1, -1)*r_prev(-1, 1-l)
            r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
            ! m = l, n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 1)*r_prev(n, l-1) - r1(0, -1)*r_prev(n, 1-l)
                v = r1(1, 1)*r_prev(n-1, l-1) - r1(1, -1)*r_prev(n-1, 1-l) - &
                    & r1(-1, 1)*r_prev(1-n, l-1) + r1(-1, -1)*r_prev(1-n, 1-l)
                w = r1(1, 1)*r_prev(n+1, l-1) - r1(1, -1)*r_prev(n+1, 1-l) + &
                    & r1(-1, 1)*r_prev(-n-1, l-1) - r1(-1, -1)*r_prev(-n-1, 1-l)
                r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, l) = r(n, l) / scal_uvw_m(l)
                !dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
                u = r1(0, 1)*r_prev(-n, l-1) - r1(0, -1)*r_prev(-n, 1-l)
                v = r1(1, 1)*r_prev(1-n, l-1) - r1(1, -1)*r_prev(1-n, 1-l) + &
                    & r1(-1, 1)*r_prev(n-1, l-1) - r1(-1, -1)*r_prev(n-1, 1-l)
                w = r1(1, 1)*r_prev(-n-1, l-1) - r1(1, -1)*r_prev(-n-1, 1-l) - &
                    & r1(-1, 1)*r_prev(n+1, l-1) + r1(-1, -1)*r_prev(n+1, 1-l)
                r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, l) = r(-n, l) / scal_uvw_m(l)
                !dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
            end do
            dst(ind+l) = alpha * dot_product(src(ind-l:ind+l), r(-l:l, l))
            ! m = -l, n = l
            v = r1(1, 1)*r_prev(l-1, 1-l) + r1(1, -1)*r_prev(l-1, l-1) - &
                & r1(-1, 1)*r_prev(1-l, 1-l) - r1(-1, -1)*r_prev(1-l, l-1)
            r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind-l) = src(ind+l)*r(l, -l)
            ! m = -l, n = -l
            v = r1(1, 1)*r_prev(1-l, 1-l) + r1(1, -1)*r_prev(1-l, l-1) + &
                & r1(-1, 1)*r_prev(l-1, 1-l) + r1(-1, -1)*r_prev(l-1, l-1)
            r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
            ! m = -l, n = l-1
            u = r1(0, 1)*r_prev(l-1, 1-l) + r1(0, -1)*r_prev(l-1, l-1)
            v = r1(1, 1)*r_prev(l-2, 1-l) + r1(1, -1)*r_prev(l-2, l-1) - &
                & r1(-1, 1)*r_prev(2-l, 1-l) - r1(-1, -1)*r_prev(2-l, l-1)
            r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
            ! m = -l, n = 1-l
            u = r1(0, 1)*r_prev(1-l, 1-l) + r1(0, -1)*r_prev(1-l, l-1)
            v = r1(1, 1)*r_prev(2-l, 1-l) + r1(1, -1)*r_prev(2-l, l-1) + &
                & r1(-1, 1)*r_prev(l-2, 1-l) + r1(-1, -1)*r_prev(l-2, l-1)
            r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
            ! m = -l, n = 1
            u = r1(0, 1)*r_prev(1, 1-l) + r1(0, -1)*r_prev(1, l-1)
            v = r1(1, 1)*r_prev(0, 1-l) + r1(1, -1)*r_prev(0, l-1)
            w = r1(1, 1)*r_prev(2, 1-l) + r1(1, -1)*r_prev(2, l-1) + &
                & r1(-1, 1)*r_prev(-2, 1-l) + r1(-1, -1)*r_prev(-2, l-1)
            r(1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, -l) = r(1, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
            ! m = -l, n = -1
            u = r1(0, 1)*r_prev(-1, 1-l) + r1(0, -1)*r_prev(-1, l-1)
            v = r1(-1, 1)*r_prev(0, 1-l) + r1(-1, -1)*r_prev(0, l-1)
            w = r1(1, 1)*r_prev(-2, 1-l) + r1(1, -1)*r_prev(-2, l-1) - &
                & r1(-1, 1)*r_prev(2, 1-l) - r1(-1, -1)*r_prev(2, l-1)
            r(-1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
            ! m = -l, n = 0
            u = r1(0, 1)*r_prev(0, 1-l) + r1(0, -1)*r_prev(0, l-1)
            v = r1(1, 1)*r_prev(1, 1-l) + r1(1, -1)*r_prev(1, l-1) + &
                & r1(-1, 1)*r_prev(-1, 1-l) + r1(-1, -1)*r_prev(-1, l-1)
            r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
            ! m = -l, n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 1)*r_prev(n, 1-l) + r1(0, -1)*r_prev(n, l-1)
                v = r1(1, 1)*r_prev(n-1, 1-l) + r1(1, -1)*r_prev(n-1, l-1) - &
                    & r1(-1, 1)*r_prev(1-n, 1-l) - r1(-1, -1)*r_prev(1-n, l-1)
                w = r1(1, 1)*r_prev(n+1, 1-l) + r1(1, -1)*r_prev(n+1, l-1) + &
                    & r1(-1, 1)*r_prev(-n-1, 1-l) + r1(-1, -1)*r_prev(-n-1, l-1)
                r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, -l) = r(n, -l) / scal_uvw_m(l)
                !dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
                u = r1(0, 1)*r_prev(-n, 1-l) + r1(0, -1)*r_prev(-n, l-1)
                v = r1(1, 1)*r_prev(1-n, 1-l) + r1(1, -1)*r_prev(1-n, l-1) + &
                    & r1(-1, 1)*r_prev(n-1, 1-l) + r1(-1, -1)*r_prev(n-1, l-1)
                w = r1(1, 1)*r_prev(-n-1, 1-l) + r1(1, -1)*r_prev(-n-1, l-1) - &
                    & r1(-1, 1)*r_prev(n+1, 1-l) - r1(-1, -1)*r_prev(n+1, l-1)
                r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
                !dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
            end do
            dst(ind-l) = alpha * dot_product(src(ind-l:ind+l), r(-l:l, -l))
            ! Now deal with m=1-l..l-1
            do m = 1-l, l-1
                ! n = l
                v = r1(1, 0)*r_prev(l-1, m) - r1(-1, 0)*r_prev(1-l, m)
                r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
                !dst(ind+m) = src(ind+l) * r(l, m)
                ! n = -l
                v = r1(1, 0)*r_prev(1-l, m) + r1(-1, 0)*r_prev(l-1, m)
                r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
                ! n = l-1
                u = r1(0, 0) * r_prev(l-1, m)
                v = r1(1, 0)*r_prev(l-2, m) - r1(-1, 0)*r_prev(2-l, m)
                r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
                r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
                ! n = 1-l
                u = r1(0, 0) * r_prev(1-l, m)
                v = r1(1, 0)*r_prev(2-l, m) + r1(-1, 0)*r_prev(l-2, m)
                r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
                r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
                ! n = 0
                u = r1(0, 0) * r_prev(0, m)
                v = r1(1, 0)*r_prev(1, m) + r1(-1, 0)*r_prev(-1, m)
                r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
                r(0, m) = r(0, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
                ! n = 1
                u = r1(0, 0) * r_prev(1, m)
                v = r1(1, 0) * r_prev(0, m)
                w = r1(1, 0)*r_prev(2, m) + r1(-1, 0)*r_prev(-2, m)
                r(1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
                r(1, m) = r(1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
                ! n = -1
                u = r1(0, 0) * r_prev(-1, m)
                v = r1(-1, 0) * r_prev(0, m)
                w = r1(1, 0)*r_prev(-2, m) - r1(-1, 0)*r_prev(2, m)
                r(-1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
                r(-1, m) = r(-1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
                ! n = 2-l..l-2, n != -1,0,1
                do n = 2, l-2
                    u = r1(0, 0) * r_prev(n, m)
                    v = r1(1, 0)*r_prev(n-1, m) - r1(-1, 0)*r_prev(1-n, m)
                    w = r1(1, 0)*r_prev(n+1, m) + r1(-1, 0)*r_prev(-1-n, m)
                    r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                    r(n, m) = r(n, m) / scal_uvw_m(m)
                    !dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                    u = r1(0, 0) * r_prev(-n, m)
                    v = r1(1, 0)*r_prev(1-n, m) + r1(-1, 0)*r_prev(n-1, m)
                    w = r1(1, 0)*r_prev(-n-1, m) - r1(-1, 0)*r_prev(n+1, m)
                    r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                    r(-n, m) = r(-n, m) / scal_uvw_m(m)
                    !dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
                end do
                dst(ind+m) = alpha * dot_product(src(ind-l:ind+l), r(-l:l, m))
            end do
        end do
    else
        ! Compute rotations/reflections
        ! l = 0
        dst(1) = beta*dst(1) + alpha*src(1)
        if (p .eq. 0) then
            return
        end if
        ! l = 1
        dst(2) = beta*dst(2) + alpha*dot_product(r1(-1:1,-1), src(2:4))
        dst(3) = beta*dst(3) + alpha*dot_product(r1(-1:1,0), src(2:4))
        dst(4) = beta*dst(4) + alpha*dot_product(r1(-1:1,1), src(2:4))
        if (p .eq. 1) then
            return
        end if
        ! Set pointers
        n = 2*p + 1
        l = n ** 2
        r_prev(-p:p, -p:p) => work(1:l)
        m = 2*l
        r(-p:p, -p:p) => work(l+1:m)
        l = m + n
        scal_uvw_m(-p:p) => work(m+1:l)
        m = l + n
        scal_u_n(-p:p) => work(l+1:m)
        l = m + n
        scal_v_n(-p:p) => work(m+1:l)
        m = l + n
        scal_w_n(-p:p) => work(l+1:m)
        ! l = 2
        r(2, 2) = (r1(1, 1)*r1(1, 1) - r1(1, -1)*r1(1, -1) - &
            & r1(-1, 1)*r1(-1, 1) + r1(-1, -1)*r1(-1, -1)) / two
        r(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
        r(2, 0) = sqrt3 / two * (r1(1, 0)*r1(1, 0) - r1(-1, 0)*r1(-1, 0))
        r(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
        r(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
        r(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
        r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
        r(1, 0) = sqrt3 * r1(1, 0) * r1(0, 0)
        r(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
        r(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
        r(0, 2) = sqrt3 / two * (r1(0, 1)*r1(0, 1) - r1(0, -1)*r1(0, -1))
        r(0, 1) = sqrt3 * r1(0, 1) * r1(0, 0)
        r(0, 0) = (three*r1(0, 0)*r1(0, 0)-one) / two
        r(0, -1) = sqrt3 * r1(0, -1) * r1(0, 0)
        r(0, -2) = sqrt3 * r1(0, 1) * r1(0, -1)
        r(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
        r(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
        r(-1, 0) = sqrt3 * r1(-1, 0) * r1(0, 0)
        r(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(0, -1)*r1(-1, 0)
        r(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
        r(-2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
        r(-2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
        r(-2, 0) = sqrt3 * r1(1, 0) * r1(-1, 0)
        r(-2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
        r(-2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
        do m = -2, 2
            dst(7+m) = beta*dst(7+m) + alpha*dot_product(r(-2:2, m), src(5:9))
        end do
        ! l > 2
        do l = 3, p
            ! Swap previous and current rotation matrices
            r_swap => r_prev
            r_prev => r
            r => r_swap
            ! Prepare scalar factors
            scal_uvw_m(0) = dble(l)
            do m = 1, l-1
                u = sqrt(dble(l*l-m*m))
                scal_uvw_m(m) = u
                scal_uvw_m(-m) = u
            end do
            u = two * dble(l)
            u = sqrt(dble(u*(u-one)))
            scal_uvw_m(l) = u
            scal_uvw_m(-l) = u
            scal_u_n(0) = dble(l)
            scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
            scal_w_n(0) = zero
            do n = 1, l-2
                u = sqrt(dble(l*l-n*n))
                scal_u_n(n) = u
                scal_u_n(-n) = u
                v = dble(l+n)
                v = sqrt(v*(v-one)) / two
                scal_v_n(n) = v
                scal_v_n(-n) = v
                w = dble(l-n)
                w = -sqrt(w*(w-one)) / two
                scal_w_n(n) = w
                scal_w_n(-n) = w
            end do
            u = sqrt(dble(2*l-1))
            scal_u_n(l-1) = u
            scal_u_n(1-l) = u
            scal_u_n(l) = zero
            scal_u_n(-l) = zero
            v = sqrt(dble((2*l-1)*(2*l-2))) / two
            scal_v_n(l-1) = v
            scal_v_n(1-l) = v
            v = sqrt(dble(2*l*(2*l-1))) / two
            scal_v_n(l) = v
            scal_v_n(-l) = v
            scal_w_n(l-1) = zero
            scal_w_n(l) = zero
            scal_w_n(-l) = zero
            scal_w_n(1-l) = zero
            ind = l*l + l + 1
            ! m = l, n = l
            v = r1(1, 1)*r_prev(l-1, l-1) - r1(1, -1)*r_prev(l-1, 1-l) - &
                & r1(-1, 1)*r_prev(1-l, l-1) + r1(-1, -1)*r_prev(1-l, 1-l)
            r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind+l) = src(ind+l) * r(l, l)
            ! m = l, n = -l
            v = r1(1, 1)*r_prev(1-l, l-1) - r1(1, -1)*r_prev(1-l, 1-l) + &
                & r1(-1, 1)*r_prev(l-1, l-1) - r1(-1, -1)*r_prev(l-1, 1-l)
            r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
            ! m = l, n = l-1
            u = r1(0, 1)*r_prev(l-1, l-1) - r1(0, -1)*r_prev(l-1, 1-l)
            v = r1(1, 1)*r_prev(l-2, l-1) - r1(1, -1)*r_prev(l-2, 1-l) - &
                & r1(-1, 1)*r_prev(2-l, l-1) + r1(-1, -1)*r_prev(2-l, 1-l)
            r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
            ! m = l, n = 1-l
            u = r1(0, 1)*r_prev(1-l, l-1) - r1(0, -1)*r_prev(1-l, 1-l)
            v = r1(1, 1)*r_prev(2-l, l-1) - r1(1, -1)*r_prev(2-l, 1-l) + &
                & r1(-1, 1)*r_prev(l-2, l-1) - r1(-1, -1)*r_prev(l-2, 1-l)
            r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
            ! m = l, n = 1
            u = r1(0, 1)*r_prev(1, l-1) - r1(0, -1)*r_prev(1, 1-l)
            v = r1(1, 1)*r_prev(0, l-1) - r1(1, -1)*r_prev(0, 1-l)
            w = r1(1, 1)*r_prev(2, l-1) - r1(1, -1)*r_prev(2, 1-l) + &
                & r1(-1, 1)*r_prev(-2, l-1) - r1(-1, -1)*r_prev(-2, 1-l)
            r(1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, l) = r(1, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
            ! m = l, n = -1
            u = r1(0, 1)*r_prev(-1, l-1) - r1(0, -1)*r_prev(-1, 1-l)
            v = r1(-1, 1)*r_prev(0, l-1) - r1(-1, -1)*r_prev(0, 1-l)
            w = r1(1, 1)*r_prev(-2, l-1) - r1(1, -1)*r_prev(-2, 1-l) - &
                & r1(-1, 1)*r_prev(2, l-1) + r1(-1, -1)*r_prev(2, 1-l)
            r(-1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, l) = r(-1, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
            ! m = l, n = 0
            u = r1(0, 1)*r_prev(0, l-1) - r1(0, -1)*r_prev(0, 1-l)
            v = r1(1, 1)*r_prev(1, l-1) - r1(1, -1)*r_prev(1, 1-l) + &
                & r1(-1, 1)*r_prev(-1, l-1) - r1(-1, -1)*r_prev(-1, 1-l)
            r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
            ! m = l, n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 1)*r_prev(n, l-1) - r1(0, -1)*r_prev(n, 1-l)
                v = r1(1, 1)*r_prev(n-1, l-1) - r1(1, -1)*r_prev(n-1, 1-l) - &
                    & r1(-1, 1)*r_prev(1-n, l-1) + r1(-1, -1)*r_prev(1-n, 1-l)
                w = r1(1, 1)*r_prev(n+1, l-1) - r1(1, -1)*r_prev(n+1, 1-l) + &
                    & r1(-1, 1)*r_prev(-n-1, l-1) - r1(-1, -1)*r_prev(-n-1, 1-l)
                r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, l) = r(n, l) / scal_uvw_m(l)
                !dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
                u = r1(0, 1)*r_prev(-n, l-1) - r1(0, -1)*r_prev(-n, 1-l)
                v = r1(1, 1)*r_prev(1-n, l-1) - r1(1, -1)*r_prev(1-n, 1-l) + &
                    & r1(-1, 1)*r_prev(n-1, l-1) - r1(-1, -1)*r_prev(n-1, 1-l)
                w = r1(1, 1)*r_prev(-n-1, l-1) - r1(1, -1)*r_prev(-n-1, 1-l) - &
                    & r1(-1, 1)*r_prev(n+1, l-1) + r1(-1, -1)*r_prev(n+1, 1-l)
                r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, l) = r(-n, l) / scal_uvw_m(l)
                !dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
            end do
            dst(ind+l) = beta*dst(ind+l) + &
                & alpha*dot_product(src(ind-l:ind+l), r(-l:l, l))
            ! m = -l, n = l
            v = r1(1, 1)*r_prev(l-1, 1-l) + r1(1, -1)*r_prev(l-1, l-1) - &
                & r1(-1, 1)*r_prev(1-l, 1-l) - r1(-1, -1)*r_prev(1-l, l-1)
            r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind-l) = src(ind+l)*r(l, -l)
            ! m = -l, n = -l
            v = r1(1, 1)*r_prev(1-l, 1-l) + r1(1, -1)*r_prev(1-l, l-1) + &
                & r1(-1, 1)*r_prev(l-1, 1-l) + r1(-1, -1)*r_prev(l-1, l-1)
            r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
            ! m = -l, n = l-1
            u = r1(0, 1)*r_prev(l-1, 1-l) + r1(0, -1)*r_prev(l-1, l-1)
            v = r1(1, 1)*r_prev(l-2, 1-l) + r1(1, -1)*r_prev(l-2, l-1) - &
                & r1(-1, 1)*r_prev(2-l, 1-l) - r1(-1, -1)*r_prev(2-l, l-1)
            r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
            ! m = -l, n = 1-l
            u = r1(0, 1)*r_prev(1-l, 1-l) + r1(0, -1)*r_prev(1-l, l-1)
            v = r1(1, 1)*r_prev(2-l, 1-l) + r1(1, -1)*r_prev(2-l, l-1) + &
                & r1(-1, 1)*r_prev(l-2, 1-l) + r1(-1, -1)*r_prev(l-2, l-1)
            r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
            ! m = -l, n = 1
            u = r1(0, 1)*r_prev(1, 1-l) + r1(0, -1)*r_prev(1, l-1)
            v = r1(1, 1)*r_prev(0, 1-l) + r1(1, -1)*r_prev(0, l-1)
            w = r1(1, 1)*r_prev(2, 1-l) + r1(1, -1)*r_prev(2, l-1) + &
                & r1(-1, 1)*r_prev(-2, 1-l) + r1(-1, -1)*r_prev(-2, l-1)
            r(1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, -l) = r(1, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
            ! m = -l, n = -1
            u = r1(0, 1)*r_prev(-1, 1-l) + r1(0, -1)*r_prev(-1, l-1)
            v = r1(-1, 1)*r_prev(0, 1-l) + r1(-1, -1)*r_prev(0, l-1)
            w = r1(1, 1)*r_prev(-2, 1-l) + r1(1, -1)*r_prev(-2, l-1) - &
                & r1(-1, 1)*r_prev(2, 1-l) - r1(-1, -1)*r_prev(2, l-1)
            r(-1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
            ! m = -l, n = 0
            u = r1(0, 1)*r_prev(0, 1-l) + r1(0, -1)*r_prev(0, l-1)
            v = r1(1, 1)*r_prev(1, 1-l) + r1(1, -1)*r_prev(1, l-1) + &
                & r1(-1, 1)*r_prev(-1, 1-l) + r1(-1, -1)*r_prev(-1, l-1)
            r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
            ! m = -l, n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 1)*r_prev(n, 1-l) + r1(0, -1)*r_prev(n, l-1)
                v = r1(1, 1)*r_prev(n-1, 1-l) + r1(1, -1)*r_prev(n-1, l-1) - &
                    & r1(-1, 1)*r_prev(1-n, 1-l) - r1(-1, -1)*r_prev(1-n, l-1)
                w = r1(1, 1)*r_prev(n+1, 1-l) + r1(1, -1)*r_prev(n+1, l-1) + &
                    & r1(-1, 1)*r_prev(-n-1, 1-l) + r1(-1, -1)*r_prev(-n-1, l-1)
                r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, -l) = r(n, -l) / scal_uvw_m(l)
                !dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
                u = r1(0, 1)*r_prev(-n, 1-l) + r1(0, -1)*r_prev(-n, l-1)
                v = r1(1, 1)*r_prev(1-n, 1-l) + r1(1, -1)*r_prev(1-n, l-1) + &
                    & r1(-1, 1)*r_prev(n-1, 1-l) + r1(-1, -1)*r_prev(n-1, l-1)
                w = r1(1, 1)*r_prev(-n-1, 1-l) + r1(1, -1)*r_prev(-n-1, l-1) - &
                    & r1(-1, 1)*r_prev(n+1, 1-l) - r1(-1, -1)*r_prev(n+1, l-1)
                r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
                !dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
            end do
            dst(ind-l) = beta*dst(ind-l) + &
                & alpha*dot_product(src(ind-l:ind+l), r(-l:l, -l))
            ! Now deal with m=1-l..l-1
            do m = 1-l, l-1
                ! n = l
                v = r1(1, 0)*r_prev(l-1, m) - r1(-1, 0)*r_prev(1-l, m)
                r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
                !dst(ind+m) = src(ind+l) * r(l, m)
                ! n = -l
                v = r1(1, 0)*r_prev(1-l, m) + r1(-1, 0)*r_prev(l-1, m)
                r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
                ! n = l-1
                u = r1(0, 0) * r_prev(l-1, m)
                v = r1(1, 0)*r_prev(l-2, m) - r1(-1, 0)*r_prev(2-l, m)
                r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
                r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
                ! n = 1-l
                u = r1(0, 0) * r_prev(1-l, m)
                v = r1(1, 0)*r_prev(2-l, m) + r1(-1, 0)*r_prev(l-2, m)
                r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
                r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
                ! n = 0
                u = r1(0, 0) * r_prev(0, m)
                v = r1(1, 0)*r_prev(1, m) + r1(-1, 0)*r_prev(-1, m)
                r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
                r(0, m) = r(0, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
                ! n = 1
                u = r1(0, 0) * r_prev(1, m)
                v = r1(1, 0) * r_prev(0, m)
                w = r1(1, 0)*r_prev(2, m) + r1(-1, 0)*r_prev(-2, m)
                r(1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
                r(1, m) = r(1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
                ! n = -1
                u = r1(0, 0) * r_prev(-1, m)
                v = r1(-1, 0) * r_prev(0, m)
                w = r1(1, 0)*r_prev(-2, m) - r1(-1, 0)*r_prev(2, m)
                r(-1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
                r(-1, m) = r(-1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
                ! n = 2-l..l-2, n != -1,0,1
                do n = 2, l-2
                    u = r1(0, 0) * r_prev(n, m)
                    v = r1(1, 0)*r_prev(n-1, m) - r1(-1, 0)*r_prev(1-n, m)
                    w = r1(1, 0)*r_prev(n+1, m) + r1(-1, 0)*r_prev(-1-n, m)
                    r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                    r(n, m) = r(n, m) / scal_uvw_m(m)
                    !dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                    u = r1(0, 0) * r_prev(-n, m)
                    v = r1(1, 0)*r_prev(1-n, m) + r1(-1, 0)*r_prev(n-1, m)
                    w = r1(1, 0)*r_prev(-n-1, m) - r1(-1, 0)*r_prev(n+1, m)
                    r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                    r(-n, m) = r(-n, m) / scal_uvw_m(m)
                    !dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
                end do
                dst(ind+m) = beta*dst(ind+m) + &
                    & alpha*dot_product(src(ind-l:ind+l), r(-l:l, m))
            end do
        end do
    end if
end subroutine fmm_sph_transform_work

!> Save matrix of transformation of spherical harmonics
!!
!! This procedure is based on @ref fmm_sph_transform_out but stores
!! corresponding transformation matrix in a data-efficient sparse way. For a
!! full description of the input parameter `r1` please take a look at @ref
!! fmm_sph_transform.
!!
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] r1: Transformation from new to old cartesian coordinates
!! @param[out] mat: Transformation matrix in a sparse storage
!!
!! @sa fmm_sph_transform, fmm_sph_transform_use_mat
subroutine fmm_sph_transform_get_mat(p, r1, mat)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1(-1:1, -1:1)
    ! Output
    real(dp), intent(out) :: mat((p+1)*(2*p+1)*(2*p+3)/3)
    ! Local variables
    real(dp) :: u, v, w, r_prev(-p:p, -p:p), r(-p:p, -p:p)
    real(dp) :: scal_uvw_m(-p:p), scal_u_n(-p:p), scal_v_n(-p:p)
    real(dp) :: scal_w_n(-p:p)
    integer :: l, m, n, ind
    ! Compute rotations/reflections
    ! l = 0
    !dst(1) = src(1)
    mat(1) = 1
    if (p .eq. 0) then
        return
    end if
    ! l = 1
    do m = -1, 1
        !dst(3+m) = 0
        ind = 3*m + 6
        do n = -1, 1
            !dst(3+m) = dst(3+m) + r1(n, m)*src(3+n)
            mat(ind+n) = r1(n, m)
        end do
    end do
    if (p .eq. 1) then
        return
    end if
    ! l = 2
    r(2, 2) = (r1(1, 1)*r1(1, 1) - r1(1, -1)*r1(1, -1) - &
        & r1(-1, 1)*r1(-1, 1) + r1(-1, -1)*r1(-1, -1)) / 2
    r(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
    r(2, 0) = sqrt(dble(3)) / 2 * (r1(1, 0)*r1(1, 0) - r1(-1, 0)*r1(-1, 0))
    r(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
    r(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
    r(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
    r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
    r(1, 0) = sqrt(dble(3)) * r1(1, 0) * r1(0, 0)
    r(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
    r(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
    r(0, 2) = sqrt(dble(3)) / 2 * (r1(0, 1)*r1(0, 1) - r1(0, -1)*r1(0, -1))
    r(0, 1) = sqrt(dble(3)) * r1(0, 1) * r1(0, 0)
    r(0, 0) = (3*r1(0, 0)*r1(0, 0)-1) / 2
    r(0, -1) = sqrt(dble(3)) * r1(0, -1) * r1(0, 0)
    r(0, -2) = sqrt(dble(3)) * r1(0, 1) * r1(0, -1)
    r(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
    r(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
    r(-1, 0) = sqrt(dble(3)) * r1(-1, 0) * r1(0, 0)
    r(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(0, -1)*r1(-1, 0)
    r(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
    r(-2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
    r(-2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
    r(-2, 0) = sqrt(dble(3)) * r1(1, 0) * r1(-1, 0)
    r(-2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
    r(-2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
    do m = -2, 2
        !dst(7+m) = 0
        ind = 5*m + 23
        do n = -2, 2
            !dst(7+m) = dst(7+m) + src(7+n)*r(n, m)
            mat(ind+n) = r(n, m)
        end do
    end do
    ! l > 2
    do l = 3, p
        ! Prepare previous rotation matrix
        r_prev(1-l:l-1, 1-l:l-1) = r(1-l:l-1, 1-l:l-1)
        ! Prepare scalar factors
        scal_uvw_m(0) = l
        do m = 1, l-1
        u = sqrt(dble(l*l-m*m))
        scal_uvw_m(m) = u
        scal_uvw_m(-m) = u
        end do
        u = 2 * l
        u = sqrt(dble(u*(u-1)))
        scal_uvw_m(l) = u
        scal_uvw_m(-l) = u
        scal_u_n(0) = l
        scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
        scal_w_n(0) = 0
        do n = 1, l-2
        u = sqrt(dble(l*l-n*n))
        scal_u_n(n) = u
        scal_u_n(-n) = u
        v = l + n
        v = sqrt(v*(v-1)) / 2
        scal_v_n(n) = v
        scal_v_n(-n) = v
        w = l - n
        w = -sqrt(w*(w-1)) / 2
        scal_w_n(n) = w
        scal_w_n(-n) = w
        end do
        u = sqrt(dble(2*l-1))
        scal_u_n(l-1) = u
        scal_u_n(1-l) = u
        scal_u_n(l) = 0
        scal_u_n(-l) = 0
        v = sqrt(dble((2*l-1)*(2*l-2))) / 2
        scal_v_n(l-1) = v
        scal_v_n(1-l) = v
        v = sqrt(dble(2*l*(2*l-1))) / 2
        scal_v_n(l) = v
        scal_v_n(-l) = v
        scal_w_n(l-1) = 0
        scal_w_n(l) = 0
        scal_w_n(-l) = 0
        scal_w_n(1-l) = 0
        ind = l*l + l + 1
        ! m = l, n = l
        v = r1(1, 1)*r_prev(l-1, l-1) - r1(1, -1)*r_prev(l-1, 1-l) - &
            & r1(-1, 1)*r_prev(1-l, l-1) + r1(-1, -1)*r_prev(1-l, 1-l)
        r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        !dst(ind+l) = src(ind+l) * r(l, l)
        ! m = l, n = -l
        v = r1(1, 1)*r_prev(1-l, l-1) - r1(1, -1)*r_prev(1-l, 1-l) + &
            & r1(-1, 1)*r_prev(l-1, l-1) - r1(-1, -1)*r_prev(l-1, 1-l)
        r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
        ! m = l, n = l-1
        u = r1(0, 1)*r_prev(l-1, l-1) - r1(0, -1)*r_prev(l-1, 1-l)
        v = r1(1, 1)*r_prev(l-2, l-1) - r1(1, -1)*r_prev(l-2, 1-l) - &
            & r1(-1, 1)*r_prev(2-l, l-1) + r1(-1, -1)*r_prev(2-l, 1-l)
        r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
        ! m = l, n = 1-l
        u = r1(0, 1)*r_prev(1-l, l-1) - r1(0, -1)*r_prev(1-l, 1-l)
        v = r1(1, 1)*r_prev(2-l, l-1) - r1(1, -1)*r_prev(2-l, 1-l) + &
            & r1(-1, 1)*r_prev(l-2, l-1) - r1(-1, -1)*r_prev(l-2, 1-l)
        r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
        ! m = l, n = 1
        u = r1(0, 1)*r_prev(1, l-1) - r1(0, -1)*r_prev(1, 1-l)
        v = r1(1, 1)*r_prev(0, l-1) - r1(1, -1)*r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(2, l-1) - r1(1, -1)*r_prev(2, 1-l) + &
            & r1(-1, 1)*r_prev(-2, l-1) - r1(-1, -1)*r_prev(-2, 1-l)
        r(1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, l) = r(1, l) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
        ! m = l, n = -1
        u = r1(0, 1)*r_prev(-1, l-1) - r1(0, -1)*r_prev(-1, 1-l)
        v = r1(-1, 1)*r_prev(0, l-1) - r1(-1, -1)*r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(-2, l-1) - r1(1, -1)*r_prev(-2, 1-l) - &
            & r1(-1, 1)*r_prev(2, l-1) + r1(-1, -1)*r_prev(2, 1-l)
        r(-1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, l) = r(-1, l) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
        ! m = l, n = 0
        u = r1(0, 1)*r_prev(0, l-1) - r1(0, -1)*r_prev(0, 1-l)
        v = r1(1, 1)*r_prev(1, l-1) - r1(1, -1)*r_prev(1, 1-l) + &
            & r1(-1, 1)*r_prev(-1, l-1) - r1(-1, -1)*r_prev(-1, 1-l)
        r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
        ! m = l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1)*r_prev(n, l-1) - r1(0, -1)*r_prev(n, 1-l)
            v = r1(1, 1)*r_prev(n-1, l-1) - r1(1, -1)*r_prev(n-1, 1-l) - &
                & r1(-1, 1)*r_prev(1-n, l-1) + r1(-1, -1)*r_prev(1-n, 1-l)
            w = r1(1, 1)*r_prev(n+1, l-1) - r1(1, -1)*r_prev(n+1, 1-l) + &
                & r1(-1, 1)*r_prev(-n-1, l-1) - r1(-1, -1)*r_prev(-n-1, 1-l)
            r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, l) = r(n, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
            u = r1(0, 1)*r_prev(-n, l-1) - r1(0, -1)*r_prev(-n, 1-l)
            v = r1(1, 1)*r_prev(1-n, l-1) - r1(1, -1)*r_prev(1-n, 1-l) + &
                & r1(-1, 1)*r_prev(n-1, l-1) - r1(-1, -1)*r_prev(n-1, 1-l)
            w = r1(1, 1)*r_prev(-n-1, l-1) - r1(1, -1)*r_prev(-n-1, 1-l) - &
                & r1(-1, 1)*r_prev(n+1, l-1) + r1(-1, -1)*r_prev(n+1, 1-l)
            r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, l) = r(-n, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
        end do
        ! m = -l, n = l
        v = r1(1, 1)*r_prev(l-1, 1-l) + r1(1, -1)*r_prev(l-1, l-1) - &
            & r1(-1, 1)*r_prev(1-l, 1-l) - r1(-1, -1)*r_prev(1-l, l-1)
        r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        !dst(ind-l) = src(ind+l)*r(l, -l)
        ! m = -l, n = -l
        v = r1(1, 1)*r_prev(1-l, 1-l) + r1(1, -1)*r_prev(1-l, l-1) + &
            & r1(-1, 1)*r_prev(l-1, 1-l) + r1(-1, -1)*r_prev(l-1, l-1)
        r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
        ! m = -l, n = l-1
        u = r1(0, 1)*r_prev(l-1, 1-l) + r1(0, -1)*r_prev(l-1, l-1)
        v = r1(1, 1)*r_prev(l-2, 1-l) + r1(1, -1)*r_prev(l-2, l-1) - &
            & r1(-1, 1)*r_prev(2-l, 1-l) - r1(-1, -1)*r_prev(2-l, l-1)
        r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
        ! m = -l, n = 1-l
        u = r1(0, 1)*r_prev(1-l, 1-l) + r1(0, -1)*r_prev(1-l, l-1)
        v = r1(1, 1)*r_prev(2-l, 1-l) + r1(1, -1)*r_prev(2-l, l-1) + &
            & r1(-1, 1)*r_prev(l-2, 1-l) + r1(-1, -1)*r_prev(l-2, l-1)
        r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
        ! m = -l, n = 1
        u = r1(0, 1)*r_prev(1, 1-l) + r1(0, -1)*r_prev(1, l-1)
        v = r1(1, 1)*r_prev(0, 1-l) + r1(1, -1)*r_prev(0, l-1)
        w = r1(1, 1)*r_prev(2, 1-l) + r1(1, -1)*r_prev(2, l-1) + &
            & r1(-1, 1)*r_prev(-2, 1-l) + r1(-1, -1)*r_prev(-2, l-1)
        r(1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, -l) = r(1, -l) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
        ! m = -l, n = -1
        u = r1(0, 1)*r_prev(-1, 1-l) + r1(0, -1)*r_prev(-1, l-1)
        v = r1(-1, 1)*r_prev(0, 1-l) + r1(-1, -1)*r_prev(0, l-1)
        w = r1(1, 1)*r_prev(-2, 1-l) + r1(1, -1)*r_prev(-2, l-1) - &
            & r1(-1, 1)*r_prev(2, 1-l) - r1(-1, -1)*r_prev(2, l-1)
        r(-1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
        ! m = -l, n = 0
        u = r1(0, 1)*r_prev(0, 1-l) + r1(0, -1)*r_prev(0, l-1)
        v = r1(1, 1)*r_prev(1, 1-l) + r1(1, -1)*r_prev(1, l-1) + &
            & r1(-1, 1)*r_prev(-1, 1-l) + r1(-1, -1)*r_prev(-1, l-1)
        r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
        ! m = -l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1)*r_prev(n, 1-l) + r1(0, -1)*r_prev(n, l-1)
            v = r1(1, 1)*r_prev(n-1, 1-l) + r1(1, -1)*r_prev(n-1, l-1) - &
                & r1(-1, 1)*r_prev(1-n, 1-l) - r1(-1, -1)*r_prev(1-n, l-1)
            w = r1(1, 1)*r_prev(n+1, 1-l) + r1(1, -1)*r_prev(n+1, l-1) + &
                & r1(-1, 1)*r_prev(-n-1, 1-l) + r1(-1, -1)*r_prev(-n-1, l-1)
            r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, -l) = r(n, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
            u = r1(0, 1)*r_prev(-n, 1-l) + r1(0, -1)*r_prev(-n, l-1)
            v = r1(1, 1)*r_prev(1-n, 1-l) + r1(1, -1)*r_prev(1-n, l-1) + &
                & r1(-1, 1)*r_prev(n-1, 1-l) + r1(-1, -1)*r_prev(n-1, l-1)
            w = r1(1, 1)*r_prev(-n-1, 1-l) + r1(1, -1)*r_prev(-n-1, l-1) - &
                & r1(-1, 1)*r_prev(n+1, 1-l) - r1(-1, -1)*r_prev(n+1, l-1)
            r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
        end do
        ! Now deal with m=1-l..l-1
        do m = 1-l, l-1
            ! n = l
            v = r1(1, 0)*r_prev(l-1, m) - r1(-1, 0)*r_prev(1-l, m)
            r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            !dst(ind+m) = src(ind+l) * r(l, m)
            ! n = -l
            v = r1(1, 0)*r_prev(1-l, m) + r1(-1, 0)*r_prev(l-1, m)
            r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
            ! n = l-1
            u = r1(0, 0) * r_prev(l-1, m)
            v = r1(1, 0)*r_prev(l-2, m) - r1(-1, 0)*r_prev(2-l, m)
            r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
            ! n = 1-l
            u = r1(0, 0) * r_prev(1-l, m)
            v = r1(1, 0)*r_prev(2-l, m) + r1(-1, 0)*r_prev(l-2, m)
            r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
            ! n = 0
            u = r1(0, 0) * r_prev(0, m)
            v = r1(1, 0)*r_prev(1, m) + r1(-1, 0)*r_prev(-1, m)
            r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
            r(0, m) = r(0, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
            ! n = 1
            u = r1(0, 0) * r_prev(1, m)
            v = r1(1, 0) * r_prev(0, m)
            w = r1(1, 0)*r_prev(2, m) + r1(-1, 0)*r_prev(-2, m)
            r(1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, m) = r(1, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
            ! n = -1
            u = r1(0, 0) * r_prev(-1, m)
            v = r1(-1, 0) * r_prev(0, m)
            w = r1(1, 0)*r_prev(-2, m) - r1(-1, 0)*r_prev(2, m)
            r(-1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, m) = r(-1, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
            ! n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 0) * r_prev(n, m)
                v = r1(1, 0)*r_prev(n-1, m) - r1(-1, 0)*r_prev(1-n, m)
                w = r1(1, 0)*r_prev(n+1, m) + r1(-1, 0)*r_prev(-1-n, m)
                r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, m) = r(n, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                u = r1(0, 0) * r_prev(-n, m)
                v = r1(1, 0)*r_prev(1-n, m) + r1(-1, 0)*r_prev(n-1, m)
                w = r1(1, 0)*r_prev(-n-1, m) - r1(-1, 0)*r_prev(n+1, m)
                r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, m) = r(-n, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
            end do
        end do
        do m = -l, l
            ind = 2*l*(l+1)*(2*l+1)/3 + l + 1
            do n = -l, l
                mat(ind+m*(2*l+1)+n) = r(n, m)
            end do
        end do
    end do
end subroutine fmm_sph_transform_get_mat

!> Apply matrix of transformation of spherical harmonics
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha \mathrm{mat} \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! \mathrm{mat} \f$ is a matrix of transformation of spherical harmonics
!! precomputed by @ref fmm_sph_transform_get_mat.
!!
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] mat: Transformation matrix in a sparse storage
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multiplier for `dst`
!! @param[inout] dst: Coefficients of transformed spherical harmonics
!!
!! @sa fmm_sph_transform, fmm_sph_transform_get_mat
subroutine fmm_sph_transform_use_mat(p, mat, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(2*p+1)*(2*p+3)/3), alpha, &
        & src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Local variables
    integer :: l, m, n, ind, indl
    ! Init dst
    if (beta .eq. zero) then
        dst = zero
    else
        dst = beta * dst
    end if
    ! l = 0
    dst(1) = dst(1) + alpha*src(1)
    if(p .eq. 0) then
        return
    end if
    do l = 1, p
        ! magical value for the offset to the current reflection matrix
        ind = 2*l*(l+1)*(2*l+1)/3 + l + 1
        ! offset for current spherical harmonics
        indl = l*l + l + 1
        do m = -l, l
            do n = -l, l
                dst(indl+m) = dst(indl+m) + &
                    & alpha*mat(ind+m*(2*l+1)+n)*src(indl+n)
            end do
        end do
    end do
end subroutine fmm_sph_transform_use_mat

!> Rotate spherical harmonics around OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Call corresponding work routine
    call fmm_sph_rotate_oz_work(p, vcos, vsin, alpha, src, beta, dst)
end subroutine fmm_sph_rotate_oz

!> Rotate spherical harmonics around OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz_work(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Local variables
    integer :: l, m, ind
    real(dp) :: v1, v2, v3, v4
    ! In case alpha is zero just scale output
    if (alpha .eq. zero) then
        ! Set output to zero if beta is also zero
        if (beta .eq. zero) then
            dst = zero
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-zero
    ! In case beta is zero output is just overwritten without being read
    if (beta .eq. zero) then
        ! l = 0
        dst(1) = alpha*src(1)
        ! l > 0
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = alpha*src(ind)
            ! m != 0
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = alpha * (v1*v3-v2*v4)
                ! m < 0
                dst(ind-m) = alpha * (v1*v4+v2*v3)
            end do
        end do
    else
        ! l = 0
        dst(1) = beta*dst(1) + alpha*src(1)
        ! l > 0
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = beta*dst(ind) + alpha*src(ind)
            ! m != 0
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = beta*dst(ind+m) + alpha*(v1*v3-v2*v4)
                ! m < 0
                dst(ind-m) = beta*dst(ind-m) + alpha*(v1*v4+v2*v3)
            end do
        end do
    end if
end subroutine fmm_sph_rotate_oz_work

!> Rotate spherical harmonics around OZ axis in an opposite direction
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz_adj(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Call corresponding work routine
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src, beta, dst)
end subroutine fmm_sph_rotate_oz_adj

!> Rotate spherical harmonics around OZ axis in an opposite direction
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Local variables
    integer :: l, m, ind
    real(dp) :: v1, v2, v3, v4
    ! In case alpha is zero just scale output
    if (alpha .eq. zero) then
        ! Set output to zero if beta is also zero
        if (beta .eq. zero) then
            dst = zero
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-zero
    ! In case beta is zero output is just overwritten without being read
    if (beta .eq. zero) then
        ! l = 0
        dst(1) = alpha*src(1)
        ! l > 0
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = alpha*src(ind)
            ! m != 0
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = alpha * (v1*v3+v2*v4)
                ! m < 0
                dst(ind-m) = alpha * (v2*v3-v1*v4)
            end do
        end do
    else
        ! l = 0
        dst(1) = beta*dst(1) + alpha*src(1)
        ! l > 0
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = beta*dst(ind) + alpha*src(ind)
            ! m != 0
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = beta*dst(ind+m) + alpha*(v1*v3+v2*v4)
                ! m < 0
                dst(ind-m) = beta*dst(ind-m) + alpha*(v2*v3-v1*v4)
            end do
        end do
    end if
end subroutine fmm_sph_rotate_oz_adj_work

!> Transform spherical harmonics in the OXZ plane
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of transformation of coordinates in the OXZ plane (y
!! coordinate remains the same) presented by 2-by-2 matrix `r1xz`.
!!
!! Based on @ref fmm_sph_transform
!! by assuming `r1(-1, 0) = r1(-1, 1) = r1(0, -1) = r1(1, -1) = 0` and
!! `r1(-1, -1) = 1`, which corresponds to the following transformation matrix:
!! \f[
!!      R_1 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & a & b \\ 0 & c & d
!!      \end{bmatrix},
!! \f]
!! where unkown elements represent input `r1xz` 2x2 array:
!! \f[
!!      R_1^{xz} = \begin{bmatrix} a & b \\ c & d \end{bmatrix}.
!! \f]
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] r1xz: Transformation from new to old coordinates in the OXZ plane
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of transformed spherical harmonics
subroutine fmm_sph_rotate_oxz(p, ctheta, stheta, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: ctheta, stheta, alpha, src((p+1)**2), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)**2)
    ! Temporary workspace
    real(dp) :: work(2*(2*p+1)*(2*p+3)+p)
    ! Call corresponding work routine
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, alpha, src, beta, dst, &
        & work)
end subroutine fmm_sph_rotate_oxz

!> Transform spherical harmonics in the OXZ plane
!!
!! This function implements @ref fmm_sph_rotate_oxz with predefined values
!! of parameters \p alpha=one and \p beta=zero.
!! 
!! @param[in] p: maximum order of spherical harmonics
!! @param[in] r1xz: 2D transformation matrix in the OXZ plane
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[out] dst: coefficients of rotated spherical harmonics
!! @param[out] work: Temporary workspace of a size (2*(2*p+1)*(2*p+3))
subroutine fmm_sph_rotate_oxz_work(p, ctheta, stheta, alpha, src, beta, dst, &
        & work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: ctheta, stheta, alpha, src((p+1)**2), beta
    ! Output
    real(dp), intent(out) :: dst((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(4*p*p+13*p+4)
    ! Local variables
    real(dp) :: u, v, w, fl, fl2, tmp1, tmp2, vu(2), vv(2), vw(2), &
        ctheta2, stheta2, cstheta
    integer :: l, m, n, ind, k
    ! Pointers for a workspace
    real(dp), pointer :: r(:, :, :), r_prev(:, :, :), scal_uvw_m(:), &
        & scal_u_n(:), scal_v_n(:), scal_w_n(:), r_swap(:, :, :), vsqr(:)
    !! Spherical harmonics Y_l^m with negative m transform into harmonics Y_l^m
    !! with the same l and negative m, while harmonics Y_l^m with non-negative
    !! m transform into harmonics Y_l^m with the same l and non-negative m.
    !! Transformations for non-negative m will be stored in r(1, :, :) and for
    !! negative m will be in r(2, :, :)
    ! In case alpha is zero just scale output
    if (alpha .eq. zero) then
        ! Set output to zero if beta is also zero
        if (beta .eq. zero) then
            dst = zero
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-zero
    ! In case beta is zero output is just overwritten without being read
    if (beta .eq. zero) then
        ! Compute rotations/reflections
        ! l = 0
        dst(1) = alpha * src(1)
        if (p .eq. 0) then
            return
        end if
        ! l = 1
        dst(2) = alpha * src(2)
        dst(3) = alpha * (src(3)*ctheta - src(4)*stheta)
        dst(4) = alpha * (src(3)*stheta + src(4)*ctheta)
        if (p .eq. 1) then
            return
        end if
        ! Set pointers
        l = 2 * (p+1) * (p+1)
        r(1:2, 0:p, 0:p) => work(1:l)
        m = 2 * l
        r_prev(1:2, 0:p, 0:p) => work(l+1:m)
        l = m + p + 1
        scal_uvw_m(0:p) => work(m+1:l)
        m = l + p
        scal_u_n(0:p-1) => work(l+1:m)
        l = m + p + 1
        scal_v_n(0:p) => work(m+1:l)
        m = l + p - 2
        scal_w_n(1:p-2) => work(l+1:m)
        l = m + p
        vsqr(1:p) => work(m+1:l)
        ! l = 2, m >= 0
        ctheta2 = ctheta * ctheta
        cstheta = ctheta * stheta
        stheta2 = stheta * stheta
        r(1, 2, 2) = (ctheta2 + one) / two
        r(1, 1, 2) = cstheta
        r(1, 0, 2) = sqrt3 / two * stheta2
        dst(9) = alpha * (src(9)*r(1, 2, 2) + src(8)*r(1, 1, 2) + &
            & src(7)*r(1, 0, 2))
        r(1, 2, 1) = -cstheta
        r(1, 1, 1) = ctheta2 - stheta2
        r(1, 0, 1) = sqrt3 * cstheta
        dst(8) = alpha * (src(9)*r(1, 2, 1) + src(8)*r(1, 1, 1) + &
            & src(7)*r(1, 0, 1))
        r(1, 2, 0) = sqrt3 / two * stheta2
        r(1, 1, 0) = -sqrt3 * cstheta
        r(1, 0, 0) = (three*ctheta2-one) / two
        dst(7) = alpha * (src(9)*r(1, 2, 0) + src(8)*r(1, 1, 0) + &
            & src(7)*r(1, 0, 0))
        ! l = 2,  m < 0
        r(2, 1, 1) = ctheta
        r(2, 2, 1) = -stheta
        dst(6) = alpha * (src(6)*r(2, 1, 1) + src(5)*r(2, 2, 1))
        r(2, 1, 2) = stheta
        r(2, 2, 2) = ctheta
        dst(5) = alpha * (src(6)*r(2, 1, 2) + src(5)*r(2, 2, 2))
        ! l > 2
        vsqr(1) = one
        vsqr(2) = four
        do l = 3, p
            ! Swap previous and current rotation matrices
            r_swap => r_prev
            r_prev => r
            r => r_swap
            ! Prepare scalar factors
            fl = dble(l)
            fl2 = fl * fl
            vsqr(l) = fl2
            scal_uvw_m(0) = fl
            do m = 1, l-1
                u = sqrt(fl2 - vsqr(m))
                scal_uvw_m(m) = u
            end do
            u = two * dble(l)
            u = sqrt(dble(u*(u-one)))
            scal_uvw_m(l) = u
            scal_u_n(0) = dble(l)
            scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
            do n = 1, l-2
                u = sqrt(fl2-vsqr(n))
                scal_u_n(n) = u
            end do
            do n = 1, l-2
                v = dble(l+n)
                v = sqrt(v*v-v) / two
                scal_v_n(n) = v
                w = dble(l-n)
                w = -sqrt(w*w-w) / two
                scal_w_n(n) = w
            end do
            u = sqrt(dble(2*l-1))
            scal_u_n(l-1) = u
            v = sqrt(dble((2*l-1)*(2*l-2))) / two
            scal_v_n(l-1) = v
            v = sqrt(dble(2*l*(2*l-1))) / two
            scal_v_n(l) = v
            ind = l*l + l + 1
            ! m = l, n = l and m = -l, n = - l
            vv = ctheta*r_prev(:, l-1, l-1) + r_prev(2:1:-1, l-1, l-1)
            r(:, l, l) = vv * scal_v_n(l) / scal_uvw_m(l)
            tmp1 = src(ind+l) * r(1, l, l)
            tmp2 = src(ind-l) * r(2, l, l)
            ! m = l, n = l-1 and m = -l, n = 1-l
            vu = stheta * r_prev(:, l-1, l-1)
            vv = ctheta*r_prev(:, l-2, l-1) + r_prev(2:1:-1, l-2, l-1)
            r(:, l-1, l) = (vu*scal_u_n(l-1)+vv*scal_v_n(l-1)) / scal_uvw_m(l)
            tmp1 = tmp1 + src(ind+l-1)*r(1, l-1, l)
            tmp2 = tmp2 + src(ind-l+1)*r(2, l-1, l)
            ! m = l, n = 1 and m = -l, n = -1
            vu = stheta * r_prev(:, 1, l-1)
            vv(1) = ctheta * r_prev(1, 0, l-1)
            vv(2) = r_prev(1, 0, l-1)
            vw = ctheta*r_prev(:, 2, l-1) - r_prev(2:1:-1, 2, l-1)
            r(:, 1, l) = vu*scal_u_n(1) + vw*scal_w_n(1) + sqrt2*scal_v_n(1)*vv
            r(:, 1, l) = r(:, 1, l) / scal_uvw_m(l)
            tmp1 = tmp1 + src(ind+1)*r(1, 1, l)
            tmp2 = tmp2 + src(ind-1)*r(2, 1, l)
            ! m = l, n = 0
            u = stheta * r_prev(1, 0, l-1)
            v = ctheta*r_prev(1, 1, l-1) - r_prev(2, 1, l-1)
            r(1, 0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
            tmp1 = tmp1 + src(ind)*r(1, 0, l)
            ! m = l, n = 2..l-2 and m = -l, n = 2-l..-2
            do n = 2, l-2
                vu = stheta * r_prev(:, n, l-1)
                vv = ctheta*r_prev(:, n-1, l-1) + r_prev(2:1:-1, n-1, l-1)
                vw = ctheta*r_prev(:, n+1, l-1) - r_prev(2:1:-1, n+1, l-1)
                vu = vu*scal_u_n(n) + vv*scal_v_n(n) + vw*scal_w_n(n)
                r(:, n, l) = vu / scal_uvw_m(l)
                tmp1 = tmp1 + src(ind+n)*r(1, n, l)
                tmp2 = tmp2 + src(ind-n)*r(2, n, l)
            end do
            dst(ind+l) = alpha * tmp1
            dst(ind-l) = alpha * tmp2
            ! Now deal with m = 0
            ! n = l and n = -l
            v = -stheta * r_prev(1, l-1, 0)
            u = scal_v_n(l) / scal_uvw_m(0)
            r(1, l, 0) = v * u
            tmp1 = src(ind+l) * r(1, l, 0)
            ! n = l-1
            u = ctheta * r_prev(1, l-1, 0)
            v = -stheta * r_prev(1, l-2, 0)
            w = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(1, l-1, 0) = w / scal_uvw_m(0)
            tmp1 = tmp1 + src(ind+l-1)*r(1, l-1, 0)
            ! n = 0
            u = ctheta * r_prev(1, 0, 0)
            v = -stheta * r_prev(1, 1, 0)
            w = u*scal_u_n(0) + v*scal_v_n(0)
            r(1, 0, 0) = w / scal_uvw_m(0)
            tmp1 = tmp1 + src(ind)*r(1, 0, 0)
            ! n = 1
            v = sqrt2*scal_v_n(1)*r_prev(1, 0, 0) + &
                & scal_w_n(1)*r_prev(1, 2, 0)
            u = ctheta * r_prev(1, 1, 0)
            w = scal_u_n(1)*u - stheta*v
            r(1, 1, 0) = w / scal_uvw_m(0)
            tmp1 = tmp1 + src(ind+1)*r(1, 1, 0)
            ! n = 2..l-2
            do n = 2, l-2
                v = scal_v_n(n)*r_prev(1, n-1, 0) + &
                    & scal_w_n(n)*r_prev(1, n+1, 0)
                u = ctheta * r_prev(1, n, 0)
                w = scal_u_n(n)*u - stheta*v
                r(1, n, 0) = w / scal_uvw_m(0)
                tmp1 = tmp1 + src(ind+n)*r(1, n, 0)
            end do
            dst(ind) = alpha * tmp1
            ! Now deal with m=1..l-1 and m=1-l..-1
            do m = 1, l-1
                ! n = l and n = -l
                vv = -stheta * r_prev(:, l-1, m)
                u = scal_v_n(l) / scal_uvw_m(m)
                r(:, l, m) = vv * u
                tmp1 = src(ind+l) * r(1, l, m)
                tmp2 = src(ind-l) * r(2, l, m)
                ! n = l-1 and n = 1-l
                vu = ctheta * r_prev(:, l-1, m)
                vv = -stheta * r_prev(:, l-2, m)
                vw = vu*scal_u_n(l-1) + vv*scal_v_n(l-1)
                r(:, l-1, m) = vw / scal_uvw_m(m)
                tmp1 = tmp1 + src(ind+l-1)*r(1, l-1, m)
                tmp2 = tmp2 + src(ind-l+1)*r(2, l-1, m)
                ! n = 0
                u = ctheta * r_prev(1, 0, m)
                v = -stheta * r_prev(1, 1, m)
                w = u*scal_u_n(0) + v*scal_v_n(0)
                r(1, 0, m) = w / scal_uvw_m(m)
                tmp1 = tmp1 + src(ind)*r(1, 0, m)
                ! n = 1
                v = sqrt2*scal_v_n(1)*r_prev(1, 0, m) + &
                    & scal_w_n(1)*r_prev(1, 2, m)
                u = ctheta * r_prev(1, 1, m)
                w = scal_u_n(1)*u - stheta*v
                r(1, 1, m) = w / scal_uvw_m(m)
                tmp1 = tmp1 + src(ind+1)*r(1, 1, m)
                ! n = -1
                u = ctheta * r_prev(2, 1, m)
                w = -stheta * r_prev(2, 2, m)
                v = u*scal_u_n(1) + w*scal_w_n(1)
                r(2, 1, m) = v / scal_uvw_m(m)
                tmp2 = tmp2 + src(ind-1)*r(2, 1, m)
                ! n = 2..l-2 and n = 2-l..-2
                do n = 2, l-2
                    vv = scal_v_n(n)*r_prev(:, n-1, m) + &
                        & scal_w_n(n)*r_prev(:, n+1, m)
                    vu = ctheta * r_prev(:, n, m)
                    vw = scal_u_n(n)*vu - stheta*vv
                    r(:, n, m) = vw / scal_uvw_m(m)
                    tmp1 = tmp1 + src(ind+n)*r(1, n, m)
                    tmp2 = tmp2 + src(ind-n)*r(2, n, m)
                end do
                dst(ind+m) = alpha * tmp1
                dst(ind-m) = alpha * tmp2
            end do
        end do
    else
        stop "Not Implemented"
    end if
end subroutine fmm_sph_rotate_oxz_work

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_ztranslate(z, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(p+1))
    ! Call corresponding work routine
    call fmm_m2m_ztranslate_work(z, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
end subroutine fmm_m2m_ztranslate

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_m2m_ztranslate_work(z, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, tmp3, res1, res2, pow_r1
    integer :: j, k, n, indj, indjn, indjk1, indjk2
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r2(:)
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_m = zero
        else
            dst_m = beta * dst_m
        end if
        return
    end if
    ! Now alpha is non-zero
    ! If harmonics have different centers
    if (z .ne. 0) then
        ! Prepare pointers
        pow_r2(1:p+1) => work(1:p+1)
        ! Get ratios r1 and r2
        r1 = src_r / dst_r
        r2 = z / dst_r
        ! Get powers of ratio r2/r1
        r2 = r2 / r1
        pow_r1 = r1
        pow_r2(1) = one
        do j = 2, p+1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Do actual M2M
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            do j = 0, p
                ! Offset for dst_m
                indj = j*j + j + 1
                ! k = 0
                tmp1 = alpha * pow_r1 * vscales(indj)
                pow_r1 = pow_r1 * r1
                tmp2 = tmp1
                res1 = zero
                ! Offset for vcnk
                indjk1 = j*(j+1)/2 + 1
                do n = 0, j
                    ! Offset for src_m
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp3 = pow_r2(n+1) / &
                        & vscales(indjn) * vcnk(indjk1+n)**2
                    res1 = res1 + tmp3*src_m(indjn)
                end do
                dst_m(indj) = tmp2 * res1
                ! k != 0
                do k = 1, j
                    tmp2 = tmp1
                    res1 = zero
                    res2 = zero
                    ! Offsets for vcnk
                    indjk1 = (j-k)*(j-k+1)/2 + 1
                    indjk2 = (j+k)*(j+k+1)/2 + 1
                    do n = 0, j-k
                        ! Offset for src_m
                        indjn = (j-n)**2 + (j-n) + 1
                        tmp3 = pow_r2(n+1) / &
                            & vscales(indjn) * vcnk(indjk1+n) * &
                            & vcnk(indjk2+n)
                        res1 = res1 + tmp3*src_m(indjn+k)
                        res2 = res2 + tmp3*src_m(indjn-k)
                    end do
                    dst_m(indj+k) = tmp2 * res1
                    dst_m(indj-k) = tmp2 * res2
                end do
            end do
        ! Update output if beta is non-zero
        else
            do j = 0, p
                ! Offset for dst_m
                indj = j*j + j + 1
                ! k = 0
                tmp1 = alpha * pow_r1 * vscales(indj)
                pow_r1 = pow_r1 * r1
                tmp2 = tmp1
                res1 = zero
                ! Offset for vcnk
                indjk1 = j * (j+1) /2 + 1
                do n = 0, j
                    ! Offset for src_m
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp3 = pow_r2(n+1) / &
                        & vscales(indjn) * vcnk(indjk1+n)**2
                    res1 = res1 + tmp3*src_m(indjn)
                end do
                dst_m(indj) = beta*dst_m(indj) + tmp2*res1
                ! k != 0
                do k = 1, j
                    tmp2 = tmp1
                    res1 = zero
                    res2 = zero
                    ! Offsets for vcnk
                    indjk1 = (j-k)*(j-k+1)/2 + 1
                    indjk2 = (j+k)*(j+k+1)/2 + 1
                    do n = 0, j-k
                        ! Offset for src_m
                        indjn = (j-n)**2 + (j-n) + 1
                        tmp3 = pow_r2(n+1) / &
                            & vscales(indjn) * vcnk(indjk1+n) * &
                            & vcnk(indjk2+n)
                        res1 = res1 + tmp3*src_m(indjn+k)
                        res2 = res2 + tmp3*src_m(indjn-k)
                    end do
                    dst_m(indj+k) = beta*dst_m(indj+k) + tmp2*res1
                    dst_m(indj-k) = beta*dst_m(indj-k) + tmp2*res2
                end do
            end do
        end if
    ! If harmonics are located at the same point
    else
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_m(k) = tmp1 * src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        ! Update output if beta is non-zero
        else
            r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_m(k) = beta*dst_m(k) + tmp1*src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        end if
    end if
end subroutine fmm_m2m_ztranslate_work

!> Adjoint M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of
!! multipole-to-multipole translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_ztranslate_adj(z, src_r, dst_r, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(p+1))
    ! Call corresponding work routine
    call fmm_m2m_ztranslate_adj_work(z, src_r, dst_r, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m, work)
end subroutine fmm_m2m_ztranslate_adj

!> Adjoint M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of
!! multipole-to-multipole translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_m2m_ztranslate_adj_work(z, src_r, dst_r, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m, work)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, res1, res2
    integer :: j, k, n, indj, indn, indjk1, indjk2
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r1(:), pow_r2(:)
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_m = zero
        else
            dst_m = beta * dst_m
        end if
        return
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        ! Prepare pointers
        pow_r1(1:p+1) => work(1:p+1)
        pow_r2(1:p+1) => work(p+2:2*(p+1))
        ! Get powers of r1 and r2
        r1 = dst_r / src_r
        r2 = -z / src_r
        pow_r1(1) = r1
        pow_r2(1) = one
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Do actual adjoint M2M
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            do n = 0, p
                ! Offset for dst_m
                indn = n*n + n + 1
                tmp1 = alpha * pow_r1(n+1) / vscales(indn)
                ! k = 0
                res1 = zero
                do j = n, p
                    ! Offset for src_m
                    indj = j*j + j + 1
                    ! Offsets for vcnk
                    indjk1 = j*(j+1)/2 + 1
                    tmp2 = tmp1 * vscales(indj) * pow_r2(j-n+1) * &
                        & vcnk(indjk1+j-n)**2
                    res1 = res1 + tmp2*src_m(indj)
                end do
                dst_m(indn) = res1
                ! k != 1
                do k = 1, n
                    res1 = zero
                    res2 = zero
                    do j = n, p
                        ! Offset for src_m
                        indj = j*j + j + 1
                        ! Offsets for vcnk
                        indjk1 = (j-k)*(j-k+1)/2 + 1
                        indjk2 = (j+k)*(j+k+1)/2 + 1
                        tmp2 = tmp1 * vscales(indj) * pow_r2(j-n+1) * &
                            & vcnk(indjk1+j-n) * vcnk(indjk2+j-n)
                        res1 = res1 + tmp2*src_m(indj+k)
                        res2 = res2 + tmp2*src_m(indj-k)
                    end do
                    dst_m(indn+k) = res1
                    dst_m(indn-k) = res2
                end do
            end do
        ! Update output if beta is non-zero
        else
            do n = 0, p
                ! Offset for dst_m
                indn = n*n + n + 1
                tmp1 = alpha * pow_r1(n+1) / vscales(indn)
                ! k = 0
                res1 = zero
                do j = n, p
                    ! Offset for src_m
                    indj = j*j + j + 1
                    ! Offsets for vcnk
                    indjk1 = j*(j+1)/2 + 1
                    tmp2 = tmp1 * vscales(indj) * pow_r2(j-n+1) * &
                        & vcnk(indjk1+j-n)**2
                    res1 = res1 + tmp2*src_m(indj)
                end do
                dst_m(indn) = beta*dst_m(indn) + res1
                ! k != 1
                do k = 1, n
                    res1 = zero
                    res2 = zero
                    do j = n, p
                        ! Offset for src_m
                        indj = j*j + j + 1
                        ! Offsets for vcnk
                        indjk1 = (j-k)*(j-k+1)/2 + 1
                        indjk2 = (j+k)*(j+k+1)/2 + 1
                        tmp2 = tmp1 * vscales(indj) * pow_r2(j-n+1) * &
                            & vcnk(indjk1+j-n) * vcnk(indjk2+j-n)
                        res1 = res1 + tmp2*src_m(indj+k)
                        res2 = res2 + tmp2*src_m(indj-k)
                    end do
                    dst_m(indn+k) = beta*dst_m(indn+k) + res1
                    dst_m(indn-k) = beta*dst_m(indn-k) + res2
                end do
            end do
        end if
    ! If harmonics are located at the same point
    else
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            r1 = dst_r / src_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                res1 = zero
                do k = indj-j, indj+j
                    dst_m(k) = tmp1 * src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        ! Update output if beta is non-zero
        else
            r1 = dst_r / src_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                res1 = zero
                do k = indj-j, indj+j
                    dst_m(k) = beta*dst_m(k) + tmp1*src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        end if
    end if
end subroutine fmm_m2m_ztranslate_adj_work

!> Save matrix of M2M translation along OZ axis
!!
!! In a case input `z` is zero no translation matrix is computed, as @ref
!! fmm_m2m_scale shall be used in this case without any precomputed matrices
!!
!!
!! @param[in] z: the OZ coordinate of the radius-vector from new to old centers
!!      of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] mat: Translation matrix for spherical harmonics
subroutine fmm_m2m_ztranslate_get_mat(z, src_r, dst_r, p, vscales, vfact, mat)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), vfact(2*p+1)
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: mat((p+1)*(p+2)*(p+3)/6)
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn, indjn, indmat
    if (z .eq. zero) then
        ! Fill matrix with NaNs
        tmp1 = zero
        mat = tmp1 / tmp1
        return
    end if
    r1 = src_r / dst_r
    r2 = z / dst_r
    pow_r1(1) = r1
    pow_r2(1) = one
    do j = 2, p+1
        pow_r1(j) = pow_r1(j-1) * r1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * vfact(j-k+1) * vfact(j+k+1)
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                mat(indmat) = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                    & vscales(indjn) / vfact(n+1) / vfact(n+1) / &
                    & vfact(j-n-k+1) / vfact(j-n+k+1)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2m_ztranslate_get_mat

!> Apply matrix of M2M translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is an matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_m2m_ztranslate_get_mat in this case as @ref fmm_m2m_scale is intended to
!! treat the case.
!!
!!
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_ztranslate_use_mat(p, mat, alpha, src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(p+2)*(p+3)/6), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    integer :: indmat, j, indj, n, indn, indjn, k
    real(dp) :: tmp1, tmp2
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! Init output properly
        if (beta .eq. zero) then
            dst_m(indj-j:indj+j) = zero
        else
            dst_m(indj-j:indj+j) = beta * dst_m(indj-j:indj+j)
        end if
        ! k = 0
        tmp1 = zero
        do n = 0, j
            indn = n*n + n + 1
            indjn = (j-n)**2 + (j-n) + 1
            tmp1 = tmp1 + mat(indmat)*src_m(indjn)
            indmat = indmat + 1
        end do
        dst_m(indj) = dst_m(indj) + alpha*tmp1
        ! k > 0
        do k = 1, j
            tmp1 = zero
            tmp2 = zero
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                tmp1 = tmp1 + mat(indmat)*src_m(indjn+k)
                tmp2 = tmp2 + mat(indmat)*src_m(indjn-k)
                indmat = indmat + 1
            end do
            dst_m(indj+k) = dst_m(indj+k) + alpha*tmp1
            dst_m(indj-k) = dst_m(indj-k) + alpha*tmp2
        end do
    end do
end subroutine fmm_m2m_ztranslate_use_mat

!> Apply adjoint matrix of M2M translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of a
!! multipole-to-multipole translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_m2m_ztranslate_get_mat in this case as @ref fmm_m2m_scale is intended to
!! treat the case.
!!
!!
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
!!
!!
!! TODO: improve performance
subroutine fmm_m2m_ztranslate_use_mat_adj(p, mat, alpha, src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(p+2)*(p+3)/6), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    integer :: indmat, j, indj, n, indn, indjn, k
    real(dp) :: tmp1, tmp2
    ! Init output properly
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = 0, j
            indn = n*n + n + 1
            indjn = (j-n)**2 + (j-n) + 1
            dst_m(indjn) = dst_m(indjn) + alpha*mat(indmat)*src_m(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                dst_m(indjn+k) = dst_m(indjn+k) + &
                    & alpha*mat(indmat)*src_m(indj+k)
                dst_m(indjn-k) = dst_m(indjn-k) + &
                    & alpha*mat(indmat)*src_m(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2m_ztranslate_use_mat_adj

!> Scale M2M, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation when harmonics are centered in the same point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_scale(src_r, dst_r, p, alpha, src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: src_r, dst_r, alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    r1 = src_r / dst_r
    tmp1 = alpha * r1
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_m(k) = dst_m(k) + src_m(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_m2m_scale

!> Adjoint scale M2M, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of a
!! multipole-to-multipole translation when harmonics are centered in the same
!! point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_scale_adj(src_r, dst_r, p, alpha, src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: src_r, dst_r, alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    r1 = dst_r / src_r
    tmp1 = alpha * r1
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_m(k) = dst_m(k) + src_m(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_m2m_scale_adj

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_rotation(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p + 19*p + 8)
    ! Call corresponding work routine
    call fmm_m2m_rotation_work(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
end subroutine fmm_m2m_rotation

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_m2m_rotation_work(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_m(:), tmp_m2(:), vcos(:), vsin(:)
    ! Convert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        ! Workspace here is 2*(p+1)
        call fmm_m2m_ztranslate_work(c(3), src_r, dst_r, p, vscales, vcnk, &
            & alpha, src_m, beta, dst_m, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_m(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_m2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_m, zero, tmp_m)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_m, zero, &
        & tmp_m2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_m2m_ztranslate_work(rho, src_r, dst_r, p, vscales, vcnk, one, &
        & tmp_m2, zero, tmp_m, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_m, zero, tmp_m2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_rotation_work

!> Adjoint M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is aa adjoint matrix of a
!! multipole-to-multipole translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] src_m: expansion in old harmonics
!! @param[inout] dst_m: expansion in new harmonics
subroutine fmm_m2m_rotation_adj(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p + 19*p + 8)
    ! Call corresponding work routine
    call fmm_m2m_rotation_adj_work(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
end subroutine fmm_m2m_rotation_adj

!> Adjoint M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is aa adjoint matrix of a
!! multipole-to-multipole translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] src_m: expansion in old harmonics
!! @param[inout] dst_m: expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_m2m_rotation_adj_work(c, src_r, dst_r, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_m(:), tmp_m2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. 0) then
        ! Workspace here is 2*(p+1)
        call fmm_m2m_ztranslate_adj_work(c(3), src_r, dst_r, p, vscales, &
            & vcnk, alpha, src_m, beta, dst_m, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_m(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_m2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_m, zero, tmp_m)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_m, zero, &
        & tmp_m2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_m2m_ztranslate_adj_work(rho, src_r, dst_r, p, vscales, vcnk, &
        & one, tmp_m2, zero, tmp_m, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_m, zero, tmp_m2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_rotation_adj_work

!> Compute reflection that aligns given vector along OZ axis
!!
!! Computes stable reflection, just as in Householder QR algorithm. This
!! procedure is for testing purposes here, as it is already inplace in
!! performance-oriented implementations.
!!
!! @param[in] c: Input vector to align along OZ axis
!! @param[out] z: OZ coordinate of the reflected vector `c`
!! @param[out] mat: Reflection matrix from new (y,z,x) to old (y,z,x)
subroutine coord_reflect_get_mat(c, z, mat)
    ! Input
    real(dp), intent(in) :: c(3)
    ! Output
    real(dp), intent(out) :: z, mat(3, 3)
    ! Local variables
    real(dp) :: max123, ssq123, c1(3), c1_norm
    integer :: m, n
    ! Compute rho (norm of the input vector c)
    max123 = abs(c(1))
    ssq123 = one
    if (c(2) .ne. zero) then
        if (abs(c(2)) .gt. max123) then
            ssq123 = one + ssq123*(max123/c(2))**2
            max123 = abs(c(2))
        else
            ssq123 = ssq123 + (c(2)/max123)**2
        end if
    end if
    if (c(3) .ne. zero) then
        if (abs(c(3)) .gt. max123) then
            ssq123 = one + ssq123*(max123/c(3))**2
            max123 = abs(c(3))
        else
            ssq123 = ssq123 + (c(3)/max123)**2
        end if
    end if
    z = max123 * sqrt(ssq123)
    ! Reorder (x,y,z) -> (y,z,x) and set proper reflection normal vector c1
    if (c(3) .ge. zero) z = -z
    c1(1) = c(2)
    c1(2) = c(3) - z ! -z and c(3) have the same sign
    c1(3) = c(1)
    ! Normalize vector c1. We know for sure c1(2) is maximum and is non-zero
    max123 = abs(c1(2))
    ssq123 = one + (c1(1)/c1(2))**2 + (c1(3)/c1(2))**2
    c1_norm = max123 * sqrt(ssq123)
    c1 = c1 / c1_norm
    mat = zero
    do m = 1, 3
        mat(m, m) = 1
        do n = 1, 3
            mat(n, m) = mat(n, m) - two*c1(n)*c1(m)
        end do
    end do
end subroutine coord_reflect_get_mat

!> Save matrices of M2M operation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[out] transform_mat: Matrix of reflection
!! @param[out] ztranslate_mat: Matrix of OZ translation
subroutine fmm_m2m_reflection_get_mat(c, src_r, dst_r, p, vscales, vfact, &
        & transform_mat, ztranslate_mat)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1)
    integer, intent(in) :: p
    ! Outputs
    real(dp), intent(out) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(dp), intent(out) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    ! Local variables
    real(dp) :: rho, r1(3, 3)
    ! If no need for transformation, just generate translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If even z coordinate is zero do nothing, as this case requires a
        ! simple scaling
        if (c(3) .ne. zero) then
            call fmm_m2m_ztranslate_get_mat(c(3), src_r, dst_r, p, vscales, &
                & vfact, ztranslate_mat)
        end if
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform_get_mat(p, r1, transform_mat)
    call fmm_m2m_ztranslate_get_mat(rho, src_r, dst_r, p, vscales, vfact, &
        & ztranslate_mat)
end subroutine fmm_m2m_reflection_get_mat

!> Apply matrices of M2M translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion of new harmonics
subroutine fmm_m2m_reflection_use_mat(c, src_r, dst_r, p, transform_mat, &
        & ztranslate_mat, alpha, src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, &
        & transform_mat((p+1)*(2*p+1)*(2*p+3)/3), &
        & ztranslate_mat((p+1)*(p+2)*(p+3)/6), &
        & alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: tmp_m((p+1)*(p+1)), tmp_m2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same just scale
        if (c(3) .eq. 0) then
            call fmm_m2m_scale(src_r, dst_r, p, alpha, src_m, beta, dst_m)
        ! Otherwise apply ztranslate matrix
        else
            call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, alpha, src_m, &
                & beta, dst_m)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(p, transform_mat, alpha, src_m, zero, tmp_m)
    call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, tmp_m, zero, &
        & tmp_m2)
    call fmm_sph_transform_use_mat(p, transform_mat, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_reflection_use_mat

!> Adjoint apply matrices of M2M translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion of new harmonics
subroutine fmm_m2m_reflection_use_mat_adj(c, src_r, dst_r, p, transform_mat, &
        & ztranslate_mat, alpha, src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, &
        & transform_mat((p+1)*(2*p+1)*(2*p+3)/3), &
        & ztranslate_mat((p+1)*(p+2)*(p+3)/6), &
        & alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: tmp_m((p+1)*(p+1)), tmp_m2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same just scale
        if (c(3) .eq. 0) then
            call fmm_m2m_scale_adj(src_r, dst_r, p, alpha, src_m, beta, dst_m)
        ! Otherwise apply ztranslate matrix
        else
            call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, alpha, &
                & src_m, beta, dst_m)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(p, transform_mat, alpha, src_m, zero, tmp_m)
    call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, one, tmp_m, zero, &
        & tmp_m2)
    call fmm_sph_transform_use_mat(p, transform_mat, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_reflection_use_mat_adj

!> Direct L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of local-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_ztranslate(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(p+1))
    ! Call corresponding work routine
    call fmm_l2l_ztranslate_work(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
end subroutine fmm_l2l_ztranslate

!> Direct L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of local-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_l2l_ztranslate_work(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2
    integer :: j, k, n, indj, indn
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r1(:), pow_r2(:)
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        ! Prepare pointers
        pow_r1(1:p+1) => work(1:p+1)
        pow_r2(1:p+1) => work(p+2:2*(p+1))
        ! Get powers of r1 and r2
        r1 = z / src_r
        r2 = dst_r / z
        pow_r1(1) = 1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Do actual L2L
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = alpha * pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1) * &
                    & vscales(indj)
                do n = j, p
                    indn = n*n + n + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) * &
                        & vfact(n-k+1) * vfact(n+k+1) / vfact(n-j+1) / &
                        & vfact(n-j+1)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    if (k .eq. 0) then
                        dst_l(indj) = dst_l(indj) + tmp2*src_l(indn)
                    else
                        dst_l(indj+k) = dst_l(indj+k) + tmp2*src_l(indn+k)
                        dst_l(indj-k) = dst_l(indj-k) + tmp2*src_l(indn-k)
                    end if
                end do
            end do
        end do
    ! If harmonics are located at the same point
    else
        r1 = dst_r / src_r
        tmp1 = alpha
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_l2l_ztranslate_work

!> Adjoint L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of
!! local-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_ztranslate_adj(z, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(p+1))
    ! Call corresponding work routine
    call fmm_l2l_ztranslate_adj_work(z, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l, work)
end subroutine fmm_l2l_ztranslate_adj

!> Adjoint L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of
!! local-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_l2l_ztranslate_adj_work(z, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l, work)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2
    integer :: j, k, n, indj, indn
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r1(:), pow_r2(:)
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        ! Prepare pointers
        pow_r1(1:p+1) => work(1:p+1)
        pow_r2(1:p+1) => work(p+2:2*(p+1))
        ! Get powers of r1 and r2
        r1 = -z / dst_r
        r2 = -src_r / z
        pow_r1(1) = one
        pow_r2(1) = one
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = alpha * pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1) * &
                    & vscales(indj)
                do n = j, p
                    indn = n*n + n + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) * &
                        & vfact(n-k+1) * vfact(n+k+1) / vfact(n-j+1) / &
                        & vfact(n-j+1)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    if (k .eq. 0) then
                        dst_l(indn) = dst_l(indn) + tmp2*src_l(indj)
                    else
                        dst_l(indn+k) = dst_l(indn+k) + tmp2*src_l(indj+k)
                        dst_l(indn-k) = dst_l(indn-k) + tmp2*src_l(indj-k)
                    end if
                end do
            end do
        end do
    ! If harmonics are located at the same point
    else
        r1 = src_r / dst_r
        tmp1 = alpha
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_l2l_ztranslate_adj_work

!> Save matrix of L2L translation along OZ axis
!!
!! In a case input `z` is zero no translation matrix is computed, as @ref
!! fmm_l2l_scale shall be used in this case without any precomputed matrices
!!
!!
!! @param[in] z: the OZ coordinate of the radius-vector from new to old centers
!!      of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] mat: Translation matrix for spherical harmonics
subroutine fmm_l2l_ztranslate_get_mat(z, src_r, dst_r, p, vscales, vfact, mat)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), vfact(2*p+1)
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: mat((p+1)*(p+2)*(p+3)/6)
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn, indmat
    if (z .eq. zero) then
        ! Fill matrix with NaNs
        tmp1 = zero
        mat = tmp1 / tmp1
        return
    end if
    r1 = z / src_r
    r2 = dst_r / z
    pow_r1(1) = one
    pow_r2(1) = one
    do j = 2, p+1
        pow_r1(j) = pow_r1(j-1) * r1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1) * vscales(indj)
            do n = j, p
                indn = n*n + n + 1
                mat(indmat) = tmp1 * pow_r1(n+1) / vscales(indn) * &
                    & vfact(n-k+1) * vfact(n+k+1) / vfact(n-j+1) / &
                    & vfact(n-j+1)
                if (mod(n+j, 2) .eq. 1) then
                    mat(indmat) = -mat(indmat)
                end if
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_l2l_ztranslate_get_mat

!> Apply matrix of L2L translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of local-to-local
!! translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_l2l_ztranslate_get_mat in this case as @ref fmm_l2l_scale is intended to
!! treat the case.
!!
!!
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Multipole expansion in new harmonics
subroutine fmm_l2l_ztranslate_use_mat(p, mat, alpha, src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(p+2)*(p+3)/6), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    integer :: j, k, n, indj, indn, indmat
    real(dp) :: tmp1, tmp2
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! Init output properly
        if (beta .eq. zero) then
            dst_l(indj-j:indj+j) = zero
        else
            dst_l(indj-j:indj+j) = beta * dst_l(indj-j:indj+j)
        end if
        ! k = 0
        tmp1 = zero
        do n = j, p
            indn = n*n + n + 1
            tmp1 = tmp1 + mat(indmat)*src_l(indn)
            indmat = indmat + 1
        end do
        dst_l(indj) = dst_l(indj) + alpha*tmp1
        ! k > 0
        do k = 1, j
            tmp1 = zero
            tmp2 = zero
            do n = j, p
                indn = n*n + n + 1
                tmp1 = tmp1 + mat(indmat)*src_l(indn+k)
                tmp2 = tmp2 + mat(indmat)*src_l(indn-k)
                indmat = indmat + 1
            end do
            dst_l(indj+k) = dst_l(indj+k) + alpha*tmp1
            dst_l(indj-k) = dst_l(indj-k) + alpha*tmp2
        end do
    end do
end subroutine fmm_l2l_ztranslate_use_mat

!> Apply adjoint matrix of L2L translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_l2l_ztranslate_get_mat in this case as @ref fmm_l2l_scale_adj is
!! intended to treat the case.
!!
!!
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
!!
!!
!! TODO: improve performance
subroutine fmm_l2l_ztranslate_use_mat_adj(p, mat, alpha, src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(p+2)*(p+3)/6), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    integer :: j, k, n, indj, indn, indmat
    real(dp) :: tmp1, tmp2
    ! Init output properly
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = j, p
            indn = n*n + n + 1
            dst_l(indn) = dst_l(indn) + alpha*mat(indmat)*src_l(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = j, p
                indn = n*n + n + 1
                dst_l(indn+k) = dst_l(indn+k) + alpha*mat(indmat)*src_l(indj+k)
                dst_l(indn-k) = dst_l(indn-k) + alpha*mat(indmat)*src_l(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_l2l_ztranslate_use_mat_adj

!> Scale L2L, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of a local-to-local
!! translation when harmonics are centered in the same point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_l: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_l: Multipole expansion in new harmonics
subroutine fmm_l2l_scale(src_r, dst_r, p, alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: src_r, dst_r, alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    r1 = dst_r / src_r
    tmp1 = alpha
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_l(k) = dst_l(k) + src_l(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_l2l_scale

!> Adjoint scale L2L, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation when harmonics are centered in the same point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_l: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_l: Multipole expansion in new harmonics
subroutine fmm_l2l_scale_adj(src_r, dst_r, p, alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: src_r, dst_r, alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    r1 = src_r / dst_r
    tmp1 = alpha
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_l(k) = dst_l(k) + src_l(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_l2l_scale_adj

!> Direct L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of a local-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_rotation(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p+19*p+8)
    ! Call corresponding work routine
    call fmm_l2l_rotation_work(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
end subroutine fmm_l2l_rotation

!> Direct L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of a local-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_l2l_rotation_work(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_l(:), tmp_l2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        ! Workspace here is 2*(p+1)
        call fmm_l2l_ztranslate_work(c(3), src_r, dst_r, p, vscales, vfact, &
            & alpha, src_l, beta, dst_l, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_l(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_l2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_l, zero, tmp_l)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_l, zero, &
        & tmp_l2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_l2l_ztranslate_work(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_l2, zero, tmp_l, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_l, zero, tmp_l2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_rotation_work

!> Adjoint L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_rotation_adj(c, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p + 19*p + 8)
    ! Call corresponding work routine
    call fmm_l2l_rotation_adj_work(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
end subroutine fmm_l2l_rotation_adj

!> Adjoint L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_l2l_rotation_adj_work(c, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_l(:), tmp_l2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        ! Workspace here is 2*(p+1)
        call fmm_l2l_ztranslate_adj_work(c(3), src_r, dst_r, p, vscales, &
            & vfact, alpha, src_l, beta, dst_l, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_l(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_l2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_l, zero, tmp_l)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_l, zero, &
        & tmp_l2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_l2l_ztranslate_adj_work(rho, src_r, dst_r, p, vscales, vfact, &
        & one, tmp_l2, zero, tmp_l, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_l, zero, tmp_l2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_rotation_adj_work

!> Save matrices of L2L operation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[out] transform_mat: Matrix of reflection
!! @param[out] ztranslate_mat: Matrix of OZ translation
subroutine fmm_l2l_reflection_get_mat(c, src_r, dst_r, p, vscales, vfact, &
        & transform_mat, ztranslate_mat)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1)
    integer, intent(in) :: p
    ! Outputs
    real(dp), intent(out) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/6)
    real(dp), intent(out) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    ! Local variables
    real(dp) :: rho, r1(3, 3)
    ! If no need for transformation, just generate translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If even z coordinate is zero do nothing, as this case requires a
        ! simple scaling
        if (c(3) .ne. zero) then
            call fmm_l2l_ztranslate_get_mat(c(3), src_r, dst_r, p, vscales, &
                & vfact, ztranslate_mat)
        end if
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform_get_mat(p, r1, transform_mat)
    call fmm_l2l_ztranslate_get_mat(rho, src_r, dst_r, p, vscales, vfact, &
        & ztranslate_mat)
end subroutine fmm_l2l_reflection_get_mat

!> Apply matrices of L2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Multipole expansion of new harmonics
subroutine fmm_l2l_reflection_use_mat(c, src_r, dst_r, p, transform_mat, &
        & ztranslate_mat, alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, &
        & transform_mat((p+1)*(2*p+1)*(2*p+3)/3), &
        & ztranslate_mat((p+1)*(p+2)*(p+3)/6), &
        & alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: tmp_l((p+1)*(p+1)), tmp_l2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same just scale
        if (c(3) .eq. 0) then
            call fmm_l2l_scale(src_r, dst_r, p, alpha, src_l, beta, dst_l)
        ! Otherwise apply ztranslate matrix
        else
            call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, alpha, src_l, &
                & beta, dst_l)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(p, transform_mat, alpha, src_l, zero, tmp_l)
    call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, one, tmp_l, zero, &
        & tmp_l2)
    call fmm_sph_transform_use_mat(p, transform_mat, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_reflection_use_mat

!> Adjoint apply matrices of L2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Multipole expansion of new harmonics
subroutine fmm_l2l_reflection_use_mat_adj(c, src_r, dst_r, p, transform_mat, &
        & ztranslate_mat, alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, &
        & transform_mat((p+1)*(2*p+1)*(2*p+3)/3), &
        & ztranslate_mat((p+1)*(p+2)*(p+3)/6), &
        & alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: tmp_l((p+1)*(p+1)), tmp_l2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same just scale
        if (c(3) .eq. 0) then
            call fmm_l2l_scale_adj(src_r, dst_r, p, alpha, src_l, beta, dst_l)
        ! Otherwise apply ztranslate matrix
        else
            call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, alpha, &
                & src_l, beta, dst_l)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(p, transform_mat, alpha, src_l, zero, tmp_l)
    call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, one, tmp_l, zero, &
        & tmp_l2)
    call fmm_sph_transform_use_mat(p, transform_mat, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_reflection_use_mat_adj

!> Direct M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of multipole-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (multipole) harmonics
!! @param[in] dst_r: Radius of new (local) harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old (multipole) harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new (local) harmonics
subroutine fmm_m2l_ztranslate(z, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
        & beta
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Temporary workspace
    real(dp) :: work((pm+2)*(pm+1))
    ! Call corresponding work routine
    call fmm_m2l_ztranslate_work(z, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
end subroutine fmm_m2l_ztranslate

!> Direct M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of multipole-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (multipole) harmonics
!! @param[in] dst_r: Radius of new (local) harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old (multipole) harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new (local) harmonics
!! @param[out] work: Temporary workspace of a size (pm+2)*(pm+1)
subroutine fmm_m2l_ztranslate_work(z, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
        & beta
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work((pm+2)*(pm+1))
    ! Local variables
    real(dp) :: tmp1, tmp2, r1, r2, res1, res2, pow_r2
    integer :: j, k, n, indj, indjn, indk1, indk2
    ! Pointers for temporary values of powers
    real(dp), pointer :: src_m2(:), pow_r1(:)
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_l = zero
        else
            dst_l = beta * dst_l
        end if
        return
    end if
    ! Now alpha is non-zero
    ! z cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (z .eq. zero) then
        return
    end if
    ! Prepare pointers
    n = (pm+1) ** 2
    src_m2(1:n) => work(1:n)
    pow_r1(1:pm+1) => work(n+1:n+pm+1)
    ! Get powers of r1 and r2
    r1 = src_r / z
    r2 = dst_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    tmp1 = abs(r1)
    do j = 0, pm
        indj = j*j + j + 1
        pow_r1(j+1) = tmp1 / vscales(indj)
        tmp1 = tmp1 * r1
    end do
    pow_r2 = one
    ! Reorder source harmonics from (degree, order) to (order, degree)
    ! Zero order k=0 at first
    do j = 0, pm
        indj = j*j + j + 1
        src_m2(j+1) = pow_r1(j+1) * src_m(indj)
    end do
    ! Non-zero orders next, a positive k followed by a negative -k
    indk1 = pm + 2
    do k = 1, pm
        n = pm - k + 1
        indk2 = indk1 + n
        do j = k, pm
            indj = j*j + j + 1
            src_m2(indk1+j-k) = pow_r1(j+1) * src_m(indj+k)
            src_m2(indk2+j-k) = pow_r1(j+1) * src_m(indj-k)
        end do
        indk1 = indk2 + n
    end do
    ! Do actual M2L
    ! Overwrite output if beta is zero
    if (beta .eq. zero) then
        do j = 0, pl
            ! Offset for dst_l
            indj = j*j + j + 1
            ! k = 0
            tmp1 = alpha * vscales(indj) * pow_r2
            pow_r2 = pow_r2 * r2
            res1 = zero
            do n = 0, pm
                res1 = res1 + m2l_ztranslate_coef(n+1, 1, j+1)*src_m2(n+1)
            end do
            dst_l(indj) = tmp1 * res1
            ! k != 0
            do k = 1, j
                ! Offsets for src_m2
                indk1 = pm + 2 + (2*pm-k+2)*(k-1)
                indk2 = indk1 + pm - k + 1
                res1 = zero
                res2 = zero
                do n = k, pm
                    res1 = res1 + &
                        & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                        & src_m2(indk1+n-k)
                    res2 = res2 + &
                        & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                        & src_m2(indk2+n-k)
                end do
                dst_l(indj+k) = tmp1 * res1
                dst_l(indj-k) = tmp1 * res2
            end do
        end do
    else
        do j = 0, pl
            ! Offset for dst_l
            indj = j*j + j + 1
            ! k = 0
            tmp1 = alpha * vscales(indj) * pow_r2
            pow_r2 = pow_r2 * r2
            res1 = zero
            do n = 0, pm
                res1 = res1 + m2l_ztranslate_coef(n+1, 1, j+1)*src_m2(n+1)
            end do
            dst_l(indj) = beta*dst_l(indj) + tmp1*res1
            ! k != 0
            do k = 1, j
                ! Offsets for src_m2
                indk1 = pm + 2 + (2*pm-k+2)*(k-1)
                indk2 = indk1 + pm - k + 1
                res1 = zero
                res2 = zero
                do n = k, pm
                    res1 = res1 + &
                        & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                        & src_m2(indk1+n-k)
                    res2 = res2 + &
                        & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                        & src_m2(indk2+n-k)
                end do
                dst_l(indj+k) = beta*dst_l(indj+k) + tmp1*res1
                dst_l(indj-k) = beta*dst_l(indj-k) + tmp1*res2
            end do
        end do
    end if
end subroutine fmm_m2l_ztranslate_work

!> Adjoint M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of
!! multipole-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (local) harmonics
!! @param[in] dst_r: Radius of new (multipole) harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] m2l_ztranslate_adj_coef:
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old (local) harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new (multipole) harmonics
subroutine fmm_m2l_ztranslate_adj(z, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), alpha, &
        & src_l((pl+1)*(pl+1)), beta
    integer, intent(in) :: pl, pm
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Temporary workspace
    real(dp) :: work((pl+2)*(pl+1))
    ! Call corresponding work routine
    call fmm_m2l_ztranslate_adj_work(z, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
end subroutine fmm_m2l_ztranslate_adj

!> Adjoint M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of
!! multipole-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (local) harmonics
!! @param[in] dst_r: Radius of new (multipole) harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] m2l_ztranslate_adj_coef:
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old (local) harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new (multipole) harmonics
subroutine fmm_m2l_ztranslate_adj_work(z, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), alpha, &
        & src_l((pl+1)*(pl+1)), beta
    integer, intent(in) :: pl, pm
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work((pl+2)*(pl+1))
    ! Local variables
    real(dp) :: tmp1, tmp2, tmp3, r1, r2, pow_r1, res1, res2
    integer :: j, k, n, indj, indn, indk1, indk2
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r2(:), src_l2(:)
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_m = zero
        else
            dst_m = beta * dst_m
        end if
        return
    end if
    ! Now alpha is non-zero
    ! z cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (z .eq. zero) then
        return
    end if
    ! Prepare pointers
    n = (pl+1) ** 2
    src_l2(1:n) => work(1:n)
    pow_r2(1:pl+1) => work(n+1:n+pl+1)
    ! Get powers of r1 and r2
    r1 = -dst_r / z
    r2 = -src_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    tmp1 = one
    do j = 0, pl
        indj = j*j + j + 1
        pow_r2(j+1) = tmp1 * vscales(indj)
        tmp1 = tmp1 * r2
    end do
    pow_r1 = abs(r1)
    ! Reorder source harmonics from (degree, order) to (order, degree)
    ! Zero order k=0 at first
    do j = 0, pl
        indj = j*j + j + 1
        src_l2(j+1) = pow_r2(j+1) * src_l(indj)
    end do
    ! Non-zero orders next, a positive k followed by a negative -k
    indk1 = pl + 2
    do k = 1, pl
        n = pl - k + 1
        indk2 = indk1 + n
        do j = k, pl
            indj = j*j + j + 1
            src_l2(indk1+j-k) = pow_r2(j+1) * src_l(indj+k)
            src_l2(indk2+j-k) = pow_r2(j+1) * src_l(indj-k)
        end do
        indk1 = indk2 + n
    end do
    ! Do actual adjoint M2L
    ! Overwrite output if beta is zero
    if (beta .eq. zero) then
        do n = 0, pm
            indn = n*n + n + 1
            ! k = 0
            tmp1 = alpha * pow_r1 / vscales(indn)
            pow_r1 = pow_r1 * r1
            res1 = zero
            do j = 0, pl
                indj = j*j + j + 1
                res1 = res1 + m2l_ztranslate_adj_coef(j+1, 1, n+1)*src_l2(j+1)
            end do
            dst_m(indn) = tmp1 * res1
            ! k != 0
            do k = 1, n
                ! Offsets for src_l2
                indk1 = pl + 2 + (2*pl-k+2)*(k-1)
                indk2 = indk1 + pl - k + 1
                res1 = zero
                res2 = zero
                do j = k, pl
                    indj = j*j + j + 1
                    res1 = res1 + &
                        & m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1)* &
                        & src_l2(indk1+j-k)
                    res2 = res2 + &
                        & m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1)* &
                        & src_l2(indk2+j-k)
                end do
                dst_m(indn+k) = tmp1 * res1
                dst_m(indn-k) = tmp1 * res2
            end do
        end do
    ! Update output if beta is non-zero
    else
        do n = 0, pm
            indn = n*n + n + 1
            ! k = 0
            tmp1 = alpha * pow_r1 / vscales(indn)
            pow_r1 = pow_r1 * r1
            res1 = zero
            do j = 0, pl
                indj = j*j + j + 1
                res1 = res1 + m2l_ztranslate_adj_coef(j+1, 1, n+1)*src_l2(j+1)
            end do
            dst_m(indn) = beta*dst_m(indn) + tmp1*res1
            ! k != 0
            do k = 1, n
                ! Offsets for src_l2
                indk1 = pl + 2 + (2*pl-k+2)*(k-1)
                indk2 = indk1 + pl - k + 1
                res1 = zero
                res2 = zero
                do j = k, pl
                    indj = j*j + j + 1
                    res1 = res1 + &
                        & m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1)* &
                        & src_l2(indk1+j-k)
                    res2 = res2 + &
                        & m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1)* &
                        & src_l2(indk2+j-k)
                end do
                dst_m(indn+k) = beta*dst_m(indn+k) + tmp1*res1
                dst_m(indn-k) = beta*dst_m(indn-k) + tmp1*res2
            end do
        end do
    end if
end subroutine fmm_m2l_ztranslate_adj_work

!> Save matrix of M2L translation along OZ axis
!!
!! In a case input `z` is zero no translation matrix is computed, as @ref
!! fmm_m2m_scale shall be used in this case without any precomputed matrices
!!
!!
!! @param[in] z: the OZ coordinate of the radius-vector from new to old centers
!!      of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] mat: Translation matrix for spherical harmonics
subroutine fmm_m2l_ztranslate_get_mat(z, src_r, dst_r, pm, pl, vscales, &
        & vfact, mat)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & vfact(2*(pm+pl)+1)
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(out) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    ! Local variables
    real(dp) :: r1, r2, tmp1, pow_r1(pm+1), pow_r2(pl+1)
    integer :: j, k, n, indj, indn, indjn, indmat
    if (z .eq. zero) then
        ! M and L cannot be located in the same place
        return
    end if
    r1 = src_r / z
    r2 = dst_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    pow_r1(1) = abs(r1)
    pow_r2(1) = one
    do j = 2, pm+1
        pow_r1(j) = pow_r1(j-1) * r1
    end do
    do j = 2, pl+1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1)
            do n = k, pm
                indn = n*n + n + 1
                mat(indmat) = tmp1 * pow_r1(n+1) / vscales(indn) / &
                    & vfact(n-k+1) / vfact(n+k+1) * vfact(j+n+1) * vfact(j+n+1)
                if (mod(n+k, 2) .eq. 1) then
                    mat(indmat) = -mat(indmat)
                end if
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2l_ztranslate_get_mat

!> Apply matrix of M2L translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_m2l_ztranslate_get_mat in this case as @ref fmm_m2l_scale is intended to
!! treat the case.
!!
!!
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_m2l_ztranslate_use_mat(pm, pl, mat, alpha, src_m, beta, dst_l)
    ! Inputs
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6), alpha, src_m((pm+1)*(pm+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Local variables
    integer :: j, k, n, indj, indn, indjn, indmat
    real(dp) :: tmp1, tmp2
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        ! Init output properly
        if (beta .eq. zero) then
            dst_l(indj-j:indj+j) = zero
        else
            dst_l(indj-j:indj+j) = beta * dst_l(indj-j:indj+j)
        end if
        ! k = 0
        tmp1 = zero
        do n = 0, pm
            indn = n*n + n + 1
            tmp1 = tmp1 + mat(indmat)*src_m(indn)
            indmat = indmat + 1
        end do
        dst_l(indj) = dst_l(indj) + alpha*tmp1
        ! k > 0
        do k = 1, j
            tmp1 = zero
            tmp2 = zero
            do n = k, pm
                indn = n*n + n + 1
                tmp1 = tmp1 + mat(indmat)*src_m(indn+k)
                tmp2 = tmp2 + mat(indmat)*src_m(indn-k)
                indmat = indmat + 1
            end do
            dst_l(indj+k) = dst_l(indj+k) + alpha*tmp1
            dst_l(indj-k) = dst_l(indj-k) + alpha*tmp2
        end do
    end do
end subroutine fmm_m2l_ztranslate_use_mat

!> Apply adjoint matrix of M2L translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of a
!! multipole-to-local translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_m2l_ztranslate_get_mat in this case as @ref fmm_m2l_scale_adj is
!! intended to treat the case.
!!
!!
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_m: Multipole expansion in new harmonics
!!
!!
!! TODO: improve performance
subroutine fmm_m2l_ztranslate_use_mat_adj(pl, pm, mat, alpha, src_l, beta, &
        & dst_m)
    ! Inputs
    integer, intent(in) :: pl, pm
    real(dp), intent(in) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6), alpha, src_l((pl+1)*(pl+1)), &
        & beta
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Local variables
    integer :: j, k, n, indj, indn, indjn, indmat
    ! Init output properly
    if (beta .eq. zero) then
        dst_m(:) = zero
    else
        dst_m(:) = beta * dst_m(:)
    end if
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        ! k = 0
        do n = 0, pm
            indn = n*n + n + 1
            dst_m(indn) = dst_m(indn) + alpha*mat(indmat)*src_l(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = k, pm
                indn = n*n + n + 1
                dst_m(indn+k) = dst_m(indn+k) + alpha*mat(indmat)*src_l(indj+k)
                dst_m(indn-k) = dst_m(indn-k) + alpha*mat(indmat)*src_l(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2l_ztranslate_use_mat_adj

!> Save matrices of M2L operation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[out] transform_mat: Matrix of reflection
!! @param[out] ztranslate_mat: Matrix of OZ translation
subroutine fmm_m2l_reflection_get_mat(c, src_r, dst_r, pm, pl, vscales, &
        & vfact, transform_mat, ztranslate_mat)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & vfact(2*(pm+pl)+1)
    integer, intent(in) :: pm, pl
    ! Outputs
    real(dp), intent(out) :: transform_mat((max(pm,pl)+1)*(2*max(pm,pl)+1)* &
        & (2*max(pm,pl)+3)/3)
    real(dp), intent(out) :: ztranslate_mat((min(pm,pl)+1)*(min(pm,pl)+2)* &
        & (3*max(pm,pl)+3-min(pm,pl))/6)
    ! Local variables
    real(dp) :: rho, r1(3, 3)
    ! If no need for transformation, just generate translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If M and L harmonics are at the same place ignore M2L completely
        if (c(3) .ne. zero) then
            call fmm_m2l_ztranslate_get_mat(c(3), src_r, dst_r, pm, pl, &
                & vscales, vfact, ztranslate_mat)
        end if
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform_get_mat(max(pm,pl), r1, transform_mat)
    call fmm_m2l_ztranslate_get_mat(rho, src_r, dst_r, pm, pl, vscales, &
        & vfact, ztranslate_mat)
end subroutine fmm_m2l_reflection_get_mat

!> Apply matrices of M2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Multipole expansion of new harmonics
subroutine fmm_m2l_reflection_use_mat(c, src_r, dst_r, pm, pl, transform_mat, &
        & ztranslate_mat, alpha, src_m, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, alpha, src_m((pm+1)*(pm+1)), &
        & beta
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: &
        & transform_mat((max(pm,pl)+1)*(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3), &
        & ztranslate_mat((min(pm,pl)+1)*(min(pm,pl)+2)* &
        & (3*max(pm,pl)+3-min(pm,pl))/6)
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Local variables
    real(dp) :: tmp_m((pm+1)*(pm+1)), tmp_l((pl+1)*(pl+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! Only apply ztranslate matrix if c is non-zero
        if (c(3) .ne. 0) then
            call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, alpha, &
                & src_m, beta, dst_l)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(pm, transform_mat, alpha, src_m, zero, &
        & tmp_m)
    call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, tmp_m, zero, &
        & tmp_l)
    call fmm_sph_transform_use_mat(pl, transform_mat, one, tmp_l, beta, dst_l)
end subroutine fmm_m2l_reflection_use_mat

!> Adjoint apply matrices of M2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Local expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion of new harmonics
subroutine fmm_m2l_reflection_use_mat_adj(c, src_r, dst_r, pl, pm, &
        & transform_mat, ztranslate_mat, alpha, src_l, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, alpha, src_l((pl+1)*(pl+1)), &
        & beta
    integer, intent(in) :: pl, pm
    real(dp), intent(in) :: &
        & transform_mat((max(pm,pl)+1)*(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3), &
        & ztranslate_mat((min(pm,pl)+1)*(min(pm,pl)+2)* &
        & (3*max(pm,pl)+3-min(pm,pl))/6)
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Local variables
    real(dp) :: tmp_m((pm+1)*(pm+1)), tmp_l((pl+1)*(pl+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same ignore M2L completely
        if (c(3) .ne. 0) then
            call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, &
                & alpha, src_l, beta, dst_m)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(pl, transform_mat, alpha, src_l, zero, &
        & tmp_l)
    call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, one, tmp_l, &
        & zero, tmp_m)
    call fmm_sph_transform_use_mat(pm, transform_mat, one, tmp_m, beta, dst_m)
end subroutine fmm_m2l_reflection_use_mat_adj

!> Direct M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of a multipole-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_m2l_rotation(c, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
        & beta
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Temporary workspace
    real(dp) :: work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8)
    ! Call corresponding work routine
    call fmm_m2l_rotation_work(c, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
end subroutine fmm_m2l_rotation

!> Direct M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of a multipole-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8 where p is a
!!      maximum of pm and pl
subroutine fmm_m2l_rotation_work(c, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
        & beta
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Temporary workspace
    real(dp), intent(out), target :: &
        & work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n, p
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_ml(:), tmp_ml2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        ! Workspace here is (pm+2)*(pm+1)
        call fmm_m2l_ztranslate_work(c(3), src_r, dst_r, pm, pl, vscales, &
            & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
        return
    end if
    ! Prepare pointers
    p = max(pm, pl)
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_ml(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_ml2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(pm, vcos, vsin, alpha, src_m, zero, tmp_ml)
    ! Perform rotation in the OXZ plane, work size is 4*pm*pm+13*pm+4
    call fmm_sph_rotate_oxz_work(pm, ctheta, -stheta, one, tmp_ml, zero, &
        & tmp_ml2, work)
    ! OZ translation, workspace here is (pm+2)*(pm+1)
    call fmm_m2l_ztranslate_work(rho, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, one, tmp_ml2, zero, tmp_ml, work)
    ! Backward rotation in the OXZ plane, work size is 4*pl*pl+13*pl+4
    call fmm_sph_rotate_oxz_work(pl, ctheta, stheta, one, tmp_ml, zero, &
        & tmp_ml2, work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(pl, vcos, vsin, one, tmp_ml2, beta, dst_l)
end subroutine fmm_m2l_rotation_work

!> Adjoint M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of a
!! multipole-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] m2l_ztranslate_adj_oef:
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: expansion in new harmonics
subroutine fmm_m2l_rotation_adj(c, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), alpha, &
        & src_l((pl+1)*(pl+1)), beta
    integer, intent(in) :: pl, pm
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Temporary workspace
    real(dp) :: work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8)
    ! Call corresponding work routine
    call fmm_m2l_rotation_adj_work(c, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
end subroutine fmm_m2l_rotation_adj

!> Adjoint M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of a
!! multipole-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] m2l_ztranslate_adj_coef:
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8 where p is a
!!      maximum of pm and pl
subroutine fmm_m2l_rotation_adj_work(c, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), alpha, &
        & src_l((pl+1)*(pl+1)), beta
    integer, intent(in) :: pl, pm
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Temporary workspace
    real(dp), intent(out), target :: &
        & work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n, p
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_ml(:), tmp_ml2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. 0) then
        ! Workspace here is (pl+2)*(pl+1)
        call fmm_m2l_ztranslate_adj_work(c(3), src_r, dst_r, pl, pm, vscales, &
            & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
        return
    end if
    ! Prepare pointers
    p = max(pm, pl)
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_ml(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_ml2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(pl, vcos, vsin, alpha, src_l, zero, tmp_ml)
    ! Perform rotation in the OXZ plane, work size is 4*pl*pl+13*pl+4
    call fmm_sph_rotate_oxz_work(pl, ctheta, -stheta, one, tmp_ml, zero, &
        & tmp_ml2, work)
    ! OZ translation, workspace here is (pl+2)*(pl+1)
    call fmm_m2l_ztranslate_adj_work(rho, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, one, tmp_ml2, zero, tmp_ml, work)
    ! Backward rotation in the OXZ plane, work size is 4*pm*pm+13*pm+4
    call fmm_sph_rotate_oxz_work(pm, ctheta, stheta, one, tmp_ml, zero, &
        & tmp_ml2, work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(pm, vcos, vsin, one, tmp_ml2, beta, dst_m)
end subroutine fmm_m2l_rotation_adj_work

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
        & cnode, rnode, snode)
    ! Inputs
    integer, intent(in) :: nsph
    real(dp), intent(in) :: csph(3, nsph), rsph(nsph)
    ! Outputs
    integer, intent(out) :: order(nsph), cluster(2, 2*nsph-1), &
        & children(2, 2*nsph-1), parent(2*nsph-1), snode(nsph)
    real(dp), intent(out) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
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
            call tree_rib_node_bisect(nsph, csph, n, order(s:e), div)
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
subroutine tree_rib_node_bisect(nsph, csph, n, order, div)
    ! Inputs
    integer, intent(in) :: nsph, n
    real(dp), intent(in) :: csph(3, nsph)
    ! Outputs
    integer, intent(inout) :: order(n)
    integer, intent(out) :: div
    ! Local variables
    real(dp) :: c(3), tmp_csph(3, n), s(3)
    real(dp), allocatable :: work(:)
    external :: dgesvd
    integer :: i, l, r, lwork, info, tmp_order(n), istat
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
    lwork = s(1)
    allocate(work(lwork), stat=istat)
    if (istat .ne. 0) stop "allocation failed"
    ! Get right singular vectors
    call dgesvd('N', 'O', 3, n, tmp_csph, 3, s, tmp_csph, 3, tmp_csph, 3, &
        & work, lwork, info)
    if (info .ne. 0) stop "dgesvd did not converge"
    deallocate(work, stat=istat)
    if (istat .ne. 0) stop "deallocation failed"
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
            !write(*,*) "FAR:", j
            work(3, iwork) = 1
        else if (npairs .eq. 1) then
            ! Mark as near admissible pair if both nodes are leaves
            !write(*,*) "NEAR:", j
            work(3, iwork) = 2
        else if (jwork+npairs .gt. lwork) then
            ! Exit procedure, since work array was too small
            !write(*,*) "SMALL LWORK"
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
            !write(*,*) "NON:", j
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

! Get near and far admissible pairs from work array of tree_get_farnear_work
! Works only for binary tree
subroutine tree_get_farnear(jwork, lwork, work, n, nnfar, nfar, sfar, far, &
        & nnnear, nnear, snear, near)
! Parameters:
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
    integer :: cfar(n+1), cnear(n+1)
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
end subroutine tree_get_farnear

!> Transfer multipole coefficients over a tree
subroutine tree_m2m_rotation(ddx_data, node_m)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*ddx_data % pm**2 + 19*ddx_data % pm + 8)
    ! Call corresponding work routine
    call tree_m2m_rotation_work(ddx_data, node_m, work)
end subroutine tree_m2m_rotation

!> Transfer multipole coefficients over a tree
subroutine tree_m2m_rotation_work(ddx_data, node_m, work)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Temporary workspace
    real(dp), intent(out) :: work(6*ddx_data % pm**2 + 19*ddx_data % pm + 8)
    ! Local variables
    integer :: i, j
    real(dp) :: c1(3), c(3), r1, r
    ! Bottom-to-top pass
    do i = ddx_data % nclusters, 1, -1
        ! Leaf node does not need any update
        if (ddx_data % children(1, i) == 0) cycle
        c = ddx_data % cnode(:, i)
        r = ddx_data % rnode(i)
        ! First child initializes output
        j = ddx_data % children(1, i)
        c1 = ddx_data % cnode(:, j)
        r1 = ddx_data % rnode(j)
        call fmm_m2m_rotation_work(c1-c, r1, r, &
            & ddx_data % pm, &
            & ddx_data % vscales, &
            & ddx_data % vcnk, one, &
            & node_m(:, j), zero, node_m(:, i), work)
        ! All other children update the same output
        do j = ddx_data % children(1, i)+1, ddx_data % children(2, i)
            c1 = ddx_data % cnode(:, j)
            r1 = ddx_data % rnode(j)
            call fmm_m2m_rotation_work(c1-c, r1, r, ddx_data % pm, &
                & ddx_data % vscales, ddx_data % vcnk, one, &
                & node_m(:, j), one, node_m(:, i), work)
        end do
    end do
end subroutine tree_m2m_rotation_work

!> Adjoint transfer multipole coefficients over a tree
subroutine tree_m2m_rotation_adj(ddx_data, node_m)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*ddx_data % pm**2 + 19*ddx_data % pm + 8)
    ! Call corresponding work routine
    call tree_m2m_rotation_adj_work(ddx_data, node_m, work)
end subroutine tree_m2m_rotation_adj

!> Adjoint transfer multipole coefficients over a tree
subroutine tree_m2m_rotation_adj_work(ddx_data, node_m, work)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Temporary workspace
    real(dp), intent(out) :: work(6*ddx_data % pm**2 + 19*ddx_data % pm + 8)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Top-to-bottom pass
    do i = 2, ddx_data % nclusters
        j = ddx_data % parent(i)
        c = ddx_data % cnode(:, j)
        r = ddx_data % rnode(j)
        c1 = ddx_data % cnode(:, i)
        r1 = ddx_data % rnode(i)
        call fmm_m2m_rotation_adj_work(c-c1, r, r1, ddx_data % pm, &
            & ddx_data % vscales, ddx_data % vcnk, one, node_m(:, j), one, &
            & node_m(:, i), work)
    end do
end subroutine tree_m2m_rotation_adj_work

!> Save M2M matrices for entire tree
subroutine tree_m2m_reflection_get_mat(ddx_data)
    ! Input/output
    type(ddx_type), intent(inout) :: ddx_data
    ! Local variables
    integer :: i, j, istatus
    real(dp) :: c1(3), c(3), r1, r
    ! For each non-root node define m2m matrices
    do i = 2, ddx_data % nclusters
        j = ddx_data % parent(i)
        c = ddx_data % cnode(:, j)
        r = ddx_data % rnode(j)
        c1 = ddx_data % cnode(:, i)
        r1 = ddx_data % rnode(i)
        call fmm_m2m_reflection_get_mat(c1-c, r1, r, ddx_data % pm, &
            & ddx_data % vscales, ddx_data % vfact, &
            & ddx_data % m2m_reflect_mat(:, i-1), &
            & ddx_data % m2m_ztranslate_mat(:, i-1))
    end do
end subroutine tree_m2m_reflection_get_mat

!> Transfer multipole coefficients over a tree
subroutine tree_m2m_reflection_use_mat(ddx_data, node_m)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Local variables
    integer :: i, j
    real(dp) :: c1(3), c(3), r1, r
    ! Bottom-to-top pass
    do i = ddx_data % nclusters, 1, -1
        ! Leaf node does not need any update
        if (ddx_data % children(1, i) == 0) cycle
        c = ddx_data % cnode(:, i)
        r = ddx_data % rnode(i)
        ! First child initializes output
        j = ddx_data % children(1, i)
        c1 = ddx_data % cnode(:, j)
        r1 = ddx_data % rnode(j)
        call fmm_m2m_reflection_use_mat(c1-c, r1, r, ddx_data % pm, &
            & ddx_data % m2m_reflect_mat(:, j-1), &
            & ddx_data % m2m_ztranslate_mat(:, j-1), one, &
            & node_m(:, j), zero, node_m(:, i))
        ! All other children update the same output
        do j = ddx_data % children(1, i)+1, ddx_data % children(2, i)
            c1 = ddx_data % cnode(:, j)
            r1 = ddx_data % rnode(j)
            call fmm_m2m_reflection_use_mat(c1-c, r1, r, ddx_data % pm, &
                & ddx_data % m2m_reflect_mat(:, j-1), &
                & ddx_data % m2m_ztranslate_mat(:, j-1), one, &
                & node_m(:, j), one, node_m(:, i))
        end do
    end do
end subroutine tree_m2m_reflection_use_mat

!> Transfer multipole coefficients over a tree
subroutine tree_m2m_reflection_use_mat_adj(ddx_data, node_m)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Local variables
    integer :: i, j
    real(dp) :: c1(3), c(3), r1, r
    ! Top-to-bottom pass
    do i = 2, ddx_data % nclusters
        j = ddx_data % parent(i)
        c = ddx_data % cnode(:, j)
        r = ddx_data % rnode(j)
        c1 = ddx_data % cnode(:, i)
        r1 = ddx_data % rnode(i)
        call fmm_m2m_reflection_use_mat_adj(c-c1, r, r1, ddx_data % pm, &
            & ddx_data % m2m_reflect_mat(:, i-1), &
            & ddx_data % m2m_ztranslate_mat(:, i-1), &
            & one, node_m(:, j), one, node_m(:, i))
    end do
end subroutine tree_m2m_reflection_use_mat_adj

!> Transfer local coefficients over a tree
subroutine tree_l2l_rotation(ddx_data, node_l)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*ddx_data % pl**2 + 19*ddx_data % pl + 8)
    ! Call corresponding work routine
    call tree_l2l_rotation_work(ddx_data, node_l, work)
end subroutine tree_l2l_rotation

!> Transfer local coefficients over a tree
subroutine tree_l2l_rotation_work(ddx_data, node_l, work)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Temporary workspace
    real(dp), intent(out) :: work(6*ddx_data % pl**2 + 19*ddx_data % pl + 8)
    ! Local variables
    integer :: i, j
    real(dp) :: c1(3), c(3), r1, r
    ! Top-to-bottom pass
    do i = 2, ddx_data % nclusters
        j = ddx_data % parent(i)
        c = ddx_data % cnode(:, j)
        r = ddx_data % rnode(j)
        c1 = ddx_data % cnode(:, i)
        r1 = ddx_data % rnode(i)
        call fmm_l2l_rotation_work(c-c1, r, r1, ddx_data % pl, &
            & ddx_data % vscales, ddx_data % vfact, one, &
            & node_l(:, j), one, node_l(:, i), work)
    end do
end subroutine tree_l2l_rotation_work

!> Adjoint transfer local coefficients over a tree
subroutine tree_l2l_rotation_adj(ddx_data, node_l)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*ddx_data % pl**2 + 19*ddx_data % pl + 8)
    ! Call corresponding work routine
    call tree_l2l_rotation_adj_work(ddx_data, node_l, work)
end subroutine tree_l2l_rotation_adj

!> Adjoint transfer local coefficients over a tree
subroutine tree_l2l_rotation_adj_work(ddx_data, node_l, work)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Temporary workspace
    real(dp), intent(out) :: work(6*ddx_data % pl**2 + 19*ddx_data % pl + 8)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Bottom-to-top pass
    do i = ddx_data % nclusters, 1, -1
        ! Leaf node does not need any update
        if (ddx_data % children(1, i) == 0) cycle
        c = ddx_data % cnode(:, i)
        r = ddx_data % rnode(i)
        ! First child initializes output
        j = ddx_data % children(1, i)
        c1 = ddx_data % cnode(:, j)
        r1 = ddx_data % rnode(j)
        call fmm_l2l_rotation_adj_work(c1-c, r1, r, ddx_data % pl, &
            & ddx_data % vscales, ddx_data % vfact, one, &
            & node_l(:, j), zero, node_l(:, i), work)
        ! All other children update the same output
        do j = ddx_data % children(1, i)+1, ddx_data % children(2, i)
            c1 = ddx_data % cnode(:, j)
            r1 = ddx_data % rnode(j)
            call fmm_l2l_rotation_adj_work(c1-c, r1, r, ddx_data % pl, &
                & ddx_data % vscales, ddx_data % vfact, one, &
                & node_l(:, j), one, node_l(:, i), work)
        end do
    end do
end subroutine tree_l2l_rotation_adj_work

!> Save L2L matrices for entire tree
subroutine tree_l2l_reflection_get_mat(ddx_data)
    ! Input/output
    type(ddx_type), intent(inout) :: ddx_data
    ! Local variables
    integer :: i, j, istatus
    real(dp) :: c1(3), c(3), r1, r
    ! For each non-root node define m2m matrices
    do i = 2, ddx_data % nclusters
        j = ddx_data % parent(i)
        c = ddx_data % cnode(:, j)
        r = ddx_data % rnode(j)
        c1 = ddx_data % cnode(:, i)
        r1 = ddx_data % rnode(i)
        call fmm_l2l_reflection_get_mat(c-c1, r, r1, ddx_data % pl, &
            & ddx_data % vscales, ddx_data % vfact, &
            & ddx_data % l2l_reflect_mat(:, i-1), &
            & ddx_data % l2l_ztranslate_mat(:, i-1))
    end do
end subroutine tree_l2l_reflection_get_mat

!> Transfer local coefficients over a tree
subroutine tree_l2l_reflection_use_mat(ddx_data, node_l)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Local variables
    integer :: i, j
    real(dp) :: c1(3), c(3), r1, r
    ! Top-to-bottom pass
    do i = 2, ddx_data % nclusters
        j = ddx_data % parent(i)
        c = ddx_data % cnode(:, j)
        r = ddx_data % rnode(j)
        c1 = ddx_data % cnode(:, i)
        r1 = ddx_data % rnode(i)
        call fmm_l2l_reflection_use_mat(c-c1, r, r1, ddx_data % pl, &
            & ddx_data % l2l_reflect_mat(:, i-1), &
            & ddx_data % l2l_ztranslate_mat(:, i-1), one, &
            & node_l(:, j), one, node_l(:, i))
    end do
end subroutine tree_l2l_reflection_use_mat

!> Adjoint Transfer local coefficients over a tree
subroutine tree_l2l_reflection_use_mat_adj(ddx_data, node_l)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Local variables
    integer :: i, j
    real(dp) :: c1(3), c(3), r1, r
    ! Bottom-to-top pass
    do i = ddx_data % nclusters, 1, -1
        ! Leaf node does not need any update
        if (ddx_data % children(1, i) == 0) cycle
        c = ddx_data % cnode(:, i)
        r = ddx_data % rnode(i)
        ! First child initializes output
        j = ddx_data % children(1, i)
        c1 = ddx_data % cnode(:, j)
        r1 = ddx_data % rnode(j)
        call fmm_l2l_reflection_use_mat_adj(c1-c, r1, r, ddx_data % pl, &
            & ddx_data % l2l_reflect_mat(:, j-1), &
            & ddx_data % l2l_ztranslate_mat(:, j-1), &
            & one, node_l(:, j), zero, node_l(:, i))
        ! All other children update the same output
        do j = ddx_data % children(1, i)+1, ddx_data % children(2, i)
            c1 = ddx_data % cnode(:, j)
            r1 = ddx_data % rnode(j)
            call fmm_l2l_reflection_use_mat_adj(c1-c, r1, r, ddx_data % pl, &
                & ddx_data % l2l_reflect_mat(:, j-1), &
                & ddx_data % l2l_ztranslate_mat(:, j-1), &
                & one, node_l(:, j), one, node_l(:, i))
        end do
    end do
end subroutine tree_l2l_reflection_use_mat_adj

!> Transfer multipole local coefficients into local over a tree
subroutine tree_m2l_rotation(ddx_data, node_m, node_l)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Output
    real(dp), intent(out) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*max(ddx_data % pm, ddx_data % pl)**2 + &
        & 19*max(ddx_data % pm, ddx_data % pl) + 8)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Any order of this cycle is OK
    do i = 1, ddx_data % nclusters
        ! If no far admissible pairs just set output to zero
        if (ddx_data % nfar(i) .eq. 0) then
            node_l(:, i) = zero
            cycle
        end if
        c = ddx_data % cnode(:, i)
        r = ddx_data % rnode(i)
        ! Use the first far admissible pair to initialize output
        k = ddx_data % far(ddx_data % sfar(i))
        c1 = ddx_data % cnode(:, k)
        r1 = ddx_data % rnode(k)
        call fmm_m2l_rotation_work(c1-c, r1, r, ddx_data % pm, ddx_data % pl, &
            & ddx_data % vscales, ddx_data % m2l_ztranslate_coef, one, &
            & node_m(:, k), zero, node_l(:, i), work)
        do j = ddx_data % sfar(i)+1, ddx_data % sfar(i+1)-1
            k = ddx_data % far(j)
            c1 = ddx_data % cnode(:, k)
            r1 = ddx_data % rnode(k)
            call fmm_m2l_rotation_work(c1-c, r1, r, ddx_data % pm, &
                & ddx_data % pl, ddx_data % vscales, &
                & ddx_data % m2l_ztranslate_coef, one, node_m(:, k), one, &
                & node_l(:, i), work)
        end do
    end do
end subroutine tree_m2l_rotation

!> Adjoint transfer multipole local coefficients into local over a tree
subroutine tree_m2l_rotation_adj(ddx_data, node_l, node_m)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Output
    real(dp), intent(out) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Any order of this cycle is OK
    node_m = zero
    do i = 1, ddx_data % nclusters
        ! If no far admissible pairs just set output to zero
        if (ddx_data % nfar(i) .eq. 0) then
            cycle
        end if
        c = ddx_data % cnode(:, i)
        r = ddx_data % rnode(i)
        ! Use the first far admissible pair to initialize output
        k = ddx_data % far(ddx_data % sfar(i))
        c1 = ddx_data % cnode(:, k)
        r1 = ddx_data % rnode(k)
        call fmm_m2l_rotation_adj(c-c1, r, r1, ddx_data % pl, ddx_data % pm, &
            & ddx_data % vscales, ddx_data % m2l_ztranslate_adj_coef, one, &
            & node_l(:, i), one, node_m(:, k))
        do j = ddx_data % sfar(i)+1, ddx_data % sfar(i+1)-1
            k = ddx_data % far(j)
            c1 = ddx_data % cnode(:, k)
            r1 = ddx_data % rnode(k)
            call fmm_m2l_rotation_adj(c-c1, r, r1, ddx_data % pl, ddx_data % pm, &
                & ddx_data % vscales, ddx_data % m2l_ztranslate_adj_coef, one, &
                & node_l(:, i), one, node_m(:, k))
        end do
    end do
end subroutine tree_m2l_rotation_adj

!> Precompute M2L translations for all nodes
subroutine tree_m2l_reflection_get_mat(ddx_data)
    ! Input/output
    type(ddx_type), intent(inout) :: ddx_data
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Any order of this cycle is OK
    do i = 1, ddx_data % nclusters
        c = ddx_data % cnode(:, i)
        r = ddx_data % rnode(i)
        do j = ddx_data % sfar(i), ddx_data % sfar(i+1)-1
            k = ddx_data % far(j)
            c1 = ddx_data % cnode(:, k)
            r1 = ddx_data % rnode(k)
            call fmm_m2l_reflection_get_mat(c1-c, r1, r, ddx_data % pm, &
                & ddx_data % pl, ddx_data % vscales, ddx_data % vfact, &
                & ddx_data % m2l_reflect_mat(:, j), &
                & ddx_data % m2l_ztranslate_mat(:, j))
        end do
    end do
end subroutine tree_m2l_reflection_get_mat

!> Transfer multipole local coefficients into local over a tree
subroutine tree_m2l_reflection_use_mat(ddx_data, node_m, node_l)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Output
    real(dp), intent(out) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Any order of this cycle is OK
    do i = 1, ddx_data % nclusters
        ! If no far admissible pairs just set output to zero
        if (ddx_data % nfar(i) .eq. 0) then
            node_l(:, i) = zero
            cycle
        end if
        c = ddx_data % cnode(:, i)
        r = ddx_data % rnode(i)
        ! Use the first far admissible pair to initialize output
        j = ddx_data % sfar(i)
        k = ddx_data % far(j)
        c1 = ddx_data % cnode(:, k)
        r1 = ddx_data % rnode(k)
        call fmm_m2l_reflection_use_mat(c1-c, r1, r, ddx_data % pm, &
            & ddx_data % pl, ddx_data % m2l_reflect_mat(:, j), &
            & ddx_data % m2l_ztranslate_mat(:, j), one, node_m(:, k), zero, &
            & node_l(:, i))
        do j = ddx_data % sfar(i)+1, ddx_data % sfar(i+1)-1
            k = ddx_data % far(j)
            c1 = ddx_data % cnode(:, k)
            r1 = ddx_data % rnode(k)
            call fmm_m2l_reflection_use_mat(c1-c, r1, r, ddx_data % pm, &
                & ddx_data % pl, ddx_data % m2l_reflect_mat(:, j), &
                & ddx_data % m2l_ztranslate_mat(:, j), one, node_m(:, k), one, &
                & node_l(:, i))
        end do
    end do
end subroutine tree_m2l_reflection_use_mat

!> Transfer multipole local coefficients into local over a tree
subroutine tree_m2l_reflection_use_mat_adj(ddx_data, node_l, node_m)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Output
    real(dp), intent(out) :: node_m((ddx_data % pm+1)**2, ddx_data % nclusters)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Any order of this cycle is OK
    node_m = zero
    do i = 1, ddx_data % nclusters
        ! If no far admissible pairs just set output to zero
        if (ddx_data % nfar(i) .eq. 0) then
            cycle
        end if
        c = ddx_data % cnode(:, i)
        r = ddx_data % rnode(i)
        ! Use the first far admissible pair to initialize output
        j = ddx_data % sfar(i)
        k = ddx_data % far(j)
        c1 = ddx_data % cnode(:, k)
        r1 = ddx_data % rnode(k)
        call fmm_m2l_reflection_use_mat_adj(c-c1, r, r1, ddx_data % pl, &
            & ddx_data % pm, ddx_data % m2l_reflect_mat(:, j), &
            & ddx_data % m2l_ztranslate_mat(:, j), one, node_l(:, i), one, &
            & node_m(:, k))
        do j = ddx_data % sfar(i)+1, ddx_data % sfar(i+1)-1
            k = ddx_data % far(j)
            c1 = ddx_data % cnode(:, k)
            r1 = ddx_data % rnode(k)
            call fmm_m2l_reflection_use_mat_adj(c-c1, r, r1, ddx_data % pl, &
                & ddx_data % pm, ddx_data % m2l_reflect_mat(:, j), &
                & ddx_data % m2l_ztranslate_mat(:, j), one, node_l(:, i), one, &
                & node_m(:, k))
        end do
    end do
end subroutine tree_m2l_reflection_use_mat_adj

subroutine tree_l2p(ddx_data, alpha, node_l, beta, grid_v)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters), &
        & alpha, beta
    ! Output
    real(dp), intent(inout) :: grid_v(ddx_data % ngrid, ddx_data % nsph)
    ! Local variables
    real(dp) :: sph_l((ddx_data % pl+1)**2, ddx_data % nsph), c(3)
    integer :: isph
    external :: dgemm
    ! Init output
    if (beta .eq. zero) then
        grid_v = zero
    else
        grid_v = beta * grid_v
    end if
    ! Get data from all clusters to spheres
    do isph = 1, ddx_data % nsph
        sph_l(:, isph) = node_l(:, ddx_data % snode(isph))
    end do
    ! Get values at grid points
    call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
        & (ddx_data % pl+1)**2, alpha, ddx_data % l2grid, &
        & ddx_data % vgrid_nbasis, sph_l, (ddx_data % pl+1)**2, beta, grid_v, &
        & ddx_data % ngrid)
end subroutine tree_l2p

subroutine tree_l2p_adj(ddx_data, alpha, grid_v, beta, node_l)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: grid_v(ddx_data % ngrid, ddx_data % nsph), alpha, &
        & beta
    ! Output
    real(dp), intent(inout) :: node_l((ddx_data % pl+1)**2, &
        & ddx_data % nclusters)
    ! Local variables
    real(dp) :: sph_l((ddx_data % pl+1)**2, ddx_data % nsph), c(3)
    integer :: isph, inode
    external :: dgemm
    ! Init output
    if (beta .eq. zero) then
        node_l = zero
    else
        node_l = beta * node_l
    end if
    ! Get weights of spherical harmonics at each sphere
    call dgemm('N', 'N', (ddx_data % pl+1)**2, ddx_data % nsph, &
        & ddx_data % ngrid, one, ddx_data % l2grid, ddx_data % vgrid_nbasis, &
        & grid_v, ddx_data % ngrid, zero, sph_l, &
        & (ddx_data % pl+1)**2)
    ! Get data from all clusters to spheres
    do isph = 1, ddx_data % nsph
        inode = ddx_data % snode(isph)
        node_l(:, inode) = node_l(:, inode) + alpha*sph_l(:, isph)
    end do
end subroutine tree_l2p_adj

subroutine tree_m2p(ddx_data, p, alpha, sph_m, beta, grid_v)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    integer, intent(in) :: p
    real(dp), intent(in) :: sph_m((p+1)**2, ddx_data % nsph), alpha, beta
    ! Output
    real(dp), intent(inout) :: grid_v(ddx_data % ngrid, ddx_data % nsph)
    ! Local variables
    integer :: isph, inode, jnear, jnode, jsph, igrid
    real(dp) :: c(3)
    ! Temporary workspace
    real(dp) :: work((p+1)*(p+1)+3*p)
    ! Init output
    if (beta .eq. zero) then
        grid_v = zero
    else
        grid_v = beta * grid_v
    end if
    ! Cycle over all spheres
    do isph = 1, ddx_data % nsph
        ! Cycle over all near-field admissible pairs of spheres
        inode = ddx_data % snode(isph)
        do jnear = ddx_data % snear(inode), ddx_data % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = ddx_data % near(jnear)
            jsph = ddx_data % order(ddx_data % cluster(1, jnode))
            ! Ignore self-interaction
            if(isph .eq. jsph) cycle
            ! Accumulate interaction for external grid points only
            do igrid = 1, ddx_data % ngrid
                if(ddx_data % ui(igrid, isph) .eq. zero) cycle
                c = ddx_data % cgrid(:, igrid)*ddx_data % rsph(isph) - &
                    & ddx_data % csph(:, jsph) + ddx_data % csph(:, isph)
                call fmm_m2p_work(c, ddx_data % rsph(jsph), p, &
                    & ddx_data % vscales_rel, alpha, sph_m(:, jsph), one, &
                    & grid_v(igrid, isph), work)
            end do
        end do
    end do
end subroutine tree_m2p

subroutine tree_m2p_adj(ddx_data, p, alpha, grid_v, beta, sph_m)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    integer, intent(in) :: p
    real(dp), intent(in) :: grid_v(ddx_data % ngrid, ddx_data % nsph), alpha, &
        & beta
    ! Output
    real(dp), intent(inout) :: sph_m((p+1)**2, ddx_data % nsph)
    ! Temporary workspace
    real(dp) :: work((p+1)*(p+1)+3*p)
    ! Local variables
    integer :: isph, inode, jnear, jnode, jsph, igrid
    real(dp) :: c(3)
    ! Init output
    if (beta .eq. zero) then
        sph_m = zero
    else
        sph_m = beta * sph_m
    end if
    ! Cycle over all spheres
    do isph = 1, ddx_data % nsph
        ! Cycle over all near-field admissible pairs of spheres
        inode = ddx_data % snode(isph)
        do jnear = ddx_data % snear(inode), ddx_data % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = ddx_data % near(jnear)
            jsph = ddx_data % order(ddx_data % cluster(1, jnode))
            ! Ignore self-interaction
            if(isph .eq. jsph) cycle
            ! Accumulate interaction for external grid points only
            do igrid = 1, ddx_data % ngrid
                if(ddx_data % ui(igrid, isph) .eq. zero) cycle
                c = ddx_data % cgrid(:, igrid)*ddx_data % rsph(isph) - &
                    & ddx_data % csph(:, jsph) + ddx_data % csph(:, isph)
                call fmm_m2p_adj_work(c, alpha*grid_v(igrid, isph), &
                    & ddx_data % rsph(jsph), p, ddx_data % vscales_rel, one, &
                    & sph_m(:, jsph), work)
            end do
        end do
    end do
end subroutine tree_m2p_adj

subroutine tree_m2p_get_mat(ddx_data)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    ! Local variables
    integer :: isph, inode, jnear, jnode, jsph, i, j, icav, igrid, m2p_column
    real(dp) :: c(3)
    ! Cycle over all spheres
    m2p_column = 1
    do isph = 1, ddx_data % nsph
        ! Cycle over all near-field admissible pairs of spheres
        inode = ddx_data % snode(isph)
        do jnear = ddx_data % snear(inode), ddx_data % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = ddx_data % near(jnear)
            jsph = ddx_data % order(ddx_data % cluster(1, jnode))
            ! Ignore self-interaction
            if (isph .eq. jsph) cycle
            ! M2P marices from spherical harmonics of jsph sphere to all
            ! external grid points of isph sphere
            do icav = ddx_data % icav_ia(isph), ddx_data % icav_ia(isph+1)-1
                igrid = ddx_data % icav_ja(icav)
                c = ddx_data % cgrid(:, igrid)*ddx_data % rsph(isph) - &
                    & ddx_data % csph(:, jsph) + ddx_data % csph(:, isph)
                call fmm_m2p_mat(c, ddx_data % rsph(jsph), ddx_data % m2p_lmax, &
                    & ddx_data % vscales, ddx_data % m2p_mat(:, m2p_column))
                m2p_column = m2p_column + 1
            end do
        end do
    end do
end subroutine tree_m2p_get_mat

subroutine tree_m2p_use_mat(ddx_data, p, alpha, sph_m, beta, grid_v)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    integer, intent(in) :: p
    real(dp), intent(in) :: sph_m((p+1)**2, ddx_data % nsph), alpha, beta
    ! Output
    real(dp), intent(inout) :: grid_v(ddx_data % ngrid, ddx_data % nsph)
    ! Local variables
    integer :: isph, inode, jnear, jnode, jsph, icav_sph, igrid, m2p_column
    real(dp) :: c(3), x(ddx_data % ngrid)
    ! Init output
    if (beta .eq. zero) then
        grid_v = zero
    else
        grid_v = beta * grid_v
    end if
    ! Cycle over all spheres
    m2p_column = 1
    do isph = 1, ddx_data % nsph
        x = zero
        ! Cycle over all near-field admissible pairs of spheres
        inode = ddx_data % snode(isph)
        do jnear = ddx_data % snear(inode), ddx_data % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = ddx_data % near(jnear)
            jsph = ddx_data % order(ddx_data % cluster(1, jnode))
            ! Ignore self-interaction
            if(isph .eq. jsph) cycle
            ! Accumulate interaction for external grid points only
            call dgemv('T', (p+1)**2, ddx_data % ncav_sph(isph), one, &
                & ddx_data % m2p_mat(1, m2p_column), ddx_data % m2p_nbasis, &
                & sph_m(1, jsph), 1, one, x(1), 1)
            m2p_column = m2p_column + ddx_data % ncav_sph(isph)
        end do
        ! Move sparsely stored x into full grid_v
        do icav_sph = 1, ddx_data % ncav_sph(isph)
            igrid = ddx_data % icav_ja(ddx_data % icav_ia(isph)+icav_sph-1)
            grid_v(igrid, isph) = grid_v(igrid, isph) + alpha*x(icav_sph)
        end do
    end do
end subroutine tree_m2p_use_mat

subroutine tree_m2p_use_mat_adj(ddx_data, p, alpha, grid_v, beta, sph_m)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    integer, intent(in) :: p
    real(dp), intent(in) :: grid_v(ddx_data % ngrid, ddx_data % nsph), alpha, &
        & beta
    ! Output
    real(dp), intent(inout) :: sph_m((p+1)**2, ddx_data % nsph)
    ! Local variables
    integer :: isph, inode, jnear, jnode, jsph, icav_sph, igrid, m2p_column
    real(dp) :: c(3), x(ddx_data % ngrid)
    ! Init output
    if (beta .eq. zero) then
        sph_m = zero
    else
        sph_m = beta * sph_m
    end if
    ! Cycle over all spheres
    m2p_column = 1
    do isph = 1, ddx_data % nsph
        ! Move full grid_v into sparsely stored x
        do icav_sph = 1, ddx_data % ncav_sph(isph)
            igrid = ddx_data % icav_ja(ddx_data % icav_ia(isph)+icav_sph-1)
            x(icav_sph) = alpha * grid_v(igrid, isph)
        end do
        ! Cycle over all near-field admissible pairs of spheres
        inode = ddx_data % snode(isph)
        do jnear = ddx_data % snear(inode), ddx_data % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = ddx_data % near(jnear)
            jsph = ddx_data % order(ddx_data % cluster(1, jnode))
            ! Ignore self-interaction
            if(isph .eq. jsph) cycle
            ! Accumulate interaction for external grid points only
            call dgemv('N', (p+1)**2, ddx_data % ncav_sph(isph), one, &
                & ddx_data % m2p_mat(1, m2p_column), ddx_data % m2p_nbasis, &
                & x(1), 1, one, sph_m(1, jsph), 1)
            m2p_column = m2p_column + ddx_data % ncav_sph(isph)
        end do
    end do
end subroutine tree_m2p_use_mat_adj

subroutine fdoka(ddx_data, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx )
    type(ddx_type), intent(in) :: ddx_data
      integer,                         intent(in)    :: isph
      real(dp),  dimension(ddx_data % nbasis, ddx_data % nsph), intent(in)    :: sigma
      real(dp),  dimension(ddx_data % ngrid),       intent(in)    :: xi
      real(dp),  dimension(ddx_data % nbasis),      intent(inout) :: basloc, vplm
      real(dp),  dimension(3, ddx_data % nbasis),    intent(inout) :: dbsloc
      real(dp),  dimension(ddx_data % lmax+1),      intent(inout) :: vcos, vsin
      real(dp),  dimension(3),           intent(inout) :: fx
!
      integer :: ig, ij, jsph, l, ind, m
      real(dp)  :: vvij, tij, xij, oij, t, fac, fl, f1, f2, f3, beta, tlow, thigh
      real(dp)  :: vij(3), sij(3), alp(3), va(3)
      real(dp), external :: dnrm2
!      
!-----------------------------------------------------------------------------------
!    
      tlow  = one - pt5*(one - ddx_data % se)*ddx_data % eta
      thigh = one + pt5*(one + ddx_data % se)*ddx_data % eta
!    
      do ig = 1, ddx_data % ngrid
        va = zero
        do ij = ddx_data % inl(isph), ddx_data % inl(isph+1) - 1
          jsph = ddx_data % nl(ij)
          vij  = ddx_data % csph(:,isph) + &
              & ddx_data % rsph(isph)*ddx_data % cgrid(:,ig) - &
              & ddx_data % csph(:,jsph)
          !vvij = sqrt(dot_product(vij,vij))
          vvij = dnrm2(3, vij, 1)
          tij  = vvij/ddx_data % rsph(jsph)
!    
          if (tij.ge.thigh) cycle
!    
          sij  = vij/vvij
          !call dbasis(sij,basloc,dbsloc,vplm,vcos,vsin)
          call dbasis(ddx_data, sij, basloc, dbsloc, vplm, vcos, vsin)
          alp  = zero
          t    = one
          do l = 1, ddx_data % lmax
            ind = l*l + l + 1
            fl  = dble(l)
            fac = t/(ddx_data % vscales(ind)**2)
            do m = -l, l
              f2 = fac*sigma(ind+m,jsph)
              f1 = f2*fl*basloc(ind+m)
              alp(:) = alp(:) + f1*sij(:) + f2*dbsloc(:,ind+m)
            end do
            t = t*tij
          end do
          beta = intmlp(ddx_data, tij,sigma(:,jsph),basloc)
          xij = fsw(tij,ddx_data % se,ddx_data % eta)
          if (ddx_data % fi(ig,isph).gt.one) then
            oij = xij/ddx_data % fi(ig,isph)
            f2  = -oij/ddx_data % fi(ig,isph)
          else
            oij = xij
            f2  = zero
          end if
          f1 = oij/ddx_data % rsph(jsph)
          va(:) = va(:) + f1*alp(:) + beta*f2*ddx_data % zi(:,ig,isph)
          if (tij .gt. tlow) then
            f3 = beta*dfsw(tij,ddx_data % se,ddx_data % eta)/ddx_data % rsph(jsph)
            if (ddx_data % fi(ig,isph).gt.one) f3 = f3/ddx_data % fi(ig,isph)
            va(:) = va(:) + f3*sij(:)
          end if
        end do
        fx = fx - ddx_data % wgrid(ig)*xi(ig)*va(:)
      end do
!      
      return
!      
!      
end subroutine fdoka
!-----------------------------------------------------------------------------------
!
!      
!      
!      
!-----------------------------------------------------------------------------------
subroutine fdokb(ddx_data, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx )
    type(ddx_type), intent(in) :: ddx_data
      integer,                         intent(in)    :: isph
      real(dp),  dimension(ddx_data % nbasis, ddx_data % nsph), intent(in)    :: sigma
      real(dp),  dimension(ddx_data % ngrid, ddx_data % nsph),  intent(in)    :: xi
      real(dp),  dimension(ddx_data % nbasis),      intent(inout) :: basloc, vplm
      real(dp),  dimension(3, ddx_data % nbasis),    intent(inout) :: dbsloc
      real(dp),  dimension(ddx_data % lmax+1),      intent(inout) :: vcos, vsin
      real(dp),  dimension(3),           intent(inout) :: fx
!
      integer :: ig, ji, jsph, l, ind, m, jk, ksph
      logical :: proc
      real(dp)  :: vvji, tji, xji, oji, t, fac, fl, f1, f2, beta, di, tlow, thigh
      real(dp)  :: b, g1, g2, vvjk, tjk, f, xjk
      real(dp)  :: vji(3), sji(3), alp(3), vb(3), vjk(3), sjk(3), vc(3)
      real(dp) :: rho, ctheta, stheta, cphi, sphi
      real(dp), external :: dnrm2
!
!-----------------------------------------------------------------------------------
!
      tlow  = one - pt5*(one - ddx_data % se)*ddx_data % eta
      thigh = one + pt5*(one + ddx_data % se)*ddx_data % eta
!
      do ig = 1, ddx_data % ngrid
        vb = zero
        vc = zero
        do ji = ddx_data % inl(isph), ddx_data % inl(isph+1) - 1
          jsph = ddx_data % nl(ji)
          vji  = ddx_data % csph(:,jsph) + &
              & ddx_data % rsph(jsph)*ddx_data % cgrid(:,ig) - &
              & ddx_data % csph(:,isph)
          !vvji = sqrt(dot_product(vji,vji))
          vvji = dnrm2(3, vji, 1)
          tji  = vvji/ddx_data % rsph(isph)
!
          if (tji.gt.thigh) cycle
!
          sji  = vji/vvji
          !call dbasis(sji,basloc,dbsloc,vplm,vcos,vsin)
          call dbasis(ddx_data, sji, basloc, dbsloc, vplm, vcos, vsin)
!
          alp = zero
          t   = one
          do l = 1, ddx_data % lmax
            ind = l*l + l + 1
            fl  = dble(l)
            fac = t/(ddx_data % vscales(ind)**2)
            do m = -l, l
              f2 = fac*sigma(ind+m,isph)
              f1 = f2*fl*basloc(ind+m)
              alp = alp + f1*sji + f2*dbsloc(:,ind+m)
            end do
            t = t*tji
          end do
          xji = fsw(tji,ddx_data % se,ddx_data % eta)
          if (ddx_data % fi(ig,jsph).gt.one) then
            oji = xji/ddx_data % fi(ig,jsph)
          else
            oji = xji
          end if
          f1 = oji/ddx_data % rsph(isph)
          vb = vb + f1*alp*xi(ig,jsph)
          if (tji .gt. tlow) then
            beta = intmlp(ddx_data, tji, sigma(:,isph), basloc)
            if (ddx_data % fi(ig,jsph) .gt. one) then
              di  = one/ddx_data % fi(ig,jsph)
              fac = di*xji
              proc = .false.
              b    = zero
              do jk = ddx_data % inl(jsph), ddx_data % inl(jsph+1) - 1
                ksph = ddx_data % nl(jk)
                vjk  = ddx_data % csph(:,jsph) + &
                    & ddx_data % rsph(jsph)*ddx_data % cgrid(:,ig) - &
                    & ddx_data % csph(:,ksph)
                !vvjk = sqrt(dot_product(vjk,vjk))
                vvjk = dnrm2(3, vjk, 1)
                tjk  = vvjk/ddx_data % rsph(ksph)
                if (ksph.ne.isph) then
                  if (tjk .le. thigh) then
                    proc = .true.
                    sjk  = vjk/vvjk
                    !call ylmbas(sjk,basloc,vplm,vcos,vsin)
                    call ylmbas(sjk, rho, ctheta, stheta, cphi, sphi, &
                        & ddx_data % lmax, ddx_data % vscales, basloc, vplm, &
                        & vcos, vsin)
                    g1  = intmlp(ddx_data, tjk, sigma(:,ksph), basloc)
                    xjk = fsw(tjk, ddx_data % se, ddx_data % eta)
                    b   = b + g1*xjk
                  end if
                end if
              end do
              if (proc) then
                g1 = di*di*dfsw(tji,ddx_data % se,ddx_data % eta)/ddx_data % rsph(isph)
                g2 = g1*xi(ig,jsph)*b
                vc = vc + g2*sji
              end if
            else
              di  = one
              fac = zero
            end if
            f2 = (one-fac)*di*dfsw(tji,ddx_data % se,ddx_data % eta)/ddx_data % rsph(isph)
            vb = vb + f2*xi(ig,jsph)*beta*sji
          end if 
        end do
        fx = fx + ddx_data % wgrid(ig)*(vb - vc)
      end do
      return
  end subroutine fdokb
!-----------------------------------------------------------------------------------
!
!
!
!
!-----------------------------------------------------------------------------------
subroutine fdoga(ddx_data, isph, xi, phi, fx )
    type(ddx_type), intent(in) :: ddx_data
      integer,                        intent(in)    :: isph
      real(dp),  dimension(ddx_data % ngrid, ddx_data % nsph), intent(in)    :: xi, phi
      real(dp),  dimension(3),          intent(inout) :: fx
!
      integer :: ig, ji, jsph
      real(dp)  :: vvji, tji, fac, swthr
      real(dp)  :: alp(3), vji(3), sji(3)
      real(dp), external :: dnrm2
!
!-----------------------------------------------------------------------------------
!
      do ig = 1, ddx_data % ngrid
        alp = zero
        if (ddx_data % ui(ig,isph) .gt. zero .and. ddx_data % ui(ig,isph).lt.one) then
          alp = alp + phi(ig,isph)*xi(ig,isph)*ddx_data % zi(:,ig,isph)
        end if
        do ji = ddx_data % inl(isph), ddx_data % inl(isph+1) - 1
          jsph  = ddx_data % nl(ji)
          vji   = ddx_data % csph(:,jsph) + &
              & ddx_data % rsph(jsph)*ddx_data % cgrid(:,ig) - &
              & ddx_data % csph(:,isph)
          !vvji  = sqrt(dot_product(vji,vji))
          vvji = dnrm2(3, vji, 1)
          tji   = vvji/ddx_data % rsph(isph)
          swthr = one + (ddx_data % se + 1.d0)*ddx_data % eta / 2.d0
          if (tji.lt.swthr .and. tji.gt.swthr-ddx_data % eta .and. ddx_data % ui(ig,jsph).gt.zero) then
            sji = vji/vvji
            fac = - dfsw(tji,ddx_data % se,ddx_data % eta)/ddx_data % rsph(isph)
            alp = alp + fac*phi(ig,jsph)*xi(ig,jsph)*sji
          end if
        end do
        fx = fx - ddx_data % wgrid(ig)*alp
      end do
!
      return 
!
!
end subroutine fdoga

subroutine efld(nsrc,src,csrc,ntrg,ctrg,ef)
integer,                    intent(in)    :: nsrc, ntrg
real*8,  dimension(nsrc),   intent(in)    :: src
real*8,  dimension(3,nsrc), intent(in)    :: csrc
real*8,  dimension(3,ntrg), intent(in)    :: ctrg
real*8,  dimension(3,ntrg), intent(inout) :: ef
!
integer :: i, j
real*8  :: dx, dy, dz, r2, rr, r3, f, e(3)
real*8, parameter :: zero=0.0d0
!
ef = zero
do j = 1, ntrg
  e = zero
  do i = 1, nsrc
    dx   = ctrg(1,j) - csrc(1,i)
    dy   = ctrg(2,j) - csrc(2,i)
    dz   = ctrg(3,j) - csrc(3,i)
    r2   = dx*dx + dy*dy + dz*dz
    rr   = sqrt(r2)
    r3   = r2*rr
    f    = src(i)/r3
    e(1) = e(1) + f*dx
    e(2) = e(2) + f*dy
    e(3) = e(3) + f*dz
  end do
  ef(:,j) = e
end do
end subroutine efld

subroutine tree_grad_m2m(ddx_data, sph_m, sph_m_grad, work)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: sph_m(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: sph_m_grad((ddx_data % lmax+2)**2, 3, &
        & ddx_data % nsph)
    ! Temporary workspace
    real(dp), intent(out) :: work((ddx_data % lmax+2)**2, ddx_data % nsph)
    ! Local variables
    integer :: isph, l, indi, indj, m
    real(dp) :: tmp1, tmp2
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    ! Set coordinate transformations
    zx_coord_transform = zero
    zx_coord_transform(3, 2) = one
    zx_coord_transform(2, 3) = one
    zx_coord_transform(1, 1) = one
    zy_coord_transform = zero
    zy_coord_transform(1, 2) = one
    zy_coord_transform(2, 1) = one
    zy_coord_transform(3, 3) = one
    ! At first reflect harmonics of a degree up to lmax
    sph_m_grad(1:ddx_data % nbasis, 3, :) = sph_m
    do isph = 1, ddx_data % nsph
        call fmm_sph_transform(ddx_data % lmax, zx_coord_transform, one, &
            & sph_m(:, isph), zero, sph_m_grad(1:ddx_data % nbasis, 1, isph))
        call fmm_sph_transform(ddx_data % lmax, zy_coord_transform, one, &
            & sph_m(:, isph), zero, sph_m_grad(1:ddx_data % nbasis, 2, isph))
    end do
    ! Derivative of M2M translation over OZ axis at the origin consists of 2
    ! steps:
    !   1) increase degree l and scale by sqrt((2*l+1)*(l*l-m*m)) / sqrt(2*l-1)
    !   2) scale by 1/rsph(isph)
    do l = ddx_data % lmax+1, 1, -1
        indi = l*l + l + 1
        indj = indi - 2*l
        tmp1 = sqrt(dble(2*l+1)) / sqrt(dble(2*l-1))
        do m = 1-l, l-1
            tmp2 = sqrt(dble(l*l-m*m)) * tmp1
            sph_m_grad(indi+m, :, :) = tmp2 * sph_m_grad(indj+m, :, :)
        end do
        sph_m_grad(indi+l, :, :) = zero
        sph_m_grad(indi-l, :, :) = zero
    end do
    sph_m_grad(1, :, :) = zero
    ! Scale by 1/rsph(isph) and rotate harmonics of degree up to lmax+1 back to
    ! the initial axis. Coefficient of 0-th degree is zero so we ignore it.
    do isph = 1, ddx_data % nsph
        sph_m_grad(:, 3, isph) = sph_m_grad(:, 3, isph) / ddx_data % rsph(isph)
        work(:, isph) = sph_m_grad(:, 1, isph) / ddx_data % rsph(isph)
        call fmm_sph_transform(ddx_data % lmax+1, zx_coord_transform, one, &
            & work(:, isph), zero, sph_m_grad(:, 1, isph))
        work(:, isph) = sph_m_grad(:, 2, isph) / ddx_data % rsph(isph)
        call fmm_sph_transform(ddx_data % lmax+1, zy_coord_transform, one, &
            & work(:, isph), zero, sph_m_grad(:, 2, isph))
    end do
end subroutine tree_grad_m2m

subroutine tree_grad_l2l(ddx_data, node_l, sph_l_grad, work)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: node_l((ddx_data % pl+1)**2, ddx_data % nclusters)
    ! Output
    real(dp), intent(out) :: sph_l_grad((ddx_data % pl+1)**2, 3, ddx_data % nsph)
    ! Temporary workspace
    real(dp), intent(out) :: work((ddx_data % pl+1)**2, ddx_data % nsph)
    ! Local variables
    integer :: isph, inode, l, indi, indj, m
    real(dp) :: tmp1, tmp2
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    ! Gradient of L2L reduces degree by 1, so exit if degree of harmonics is 0
    ! or -1 (which means no FMM at all)
    if (ddx_data % pl .le. 0) return
    ! Set coordinate transformations
    zx_coord_transform = zero
    zx_coord_transform(3, 2) = one
    zx_coord_transform(2, 3) = one
    zx_coord_transform(1, 1) = one
    zy_coord_transform = zero
    zy_coord_transform(1, 2) = one
    zy_coord_transform(2, 1) = one
    zy_coord_transform(3, 3) = one
    ! At first reflect harmonics of a degree up to pl
    do isph = 1, ddx_data % nsph
        inode = ddx_data % snode(isph)
        sph_l_grad(:, 3, isph) = node_l(:, inode)
        call fmm_sph_transform(ddx_data % pl, zx_coord_transform, one, &
            & node_l(:, inode), zero, sph_l_grad(:, 1, isph))
        call fmm_sph_transform(ddx_data % pl, zy_coord_transform, one, &
            & node_l(:, inode), zero, sph_l_grad(:, 2, isph))
    end do
    ! Derivative of L2L translation over OZ axis at the origin consists of 2
    ! steps:
    !   1) decrease degree l and scale by sqrt((2*l-1)*(l*l-m*m)) / sqrt(2*l+1)
    !   2) scale by 1/rsph(isph)
    do l = 1, ddx_data % pl
        indi = l*l + l + 1
        indj = indi - 2*l
        tmp1 = -sqrt(dble(2*l-1)) / sqrt(dble(2*l+1))
        do m = 1-l, l-1
            tmp2 = sqrt(dble(l*l-m*m)) * tmp1
            sph_l_grad(indj+m, :, :) = tmp2 * sph_l_grad(indi+m, :, :)
        end do
    end do
    ! Scale by 1/rsph(isph) and rotate harmonics of degree up to pl-1 back to
    ! the initial axis. Coefficient of pl-th degree is zero so we ignore it.
    do isph = 1, ddx_data % nsph
        sph_l_grad(1:ddx_data % pl**2, 3, isph) = &
            & sph_l_grad(1:ddx_data % pl**2, 3, isph) / ddx_data % rsph(isph)
        work(1:ddx_data % pl**2, isph) = &
            & sph_l_grad(1:ddx_data % pl**2, 1, isph) / ddx_data % rsph(isph)
        call fmm_sph_transform(ddx_data % pl-1, zx_coord_transform, one, &
            & work(1:ddx_data % pl**2, isph), zero, &
            & sph_l_grad(1:ddx_data % pl**2, 1, isph))
        work(1:ddx_data % pl**2, isph) = &
            & sph_l_grad(1:ddx_data % pl**2, 2, isph) / ddx_data % rsph(isph)
        call fmm_sph_transform(ddx_data % pl-1, zy_coord_transform, one, &
            & work(1:ddx_data % pl**2, isph), zero, &
            & sph_l_grad(1:ddx_data % pl**2, 2, isph))
    end do
    ! Set degree pl to zero to avoid problems if user actually uses it
    l = ddx_data % pl
    indi = l*l + l + 1
    sph_l_grad(indi-l:indi+l, :, :) = zero
end subroutine tree_grad_l2l

end module ddx_core

