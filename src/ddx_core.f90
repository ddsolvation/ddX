!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_core.f90
!! Core routines and parameters of the entire ddX software
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

!> Core routines and parameters of the ddX software
module ddx_core
! Get ddx_params_type and all compile-time definitions
use ddx_parameters
! Get ddx_constants_type and all run-time constants
use ddx_constants
! Get ddx_workspace_type for temporary buffers
use ddx_workspace
! Get harmonics-related functions
use ddx_harmonics
! Enable OpenMP
use omp_lib
implicit none

!> @defgroup Fortran_interface_core Fortran interface: core routines

!> This defined type contains the primal and adjoint RHSs, the solution of
!! the primal and adjoint linear systems, useful intermediates for the
!! computation of the forces, and the information about the convergence of
!! the linear system solver (time, number of iterations, relative difference
!! at each iteration).
type ddx_state_type
    !! High-level entities to be stored and accessed by users

    !!
    !! Quantities common to all models.
    !!
    !> Representation of the solute density in spherical harmonics
    !! (\f$ \Psi \f$). It is used as RHS for the adjoint linear system.
    !! Dimension (nbasis, nsph).
    real(dp), allocatable :: psi(:, :)
    !> Electric potential at the cavity points. It is used to construct
    !! the RHS for the primal linear system. Dimension (ncav).
    real(dp), allocatable :: phi_cav(:)
    !> Potential at all the grid points. Dimension (ngrid, nsph).
    real(dp), allocatable :: phi_grid(:, :)
    !> Zeta intermediate for the forces. Dimension (ncav).
    real(dp), allocatable :: zeta(:)
    !> Error flag, nonzero in case of an error.
    integer :: error_flag
    !> Last error message.
    character(len=255) :: error_message

    !!
    !! ddCOSMO quantities (used also by ddPCM).
    !!
    !> Variable \f$ \Phi \f$ of a dimension (nbasis, nsph).
    !! Allocated and used by the COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: phi(:, :)
    !> Solution of the ddCOSMO system of a dimension (nbasis, nsph).
    !! Allocated and used by the COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: xs(:, :)
    !> Number of iterations to solve the primal ddCOSMO linear system.
    !! Used by the COSMO (model=1) and PCM (model=2) models.
    integer :: xs_niter
    !> Relative error of the iterative solver at each iteration of the primal
    !! ddCOSMO linear system, dimension (maxiter).
    !! Allocated and used by the COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: xs_rel_diff(:)
    !> Time to solve the primal ddCOSMO linear system.
    !! Used by the COSMO (model=1) and PCM (model=2) models.
    real(dp) :: xs_time
    !> Values of s at grid points. Dimension is (ngrid, nsph).
    !! Allocated and used by the COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: sgrid(:, :)
    !> Solution of the adjoint ddCOSMO linear system, dimension (nbasis, nsph).
    !! Allocated and used by the COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: s(:, :)
    !> Number of iterations required to solve the adjoint ddCOSMO linear
    !! system. Used by the COSMO (model=1) and PCM (model=2) models.
    integer :: s_niter
    !> Relative error of the iterative solver at each iteration of the adjoint
    !! ddCOSMO linear system, dimension (maxiter).
    !! Allocated and used by the COSMO (model=1) and PCM (model=2) models.
    real(dp), allocatable :: s_rel_diff(:)
    !> Time to solve the adjoint ddCOSMO linear system.
    !! Used by the COSMO (model=1) and PCM (model=2) models.
    real(dp) :: s_time

    !!
    !! ddPCM specific quantities.
    !!
    !> Variable \f$ \Phi_\infty \f$ of a dimension (nbasis, nsph).
    !! Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: phiinf(:, :)
    !> Variable \f$ \Phi_\varepsilon \f$ of a dimension (nbasis, nsph).
    !! Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: phieps(:, :)
    !> Number of iterations to solve the primal ddPCM linear system system.
    !! Allocated and used only by PCM (model=2) model.
    integer :: phieps_niter
    !> Relative error of the iterative solver at each iteration of the primal
    !! ddPCM linear system, dimension (maxiter).
    !! Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: phieps_rel_diff(:)
    !> Time to solve the primal ddPCM linear system.
    !! Used only by PCM (model=2) model.
    real(dp) :: phieps_time
    !> Shortcut of \f$ \Phi_\varepsilon - \Phi \f$ for the computation of
    !! the ddPCM forces.
    !! Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: g(:, :)
    !> Solution of the adjoint ddPCM linear system, dimension (nbasis, nsph).
    !! Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: y(:, :)
    !> Number of iteration to solve the adjoint ddPCM linear system.
    !! Used only by PCM (model=2) model.
    integer :: y_niter
    !> Relative error of the iterative solver at each iteration of the adjoint
    !! ddPCM linear system. Allocated and used only by the PCM (model=2) model.
    real(dp), allocatable :: y_rel_diff(:)
    !> Time to solve adjoint ddPCM system
    real(dp) :: y_time
    !> Solution of the adjoint ddPCM linear system, dimension (nbasis, nsph).
    !! Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: ygrid(:, :)
    !> Effective total adjoint solution of the ddPCM model, defined as
    !! \f$ Q := S - \frac{4\pi}{\varepsilon-1}Y \f$. Dimension (nbasis, nsph).
    !! Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: q(:, :)
    !> Values of Q at grid points. Dimension (ngrid, nsph).
    !! Allocated and used only by PCM (model=2) model.
    real(dp), allocatable :: qgrid(:, :)

    !!
    !! ddLPB quantities
    !!
    !> Solution to the ddLPB linear system. Dimension (nbasis, nsph, 2).
    !! Allocated and used only by the LPB (model=3) model.
    real(dp), allocatable :: x_lpb(:,:,:)
    !> Solution to the ddLPB adjoint linear system.
    !! Dimension (nbasis, nsph, 2).
    !! Allocated and used only by the LPB (model=3) model.
    real(dp), allocatable :: x_adj_lpb(:,:,:)
    !> g RHS for ddLPB.
    !! Allocated and used only by the LPB (model=3) model.
    real(dp), allocatable :: g_lpb(:,:)
    !> f RHS for ddLPB.
    !! Allocated and used only by the LPB (model=3) model.
    real(dp), allocatable :: f_lpb(:,:)
    !> Number of iterations required to solve the primal ddLPB linear system.
    !! Allocated and used only by the LPB (model=3) model.
    integer :: x_lpb_niter
    !> Number of iterations required to solve the adjoint ddLPB linear ststem.
    !! Allocated and used only by the LPB (model=3) model.
    integer :: x_adj_lpb_niter
    !> Time required to solve the primal ddLPB linear system.
    !! Allocated and used only by the LPB (model=3) model.
    real(dp) :: x_lpb_time
    !> Time required to solve the adjoint ddLPB linear system.
    !! Allocated and used only by the LPB (model=3) model.
    real(dp) :: x_adj_lpb_time
    !> Relative error of the iterative solver at each iteration of the primal
    !! ddLPB linear system, dimension (maxiter).
    !! Allocated and used only by the LPB (model=3) model.
    real(dp), allocatable :: x_lpb_rel_diff(:)
    !> Relative error of the iterative solver at each iteration of the adjoint
    !! ddLPB linear system, dimension (maxiter).
    !! Allocated and used only by the LPB (model=3) model.
    real(dp), allocatable :: x_adj_lpb_rel_diff(:)

end type ddx_state_type

!> Main ddX type that stores all the required information.
!! Container for the params, contants and workspace derived types.
type ddx_type
    !! New types inside the old one for an easier shift to the new design
    type(ddx_params_type) :: params
    type(ddx_constants_type) :: constants
    type(ddx_workspace_type) :: workspace
    !> Flag if there were an error
    integer :: error_flag = 0
    !> Last error message
    character(len=255) :: error_message
end type ddx_type

contains

!------------------------------------------------------------------------------
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
!! @param[in] ngrid: Number of Lebedev grid points `ngrid` >= 0
!! @param[in] force: 1 if forces are required and 0 otherwise
!! @param[in] fmm: 1 to use FMM acceleration and 0 otherwise
!! @param[in] pm: Maximal degree of multipole spherical harmonics. Ignored in
!!      the case `fmm=0`. Value -1 means no far-field FMM interactions are
!!      computed. `pm` >= -1.
!! @param[in] pl: Maximal degree of local spherical harmonics. Ignored in
!!      the case `fmm=0`. Value -1 means no far-field FMM interactions are
!!      computed. `pl` >= -1.
!! @param[in] se: Shift of characteristic function. -1 for interior, 0 for
!!      centered and 1 for outer regularization
!! @param[in] eta: Regularization parameter. 0 < eta <= 1.
!! @param[in] eps: Relative dielectric permittivity. eps > 1.
!! @param[in] kappa: Debye-H\"{u}ckel parameter
!! @param[in] matvecmem: handling of the sparse matrices. 1 for precomputin
!!      them and keeping them in memory, 0 for assembling the matrix-vector
!!      product on-the-fly. 
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
!------------------------------------------------------------------------------
subroutine ddinit(nsph, charge, x, y, z, rvdw, model, lmax, ngrid, force, &
        & fmm, pm, pl, se, eta, eps, kappa, &
        & matvecmem, maxiter, jacobi_ndiis, &
        & nproc, output_filename, ddx_data)
    ! Inputs
    integer, intent(in) :: nsph, model, lmax, force, fmm, pm, pl, &
        & matvecmem, maxiter, jacobi_ndiis, &
        & ngrid
    real(dp), intent(in):: charge(nsph), x(nsph), y(nsph), z(nsph), &
        & rvdw(nsph), se, eta, eps, kappa
    character(len=255), intent(in) :: output_filename
    ! Output
    type(ddx_type), target, intent(out) :: ddx_data
    ! Inouts
    integer, intent(inout) :: nproc
    ! Local variables
    real(dp), allocatable :: csph(:, :)
    integer :: istatus
    ! set the object in non error state
    ddx_data % error_flag = 0
    ddx_data % error_message = ''
    ! Init underlying objects
    allocate(csph(3, nsph), stat=istatus)
    if (istatus.ne.0) then
        write(ddx_data % error_message, *) 'Allocation failed in ddinit'
        ddx_data % error_flag = 1
        return
    end if
    csph(1, :) = x
    csph(2, :) = y
    csph(3, :) = z
    call params_init(model, force, eps, kappa, eta, se, lmax, ngrid, &
        & matvecmem, maxiter, jacobi_ndiis, &
        & fmm, pm, pl, nproc, nsph, charge, &
        & csph, rvdw, &
        & output_filename, ddx_data % params)
    if (ddx_data % params % error_flag .ne. 0) then
        ddx_data % error_flag = ddx_data % params % error_flag
        ddx_data % error_message = ddx_data % params % error_message
        return
    end if
    call constants_init(ddx_data % params, ddx_data % constants)
    if (ddx_data % constants % error_flag .ne. 0) then
        ddx_data % error_flag = ddx_data % constants % error_flag
        ddx_data % error_message = ddx_data % constants % error_message
        return
    end if
    call workspace_init(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace)
    if (ddx_data % workspace % error_flag .ne. 0) then
        ddx_data % error_flag = ddx_data % workspace % error_flag
        ddx_data % error_message = ddx_data % workspace % error_message
        return
    end if
    deallocate(csph, stat=istatus)
    if (istatus.ne.0) then
        write(ddx_data % error_message, *) 'Deallocation failed in ddinit'
        ddx_data % error_flag = 1
        return
    end if
end subroutine ddinit

!> Initialize the ddx_state object
!> @ingroup Fortran_interface_core
!!
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] state: ddx state (contains solutions and RHSs)
!!
subroutine ddx_init_state(params, constants, state)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(out) :: state
    integer :: istatus

    state % error_flag = 0
    state % error_message = ''

    allocate(state % psi(constants % nbasis, params % nsph), stat=istatus)
    if (istatus .ne. 0) then
        state % error_flag = 1
        state % error_message = "ddinit: `psi` allocation failed"
        return
    end if
    allocate(state % phi_cav(constants % ncav), stat=istatus)
    if (istatus .ne. 0) then
        state % error_flag = 1
        state % error_message = "ddinit: `phi_cav` allocation failed"
        return
    end if

    ! COSMO model
    if (params % model .eq. 1) then
        allocate(state % phi_grid(params % ngrid, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `phi_grid` " // &
                & "allocation failed"
            return
        end if
        allocate(state % phi(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `phi` " // &
                & "allocation failed"
            return
        end if
        allocate(state % xs(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `xs` " // &
                & "allocation failed"
            return
        end if
        allocate(state % xs_rel_diff(params % maxiter), &
            & stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `xs_rel_diff` " // &
                & "allocation failed"
            return
        end if
        allocate(state % s(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `s` " // &
            & "allocation failed"
            return
        end if
        allocate(state % s_rel_diff(params % maxiter), &
            & stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `s_rel_diff` " // &
                & "allocation failed"
            return
        end if
        allocate(state % sgrid(params % ngrid, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `sgrid` " // &
                & "allocation failed"
            return
        end if
        allocate(state % zeta(constants % ncav), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `zeta` " // &
                & "allocation failed"
            return
        end if
    ! PCM model
    else if (params % model .eq. 2) then
        allocate(state % phi_grid(params % ngrid, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `phi_grid` " // &
                & "allocation failed"
            return
        end if
        allocate(state % phi(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `phi` " // &
                & "allocation failed"
            return
        end if
        allocate(state % phiinf(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `phiinf` " // &
                & "allocation failed"
            return
        end if
        allocate(state % phieps(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `phieps` " // &
                & "allocation failed"
            return
        end if
        allocate(state % phieps_rel_diff(params % maxiter), &
            & stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `xs_rel_diff` " // &
                & "allocation failed"
            return
        end if
        allocate(state % xs(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `xs` " // &
                & "allocation failed"
            return
        end if
        allocate(state % xs_rel_diff(params % maxiter), &
            & stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `xs_rel_diff` " // &
                & "allocation failed"
            return
        end if
        allocate(state % s(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `s` " // &
                & "allocation failed"
            return
        end if
        allocate(state % s_rel_diff(params % maxiter), &
            & stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `xs_rel_diff` " // &
                & "allocation failed"
            return
        end if
        allocate(state % sgrid(params % ngrid, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `sgrid` " // &
                & "allocation failed"
            return
        end if
        allocate(state % y(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `y` " // &
                & "allocation failed"
            return
        end if
        allocate(state % y_rel_diff(params % maxiter), &
            & stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `y_rel_diff` " // &
                & "allocation failed"
            return
        end if
        allocate(state % ygrid(params % ngrid, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `ygrid` " // &
                & "allocation failed"
            return
        end if
        allocate(state % g(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `g` " // &
                & "allocation failed"
            return
        end if
        allocate(state % q(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `q` " // &
                & "allocation failed"
            return
        end if
        allocate(state % qgrid(params % ngrid, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `qgrid` " // &
                & "allocation failed"
            return
        end if
        allocate(state % zeta(constants % ncav), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `zeta` " // &
                & "allocation failed"
            return
        end if
    ! LPB model
    else if (params % model .eq. 3) then
        allocate(state % phi_grid(params % ngrid, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `phi_grid` " // &
                & "allocation failed"
            return
        end if
        allocate(state % phi(constants % nbasis, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `phi` " // &
                & "allocation failed"
            return
        end if
        allocate(state % zeta(constants % ncav), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `zeta` " // &
                & "allocation failed"
            !write(*, *) "Error in allocation of M2P matrices"
            return
        end if
        allocate(state % x_lpb_rel_diff(params % maxiter), &
            & stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `x_lpb_rel_diff` " // &
                & "allocation failed"
            return
        end if
        allocate(state % x_lpb(constants % nbasis, &
            & params % nsph, 2), &
            & stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `x_lpb` " // &
                & "allocation failed"
            return
        end if
        allocate(state % x_adj_lpb(constants % nbasis, &
            & params % nsph, 2), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `x_adj_lpb` " // &
                & "allocation failed"
            return
        end if
        allocate(state % x_adj_lpb_rel_diff(params % maxiter), &
            & stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `x_adj_lpb_rel_diff` " // &
                & "allocation failed"
            return
        end if
        allocate(state % g_lpb(params % ngrid, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `g_lpb` " // &
                & "allocation failed"
            return
        end if
        allocate(state % f_lpb(params % ngrid, &
            & params % nsph), stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "ddinit: `f_lpb` " // &
                & "allocation failed"
            return
        end if
    end if
end subroutine ddx_init_state

!------------------------------------------------------------------------------
!> Read configuration from a file
!!
!! @param[in] fname: Filename containing all the required info
!! @param[out] ddx_data: Object containing all inputs
!! @param[out] info: flag of succesfull exit
!!      = 0: Succesfull exit
!!      < 0: If info=-i then i-th argument had an illegal value
!!      > 0: Allocation of a buffer for the output ddx_data failed
!------------------------------------------------------------------------------
subroutine ddfromfile(fname, ddx_data, tol)
    ! Input
    character(len=*), intent(in) :: fname
    ! Outputs
    type(ddx_type), intent(out) :: ddx_data
    real(dp), intent(out) :: tol
    ! Local variables
    integer :: nproc, model, lmax, ngrid, force, fmm, pm, pl, &
        & nsph, i, matvecmem, maxiter, jacobi_ndiis, &
        & istatus
    real(dp) :: eps, se, eta, kappa
    real(dp), allocatable :: charge(:), x(:), y(:), z(:), rvdw(:)
    character(len=255) :: output_filename
    !! Read all the parameters from the file
    ! Open a configuration file
    open(unit=100, file=fname, form='formatted', access='sequential')
    ! Printing flag
    read(100, *) output_filename
    ! Number of OpenMP threads to be used
    read(100, *) nproc
    if(nproc .lt. 0) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 2nd line of a config file ", fname, &
            & ": `nproc` must be a positive integer value."
        ddx_data % error_flag = 1
        return
    end if
    ! Model to be used: 1 for COSMO, 2 for PCM and 3 for LPB
    read(100, *) model
    if((model .lt. 1) .or. (model .gt. 3)) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 3rd line of a config file ", fname, &
            & ": `model` must be an integer of a value 1, 2 or 3."
        ddx_data % error_flag = 1
        return
    end if
    ! Max degree of modeling spherical harmonics
    read(100, *) lmax
    if(lmax .lt. 0) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 4th line of a config file ", fname, &
            & ": `lmax` must be a non-negative integer value."
        ddx_data % error_flag = 1
        return
    end if
    ! Approximate number of Lebedev points
    read(100, *) ngrid
    if(ngrid .lt. 0) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 5th line of a config file ", fname, &
            & ": `ngrid` must be a non-negative integer value."
        ddx_data % error_flag = 1
        return
    end if
    ! Dielectric permittivity constant of the solvent
    read(100, *) eps
    if(eps .lt. zero) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 6th line of a config file ", fname, &
            & ": `eps` must be a non-negative floating point value."
        ddx_data % error_flag = 1
        return
    end if
    ! Shift of the regularized characteristic function
    read(100, *) se
    if((se .lt. -one) .or. (se .gt. one)) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 7th line of a config file ", fname, &
            & ": `se` must be a floating point value in a range [-1, 1]."
        ddx_data % error_flag = 1
        return
    end if
    ! Regularization parameter
    read(100, *) eta
    if((eta .lt. zero) .or. (eta .gt. one)) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 8th line of a config file ", fname, &
            & ": `eta` must be a floating point value in a range [0, 1]."
        ddx_data % error_flag = 1
        return
    end if
    ! Debye H\"{u}ckel parameter
    read(100, *) kappa
    if(kappa .lt. zero) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 9th line of a config file ", fname, &
            & ": `kappa` must be a non-negative floating point value."
        ddx_data % error_flag = 1
        return
    end if
    ! whether the (sparse) matrices are precomputed and kept in memory (1) or not (0).
    read(100, *) matvecmem 
    if((matvecmem.lt. 0) .or. (matvecmem .gt. 1)) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 10th line of a config file ", fname, &
            & ": `matvecmem` must be an integer value of a value 0 or 1."
        ddx_data % error_flag = 1
        return
    end if
    ! Relative convergence threshold for the iterative solver
    read(100, *) tol
    if((tol .lt. 1d-14) .or. (tol .gt. one)) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 12th line of a config file ", fname, &
            & ": `tol` must be a floating point value in a range [1d-14, 1]."
        ddx_data % error_flag = 1
        return
    end if
    ! Maximum number of iterations for the iterative solver
    read(100, *) maxiter
    if((maxiter .le. 0)) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 13th line of a config file ", fname, &
            & ": `maxiter` must be a positive integer value."
        ddx_data % error_flag = 1
        return
    end if
    ! Number of extrapolation points for Jacobi/DIIS solver
    read(100, *) jacobi_ndiis
    if((jacobi_ndiis .lt. 0)) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 14th line of a config file ", fname, &
            & ": `jacobi_ndiis` must be a non-negative integer value."
        ddx_data % error_flag = 1
        return
    end if
    ! Whether to compute (1) or not (0) forces as analytical gradients
    read(100, *) force
    if((force .lt. 0) .or. (force .gt. 1)) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 17th line of a config file ", fname, &
            & ": `force` must be an integer value of a value 0 or 1."
        ddx_data % error_flag = 1
        return
    end if
    ! Whether to use (1) or not (0) the FMM to accelerate computations
    read(100, *) fmm
    if((fmm .lt. 0) .or. (fmm .gt. 1)) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 18th line of a config file ", fname, &
            & ": `fmm` must be an integer value of a value 0 or 1."
        ddx_data % error_flag = 1
        return
    end if
    ! Max degree of multipole spherical harmonics for the FMM
    read(100, *) pm
    if(pm .lt. 0) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 19th line of a config file ", fname, &
            & ": `pm` must be a non-negative integer value."
        ddx_data % error_flag = 1
        return
    end if
    ! Max degree of local spherical harmonics for the FMM
    read(100, *) pl
    if(pl .lt. 0) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 20th line of a config file ", fname, &
            & ": `pl` must be a non-negative integer value."
        ddx_data % error_flag = 1
        return
    end if
    ! Number of input spheres
    read(100, *) nsph
    if(nsph .le. 0) then
        write(ddx_data % error_message, "(3A)") &
            & "Error on the 21th line of a config file ", fname, &
            & ": `nsph` must be a positive integer value."
        ddx_data % error_flag = 1
        return
    end if
    ! Coordinates, radii and charges
    allocate(charge(nsph), x(nsph), y(nsph), z(nsph), rvdw(nsph), stat=istatus)
    if(istatus .ne. 0) then
        write(ddx_data % error_message, "(2A)") &
            & "Could not allocate space for coordinates, radii ", &
            & "and charges of atoms"
        ddx_data % error_flag = 1
        return
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

    ! adjust ngrid
    call closest_supported_lebedev_grid(ngrid)

    !! Initialize ddx_data object
    call ddinit(nsph, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, &
        & pm, pl, se, eta, eps, kappa, matvecmem, &
        & maxiter, jacobi_ndiis, nproc, output_filename, ddx_data)
    !! Clean local temporary data
    deallocate(charge, x, y, z, rvdw, stat=istatus)
    if(istatus .ne. 0) then
        write(ddx_data % error_message, "(2A)") &
            & "Could not deallocate space for coordinates, ", &
            & "radii and charges of atoms"
        ddx_data % error_flag = 1
        return
    end if
end subroutine ddfromfile

!------------------------------------------------------------------------------
!> Deallocate object with corresponding data
!!
!! @param[inout] ddx_data: object to deallocate
!------------------------------------------------------------------------------
subroutine ddfree(ddx_data)
    ! Input/output
    type(ddx_type), intent(inout) :: ddx_data
    ! Local variables
    call workspace_free(ddx_data % workspace)
    if (ddx_data % workspace % error_flag .ne. 0) then
        write(ddx_data % error_message, *) "workspace_free failed!"
        ddx_data % error_flag = 1
        return
    end if
    call constants_free(ddx_data % constants)
    if (ddx_data % workspace % error_flag .ne. 0) then
        write(ddx_data % error_message, *) "constants_free failed!"
        ddx_data % error_flag = 1
        return
    end if
    call params_free(ddx_data % params)
    if (ddx_data % workspace % error_flag .ne. 0) then
        write(ddx_data % error_message, *) "params_free failed!"
        ddx_data % error_flag = 1
        return
    end if
end subroutine ddfree

!> Deallocate the ddx_state object
!> @ingroup Fortran_interface_core
!!
!! @param[inout] state: ddx state (contains solutions and RHSs)
!!
subroutine ddx_free_state(state)
    implicit none
    type(ddx_state_type), intent(inout) :: state
    integer :: istatus

    if (allocated(state % phi_cav)) then
        deallocate(state % phi_cav, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`phi_cav` deallocation failed!"
            return
        endif
    end if
    if (allocated(state % psi)) then
        deallocate(state % psi, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`psi` deallocation failed!"
            return
        endif
    end if
    if (allocated(state % phi_grid)) then
        deallocate(state % phi_grid, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`phi_grid` deallocation failed!"
            return
        endif
    end if
    if (allocated(state % phi)) then
        deallocate(state % phi, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`phi` deallocation failed!"
        endif
    end if
    if (allocated(state % phiinf)) then
        deallocate(state % phiinf, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`phiinf` deallocation failed!"
        endif
    end if
    if (allocated(state % phieps)) then
        deallocate(state % phieps, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`phieps` deallocation failed!"
        endif
    end if
    if (allocated(state % phieps_rel_diff)) then
        deallocate(state % phieps_rel_diff, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`phieps_rel_diff` deallocation failed!"
        endif
    end if
    if (allocated(state % xs)) then
        deallocate(state % xs, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`xs` deallocation failed!"
        endif
    end if
    if (allocated(state % xs_rel_diff)) then
        deallocate(state % xs_rel_diff, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`xs_rel_diff` deallocation failed!"
        endif
    end if
    if (allocated(state % s)) then
        deallocate(state % s, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`s` deallocation failed!"
        endif
    end if
    if (allocated(state % s_rel_diff)) then
        deallocate(state % s_rel_diff, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`s_rel_diff` deallocation failed!"
        endif
    end if
    if (allocated(state % sgrid)) then
        deallocate(state % sgrid, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`sgrid` deallocation failed!"
        endif
    end if
    if (allocated(state % y)) then
        deallocate(state % y, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`y` deallocation failed!"
        endif
    end if
    if (allocated(state % y_rel_diff)) then
        deallocate(state % y_rel_diff, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`y_rel_diff` deallocation failed!"
        endif
    end if
    if (allocated(state % ygrid)) then
        deallocate(state % ygrid, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`ygrid` deallocation failed!"
        endif
    end if
    if (allocated(state % g)) then
        deallocate(state % g, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`g` deallocation failed!"
        endif
    end if
    if (allocated(state % q)) then
        deallocate(state % q, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`q` deallocation failed!"
        endif
    end if
    if (allocated(state % qgrid)) then
        deallocate(state % qgrid, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`qgrid` deallocation failed!"
        endif
    end if
    if (allocated(state % zeta)) then
        deallocate(state % zeta, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`zeta` deallocation failed!"
        endif
    end if
    if (allocated(state % x_lpb)) then
        deallocate(state % x_lpb, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`x_lpb` deallocation failed!"
        endif
    end if
    if (allocated(state % x_lpb_rel_diff)) then
        deallocate(state % x_lpb_rel_diff, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`x_lpb_rel_diff` deallocation failed!"
        end if
    end if
    if (allocated(state % x_adj_lpb)) then
        deallocate(state % x_adj_lpb, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "x_adj_lpb deallocation failed!"
        endif
    end if
    if (allocated(state % x_adj_lpb_rel_diff)) then
        deallocate(state % x_adj_lpb_rel_diff, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`x_adj_lpb_rel_diff` deallocation failed!"
        end if
    end if
    if (allocated(state % g_lpb)) then
        deallocate(state % g_lpb, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`g_lpb` deallocation failed!"
        endif
    end if
    if (allocated(state % f_lpb)) then
        deallocate(state % f_lpb, stat=istatus)
        if (istatus .ne. 0) then
            state % error_flag = 1
            state % error_message = "`f_lpb` deallocation failed!"
        endif
    end if
end subroutine ddx_free_state

!------------------------------------------------------------------------------
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
subroutine print_spherical(iunit, label, nbasis, lmax, ncol, icol, x)
    ! Inputs
    character (len=*), intent(in) :: label
    integer, intent(in) :: nbasis, lmax, ncol, icol, iunit
    real(dp), intent(in) :: x(nbasis, ncol)
    ! Local variables
    integer :: l, m, ind, noff, nprt, ic, j
    ! Print header:
    if (ncol .eq. 1) then
        write (iunit,'(3x,a,1x,"(column ",i4")")') label, icol
    else
        write (iunit,'(3x,a)') label
    endif
    ! Print entries:
    if (ncol .eq. 1) then
        do l = 0, lmax
            ind = l*l + l + 1
            do m = -l, l
                write(iunit,1000) l, m, x(ind+m, 1)
            end do
        end do
    else
        noff = mod(ncol, 5)
        nprt = max(ncol-noff, 0)
        do ic = 1, nprt, 5
            write(iunit,1010) (j, j = ic, ic+4)
            do l = 0, lmax
                ind = l*l + l + 1
                do m = -l, l
                    write(iunit,1020) l, m, x(ind+m, ic:ic+4)
                end do
            end do
        end do
        write (iunit,1010) (j, j = nprt+1, nprt+noff)
        do l = 0, lmax
            ind = l*l + l + 1
            do m = -l, l
                write(iunit,1020) l, m, x(ind+m, nprt+1:nprt+noff)
            end do
        end do
    end if
    1000 format(1x,i3,i4,f14.8)
    1010 format(8x,5i14)
    1020 format(1x,i3,i4,5f14.8)
end subroutine print_spherical

!------------------------------------------------------------------------------
!> Print array of quadrature points
!!
!! Prints (ngrid, ncol) array
!!
!! @param[in] label: Label to print
!! @param[in] ngrid: Number of rows of input x. ngrid >= 0
!! @param[in] ncol: Number of columns to print
!! @param[in] icol: This number is only for printing purposes
!! @param[in] x: Actual data to print
subroutine print_nodes(iunit, label, ngrid, ncol, icol, x)
    ! Inputs
    character (len=*), intent(in) :: label
    integer, intent(in) :: ngrid, ncol, icol, iunit
    real(dp), intent(in) :: x(ngrid, ncol)
    ! Local variables
    integer :: ig, noff, nprt, ic, j
    ! Print header :
    if (ncol .eq. 1) then
        write (iunit,'(3x,a,1x,"(column ",i4")")') label, icol
    else
        write (iunit,'(3x,a)') label
    endif
    ! Print entries :
    if (ncol .eq. 1) then
        do ig = 1, ngrid
            write(iunit,1000) ig, x(ig, 1)
        enddo
    else
        noff = mod(ncol, 5)
        nprt = max(ncol-noff, 0)
        do ic = 1, nprt, 5
            write(iunit,1010) (j, j = ic, ic+4)
            do ig = 1, ngrid
                write(iunit,1020) ig, x(ig, ic:ic+4)
            end do
        end do
        write (iunit,1010) (j, j = nprt+1, nprt+noff)
        do ig = 1, ngrid
            write(iunit,1020) ig, x(ig, nprt+1:nprt+noff)
        end do
    end if
    !
    1000 format(1x,i5,f14.8)
    1010 format(6x,5i14)
    1020 format(1x,i5,5f14.8)
    !
end subroutine print_nodes

!------------------------------------------------------------------------------
!> Print dd Solution vector
!!
!! @param[in] label : Label to print
!! @param[in] vector: Vector to print
!------------------------------------------------------------------------------
subroutine print_ddvector(ddx_data, label, vector)
    ! Inputs
    type(ddx_type), intent(in)  :: ddx_data
    character(len=*) :: label
    real(dp) :: vector(ddx_data % constants % nbasis, ddx_data % params % nsph)
    ! Local variables
    integer :: isph, lm
    write(6,*) label
    do isph = 1, ddx_data % params % nsph
      do lm = 1, ddx_data % constants % nbasis
        write(6,'(F15.8)') vector(lm,isph)
      end do
    end do
    return
end subroutine print_ddvector


!------------------------------------------------------------------------------
!> Integrate against spherical harmonics
!!
!! Integrate by Lebedev spherical quadrature. This function can be simply
!! substituted by a matrix-vector product.
!!
!! @param[in] nsph: Number of spheres. `nsph` >= 0.
!! @param[in] nbasis: Number of spherical harmonics. `nbasis` >= 0.
!! @param[in] ngrid: Number of Lebedev grid points. `ngrid` >= 0.
!! @param[in] vwgrid: Values of spherical harmonics at Lebedev grid points,
!!      multiplied by weights of grid points. Dimension is (ldvwgrid, ngrid).
!! @param[in] ldvwgrid: Leading dimension of `vwgrid`.
!! @param[in] x_grid: Input values at grid points of the sphere. Dimension is
!!      (ngrid, nsph).
!! @param[out] x_lm: Output spherical harmonics. Dimension is (nbasis, nsph).
subroutine ddintegrate(nsph, nbasis, ngrid, vwgrid, ldvwgrid, x_grid, x_lm)
    !! Inputs
    integer, intent(in) :: nsph, nbasis, ngrid, ldvwgrid
    real(dp), intent(in) :: vwgrid(ldvwgrid, ngrid)
    real(dp), intent(in) :: x_grid(ngrid, nsph)
    !! Output
    real(dp), intent(out) :: x_lm(nbasis, nsph)
    !! Just call a single dgemm to do the job
    call dgemm('N', 'N', nbasis, nsph, ngrid, one, vwgrid, ldvwgrid, x_grid, &
        & ngrid, zero, x_lm, nbasis)
end subroutine ddintegrate

!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
subroutine dbasis(params, constants, x, basloc, dbsloc, vplm, vcos, vsin)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), dimension(3),        intent(in)    :: x
    real(dp), dimension(constants % nbasis),   intent(inout) :: basloc, vplm
    real(dp), dimension(3, constants % nbasis), intent(inout) :: dbsloc
    real(dp), dimension(params % lmax+1),   intent(inout) :: vcos, vsin
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
    call polleg_work( cthe, sthe, params % lmax, vplm, vcos )
    !
    !     evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is 
    !     pointless if z = 1, as the only non vanishing terms will be the 
    !     ones with m=0.
    !
    !     not ( NORTH or SOUTH pole )
    if ( sthe.ne.zero ) then
        call trgev( cphi, sphi, params % lmax, vcos, vsin )
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
    do l = 0, params % lmax
        ind = l*l + l + 1
        ! m = 0
        fln = constants % vscales(ind)   
        basloc(ind) = fln*vplm(ind)
        if (l.gt.0) then
            dbsloc(:,ind) = fln*vplm(ind+1)*et(:)
        else
            dbsloc(:,ind) = zero
        end if
        !dir$ simd
        do m = 1, l
            fln = constants % vscales(ind+m)
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
!> TODO
!------------------------------------------------------------------------------
real(dp) function intmlp(params, constants, t, sigma, basloc )
!  
      type(ddx_params_type), intent(in) :: params
      type(ddx_constants_type), intent(in) :: constants
      real(dp), intent(in) :: t
      real(dp), dimension(constants % nbasis), intent(in) :: sigma, basloc
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
      do l = 0, params % lmax
!      
        ind = l*l + l + 1
!
!       update factor 4pi / (2l+1) * t^l
        fac = tt / constants % vscales(ind)**2
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

!------------------------------------------------------------------------------
!> Weigh potential at cavity points by characteristic function
!> TODO use cavity points in CSR format
!------------------------------------------------------------------------------
subroutine wghpot(ncav, phi_cav, nsph, ngrid, ui, phi_grid, g)
    !! Inputs
    integer, intent(in) :: ncav, nsph, ngrid
    real(dp), intent(in) :: phi_cav(ncav), ui(ngrid, nsph)
    !! Outputs
    real(dp), intent(out) :: g(ngrid, nsph), phi_grid(ngrid, nsph)
    !! Local variables
    integer isph, igrid, icav
    !! Code
    ! Initialize
    icav = 0 
    g = zero
    phi_grid = zero
    ! Loop over spheres
    do isph = 1, nsph
        ! Loop over points
        do igrid = 1, ngrid
            ! Non-zero contribution from point
            if (ui(igrid, isph) .ne. zero) then
                ! Advance cavity point counter
                icav = icav + 1
                phi_grid(igrid, isph) = phi_cav(icav)
                ! Weigh by (negative) characteristic function
                g(igrid, isph) = -ui(igrid, isph) * phi_cav(icav)
            endif
        enddo
    enddo
end subroutine wghpot

!------------------------------------------------------------------------------------------------
!> Compute the local Sobolev H^(-1/2)-norm on one sphere of u
!------------------------------------------------------------------------------
subroutine hsnorm(lmax, nbasis, u, unorm)
    integer, intent(in) :: lmax, nbasis
    real(dp), dimension(nbasis), intent(in) :: u
    real(dp), intent(inout) :: unorm
    integer :: l, m, ind
    real(dp)  :: fac
    ! initialize
    unorm = zero
    do l = 0, lmax
        ind = l*l + l + 1
        fac = one/(one + dble(l))
        do m = -l, l
            unorm = unorm + fac*u(ind+m)*u(ind+m)
        end do
    end do
!   the much neglected square root
    unorm = sqrt(unorm)
end subroutine hsnorm

!------------------------------------------------------------------------------
!> Compute the global  Sobolev H^(-1/2)-norm of x
!-------------------------------------------------------------------------------
real(dp) function hnorm(lmax, nbasis, nsph, x)
    integer, intent(in) :: lmax, nbasis, nsph
    real(dp),  dimension(nbasis, nsph), intent(in) :: x
    integer                                        :: isph, istatus
    real(dp)                                       :: vrms, fac
!
    vrms = 0.0_dp
    !$omp parallel do default(none) shared(nsph,lmax,nbasis,x) &
    !$omp private(isph,fac) schedule(dynamic) reduction(+:vrms)
    do isph = 1, nsph
        call hsnorm(lmax, nbasis, x(:,isph), fac)
        vrms = vrms + fac*fac
    enddo
!   call rmsvec(nsph, u, vrms, vmax)
    hnorm = sqrt(vrms/dble(nsph))
end function hnorm

!------------------------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
real(dp) function rmsnorm(lmax, nbasis, nsph, x)
    integer,                            intent(in) :: lmax, nbasis, nsph
    real(dp),  dimension(nbasis, nsph), intent(in) :: x
!
    integer  :: n
    real(dp) :: vrms, vmax
!
    n = nbasis*nsph
    call rmsvec(n,x,vrms,vmax)
!
    rmsnorm = vrms
  
end function rmsnorm

!------------------------------------------------------------------------------
!> compute root-mean-square and max norm
!------------------------------------------------------------------------------
subroutine rmsvec( n, v, vrms, vmax )
    implicit none
    integer,               intent(in)    :: n
    real(dp),  dimension(n), intent(in)    :: v
    real(dp),                intent(inout) :: vrms, vmax
    integer :: i
    real(dp), parameter :: zero=0.0d0

    vrms = zero
    vmax = zero
    do i = 1,n
      vmax = max(vmax,abs(v(i)))
      vrms = vrms + v(i)*v(i)
    end do
    ! the much neglected square root
    vrms = sqrt(vrms/dble(n))
endsubroutine rmsvec

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine adjrhs(params, constants, isph, xi, vlm, work)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: isph
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: xi
    real(dp), dimension(constants % nbasis), intent(inout) :: vlm
    real(dp), dimension(params % lmax+1), intent(inout) :: work

    integer :: ij, jsph, ig, l, ind, m
    real(dp)  :: vji(3), vvji, tji, sji(3), xji, oji, fac, ffac, t
    real(dp) :: rho, ctheta, stheta, cphi, sphi

    do ij = constants % inl(isph), constants % inl(isph+1)-1
      jsph = constants % nl(ij)
      do ig = 1, params % ngrid
        vji  = params % csph(:,jsph) + params % rsph(jsph)* &
            & constants % cgrid(:,ig) - params % csph(:,isph)
        vvji = sqrt(dot_product(vji,vji))
        tji  = vvji/params % rsph(isph)
        if ( tji.lt.( one + (params % se+one)/two*params % eta ) ) then
          xji = fsw( tji, params % se, params % eta )
          if ( constants % fi(ig,jsph).gt.one ) then
            oji = xji/ constants % fi(ig,jsph)
          else
            oji = xji
          endif
          fac = constants % wgrid(ig) * xi(ig,jsph) * oji
          call fmm_l2p_adj_work(vji, fac, params % rsph(isph), &
              & params % lmax, constants % vscales_rel, one, vlm, work)
        endif
      enddo
    enddo
end subroutine adjrhs

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine calcv(params, constants, isph, pot, sigma, work)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: isph
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: sigma
    real(dp), dimension(params % ngrid), intent(inout) :: pot
    real(dp), dimension(params % lmax+1), intent(inout) :: work

    integer :: its, ij, jsph
    real(dp) :: vij(3), sij(3)
    real(dp) :: vvij, tij, xij, oij, stslm, stslm2, stslm3, &
        & thigh, rho, ctheta, stheta, cphi, sphi

    thigh = one + pt5*(params % se + one)*params % eta
    pot(:) = zero
    ! loop over grid points
    do its = 1, params % ngrid
        ! contribution from integration point present
        if (constants % ui(its,isph).lt.one) then
            ! loop over neighbors of i-sphere
            do ij = constants % inl(isph), constants % inl(isph+1)-1
                jsph = constants % nl(ij)
                ! compute t_n^ij = | r_i + \rho_i s_n - r_j | / \rho_j
                vij  = params % csph(:,isph) + params % rsph(isph)* &
                    & constants % cgrid(:,its) - params % csph(:,jsph)
                vvij = sqrt( dot_product( vij, vij ) )
                tij  = vvij / params % rsph(jsph)
                ! point is INSIDE j-sphere
                if (tij.lt.thigh .and. tij.gt.zero) then
                    xij = fsw(tij, params % se, params % eta)
                    if (constants % fi(its,isph).gt.one) then
                        oij = xij / constants % fi(its,isph)
                    else
                        oij = xij
                    end if
                    call fmm_l2p_work(vij, params % rsph(jsph), params % lmax, &
                        & constants % vscales_rel, oij, sigma(:, jsph), one, &
                        & pot(its), work)
                end if
            end do
        end if
    end do
end subroutine calcv

!------------------------------------------------------------------------------
!> Evaluate values of spherical harmonics at Lebedev grid points
!------------------------------------------------------------------------------
subroutine ddeval_grid(params, constants, alpha, x_sph, beta, x_grid)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: alpha, x_sph(constants % nbasis, params % nsph), &
        & beta
    !! Output
    real(dp), intent(inout) :: x_grid(params % ngrid, params % nsph)
    !! The code
    call ddeval_grid_work(constants % nbasis, params % ngrid, params % nsph, &
        & constants % vgrid, constants % vgrid_nbasis, alpha, x_sph, beta, &
        & x_grid)
end subroutine ddeval_grid

!------------------------------------------------------------------------------
!> Evaluate values of spherical harmonics at Lebedev grid points
!------------------------------------------------------------------------------
subroutine ddeval_grid_work(nbasis, ngrid, nsph, vgrid, ldvgrid, alpha, &
        & x_sph, beta, x_grid)
    !! Inputs
    integer, intent(in) :: nbasis, ngrid, nsph, ldvgrid
    real(dp), intent(in) :: vgrid(ldvgrid, ngrid), alpha, &
        & x_sph(nbasis, nsph), beta
    !! Output
    real(dp), intent(inout) :: x_grid(ngrid, nsph)
    !! Local variables
    external :: dgemm
    !! The code
    call dgemm('T', 'N', ngrid, nsph, nbasis, alpha, vgrid, ldvgrid, x_sph, &
        & nbasis, beta, x_grid, ngrid)
end subroutine ddeval_grid_work

!> Integrate values at grid points into spherical harmonics
subroutine ddintegrate_sph(params, constants, alpha, x_grid, beta, x_sph)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: alpha, x_grid(params % ngrid, params % nsph), &
        & beta
    !! Output
    real(dp), intent(inout) :: x_sph(constants % nbasis, params % nsph)
    !! The code
    call ddintegrate_sph_work(constants % nbasis, params % ngrid, &
        & params % nsph, constants % vwgrid, constants % vgrid_nbasis, alpha, &
        & x_grid, beta, x_sph)
end subroutine ddintegrate_sph

!> Integrate values at grid points into spherical harmonics
subroutine ddintegrate_sph_work(nbasis, ngrid, nsph, vwgrid, ldvwgrid, alpha, &
        & x_grid, beta, x_sph)
    !! Inputs
    integer, intent(in) :: nbasis, ngrid, nsph, ldvwgrid
    real(dp), intent(in) :: vwgrid(ldvwgrid, ngrid), alpha, &
        & x_grid(ngrid, nsph), beta
    !! Outputs
    real(dp), intent(inout) :: x_sph(nbasis, nsph)
    !! Local variables
    !! The code
    call dgemm('N', 'N', nbasis, nsph, ngrid, alpha, vwgrid, ldvwgrid, &
        & x_grid, ngrid, beta, x_sph, nbasis)
end subroutine ddintegrate_sph_work

!------------------------------------------------------------------------------
!> Unwrap values at cavity points into values at all grid points
!------------------------------------------------------------------------------
subroutine ddcav_to_grid(params, constants, x_cav, x_grid)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x_cav(constants % ncav)
    !! Output
    real(dp), intent(inout) :: x_grid(params % ngrid, params % nsph)
    !! Local variables
    character(len=255) :: string
    !! The code
    call ddcav_to_grid_work(params % ngrid, params % nsph, constants % ncav, &
        & constants % icav_ia, constants % icav_ja, x_cav, x_grid)
end subroutine ddcav_to_grid

!------------------------------------------------------------------------------
!> Unwrap values at cavity points into values at all grid points
!------------------------------------------------------------------------------
subroutine ddcav_to_grid_work(ngrid, nsph, ncav, icav_ia, icav_ja, x_cav, &
        & x_grid)
    !! Inputs
    integer, intent(in) :: ngrid, nsph, ncav
    integer, intent(in) :: icav_ia(nsph+1), icav_ja(ncav)
    real(dp), intent(in) :: x_cav(ncav)
    !! Output
    real(dp), intent(out) :: x_grid(ngrid, nsph)
    !! Local variables
    integer :: isph, icav, igrid, igrid_old
    !! The code
    do isph = 1, nsph
        igrid_old = 0
        igrid = 0
        do icav = icav_ia(isph), icav_ia(isph+1)-1
            igrid = icav_ja(icav)
            x_grid(igrid_old+1:igrid-1, isph) = zero
            x_grid(igrid, isph) = x_cav(icav)
            igrid_old = igrid
        end do
        x_grid(igrid+1:ngrid, isph) = zero
    end do
end subroutine ddcav_to_grid_work

!> Integrate by a characteristic function at Lebedev grid points
!! \xi(n,i) = sum w_n U_n^i Y_l^m(s_n) [S_i]_l^m
!!            l,m
!!
subroutine ddproject_cav(params, constants, s, xi)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in)  :: s(constants%nbasis, params%nsph)
    real(dp), intent(out) :: xi(constants%ncav)
    integer :: its, isph, ii
    ii = 0
    do isph = 1, params%nsph
        do its = 1, params%ngrid
            if (constants%ui(its, isph) .gt. zero) then
                ii     = ii + 1
                xi(ii) = constants%ui(its, isph) &
                    &* dot_product(constants%vwgrid(:, its), s(:, isph))
            end if
        end do
    end do
end subroutine ddproject_cav

!------------------------------------------------------------------------------
!> Transfer multipole coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_m2m_rotation(params, constants, node_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pm**2 + 19*params % pm + 8)
    ! Call corresponding work routine
    call tree_m2m_rotation_work(params, constants, node_m, work)
end subroutine tree_m2m_rotation

!------------------------------------------------------------------------------
!> Transfer multipole coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_m2m_rotation_work(params, constants, node_m, work)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp), intent(out) :: work(6*params % pm**2 + 19*params % pm + 8)
    ! Local variables
    integer :: i, j
    real(dp) :: c1(3), c(3), r1, r
    ! Bottom-to-top pass
    !!$omp parallel do default(none) shared(constants,params,node_m) &
    !!$omp private(i,j,c1,c,r1,r,work)
    do i = constants % nclusters, 1, -1
        ! Leaf node does not need any update
        if (constants % children(1, i) == 0) cycle
        c = constants % cnode(:, i)
        r = constants % rnode(i)
        ! First child initializes output
        j = constants % children(1, i)
        c1 = constants % cnode(:, j)
        r1 = constants % rnode(j)
        c1 = c1 - c
        call fmm_m2m_rotation_work(c1, r1, r, &
            & params % pm, &
            & constants % vscales, &
            & constants % vcnk, one, &
            & node_m(:, j), zero, node_m(:, i), work)
        ! All other children update the same output
        do j = constants % children(1, i)+1, constants % children(2, i)
            c1 = constants % cnode(:, j)
            r1 = constants % rnode(j)
            c1 = c1 - c
            call fmm_m2m_rotation_work(c1, r1, r, params % pm, &
                & constants % vscales, constants % vcnk, one, &
                & node_m(:, j), one, node_m(:, i), work)
        end do
    end do
end subroutine tree_m2m_rotation_work

!------------------------------------------------------------------------------
!> Transfer multipole coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_m2m_bessel_rotation(params, constants, node_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Temporary workspace
    !real(dp) :: work(6*params % pm**2 + 19*params % pm + 8)
    ! Call corresponding work routine
    call tree_m2m_bessel_rotation_work(params, constants, node_m)
end subroutine tree_m2m_bessel_rotation

!------------------------------------------------------------------------------
!> Transfer multipole coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_m2m_bessel_rotation_work(params, constants, node_m)
    use complex_bessel
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pm**2 + 19*params % pm + 8)
    complex(dp) :: work_complex(2*params % pm+1)
    ! Local variables
    integer :: i, j
    real(dp) :: c1(3), c(3), r1, r
    ! Bottom-to-top pass
    do i = constants % nclusters, 1, -1
        ! Leaf node does not need any update
        if (constants % children(1, i) == 0) cycle
        c = constants % cnode(:, i)
        r = constants % rnode(i)
        ! First child initializes output
        j = constants % children(1, i)
        c1 = constants % cnode(:, j)
        r1 = constants % rnode(j)
!        call fmm_m2m_bessel_rotation(c1-c, r1, r, params % kappa, &
!            & params % pm, &
!            & constants % vscales, &
!            & constants % vcnk, one, &
!            & node_m(:, j), zero, node_m(:, i))
        c1 = params % kappa*(c1 - c)
        call fmm_m2m_bessel_rotation_work(c1, &
            & constants % SK_rnode(:, j), constants % SK_rnode(:, i), &
            & params % pm, &
            & constants % vscales, &
            & constants % vcnk, one, &
            & node_m(:, j), zero, node_m(:, i), work, work_complex)
        ! All other children update the same output
        do j = constants % children(1, i)+1, constants % children(2, i)
            c1 = constants % cnode(:, j)
            r1 = constants % rnode(j)
!            call fmm_m2m_bessel_rotation(c1-c, r1, r, params % kappa, &
!                & params % pm, &
!                & constants % vscales, constants % vcnk, one, &
!                & node_m(:, j), one, node_m(:, i))
            c1 = params % kappa*(c1 - c)
            call fmm_m2m_bessel_rotation_work(c1, &
                & constants % SK_rnode(:, j), constants % SK_rnode(:, i), &
                & params % pm, &
                & constants % vscales, &
                & constants % vcnk, one, &
                & node_m(:, j), one, node_m(:, i), work, work_complex)
        end do
    end do
end subroutine tree_m2m_bessel_rotation_work

!------------------------------------------------------------------------------
!> Adjoint transfer multipole coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_m2m_rotation_adj(params, constants, node_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pm**2 + 19*params % pm + 8)
    ! Call corresponding work routine
    call tree_m2m_rotation_adj_work(params, constants, node_m, work)
end subroutine tree_m2m_rotation_adj

!------------------------------------------------------------------------------
!> Adjoint transfer multipole coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_m2m_rotation_adj_work(params, constants, node_m, work)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp), intent(out) :: work(6*params % pm**2 + 19*params % pm + 8)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Top-to-bottom pass
    do i = 2, constants % nclusters
        j = constants % parent(i)
        c = constants % cnode(:, j)
        r = constants % rnode(j)
        c1 = constants % cnode(:, i)
        r1 = constants % rnode(i)
        c1 = c - c1
        call fmm_m2m_rotation_adj_work(c1, r, r1, params % pm, &
            & constants % vscales, constants % vcnk, one, node_m(:, j), one, &
            & node_m(:, i), work)
    end do
end subroutine tree_m2m_rotation_adj_work
!------------------------------------------------------------------------------
!> Adjoint transfer multipole coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_m2m_bessel_rotation_adj(params, constants, node_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Temporary workspace
    !real(dp), intent(out) :: work(6*params % pm**2 + 19*params % pm + 8)
    call tree_m2m_bessel_rotation_adj_work(params, constants, node_m)
end subroutine tree_m2m_bessel_rotation_adj

!------------------------------------------------------------------------------
!> Adjoint transfer multipole coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_m2m_bessel_rotation_adj_work(params, constants, node_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pm**2 + 19*params % pm + 8)
    complex(dp) :: work_complex(2*params % pm+1)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Top-to-bottom pass
    do i = 2, constants % nclusters
        j = constants % parent(i)
        c = constants % cnode(:, j)
        c1 = constants % cnode(:, i)
        c1 = params % kappa*(c1 - c)
        ! Little bit confusion about the i and j indices
        call fmm_m2m_bessel_rotation_adj_work(c1, constants % SK_rnode(:, i), &
            & constants % SK_rnode(:, j),&
            & params % pm, constants % vscales, constants % vcnk, one, &
            & node_m(:, j), one, node_m(:, i), work, work_complex)
    end do
end subroutine tree_m2m_bessel_rotation_adj_work

!------------------------------------------------------------------------------
!> Transfer local coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_l2l_rotation(params, constants, node_l)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pl**2 + 19*params % pl + 8)
    ! Call corresponding work routine
    call tree_l2l_rotation_work(params, constants, node_l, work)
end subroutine tree_l2l_rotation

!------------------------------------------------------------------------------
!> Transfer local coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_l2l_rotation_work(params, constants, node_l, work)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp), intent(out) :: work(6*params % pl**2 + 19*params % pl + 8)
    ! Local variables
    integer :: i, j
    real(dp) :: c1(3), c(3), r1, r
    ! Top-to-bottom pass
    !!$omp parallel do default(none) shared(constants,params,node_l) &
    !!$omp private(i,j,c1,c,r1,r,work)
    do i = 2, constants % nclusters
        j = constants % parent(i)
        c = constants % cnode(:, j)
        r = constants % rnode(j)
        c1 = constants % cnode(:, i)
        r1 = constants % rnode(i)
        c1 = c - c1
        call fmm_l2l_rotation_work(c1, r, r1, params % pl, &
            & constants % vscales, constants % vfact, one, &
            & node_l(:, j), one, node_l(:, i), work)
    end do
end subroutine tree_l2l_rotation_work

!------------------------------------------------------------------------------
!> Transfer local coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_l2l_bessel_rotation(params, constants, node_l)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Temporary workspace
    !real(dp) :: work(6*params % pl**2 + 19*params % pl + 8)
    ! Call corresponding work routine
    call tree_l2l_bessel_rotation_work(params, constants, node_l)
end subroutine tree_l2l_bessel_rotation

!------------------------------------------------------------------------------
!> Transfer local coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_l2l_bessel_rotation_work(params, constants, node_l)
    use complex_bessel
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pm**2 + 19*params % pm + 8)
    complex(dp) :: work_complex(2*params % pm+1)
    ! Local variables
    integer :: i, j
    real(dp) :: c_child(3), c_parent(3), c_diff(3)
    ! Top-to-bottom pass
    do i = 2, constants % nclusters
        j = constants % parent(i)
        c_child = constants % cnode(:, j)
        c_parent = constants % cnode(:, i)
        c_diff = params % kappa*(c_child - c_parent)
        call fmm_l2l_bessel_rotation_work(c_diff, &
            & constants % si_rnode(:, j), constants % si_rnode(:, i), &
            & params % pl, constants % vscales, constants % vfact, one, &
            & node_l(:, j), one, node_l(:, i), work, work_complex)
    end do
end subroutine tree_l2l_bessel_rotation_work

!------------------------------------------------------------------------------
!> Adjoint transfer local coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_l2l_rotation_adj(params, constants, node_l)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pl**2 + 19*params % pl + 8)
    ! Call corresponding work routine
    call tree_l2l_rotation_adj_work(params, constants, node_l, work)
end subroutine tree_l2l_rotation_adj

!------------------------------------------------------------------------------
!> Adjoint transfer local coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_l2l_rotation_adj_work(params, constants, node_l, work)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp), intent(out) :: work(6*params % pl**2 + 19*params % pl + 8)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Bottom-to-top pass
    do i = constants % nclusters, 1, -1
        ! Leaf node does not need any update
        if (constants % children(1, i) == 0) cycle
        c = constants % cnode(:, i)
        r = constants % rnode(i)
        ! First child initializes output
        j = constants % children(1, i)
        c1 = constants % cnode(:, j)
        r1 = constants % rnode(j)
        call fmm_l2l_rotation_adj_work(c1-c, r1, r, params % pl, &
            & constants % vscales, constants % vfact, one, &
            & node_l(:, j), zero, node_l(:, i), work)
        ! All other children update the same output
        do j = constants % children(1, i)+1, constants % children(2, i)
            c1 = constants % cnode(:, j)
            r1 = constants % rnode(j)
            call fmm_l2l_rotation_adj_work(c1-c, r1, r, params % pl, &
                & constants % vscales, constants % vfact, one, &
                & node_l(:, j), one, node_l(:, i), work)
        end do
    end do
end subroutine tree_l2l_rotation_adj_work

!------------------------------------------------------------------------------
!> Adjoint transfer local coefficients over a tree
!------------------------------------------------------------------------------
subroutine tree_l2l_bessel_rotation_adj(params, constants, node_l)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    ! Output
    real(dp), intent(inout) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pl**2 + 19*params % pl + 8)
    complex(dp) :: work_complex(2*params % pl+1)
    ! Local variables
    integer :: i, j
    real(dp) :: c_parent(3), c_child(3), c_diff(3)
    ! Bottom-to-top pass
    do i = constants % nclusters, 1, -1
        c_child = constants % cnode(:, i)
        j = constants % parent(i)
        if (j == 0) cycle
        c_parent = constants % cnode(:, j)
        c_diff = params % kappa*(c_parent - c_child)
        call fmm_l2l_bessel_rotation_adj_work(c_diff, &
            & constants % si_rnode(:, i), constants % si_rnode(:, j), &
            & params % pl, constants % vscales, constants % vfact, one, &
            & node_l(:, i), one, node_l(:, j), work, work_complex)
    end do
end subroutine tree_l2l_bessel_rotation_adj

!------------------------------------------------------------------------------
!> Transfer multipole local coefficients into local over a tree
!------------------------------------------------------------------------------
subroutine tree_m2l_rotation(params, constants, node_m, node_l)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Output
    real(dp), intent(out) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*max(params % pm, params % pl)**2 + &
        & 19*max(params % pm, params % pl) + 8)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Any order of this cycle is OK
    !$omp parallel do default(none) shared(constants,params,node_m,node_l) &
    !$omp private(i,c,r,k,c1,r1,work) schedule(dynamic)
    do i = 1, constants % nclusters
        ! If no far admissible pairs just set output to zero
        if (constants % nfar(i) .eq. 0) then
            node_l(:, i) = zero
            cycle
        end if
        c = constants % cnode(:, i)
        r = constants % rnode(i)
        ! Use the first far admissible pair to initialize output
        k = constants % far(constants % sfar(i))
        c1 = constants % cnode(:, k)
        r1 = constants % rnode(k)
        c1 = c1 - c
        call fmm_m2l_rotation_work(c1, r1, r, params % pm, params % pl, &
            & constants % vscales, constants % m2l_ztranslate_coef, one, &
            & node_m(:, k), zero, node_l(:, i), work)
        do j = constants % sfar(i)+1, constants % sfar(i+1)-1
            k = constants % far(j)
            c1 = constants % cnode(:, k)
            r1 = constants % rnode(k)
            c1 = c1 - c
            call fmm_m2l_rotation_work(c1, r1, r, params % pm, &
                & params % pl, constants % vscales, &
                & constants % m2l_ztranslate_coef, one, node_m(:, k), one, &
                & node_l(:, i), work)
        end do
    end do
end subroutine tree_m2l_rotation

!------------------------------------------------------------------------------
!> Transfer multipole local coefficients into local over a tree
!------------------------------------------------------------------------------
subroutine tree_m2l_bessel_rotation(params, constants, node_m, node_l)
    use complex_bessel
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Output
    real(dp), intent(out) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pm**2 + 19*params % pm + 8)
    complex(dp) :: work_complex(2*params % pm+1)
    ! Local variables
    integer :: i, j, k, NZ, ierr
    real(dp) :: c1(3), c(3), r1, r
    ! Any order of this cycle is OK
    !$omp parallel do default(none) shared(constants,params,node_m,node_l) &
    !$omp private(i,c,r,k,c1,r1,work,work_complex) schedule(dynamic)
    do i = 1, constants % nclusters
        ! If no far admissible pairs just set output to zero
        if (constants % nfar(i) .eq. 0) then
            node_l(:, i) = zero
            cycle
        end if
        c = constants % cnode(:, i)
        r = constants % rnode(i)
        ! Use the first far admissible pair to initialize output
        k = constants % far(constants % sfar(i))
        c1 = constants % cnode(:, k)
        r1 = constants % rnode(k)
!        call fmm_m2l_bessel_rotation(c1-c, r1, r, params % kappa, &
!            & params % pm, &
!            & constants % vscales, constants % vcnk, one, &
!            & node_m(:, k), zero, node_l(:, i))
        c1 = params % kappa*(c1 - c)
        call fmm_m2l_bessel_rotation_work(c1, &
            & constants % SK_rnode(:, k), constants % SI_rnode(:, i), &
            & params % pm, &
            & constants % vscales, constants % vcnk, one, &
            & node_m(:, k), zero, node_l(:, i), work, work_complex)
        do j = constants % sfar(i)+1, constants % sfar(i+1)-1
            k = constants % far(j)
            c1 = constants % cnode(:, k)
            r1 = constants % rnode(k)
!            call fmm_m2l_bessel_rotation(c1-c, r1, r, params % kappa, &
!                & params % pm, &
!                & constants % vscales, &
!                & constants % vcnk, one, node_m(:, k), one, &
!                & node_l(:, i))
            c1 = params % kappa*(c1 - c)
            call fmm_m2l_bessel_rotation_work(c1, constants % SK_rnode(:, k), &
                & constants % SI_rnode(:, i), params % pm, &
                & constants % vscales, constants % vcnk, one, &
                & node_m(:, k), one, node_l(:, i), work, work_complex)
        end do
    end do
end subroutine tree_m2l_bessel_rotation
!------------------------------------------------------------------------------
!> Adjoint transfer multipole local coefficients into local over a tree
!------------------------------------------------------------------------------
subroutine tree_m2l_bessel_rotation_adj(params, constants, node_l, node_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Output
    real(dp), intent(inout) :: node_m((params % pm+1)**2, constants % nclusters)
    !Call corresponding work routine
    call tree_m2l_bessel_rotation_adj_work(params, constants, node_l, node_m)
end subroutine tree_m2l_bessel_rotation_adj


!------------------------------------------------------------------------------
!> Adjoint transfer multipole local coefficients into local over a tree
!------------------------------------------------------------------------------
subroutine tree_m2l_bessel_rotation_adj_work(params, constants, node_l, node_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Output
    real(dp), intent(out) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Temporary workspace
    real(dp) :: work(6*params % pm**2 + 19*params % pm + 8)
    complex(dp) :: work_complex(2*params % pm+1)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Any order of this cycle is OK
    node_m = zero
    !!$omp parallel do shared(constants,params,node_l,node_m) &
    !!$omp private(i,c,r,k,c1,work,work_complex,j) schedule(dynamic)
    do i = 1, constants % nclusters
        ! If no far admissible pairs just set output to zero
        if (constants % nfar(i) .eq. 0) then
            cycle
        end if
        c = constants % cnode(:, i)
        r = constants % rnode(i)
        ! Use the first far admissible pair to initialize output
        k = constants % far(constants % sfar(i))
        c1 = constants % cnode(:, k)
        c1 = params % kappa*(c1 - c)
        call fmm_m2l_bessel_rotation_adj_work(c1, constants % SI_rnode(:, i), &
            & constants % SK_rnode(:, k), params % pm, &
            & constants % vscales, constants % m2l_ztranslate_adj_coef, one, &
            & node_l(:, i), one, node_m(:, k), work, work_complex)
        do j = constants % sfar(i)+1, constants % sfar(i+1)-1
            k = constants % far(j)
            c1 = constants % cnode(:, k)
            c1 = params % kappa*(c1 - c)
            call fmm_m2l_bessel_rotation_adj_work(c1, constants % SI_rnode(:, i), &
                & constants % SK_rnode(:, k), params % pm, &
                & constants % vscales, constants % m2l_ztranslate_adj_coef, one, &
                & node_l(:, i), one, node_m(:, k), work, work_complex)
        end do
    end do
end subroutine tree_m2l_bessel_rotation_adj_work

!------------------------------------------------------------------------------
!> Adjoint transfer multipole local coefficients into local over a tree
!------------------------------------------------------------------------------
subroutine tree_m2l_rotation_adj(params, constants, node_l, node_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Output
    real(dp), intent(out) :: node_m((params % pm+1)**2, constants % nclusters)
    ! Local variables
    integer :: i, j, k
    real(dp) :: c1(3), c(3), r1, r
    ! Any order of this cycle is OK
    node_m = zero
    do i = 1, constants % nclusters
        ! If no far admissible pairs just set output to zero
        if (constants % nfar(i) .eq. 0) then
            cycle
        end if
        c = constants % cnode(:, i)
        r = constants % rnode(i)
        ! Use the first far admissible pair to initialize output
        k = constants % far(constants % sfar(i))
        c1 = constants % cnode(:, k)
        r1 = constants % rnode(k)
        call fmm_m2l_rotation_adj(c-c1, r, r1, params % pl, params % pm, &
            & constants % vscales, constants % m2l_ztranslate_adj_coef, one, &
            & node_l(:, i), one, node_m(:, k))
        do j = constants % sfar(i)+1, constants % sfar(i+1)-1
            k = constants % far(j)
            c1 = constants % cnode(:, k)
            r1 = constants % rnode(k)
            call fmm_m2l_rotation_adj(c-c1, r, r1, params % pl, params % pm, &
                & constants % vscales, constants % m2l_ztranslate_adj_coef, one, &
                & node_l(:, i), one, node_m(:, k))
        end do
    end do
end subroutine tree_m2l_rotation_adj

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_l2p(params, constants, alpha, node_l, beta, grid_v, sph_l)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: node_l((params % pl+1)**2, constants % nclusters), &
        & alpha, beta
    ! Output
    real(dp), intent(inout) :: grid_v(params % ngrid, params % nsph)
    ! Scratch
    real(dp), intent(out) :: sph_l((params % pl+1)**2, params % nsph)
    ! Local variables
    real(dp) :: c(3)
    integer :: isph
    external :: dgemm


    ! Init output
    if (beta .eq. zero) then
        grid_v = zero
    else
        grid_v = beta * grid_v
    end if
    ! Get data from all clusters to spheres
    !$omp parallel do default(none) shared(params,constants,node_l,sph_l) &
    !$omp private(isph) schedule(dynamic)
    do isph = 1, params % nsph
        sph_l(:, isph) = node_l(:, constants % snode(isph))
    end do
    ! Get values at grid points
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & (params % pl+1)**2, alpha, constants % vgrid2, &
        & constants % vgrid_nbasis, sph_l, (params % pl+1)**2, beta, grid_v, &
        & params % ngrid)

end subroutine tree_l2p

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_l2p_bessel(params, constants, alpha, node_l, beta, grid_v)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: node_l((params % pl+1)**2, constants % nclusters), &
        & alpha, beta
    ! Output
    real(dp), intent(inout) :: grid_v(params % ngrid, params % nsph)
    ! Local variables
    real(dp) :: sph_l((params % pl+1)**2, params % nsph), c(3)
    integer :: isph
    external :: dgemm
    ! Init output
    if (beta .eq. zero) then
        grid_v = zero
    else
        grid_v = beta * grid_v
    end if
    ! Get data from all clusters to spheres
    !$omp parallel do default(none) shared(params,constants,node_l,sph_l) &
    !$omp private(isph) schedule(dynamic)
    do isph = 1, params % nsph
        sph_l(:, isph) = node_l(:, constants % snode(isph))
    end do
    ! Get values at grid points
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & (params % pl+1)**2, alpha, constants % vgrid, &
        & constants % vgrid_nbasis, sph_l, (params % pl+1)**2, beta, grid_v, &
        & params % ngrid)
end subroutine tree_l2p_bessel

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_l2p_adj(params, constants, alpha, grid_v, beta, node_l, sph_l)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: grid_v(params % ngrid, params % nsph), alpha, &
        & beta
    ! Output
    real(dp), intent(inout) :: node_l((params % pl+1)**2, &
        & constants % nclusters)
    ! Scractch
    real(dp), intent(out) :: sph_l((params % pl+1)**2, params % nsph)
    ! Local variables
    real(dp) :: c(3)
    integer :: isph, inode
    external :: dgemm
    ! Init output
    if (beta .eq. zero) then
        node_l = zero
    else
        node_l = beta * node_l
    end if
    ! Get weights of spherical harmonics at each sphere
    call dgemm('N', 'N', (params % pl+1)**2, params % nsph, &
        & params % ngrid, one, constants % vgrid2, constants % vgrid_nbasis, &
        & grid_v, params % ngrid, zero, sph_l, &
        & (params % pl+1)**2)
    ! Get data from all clusters to spheres
    !$omp parallel do default(none) shared(params,constants,node_l,sph_l, &
    !$omp alpha) private(isph,inode) schedule(dynamic)
    do isph = 1, params % nsph
        inode = constants % snode(isph)
        node_l(:, inode) = node_l(:, inode) + alpha*sph_l(:, isph)
    end do
end subroutine tree_l2p_adj

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_l2p_bessel_adj(params, constants, alpha, grid_v, beta, node_l)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: grid_v(params % ngrid, params % nsph), alpha, &
        & beta
    ! Output
    real(dp), intent(inout) :: node_l((params % pl+1)**2, &
        & constants % nclusters)
    ! Local variables
    real(dp) :: sph_l((params % pl+1)**2, params % nsph), c(3)
    integer :: isph, inode
    external :: dgemm
    ! Init output
    if (beta .eq. zero) then
        node_l = zero
    else
        node_l = beta * node_l
    end if
    ! Get weights of spherical harmonics at each sphere
    call dgemm('N', 'N', (params % pl+1)**2, params % nsph, &
        & params % ngrid, one, constants % vgrid, constants % vgrid_nbasis, &
        & grid_v, params % ngrid, zero, sph_l, &
        & (params % pl+1)**2)
    ! Get data from all clusters to spheres
    !$omp parallel do default(none) shared(params,constants,node_l,sph_l, &
    !$omp alpha) private(isph,inode) schedule(dynamic)
    do isph = 1, params % nsph
        inode = constants % snode(isph)
        node_l(:, inode) = node_l(:, inode) + alpha*sph_l(:, isph)
    end do
end subroutine tree_l2p_bessel_adj

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_m2p(params, constants, p, alpha, sph_m, beta, grid_v)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: p
    real(dp), intent(in) :: sph_m((p+1)**2, params % nsph), alpha, beta
    ! Output
    real(dp), intent(inout) :: grid_v(params % ngrid, params % nsph)
    ! Local variables
    integer :: isph, inode, jnear, jnode, jsph, igrid
    real(dp) :: c(3)
    ! Temporary workspace
    real(dp) :: work(p+1)
    ! Init output
    if (beta .eq. zero) then
        grid_v = zero
    else
        grid_v = beta * grid_v
    end if
    ! Cycle over all spheres
    !$omp parallel do default(none) shared(params,constants,grid_v,p, &
    !$omp alpha,sph_m), private(isph,inode,jnear,jnode,jsph,igrid,c,work) &
    !$omp schedule(dynamic)
    do isph = 1, params % nsph
        ! Cycle over all near-field admissible pairs of spheres
        inode = constants % snode(isph)
        do jnear = constants % snear(inode), constants % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = constants % near(jnear)
            jsph = constants % order(constants % cluster(1, jnode))
            ! Ignore self-interaction
            if(isph .eq. jsph) cycle
            ! Accumulate interaction for external grid points only
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
                c = constants % cgrid(:, igrid)*params % rsph(isph) - &
                    & params % csph(:, jsph) + params % csph(:, isph)
                call fmm_m2p_work(c, params % rsph(jsph), p, &
                    & constants % vscales_rel, alpha, sph_m(:, jsph), one, &
                    & grid_v(igrid, isph), work)
            end do
        end do
    end do
end subroutine tree_m2p

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_m2p_bessel(params, constants, p, alpha, sph_p, sph_m, beta, grid_v)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: p, sph_p
    real(dp), intent(in) :: sph_m((sph_p+1)**2, params % nsph), alpha, beta
    ! Output
    real(dp), intent(inout) :: grid_v(params % ngrid, params % nsph)
    ! Local variables
    integer :: isph, inode, jnear, jnode, jsph, igrid
    real(dp) :: c(3)
    ! Temporary workspace
    real(dp) :: work(p+1)
    complex(dp) :: work_complex(p+1)
    ! Init output
    if (beta .eq. zero) then
        grid_v = zero
    else
        grid_v = beta * grid_v
    end if
    ! Cycle over all spheres
    !$omp parallel do default(none) shared(params,constants,grid_v,p, &
    !$omp alpha,sph_m), private(isph,inode,jnear,jnode,jsph,igrid,c,work, &
    !$omp work_complex) schedule(dynamic)
    do isph = 1, params % nsph
        ! Cycle over all near-field admissible pairs of spheres
        inode = constants % snode(isph)
        do jnear = constants % snear(inode), constants % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = constants % near(jnear)
            jsph = constants % order(constants % cluster(1, jnode))
            ! Ignore self-interaction
            !if(isph .eq. jsph) cycle
            ! Accumulate interaction for external grid points only
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
                c = constants % cgrid(:, igrid)*params % rsph(isph) - &
                    & params % csph(:, jsph) + params % csph(:, isph)
                c = c * params % kappa
                call fmm_m2p_bessel_work(c, p, constants % vscales, &
                    & constants % SK_ri(:, jsph), alpha, sph_m(:, jsph), one, &
                    & grid_v(igrid, isph), work_complex, work)
            end do
        end do
    end do
end subroutine tree_m2p_bessel

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_m2p_adj(params, constants, p, alpha, grid_v, beta, sph_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: p
    real(dp), intent(in) :: grid_v(params % ngrid, params % nsph), alpha, &
        & beta
    ! Output
    real(dp), intent(inout) :: sph_m((p+1)**2, params % nsph)
    ! Temporary workspace
    real(dp) :: work(p+1)
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
    !!$omp parallel do default(none) shared(params,constants,grid_v,p, &
    !!$omp alpha,sph_m), private(isph,inode,jnear,jnode,jsph,igrid,c,work) &
    !!$omp schedule(dynamic)
    do isph = 1, params % nsph
        ! Cycle over all near-field admissible pairs of spheres
        inode = constants % snode(isph)
        do jnear = constants % snear(inode), constants % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = constants % near(jnear)
            jsph = constants % order(constants % cluster(1, jnode))
            ! Ignore self-interaction
            if(isph .eq. jsph) cycle
            ! Accumulate interaction for external grid points only
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
                c = constants % cgrid(:, igrid)*params % rsph(isph) - &
                    & params % csph(:, jsph) + params % csph(:, isph)
                call fmm_m2p_adj_work(c, alpha*grid_v(igrid, isph), &
                    & params % rsph(jsph), p, constants % vscales_rel, one, &
                    & sph_m(:, jsph), work)
            end do
        end do
    end do
end subroutine tree_m2p_adj

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_m2p_bessel_adj(params, constants, p, alpha, grid_v, beta, sph_p, &
        & sph_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: p, sph_p
    real(dp), intent(in) :: grid_v(params % ngrid, params % nsph), alpha, &
        & beta
    ! Output
    real(dp), intent(inout) :: sph_m((sph_p+1)**2, params % nsph)
    ! Temporary workspace
    real(dp) :: work(p+1)
    ! Local variables
    !real(dp), allocatable :: tmp(:, :, :)
    integer :: isph, inode, jnear, jnode, jsph, igrid, info, iproc
    real(dp) :: c(3), w
    ! Init output
    if (beta .eq. zero) then
        sph_m = zero
    else
        sph_m = beta * sph_m
    end if
    ! Cycle over all spheres
    !!$omp parallel do default(none) shared(params,constants,p,sph_m,tmp, &
    !!$omp alpha,grid_v) private(isph,inode,jnear,jnode,jsph,igrid,c, &
    !!$omp iproc) schedule(dynamic)
    do isph = 1, params % nsph
        ! Cycle over all near-field admissible pairs of spheres
        inode = constants % snode(isph)
        ! iproc = omp_get_thread_num() + 1
        do jnear = constants % snear(inode), constants % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = constants % near(jnear)
            jsph = constants % order(constants % cluster(1, jnode))
            ! Ignore self-interaction
            !if(isph .eq. jsph) cycle
            ! Accumulate interaction for external grid points only
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
                c = constants % cgrid(:, igrid)*params % rsph(isph) - &
                    & params % csph(:, jsph) + params % csph(:, isph)
                call fmm_m2p_bessel_adj(c, alpha*grid_v(igrid, isph), &
                    & params % rsph(jsph), params % kappa, p, &
                    & constants % vscales, one, sph_m(:, jsph))
            end do
        end do
    end do
end subroutine tree_m2p_bessel_adj

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_m2p_bessel_nodiag_adj(params, constants, p, alpha, grid_v, beta, sph_p, &
        & sph_m)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: p, sph_p
    real(dp), intent(in) :: grid_v(params % ngrid, params % nsph), alpha, &
        & beta
    ! Output
    real(dp), intent(inout) :: sph_m((sph_p+1)**2, params % nsph)
    ! Temporary workspace
    real(dp) :: work(p+1)
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
    do isph = 1, params % nsph
        ! Cycle over all near-field admissible pairs of spheres
        inode = constants % snode(isph)
        do jnear = constants % snear(inode), constants % snear(inode+1)-1
            ! Near-field interactions are possible only between leaf nodes,
            ! which must contain only a single input sphere
            jnode = constants % near(jnear)
            jsph = constants % order(constants % cluster(1, jnode))
            ! Ignore self-interaction
            if(isph .eq. jsph) cycle
            ! Accumulate interaction for external grid points only
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
                c = constants % cgrid(:, igrid)*params % rsph(isph) - &
                    & params % csph(:, jsph) + params % csph(:, isph)
                call fmm_m2p_bessel_adj(c, alpha*grid_v(igrid, isph), &
                    & params % rsph(jsph), params % kappa, p, constants % vscales, one, &
                    & sph_m(:, jsph))
            end do
        end do
    end do
end subroutine tree_m2p_bessel_nodiag_adj

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine efld(ncav,zeta,ccav,nsph,csph,force)
integer,                    intent(in)    :: ncav, nsph
real*8,  dimension(ncav),   intent(in)    :: zeta
real*8,  dimension(3,ncav), intent(in)    :: ccav
real*8,  dimension(3,nsph), intent(in)    :: csph
real*8,  dimension(3,nsph), intent(inout) :: force
!
integer :: icav, isph
real*8  :: dx, dy, dz, rijn2, rijn, rijn3, f, sum_int(3)
real*8, parameter :: zero=0.0d0
!
force = zero
do isph = 1, nsph
  sum_int = zero
  do icav = 1, ncav
    dx   = csph(1,isph) - ccav(1,icav)
    dy   = csph(2,isph) - ccav(2,icav)
    dz   = csph(3,isph) - ccav(3,icav)
    rijn2   = dx*dx + dy*dy + dz*dz
    rijn   = sqrt(rijn2)
    rijn3   = rijn2*rijn
    f    = zeta(icav)/rijn3
    sum_int(1) = sum_int(1) + f*dx
    sum_int(2) = sum_int(2) + f*dy
    sum_int(3) = sum_int(3) + f*dz
  end do
  force(:,isph) = sum_int
end do
end subroutine efld

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_grad_m2m(params, constants, sph_m, sph_m_grad, work)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: sph_m(constants % nbasis, params % nsph)
    ! Output
    real(dp), intent(inout) :: sph_m_grad((params % lmax+2)**2, 3, &
        & params % nsph)
    ! Temporary workspace
    real(dp), intent(inout) :: work((params % lmax+2)**2, params % nsph)
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
    sph_m_grad(1:constants % nbasis, 3, :) = sph_m
    do isph = 1, params % nsph
        call fmm_sph_transform(params % lmax, zx_coord_transform, one, &
            & sph_m(:, isph), zero, sph_m_grad(1:constants % nbasis, 1, isph))
        call fmm_sph_transform(params % lmax, zy_coord_transform, one, &
            & sph_m(:, isph), zero, sph_m_grad(1:constants % nbasis, 2, isph))
    end do
    ! Derivative of M2M translation over OZ axis at the origin consists of 2
    ! steps:
    !   1) increase degree l and scale by sqrt((2*l+1)*(l*l-m*m)) / sqrt(2*l-1)
    !   2) scale by 1/rsph(isph)
    do l = params % lmax+1, 1, -1
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
    do isph = 1, params % nsph
        sph_m_grad(:, 3, isph) = sph_m_grad(:, 3, isph) / params % rsph(isph)
        work(:, isph) = sph_m_grad(:, 1, isph) / params % rsph(isph)
        call fmm_sph_transform(params % lmax+1, zx_coord_transform, one, &
            & work(:, isph), zero, sph_m_grad(:, 1, isph))
        work(:, isph) = sph_m_grad(:, 2, isph) / params % rsph(isph)
        call fmm_sph_transform(params % lmax+1, zy_coord_transform, one, &
            & work(:, isph), zero, sph_m_grad(:, 2, isph))
    end do
end subroutine tree_grad_m2m

!------------------------------------------------------------------------------
!> TODO
!------------------------------------------------------------------------------
subroutine tree_grad_l2l(params, constants, node_l, sph_l_grad, work)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: node_l((params % pl+1)**2, constants % nclusters)
    ! Output
    real(dp), intent(out) :: sph_l_grad((params % pl+1)**2, 3, params % nsph)
    ! Temporary workspace
    real(dp), intent(out) :: work((params % pl+1)**2, params % nsph)
    ! Local variables
    integer :: isph, inode, l, indi, indj, m
    real(dp) :: tmp1, tmp2
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    ! Gradient of L2L reduces degree by 1, so exit if degree of harmonics is 0
    ! or -1 (which means no FMM at all)
    if (params % pl .le. 0) return
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
    do isph = 1, params % nsph
        inode = constants % snode(isph)
        sph_l_grad(:, 3, isph) = node_l(:, inode)
        call fmm_sph_transform(params % pl, zx_coord_transform, one, &
            & node_l(:, inode), zero, sph_l_grad(:, 1, isph))
        call fmm_sph_transform(params % pl, zy_coord_transform, one, &
            & node_l(:, inode), zero, sph_l_grad(:, 2, isph))
    end do
    ! Derivative of L2L translation over OZ axis at the origin consists of 2
    ! steps:
    !   1) decrease degree l and scale by sqrt((2*l-1)*(l*l-m*m)) / sqrt(2*l+1)
    !   2) scale by 1/rsph(isph)
    do l = 1, params % pl
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
    do isph = 1, params % nsph
        sph_l_grad(1:params % pl**2, 3, isph) = &
            & sph_l_grad(1:params % pl**2, 3, isph) / params % rsph(isph)
        work(1:params % pl**2, isph) = &
            & sph_l_grad(1:params % pl**2, 1, isph) / params % rsph(isph)
        call fmm_sph_transform(params % pl-1, zx_coord_transform, one, &
            & work(1:params % pl**2, isph), zero, &
            & sph_l_grad(1:params % pl**2, 1, isph))
        work(1:params % pl**2, isph) = &
            & sph_l_grad(1:params % pl**2, 2, isph) / params % rsph(isph)
        call fmm_sph_transform(params % pl-1, zy_coord_transform, one, &
            & work(1:params % pl**2, isph), zero, &
            & sph_l_grad(1:params % pl**2, 2, isph))
    end do
    ! Set degree pl to zero to avoid problems if user actually uses it
    l = params % pl
    indi = l*l + l + 1
    sph_l_grad(indi-l:indi+l, :, :) = zero
end subroutine tree_grad_l2l

!> Print the ddX logo
!> @ingroup Fortran_interface_core
!!
!! @param[out] string: container for the logo
!!
subroutine get_banner(string)
    implicit none
    character (len=2047), intent(out) :: string
    character (len=10) :: vstr
    write(vstr, *) "0.1.3"
    write(string, *) &
        & " +----------------------------------------------------------------+", &
        & NEW_LINE('a'), &
        & "  |                                                                |", &
        & NEW_LINE('a'), &
        & "  |                        888      888 Y8b    d8Y                 |", &
        & NEW_LINE('a'), &
        & "  |                        888      888  Y8b  d8Y                  |", &
        & NEW_LINE('a'), &
        & "  |                        888      888   Y8888Y                   |", &
        & NEW_LINE('a'), &
        & "  |                    .d88888  .d88888    Y88Y                    |", &
        & NEW_LINE('a'), &
        & "  |                   d88  888 d88  888    d88b                    |", &
        & NEW_LINE('a'), &
        & "  |                   888  888 888  888   d8888b                   |", &
        & NEW_LINE('a'), &
        & "  |                   Y88b 888 Y88b 888  d8Y  Y8b                  |", &
        & NEW_LINE('a'), &
        & "  |                     Y88888   Y88888 d8Y    Y8b                 |", &
        & NEW_LINE('a'), &
        & "  |                                                                |", &
        & NEW_LINE('a'), &
        & "  |                 https://ddsolvation.github.io/ddX/             |", &
        & NEW_LINE('a'), &
        & "  |                         Version:", vstr, "                     |", &
        & NEW_LINE('a'), &
        & "  |                                                                |", &
        & NEW_LINE('a'), &
        & "  +----------------------------------------------------------------+"
end subroutine get_banner

!> Transform a function defined at the exposed cavity points (cav) to
!> a spherical harmonics expansion. Note that the function is also
!> multiplied by the characteristic function U.
!!
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] property_cav: property defined at the exposed cavity points,
!!      size (ncav)
!! @param[out] property_sph: property as a spherical harmonics expansion,
!!      size (nbasis, nsph)
subroutine cav_to_spherical(params, constants, workspace, property_cav, &
        & property_sph)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: property_cav(constants % ncav)
    real(dp), intent(out) :: property_sph(constants % nbasis, params % nsph)

    ! multiply by the characteristic function U
    workspace % tmp_cav = property_cav * constants % ui_cav

    ! extend the function to the sphere intersection with zeros
    call ddcav_to_grid_work(params % ngrid, params % nsph, constants % ncav, &
        & constants % icav_ia, constants % icav_ja, workspace % tmp_cav, &
        & workspace % tmp_grid)

    ! integrate against spherical harmonics
    call ddintegrate(params % nsph, constants % nbasis, &
        & params % ngrid, constants % vwgrid, constants % vgrid_nbasis, &
        & workspace % tmp_grid, property_sph)
end subroutine cav_to_spherical

end module ddx_core
