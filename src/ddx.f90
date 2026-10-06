!> High-level module of the ddX software
module ddx
! Get ddcosmo-module
use ddx_cosmo
! Get ddpcm-module
use ddx_pcm
! Get ddlpb-module
use ddx_lpb

use ddx_defaults
implicit none

contains

!> @defgroup Fortran_interface_core Fortran interface: core routines

!> Wrapper to the initialization routine which supports optional arguments.
!! This will make future maintenance easier.
!!
!> @ingroup Fortran_interface_core
!! @param[in] model: 1 for COSMO, 2 for PCM and 3 for LPB
!! @param[in] nsph: number of spheres n > 0
!! @param[in] coords: coordinates of the spheres, size (3, nsph)
!! @param[in] radii: radii of the spheres, size (nsph)
!! @param[in] eps: relative dielectric permittivity eps > 1
!! @param[out] ddx_data: Container for ddX data structures
!! @param[inout] ddx_error: ddX error
!!
!! @param[in,optional] force: 1 if forces are required and 0 otherwise
!! @param[in,optional] kappa: Debye-H\"{u}ckel parameter
!! @param[in,optional] eta: regularization parameter 0 < eta <= 1
!! @param[in,optional] shift: shift of characteristic function -1 for interior,
!!                     0 for centered and 1 for outer regularization
!! @param[in,optional] lmax: Maximal degree of modeling spherical harmonics,
!!                     `lmax` >= 0
!! @param[in,optional] ngrid: Number of Lebedev grid points `ngrid` >= 0
!! @param[in,optional] incore: handling of the sparse matrices, 1 for
!!                     precomputing them and keeping them in memory,
!!                     0 for assembling the matrix-vector products on-the-fly
!! @param[in,optional] maxiter: Maximum number of iterations for the
!!                     iterative solver,  maxiter > 0
!! @param[in,optional] jacobi_ndiis: Number of extrapolation points for
!!                     the  Jacobi/DIIS solver, ndiis >= 0
!! @param[in,optional] enable_fmm: 1 to use FMM acceleration and 0 otherwise
!! @param[in,optional] pm: Maximal degree of multipole spherical harmonics.
!!                     Ignored in the case `fmm=0`. Value -1 means no
!!                     far-field FMM interactions are computed, `pm` >= -1
!! @param[in,optional] pl: Maximal degree of local spherical harmonics.
!!                     Ignored in the case `fmm=0`. Value -1 means no
!!                     far-field FMM interactions are computed, `pl` >= -1
!! @param[in,optional] nproc: Number of OpenMP threads, nproc >= 0.
!! @param[in,optional] logfile: file name for log information.
!! @param[in,optional] switching: kind of switching, 0 legacy, 1 new version.
!!
subroutine ddinit(model, nsph, coords, radii, eps, ddx_data, ddx_error, &
        & force, kappa, eta, shift, lmax, ngrid, incore, maxiter, &
        & jacobi_ndiis, enable_fmm, pm, pl, nproc, logfile, adjoint, &
        & eps_int, switching)

    ! mandatory arguments
    integer, intent(in) :: model, nsph
    real(dp), intent(in) :: coords(3, nsph), radii(nsph)
    type(ddx_type), target, intent(out) :: ddx_data
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: eps

    ! optional arguments
    integer, intent(in), optional :: force, adjoint, lmax, ngrid, incore, &
        & maxiter, jacobi_ndiis, enable_fmm, pm, pl, nproc, switching
    real(dp), intent(in), optional :: kappa, eta, shift, eps_int
    character(len=255), intent(in), optional :: logfile

    ! local copies of the optional arguments with default values
    integer :: local_force = default_force
    integer :: local_adjoint = default_adjoint
    integer :: local_lmax = default_lmax
    integer :: local_ngrid = default_ngrid
    integer :: local_incore = default_incore
    integer :: local_maxiter = default_maxiter
    integer :: local_jacobi_ndiis = default_jacobi_ndiis
    integer :: local_enable_fmm = default_enable_fmm
    integer :: local_pm = default_pm
    integer :: local_pl = default_pl
    integer :: local_nproc = default_nproc
    integer :: local_switching = default_switching

    real(dp) :: local_kappa = default_kappa
    real(dp) :: local_eta = default_eta
    real(dp) :: local_shift
    real(dp) :: local_eps_int = default_eps_int

    character(len=255) :: local_logfile = default_logfile

    ! arrays for x, y, z coordinates
    real(dp), allocatable :: x(:), y(:), z(:)
    integer :: info

    ! The default switching shift depends on the model
    if (model.lt.1 .or. model.gt.3) then
        call update_error(ddx_error, "ddinit: wrong value of `model`")
        return
    end if
    local_shift = default_se(model)

    ! Update local variables with provided optional arguments
    if (present(force)) local_force = force
    if (present(lmax)) local_lmax = lmax
    if (present(ngrid)) local_ngrid = ngrid
    if (present(incore)) local_incore = incore
    if (present(maxiter)) local_maxiter = maxiter
    if (present(jacobi_ndiis)) local_jacobi_ndiis = jacobi_ndiis
    if (present(enable_fmm)) local_enable_fmm = enable_fmm
    if (present(pm)) local_pm = pm
    if (present(pl)) local_pl = pl
    if (present(nproc)) local_nproc = nproc
    if (present(kappa)) local_kappa = kappa
    if (present(eta)) local_eta = eta
    if (present(logfile)) local_logfile = logfile
    if (present(switching)) local_switching = switching
    if (present(shift)) local_shift = shift

    ! this are not yet supported, but they will probably
    if (present(adjoint)) then
        local_adjoint = adjoint
        call update_error(ddx_error, &
            & "ddinit: adjoint argument is not yet supported")
        return
    end if
    if (present(eps_int)) then
        local_eps_int = eps_int
        call update_error(ddx_error, &
            & "ddinit: eps_int argument is not yet supported")
        return
    end if

    allocate(x(nsph), y(nsph), z(nsph), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "ddinit: allocation failed")
        return
    end if
    x = coords(1,:)
    y = coords(2,:)
    z = coords(3,:)

    call allocate_model(nsph, x, y, z, radii, model, local_lmax, local_ngrid, &
        & local_force, local_enable_fmm, local_pm, local_pl, local_shift, &
        & local_eta, eps, local_kappa, local_incore, local_maxiter, &
        & local_jacobi_ndiis, local_nproc, local_logfile, local_switching, &
        & ddx_data, ddx_error)

    deallocate(x, y, z, stat=info)
    if (info.ne.0) then
        call update_error(ddx_error, "ddinit: deallocation failed")
        return
    end if
end subroutine ddinit

!> Main solver routine
!!
!! Solves the solvation problem, computes the energy, and if required
!! computes the forces.
!!
!> @ingroup Fortran_interface_core
!! @param[in] ddx_data: ddX object with all input information
!! @param[inout] state: ddx state (contains RHSs and solutions)
!! @param[in] electrostatics: electrostatic property container
!! @param[in] psi: RHS of the adjoint problem
!! @param[in] tol: tolerance for the linear system solvers
!! @param[out] esolv: solvation energy
!! @param[inout] ddx_error: ddX error
!! @param[out] force: Analytical forces (optional argument, only if
!!             required)
!! @param[in] read_guess: optional argument, if true read the guess
!!            from the state object
!!
subroutine ddrun(ddx_data, state, electrostatics, psi, tol, esolv, &
        & ddx_error, force, read_guess)
    type(ddx_type), intent(inout) :: ddx_data
    type(ddx_state_type), intent(inout) :: state
    type(ddx_electrostatics_type), intent(in) :: electrostatics
    real(dp), intent(in) :: psi(ddx_data % constants % nbasis, &
        & ddx_data % params % nsph)
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: esolv
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(out), optional :: force(3, ddx_data % params % nsph)
    logical, intent(in), optional :: read_guess
    ! local variables
    logical :: do_guess

    ! decide if the guess has to be read or must be done
    if (present(read_guess)) then
        do_guess = .not.read_guess
    else
        do_guess = .true.
    end if

    ! if the forces are to be computed, but the array is not passed, raise
    ! an error
    if ((.not.present(force)) .and. (ddx_data % params % force .eq. 1)) then
        call update_error(ddx_error, &
            & "ddrun: forces are to be computed, but the optional force" // &
            & " array has not been passed.")
        return
    end if

    call setup(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, electrostatics, psi, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddrun: setup returned an error, exiting")
        return
    end if

    ! solve the primal linear system
    if (do_guess) then
        call fill_guess(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, state, tol, ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddrun: fill_guess returned an error, exiting")
            return
        end if
    end if
    call solve(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, tol, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddrun: solve returned an error, exiting")
        return
    end if

    ! compute the energy
    call energy(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, esolv, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddrun: energy returned an error, exiting")
        return
    end if

    ! solve the primal linear system
    if (ddx_data % params % force .eq. 1) then
        if (do_guess) then
            call fill_guess_adjoint(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, tol, ddx_error)
            if (ddx_error % flag .ne. 0) then
                call update_error(ddx_error, &
                    & "ddrun: fill_guess_adjoint returned an error, exiting")
                return
            end if
        end if
        call solve_adjoint(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, state, tol, ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddrun: solve_adjoint returned an error, exiting")
            return
        end if
    end if

    ! compute the forces
    if (ddx_data % params % force .eq. 1) then
        force = zero
        call solvation_force_terms(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, state, electrostatics, force, ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddrun: solvation_force_terms returned an error, exiting")
            return
        end if
    end if

end subroutine ddrun

!> Setup the state for the different models
!!
!> @ingroup Fortran_interface_core
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] electrostatics: electrostatic property container
!! @param[in] psi: Representation of the solute potential in spherical
!!     harmonics, size (nbasis, nsph)
!! @param[inout] ddx_error: ddX error
!!
subroutine setup(params, constants, workspace, state, electrostatics, &
        & psi, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_electrostatics_type), intent(in) :: electrostatics
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error

    if (params % model .eq. 1) then
        call cosmo_setup(params, constants, workspace, state, &
            & electrostatics % phi_cav, psi, ddx_error)
    else if (params % model .eq. 2) then
        call pcm_setup(params, constants, workspace, state, &
            & electrostatics % phi_cav, psi, ddx_error)
    else if (params % model .eq. 3) then
        call lpb_setup(params, constants, workspace, state, &
            & electrostatics % phi_cav, electrostatics % e_cav, &
            & psi, ddx_error)
    else
        call update_error(ddx_error, "Unknow model in setup.")
        return
    end if

end subroutine setup

!> Do a guess for the primal linear system for the different models
!!
!> @ingroup Fortran_interface_core
!! @param[in] params: User specified parameters
!! @param[inout] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: tolerance
!! @param[inout] ddx_error: ddX error
!!
subroutine fill_guess(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error

    if (.not.state % rhs_done) then
        call update_error(ddx_error, &
            & "In fill_guess, the RHS is not initialized")
        return
    end if

    if (params % model .eq. 1) then
        call cosmo_guess(params, constants, workspace, state, ddx_error)
    else if (params % model .eq. 2) then
        call pcm_guess(params, constants, workspace, state, ddx_error)
    else if (params % model .eq. 3) then
        call lpb_guess(params, constants, workspace, state, tol, ddx_error)
    else
        call update_error(ddx_error, "Unknow model in fill_guess.")
        return
    end if

    if (ddx_error % flag .eq. 0) then
        state % guess_done = .true.
        state % solved = .false.
    end if

end subroutine fill_guess

!> Do a guess for the adjoint linear system for the different models
!!
!> @ingroup Fortran_interface_core
!! @param[in] params: User specified parameters
!! @param[inout] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] tol: tolerance
!! @param[inout] ddx_error: ddX error
!!
subroutine fill_guess_adjoint(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error

    if (.not.state % adjoint_rhs_done) then
        call update_error(ddx_error, &
            & "In fill_guess_adjoint, the adjoint RHS is not initialized")
        return
    end if

    if (params % model .eq. 1) then
        call cosmo_guess_adjoint(params, constants, workspace, state, ddx_error)
    else if (params % model .eq. 2) then
        call pcm_guess_adjoint(params, constants, workspace, state, ddx_error)
    else if (params % model .eq. 3) then
        call lpb_guess_adjoint(params, constants, workspace, state, tol, ddx_error)
    else
        call update_error(ddx_error, "Unknow model in fill_guess_adjoint.")
        return
    end if

    if (ddx_error % flag .eq. 0) then
        state % adjoint_guess_done = .true.
        state % adjoint_solved = .false.
    end if

end subroutine fill_guess_adjoint

!> Solve the primal linear system for the different models
!!
!> @ingroup Fortran_interface_core
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!! @param[inout] ddx_error: ddX error
!!
subroutine solve(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error

    if (.not.state % rhs_done) then
        call update_error(ddx_error, "In solve, the RHS is not initialized")
        return
    end if
    if (.not.(state % solved .or. state % guess_done)) then
        call update_error(ddx_error, &
            & "In solve, no guess or previous solution is provided")
        return
    end if

    if (params % model .eq. 1) then
        call cosmo_solve(params, constants, workspace, state, tol, ddx_error)
    else if (params % model .eq. 2) then
        call pcm_solve(params, constants, workspace, state, tol, ddx_error)
    else if (params % model .eq. 3) then
        call lpb_solve(params, constants, workspace, state, tol, ddx_error)
    else
        call update_error(ddx_error, "Unknow model in solve.")
        return
    end if

    if (ddx_error % flag .eq. 0) then
        state % solved = .true.
    end if

end subroutine solve

!> Solve the adjoint linear system for the different models
!!
!> @ingroup Fortran_interface_core
!! @param[in] params       : General options
!! @param[in] constants    : Precomputed constants
!! @param[inout] workspace : Preallocated workspaces
!! @param[inout] state     : Solutions, guesses and relevant quantities
!! @param[in] tol          : Tolerance for the iterative solvers
!! @param[inout] ddx_error: ddX error
!!
subroutine solve_adjoint(params, constants, workspace, state, tol, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: tol
    type(ddx_error_type), intent(inout) :: ddx_error

    if (.not.state % adjoint_rhs_done) then
        call update_error(ddx_error, &
            & "In solve_adjoint, the RHS is not initialized")
        return
    end if
    if (.not.(state % adjoint_solved .or. state % adjoint_guess_done)) then
        call update_error(ddx_error, &
            & "In solve_adjoint, no guess or previous solution is provided")
        return
    end if

    if (params % model .eq. 1) then
        call cosmo_solve_adjoint(params, constants, workspace, state, tol, ddx_error)
    else if (params % model .eq. 2) then
        call pcm_solve_adjoint(params, constants, workspace, state, tol, ddx_error)
    else if (params % model .eq. 3) then
        call lpb_solve_adjoint(params, constants, workspace, state, tol, ddx_error)
    else
        call update_error(ddx_error, "Unknow model in solve_adjoint.")
        return
    end if

    if (ddx_error % flag .eq. 0) then
        state % adjoint_solved = .true.
    end if

end subroutine solve_adjoint

!> Compute the energy for the different models
!!
!> @ingroup Fortran_interface_core
!! @param[in] params: General options
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[in] state: ddx state (contains solutions and RHSs)
!! @param[out] solvation_energy: resulting energy
!! @param[inout] ddx_error: ddX error
!!
subroutine energy(params, constants, workspace, state, solvation_energy, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(in) :: workspace
    type(ddx_state_type), intent(in) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(out) :: solvation_energy

    ! dummy operation on unused interface arguments
    if (allocated(workspace % tmp_pot)) continue

    if (.not.state % solved) then
        call update_error(ddx_error, &
            & "In energy, the solution is not available.")
        return
    end if
    if (.not.state % adjoint_rhs_done) then
        call update_error(ddx_error, &
            & "In energy, the adjoint RHS is not initialized")
        return
    end if


    if (params % model .eq. 1) then
        call cosmo_energy(constants, state, solvation_energy, ddx_error)
    else if (params % model .eq. 2) then
        call pcm_energy(constants, state, solvation_energy, ddx_error)
    else if (params % model .eq. 3) then
        call lpb_energy(constants, state, solvation_energy, ddx_error)
    else
        call update_error(ddx_error, "Unknow model in energy.")
        return
    end if

end subroutine energy

!> Compute the solvation terms of the forces (solute aspecific) for the
!! different models. This must be summed to the solute specific term to get
!! the full forces
!!
!> @ingroup Fortran_interface_core
!! @param[in] params: General options
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: Solutions and relevant quantities
!! @param[in] electrostatics: Electrostatic properties container.
!! @param[out] force: Geometrical contribution to the forces
!! @param[inout] ddx_error: ddX error
!!
subroutine solvation_force_terms(params, constants, workspace, &
        & state, electrostatics, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_electrostatics_type), intent(in) :: electrostatics
    real(dp), intent(out) :: force(3, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error

    if (.not.state % solved) then
        call update_error(ddx_error, &
            & "In solvation_force_terms, the solution is not available.")
        return
    end if
    if (.not.state % adjoint_solved) then
        call update_error(ddx_error, &
            & "In solvation_force_terms, the adjoint solution is not available.")
        return
    end if
    if (.not.state % rhs_done) then
        call update_error(ddx_error, &
            & "In solvation_force_terms, the RHS is not initialized")
        return
    end if
    if (.not.state % adjoint_rhs_done) then
        call update_error(ddx_error, &
            & "In solvation_force_terms, the adjoint RHS is not initialized")
        return
    end if

    if (params % model .eq. 1) then
        call cosmo_solvation_force_terms(params, constants, workspace, &
            & state, electrostatics % e_cav, force, ddx_error)
    else if (params % model .eq. 2) then
        call pcm_solvation_force_terms(params, constants, workspace, &
            & state, electrostatics % e_cav, force, ddx_error)
    else if (params % model .eq. 3) then
        call lpb_solvation_force_terms(params, constants, workspace, &
            & state, electrostatics % g_cav, force, ddx_error)
    else
        call update_error(ddx_error, "Unknow model in solvation_force_terms.")
        return
    end if

end subroutine solvation_force_terms

end module ddx
