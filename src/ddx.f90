!> High-level module of the ddX software
module ddx
! Get ddcosmo-module
use ddx_cosmo
! Get ddpcm-module
use ddx_pcm
! Get ddlpb-module
use ddx_lpb
implicit none

contains

!> @defgroup Fortran_interface_core Fortran interface: core routines

!> Read the configuration from ddX input file and return a ddx_data
!! structure
!!
!> @ingroup Fortran_interface_core
!! @param[in] fname: Filename containing all the required info
!! @param[out] ddx_data: Object containing all inputs
!! @param[out] tol: tolerance for iterative solvers
!! @param[out] charges: charge array, size(nsph)
!! @param[inout] ddx_error: ddX error
!!
subroutine ddfromfile(fname, ddx_data, tol, charges, ddx_error)
    implicit none
    character(len=*), intent(in) :: fname
    type(ddx_type), intent(out) :: ddx_data
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(out) :: tol
    real(dp), allocatable, intent(out) :: charges(:)
    ! Local variables
    integer :: nproc, model, lmax, ngrid, force, fmm, pm, pl, &
        & nsph, i, matvecmem, maxiter, jacobi_ndiis, &
        & istatus
    real(dp) :: eps, se, eta, kappa
    real(dp), allocatable :: csph(:, :), radii(:)
    character(len=255) :: output_filename
    !! Read all the parameters from the file
    ! Open a configuration file
    open(unit=100, file=fname, form='formatted', access='sequential')
    ! Printing flag
    read(100, *) output_filename
    ! Number of OpenMP threads to be used
    read(100, *) nproc
    if(nproc .lt. 0) then
        call update_error(ddx_error, "Error on the 2nd line of a config " // &
            & "file " // trim(fname) // ": `nproc` must be a positive " // &
            & "integer value.")
    end if
    ! Model to be used: 1 for COSMO, 2 for PCM and 3 for LPB
    read(100, *) model
    if((model .lt. 1) .or. (model .gt. 3)) then
        call update_error(ddx_error, "Error on the 3rd line of a config file " // &
            & trim(fname) // ": `model` must be an integer of a value " // &
            & "1, 2 or 3.")
    end if
    ! Max degree of modeling spherical harmonics
    read(100, *) lmax
    if(lmax .lt. 0) then
        call update_error(ddx_error, "Error on the 4th line of a config file " // &
            & trim(fname) // ": `lmax` must be a non-negative integer value.")
    end if
    ! Approximate number of Lebedev points
    read(100, *) ngrid
    if(ngrid .lt. 0) then
        call update_error(ddx_error, "Error on the 5th line of a config file " // &
            & trim(fname) // ": `ngrid` must be a non-negative integer value.")
    end if
    ! Dielectric permittivity constant of the solvent
    read(100, *) eps
    if(eps .lt. zero) then
        call update_error(ddx_error, "Error on the 6th line of a config file " // &
            & trim(fname) // ": `eps` must be a non-negative floating " // &
            & "point value.")
    end if
    ! Shift of the regularized characteristic function
    read(100, *) se
    if((se .lt. -one) .or. (se .gt. one)) then
        call update_error(ddx_error, "Error on the 7th line of a config file " // &
            & trim(fname) // ": `se` must be a floating point value in a " // &
            & " range [-1, 1].")
    end if
    ! Regularization parameter
    read(100, *) eta
    if((eta .lt. zero) .or. (eta .gt. one)) then
        call update_error(ddx_error, "Error on the 8th line of a config file " // &
            & trim(fname) // ": `eta` must be a floating point value " // &
            & "in a range [0, 1].")
    end if
    ! Debye H\"{u}ckel parameter
    read(100, *) kappa
    if(kappa .lt. zero) then
        call update_error(ddx_error, "Error on the 9th line of a config file " // &
            & trim(fname) // ": `kappa` must be a non-negative floating " // &
            & "point value.")
    end if
    ! whether the (sparse) matrices are precomputed and kept in memory (1)
    ! or not (0).
    read(100, *) matvecmem 
    if((matvecmem.lt. 0) .or. (matvecmem .gt. 1)) then
        call update_error(ddx_error, "Error on the 10th line of a config " // &
            & "file " // trim(fname) // ": `matvecmem` must be an " // &
            & "integer value of a value 0 or 1.")
    end if
    ! Relative convergence threshold for the iterative solver
    read(100, *) tol
    if((tol .lt. 1d-14) .or. (tol .gt. one)) then
        call update_error(ddx_error, "Error on the 12th line of a config " // &
            & "file " // trim(fname) // ": `tol` must be a floating " // &
            & "point value in a range [1d-14, 1].")
    end if
    ! Maximum number of iterations for the iterative solver
    read(100, *) maxiter
    if((maxiter .le. 0)) then
        call update_error(ddx_error, "Error on the 13th line of a config " // &
            & "file " // trim(fname) // ": `maxiter` must be a positive " // &
            & " integer value.")
    end if
    ! Number of extrapolation points for Jacobi/DIIS solver
    read(100, *) jacobi_ndiis
    if((jacobi_ndiis .lt. 0)) then
        call update_error(ddx_error, "Error on the 14th line of a config " // &
            & "file " // trim(fname) // ": `jacobi_ndiis` must be a " // &
            & "non-negative integer value.")
    end if
    ! Whether to compute (1) or not (0) forces as analytical gradients
    read(100, *) force
    if((force .lt. 0) .or. (force .gt. 1)) then
        call update_error(ddx_error, "Error on the 17th line of a config " // &
            & "file " // trim(fname) // ": `force` must be an integer " // &
            "value of a value 0 or 1.")
    end if
    ! Whether to use (1) or not (0) the FMM to accelerate computations
    read(100, *) fmm
    if((fmm .lt. 0) .or. (fmm .gt. 1)) then
        call update_error(ddx_error, "Error on the 18th line of a config " // &
            & "file " // trim(fname) // ": `fmm` must be an integer " // &
            & "value of a value 0 or 1.")
    end if
    ! Max degree of multipole spherical harmonics for the FMM
    read(100, *) pm
    if(pm .lt. 0) then
        call update_error(ddx_error, "Error on the 19th line of a config " // &
            & "file " // trim(fname) // ": `pm` must be a non-negative " // &
            & "integer value.")
    end if
    ! Max degree of local spherical harmonics for the FMM
    read(100, *) pl
    if(pl .lt. 0) then
        call update_error(ddx_error, "Error on the 20th line of a config " // &
            & "file " // trim(fname) // ": `pl` must be a non-negative " // &
            & "integer value.")
    end if
    ! Number of input spheres
    read(100, *) nsph
    if(nsph .le. 0) then
        call update_error(ddx_error, "Error on the 21th line of a config " // &
            & "file " // trim(fname) // ": `nsph` must be a positive " // &
            & "integer value.")
    end if

    ! return in case of errors in the parameters
    if (ddx_error % flag .ne. 0) return

    ! Coordinates, radii and charges
    allocate(charges(nsph), csph(3, nsph), radii(nsph), stat=istatus)
    if(istatus .ne. 0) then
        call update_error(ddx_error, "Could not allocate space for " // &
            & "coordinates, radii and charges of atoms.")
        return
    end if
    do i = 1, nsph
        read(100, *) charges(i), csph(1, i), csph(2, i), csph(3, i), radii(i)
    end do
    ! Finish reading
    close(100)
    !! Convert Angstrom input into Bohr units
    csph = csph * tobohr
    radii = radii * tobohr
    kappa = kappa / tobohr

    ! adjust ngrid
    call closest_supported_lebedev_grid(ngrid)

    !! Initialize ddx_data object
    call ddinit(model, nsph, csph, radii, eps, ddx_data, ddx_error, &
        & force=force, kappa=kappa, eta=eta, shift=se, lmax=lmax, &
        & ngrid=ngrid, incore=matvecmem, maxiter=maxiter, &
        & jacobi_ndiis=jacobi_ndiis, enable_fmm=fmm, pm=pm, pl=pl, &
        & nproc=nproc, logfile=output_filename)

    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddinit returned an error, exiting")
        return
    end if
    !! Clean local temporary data
    deallocate(radii, csph, stat=istatus)
    if(istatus .ne. 0) then
        call update_error(ddx_error, "Could not deallocate space for " // &
            & "coordinates, radii and charges of atoms")
        return
    end if
end subroutine ddfromfile

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
!!
subroutine ddinit(model, nsph, coords, radii, eps, ddx_data, ddx_error, &
        & force, kappa, eta, shift, lmax, ngrid, incore, maxiter, &
        & jacobi_ndiis, enable_fmm, pm, pl, nproc, logfile, adjoint, eps_int)

    ! mandatory arguments
    integer, intent(in) :: model, nsph
    real(dp), intent(in) :: coords(3, nsph), radii(nsph)
    type(ddx_type), target, intent(out) :: ddx_data
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: eps

    ! optional arguments
    integer, intent(in), optional :: force, adjoint, lmax, ngrid, incore, &
        & maxiter, jacobi_ndiis, enable_fmm, pm, pl, nproc
    real(dp), intent(in), optional :: kappa, eta, shift, eps_int
    character(len=255), intent(in), optional :: logfile

    ! local copies of the optional arguments with dummy default values
    integer :: local_force = 0
    integer :: local_adjoint = 0
    integer :: local_lmax = 6
    integer :: local_ngrid = 302
    integer :: local_incore = 0
    integer :: local_maxiter = 100
    integer :: local_jacobi_ndiis = 20
    integer :: local_enable_fmm = 1
    integer :: local_pm = 8
    integer :: local_pl = 8
    integer :: local_nproc = 1
    real(dp) :: local_kappa = 0.0d0
    real(dp) :: local_eta = 0.1d0
    real(dp) :: local_shift
    real(dp) :: local_eps_int = 1.0d0
    character(len=255) :: local_logfile = ""

    ! arrays for x, y, z coordinates
    real(dp), allocatable :: x(:), y(:), z(:)
    integer :: info

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

    ! for the shift eta the default value depends on the model
    ! ddCOSMO has an interal shift, ddPCM and ddLPB a symmetric shift
    if (present(shift)) then
        local_shift = shift
    else
        if (model.eq.1) then
            local_shift = -one
        else
            local_shift = zero
        end if
    end if

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
        & local_jacobi_ndiis, local_nproc, local_logfile, ddx_data, ddx_error)

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

    call time_push()
    call setup(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, electrostatics, psi, ddx_error)
    call time_pull("setup")
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddrun: setup returned an error, exiting")
        return
    end if

    ! solve the primal linear system
    if (do_guess) then
        call time_push()
        call fill_guess(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, state, tol, ddx_error)
        call time_pull("guess")
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddrun: fill_guess returned an error, exiting")
            return
        end if
    end if
    call time_push()
    call solve(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, tol, ddx_error)
    call time_pull("solve")
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddrun: solve returned an error, exiting")
        return
    end if

    ! compute the energy
    call time_push()
    call energy(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, state, esolv, ddx_error)
    call time_pull("energy")
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddrun: energy returned an error, exiting")
        return
    end if

    ! solve the primal linear system
    if (ddx_data % params % force .eq. 1) then
        if (do_guess) then
            call time_push()
            call fill_guess_adjoint(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, tol, ddx_error)
            call time_pull("guess adjoint")
            if (ddx_error % flag .ne. 0) then
                call update_error(ddx_error, &
                    & "ddrun: fill_guess_adjoint returned an error, exiting")
                return
            end if
        end if
        call time_push()
        call solve_adjoint(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, state, tol, ddx_error)
        call time_pull("solve adjoint")
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddrun: solve_adjoint returned an error, exiting")
            return
        end if
    end if

    ! compute the forces
    if (ddx_data % params % force .eq. 1) then
        force = zero
        call time_push()
        call solvation_force_terms(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, state, electrostatics, force, ddx_error)
        call time_pull("solvation forces")
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
