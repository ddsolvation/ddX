module ddx_cinterface
    use, intrinsic :: iso_c_binding
    use ddx_core
    use ddx_definitions
    use ddx_errors
    use ddx_constants
    use ddx_parameters
    use ddx_workspace
    use ddx_multipolar_solutes
    !
    use ddx_cosmo
    use ddx_pcm
    use ddx_lpb
    implicit none

    type ddx_setup
        type(ddx_params_type)    :: params
        type(ddx_constants_type) :: constants
        type(ddx_workspace_type) :: workspace
    end type ddx_setup

contains
!
! Generic stuff
!
subroutine ddx_get_banner(message, maxlen) bind(C)
    integer(c_int), intent(in), value :: maxlen
    character(len=1, kind=C_char), intent(out) :: message(maxlen)
    integer :: length, i
    character (len=4095) :: header
    call get_banner(header)
    message(maxlen) = c_null_char
    length = min(maxlen-1, 4095)
    do i = length, 1, -1
        if (header(i:i) .eq. ' ') then
            length = i-1
        else
            exit
        endif
    enddo
    message(length + 1) = c_null_char
    do i = 1, length
        message(i) = header(i:i)
    enddo
end

function ddx_supported_lebedev_grids(n, grids) result(c_ngrids) bind(C)
    integer(c_int), intent(in), value ::  n
    integer(c_int), intent(out) :: grids(n)
    integer(c_int) :: c_ngrids
    c_ngrids = min(n, nllg)
    grids(1:c_ngrids) = ng0(1:c_ngrids)
end

! With reference to a atomic sphere `sphere` of radius `r` centred at `a` compute:
!       4π     |x - a|^l
!     ------  ----------- Y_l^m(|x - a|)
!     2l + 1       r^l
! lmax should be identical to the value stored inside c_ddx or less.
subroutine ddx_scaled_ylm(c_ddx, lmax, x, sphere, ylm) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: lmax
    real(c_double), intent(in) :: x(3)
    integer(c_int), intent(in), value :: sphere
    real(c_double), intent(out) :: ylm((lmax+1)**2)
    real(dp) :: delta(3)
    real(dp) :: ratio, rho, ctheta, stheta, cphi, sphi, dnorm
    real(dp) :: vplm((lmax+1)**2), vcos(lmax+1), vsin(lmax+1)
    double precision, external :: dnrm2
    integer :: ind, m, l
    type(ddx_setup), pointer :: ddx
    call c_f_pointer(c_ddx, ddx)

    associate(rsph => ddx%params%rsph, lmax => ddx%params%lmax, vscales => ddx%constants%vscales)
        delta = x(:) - ddx%params%csph(:, sphere)
        dnorm = dnrm2(3, delta, 1)
        delta = delta / dnorm
        call ylmbas(delta, rho, ctheta, stheta, cphi, sphi, lmax, vscales, ylm, vplm, vcos, vsin)
        do l = 0, lmax
            ratio = ddx%constants%v4pi2lp1(l+1) * (dnorm / rsph(sphere))**l
            do m = -l, l
                ind = l*(l+1) + m + 1
                ylm(ind) = ylm(ind) * ratio
            enddo
        enddo
    end associate
end

!
! Error
!

function ddx_allocate_error() result(c_error) bind(C)
    type(c_ptr) :: c_error
    type(ddx_error_type), pointer :: error
    allocate(error)
    c_error = c_loc(error)
end function

function ddx_get_error_flag(c_error) result(has_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(ddx_error_type), pointer :: error
    integer(c_int) :: has_error
    call c_f_pointer(c_error, error)
    has_error = error%flag
end

subroutine ddx_get_error_message(c_error, message, maxlen) bind(C)
    type(c_ptr), intent(in), value :: c_error
    integer(c_int), intent(in), value :: maxlen
    character(len=1, kind=C_char), intent(out) :: message(maxlen)
    character(len=2000) :: error_message
    type(ddx_error_type), pointer :: error
    integer :: length, i
    call c_f_pointer(c_error, error)
    if (error%flag .ne. 0) then
        error_message = error%message
    else
        error_message = ''
    endif

    ! Convert to C message
    message(maxlen) = c_null_char
    length = min(maxlen-1, error%message_length-1)
    do i = length, 1, -1
        if (error_message(i:i) .eq. ' ') then
            length = i-1
        else
            exit
        endif
    enddo
    message(length + 1) = c_null_char
    do i = 1, length
        message(i) = error_message(i:i)
    enddo
end

!
! Setup object
!
function ddx_allocate_model(model, enable_force, solvent_epsilon, solvent_kappa, eta, se, lmax, &
        & n_lebedev, incore, maxiter, jacobi_n_diis, enable_fmm, fmm_multipole_lmax, fmm_local_lmax, &
        & n_proc, n_spheres, sphere_centres, sphere_radii, length_logfile, c_logfile, c_error) result(c_ddx) bind(C)
    integer(c_int), intent(in), value :: model, enable_force, lmax, n_lebedev, maxiter, &
        & incore, jacobi_n_diis, enable_fmm, fmm_multipole_lmax, fmm_local_lmax, n_proc, &
        & n_spheres, length_logfile
    real(c_double), intent(in) :: sphere_centres(3, n_spheres), sphere_radii(n_spheres)
    real(c_double), intent(in), value :: eta, se, solvent_epsilon, solvent_kappa
    type(c_ptr), intent(in), value :: c_error
    type(ddx_error_type), pointer :: error
    !type(c_funptr), value :: printfctn
    type(c_ptr) :: c_ddx
    integer :: passproc, i
    type(ddx_setup), pointer :: ddx
    character(len=1, kind=C_char), intent(in) :: c_logfile(length_logfile)
    character(len=255) :: logfile

    call c_f_pointer(c_error, error)

    ! Allocate DDX object
    allocate(ddx)
    c_ddx = c_loc(ddx)

    ! Convert to Fortran objects
    passproc = n_proc
    logfile = ''
    do i = 1, min(length_logfile, 255)
        if (c_logfile(i) .ne. c_null_char) then
            logfile(i:i) = c_logfile(i)
        endif
    enddo

    call params_init(model, enable_force, solvent_epsilon, solvent_kappa, eta, se, lmax, &
        & n_lebedev, incore, maxiter, jacobi_n_diis, enable_fmm, &
        & fmm_multipole_lmax, fmm_local_lmax, passproc, n_spheres, &
        & sphere_centres, sphere_radii, logfile, ddx%params, error)
    if (error%flag .ne. 0) then
        return
    endif
    call constants_init(ddx%params, ddx%constants, error)
    if (error%flag .ne. 0) then
        return
    endif
    call workspace_init(ddx%params, ddx%constants, ddx%workspace, error)
    if (error%flag .ne. 0) then
        return
    endif
end function

subroutine ddx_deallocate_model(c_ddx, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_error
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call params_free(ddx%params, error)
    call constants_free(ddx%constants, error)
    call workspace_free(ddx%workspace, error)
    deallocate(ddx)
end

subroutine ddx_get_logfile(c_ddx, message, maxlen) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: maxlen
    character(len=1, kind=C_char), intent(out) :: message(maxlen)
    type(ddx_setup), pointer :: ddx
    integer :: length, i
    call c_f_pointer(c_ddx, ddx)
    ! Convert to C message
    message(maxlen) = c_null_char
    length = min(maxlen-1, 255)
    do i = length, 1, -1
        if (ddx%params%output_filename(i:i) .eq. ' ') then
            length = i-1
        else
            exit
        endif
    enddo
    message(length + 1) = c_null_char
    do i = 1, length
        message(i) = ddx%params%output_filename(i:i)
    enddo
end

! Generated block, see scripts/generate_cinterface.py
function ddx_get_enable_fmm(c_ddx) result(c_fmm) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_fmm
    call c_f_pointer(c_ddx, ddx)
    c_fmm = ddx % params % fmm
end function

function ddx_get_enable_force(c_ddx) result(c_force) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_force
    call c_f_pointer(c_ddx, ddx)
    c_force = ddx % params % force
end function

function ddx_get_jacobi_n_diis(c_ddx) result(c_jacobi_ndiis) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_jacobi_ndiis
    call c_f_pointer(c_ddx, ddx)
    c_jacobi_ndiis = ddx % params % jacobi_ndiis
end function

function ddx_get_lmax(c_ddx) result(c_lmax) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_lmax
    call c_f_pointer(c_ddx, ddx)
    c_lmax = ddx % params % lmax
end function

function ddx_get_incore(c_ddx) result(c_matvecmem) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_matvecmem
    call c_f_pointer(c_ddx, ddx)
    c_matvecmem = ddx % params % matvecmem
end function

function ddx_get_maxiter(c_ddx) result(c_maxiter) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_maxiter
    call c_f_pointer(c_ddx, ddx)
    c_maxiter = ddx % params % maxiter
end function

function ddx_get_model(c_ddx) result(c_model) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_model
    call c_f_pointer(c_ddx, ddx)
    c_model = ddx % params % model
end function

function ddx_get_n_lebedev(c_ddx) result(c_ngrid) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_ngrid
    call c_f_pointer(c_ddx, ddx)
    c_ngrid = ddx % params % ngrid
end function

function ddx_get_n_spheres(c_ddx) result(c_nsph) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_nsph
    call c_f_pointer(c_ddx, ddx)
    c_nsph = ddx % params % nsph
end function

function ddx_get_n_proc(c_ddx) result(c_nproc) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_nproc
    call c_f_pointer(c_ddx, ddx)
    c_nproc = ddx % params % nproc
end function

function ddx_get_fmm_local_lmax(c_ddx) result(c_pl) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_pl
    call c_f_pointer(c_ddx, ddx)
    c_pl = ddx % params % pl
end function

function ddx_get_fmm_multipole_lmax(c_ddx) result(c_pm) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_pm
    call c_f_pointer(c_ddx, ddx)
    c_pm = ddx % params % pm
end function

function ddx_get_solvent_epsilon(c_ddx) result(c_eps) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    real(c_double) :: c_eps
    call c_f_pointer(c_ddx, ddx)
    c_eps = ddx % params % eps
end function

function ddx_get_eta(c_ddx) result(c_eta) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    real(c_double) :: c_eta
    call c_f_pointer(c_ddx, ddx)
    c_eta = ddx % params % eta
end function

function ddx_get_solvent_kappa(c_ddx) result(c_kappa) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    real(c_double) :: c_kappa
    call c_f_pointer(c_ddx, ddx)
    c_kappa = ddx % params % kappa
end function

function ddx_get_shift(c_ddx) result(c_se) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    real(c_double) :: c_se
    call c_f_pointer(c_ddx, ddx)
    c_se = ddx % params % se
end function

subroutine ddx_get_sphere_centres(c_ddx, nsph, c_csph) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: nsph
    real(c_double), intent(out) :: c_csph(3, nsph)
    type(ddx_setup), pointer :: ddx
    call c_f_pointer(c_ddx, ddx)
    c_csph(:, :) = ddx % params % csph(:, :)
end subroutine

subroutine ddx_get_sphere_radii(c_ddx, nsph, c_rsph) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: nsph
    real(c_double), intent(out) :: c_rsph(nsph)
    type(ddx_setup), pointer :: ddx
    call c_f_pointer(c_ddx, ddx)
    c_rsph(:) = ddx % params % rsph(:)
end subroutine

function ddx_get_n_basis(c_ddx) result(c_nbasis) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_nbasis
    call c_f_pointer(c_ddx, ddx)
    c_nbasis = ddx % constants % nbasis
end function

function ddx_get_n_cav(c_ddx) result(c_ncav) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int) :: c_ncav
    call c_f_pointer(c_ddx, ddx)
    c_ncav = ddx % constants % ncav
end function

subroutine ddx_get_cavity(c_ddx, ncav, c_ccav) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: ncav
    real(c_double), intent(out) :: c_ccav(3, ncav)
    type(ddx_setup), pointer :: ddx
    call c_f_pointer(c_ddx, ddx)
    c_ccav(:, :) = ddx % constants % ccav(:, :)
end subroutine
! end generated block

!
! State object
!
function ddx_allocate_state(c_ddx, c_error) result(c_state) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_error
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(c_ptr) :: c_state
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    allocate(state)
    call ddx_init_state(ddx%params, ddx%constants, state, error)
    c_state = c_loc(state)
end function

subroutine ddx_deallocate_state(c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_state, c_error
    type(ddx_state_type), pointer :: state
    type(ddx_error_type), pointer :: error
    call c_f_pointer(c_state, state)
    call c_f_pointer(c_error, error)
    call ddx_free_state(state, error)
    deallocate(state)
end subroutine

subroutine ddx_get_x(c_state, nbasis, nsph, x) bind(C)
    type(c_ptr), intent(in), value :: c_state
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nsph, nbasis
    real(c_double), intent(out)  :: x(nbasis, nsph)
    call c_f_pointer(c_state, state)

    if (allocated(state%x_lpb)) then  ! Case for ddLPB
        x(:, :) = state%x_lpb(:, :, 1)  ! Use only the first of the two solutions
    else
        x(:, :) = state%xs(:, :)
    endif
end subroutine

function ddx_get_x_niter(c_state) bind(C) result(c_niter)
    type(c_ptr), intent(in), value :: c_state
    type(ddx_state_type), pointer :: state
    integer(c_int) :: c_niter
    call c_f_pointer(c_state, state)
    if (allocated(state%x_lpb)) then  ! Case for ddLPB
        c_niter = state%x_lpb_niter
    else
        c_niter = state%xs_niter
    endif
end function

subroutine ddx_get_s(c_state, nbasis, nsph, s) bind(C)
    type(c_ptr), intent(in), value :: c_state
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nsph, nbasis
    real(c_double), intent(out)  :: s(nbasis, nsph)
    call c_f_pointer(c_state, state)
    if (allocated(state%x_adj_lpb)) then  ! Case for ddLPB
        s(:, :) = state%x_adj_lpb(:, :, 1)  ! Use only the first of the two solutions
    else
        s(:, :) = state%s(:, :)
    endif
end subroutine

function ddx_get_s_niter(c_state) bind(C) result(c_niter)
    type(c_ptr), intent(in), value :: c_state
    type(ddx_state_type), pointer :: state
    integer(c_int) :: c_niter
    call c_f_pointer(c_state, state)
    if (allocated(state%x_adj_lpb)) then  ! Case for ddLPB
        c_niter = state%x_adj_lpb_niter
    else
        c_niter = state%s_niter
    endif
end function

subroutine ddx_get_xi(c_state, c_ddx, ncav, xi) bind(C)
    type(c_ptr), intent(in), value :: c_state, c_ddx
    type(ddx_state_type), pointer :: state
    type(ddx_setup), pointer :: ddx
    integer(c_int), intent(in), value :: ncav
    real(c_double), intent(out)  :: xi(ncav)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddproject_cav(ddx%params, ddx%constants, state%q, xi)
end subroutine

subroutine ddx_get_zeta_dip(c_state, c_ddx, ncav, zeta_dip) bind(C)
    type(c_ptr), intent(in), value :: c_state, c_ddx
    type(ddx_state_type), pointer :: state
    type(ddx_setup), pointer :: ddx
    integer(c_int), intent(in), value :: ncav
    real(c_double), intent(out)  :: zeta_dip(3, ncav)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    if (allocated(state%zeta_dip)) then
        zeta_dip(:, :) = state%zeta_dip(:, :)
    else
        zeta_dip(:, :) = 0d0
    endif
end subroutine

!
! Cosmo
!
! Setup the problem in the state
subroutine ddx_cosmo_setup(c_ddx, c_state, ncav, nbasis, nsph, psi, phi_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    integer(c_int), intent(in), value :: ncav, nbasis, nsph
    real(c_double), intent(in) :: phi_cav(ncav), psi(nbasis, nsph)
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddcosmo_setup(ddx%params, ddx%constants, ddx%workspace, state, phi_cav, psi)
end subroutine

! Put a guess for the problem into the state (optional)
subroutine ddx_cosmo_guess(c_ddx, c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddcosmo_guess(ddx%params, ddx%constants, ddx%workspace, state, error)
end

! Solve the problem
subroutine ddx_cosmo_solve(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddcosmo_solve(ddx%params, ddx%constants, ddx%workspace, state, tol, error)
end subroutine

! Put a guess for the adjoint problem into the state (optional)
subroutine ddx_cosmo_guess_adjoint(c_ddx, c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddcosmo_guess_adjoint(ddx%params, ddx%constants, ddx%workspace, state, error)
end

! Solve the adjoint problem inside the state
subroutine ddx_cosmo_solve_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddcosmo_solve_adjoint(ddx%params, ddx%constants, ddx%workspace, state, tol, error)
end

! Compute the ddCOSMO energy
function ddx_cosmo_energy(c_ddx, c_state, c_error) result(c_energy) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    real(c_double) :: c_energy
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddcosmo_energy(ddx%constants, state, c_energy)
end function

! Compute the forces
subroutine ddx_cosmo_solvation_force_terms(c_ddx, c_state, nsph, forces, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nsph
    real(c_double), intent(out) :: forces(3, nsph)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddcosmo_solvation_force_terms(ddx%params, ddx%constants, ddx%workspace, state, forces)
end

!
! PCM
!
subroutine ddx_pcm_setup(c_ddx, c_state, ncav, nbasis, nsph, psi, phi_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    integer(c_int), intent(in), value :: ncav, nbasis, nsph
    real(c_double), intent(in) :: phi_cav(ncav), psi(nbasis, nsph)
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddpcm_setup(ddx%params, ddx%constants, ddx%workspace, state, phi_cav, psi)
end subroutine

subroutine ddx_pcm_guess(c_ddx, c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddpcm_guess(ddx%params, ddx%constants, ddx%workspace, state, error)
end

subroutine ddx_pcm_solve(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddpcm_solve(ddx%params, ddx%constants, ddx%workspace, state, tol, error)
end subroutine

subroutine ddx_pcm_guess_adjoint(c_ddx, c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddpcm_guess_adjoint(ddx%params, ddx%constants, ddx%workspace, state, error)
end

subroutine ddx_pcm_solve_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddpcm_solve_adjoint(ddx%params, ddx%constants, ddx%workspace, state, tol, error)
end

function ddx_pcm_energy(c_ddx, c_state, c_error) result(c_energy) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    real(c_double) :: c_energy
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddpcm_energy(ddx%constants, state, c_energy)
end function

subroutine ddx_pcm_solvation_force_terms(c_ddx, c_state, nsph, forces, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nsph
    real(c_double), intent(out) :: forces(3, nsph)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddpcm_solvation_force_terms(ddx%params, ddx%constants, ddx%workspace, state, forces)
end

!
! LPB
!
subroutine ddx_lpb_setup(c_ddx, c_state, ncav, nbasis, nsph, psi, phi_cav, gradphi_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    integer(c_int), intent(in), value :: ncav, nbasis, nsph
    real(c_double), intent(in) :: phi_cav(ncav), psi(nbasis, nsph), gradphi_cav(3, ncav)
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddlpb_setup(ddx%params, ddx%constants, ddx%workspace, state, phi_cav, gradphi_cav, psi, error)
end subroutine

subroutine ddx_lpb_guess(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    real(c_double), intent(in), value :: tol
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddlpb_guess(ddx%params, ddx%constants, ddx%workspace, state, tol, error)
end

subroutine ddx_lpb_solve(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddlpb_solve(ddx%params, ddx%constants, ddx%workspace, state, tol, error)
end subroutine

subroutine ddx_lpb_guess_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    real(c_double), intent(in), value :: tol
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddlpb_guess_adjoint(ddx%params, ddx%constants, ddx%workspace, state, tol, error)
end

subroutine ddx_lpb_solve_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddlpb_solve_adjoint(ddx%params, ddx%constants, ddx%workspace, state, tol, error)
end

function ddx_lpb_energy(c_ddx, c_state, c_error) result(c_energy) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    real(c_double) :: c_energy
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    call ddlpb_energy(ddx%constants, state, c_energy)
end function

! TODO LPB force terms not yet supported in C and python interface


!
! multipolar solutes
!
subroutine ddx_multipole_electrostatics_0(c_ddx, nsph, ncav, nmultipoles, multipoles, &
        & phi_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer(c_int), intent(in), value :: nsph, ncav, nmultipoles
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    real(c_double), intent(out) :: phi_cav(ncav)
    type(c_ptr), intent(in), value :: c_error
    type(ddx_error_type), pointer :: error
    integer :: mmax
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_ddx, ddx)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    call build_phi(ddx%params, ddx%constants, ddx%workspace, multipoles, mmax, &
            & phi_cav, error)
end

subroutine ddx_multipole_electrostatics_1(c_ddx, nsph, ncav, nmultipoles, multipoles, &
        & phi_cav, e_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    integer(c_int), intent(in), value :: nsph, ncav, nmultipoles
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    real(c_double), intent(out) :: phi_cav(ncav), e_cav(3, ncav)
    integer :: mmax
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    call build_e(ddx%params, ddx%constants, ddx%workspace, multipoles, mmax, &
            & phi_cav, e_cav, error)
end

subroutine ddx_multipole_electrostatics_2(c_ddx, nsph, ncav, nmultipoles, multipoles, &
        & phi_cav, e_cav, g_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    integer(c_int), intent(in), value :: nsph, ncav, nmultipoles
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    real(c_double), intent(out) :: phi_cav(ncav), e_cav(3, ncav), g_cav(3, 3, ncav)
    integer :: mmax
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    call build_g(ddx%params, ddx%constants, ddx%workspace, multipoles, mmax, &
            & phi_cav, e_cav, g_cav, error)
end

subroutine ddx_multipole_psi(c_ddx, nbasis, nsph, nmultipoles, multipoles, psi, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    integer(c_int), intent(in), value :: nmultipoles, nsph, nbasis
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    real(c_double), intent(out) :: psi(nbasis, nsph)
    integer :: mmax
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    call build_psi(ddx%params, multipoles, mmax, psi)
end

subroutine ddx_multipole_forces(c_ddx, c_state, nsph, ncav, nmultipoles, multipoles, &
                & e_cav, forces, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_error_type), pointer :: error
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nmultipoles, nsph, ncav
    integer :: mmax
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph), &
        & e_cav(3, ncav)
    real(c_double), intent(out) :: forces(3, nsph)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_error, error)
    call c_f_pointer(c_state, state)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    call grad_phi(ddx%params, ddx%constants, ddx%workspace, state, mmax, &
        & multipoles, forces, e_cav, error)
end

end
