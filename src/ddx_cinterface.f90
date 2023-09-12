module ddx_cinterface
    use, intrinsic :: iso_c_binding
    use ddx
    implicit none

    type ddx_setup_type
        type(ddx_params_type)    :: params
        type(ddx_constants_type) :: constants
        type(ddx_workspace_type) :: workspace
    end type ddx_setup_type

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
    type(ddx_setup_type), pointer :: ddx_model
    call c_f_pointer(c_ddx, ddx_model)

    associate(rsph => ddx_model%params%rsph, lmax => ddx_model%params%lmax, vscales => ddx_model%constants%vscales)
        delta = x(:) - ddx_model%params%csph(:, sphere)
        dnorm = dnrm2(3, delta, 1)
        delta = delta / dnorm
        call ylmbas(delta, rho, ctheta, stheta, cphi, sphi, lmax, vscales, ylm, vplm, vcos, vsin)
        do l = 0, lmax
            ratio = ddx_model%constants%v4pi2lp1(l+1) * (dnorm / rsph(sphere))**l
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
    type(ddx_error_type), pointer :: ddx_error
    allocate(ddx_error)
    c_error = c_loc(ddx_error)
end function

function ddx_get_error_flag(c_error) result(has_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(ddx_error_type), pointer :: ddx_error
    integer(c_int) :: has_error
    call c_f_pointer(c_error, ddx_error)
    has_error = ddx_error%flag
end

subroutine ddx_get_error_message(c_error, message, maxlen) bind(C)
    type(c_ptr), intent(in), value :: c_error
    integer(c_int), intent(in), value :: maxlen
    character(len=1, kind=C_char), intent(out) :: message(maxlen)
    character(len=2047) :: error_message
    type(ddx_error_type), pointer :: ddx_error
    integer :: length, i
    call c_f_pointer(c_error, ddx_error)
    if (ddx_error%flag .ne. 0) then
        error_message = ddx_error%message
    else
        error_message = ''
    endif

    ! Convert to C message
    message(maxlen) = c_null_char
    length = min(maxlen-1, ddx_error%message_length-1)
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
! Electrostatics container
!
function ddx_allocate_electrostatics(c_ddx, c_error) result(c_electrostatics) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_error
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(c_ptr) :: c_electrostatics
    type(ddx_electrostatics_type), pointer :: electrostatics
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    allocate(electrostatics)
    call allocate_electrostatics(ddx_model%params, ddx_model%constants, electrostatics, ddx_error)
    c_electrostatics = c_loc(electrostatics)
end function

subroutine ddx_deallocate_electrostatics(c_electrostatics, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_electrostatics, c_error
    type(ddx_electrostatics_type), pointer :: electrostatics
    type(ddx_error_type), pointer :: ddx_error
    call c_f_pointer(c_electrostatics, electrostatics)
    call c_f_pointer(c_error, ddx_error)
    call deallocate_electrostatics(electrostatics, ddx_error)
    deallocate(electrostatics)
end

subroutine ddx_multipole_electrostatics(c_ddx, nsph, nmultipoles, multipoles, &
        & c_electrostatics, c_error) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: c_ddx, c_electrostatics, c_error
    integer(c_int), intent(in), value :: nsph, nmultipoles
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_electrostatics_type), pointer :: electrostatics
    integer :: mmax

    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_electrostatics, electrostatics)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)

    call multipole_electrostatics(ddx_model%params, ddx_model%constants, ddx_model%workspace, &
        & multipoles, mmax, electrostatics, ddx_error)
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
    type(ddx_error_type), pointer :: ddx_error
    !type(c_funptr), value :: printfctn
    type(c_ptr) :: c_ddx
    integer :: passproc, i
    type(ddx_setup_type), pointer :: ddx_model
    character(len=1, kind=C_char), intent(in) :: c_logfile(length_logfile)
    character(len=255) :: logfile

    call c_f_pointer(c_error, ddx_error)

    ! Allocate DDX object
    allocate(ddx_model)
    c_ddx = c_loc(ddx_model)

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
        & sphere_centres, sphere_radii, logfile, ddx_model%params, ddx_error)
    if (ddx_error%flag .ne. 0) then
        return
    endif
    call constants_init(ddx_model%params, ddx_model%constants, ddx_error)
    if (ddx_error%flag .ne. 0) then
        return
    endif
    call workspace_init(ddx_model%params, ddx_model%constants, ddx_model%workspace, ddx_error)
    if (ddx_error%flag .ne. 0) then
        return
    endif
end function

subroutine ddx_deallocate_model(c_ddx, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_error
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call params_free(ddx_model%params, ddx_error)
    call constants_free(ddx_model%constants, ddx_error)
    call workspace_free(ddx_model%workspace, ddx_error)
    deallocate(ddx_model)
end

subroutine ddx_get_logfile(c_ddx, message, maxlen) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: maxlen
    character(len=1, kind=C_char), intent(out) :: message(maxlen)
    type(ddx_setup_type), pointer :: ddx_model
    integer :: length, i
    call c_f_pointer(c_ddx, ddx_model)
    ! Convert to C message
    message(maxlen) = c_null_char
    length = min(maxlen-1, 255)
    do i = length, 1, -1
        if (ddx_model%params%output_filename(i:i) .eq. ' ') then
            length = i-1
        else
            exit
        endif
    enddo
    message(length + 1) = c_null_char
    do i = 1, length
        message(i) = ddx_model%params%output_filename(i:i)
    enddo
end

! Generated block, see scripts/generate_cinterface.py
function ddx_get_enable_fmm(c_ddx) result(c_fmm) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_fmm
    call c_f_pointer(c_ddx, ddx_model)
    c_fmm = ddx_model % params % fmm
end function

function ddx_get_enable_force(c_ddx) result(c_force) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_force
    call c_f_pointer(c_ddx, ddx_model)
    c_force = ddx_model % params % force
end function

function ddx_get_jacobi_n_diis(c_ddx) result(c_jacobi_ndiis) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_jacobi_ndiis
    call c_f_pointer(c_ddx, ddx_model)
    c_jacobi_ndiis = ddx_model % params % jacobi_ndiis
end function

function ddx_get_lmax(c_ddx) result(c_lmax) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_lmax
    call c_f_pointer(c_ddx, ddx_model)
    c_lmax = ddx_model % params % lmax
end function

function ddx_get_incore(c_ddx) result(c_matvecmem) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_matvecmem
    call c_f_pointer(c_ddx, ddx_model)
    c_matvecmem = ddx_model % params % matvecmem
end function

function ddx_get_maxiter(c_ddx) result(c_maxiter) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_maxiter
    call c_f_pointer(c_ddx, ddx_model)
    c_maxiter = ddx_model % params % maxiter
end function

function ddx_get_model(c_ddx) result(c_model) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_model
    call c_f_pointer(c_ddx, ddx_model)
    c_model = ddx_model % params % model
end function

function ddx_get_n_lebedev(c_ddx) result(c_ngrid) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_ngrid
    call c_f_pointer(c_ddx, ddx_model)
    c_ngrid = ddx_model % params % ngrid
end function

function ddx_get_n_spheres(c_ddx) result(c_nsph) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_nsph
    call c_f_pointer(c_ddx, ddx_model)
    c_nsph = ddx_model % params % nsph
end function

function ddx_get_n_proc(c_ddx) result(c_nproc) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_nproc
    call c_f_pointer(c_ddx, ddx_model)
    c_nproc = ddx_model % params % nproc
end function

function ddx_get_fmm_local_lmax(c_ddx) result(c_pl) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_pl
    call c_f_pointer(c_ddx, ddx_model)
    c_pl = ddx_model % params % pl
end function

function ddx_get_fmm_multipole_lmax(c_ddx) result(c_pm) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_pm
    call c_f_pointer(c_ddx, ddx_model)
    c_pm = ddx_model % params % pm
end function

function ddx_get_solvent_epsilon(c_ddx) result(c_eps) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    real(c_double) :: c_eps
    call c_f_pointer(c_ddx, ddx_model)
    c_eps = ddx_model % params % eps
end function

function ddx_get_eta(c_ddx) result(c_eta) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    real(c_double) :: c_eta
    call c_f_pointer(c_ddx, ddx_model)
    c_eta = ddx_model % params % eta
end function

function ddx_get_solvent_kappa(c_ddx) result(c_kappa) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    real(c_double) :: c_kappa
    call c_f_pointer(c_ddx, ddx_model)
    c_kappa = ddx_model % params % kappa
end function

function ddx_get_shift(c_ddx) result(c_se) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    real(c_double) :: c_se
    call c_f_pointer(c_ddx, ddx_model)
    c_se = ddx_model % params % se
end function

subroutine ddx_get_sphere_centres(c_ddx, nsph, c_csph) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: nsph
    real(c_double), intent(out) :: c_csph(3, nsph)
    type(ddx_setup_type), pointer :: ddx_model
    call c_f_pointer(c_ddx, ddx_model)
    c_csph(:, :) = ddx_model % params % csph(:, :)
end subroutine

subroutine ddx_get_sphere_radii(c_ddx, nsph, c_rsph) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: nsph
    real(c_double), intent(out) :: c_rsph(nsph)
    type(ddx_setup_type), pointer :: ddx_model
    call c_f_pointer(c_ddx, ddx_model)
    c_rsph(:) = ddx_model % params % rsph(:)
end subroutine

function ddx_get_n_basis(c_ddx) result(c_nbasis) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_nbasis
    call c_f_pointer(c_ddx, ddx_model)
    c_nbasis = ddx_model % constants % nbasis
end function

function ddx_get_n_cav(c_ddx) result(c_ncav) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int) :: c_ncav
    call c_f_pointer(c_ddx, ddx_model)
    c_ncav = ddx_model % constants % ncav
end function

subroutine ddx_get_cavity(c_ddx, ncav, c_ccav) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: ncav
    real(c_double), intent(out) :: c_ccav(3, ncav)
    type(ddx_setup_type), pointer :: ddx_model
    call c_f_pointer(c_ddx, ddx_model)
    c_ccav(:, :) = ddx_model % constants % ccav(:, :)
end subroutine
! end generated block

!
! State object
!
function ddx_allocate_state(c_ddx, c_error) result(c_state) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_error
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(c_ptr) :: c_state
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    allocate(state)
    call allocate_state(ddx_model%params, ddx_model%constants, state, ddx_error)
    c_state = c_loc(state)
end function

subroutine ddx_deallocate_state(c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_state, c_error
    type(ddx_state_type), pointer :: state
    type(ddx_error_type), pointer :: ddx_error
    call c_f_pointer(c_state, state)
    call c_f_pointer(c_error, ddx_error)
    call deallocate_state(state, ddx_error)
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
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int), intent(in), value :: ncav
    real(c_double), intent(out)  :: xi(ncav)
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_state, state)
    call ddproject_cav(ddx_model%params, ddx_model%constants, state%q, xi)
end subroutine

subroutine ddx_get_zeta_dip(c_state, c_ddx, ncav, zeta_dip) bind(C)
    type(c_ptr), intent(in), value :: c_state, c_ddx
    type(ddx_state_type), pointer :: state
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int), intent(in), value :: ncav
    real(c_double), intent(out)  :: zeta_dip(3, ncav)
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_state, state)
    if (allocated(state%zeta_dip)) then
        zeta_dip(:, :) = state%zeta_dip(:, :)
    else
        zeta_dip(:, :) = 0d0
    endif
end subroutine

!
! High level APIs (model nonspecific)
!
function ddx_ddrun(c_ddx, c_state, c_electrostatics, nbasis, nsph, &
        & psi, tol, force, read_guess, c_error) result(c_energy) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state, c_electrostatics, &
        & c_error
    integer(c_int), intent(in), value :: nbasis, nsph
    real(c_double), intent(in) :: psi(nbasis, nsph)
    real(c_double), intent(in), value :: tol
    real(c_double), intent(out) :: force(3, nsph)
    real(c_double) :: c_energy
    integer(c_int), intent(in), value :: read_guess
    type(ddx_setup_type), pointer :: ddx_model
    logical :: do_guess

    ! at least return something in case of errors
    c_energy = zero

    ! setup
    do_guess = read_guess.ne.0
    call c_f_pointer(c_ddx, ddx_model)
    call ddx_setup(c_ddx, c_state, c_electrostatics, nbasis, nsph, psi, c_error)
    if (ddx_get_error_flag(c_error) .ne. 0) return

    ! primal linear system
    if (do_guess) then
        call ddx_fill_guess(c_ddx, c_state, tol, c_error)
        if (ddx_get_error_flag(c_error) .ne. 0) return
    end if
    call ddx_solve(c_ddx, c_state, tol, c_error)
    if (ddx_get_error_flag(c_error) .ne. 0) return

    ! energy
    c_energy = ddx_energy(c_ddx, c_state, c_error)
    if (ddx_get_error_flag(c_error) .ne. 0) return

    ! adjoint linear system
    if (ddx_model%params%force .eq. 1) then
        if (do_guess) then
            call ddx_fill_guess(c_ddx, c_state, tol, c_error)
        end if
        if (ddx_get_error_flag(c_error) .ne. 0) return
        call ddx_solve_adjoint(c_ddx, c_state, tol, c_error)
        if (ddx_get_error_flag(c_error) .ne. 0) return
    end if

    ! forces
    if (ddx_model%params%force .eq. 1) then
        call ddx_solvation_force_terms(c_ddx, c_state, c_electrostatics, &
            & nsph, force, c_error)
        if (ddx_get_error_flag(c_error) .ne. 0) return
    end if

end function

subroutine ddx_setup(c_ddx, c_state, c_electrostatics, nbasis, nsph, psi, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state, c_electrostatics, &
        & c_error
    integer(c_int), intent(in), value :: nbasis, nsph
    real(c_double), intent(in) :: psi(nbasis, nsph)
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_state_type), pointer :: state
    type(ddx_electrostatics_type), pointer :: electrostatics
    type(ddx_error_type), pointer :: ddx_error
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_state, state)
    call c_f_pointer(c_electrostatics, electrostatics)
    call c_f_pointer(c_error, ddx_error)
    call setup(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, &
        & electrostatics, psi, ddx_error)
end subroutine

subroutine ddx_fill_guess(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state, c_error
    real(c_double), intent(in), value :: tol
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_state_type), pointer :: state
    type(ddx_error_type), pointer :: ddx_error
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_state, state)
    call c_f_pointer(c_error, ddx_error)
    call fill_guess(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, &
        & tol, ddx_error)
end subroutine

subroutine ddx_fill_guess_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state, c_error
    real(c_double), intent(in), value :: tol
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_state_type), pointer :: state
    type(ddx_error_type), pointer :: ddx_error
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_state, state)
    call c_f_pointer(c_error, ddx_error)
    call fill_guess_adjoint(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, &
        & tol, ddx_error)
end subroutine

subroutine ddx_solve(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state, c_error
    real(c_double), intent(in), value :: tol
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_state_type), pointer :: state
    type(ddx_error_type), pointer :: ddx_error
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_state, state)
    call c_f_pointer(c_error, ddx_error)
    call solve(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, tol, ddx_error)
end subroutine

subroutine ddx_solve_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state, c_error
    real(c_double), intent(in), value :: tol
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_state_type), pointer :: state
    type(ddx_error_type), pointer :: ddx_error
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_state, state)
    call c_f_pointer(c_error, ddx_error)
    call solve_adjoint(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, &
        & tol, ddx_error)
end subroutine

function ddx_energy(c_ddx, c_state, c_error) result(c_energy) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state, c_error
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_state_type), pointer :: state
    type(ddx_error_type), pointer :: ddx_error
    real(c_double) :: c_energy
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_state, state)
    call c_f_pointer(c_error, ddx_error)
    call energy(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, &
        & c_energy, ddx_error)
end function

subroutine ddx_solvation_force_terms(c_ddx, c_state, c_electrostatics, &
        & nsph, forces, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state, c_electrostatics, &
        & c_error
    integer(c_int), intent(in), value :: nsph
    real(c_double), intent(out) :: forces(3, nsph)
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_state_type), pointer :: state
    type(ddx_electrostatics_type), pointer :: electrostatics
    type(ddx_error_type), pointer :: ddx_error
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_state, state)
    call c_f_pointer(c_electrostatics, electrostatics)
    call c_f_pointer(c_error, ddx_error)
    call solvation_force_terms(ddx_model%params, ddx_model%constants, ddx_model%workspace, &
        & state, electrostatics, forces, ddx_error)
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
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call cosmo_setup(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, phi_cav, psi, ddx_error)
end subroutine

! Put a guess for the problem into the state (optional)
subroutine ddx_cosmo_guess(c_ddx, c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call cosmo_guess(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, ddx_error)
end

! Solve the problem
subroutine ddx_cosmo_solve(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call cosmo_solve(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, tol, ddx_error)
end subroutine

! Put a guess for the adjoint problem into the state (optional)
subroutine ddx_cosmo_guess_adjoint(c_ddx, c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call cosmo_guess_adjoint(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, ddx_error)
end

! Solve the adjoint problem inside the state
subroutine ddx_cosmo_solve_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call cosmo_solve_adjoint(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, tol, ddx_error)
end

! Compute the ddCOSMO energy
function ddx_cosmo_energy(c_ddx, c_state, c_error) result(c_energy) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    real(c_double) :: c_energy
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call cosmo_energy(ddx_model%constants, state, c_energy, ddx_error)
end function

! Compute the forces
subroutine ddx_cosmo_solvation_force_terms(c_ddx, c_state, nsph, ncav, e_cav, forces, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nsph, ncav
    real(c_double), intent(in) :: e_cav(3, ncav)
    real(c_double), intent(out) :: forces(3, nsph)
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call cosmo_solvation_force_terms(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, e_cav, forces, ddx_error)
end

!
! PCM
!
subroutine ddx_pcm_setup(c_ddx, c_state, ncav, nbasis, nsph, psi, phi_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    integer(c_int), intent(in), value :: ncav, nbasis, nsph
    real(c_double), intent(in) :: phi_cav(ncav), psi(nbasis, nsph)
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call pcm_setup(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, phi_cav, psi, ddx_error)
end subroutine

subroutine ddx_pcm_guess(c_ddx, c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call pcm_guess(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, ddx_error)
end

subroutine ddx_pcm_solve(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call pcm_solve(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, tol, ddx_error)
end subroutine

subroutine ddx_pcm_guess_adjoint(c_ddx, c_state, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call pcm_guess_adjoint(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, ddx_error)
end

subroutine ddx_pcm_solve_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call pcm_solve_adjoint(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, tol, ddx_error)
end

function ddx_pcm_energy(c_ddx, c_state, c_error) result(c_energy) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    real(c_double) :: c_energy
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call pcm_energy(ddx_model%constants, state, c_energy, ddx_error)
end function

subroutine ddx_pcm_solvation_force_terms(c_ddx, c_state, nsph, ncav, e_cav, forces, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nsph, ncav
    real(c_double), intent(in) :: e_cav(3, ncav)
    real(c_double), intent(out) :: forces(3, nsph)
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call pcm_solvation_force_terms(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, e_cav, forces, ddx_error)
end

!
! LPB
!
subroutine ddx_lpb_setup(c_ddx, c_state, ncav, nbasis, nsph, psi, phi_cav, gradphi_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    integer(c_int), intent(in), value :: ncav, nbasis, nsph
    real(c_double), intent(in) :: phi_cav(ncav), psi(nbasis, nsph), gradphi_cav(3, ncav)
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call lpb_setup(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, phi_cav, gradphi_cav, psi, ddx_error)
end subroutine

subroutine ddx_lpb_guess(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    real(c_double), intent(in), value :: tol
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call lpb_guess(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, tol, ddx_error)
end

subroutine ddx_lpb_solve(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call lpb_solve(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, tol, ddx_error)
end subroutine

subroutine ddx_lpb_guess_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    real(c_double), intent(in), value :: tol
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call lpb_guess_adjoint(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, tol, ddx_error)
end

subroutine ddx_lpb_solve_adjoint(c_ddx, c_state, tol, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    real(c_double), intent(in), value :: tol
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call lpb_solve_adjoint(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, tol, ddx_error)
end

function ddx_lpb_energy(c_ddx, c_state, c_error) result(c_energy) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    real(c_double) :: c_energy
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call lpb_energy(ddx_model%constants, state, c_energy, ddx_error)
end function

subroutine ddx_lpb_solvation_force_terms(c_ddx, c_state, nsph, ncav, g_cav, forces, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nsph, ncav
    real(c_double), intent(in) :: g_cav(3, 3, ncav)
    real(c_double), intent(out) :: forces(3, nsph)
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    call lpb_solvation_force_terms(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, g_cav, forces, ddx_error)
end

!
! multipolar solutes
!
subroutine ddx_multipole_electrostatics_0(c_ddx, nsph, ncav, nmultipoles, multipoles, &
        & phi_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    integer(c_int), intent(in), value :: nsph, ncav, nmultipoles
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    real(c_double), intent(out) :: phi_cav(ncav)
    type(c_ptr), intent(in), value :: c_error
    type(ddx_error_type), pointer :: ddx_error
    integer :: mmax
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_ddx, ddx_model)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    call multipole_electrostatics_0(ddx_model%params, ddx_model%constants, &
        & ddx_model%workspace, multipoles, mmax, phi_cav, ddx_error)
end

subroutine ddx_multipole_electrostatics_1(c_ddx, nsph, ncav, nmultipoles, multipoles, &
        & phi_cav, e_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    integer(c_int), intent(in), value :: nsph, ncav, nmultipoles
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    real(c_double), intent(out) :: phi_cav(ncav), e_cav(3, ncav)
    integer :: mmax
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    call multipole_electrostatics_1(ddx_model%params, ddx_model%constants, &
        & ddx_model%workspace, multipoles, mmax, phi_cav, e_cav, ddx_error)
end

subroutine ddx_multipole_electrostatics_2(c_ddx, nsph, ncav, nmultipoles, multipoles, &
        & phi_cav, e_cav, g_cav, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    integer(c_int), intent(in), value :: nsph, ncav, nmultipoles
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    real(c_double), intent(out) :: phi_cav(ncav), e_cav(3, ncav), g_cav(3, 3, ncav)
    integer :: mmax
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    call multipole_electrostatics_2(ddx_model%params, ddx_model%constants, &
        & ddx_model%workspace, multipoles, mmax, phi_cav, e_cav, g_cav, ddx_error)
end

subroutine ddx_multipole_psi(c_ddx, nbasis, nsph, nmultipoles, multipoles, psi, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    integer(c_int), intent(in), value :: nmultipoles, nsph, nbasis
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    real(c_double), intent(out) :: psi(nbasis, nsph)
    integer :: mmax
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    call multipole_psi(ddx_model%params, multipoles, mmax, psi)
end

subroutine ddx_multipole_force_terms(c_ddx, c_state, nsph, nmultipoles, multipoles, &
        & forces, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_error
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup_type), pointer :: ddx_model
    type(ddx_error_type), pointer :: ddx_error
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nmultipoles, nsph
    integer :: mmax
    real(c_double), intent(in) :: multipoles(nmultipoles, nsph)
    real(c_double), intent(inout) :: forces(3, nsph)
    call c_f_pointer(c_ddx, ddx_model)
    call c_f_pointer(c_error, ddx_error)
    call c_f_pointer(c_state, state)
    mmax = int(sqrt(dble(nmultipoles)) - 1d0)
    forces = zero
    call multipole_force_terms(ddx_model%params, ddx_model%constants, ddx_model%workspace, state, mmax, &
        & multipoles, forces, ddx_error)
end

end
