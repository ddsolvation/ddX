module ddx_cinterface
    use, intrinsic :: iso_c_binding
    use ddx_core
    use ddx_definitions
    use ddx_constants
    use ddx_parameters
    use ddx_workspace
    !
    use ddx_cosmo
    use ddx_pcm
    implicit none

    type ddx_setup
        type(ddx_params_type)    :: params
        type(ddx_constants_type) :: constants
        type(ddx_workspace_type) :: workspace
        character(len=255)       :: error_message  ! most recent error message (if any)
    end type ddx_setup

contains
!
! Generic stuff
!
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

subroutine ddx_nuclear_contributions(c_ddx, nsph, ncav, nbasis, phi, gradphi, psi) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: nsph, ncav, nbasis
    type(ddx_setup), pointer :: ddx
    real(c_double), intent(out) :: phi(ncav)
    real(c_double), intent(out) :: gradphi(3, ncav)
    real(c_double), intent(out) :: psi(nbasis, nsph)
    real(dp) :: hessianphi(3, 3, ncav)  ! not used, hence not allocated.
    call c_f_pointer(c_ddx, ddx)
    call mkrhs(ddx%params, ddx%constants, ddx%workspace, 1, phi, 1, gradphi, 0, &
            & hessianphi, psi)
end

!
! Setup object
!
function ddx_allocate_model(model, enable_force, solvent_epsilon, solvent_kappa, eta, lmax, &
        & n_lebedev, incore, maxiter, jacobi_n_diis, enable_fmm, fmm_multipole_lmax, fmm_local_lmax, &
        & n_proc, n_spheres, sphere_charges, sphere_centres, &
        !& sphere_radii, printfctn, info) result(c_ddx) bind(C)
        & sphere_radii, info) result(c_ddx) bind(C)
    integer(c_int), intent(in), value :: model, enable_force, lmax, n_lebedev, maxiter, &
        & incore, jacobi_n_diis, enable_fmm, fmm_multipole_lmax, fmm_local_lmax, n_proc, n_spheres
    real(c_double), intent(in) :: sphere_charges(n_spheres), sphere_centres(3, n_spheres), &
        & sphere_radii(n_spheres)
    real(c_double), intent(in), value :: eta, solvent_epsilon, solvent_kappa
    integer(c_int), intent(out) :: info
    !type(c_funptr), value :: printfctn
    type(c_ptr) :: c_ddx
    integer :: passproc, gmresr_j, gmresr_dim, itersolver
    real(dp) :: se
    type(ddx_setup), pointer :: ddx

    ! interface
    !     subroutine print_characters(string) bind(C)
    !         character(kind=c_char), dimension(*), intent(in) :: string
    !     end subroutine
    ! end interface
    ! procedure(print_characters), pointer :: cprint
    ! call c_f_procpointer(printfctn, cprint)  ! Convert C to fortran procedure pointer

    ! Allocate DDX object
    allocate(ddx)
    c_ddx = c_loc(ddx)

    passproc = n_proc
    se = 0.0        ! Hard-code centred regularisation
    itersolver = 1  ! Hard-code Jacobi DIIS
    gmresr_j   = 1  ! Dummy as GMRES not used
    gmresr_dim = 5  ! Dummy as GMRES not used

    ddx%error_message = "No error"
    call params_init(model, enable_force, solvent_epsilon, solvent_kappa, eta, se, lmax, &
        & n_lebedev, incore, itersolver, maxiter, jacobi_n_diis, gmresr_j, gmresr_dim, enable_fmm, &
        & fmm_multipole_lmax, fmm_local_lmax, passproc, n_spheres, sphere_charges, &
        & sphere_centres, sphere_radii, print_func_default, ddx%params, info)
    if (info .ne. 0) then
        ddx%error_message = ddx%params%error_message
        return
    endif
    info = 0
    call constants_init(ddx%params, ddx%constants, info)
    if (info .ne. 0) then
        ddx%error_message = ddx%constants%error_message
        return
    endif
    info = 0
    call workspace_init(ddx%params, ddx%constants, ddx%workspace, info)
    if (info .ne. 0) then
        ddx%error_message = ddx%workspace%error_message
        return
    endif
end function

subroutine ddx_deallocate_model(c_ddx) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    integer :: info
    call c_f_pointer(c_ddx, ddx)
    call params_deinit(ddx%params, info)
    call constants_free(ddx%constants, info)
    call workspace_free(ddx%workspace, info)
    deallocate(ddx)
end

subroutine ddx_get_error_message(c_ddx, message, maxlen) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: maxlen
    character(len=1, kind=C_char), intent(out) :: message(maxlen)
    type(ddx_setup), pointer :: ddx
    integer :: length, i
    character :: ch
    call c_f_pointer(c_ddx, ddx)
    message(maxlen) = c_null_char
    length = min(maxlen-1, 255)
    do i = length, 1, -1
        if (ddx%error_message(i:i) .eq. ' ') then
            length = i-1
        else
            exit
        endif
    enddo
    message(length + 1) = c_null_char
    do i = 1, length
        message(i) = ddx%error_message(i:i)
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

subroutine ddx_get_sphere_charges(c_ddx, nsph, c_charge) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: nsph
    real(c_double), intent(out) :: c_charge(nsph)
    type(ddx_setup), pointer :: ddx
    call c_f_pointer(c_ddx, ddx)
    c_charge(:) = ddx % params % charge(:)
end subroutine

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
function ddx_allocate_state(c_ddx) result(c_state) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    type(c_ptr) :: c_state
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    allocate(state)
    c_state = c_loc(state)
    call ddx_init_state(ddx%params, ddx%constants, state)
end function

subroutine ddx_deallocate_state(c_state) bind(C)
    type(c_ptr), intent(in), value :: c_state
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_state, state)
    call ddx_free_state(state)
    deallocate(state)
end subroutine

subroutine ddx_get_x(c_state, nbasis, nsph, x) bind(C)
    type(c_ptr), intent(in), value :: c_state
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nsph, nbasis
    real(c_double), intent(out)  :: x(nbasis, nsph)
    call c_f_pointer(c_state, state)
    x(:, :) = state%xs(:, :)
end subroutine

function ddx_get_x_niter(c_state) bind(C) result(c_niter)
    type(c_ptr), intent(in), value :: c_state
    type(ddx_state_type), pointer :: state
    integer(c_int) :: c_niter
    call c_f_pointer(c_state, state)
    c_niter = state%xs_niter
end function

subroutine ddx_get_s(c_state, nbasis, nsph, s) bind(C)
    type(c_ptr), intent(in), value :: c_state
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nsph, nbasis
    real(c_double), intent(out)  :: s(nbasis, nsph)
    call c_f_pointer(c_state, state)
    s(:, :) = state%s(:, :)
end subroutine

function ddx_get_s_niter(c_state) bind(C) result(c_niter)
    type(c_ptr), intent(in), value :: c_state
    type(ddx_state_type), pointer :: state
    integer(c_int) :: c_niter
    call c_f_pointer(c_state, state)
    c_niter = state%s_niter
end function

! TODO This function will change syntax later
subroutine ddx_get_xi(c_state, c_ddx, ncav, xi) bind(C)
    type(c_ptr), intent(in), value :: c_state, c_ddx
    type(ddx_state_type), pointer :: state
    type(ddx_setup), pointer :: ddx
    integer(c_int), intent(in), value :: ncav
    real(c_double), intent(out)  :: xi(ncav)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddeval_grid_work(ddx%constants%nbasis, ddx%params%ngrid, ddx%params%nsph, &
            & ddx%constants%vgrid, ddx%constants%vgrid_nbasis, one, state%s, zero, xi)
end subroutine


!
! Cosmo
!
subroutine ddx_cosmo_fill_guess(c_ddx, c_state) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddpcm_guess(ddx%params, ddx%constants, state)
end

subroutine ddx_cosmo_solve(c_ddx, c_state, ncav, phi, tol) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: ncav
    real(c_double), intent(in), value :: tol
    real(c_double), intent(in) :: phi(ncav)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddcosmo_solve(ddx%params, ddx%constants, ddx%workspace, state, phi, tol)
end subroutine

subroutine ddx_cosmo_adjoint(c_ddx, c_state, nbasis, nsph, psi, tol) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nbasis, nsph
    real(c_double), intent(in), value :: tol
    real(c_double), intent(in) :: psi(nbasis, nsph)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddcosmo_adjoint(ddx%params, ddx%constants, ddx%workspace, state, psi, tol)
end

subroutine ddx_cosmo_forces(c_ddx, c_state, nbasis, nsph, ncav, phi, gradphi, psi, forces) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: ncav, nbasis, nsph
    real(c_double), intent(in) :: phi(ncav), gradphi(3, ncav), psi(nbasis, nsph)
    real(c_double), intent(out) :: forces(3, nsph)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddcosmo_forces(ddx%params, ddx%constants, ddx%workspace, state, phi, gradphi, psi, forces)
end

!
! PCM
!
subroutine ddx_pcm_fill_guess(c_ddx, c_state) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddpcm_guess(ddx%params, ddx%constants, state)
end

subroutine ddx_pcm_solve(c_ddx, c_state, ncav, phi, tol) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: ncav
    real(c_double), intent(in), value :: tol
    real(c_double), intent(in) :: phi(ncav)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddpcm_solve(ddx%params, ddx%constants, ddx%workspace, state, phi, tol)
end subroutine

subroutine ddx_pcm_adjoint(c_ddx, c_state, nbasis, nsph, psi, tol) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: nbasis, nsph
    real(c_double), intent(in), value :: tol
    real(c_double), intent(in) :: psi(nbasis, nsph)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddpcm_adjoint(ddx%params, ddx%constants, ddx%workspace, state, psi, tol)
end

subroutine ddx_pcm_forces(c_ddx, c_state, nbasis, nsph, ncav, phi, gradphi, psi, forces) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_state
    type(ddx_setup), pointer :: ddx
    type(ddx_state_type), pointer :: state
    integer(c_int), intent(in), value :: ncav, nbasis, nsph
    real(c_double), intent(in) :: phi(ncav), gradphi(3, ncav), psi(nbasis, nsph)
    real(c_double), intent(out) :: forces(3, nsph)
    call c_f_pointer(c_ddx, ddx)
    call c_f_pointer(c_state, state)
    call ddpcm_forces(ddx%params, ddx%constants, ddx%workspace, state, phi, gradphi, psi, forces)
end

end
