!> Functionality to read input files for the standalone driver

module ddx_input
use ddx_errors, only: ddx_error_type, update_error
use ddx, only: ddinit, ddx_type, closest_supported_lebedev_grid
use ddx_definitions
use ddx_defaults
implicit none

contains

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
!> Read a ddX input file, either in the legacy line-based format or in
!! the keyword format.
!!
!! The format is detected from the first non-blank, non-comment line: a
!! legacy file starts with a single value (the output filename), a
!! keyword file starts with a `key value` pair.
!!
!! Keyword format: `key value`, `key = value` (with or without spaces).
!! `#` or `!` start a comment, keys are case-insensitive, unknown keys
!! are errors. Everything after the `atoms` leader is one sphere
!! per line:
!! charge x y z radius.
!!
!! Lengths are in Angstrom, kappa in inverse Angstrom.
!!
subroutine ddfromfile(fname, ddx_data, tol, charges, ddx_error)
    implicit none
    character(len=*), intent(in) :: fname
    type(ddx_type), intent(out) :: ddx_data
    real(dp), intent(out) :: tol
    real(dp), allocatable, intent(out) :: charges(:)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! Local variables
    character(len=512) :: line, key, val, msg, probe, error_text
    character(len=255) :: output_filename
    integer :: unit, io_info, lineno, index_comment, index_splitting, &
        & k, nsph, natoms
    integer :: nproc, model, lmax, ngrid, force, fmm, pm, pl, &
        & matvecmem, maxiter, jacobi_ndiis, switching
    real(dp) :: eps, se, eta, kappa
    real(dp), allocatable :: csph(:, :), radii(:)
    logical :: in_atoms, have_model, have_eps, have_se

    open(newunit=unit, file=fname, status="old", action="read", &
        & form="formatted", iostat=io_info)
    if (io_info .ne. 0) then
        call update_error(ddx_error, "Cannot open input file " // &
            & trim(fname))
        return
    end if
    do
        read(unit, "(a)", iostat=io_info) line
        index_comment = scan(line, "#!")
        if (index_comment .gt. 0) line(index_comment:) = " "
        if (io_info .ne. 0) exit
        line = adjustl(line)
        if (len_trim(line) .eq. 0 .or. line(1:1) .eq. "#") cycle
        exit
    end do
    index_splitting = 0
    if (io_info .eq. 0) then
        ! `key=val`, `key = val` and `key val` all count as two tokens
        probe = line
        do k = 1, len(probe)
            if (probe(k:k) .eq. achar(9) .or. probe(k:k) .eq. "=") &
                & probe(k:k) = " "
        end do
        probe = adjustl(probe)
        index_splitting = index(trim(probe), " ")
    end if
    if (index_splitting .eq. 0) then
        ! single token (or empty file): legacy format
        close(unit)
        call ddfromfile_legacy(fname, ddx_data, tol, charges, ddx_error)
        return
    end if
    rewind(unit)

    output_filename = ""
    nproc = default_nproc
    lmax = default_lmax
    ngrid = default_ngrid
    eta = default_eta
    kappa = default_kappa
    matvecmem = default_incore
    tol = default_tol
    maxiter = default_maxiter
    jacobi_ndiis = default_jacobi_ndiis
    force = default_force
    fmm = default_enable_fmm
    pm = default_pm
    pl = default_pl
    switching = default_switching
    model = 0
    eps = -one
    nsph = 0
    natoms = 0
    have_model = .false.
    have_eps = .false.
    have_se = .false.
    in_atoms = .false.

    lineno = 0
    do
        read(unit, "(a)", iostat=io_info) line
        if (io_info .ne. 0) exit
        lineno = lineno + 1
        index_comment = scan(line, "#!")
        if (index_comment .gt. 0) line(index_comment:) = " "
        do k = 1, len(line)
            if (line(k:k) .eq. achar(9)) line(k:k) = " "
        end do
        line = adjustl(line)
        if (len_trim(line) .eq. 0) cycle

        error_text = ""
        io_info = 0
        if (in_atoms) then
            natoms = natoms + 1
            read(line, *, iostat=io_info) charges(natoms), &
                & csph(1, natoms), csph(2, natoms), csph(3, natoms), &
                & radii(natoms)
            if (io_info .ne. 0) &
                & error_text = "expected `charge x y z radius`"
        else
            index_splitting = scan(line, " =")
            if (index_splitting .eq. 0) then
                key = line
                val = ""
            else
                key = line(1:index_splitting-1)
                val = adjustl(line(index_splitting+1:))
                if (val(1:1) .eq. "=") val = adjustl(val(2:))
            end if
            do k = 1, len_trim(key)
                if (key(k:k) .ge. "A" .and. key(k:k) .le. "Z") &
                    & key(k:k) = achar(iachar(key(k:k)) + 32)
            end do

            select case (key)
            case ("output", "logfile")
                output_filename = trim(val)
            case ("nproc")
                read(val, *, iostat=io_info) nproc
            case ("model")
                do k = 1, len_trim(val)
                    if (val(k:k) .ge. "A" .and. val(k:k) .le. "Z") &
                        & val(k:k) = achar(iachar(val(k:k)) + 32)
                end do
                select case (trim(val))
                case ("cosmo")
                    model = 1
                case ("pcm")
                    model = 2
                case ("lpb")
                    model = 3
                case default
                    read(val, *, iostat=io_info) model
                end select
                have_model = .true.
            case ("lmax")
                read(val, *, iostat=io_info) lmax
            case ("ngrid")
                read(val, *, iostat=io_info) ngrid
            case ("eps")
                read(val, *, iostat=io_info) eps
                have_eps = .true.
            case ("se", "shift")
                read(val, *, iostat=io_info) se
                have_se = .true.
            case ("eta")
                read(val, *, iostat=io_info) eta
            case ("kappa")
                read(val, *, iostat=io_info) kappa
            case ("matvecmem", "incore")
                read(val, *, iostat=io_info) matvecmem
            case ("tol")
                read(val, *, iostat=io_info) tol
            case ("maxiter")
                read(val, *, iostat=io_info) maxiter
            case ("jacobi_ndiis")
                read(val, *, iostat=io_info) jacobi_ndiis
            case ("force")
                read(val, *, iostat=io_info) force
            case ("fmm")
                read(val, *, iostat=io_info) fmm
            case ("pm")
                read(val, *, iostat=io_info) pm
            case ("pl")
                read(val, *, iostat=io_info) pl
            case ("switching")
                read(val, *, iostat=io_info) switching
            case ("nsph")
                read(val, *, iostat=io_info) nsph
                if (io_info .eq. 0 .and. nsph .le. 0) &
                    & error_text = "`nsph` must be positive"
            case ("atoms")
                if (nsph .le. 0) then
                    error_text = "`nsph` must be given before `atoms`"
                else
                    allocate(charges(nsph), csph(3, nsph), radii(nsph))
                    in_atoms = .true.
                end if
            case default
                error_text = "unknown key `" // trim(key) // "`"
            end select
            if (io_info .ne. 0 .and. len_trim(error_text) .eq. 0) &
                & error_text = "invalid value for `" // trim(key) // "`"
        end if

        if (len_trim(error_text) .gt. 0) then
            write(msg, "(a,a,i0,a,a)") trim(fname), ", line ", lineno, &
                & ": ", trim(error_text)
            call update_error(ddx_error, trim(msg))
            exit
        end if
        ! all atoms read: ignore whatever follows
        if (in_atoms .and. natoms .eq. nsph) exit
    end do
    close(unit)
    if (ddx_error % flag .ne. 0) return

    if (.not. have_model) &
        & call update_error(ddx_error, "missing key `model`")
    if (have_model .and. (model .lt. 1 .or. model .gt. 3)) &
        & call update_error(ddx_error, "`model` must be cosmo," // &
            & " pcm, lpb or 1..3")
    if (.not. have_eps) &
        & call update_error(ddx_error, "missing key `eps`")
    if (.not. in_atoms) then
        call update_error(ddx_error, "missing `atoms` section")
    else if (natoms .lt. nsph) then
        write(msg, "(a,i0,a,i0,a)") "`nsph` is ", nsph, " but only ", &
            & natoms, " atom lines were found"
        call update_error(ddx_error, trim(msg))
    end if
    if (nproc .lt. 0) &
        & call update_error(ddx_error, "`nproc` must be non-negative")
    if (lmax .lt. 0) &
        & call update_error(ddx_error, "`lmax` must be non-negative")
    if (ngrid .lt. 0) &
        & call update_error(ddx_error, "`ngrid` must be non-negative")
    if (have_eps .and. eps .lt. zero) &
        & call update_error(ddx_error, "`eps` must be non-negative")
    if (se .lt. -one .or. se .gt. one) &
        & call update_error(ddx_error, "`se` must be in [-1, 1]")
    if (eta .lt. zero .or. eta .gt. one) &
        & call update_error(ddx_error, "`eta` must be in [0, 1]")
    if (kappa .lt. zero) &
        & call update_error(ddx_error, "`kappa` must be non-negative")
    if (matvecmem .lt. 0 .or. matvecmem .gt. 1) &
        & call update_error(ddx_error, "`matvecmem` must be 0 or 1")
    if (tol .lt. 1d-14 .or. tol .gt. one) &
        & call update_error(ddx_error, "`tol` must be in [1d-14, 1]")
    if (maxiter .le. 0) &
        & call update_error(ddx_error, "`maxiter` must be positive")
    if (jacobi_ndiis .lt. 0) &
        & call update_error(ddx_error, "`jacobi_ndiis` must be" // &
            & " non-negative")
    if (force .lt. 0 .or. force .gt. 1) &
        & call update_error(ddx_error, "`force` must be 0 or 1")
    if (fmm .lt. 0 .or. fmm .gt. 1) &
        & call update_error(ddx_error, "`fmm` must be 0 or 1")
    if (pm .lt. 0 .or. pl .lt. 0) &
        & call update_error(ddx_error, "`pm` and `pl` must be" // &
            & " non-negative")
    if (ddx_error % flag .ne. 0) return

    ! The default switching shift depends on the model
    if (.not.have_se) se = default_se(model)
    csph = csph * tobohr
    radii = radii * tobohr
    kappa = kappa / tobohr
    call closest_supported_lebedev_grid(ngrid)

    call ddinit(model, nsph, csph, radii, eps, ddx_data, ddx_error, &
        & force=force, kappa=kappa, eta=eta, shift=se, lmax=lmax, &
        & ngrid=ngrid, incore=matvecmem, maxiter=maxiter, &
        & jacobi_ndiis=jacobi_ndiis, enable_fmm=fmm, pm=pm, pl=pl, &
        & nproc=nproc, logfile=output_filename, switching=switching)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "ddinit returned an error," // &
            & " exiting")
        return
    end if
end subroutine ddfromfile

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
subroutine ddfromfile_legacy(fname, ddx_data, tol, charges, ddx_error)
    implicit none
    character(len=*), intent(in) :: fname
    type(ddx_type), intent(out) :: ddx_data
    real(dp), intent(out) :: tol
    real(dp), allocatable, intent(out) :: charges(:)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! Local variables
    integer :: nproc, model, lmax, ngrid, force, fmm, pm, pl, &
        & nsph, i, matvecmem, maxiter, jacobi_ndiis, &
        & istatus
    real(dp) :: eps, se, eta, kappa
    real(dp), allocatable :: csph(:, :), radii(:)
    character(len=255) :: output_filename

    write(6, *) "Warning: the legacy input file is deprecated, not" // &
        & " every parameter can be set with it."

    !! Read all the parameters from the file
    ! Open a configuration file
    open(unit=100, file=fname, form="formatted", access="sequential")
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

    ! Initialize ddx_data object
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
end subroutine ddfromfile_legacy

end module ddx_input
