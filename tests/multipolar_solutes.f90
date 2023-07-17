program test_multipolar_solutes
use ddx_core
use ddx_multipolar_solutes
implicit none

integer, parameter :: natoms = 24
character(len=255) :: dummy_file_name = ''

type(ddx_type) :: nofmm, fmm
type(ddx_error_type) :: error

real(dp) :: coordinates(3, natoms), radii(natoms), charges(natoms)
real(dp) :: x(natoms), y(natoms), z(natoms)
integer :: info, nproc = 1
real(dp), allocatable :: multipoles(:, :)
real(dp) :: diff
real(dp), parameter :: threshold = 1e-10

! We manually specify two benzene molecules 20 A apart

radii = (/ &
    & 4.00253d0, &
    & 4.00253d0, &
    & 4.00253d0, &
    & 4.00253d0, &
    & 4.00253d0, &
    & 4.00253d0, &
    & 2.99956d0, &
    & 2.99956d0, &
    & 2.99956d0, &
    & 2.99956d0, &
    & 2.99956d0, &
    & 2.99956d0, &
    & 4.00253d0, &
    & 4.00253d0, &
    & 4.00253d0, &
    & 4.00253d0, &
    & 4.00253d0, &
    & 4.00253d0, &
    & 2.99956d0, &
    & 2.99956d0, &
    & 2.99956d0, &
    & 2.99956d0, &
    & 2.99956d0, &
    & 2.99956d0/)

charges = (/ &
    & -0.04192d0, &
    & -0.04192d0, &
    & -0.04198d0, &
    & -0.04192d0, &
    & -0.04192d0, &
    & -0.04198d0, &
    &  0.04193d0, &
    &  0.04193d0, &
    &  0.04197d0, &
    &  0.04193d0, &
    &  0.04193d0, &
    &  0.04197d0, &
    & -0.04192d0, &
    & -0.04192d0, &
    & -0.04198d0, &
    & -0.04192d0, &
    & -0.04192d0, &
    & -0.04198d0, &
    &  0.04193d0, &
    &  0.04193d0, &
    &  0.04197d0, &
    &  0.04193d0, &
    &  0.04193d0, &
    &  0.04197d0/)

coordinates = reshape((/ &
    &  0.00000d0,  2.29035d0,  1.32281d0, &
    &  0.00000d0,  2.29035d0, -1.32281d0, &
    &  0.00000d0,  0.00000d0, -2.64562d0, &
    &  0.00000d0, -2.29035d0, -1.32281d0, &
    &  0.00000d0, -2.29035d0,  1.32281d0, &
    &  0.00000d0,  0.00000d0,  2.64562d0, &
    &  0.00000d0,  4.05914d0,  2.34326d0, &
    &  0.00000d0,  4.05914d0, -2.34326d0, &
    &  0.00000d0,  0.00000d0, -4.68652d0, &
    &  0.00000d0, -4.05914d0, -2.34326d0, &
    &  0.00000d0, -4.05914d0,  2.34326d0, &
    &  0.00000d0,  0.00000d0,  4.68652d0, &
    & 20.00000d0,  2.29035d0,  1.32281d0, &
    & 20.00000d0,  2.29035d0, -1.32281d0, &
    & 20.00000d0,  0.00000d0, -2.64562d0, &
    & 20.00000d0, -2.29035d0, -1.32281d0, &
    & 20.00000d0, -2.29035d0,  1.32281d0, &
    & 20.00000d0,  0.00000d0,  2.64562d0, &
    & 20.00000d0,  4.05914d0,  2.34326d0, &
    & 20.00000d0,  4.05914d0, -2.34326d0, &
    & 20.00000d0,  0.00000d0, -4.68652d0, &
    & 20.00000d0, -4.05914d0, -2.34326d0, &
    & 20.00000d0, -4.05914d0,  2.34326d0, &
    & 20.00000d0,  0.00000d0,  4.68652d0/), &
    (/3, natoms/))

coordinates = coordinates * tobohr
radii = radii * tobohr

x = coordinates(1, :)
y = coordinates(2, :)
z = coordinates(3, :)

! Initialize two ddX model instances, with and without FMMs

call ddinit(natoms, x, y, z, radii, 1, 10, 302, 0, 0, 0, 0, 0.0d0, &
    & 0.1d0, 80.0d0, 0.104d0*tobohr, 0, 200, 25, nproc, &
    & dummy_file_name, nofmm, error)
call check_error(error)

call ddinit(natoms, x, y, z, radii,  1, 10, 302, 0, 1, 20, 20, 0.0d0, &
    & 0.1d0, 80.0d0, 0.104d0*tobohr, 0, 200, 25, nproc, &
    & dummy_file_name, fmm, error)
call check_error(error)

if (nofmm % constants % ncav .ne. fmm % constants % ncav) &
    & call test_error("Mismatch in geometry between FMM and no FMM.")

! Compute the multipoles from the charges

allocate(multipoles(1, natoms), stat=info)
if (info .ne. 0) call test_error("Allocation failed in test_multipolar_solutes")
multipoles(1, :) = charges / sqrt4pi

write(6, *) "Testing up to charges"
call test(fmm, nofmm, multipoles, 0, natoms, threshold, error)
call check_error(error)

deallocate(multipoles, stat=info)
if (info .ne. 0) call test_error("Deallocation failed in test_multipolar_solutes")

! Now we do charges and some random dipoles

allocate(multipoles(4, natoms), stat=info)
if (info .ne. 0) call test_error("Allocation failed in test_multipolar_solutes")
multipoles(1, :) = charges / sqrt4pi
multipoles(2, :) = 0.0d0
multipoles(3, :) = 0.1d0
multipoles(4, :) = -0.1d0

write(6, *) "Testing up to dipoles"
call test(fmm, nofmm, multipoles, 1, natoms, threshold, error)
call check_error(error)

deallocate(multipoles, stat=info)
if (info .ne. 0) call test_error("Deallocation failed in test_multipolar_solutes")

! And finally we do  charges, random dipoles and random quadrupoles

allocate(multipoles(10, natoms), stat=info)
if (info .ne. 0) call test_error("Allocation failed in test_multipolar_solutes")
multipoles(1, :) = charges / sqrt4pi
multipoles(2, :) = 0.0d0
multipoles(3, :) = 0.1d0
multipoles(4, :) = 0.1d0
multipoles(5, :) = -0.1d0
multipoles(6, :) = 0.2d0
multipoles(7, :) = -0.2d0
multipoles(8, :) = 0.3d0
multipoles(9, :) = -0.3d0
multipoles(10, :) = 0.1d0

write(6, *) "Testing up to quadrupoles"
call test(fmm, nofmm, multipoles, 2, natoms, threshold, error)
call check_error(error)

deallocate(multipoles, stat=info)
if (info .ne. 0) call test_error("Deallocation failed in test_multipolar_solutes")

! deallocate

call ddfree(nofmm)
call ddfree(fmm)

contains

subroutine test(fmm, nofmm, multipoles, mmax, nm, threshold, error)
    type(ddx_type), intent(inout) :: fmm, nofmm
    integer, intent(in) :: mmax, nm
    real(dp), intent(in) :: multipoles((mmax+1)**2, nm)
    real(dp), intent(in) :: threshold
    type(ddx_error_type), intent(inout) :: error
    real(dp), allocatable :: phi_nofmm(:), phi_fmm(:), e_nofmm(:, :), &
        & e_fmm(:, :), e_num(:, :), g_nofmm(:, :, :), g_fmm(:, :, :), &
        & g_num(:, :, :)
    integer :: info

    ! Allocate the electrostatic properties arrays

    allocate(phi_nofmm(nofmm % constants % ncav), &
        & phi_fmm(nofmm % constants % ncav), &
        & e_nofmm(3, nofmm % constants % ncav), &
        & e_fmm(3, nofmm % constants % ncav), &
        & e_num(3, nofmm % constants % ncav), &
        & g_nofmm(3, 3, nofmm % constants % ncav), &
        & g_fmm(3, 3, nofmm % constants % ncav), &
        & g_num(3, 3, nofmm % constants % ncav), &
        & stat=info)
    if (info .ne. 0) call test_error("Allocation failed in test_multipolar_solutes")

    ! Compare the potential with and without FMMs

    call build_phi(nofmm % params, nofmm % constants, nofmm % workspace, &
        & multipoles, mmax, phi_nofmm, error)

    call build_phi(fmm % params, fmm % constants, fmm % workspace, &
        & multipoles, mmax, phi_fmm, error)

    diff = norm_inf_1d(nofmm % constants % ncav, phi_nofmm, phi_fmm)
    if (diff .ge. threshold) call test_error("build_phi, phi: FMM and no FMM mismatch")

    ! Compare the field with and without FMMs, furthermore, compare the
    ! potential with the previously obtained one

    call build_e(nofmm % params, nofmm % constants, nofmm % workspace, &
        & multipoles, mmax, phi_nofmm, e_nofmm, error)

    ! Note, here phi_fmm contains the one from the previous step, so we
    ! are actually comparing build_phi with build_e

    diff = norm_inf_1d(nofmm % constants % ncav, phi_nofmm, phi_fmm)
    if (diff .ge. threshold) call test_error("build_e build_phi, phi: mismatch")

    call build_e(fmm % params, fmm % constants, fmm % workspace, &
        & multipoles, mmax, phi_fmm, e_fmm, error)

    diff = norm_inf_1d(nofmm % constants % ncav, phi_nofmm, phi_fmm)
    if (diff .ge. threshold) call test_error("build_e, phi: FMM and no FMM mismatch")

    diff = norm_inf_2d(3, nofmm % constants % ncav, e_nofmm, e_fmm)
    if (diff .ge. threshold) call test_error("build_e, e: FMM and no FMM mismatch")

    ! Compute the numerical field and check

    call numerical_field(multipoles, nofmm % params % csph, mmax, &
        & nofmm % params % nsph, nofmm % constants % ccav, &
        & nofmm % constants % ncav, e_num, error)

    diff = norm_inf_2d(3, nofmm % constants % ncav, e_nofmm, e_num)
    if (diff .ge. threshold) &
        & call test_error("The analytical field does not match the numerical field")

    ! Compare the field gradient with and without FMMs, furthermore, compare
    ! the field and the potential with the previously obtained ones

    call build_g(nofmm % params, nofmm % constants, nofmm % workspace, &
        & multipoles, mmax, phi_nofmm, e_nofmm, g_nofmm, error)

    ! Note, here phi_fmm and e_fmm contain the properties from the previous
    ! step, so we are actually comparing build_g with build_e

    diff = norm_inf_1d(nofmm % constants % ncav, phi_nofmm, phi_fmm)
    if (diff .ge. threshold) call test_error("build_g build_e, phi: mismatch")

    diff = norm_inf_2d(3, nofmm % constants % ncav, e_nofmm, e_fmm)
    if (diff .ge. threshold) call test_error("build_g build_e, e: mismatch")

    call build_g(fmm % params, fmm % constants, fmm % workspace, &
        & multipoles, mmax, phi_fmm, e_fmm, g_fmm, error)

    diff = norm_inf_1d(nofmm % constants % ncav, phi_nofmm, phi_fmm)
    if (diff .ge. threshold) call test_error("build_g, phi: FMM and no FMM mismatch")

    diff = norm_inf_2d(3, nofmm % constants % ncav, e_nofmm, e_fmm)
    if (diff .ge. threshold) call test_error("build_g, e: FMM and no FMM mismatch")

    diff = norm_inf_3d(3, 3, nofmm % constants % ncav, g_nofmm, g_fmm)
    if (diff .ge. threshold) call test_error("build_g, g: FMM and no FMM mismatch")

    ! Compute the numerical field gradient and check

    call numerical_field_gradient(multipoles, nofmm % params % csph, mmax, &
        & nofmm % params % nsph, nofmm % constants % ccav, &
        & nofmm % constants % ncav, g_num, error)

    diff = norm_inf_2d(3, nofmm % constants % ncav, g_nofmm, g_num)
    if (diff .ge. threshold) &
        & call test_error("The analytical field gradient does not match the numerical field gradient")

    ! Deallocate everything

    deallocate(phi_nofmm, phi_fmm, e_nofmm, e_fmm, e_num, g_nofmm, &
        & g_fmm, g_num, stat=info)
    if (info .ne. 0) &
        & call test_error("Deallocation failed in test_multipolar_solutes")

end subroutine test

real(dp) function norm_inf_1d(n, array_1, array_2)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: array_1(n), array_2(n)
    norm_inf_1d = maxval(abs(array_1 - array_2))
end function norm_inf_1d

real(dp) function norm_inf_2d(m, n, array_1, array_2)
    implicit none
    integer, intent(in) :: m, n
    real(dp), intent(in) :: array_1(m, n), array_2(m, n)
    norm_inf_2d = maxval(abs(array_1 - array_2))
end function norm_inf_2d

real(dp) function norm_inf_3d(l, m, n, array_1, array_2)
    implicit none
    integer, intent(in) :: l, m, n
    real(dp), intent(in) :: array_1(l, m, n), array_2(l, m, n)
    norm_inf_3d = maxval(abs(array_1 - array_2))
end function norm_inf_3d

subroutine test_error(message)
    implicit none
    character(*) :: message
    write(6,*) message
    stop 1
end subroutine test_error

subroutine numerical_field(multipoles, cm, mmax, nm, coordinates, &
        & npoints, num_field, error)
    implicit none
    integer, intent(in) :: mmax, nm, npoints
    real(dp), intent(in) :: cm(3, nm), coordinates(3, npoints), &
        & multipoles((mmax + 1)**2, nm)
    real(dp), intent(out) :: num_field(3, npoints)
    type(ddx_error_type), intent(inout) :: error
    real(dp) :: local_coords(3, npoints), v_plus(npoints), v_minus(npoints)
    real(dp), parameter :: delta = 1d-5

    local_coords = coordinates
    local_coords(1, :) = local_coords(1, :) + delta
    call build_phi_dense(multipoles, cm, mmax, nm, v_plus, local_coords, &
        & npoints, error)

    local_coords = coordinates
    local_coords(1, :) = local_coords(1, :) - delta
    call build_phi_dense(multipoles, cm, mmax, nm, v_minus, local_coords, &
        & npoints, error)

    num_field(1, :) = - (v_plus - v_minus)/(2.0d0*delta)

    local_coords = coordinates
    local_coords(2, :) = local_coords(2, :) + delta
    call build_phi_dense(multipoles, cm, mmax, nm, v_plus, local_coords, &
        & npoints, error)

    local_coords = coordinates
    local_coords(2, :) = local_coords(2, :) - delta
    call build_phi_dense(multipoles, cm, mmax, nm, v_minus, local_coords, &
        & npoints, error)

    num_field(2, :) = - (v_plus - v_minus)/(2.0d0*delta)

    local_coords = coordinates
    local_coords(3, :) = local_coords(3, :) + delta
    call build_phi_dense(multipoles, cm, mmax, nm, v_plus, local_coords, &
        & npoints, error)

    local_coords = coordinates
    local_coords(3, :) = local_coords(3, :) - delta
    call build_phi_dense(multipoles, cm, mmax, nm, v_minus, local_coords, &
        & npoints, error)

    num_field(3, :) = - (v_plus - v_minus)/(2.0d0*delta)

end subroutine numerical_field

subroutine numerical_field_gradient(multipoles, cm, mmax, nm, coordinates, &
        & npoints, num_field_gradient, error)
    implicit none
    integer, intent(in) :: mmax, nm, npoints
    real(dp), intent(in) :: cm(3, nm), coordinates(3, npoints), &
        & multipoles((mmax + 1)**2, nm)
    real(dp), intent(out) :: num_field_gradient(3, 3, npoints)
    real(dp) :: local_coords(3, npoints), e_plus(3, npoints), &
        & e_minus(3, npoints), v_scratch(npoints)
    type(ddx_error_type), intent(inout) :: error
    real(dp), parameter :: delta = 1d-5

    local_coords = coordinates
    local_coords(1, :) = local_coords(1, :) + delta
    call build_e_dense(multipoles, cm, mmax, nm, v_scratch, local_coords, &
        & npoints, e_plus, error)

    local_coords = coordinates
    local_coords(1, :) = local_coords(1, :) - delta
    call build_e_dense(multipoles, cm, mmax, nm, v_scratch, local_coords, &
        & npoints, e_minus, error)

    num_field_gradient(1, :, :) = - (e_plus - e_minus)/(2.0d0*delta)

    local_coords = coordinates
    local_coords(2, :) = local_coords(2, :) + delta
    call build_e_dense(multipoles, cm, mmax, nm, v_scratch, local_coords, &
        & npoints, e_plus, error)

    local_coords = coordinates
    local_coords(2, :) = local_coords(2, :) - delta
    call build_e_dense(multipoles, cm, mmax, nm, v_scratch, local_coords, &
        & npoints, e_minus, error)

    num_field_gradient(2, :, :) = - (e_plus - e_minus)/(2.0d0*delta)

    local_coords = coordinates
    local_coords(3, :) = local_coords(3, :) + delta
    call build_e_dense(multipoles, cm, mmax, nm, v_scratch, local_coords, &
        & npoints, e_plus, error)

    local_coords = coordinates
    local_coords(3, :) = local_coords(3, :) - delta
    call build_e_dense(multipoles, cm, mmax, nm, v_scratch, local_coords, &
        & npoints, e_minus, error)

    num_field_gradient(3, :, :) = - (e_plus - e_minus)/(2.0d0*delta)

end subroutine numerical_field_gradient

end program test_multipolar_solutes
