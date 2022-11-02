module ddx_properties
use ddx_core
implicit none

contains

!! Given a set of coordinates, partition them to the corresponding spheres and
!! return the corresponding indexes into an array. If a coordinate is outside
!! the cavity set the corresponding index to 0.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[in] coordinates: coordinates of the points, size (3, npoints)
!! @param[out] partition: partition over the spheres, size (npoints)
!! @param[in] npoints: number of points
subroutine partition_coordinates(params, constants, coordinates, partition, &
        & npoints)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: npoints
    real(dp), intent(in) :: coordinates(3, npoints)
    integer, intent(out) :: partition(npoints)
    integer :: i, isph
    real(dp) :: d(3), r2
    !TODO: make this linear scaling
    do i = 1, npoints
        partition(i) = 0
        do isph = 1, params % nsph
            d = coordinates(:, i) - params % csph(:, isph)
            r2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
            if (r2 .lt. params % rsph(isph)**2) then
                partition(i) = isph
                exit
            end if
        end do
    end do
end subroutine partition_coordinates

!! Given a set of coordinates, compute the reaction potential at them.
!! If a given point is outside the cavity return zero as the potential
!! at that point.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[in] state: ddx state
!! @param[in] coordinates: coordinates of the points, size (3, npoints)
!! @param[out] potential: reaction potential at the points, size (npoints)
!! @param[in] npoints: number of points
subroutine reaction_potential(params, constants, state, coordinates, &
        & potential, npoints)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(in), target :: state
    integer, intent(in) :: npoints
    real(dp), intent(in) :: coordinates(3, npoints)
    real(dp), intent(out) :: potential(npoints)
    integer, allocatable :: partition(:)
    integer :: info, i, isph
    real(dp) :: d(3), work(params % lmax + 1)
    real(dp), pointer :: x(:, :)

    allocate(partition(npoints), stat=info)
    if (info .ne. 0) stop 1

    call partition_coordinates(params, constants, coordinates, partition, &
        & npoints)

    if (params % model .eq. 1 .or. params % model .eq. 2) then
        x => state % xs
    else if (params % model .eq. 3) then
        x => state % x_lpb(:, :, 1)
    else
        return
    end if

    do i = 1, npoints
        isph = partition(i)
        potential(i) = zero
        if (isph .gt. 0) then
            d = coordinates(:, i) - params % csph(:, isph)
            call fmm_l2p_work(d, params % rsph(isph), params % lmax, &
                & constants % vscales_rel, one, x(:, isph), one, &
                & potential(i), work)
        end if
    end do

    deallocate(partition, stat=info)
    if (info .ne. 0) stop 1
end subroutine reaction_potential

!! Given a set of sphere indexes, compute the reaction potential at the
!! centers of the given spheres.
!! @param[in] params: ddx parameters
!! @param[in] state: ddx state
!! @param[in] indexes: sphere indexes, size (npoints)
!! @param[out] potential: reaction potential at the points, size (npoints)
!! @param[in] npoints: number of points
subroutine reaction_potential_at_centers(params, state, indexes, potential, &
        & npoints)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_state_type), intent(in), target :: state
    integer, intent(in) :: npoints
    integer, intent(in) :: indexes(npoints)
    real(dp), intent(out) :: potential(npoints)
    integer :: i, isph
    real(dp), pointer :: x(:, :)

    if (params % model .eq. 1 .or. params % model .eq. 2) then
        x => state % xs
    else if (params % model .eq. 3) then
        x => state % x_lpb(:, :, 1)
    else
        return
    end if

    do i = 1, npoints
        isph = indexes(i)
        potential(i) = x(1, isph)*sqrt4pi
    end do
end subroutine reaction_potential_at_centers

!! Given a set of sphere indexes, compute the reaction field at the
!! centers of the given spheres. This is done using finite differences.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[in] state: ddx state
!! @param[in] indexes: sphere indexes, size (npoints)
!! @param[out] field: reaction field at the points, size (3, npoints)
!! @param[in] npoints: number of points
subroutine num_reaction_field_at_centers(params, constants, state, indexes, field, &
        & npoints)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(in), target :: state
    integer, intent(in) :: npoints
    integer, intent(in) :: indexes(npoints)
    real(dp), intent(out) :: field(3, npoints)
    integer :: i, isph
    real(dp), pointer :: x(:, :)
    real(dp) :: c(3, 6), v(6), step
    step = 0.01d0

    if (params % model .eq. 1 .or. params % model .eq. 2) then
        x => state % xs
    else if (params % model .eq. 3) then
        x => state % x_lpb(:, :, 1)
    else
        return
    end if

    do i = 1, npoints
        isph = indexes(i)
        c(:, 1) = params % csph(:, isph)
        c(:, 2) = params % csph(:, isph)
        c(:, 3) = params % csph(:, isph)
        c(:, 4) = params % csph(:, isph)
        c(:, 5) = params % csph(:, isph)
        c(:, 6) = params % csph(:, isph)
        c(1, 1) = c(1, 1) - step
        c(1, 2) = c(1, 2) + step
        c(2, 3) = c(2, 3) - step
        c(2, 4) = c(2, 4) + step
        c(3, 5) = c(3, 5) - step
        c(3, 6) = c(3, 6) + step
        call reaction_potential(params, constants, state, c, v, 6)
        field(1, i) = (v(2) - v(1))/(2.0d0*step)
        field(2, i) = (v(4) - v(3))/(2.0d0*step)
        field(3, i) = (v(6) - v(5))/(2.0d0*step)
    end do
end subroutine num_reaction_field_at_centers

!! Given a set of sphere indexes, compute the reaction field at the
!! centers of the given spheres.
!! @param[in] params: ddx parameters
!! @param[in] state: ddx state
!! @param[in] indexes: sphere indexes, size (npoints)
!! @param[out] field: reaction field at the points, size (3, npoints)
!! @param[in] npoints: number of points
subroutine reaction_field_at_centers(params, state, indexes, field, npoints)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_state_type), intent(in), target :: state
    integer, intent(in) :: npoints
    integer, intent(in) :: indexes(npoints)
    real(dp), intent(out) :: field(3, npoints)
    integer :: i, isph
    real(dp), pointer :: x(:, :)
    real(dp) :: f1, f2

    if (params % model .eq. 1 .or. params % model .eq. 2) then
        x => state % xs
    else if (params % model .eq. 3) then
        x => state % x_lpb(:, :, 1)
    else
        return
    end if

    f1 = sqrt(fourpi/three)
    do i = 1, npoints
        isph = indexes(i)
        f2 = f1/params % rsph(isph)
        field(1, i) = f2*x(4, isph)
        field(2, i) = f2*x(2, isph)
        field(3, i) = f2*x(3, isph)
    end do
end subroutine reaction_field_at_centers

end module ddx_properties
