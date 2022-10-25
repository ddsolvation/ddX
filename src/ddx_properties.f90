module ddx_properties
use ddx_core
implicit none

contains

!! Given a set of coordinates, partition them to the corresponding spheres and
!! return the corresponding indexes into an array. If a coordinate is outside
!! the cavity set the corresponding index to -1.
!! @param[in] params
!! @param[in] constants
!! @param[in] coordinates
!! @param[out] partition
!! @param[in] npoints
subroutine partition_coordinates(params, constants, coordinates, partition, &
        & npoints)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: npoints
    real(dp), intent(in) :: coordinates(3, npoints)
    integer, intent(out) :: partition(npoints)
    !TODO: make this linear scaling
end subroutine partition_coordinates

subroutine reaction_potential(params, constants, state, coordinates, &
        & potential, npoints)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(in) :: state
    integer, intent(in) :: npoints
    real(dp), intent(in) :: coordinates(3, npoints)
    real(dp), intent(out) :: potential(npoints)
end subroutine reaction_potential

subroutine reaction_potential_at_centers(params, state, indexes, potential, &
        & npoints)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_state_type), intent(in) :: state
    integer, intent(in) :: npoints
    integer, intent(in) :: indexes(npoints)
    real(dp), intent(out) :: potential(npoints)
    integer :: i, isph
    real(dp), pointer :: x(:, :)

    if (params % model .eq. 1 .or. params % model .eq. 2) then
        do i = 1, npoints
            isph = indexes(i)
            potential(i) = state % xs(1, isph)*sqrt4pi
        end do
    else if (params % model .eq. 3) then
        do i = 1, npoints
            isph = indexes(i)
            potential(i) = state % x_lpb(1, isph, 1)*sqrt4pi
        end do
    end if
end subroutine reaction_potential_at_centers

subroutine reaction_field_at_centers()
end subroutine reaction_field_at_centers

end module ddx_properties
