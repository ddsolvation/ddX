module ddx_errors

!> ddX type for containing error information
type ddx_error_type
    !> Flag for error codes
    integer :: flag = 0
    !> Error message log
    integer :: max_length = 2000
    character(len=2000) :: message
    integer :: message_length = 1
end type ddx_error_type

contains

subroutine print_error(error)
    implicit none
    type(ddx_error_type), intent(in) :: error
    write(6, "(A)") error % message(0:error % message_length-2)
end subroutine print_error

subroutine check_error(error)
    implicit none
    type(ddx_error_type), intent(in) :: error
    if (error % flag .ne. 0) then
        call print_error(error)
        stop error % flag
    end if
end subroutine check_error

subroutine update_error(error, message)
    implicit none
    type(ddx_error_type), intent(inout) :: error
    character(len=*) :: message
    integer :: message_start, message_stop
    ! update the error flag by summing one
    error % flag = error % flag + 1
    message_start = error % message_length
    message_stop = message_start + len(message) + 3
    ! update the message if there is still space left
    if (message_stop .lt. error % max_length) then
        error % message(message_start:message_start) = " "
        error % message(message_start + 1:message_stop - 2) = message
        error % message(message_stop - 1:message_stop - 1) = achar(13)
        error % message(message_stop:message_stop) = achar(10)
        error % message_length = message_stop + 1
    end if
end subroutine update_error

subroutine reset_error(error)
    implicit none
    type(ddx_error_type), intent(inout) :: error
    integer :: i
    error % flag = 0
    do i = 1, error % max_length
        error % message(i:i) = " "
    end do
end subroutine reset_error

end module ddx_errors
