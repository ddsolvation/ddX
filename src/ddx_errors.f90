module ddx_errors

use ddx_definitions

!> ddX type for containing ddx_error information
type ddx_error_type
    !> Flag for ddx_error codes
    integer :: flag = 0
    !> Error message log
    integer :: max_length = 2047
    character(len=2047) :: message
    integer :: message_length = 0
end type ddx_error_type

contains

subroutine print_error(ddx_error)
    implicit none
    type(ddx_error_type), intent(in) :: ddx_error
    write(6, "(A)") ddx_error % message(1:1+ddx_error % message_length)
end subroutine print_error

subroutine check_error(ddx_error)
    implicit none
    type(ddx_error_type), intent(in) :: ddx_error
    if (ddx_error % flag .ne. 0) then
        call print_error(ddx_error)
        stop 1
    end if
end subroutine check_error

subroutine update_error(ddx_error, message)
    implicit none
    type(ddx_error_type), intent(inout) :: ddx_error
    character(len=*) :: message
    integer :: message_start, message_stop, message_length, total_length
    ! update the ddx_error flag by summing one
    ddx_error % flag  = ddx_error % flag + 1
    message_length = len(message)
    total_length = message_length + 2
    message_start = 1 + ddx_error % message_length
    message_stop = message_start + total_length
    ! update the message
    if (message_stop .lt. ddx_error % max_length) then
        ddx_error % message(message_start:message_start) = " "
        ddx_error % message(message_start+1:message_stop-2) = message(1:message_length)
        ddx_error % message(message_stop-1:message_stop-1) = achar(10)
        ddx_error % message_length = ddx_error % message_length + total_length
    else if (message_stop - ddx_error % max_length .gt. 3) then
        message_length = ddx_error % max_length - message_start
        ddx_error % message(message_start:message_start) = " "
        ddx_error % message(message_start+1:ddx_error % max_length-2) = message(1:message_length-3)
        ddx_error % message(ddx_error % max_length-1:ddx_error%max_length-1) = achar(10)
        ddx_error % message_length = ddx_error % message_length + message_length
    end if
end subroutine update_error

subroutine reset_error(ddx_error)
    implicit none
    type(ddx_error_type), intent(inout) :: ddx_error
    integer :: i
    ddx_error % flag = 0
    ddx_error % message_length = 0
    do i = 1, ddx_error % max_length
        ddx_error % message(i:i) = " "
    end do
end subroutine reset_error

end module ddx_errors
