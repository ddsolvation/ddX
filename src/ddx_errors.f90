module ddx_errors

use ddx_definitions

!> ddX type for containing error information
type ddx_error_type
    !> Flag for error codes
    integer :: flag = 0
    !> Error message log
    integer :: max_length = 2047
    character(len=2047) :: message
    integer :: message_length = 0
end type ddx_error_type

contains

subroutine print_error(error)
    implicit none
    type(ddx_error_type), intent(in) :: error
    write(6, "(A)") error % message(1:1+error % message_length)
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
    integer :: message_start, message_stop, message_length, total_length
    ! update the error flag by summing one
    error % flag  = error % flag + 1
    message_length = len(message)
    total_length = message_length + 2
    message_start = 1 + error % message_length
    message_stop = message_start + total_length
    ! update the message
    if (message_stop .lt. error % max_length) then
        error % message(message_start:message_start) = " "
        error % message(message_start+1:message_stop-2) = message(1:message_length)
        error % message(message_stop-1:message_stop-1) = achar(10)
        error % message_length = error % message_length + total_length
    else if (message_stop - error % max_length .gt. 3) then
        message_length = error % max_length - message_start
        error % message(message_start:message_start) = " "
        error % message(message_start+1:error % max_length-2) = message(1:message_length-3)
        error % message(error % max_length-1:error%max_length-1) = achar(10)
        error % message_length = error % message_length + message_length
    end if
end subroutine update_error

subroutine reset_error(error)
    implicit none
    type(ddx_error_type), intent(inout) :: error
    integer :: i
    error % flag = 0
    error % message_length = 0
    do i = 1, error % max_length
        error % message(i:i) = " "
    end do
end subroutine reset_error

end module ddx_errors
