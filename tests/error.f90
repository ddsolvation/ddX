program test_ddx_error
use ddx_errors
implicit none

type(ddx_error_type) :: error
character(len=40) :: ref_text

! we artificially set the max_length to a small value
error % max_length = 40

call update_error(error, "<1 block of 15>")

ref_text(1:16) = " <1 block of 15>"
ref_text(17:17) = achar(10)

if (error % message(1:error % message_length) .ne. ref_text(1:error%message_length)) &
    & stop 1

call update_error(error, "<2 block of 15>")
ref_text(18:33) = " <2 block of 15>"
ref_text(34:34) = achar(10)

write(6,"(A,A,A)") error % message(1:error % message_length), "<"
write(6,"(A,A,A)") ref_text(1:error % message_length), "<"

if (error % message(1:error % message_length) .ne. ref_text(1:error%message_length)) &
    & stop 2

call update_error(error, "<3 block of 15>")
ref_text(35:38) = " <3 "
ref_text(39:39) = achar(10)

if (error % message(1:error % message_length) .ne. ref_text(1:error%message_length)) &
    & stop 3

call print_error(error)
call reset_error(error)
call check_error(error)

contains

end program test_ddx_error
