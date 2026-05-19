program fpm_wrapper_standalone
implicit none
integer :: i, status, failed
integer, parameter :: ncommands = 8
character(len=1000) :: commands(ncommands)
character(len=8) :: status_string
integer :: results(ncommands)

commands = [character(len=1000) :: &
    "./tests/standalone_tests/run_test.py cosmo --fpm", &
    "./tests/standalone_tests/run_test.py cosmo_fmm --fpm", &
    "./tests/standalone_tests/run_test.py cosmo_incore --fpm", &
    "./tests/standalone_tests/run_test.py pcm --fpm", &
    "./tests/standalone_tests/run_test.py pcm_fmm --fpm", &
    "./tests/standalone_tests/run_test.py pcm_incore --fpm", &
    "./tests/standalone_tests/run_test.py lpb --fpm", &
    "./tests/standalone_tests/run_test.py lpb_fmm --fpm"]


results(:) = 0
do i = 1, ncommands
    call system(commands(i), status)
    if (status.ne.0) then
        results(i) = 1
    end if
end do

failed = sum(results)

write(6,"(A)") "OVERVIEW:"
write(6,"(A,I0,A,I0,A)") "Results: ", ncommands - failed, "/", ncommands, " passed"
do i = 1, ncommands
    if (results(i).eq.0) then
        status_string = "SUCCESS:"
    else
        status_string = "FAILURE:"
    end if
    write(6,"(A,X,A)") trim(status_string), trim(commands(i))
end do

if (failed.gt.0) call exit(1)
call exit(0)

end program fpm_wrapper_standalone
