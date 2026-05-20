program fpm_wrapper
implicit none
integer :: i, status, failed
integer, parameter :: ncommands = 9
character(len=1000) :: commands(ncommands)
character(len=8) :: status_string
integer :: results(ncommands)

commands = [character(len=1000) :: &
    "fpm run --target force -- tests/Input_force.txt", &
    "fpm run --target test_ddx_driver -- tests/data/ddpcm_force_fmm.in tests/data/ddpcm_force_fmm.out 1E-12", &
    "fpm run --target test_ddx_driver -- tests/data/ddcosmo_force_fmm.in tests/data/ddcosmo_force_fmm.out 1E-12", &
    "fpm run --target force_ddlpb -- tests/data/ddlpb_force.txt", &
    "fpm run --target ddlpb_esolv -- tests/data/ddlpb_force.txt", &
    "fpm run --target matrix_derivatives -- tests/data/ddlpb_force.txt", &
    "fpm run --target matrix_adjoint -- tests/data/ddlpb_force.txt", &
    "fpm run --target matrix_solvers -- tests/data/ddlpb_force.txt", &
    "fpm run --target test_gradients -- tests/Input_cosmo_small.txt"]


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

end program fpm_wrapper
