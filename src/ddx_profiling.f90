module ddx_profiling
    !! Unified Input/Output handling across the code.
    use ddx_definitions, only: dp
    implicit none
    private

    integer, parameter:: ntimes = 128
    integer:: tcnt = 1
    real(dp):: times(ntimes)

    public:: time_pull, time_push

    contains

    subroutine time_push()
        implicit none
        real(dp):: omp_get_wtime

        if(tcnt <= ntimes) then
            times(tcnt) = omp_get_wtime()
            tcnt = tcnt+1
        else
            write(*, *) 'time_push Cannot push another time in the buffer.'
            stop 1
        end if
    end subroutine

    subroutine time_pull(s)
        implicit none

        character(len=*), intent(in):: s
        real(dp):: elap
        character(len = 2048):: msg

        real(dp):: omp_get_wtime

        if(tcnt > 1) then
            elap = omp_get_wtime() - times(tcnt-1)
            tcnt = tcnt-1
            write(msg, "(3a, ': ', e14.6E2, ' s')") repeat('-', tcnt), '> ', s, elap
            write(*, *) '[ddX-time]', trim(msg)
        else
            write(*, *) 'time_pull Cannot pull any value.'
            stop 1
        end if
    end subroutine

end module ddx_profiling
