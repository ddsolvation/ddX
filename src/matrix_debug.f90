subroutine build_matrix(params, constants, workspace, n, matrix, matvec)
    use ddx_core
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: n
    real(dp), intent(out) :: matrix(n, n)
    external :: matvec
    real(dp), allocatable :: scr1(:), scr2(:)
    integer :: i, j, istat
    allocate(scr1(n), scr2(n), stat=istat)
    if (istat.ne.0) stop 1
    do i = 1, n
        scr1 = zero
        scr1(i) = one
        call matvec(params, constants, workspace, scr1, scr2)
        matrix(i,:) = scr2
    end do
    deallocate(scr1, scr2, stat=istat)
    if (istat.ne.0) stop 1
end subroutine build_matrix

subroutine print_matrix(string, m, n, matrix)
    use ddx_core
    implicit none
    character(*), intent(in) :: string
    integer, intent(in) :: m, n
    real(dp), intent(out) :: matrix(m, n)
    integer :: i, j
    open(file=string, unit=8, form='formatted')
    do i = 1, m
        do j = 1, n
            write(8,'(F20.10$)') matrix(i, j)
        end do
        write(8,*)
    end do
    close(8)
end subroutine print_matrix
