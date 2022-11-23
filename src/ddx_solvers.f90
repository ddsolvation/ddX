!> @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_solvers.f90
!! Iterative solvers for ddX library
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2020-12-17

!> Core routines and parameters of ddX software
module ddx_solvers
! Get the core-routines
use ddx_core

! Disable implicit types
implicit none

! Interfaces
interface
    ! Interface for the matrix-vector product function
    subroutine matvec_interface(params, constants, workspace, x, y)
        !! Add definitions for derived types
        use ddx_core
        !! Inputs
        type(ddx_params_type), intent(in) :: params
        type(ddx_constants_type), intent(in) :: constants
        real(dp), intent(in) :: x(constants % nbasis, params % nsph)
        !! Temporaries
        type(ddx_workspace_type), intent(inout) :: workspace
        ! Output
        real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    end subroutine matvec_interface

    ! Interface for the matrix-vector product function for external jacobi_diis.
    ! note that the dimension of x and y is double for ddlpb.
    subroutine matvec_interface_external(params, constants, workspace, x, y)
        !! Add definitions for derived types
        use ddx_core
        !! Inputs
        type(ddx_params_type), intent(in) :: params
        type(ddx_constants_type), intent(in) :: constants
        real(dp), intent(in) :: x(constants % nbasis, params % nsph, 2)
        !! Temporaries
        type(ddx_workspace_type), intent(inout) :: workspace
        ! Output
        real(dp), intent(out) :: y(constants % nbasis, params % nsph, 2)
    end subroutine matvec_interface_external

    ! Interface for the norm calculating function
    real(dp) function norm_interface(lmax, nbasis, nsph, x)
        !! Add definition for real(dp)
        use ddx_definitions
        !! Inputs
        integer, intent(in) :: lmax, nbasis, nsph
        real(dp), intent(in) :: x(nbasis, nsph)
    end function norm_interface
end interface

contains

!> Jacobi iterative solver
!!
!! @param[in] params: User input parameters for the ddX.
!! @param[in] constants: Precomputed constants for the ddX.
!! @param[inout] workspace: Temporary workspace for the ddX.
!! @param[in] tol: Relative error threshold to stop iterations.
!! @param[in] rhs: Right hand side.
!! @param[inout] x: On input, containts an initial guess to a solution. On
!!      exit, returns 
!! @param[inout] niter: Maximal number of iterations on input. On exit,
!!      contains a number of actually performed matrix-vector products.
!! @param[out] x_rel_diff: Relative difference between old and new values of
!!      the solution on each iteration.
!! @param[in] matvec: Routine that performs
!! @param[in] dm1vec:
!! @param[in] norm_func:
subroutine jacobi_diis(params, constants, workspace, tol, rhs, x, niter, &
        & x_rel_diff, matvec, dm1vec, norm_func)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: tol, rhs(constants % n)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Outputs
    real(dp), intent(inout) :: x(constants % n)
    integer, intent(inout) :: niter
    real(dp), intent(out) :: x_rel_diff(niter)
    !! External procedures
    procedure(matvec_interface) :: matvec, dm1vec
    procedure(norm_interface) :: norm_func
    !! Local variables
    integer :: it, nmat
    real(dp) :: diff, norm, rel_diff
    !! The code
    ! DIIS variable
    nmat = 1
    ! Iterations
    do it = 1, niter
        ! y = rhs - O x
        call matvec(params, constants, workspace, x, workspace % tmp_y)
        workspace % tmp_y = rhs - workspace % tmp_y
        ! x_new = D^-1 y
        call dm1vec(params, constants, workspace, workspace % tmp_y, &
            & workspace % tmp_x_new)
        ! DIIS extrapolation
        if (params % jacobi_ndiis .gt. 0) then
            workspace % tmp_x_diis(:, nmat) = workspace % tmp_x_new
            workspace % tmp_e_diis(:, nmat) = workspace % tmp_x_new - x
            call diis(constants % n, nmat, params % jacobi_ndiis, &
                & workspace % tmp_x_diis, workspace % tmp_e_diis, &
                & workspace % tmp_bmat, workspace % tmp_x_new)
        end if
        ! Norm of increment
        x = workspace % tmp_x_new - x
        diff = norm_func(params % lmax, constants % nbasis, params % nsph, x)
        norm = norm_func(params % lmax, constants % nbasis, params % nsph, &
            & workspace % tmp_x_new)

        ! compute rel_diff using the rule 0/0 = 0
        if (diff.eq.zero .and. norm.eq.zero) then
            rel_diff = zero
        else
            rel_diff = diff / norm
        end if

        if (params % verbose) then
            write(params % iunit, '(A,I3,X,A,E10.2)') &
                & 'Iteration:', it, 'Rel. diff:', rel_diff
            flush(params % iunit)
        end if

        x_rel_diff(it) = rel_diff
        ! Update solution
        x = workspace % tmp_x_new
        ! Check stop condition
        if (rel_diff .le. tol) then
            niter = it
            return
        end if
    end do
    workspace % error_flag = 1
    workspace % error_message = "Jacobi solver did not converge"
    return
endsubroutine jacobi_diis

!> DIIS helper routine
subroutine diis(n, nmat, ndiis, x, e, b, xnew)
    implicit none
    integer, intent(in) :: n, ndiis
    integer, intent(inout) :: nmat
    real(dp), dimension(n,ndiis), intent(inout) :: x, e
    real(dp), dimension(ndiis+1,ndiis+1), intent(inout) :: b
    real(dp), dimension(n), intent(inout) :: xnew
    integer :: nmat1, i, istatus, info
    integer :: j, k
    real(dp), allocatable :: bloc(:,:), cex(:)
    integer,  allocatable :: ipiv(:)
    real(dp), parameter :: zero = 0.0d0, one = 1.0d0
    integer,  parameter :: delta = 10

    if (nmat.ge.ndiis) then
        do j = 2, nmat - delta
            do k = 2, nmat - delta
                b(j, k) = b(j+delta, k+delta)
            end do
        end do
        do j = 1, nmat - delta
            x(:, j) = x(:, j+delta)
            e(:, j) = e(:, j+delta)
        end do
        nmat = nmat - delta
    end if
    nmat1 = nmat + 1

    allocate(bloc(nmat1, nmat1), cex(nmat1), ipiv(nmat1), stat=istatus)

    call makeb(n, nmat, ndiis, e, b)
    bloc = b(1:nmat1, 1:nmat1)
    cex = zero
    cex(1) = one
    call dgesv(nmat1, 1, bloc, nmat1, ipiv, cex, nmat1, info)
!   call gjinv(nmat1, 1, bloc, cex, ok)
    if (info.ne.0) then
        nmat = 1
        return
    end if

    xnew = zero
    do i = 1, nmat
        call daxpy(n,cex(i+1),x(:, i),1,xnew,1)
    end do
    nmat = nmat + 1

    deallocate (bloc, cex, stat=istatus)
end subroutine diis

!> DIIS helper routine
subroutine makeb(n, nmat, ndiis, e, b)
    implicit none
    integer, intent(in) :: n, nmat, ndiis
    real(dp), dimension(n,ndiis), intent(in) :: e
    real(dp), dimension(ndiis+1,ndiis+1), intent(inout) :: b

    integer :: i
    real(dp) :: bij
    real(dp), parameter :: zero = 0.0d0, one = 1.0d0

    ! 1st built
    if (nmat.eq.1) then

        !     [ 0 |  1  ]
        ! b = [ --+---- ]
        !     [ 1 | e*e ]

        b(1,1) = zero
        b(1,2) = one
        b(2,1) = one
        b(2,2) = dot_product(e(:,1),e(:,1))

    ! subsequent builts
    else

        ! first, update the lagrangian line:
        b(nmat+1,1) = one
        b(1,nmat+1) = one

        ! now, compute the new matrix elements:
        !$omp parallel do default(none) shared(nmat,e,b) &
        !$omp private(i,bij) schedule(dynamic)
        do i = 1, nmat - 1
            bij = dot_product(e(:,i), e(:,nmat))
            b(nmat+1, i+1) = bij
            b(i+1, nmat+1) = bij
        end do
        b(nmat+1,nmat+1) = dot_product(e(:,nmat),e(:,nmat))
    end if
end subroutine makeb

!> jacobi_diis_solver_external
!! WARNING: code duplication is bad, but here it makes a very specific sense.
!! in ddLPB, we are solving a linear system with Jacobi DIIS also to apply the
!! preconditioner. As calling a solver recursively is a terrible idea, we
!! create a copy of jacobi_diis that:
!! - takes as an input the size of the system to be solved (no hardcoded
!!   dependence on nsph)
!! - uses copies of the norm subroutine, so that it can also be called with 
!!   a user-given size for arrays
subroutine jacobi_diis_external(params, constants, workspace, n, tol, rhs, x, n_iter, &
          & x_rel_diff, matvec, dm1vec, norm_func)
      type(ddx_params_type),    intent(in)    :: params
      type(ddx_constants_type), intent(inout) :: constants
      type(ddx_workspace_type), intent(inout) :: workspace
      ! Inputs
      integer,                  intent(in)    :: n
      real(dp),                 intent(in)    :: tol
      real(dp),  dimension(n),  intent(in)    :: rhs
      ! Outputs
      real(dp),  dimension(n),  intent(inout) :: x
      integer,                  intent(inout) :: n_iter
      real(dp), intent(out) :: x_rel_diff(n_iter)
      external                                :: matvec, dm1vec
      procedure(norm_interface)               :: norm_func
      ! Local variables
      integer  :: it, nmat, istatus, lenb, nsph_u
      real(dp) :: diff, norm, rel_diff = zero
      logical  :: dodiis
      real(dp), allocatable :: x_new(:), y(:), x_diis(:,:), e_diis(:,:), bmat(:,:)
      ! DIIS extrapolation flag
      dodiis = (params%jacobi_ndiis.gt.0)
      ! extrapolation required
      if (dodiis) then
        ! allocate workspaces
        lenb = params % jacobi_ndiis + 1
        allocate( x_diis(n,params % jacobi_ndiis), e_diis(n,params % jacobi_ndiis), bmat(lenb,lenb) , stat=istatus )
        if (istatus .ne. 0) then
          workspace % error_flag = 1
          workspace % error_message = ' jacobi_diis: [1] failed allocation (diis)'
          return
        endif
        ! initialize the number of points for diis to one.
        ! note that nmat is updated by diis.
        nmat = 1
      endif
      ! allocate workspaces
      allocate( x_new(n), y(n) , stat=istatus )
      if (istatus .ne. 0) then
          workspace % error_flag = 1
          workspace % error_message = ' jacobi_diis: [2] failed allocation (diis)'
          return
      endif
      ! Jacobi iterations
      do it = 1, n_iter
        ! y = rhs - O x
        call matvec(params, constants, workspace, x, y )
        y = rhs - y
        ! x_new = D^-1 y
        call dm1vec(params, constants, workspace, y, x_new)
        ! DIIS extrapolation
        if (dodiis) then
          x_diis(:,nmat) = x_new
          e_diis(:,nmat) = x_new - x
          call diis(n,nmat,params % jacobi_ndiis,x_diis,e_diis,bmat,x_new)
        endif
        !increment
        x = x_new - x
        ! warning: in the following, we use a quick and dirty hack to pass
        ! to the norm function arrays double the size than in non ddlpb calculations
        ! by defining a fake number of spheres with is n/nbasis (i.e., 2*nsph in 
        ! ddlpb).
        nsph_u = n/constants % nbasis
        diff = norm_func(params % lmax, constants % nbasis, nsph_u, x)
        norm = norm_func(params % lmax, constants % nbasis, nsph_u, x_new)
        ! compute rel_diff using the rule 0/0 = 0
        if (diff.eq.zero .and. norm.eq.zero) then
            rel_diff = zero
        else
            rel_diff = diff / norm
        end if
        x_rel_diff(it) = rel_diff
        constants % inner_tol = max(rel_diff*sqrt(tol), tol/100.0d0)
        ! update
        x = x_new
        ! EXIT Jacobi loop here
        if (rel_diff .le. tol) then
            n_iter = it
            return
        end if
      enddo

      workspace % error_flag = 1
      workspace % error_message = ' jacobi_diis: [2] failed allocation (diis)'

endsubroutine jacobi_diis_external

end module ddx_solvers
