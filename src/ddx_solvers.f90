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
!! @param[out] info:
subroutine jacobi_diis(params, constants, workspace, tol, rhs, x, niter, &
        & x_rel_diff, matvec, dm1vec, norm_func, info)
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
    integer, intent(out) :: info
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

        x_rel_diff(it) = rel_diff
        ! Update solution
        x = workspace % tmp_x_new
        ! Check stop condition
        if (rel_diff .le. tol) then
            info = 0
            niter = it
            return
        end if
    end do
    info = 1
endsubroutine jacobi_diis

!subroutine diis_new(params, constants, workspace, nmat, xnew, x)
!    implicit none
!    type(ddx_params_type), intent(in) :: params
!    type(ddx_constants_type), intent(in) :: constants
!    type(ddx_workspace_type), intent(inout) :: workspace
!    real(dp), intent(inout) :: xnew(constants % n)
!    real(dp), intent(in) :: x(constants % n)
!    integer, intent(inout) :: nmat
!
!    ! add new quantities
!    workspace % tmp_x_diis(:, nmat) = xnew
!    workspace % tmp_e_diis(:, nmat) = xnew - x
!
!    ! if we run out of space, reset
!    if (nmat.ge.(2*params % jacobi_ndiis)) then
!        workspace % tmp_x_diis(:, 0:nmat) = &
!            & workspace % tmp_x_diis(:, nmat:2*nmat)
!        workspace % tmp_e_diis(:, 0:nmat) = &
!            & workspace % tmp_e_diis(:, nmat:2*nmat)
!        workspace % tmp_b_mat(0:nmat, 0:nmat) = &
!            & workspace % tmp_b_mat(nmat:2*nmat, nmat:2*nmat)
!    nmat = workspace % jacobi_ndiis
!    end if
!
!    ! build B matrix of DIIS
!    call makeb(workspace % n, nmat, 2*workspace % jacobi_ndiis, &
!        & workspace % tmp_e_diis, workspace % tmp_b_mat)
!
!end subroutine diis_new


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
    if (istatus.ne.0) then
        write(*,*) 'diis: allocation failed!'
        stop
    end if

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
    !$omp parallel do default(none) shared(nmat,x,cex) &
    !$omp private(i) schedule(dynamic) reduction(+:xnew)
    do i = 1, nmat
        xnew = xnew + cex(i+1)*x(:, i)
    end do
    nmat = nmat + 1

    deallocate (bloc, cex, stat=istatus)
    if (istatus.ne.0) then
        write(*,*) 'diis: deallocation failed!'
        stop
    end if
end subroutine diis

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

!*********************************************************
! GMRESR algorithm to solve linear system Ax = b
!
!  author:
!  m.botchev, utrecht university, december 1996
!  report bugs to botchev@cwi.nl or botchev@excite.com

! Copyright (c) 1996 by M.A.Botchev
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.

! Details to the algorithm may be found in
!  H.A. van der Vorst, !. Vuik, "GMRESR: a Family of Nested GMRES
!  Methods", Num. Lin. Alg. Appl., vol. 1(4), 369--386 (1994)

! parameter list:
! oktest is TRUE if intermediate residuals should be printed
! n      == INTEGER size of problem
! j      == INTEGER truncation parameter (work with j last vectors)
! mgmres == INTEGER dimension of the envoked GMRES
!           if mgmres.eq.0 then we get GCR algorithm, the simplest
!           version of GMRESR 
! b      == DOUBLE PRECISION righthand side vector
! x      == DOUBLE PRECISION initial guess on input,
!           (approximate) solution on output
! work   == DOUBLE PRECISION work space 1 of size n x (2*j+mgmres+2)
! eps    == DOUBLE PRECISION tolerance of stopping criterion. 
!           process is stopped as soon as 
!           (1) residual norm has been dumped by factor eps, 
!           i.e.  ||res|| / ||b|| <= eps   OR
!           (2) maximum number of iterations maxit has been performed
! stc    == CHARACTER*3
!           Determine stopping criterion (||.|| denotes the 2-norm):
!           stc='rel'    -- relative stopping crit.: ||res|| < eps*||res0||
!           stc='abs'    -- absolute stopping crit.: ||res|| < eps
! maxits == INTEGER max. no. outer_iterative_steps/truncation_length on input
!           on output it is the actual number of total iterative steps   
! resid  == DOUBLE PRECISION residual measure (depends on stopping criterion)
!           achieved on output 
! iflag  == INTEGER on output 0 - solution found within tolerance
!                             1 - no convergence within niter

subroutine gmresr(params, constants, workspace, eps, b, x, &
                & niter, resid, matvec, info)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: eps, b(constants % n)
    !! Outputs
    real(dp), intent(inout) :: x(constants % n)
    integer, intent(inout) :: niter
    !real(dp), intent(out) :: x_rel_diff(niter)
    integer, intent(out) :: info
    !! Temporary buffers
    type(ddx_workspace_type), intent(inout) :: workspace
    !! External procedures
    procedure(matvec_interface) :: matvec
    !! Local variables

!     list of variables: arrays in alphabetical order,
!     then other variables in alphabetical order

      integer i,its,nits,itsinn,k,n,j

      double precision alpha, alphai, cknorm, ckres, ddot, dnrm2, &
                       & epsinn, res0,resnor,resid

! distribute work space work(n,*) among some virtual variables;
! namely, we think of columns of work as being occupied by 
! c(n,0:j-1), u(n,0:j-1), resid(n), workgmr(n,params % gmresr_dim+1)
! therefore we define "shifts"

      integer c, u, workres, workgmre
! ----------------------------------------------------------------------------

      ! Read j and n from the parameters
      j = params % gmresr_j
      n = constants % n
!     c occupies j columns 0...j-1:
      c = 0 
!     u occupies j columns j...2*j-1:
      u = params % gmresr_j
!     resid occupies 1 column No. 2*j:
      workres = 2*j    
!     workgmre occupies gmres_dim+1 columns 2*j+1...2*j+gmres_dim+1:
      workgmre = 2*j+1 
!     so that we can access, say, to the (k)th column of the virtual
!     array c(n,0:j-1) as work(1,c+k),
!     virtual residual vector resid(n) is work(1,workres) and so on ...

! ***Furthermore, we build sequences c_k and u_k, k = 0,...,m-1
! but we store only the last j vectors in the following way:
! Assume j=3, then
! --------------------------------------------------------------
!  k    |  number of column of work | columns of work which are vectors
!       |  in which we store c_k    |  we actually store and use
!  0    |           0               |   c_0             u_0            ...
!  1    |           1               |   c_0  c_1        u_0  u_1       ...
!  2    |           2               |   c_0  c_1  c_2   u_0  u_1  u_2  ...
!  3    |           0               |   c_3  c_1  c_2   u_3  u_1  u_2  ...
!  4    |           1               |   c_3  c_4  c_2   u_3  u_4  u_2  ... 
!  5    |           2               |   c_3  c_4  c_5   u_3  u_4  u_5  ...
!  6    |           0               |   c_6  c_4  c_5   u_6  u_4  u_5  ...
! ...   |           ...             |      ...               ...
! This mapping is performed by function mod(k,j)
!

!     Reset iteration counter
      nits= 0
      its = 0
!     Calculate (initial) residual norm
      call matvec(params, constants, workspace, x, &
          & workspace % tmp_gmresr(1, workres))
      alpha = -1
!
      call daxpy(n, alpha, b, 1, workspace % tmp_gmresr(1, workres), 1)
      call dscal(n, alpha, workspace % tmp_gmresr(1, workres), 1)

!     Calculate residual norm and quit if it is zero
      res0 = dnrm2(n, b, 1)
      resnor = res0
      resid = 0

      if ( res0 .eq. 0.0d0 ) then
         info = 0
         niter = 0
         return  
      end if

      resid=resnor/res0

      if ( resid .le. eps ) then
         info = 0
         niter = 0
         return
      end if
 
!     Main iterative loop ============================
      k = -1
      do while (.true.)

!        Loop to increment dimension of projection subspace
         k=k+1

!        Number of step (not restart) to be done
         its = its + 1
!        write(*,'(A,i3)') '+++++++++++++++++ its ',its 

!        - - - - - - - - - - - - - - - - - - - - - - - - -
!        This part should deliver 
!        u(1,k) <-- invA * resid
!        where u(1,k) is the k-th column of array u(1:n,0:m) and
!        invA is some reasonable approximation to the inverse of A
!
!        If gmres_dim=0 then no inner iterations are performed 
!        to get invA, so that invA is just the identity matrix. 
!        In this case algorithm GMRESR is nothing but GCR
!
!        Otherwise for inner iterations we perform ONE restart of GMRES
!        ATTN: for this implementation it is crucial to perform only
!        ONE restart of GMRES
         if (params % gmresr_dim .eq. 0) then
!           u(1,k) := resid  
            call dcopy(n, workspace % tmp_gmresr(1, workres), 1, &
                & workspace % tmp_gmresr(1, u+mod(k,j)), 1)
            call matvec(params, constants, workspace, &
                & workspace % tmp_gmresr(1, u+mod(k,j)), &
                & workspace % tmp_gmresr(1, c+mod(k,j)))
            nits=nits+1
         else
!           Solve linear system A*u(1,k)=resid by GMRES
!           The stopping criterion for inner iterations is 
!           always absolute but it is automatically adjusted
!           not to be stricter than the stopping criterion for the 
!           outer iterations.  For example, if stop.criterion for
!           the outer iterations is relative than absolute stop.
!           criterion for the inner iterations is (eps*res0)
!           Accuracy for inner iteration:

               epsinn = eps*res0

!           After envoking gmres0 epsinn and itsinn contain actual achieved
!           accuracy and number of performed iterations respectively

            itsinn=params % gmresr_dim

            call gmres0(params, constants, workspace, n, &
                       & workspace % tmp_gmresr(1, workres), &
                       & workspace % tmp_gmresr(1, u+mod(k,j)), &
                       & workspace % tmp_gmresr(1, c+mod(k,j)), &
                       & workspace % tmp_gmresr(1, workgmre), &
                       & epsinn,itsinn,matvec)
            !write(*,*) 'GMRES ', work(1,c+mod(k,j))

            nits=nits+itsinn
         end if           
! - - - - - - - - - - - - - - - - - - - - - - - - 
      
!        Inner loop to orthogonalize 
!        c(1,k) with respect to c(1,k-j),...,c(1,k-1)
!        and to update correspondingly 
!        u(1,k) with respect to u(1,k-j),...,u(1,k-1)
!        parameter j is used only here
         do i = max0(0,k-j),k-1
            alphai = ddot(n, workspace % tmp_gmresr(1, c+mod(i,j)), 1, &
                & workspace % tmp_gmresr(1, c+mod(k,j)), 1)
            call daxpy(n, -alphai, workspace % tmp_gmresr(1, c+mod(i,j)), 1, &
                & workspace % tmp_gmresr(1, c+mod(k,j)), 1)
            call daxpy(n, -alphai, workspace % tmp_gmresr(1, u+mod(i,j)), 1, &
                & workspace % tmp_gmresr(1, u+mod(k,j)), 1)
         end do

!        Normalize c(1,k) and "normalize" u(1,k)
         cknorm = dnrm2(n, workspace % tmp_gmresr(1, c+mod(k,j)), 1)
         cknorm = 1 / cknorm
         call dscal(n, cknorm, workspace % tmp_gmresr(1, c+mod(k,j)), 1)
         call dscal(n, cknorm, workspace % tmp_gmresr(1, u+mod(k,j)), 1)

!        Update current solution and residual
         ckres = ddot(n, workspace % tmp_gmresr(1, c+mod(k,j)), 1, &
             & workspace % tmp_gmresr(1, workres), 1)
         call daxpy(n, ckres, workspace % tmp_gmresr(1, u+mod(k,j)), 1, x, 1)
         call daxpy(n, -ckres, workspace % tmp_gmresr(1, c+mod(k,j)), 1, &
             & workspace % tmp_gmresr(1, workres), 1)
         

!        call show(n,10,x,'GMRESR       ')  

!        Calculate residual norm, check convergence
         resnor = dnrm2(n, workspace % tmp_gmresr(1, workres), 1)

            resid=resnor/res0

         if ( resid .le. eps ) then
            info = 0
            niter = nits
            return
         end if
         if (its .ge. niter*j) then
            info = 1
            niter = nits
            return
         end if

!        print 11, '            ||res|| = ',resnor 
! 11     format(A,d)
! 13     format(i4,A,d)
      
      end do
! End of inifinite iterative loop =================
! End of GMRESR subroutine      
      end 
! This is the modified GMRES routine gmres0 adapted for GMRESR by 
! Mike Botchev, Utrecht University, Dec. 1996
! For detail on how to make GMRES (for GMRESR) cheaper see 
! the above-mentioned paper on GMRESR 
!*************************************************************
! This code was initially written by Youcef Saad
! then revised by Henk A. van der Vorst  
! and Mike Botchev (oct. 1996)
! ************************************************************ 
! gmres algorithm . simple version .  (may 23, 1985)
! parameter list:
! oktest == TRUE for printing intermediate results
! n      == size of problem
! im     == size of krylov subspace:  should not exceed 50 in this
!          version (can be reset in code. looking at comment below)
! rhs    == right hand side
! uu     == initial guess for vector u (see above-mentioned paper on GMRESR)
!           on input, approximate solution on output
! cc     == initial guess for vector c (see above-mentioned paper on GMRESR)
!           on input, approximate solution on output
! work0  == work space of size n x (im+1)
! eps    == tolerance for stopping criterion. process is stopped
!           as soon as ( ||.|| is the euclidean norm):
!           || current residual || <= eps  
! maxits == maximum number of iterations allowed
!           on OUTPUT: actual number of iterations performed
! ----------------------------------------------------------------
! subroutines 
! matvec      == matrix vector multiplication y <- A*x

! BLAS:
! dcopy       == y <-- x routine
! ddot        == dot product function
! dnrm2       == euclidean norm function
! daxpy       == y <-- y+ax routine
! dscal       == x <-- ax routine
! dtsrv       == to solve linear system with a triangular matrix
!*************************************************************
!-------------------------------------------------------------
! arnoldi size should not exceed 10 in this version..
! to reset modify maxdim. BUT:             ----------------
! maxdim was set to 10 because of two reasons:
! (1) it is assumed in this implementation that it is cheaper to
! make maxdim vector updates than to make 1 matrix-vector
! multiplication;
! (2) for large maxdim we may lose the orthogonality property
! on which this cheap implementation is based.
! Please keep it in mind changing maxdim
!-------------------------------------------------------------
!=============================================================================
subroutine gmres0(params, constants, workspace, n, rhs, uu, cc, work0, eps, niter, matvec)

    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    !! Temporary buffers
    type(ddx_workspace_type), intent(inout) :: workspace
      integer maxdim,maxd1,md1max
      parameter (maxdim=20, maxd1=maxdim+1, md1max=maxdim*maxd1)
    procedure(matvec_interface) :: matvec

      integer jjj,jj1
      integer i,i1,its,j,k,k1,niter,n
      double precision cc(n),coeff,coef1,dabs,ddot,dnrm2,dsqrt,eps,epsmac, &
                      & gam,rhs(n),ro,uu(n),work0(n,params % gmresr_dim+1),t     

      double precision hh(maxd1,maxdim),hh1(maxd1,maxdim),c(maxdim), &
                      & s(maxdim),rs(maxd1),rs1(maxd1)

      data (( hh(jj1,jjj), jj1=1,maxd1), jjj=1,maxdim) / md1max*0.0 / ,&
                      & epsmac / 1.d-16 / 
      its = 0

!     ----------------------------
!     Outer loop starts here.. 
!     BUT here (for GMRESR) only 1 outer loop is allowed
!     Compute initial residual vector 
!     ----------------------------
 10   continue
!        do not calculate initial residual first restart because 
!        initial guess is always zero. 
!        make initial guess zero:
         coeff = 0.0
         call dscal(n,coeff,uu,1)
!        make initial residual right-hand side:
         call dcopy(n,rhs,1,work0,1)

         ro = dnrm2 (n, work0, 1)
         if ((ro .eq. 0.0d0).or.(ro .le. eps)) then
            call matvec(params, constants, workspace, uu, cc)
            eps = ro
            niter = its 
            return
         end if

         coeff = 1 / ro
         call dscal(n,coeff,work0,1)

!        initialize 1-st term  of rhs of hessenberg system..
         rs(1) = ro
         i = 0

 4       continue
            i=i+1
            its = its + 1
            i1 = i + 1
            call  matvec(params, constants, workspace, work0(1,i), work0(1,i1))
!           -----------------------------------------
!           modified gram - schmidt...
!           -----------------------------------------
            do j=1, i
               t = ddot(n, work0(1,j),1,work0(1,i1),1)
               hh(j,i) = t
               call daxpy(n, -t, work0(1,j), 1, work0(1,i1), 1)
            end do
            t = dnrm2(n, work0(1,i1), 1)
            hh(i1,i) = t
            if (t .ne. 0.0d0)then
               t = 1 / t
               call dscal(n, t, work0(1,i1), 1)
!              save new column of hh in hh1 to reproduce vector cc later on
               call dcopy(maxd1,hh(1,i),1,hh1(1,i),1)
            endif
!           done with modified gram schmidt and arnoldi step..

!           now  update factorization of hh
            if (i .ne. 1) then
!              perform previous transformations  on i-th column of h
               do k=2,i
                  k1 = k-1
                  t = hh(k1,i)
                  hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
                  hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
               end do
            endif
            gam = dsqrt(hh(i,i)**2 + hh(i1,i)**2)
            if (gam .eq. 0.0d0) gam = epsmac

!           determine next plane rotation
            c(i) = hh(i,i)/gam
            s(i) = hh(i1,i)/gam
            rs(i1) = -s(i)*rs(i)
            rs(i) =  c(i)*rs(i)

!           determine residual norm and test for convergence-
            hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
            ro = dabs(rs(i1))
         if ((i .lt. params % gmresr_dim) .and. (ro .gt. eps))  goto 4

!        now compute solution. first solve upper triangular system.

!        rs := hh(1:i,1:i) ^-1 * rs

         call dtrsv('U','N','N',i,hh,maxd1,rs,1)   
!        done with back substitution..

!        now form linear combination to get vector uu
         do j=1, i
            t = rs(j)
            call daxpy(n, t, work0(1,j), 1, uu,1)
         end do
!        DO NOT restart outer loop EVEN when necessary (that is crucial
!        for this implementation of GMRESR):  NEVER goto 10 !  
!     if (ro .gt. eps .and. its .lt. niter) goto 10

!     Finally, reproduce vector cc as cc = A*uu = work0*hh1*rs:
!     rs := hh1(1:i1,1:i) * rs
      coeff = 1
      coef1 = 0
      call dgemv('N',i1,i,coeff,hh1,maxd1,rs,1,coef1,rs1,1)

!     now multiply Krylov basis vectors work0 by rs:
!     cc := work0*rs
      call dscal(n,coef1,cc,1)
      do j=1, i1
         t = rs1(j)
         call daxpy(n, t, work0(1,j), 1, cc,1)
      end do        

      niter=its
      eps=ro 
      return
!------------------------------- end of gmres0 ----------------------
      end
!
!     jacobi_diis_solver_external
!
!     WARNING: code duplication is bad, but here it makes a very specific sense.
!
!     in ddLPB, we are solving a linear system with Jacobi DIIS also to apply the
!     preconditioner. As calling a solver recursively is a terrible idea, we
!     create a copy of jacobi_diis that 
!
!     - takes as an input the size of the system to be solved (no hardcoded dependence
!       on nsph)
!     - uses copies of the norm subroutine, so that it can also be called with 
!       a user-given size for arrays
!
subroutine jacobi_diis_external(params, constants, workspace, n, tol, rhs, x, n_iter, &
          & x_rel_diff, matvec, dm1vec, norm_func, info)
      type(ddx_params_type),    intent(in)    :: params
      type(ddx_constants_type), intent(in)    :: constants
      type(ddx_workspace_type), intent(inout) :: workspace
      ! Inputs
      integer,                  intent(in)    :: n
      real(dp),                 intent(in)    :: tol
      real(dp),  dimension(n),  intent(in)    :: rhs
      ! Outputs
      real(dp),  dimension(n),  intent(inout) :: x
      integer,                  intent(inout) :: n_iter, info
      real(dp), intent(out) :: x_rel_diff(n_iter)
!
      external                                :: matvec, dm1vec
      procedure(norm_interface)               :: norm_func
      ! Local variables
      integer  :: it, nmat, istatus, lenb, nsph_u
      real(dp) :: diff, norm, rel_diff = zero
      logical  :: dodiis
!
      real(dp), allocatable :: x_new(:), y(:), x_diis(:,:), e_diis(:,:), bmat(:,:)
!
!---------------------------------------------------------------------------------------------
!
!
!     DIIS extrapolation flag
      dodiis = (params%jacobi_ndiis.gt.0)
!
!     extrapolation required
      if (dodiis) then
!
!       allocate workspaces
        lenb = params % jacobi_ndiis + 1
        allocate( x_diis(n,params % jacobi_ndiis), e_diis(n,params % jacobi_ndiis), bmat(lenb,lenb) , stat=istatus )
        if (istatus .ne. 0) then
          write(*,*) ' jacobi_diis: [1] failed allocation (diis)'
          stop
        endif
!        
!       initialize the number of points for diis to one.
!       note that nmat is updated by diis.
        nmat = 1
!        
      endif
!
!     allocate workspaces
      allocate( x_new(n), y(n) , stat=istatus )
      if (istatus .ne. 0) then
        write(*,*) ' jacobi_diis: [2] failed allocation (scratch)' 
        stop
      endif
!
!     Jacobi iterations
!     =================e
      do it = 1, n_iter
!
!       y = rhs - O x
        call matvec(params, constants, workspace, x, y )
        y = rhs - y
!
!       x_new = D^-1 y
        call dm1vec(params, constants, workspace, y, x_new)
!
!       DIIS extrapolation
!       ==================
        if (dodiis) then
!
          x_diis(:,nmat) = x_new
          e_diis(:,nmat) = x_new - x
!
          call diis(n,nmat,params % jacobi_ndiis,x_diis,e_diis,bmat,x_new)
!
        endif
!
!       increment
        x = x_new - x
!
!       warning: in the following, we use a quick and dirty hack to pass 
!       to the norm function arrays double the size than in non ddlpb calculations
!       by defining a fake number of spheres with is n/nbasis (i.e., 2*nsph in 
!       ddlpb).
!
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
!
!       update
        x = x_new
!
!       EXIT Jacobi loop here
!       =====================
        if (rel_diff .le. tol) then
            info = 0
            n_iter = it
            return
        end if
!
      enddo
!
!     something went wrong...
      info = 1
!
      return
!
!
endsubroutine jacobi_diis_external

end module ddx_solvers
