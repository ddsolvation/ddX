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
use ddx_core
implicit none

contains

!---------------------------------------------------------------------------------------------
! Purpose : Jacobi/DIIS solver
!
! variables:
!
!   n        : integer, input, size of the matrix
!
!   lprint   : integer, input, printing flag.
!
!   diis_max : integer, input, number of points to be used for diis extrapolation
!
!              if diis_max = 0, this is just a Jacobi solver.
!
!   norm     : integer, input, norm to be used to evaluate convergence
!              1: max |x_new - x|
!              2: rms (x_new - x)
!              3: rms (x_new - x) and max |x_new - x|
!              4: norm computed by the user-provided function u_norm(n,x)
!
!   tol      : real, input, convergence criterion. if norm = 3, convergence is 
!              achieved when rms (x_new - x) < tol and max |x_new - x| < 10*tol.
!
!   rhs      : real, dimension(n), input, right-hand side of the linear system
!
!   x        : real, dimension(n). In input, a guess of the solution (can be zero).
!              In output, the solution
!
!   n_iter   : integer, in input, the maximum number of iterations. In output,
!              the number of iterations needed to converge.
!
!   ok       : logical, output, T if the solver converged, false otherwise.
!
!   matvec   : external, subroutine to compute the required matrix-vector multiplication
!              format: subroutine matvec(n,x,y)
!
!   dm1vec   : external, subroutine to apply the inverse diagonal matrix to a vector.
!              format: subroutine dm1vec(n,x,y)
! 
!   u_norm   : external, optional function to compute the norm of a vector.
!              format: real(dp) function u_norm(n,x)
!
!---------------------------------------------------------------------------------------------
subroutine jacobi_diis(ddx_data, n, lprint, diis_max, norm, tol, rhs, x, n_iter, &
        & ok, matvec, dm1vec, u_norm)
    ! Inputs
      type(ddx_type), intent(in)       :: ddx_data
      integer,               intent(in)    :: n, diis_max, norm, lprint
      real(dp),                intent(in)    :: tol
      real(dp),  dimension(n), intent(in)    :: rhs
      ! Outputs
      real(dp),  dimension(n), intent(inout) :: x
      integer,               intent(inout) :: n_iter
      logical,               intent(inout) :: ok
      external                             :: matvec, dm1vec
      real(dp),  optional                    :: u_norm
      external                             :: u_norm
      ! Local variables
      integer :: it, nmat, istatus, lenb
      real(dp)  :: rms_norm, max_norm, tol_max, rms_norm_diff, max_norm_diff
      logical :: dodiis
!
      real(dp), allocatable :: x_new(:), y(:), x_diis(:,:), e_diis(:,:), bmat(:,:)
!
!---------------------------------------------------------------------------------------------
!
!     check inputs
      if ( (norm.eq.4) .and. (.not.present(u_norm)) ) then
        write(6,*) ' must provide a function norm(n,x) to evaluate the norm of the increment'
        stop
      endif
!
!     DIIS extrapolation flag
      dodiis =  (diis_max.ne.0)
!
!     set tolerance
      tol_max = 10.0d0 * tol
!
!     extrapolation required
      if (dodiis) then
!
!       allocate workspaces
        lenb = diis_max + 1
        allocate( x_diis(n,diis_max), e_diis(n,diis_max), bmat(lenb,lenb) , stat=istatus )
        if (istatus .ne. 0) then
          write(*,*) ' jacobi_diis: [1] failed allocation (diis)'
          stop
        endif
!        
!       an enigmatic constant
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
!     =================
      do it = 1, n_iter
!
!       y = rhs - O x
        call matvec(ddx_data, x, y )
        !call prtsph("matvec", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, y)
        y = rhs - y
        !call prtsph("matvec y", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, y)
!
!       x_new = D^-1 y
        call dm1vec(ddx_data, y, x_new)
        !call prtsph("matvec D^-1 y", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, x_new)
!
!       DIIS extrapolation
!       ==================
        if (dodiis) then
!
          x_diis(:,nmat) = x_new
          e_diis(:,nmat) = x_new - x
!
          call diis(n,nmat,diis_max,x_diis,e_diis,bmat,x_new)
!
        endif
!
!       increment
        x = x_new - x
!
!       rms/max norm of increment
        if ( norm.le.3 ) then
!
!         compute norm
          call rmsvec( n, x, rms_norm_diff, max_norm_diff )
          call rmsvec( n, x_new, rms_norm, max_norm)
!
!         check norm
          if ( norm.eq.1 ) then
!                  
            ok = (rms_norm_diff .lt. tol*rms_norm)
            
          elseif ( norm.eq.2 ) then
!                  
            ok = (max_norm_diff .lt. tol*max_norm)
!            
          else 

            ok = (rms_norm_diff .lt. tol*rms_norm) .and. (max_norm_diff .lt. tol*max_norm)
!            
          endif
!
!       user-provided norm of increment
        elseif ( norm.eq.4 ) then
!
!         just a placeholder for printing
          max_norm_diff = -1.d0
!
!         compute norm
          rms_norm_diff = u_norm(ddx_data, x )
          rms_norm = u_norm(ddx_data, x_new )
!
!         check norm
          ok = (rms_norm_diff .lt. tol*rms_norm)
!          
        endif
!
!       printing
        if ( lprint.gt.0 ) then
           if (norm.eq.1) then
             write(*,110) it, 'max', max_norm_diff/max_norm
           else if (norm.eq.2) then
             write(*,110) it, 'rms', rms_norm_diff/rms_norm
           else if (norm.eq.3) then
             write(*,100) it, rms_norm_diff/rms_norm, max_norm_diff/max_norm
           else if (norm.eq.4) then
             write(*,120) it, rms_norm_diff/rms_norm
           end if
         end if
  100   format(t3,'iter=',i4,' residual norm (rms,max): ', 2d14.4 )
  110   format(t3,'iter=',i4,' residual norm (',a,'): ', d14.4 )
  120   format(t3,'iter=',i4,' residual norm: ', d20.14 )
!
!       update
        x = x_new
!
!       EXIT Jacobi loop here
!       =====================
        if (ok) exit
!
      enddo
!
!     record number of Jacobi iterations
      n_iter = it
!
      return
!
!
endsubroutine jacobi_diis
!---------------------------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------------------------
!
subroutine diis(n,nmat,ndiis,x,e,b,xnew)
implicit none
integer,                             intent(in)    :: n, ndiis
integer,                             intent(inout) :: nmat
real(dp),  dimension(n,ndiis),         intent(inout) :: x, e
real(dp),  dimension(ndiis+1,ndiis+1), intent(inout) :: b
real(dp),  dimension(n),               intent(inout) :: xnew
!
integer :: nmat1, i, istatus
integer :: j, k
logical :: ok
!
real(dp), allocatable :: bloc(:,:), cex(:)
!
real(dp), parameter :: zero = 0.0d0, one = 1.0d0
!
!------------------------------------------------------------------------------
!
if (nmat.ge.ndiis) then
  do j = 2, nmat - 10
    do k = 2, nmat - 10
      b(j,k) = b(j+10,k+10)
    end do
  end do
  do j = 1, nmat - 10
    x(:,j) = x(:,j+10)
    e(:,j) = e(:,j+10)
  end do
  nmat = nmat - 10
end if
nmat1 = nmat + 1
allocate (bloc(nmat1,nmat1),cex(nmat1) , stat=istatus)
if ( istatus.ne.0 ) then 
  write(*,*) 'diis: allocation failed!'
  stop
endif

call makeb(n,nmat,ndiis,e,b)
bloc   = b(1:nmat1,1:nmat1)
cex    = zero
cex(1) = one
call gjinv(nmat1,1,bloc,cex,ok)
if (.not. ok) then
  nmat = 1
  return
end if
xnew = zero
do i = 1, nmat
  xnew = xnew + cex(i+1)*x(:,i)
end do
nmat = nmat + 1
deallocate (bloc,cex , stat=istatus)
if ( istatus.ne.0 ) then 
  write(*,*) 'diis: deallocation failed!'
  stop
endif
!
return
end subroutine diis
  !
subroutine makeb(n,nmat,ndiis,e,b)
implicit none
integer, intent(in) :: n, nmat, ndiis
real(dp), dimension(n,ndiis),         intent(in) :: e
real(dp), dimension(ndiis+1,ndiis+1), intent(inout) :: b
!
integer :: i
real(dp)  :: bij
real(dp), parameter :: zero = 0.0d0, one = 1.0d0
  
! 1st built
if (nmat.eq.1) then
!
!       [ 0 |  1  ]
!   b = [ --+---- ]
!       [ 1 | e*e ]
!
  b(1,1) = zero
  b(1,2) = one
  b(2,1) = one
  b(2,2) = dot_product(e(:,1),e(:,1))
!
! subsequent builts
else
!
!   first, update the lagrangian line:
  b(nmat+1,1) = one
  b(1,nmat+1) = one
!
!   now, compute the new matrix elements:
  do i = 1, nmat - 1
    bij = dot_product(e(:,i),e(:,nmat))
    b(nmat+1,i+1) = bij
    b(i+1,nmat+1) = bij
  end do
  b(nmat+1,nmat+1) = dot_product(e(:,nmat),e(:,nmat))
end if
!
return
end subroutine makeb
!
subroutine gjinv(n,nrhs,a,b,ok)
implicit none
!
integer,                    intent(in)    :: n, nrhs
logical,                    intent(inout) :: ok
real(dp),  dimension(n,n),    intent(inout) :: a
real(dp),  dimension(n,nrhs), intent(inout) :: b
!
integer :: i, j, k, irow, icol, istatus
real(dp)  :: big, dum, pinv
real(dp), parameter :: zero = 0.0d0, one = 1.0d0
!
integer, allocatable :: indxc(:), indxr(:), piv(:)
real(dp),  allocatable :: scr(:)
!
allocate (indxc(n), indxr(n), piv(n) , stat=istatus)
if ( istatus.ne.0 ) then
  write(*,*)'gjinv: allocation failed! [1]'
  stop
endif
allocate (scr(n) , stat=istatus)
if ( istatus.ne.0 ) then
  write(*,*)'gjinv: allocation failed! [2]'
  stop
endif
!
ok  = .false.
piv = 0
!
irow = 0
icol = 0
do i = 1, n
  big = zero
  do j = 1, n
    if (piv(j).ne.1) then
      do k = 1, n
        if (piv(k).eq.0) then
          if (abs(a(j,k)).gt.big) then
            big  = abs(a(j,k))
            irow = j
            icol = k
          end if
        end if
      end do
    end if
  end do
!
  piv(icol) = piv(icol) + 1
  if (piv(icol) .gt. 1) then
    write(*,1000)
    return
  end if
  if (irow.ne.icol) then
    scr         = a(irow,:)
    a(irow,:)   = a(icol,:)
    a(icol,:)   = scr  
    scr(1:nrhs) = b(irow,:)
    b(irow,:)   = b(icol,:)
    b(icol,:)   = scr(1:nrhs)       
  end if
!
  indxr(i) = irow
  indxc(i) = icol
!
  if (a(icol,icol) .eq. zero) then
    write(*,1000)
    return
  end if
!
  pinv = one/a(icol,icol)
  a(icol,icol) = one
  a(icol,:) = a(icol,:)*pinv
  b(icol,:) = b(icol,:)*pinv
!
  do j = 1, n
    if (j.ne.icol) then
      dum       = a(j,icol)
      a(j,icol) = zero
      a(j,:)    = a(j,:) - a(icol,:)*dum
      b(j,:)    = b(j,:) - b(icol,:)*dum
    end if
  end do
end do
!
do j = n, 1, -1
  if (indxr(j) .ne. indxc(j)) then
    scr           = a(:,indxr(j))
    a(:,indxr(j)) = a(:,indxc(j))
    a(:,indxc(j)) = scr
  end if
end do
!
ok = .true.
deallocate (indxr,indxc,piv,scr , stat=istatus)
if ( istatus.ne.0 ) then
  write(*,*)'gjinv: deallocation failed! [1]'
  stop
endif

return
1000 format (' warning: singular matrix in gjinv!')
end subroutine gjinv
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
!           i.e.  ||res|| / ||res0|| <= eps   OR
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
!                             1 - no convergence within maxits

subroutine gmresr(ddx_data, oktest, n, j, mgmres, b, x, work, eps, stc,&
                & maxits, resid, matvec, iflag)
! ----------------------------------------------------------
! subroutines used
!   matvec   == matrix-vector product y <- A*x

! WARNING: with respect to the original implementation of GMRES, the
!          matvec interface has been changed from matvec(x,y,n) to
!          matvec(n,x,y).

!  blas subroutines:
!   dscal
!   daxpy
!  blas functions:
!   ddot
!   dnrm2 
!**********************************************************

      external matvec

!     list of variables: arrays in alphabetical order,
!     then other variables in alphabetical order

      type(ddx_type), intent(in)  :: ddx_data
      logical oktest
      character*3 stc

      integer i,iflag,its,nits,itsinn,j,k,mgmres,n
      real(dp) maxits

      double precision b(n), x(n), work(n,0 : 2*j+mgmres+2 -1), & 
                       & alpha, alphai, cknorm, ckres, ddot, dnrm2, &
                       & eps, epsinn, res0,resnor,resid

! distribute work space work(n,*) among some virtual variables;
! namely, we think of columns of work as being occupied by 
! c(n,0:j-1), u(n,0:j-1), resid(n), workgmr(n,mgmres+1)
! therefore we define "shifts"

      integer c, u, workres, workgmre
! ----------------------------------------------------------------------------

      if((stc.NE.'rel').and.(stc.NE.'abs'))then
         write(*,*) 'Error in VACGMRESR:'
         write(*,*) 'PARAMETER STC=',stc,' SHOULD BE rel OR abs.'
         STOP
      endif

!     c occupies j columns 0...j-1:
      c = 0 
!     u occupies j columns j...2*j-1:
      u = j
!     resid occupies 1 column No. 2*j:
      workres = 2*j    
!     workgmre occupies mgmres+1 columns 2*j+1...2*j+mgmres+1:
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
      call matvec(ddx_data, n, x,work(1,workres))
      alpha = -1
!
      call daxpy(n,alpha,b,1,work(1,workres),1)
      call dscal(n,alpha,work(1,workres),1)

!     Calculate residual norm and quit if it is zero
      res0 = dnrm2(n,work(1,workres),1)
      resnor = res0
      resid = 0

      if ( res0 .eq. 0.0d0 ) then
         iflag = 0
         maxits = 0
         return  
      end if

      if (stc.eq.'abs') then
         resid=resnor
      else
         resid=resnor/res0
      endif

      if ( resid .le. eps ) then
         iflag = 0
         maxits = 0
         return
      end if 
 
!     Main iterative loop ============================
      k = -1
      do while (.true.)

         if(oktest)write(*,199)its,resid
 199     format('   its =', i4, ' resid =', d20.6)

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
!        If mgmres=0 then no inner iterations are performed 
!        to get invA, so that invA is just the identity matrix. 
!        In this case algorithm GMRESR is nothing but GCR
!
!        Otherwise for inner iterations we perform ONE restart of GMRES
!        ATTN: for this implementation it is crucial to perform only
!        ONE restart of GMRES
         if (mgmres.eq.0) then
!           u(1,k) := resid  
            call dcopy(n,work(1,workres),1,work(1,u+mod(k,j)),1)
            call matvec(ddx_data, n, work(1,u+mod(k,j)), work(1,c+mod(k,j)))
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

            if(stc.eq.'abs')then
               epsinn = eps
            else
               epsinn = eps*res0
            endif

!           After envoking gmres0 epsinn and itsinn contain actual achieved
!           accuracy and number of performed iterations respectively

            itsinn=mgmres

            call gmres0(ddx_data, oktest, n, mgmres, &
                       & work(1,workres),work(1,u+mod(k,j)), &
                       & work(1,c+mod(k,j)),work(1,workgmre), &
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
            alphai = ddot(n,work(1,c+mod(i,j)),1,work(1,c+mod(k,j)),1)
            call daxpy(n, -alphai, work(1,c+mod(i,j)), 1, &
                      & work(1,c+mod(k,j)), 1)
            call daxpy(n, -alphai, work(1,u+mod(i,j)), 1, &
                      & work(1,u+mod(k,j)), 1)
         end do

!        Normalize c(1,k) and "normalize" u(1,k)
         cknorm = dnrm2(n,work(1,c+mod(k,j)),1)
         cknorm = 1 / cknorm
         call dscal(n,cknorm,work(1,c+mod(k,j)),1)
         call dscal(n,cknorm,work(1,u+mod(k,j)),1)

!        Update current solution and residual
         ckres = ddot(n,work(1,c+mod(k,j)),1,work(1,workres),1)
         call daxpy(n, ckres,work(1,u+mod(k,j)),1,x,          1)
         call daxpy(n,-ckres,work(1,c+mod(k,j)),1,work(1,workres),1)
         

!        call show(n,10,x,'GMRESR       ')  

!        Calculate residual norm, check convergence
         resnor = dnrm2(n,work(1,workres),1)

         if (stc.eq.'abs') then
            resid=resnor
         else
            resid=resnor/res0
         endif

         if ( resid .le. eps ) then
            iflag = 0
            maxits = nits
            return
         end if
         if (its .ge. maxits*j) then
            iflag = 1
            maxits = nits
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
subroutine gmres0(ddx_data, oktest, n, im, rhs, uu, cc, work0, eps, maxits, matvec)

      type(ddx_type), intent(in)  :: ddx_data
      integer maxdim,maxd1,md1max
      parameter (maxdim=20, maxd1=maxdim+1, md1max=maxdim*maxd1)
      external matvec

      logical oktest
      integer jjj,jj1
      integer i,i1,im,its,j,k,k1,maxits,n
      double precision cc(n),coeff,coef1,dabs,ddot,dnrm2,dsqrt,eps,epsmac, &
                      & gam,rhs(n),ro,uu(n),work0(n,im+1),t     

      double precision hh(maxd1,maxdim),hh1(maxd1,maxdim),c(maxdim), &
                      & s(maxdim),rs(maxd1),rs1(maxd1)

      data (( hh(jj1,jjj), jj1=1,maxd1), jjj=1,maxdim) / md1max*0.0 / ,&
                      & epsmac / 1.d-16 / 
!-----------------------------------------------------------------------------

      if (im .gt. maxdim) then
         im = maxdim
         write (*,'(A,i2)') 'GMRES0: dimension has been reduced to ',im
         write (*,'(A)') ' => reset MAXDIM if you want it to be more'
         write (*,'(A)') ' BUT read comments near MAXDIM before'
      end if

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
            call matvec(ddx_data, n, uu, cc)
            eps = ro
            maxits = its 
            return
         end if

         coeff = 1 / ro
         call dscal(n,coeff,work0,1)

         if (oktest) write(*, 199) its, ro

!        initialize 1-st term  of rhs of hessenberg system..
         rs(1) = ro
         i = 0

 4       continue
            i=i+1
            its = its + 1
            i1 = i + 1
            call  matvec(ddx_data, n, work0(1,i), work0(1,i1))
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
            if (oktest) write(*, 199) its, ro
         if ((i .lt. im) .and. (ro .gt. eps))  goto 4

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
!     if (ro .gt. eps .and. its .lt. maxits) goto 10

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

 199  format('itsinn =', i4, ' res. norm =', d20.6)

      maxits=its
      eps=ro 
      return
!------------------------------- end of gmres0 ----------------------
      end

end module ddx_solvers
