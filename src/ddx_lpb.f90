!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_core.f90
!! Core routines and parameters of entire ddX software
!!
!! @version 1.0.0
!! @author Abhinav Jha and Michele Nottoli
!! @date 2021-02-25

module ddx_lpb
use ddx_core
use ddx_operators
use ddx_solvers
implicit none
!!
!! Logical variables for iterations, cosmo solver, and HSP solver
!!
logical :: first_out_iter
logical :: do_diag
integer :: matAB
!!
!! Hardcoded values
!!
integer :: nbasis0
integer :: lmax0
real(dp),  parameter :: epsp = 1.0d0
!!
!! Taken from Chaoyu's MATLAB code
!!
!! wij          : (?) Maybe Eq. (54) 
!! coefvec      : Intermediate value in computation of update_rhs
!! Pchi         : Pchi matrix, Eq. (87)
!! Qmat         : (?)
!! coefY        : Intermediate calculation in Q Matrix Eq. (91)
!! C_ik         : (i'_l0(r_j)/i_l0(r_j)-k'_l0(r_j)/k_l0(r_j))^{-1}
real(dp), allocatable :: wij(:,:)
real(dp), allocatable :: coefvec(:,:,:), Pchi(:,:,:), &
                         & Qmat(:,:,:), coefY(:,:,:)
real(dp), allocatable :: C_ik(:,:)
!! SI_ri        : Bessel' function of first kind
!! DI_ri        : Derivative of Bessel' function of first kind
!! SK_ri        : Bessel' function of second kind
!! DK_ri        : Derivative of Bessel' function of second kind
!! termimat     : i'_l(r_j)/i_l(r_j)
!! tol_gmres    : Tolerance of GMRES iteration
!! n_iter_gmres : Maximum number of GMRES itertation

real(dp), allocatable :: SI_ri(:,:), DI_ri(:,:), SK_ri(:,:), &
                              & DK_ri(:,:), termimat(:,:)
real(dp)              :: tol_gmres, n_iter_gmres

!! Terms related to Forces of ddLPB model
real(dp), allocatable :: diff_re_c1(:,:), diff_re_c2(:,:)
real(dp), allocatable :: diff0_c1(:,:), diff0_c2(:,:)
real(dp), allocatable :: diff_ep_c1(:), diff_ep_c2(:)
contains
  !!
  !! ddLPB calculation happens here
  !! @param[in] ddx_data : dd Data 
  !! @param[in] phi      : Boundary conditions
  !! @param[in] psi      : Electrostatic potential vector.
  !!                       Use of psi unknown
  !! @param[in] gradphi  : Gradient of phi
  !! @param[out] esolv   : Electrostatic solvation energy
  !!
  subroutine ddlpb(ddx_data, phi, gradphi, psi, esolv, force)
  ! main ddLPB
  implicit none
  ! Inputs
  type(ddx_type), intent(inout) :: ddx_data
  real(dp), intent(in)       :: phi(ddx_data % ncav), &
                                & gradphi(3,ddx_data % ncav)
  real(dp), intent(in)       :: psi(ddx_data % nbasis, ddx_data % nsph)
  ! Outputs
  real(dp), intent(out)      :: esolv, force(3, ddx_data % nsph)
  logical                    :: converged = .false.
  integer                    :: iteration = 1, istatus
  real(dp)                   :: inc, old_esolv
  !!
  !! Xr         : Reaction potential solution (Laplace equation)
  !! Xe         : Extended potential solution (HSP equation)
  !! rhs_r      : Right hand side corresponding to Laplace equation
  !! rhs_e      : Right hand side corresponding to HSP equation
  !! rhs_r_init : Initial right hand side corresponding to Laplace equation
  !! rhs_e_init : Initial right hand side corresponding to HSP equation
  !!
  real(dp), allocatable ::   Xr(:,:), Xe(:,:), rhs_r(:,:), rhs_e(:,:), &
                                  & rhs_r_init(:,:), rhs_e_init(:,:), &
                                  & Xadj_r(:,:), Xadj_e(:,:), &
                                  & Xadj_r_sgrid(:,:), Xadj_e_sgrid(:,:), &
                                  & force_r(:,:), force_e(:,:)
  !!
  !! g      : Intermediate matrix for computation of g0
  !! f      : Intermediate matrix for computation of f0
  !! g0     : Vector associated to psi_0 Eq.(77) QSM19.SISC
  !! f0     : Vector associated to partial_n_psi_0 Eq.(99) QSM19.SISC
  !! lx     : External routine from matvec.f90. Used for Jacobi solver
  !! ldm1x  : External routine from matvec.f90. Used for Jacobi solver
  !! hnorm  : External routine from matvec.f90. Used for Jacobi solver
  !!          h^-1/2 norm of the increment on each sphere
  !! ok     : Boolean to check convergence of solver
  !! n_iter : Number of iterative steps
  real(dp), allocatable :: g(:,:), f(:,:), g0(:), f0(:), phi_grid(:, :)
  integer                    :: isph
  integer                    :: i
  logical                    :: ok = .false.
  integer                    :: n_iter
  integer                    :: its
  
  ! Local variables, used in force computation
  real(dp), allocatable :: vsin(:), vcos(:), vplm(:), basloc(:), &
        & dbsloc(:, :), fx(:, :), ef(:, :)
  ! lmax0 set to minimum of 6 or given lmax.
  ! nbasis0 set to minimum of 49 or given (lmax+1)^2.
  ! Previous implementation had hard coded value 6 and 49.
  lmax0 = MIN(6, ddx_data % lmax)
  nbasis0 = MIN(49, ddx_data % nbasis)
  
  !
  ! Allocate Bessel's functions of the first kind and the second kind
  ! and their derivatives
  !
  call ddlpb_init(ddx_data)

  allocate(g(ddx_data % ngrid,ddx_data % nsph),&
           & f(ddx_data % ngrid, ddx_data % nsph), &
           & phi_grid(ddx_data % ngrid, ddx_data % nsph), &
           & g0(ddx_data % nbasis),f0(ddx_data % nbasis),&
           & rhs_r(ddx_data % nbasis, ddx_data % nsph), &
           & rhs_e(ddx_data % nbasis, ddx_data % nsph), &
           & rhs_r_init(ddx_data % nbasis, ddx_data % nsph),&
           & rhs_e_init(ddx_data % nbasis, ddx_data % nsph), &
           & Xr(ddx_data % nbasis, ddx_data % nsph),&
           & Xe(ddx_data % nbasis, ddx_data % nsph), &
           & Xadj_r(ddx_data % nbasis, ddx_data % nsph),&
           & Xadj_e(ddx_data % nbasis, ddx_data % nsph), &
           & vsin(ddx_data % lmax+1), vcos(ddx_data % lmax+1), &
           & vplm(ddx_data % nbasis), basloc(ddx_data % nbasis), &
           & dbsloc(3, ddx_data % nbasis), fx(3, ddx_data % nsph), &
           & Xadj_r_sgrid(ddx_data % ngrid, ddx_data % nsph), &
           & Xadj_e_sgrid(ddx_data % ngrid, ddx_data % nsph), &
           & force_r(3, ddx_data % nsph), &
           & force_e(3, ddx_data % nsph), &
           & stat = istatus)
  if (istatus.ne.0) write(6,*) 'ddlpb allocation failed'

  !do i = 1, ncav
  !  write(6,'(3F15.8)') phi(i), gradphi(:,i)
  !end do
  !stop

  ! Build the right hand side
  ! do i = 1, ncav
  !   write(6,'(4F20.10)') phi(i), gradphi(:,i)
  ! end do
  !
  !! wghpot: Weigh potential at cavity points. Comes from ddCOSMO
  !1         Intermediate computation of G_0 Eq.(77) QSM19.SISC
  !!
  !! @param[in]  phi : Boundary conditions (This is psi_0 Eq.(20) QSM19.SISC)
  !! @param[out] g   : Boundary conditions on solute-solvent boundary gamma_j_e
  !!
  call wghpot(ddx_data, phi, phi_grid, g)
  !!
  !! wghpot_f : Intermediate computation of F_0 Eq.(75) from QSM19.SISC
  !!
  !! @param[in]  gradphi : Gradient of psi_0
  !! @param[out] f       : Boundary conditions scaled by characteristic function
  !!
  call wghpot_f(ddx_data, gradphi,f)

  ! do isph = 1, nsph
  !   do i = 1, ngrid
  !     write(6,'(2F20.10)') g(i,isph), f(i,isph)
  !   end do
  ! end do
  
  !!
  !! Integrate Right hand side
  !! rhs_r_init: g0+f0
  !! rhs_e_init: f0
  !!
  do isph = 1, ddx_data % nsph
    !! intrhs is a subroutine in ddx_operators
    !! @param[in]  isph : Sphere number, used for output
    !! @param[in]  g    : Intermediate right side g
    !! @param[out] g0   : Integrated right side Eq.(77) in QSM19.SISC
    call intrhs(ddx_data % iprint, ddx_data % ngrid, &
                ddx_data % lmax, ddx_data % vwgrid, ddx_data % vgrid_nbasis, &
                & isph, g(:,isph), g0)
    call intrhs(ddx_data % iprint, ddx_data % ngrid, &
                ddx_data % lmax, ddx_data % vwgrid, ddx_data % vgrid_nbasis, &
                & isph,f(:,isph),f0)
    !! rhs 
    rhs_r_init(:,isph) = g0 + f0
    rhs_e_init(:,isph) = f0
  end do

  rhs_r = rhs_r_init
  rhs_e = rhs_e_init
  
  n_iter = ddx_data % maxiter
  
  first_out_iter = .true.

  do while (.not.converged)

    !! Solve the ddCOSMO step
    !! A X_r = RHS_r (= G_X+G_0) 
    !! Call Jacobi solver
    !! @param[in]      nsph*nylm : Size of matrix
    !! @param[in]      iprint    : Flag for printing
    !! @param[in]      ndiis     : Number of points to be used for 
    !!                            DIIS extrapolation. Set to 25 in ddCOSMO
    !! @param[in]      4         : Norm to be used to evaluate convergence
    !!                             4 refers to user defined norm. Here hnorm
    !! @param[in]      tol       : Convergence tolerance
    !! @param[in]      rhs_r     : Right-hand side
    !! @param[in, out] xr        : Initial guess to solution and final solution
    !! @param[in, out] n_iter    : Number of iterative steps
    !! @param[in, out] ok        : Boolean to check whether the solver converged
    !! @param[in]      lx        : External subroutine to compute matrix 
    !!                             multiplication, i.e., Lx_r, comes from matvec.f90
    !! @param[in]      ldm1x     : External subroutine to apply invert diagonal
    !!                             matrix to vector, i.e., L^{-1}x_r, comes from matvec.f90
    !! @param[in]      hnorm     : User defined norm, comes from matvec.f90
    call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, ddx_data % ndiis, &
                     & 4, ddx_data % tol, rhs_r, Xr, n_iter, ok, lx, ldm1x, hnorm)
    call convert_ddcosmo(ddx_data, 1, Xr)
    ! call print_ddvector('xr',xr)
  
    !! Solve ddLPB step
    !! B X_e = RHS_e (= F_0)
    call lpb_hsp(ddx_data, rhs_e, Xe)
    ! call print_ddvector('xe',xe)
  
    !! Update the RHS
    !! / RHS_r \ = / g + f \ - / c1 c2 \ / X_r \
    !! \ RHS_e /   \ f     /   \ c1 c2 / \ X_e /
    call update_rhs(ddx_data, rhs_r_init, rhs_e_init, rhs_r, rhs_e, Xr, Xe)
    ! call print_ddvector('rhs_r',rhs_r)
    ! call print_ddvector('rhs_e',rhs_e)

    !! Compute energy
    !! esolv = pt5*sprod(nsph*nylm,xr,psi)
    esolv = zero
    do isph = 1, ddx_data % nsph
      esolv = esolv + pt5*ddx_data % charge(isph)*Xr(1,isph)*(one/(two*sqrt(pi)))
    end do

    !! Check for convergence
    inc = abs(esolv - old_esolv)/abs(esolv)
    old_esolv = esolv
    if ((iteration.gt.1) .and. (inc.lt.ddx_data % tol)) then
      write(6,*) 'Reach tolerance.'
      converged = .true.
    end if
    write(6,*) iteration, esolv, inc
    iteration = iteration + 1

    ! to be removed
    first_out_iter = .false.
  end do
  ! Start the Force computation
  if(ddx_data % force .eq. 1) then
    ! stop "Forces not implemented in ddLPB"
    ! Call the subroutine adjoint to solve the adjoint solution
    call ddx_lpb_adjoint(ddx_data, psi, Xadj_r, Xadj_e)
    call debug_derivative(ddx_data)
    ! Call dgemm to integrate the adjoint solution on the grid points
    call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % vgrid, ddx_data % vgrid_nbasis, &
            & Xadj_r , ddx_data % nbasis, zero, Xadj_r_sgrid, &
            & ddx_data % ngrid)
    call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % vgrid, ddx_data % vgrid_nbasis, &
            & Xadj_e , ddx_data % nbasis, zero, Xadj_e_sgrid, &
            & ddx_data % ngrid)
    ! Compute A^k*Xadj_r, using Subroutine from ddCOSMO
    do isph = 1, ddx_data % nsph
      call fdoka(ddx_data, isph, Xr, Xadj_r_sgrid(:, isph), &
                & basloc, dbsloc, vplm, vcos, vsin, force_r(:,isph))
      call fdokb(ddx_data, isph, Xr, Xadj_r_sgrid, basloc, &
                & dbsloc, vplm, vcos, vsin, force_r(:, isph))
    end do
    ! Remove the factor of 4pi/2l+1
    call convert_force(ddx_data, -1, force_r)
    ! Compute B^k*Xadj_e
    do isph = 1, ddx_data % nsph
      call fdoka_b_xe(ddx_data, isph, Xe, Xadj_e_sgrid(:, isph), &
                & basloc, dbsloc, vplm, vcos, vsin, force_e(:,isph))
      call fdokb_b_xe(ddx_data, isph, Xe, Xadj_e_sgrid, &
                & basloc, dbsloc, vplm, vcos, vsin, force_e(:, isph))
    end do
    ! Add contributions with respect to C_1 and C_2 matrix
    ! It's a seperate subroutine because one wants to use convert_force before
    do isph = 1, ddx_data % nsph
      ! fdoc: Force, derivative of C
      call fdoc(ddx_data, isph, Xr, Xe, Xadj_r, Xadj_e, Xadj_r_sgrid,&
               & Xadj_e_sgrid, force_r(:, isph), force_e(:, isph))
    end do
    if (ddx_data % iprint .ge. 5) then
      write(*,*) 'Forces Laplace:'
      do isph = 1, ddx_data % nsph
        write(6,'(1x,i5,3ES25.16E3)') isph, force_r(:,isph)
      end do
      write(*,*) 'Forces HSP:'
      do isph = 1, ddx_data % nsph
        write(6,'(1x,i5,3ES25.16E3)') isph, force_e(:,isph)
      end do
  end if
    force = force_r + force_e
  endif
  deallocate(g, f, phi_grid, g0 ,f0, &
           & rhs_r, rhs_e, rhs_r_init , &
           & rhs_e_init, Xr, Xe, Xadj_r, Xadj_e, &
           & Xadj_r_sgrid, Xadj_e_sgrid, &
           & force_r, force_e, &
           & vsin, vcos, vplm, basloc,&
           & dbsloc, fx, stat = istatus)
  if (istatus.ne.0) write(6,*) 'ddlpb deallocation failed'

  return
  end subroutine ddlpb
  !!
  !! Allocate Bessel's functions of the first kind and the second kind
  !! Uses the file bessel.f90
  !! @param[out] SI_ri : Bessel's function of the first kind
  !! @param[out] DI_ri : Derivative of Bessel's function of the first kind
  !! @param[out] SK_ri : Bessel's function of the second kind
  !! @param[out] DK_ri : Derivative of Bessel's function of the second kind
  !! @param[out] NM    : Highest order computed
  !!
  subroutine ddlpb_init(ddx_data)
  use bessel
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  integer                         :: istatus, isph, NM
  allocate(SI_ri(0:ddx_data % lmax, ddx_data % nsph),&
           & DI_ri(0:ddx_data % lmax, ddx_data % nsph),&
           & SK_ri(0:ddx_data % lmax, ddx_data % nsph), &
           & DK_ri(0:ddx_data % lmax, ddx_data % nsph), &
           & diff_re_c1(ddx_data % nbasis, ddx_data % nsph), &
           & diff_re_c2(ddx_data % nbasis, ddx_data % nsph), &
           & diff0_c1(nbasis0, ddx_data % nsph), &
           & diff0_c2(nbasis0, ddx_data % nsph), &
           & diff_ep_c1(ddx_data % ncav), &
           & diff_ep_c2(ddx_data % ncav), &
           & termimat(0:ddx_data % lmax, ddx_data % nsph), stat=istatus)
  if (istatus.ne.0) then
    write(*,*)'ddinit : [1] allocation failed !'
    stop
  end if
  do isph = 1, ddx_data % nsph
    call SPHI_bessel(ddx_data % lmax,ddx_data % rsph(isph)*ddx_data % kappa,&
                     & NM,SI_ri(:,isph), &
                     & DI_ri(:,isph))
    call SPHK_bessel(ddx_data % lmax,ddx_data % rsph(isph)*ddx_data % kappa,&
                     & NM,SK_ri(:,isph), &
                     & DK_ri(:,isph))
  end do
  return
  end subroutine ddlpb_init

  !!
  !! Find intermediate F0 in the RHS of the ddLPB model given in Eq.(82)
  !! @param[in]  gradphi : Gradient of psi_0
  !! @param[out] f       : Intermediate calculation of F0
  !!
  subroutine wghpot_f(ddx_data, gradphi, f )
  use bessel
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  real(dp), dimension(3, ddx_data % ncav),       intent(in)  :: gradphi
  real(dp), dimension(ddx_data % ngrid, ddx_data % nsph),    intent(out) :: f

  integer :: isph, ig, ic, ind, ind0, jg, l, m, jsph
  real(dp) :: nderphi, sumSijn, rijn, coef_Ylm, sumSijn_pre, termi, &
      & termk, term
  real(dp), dimension(3) :: sijn, vij
  real(dp) :: rho, ctheta, stheta, cphi, sphi
  real(dp), allocatable :: SK_rijn(:), DK_rijn(:)

  integer :: l0, m0, NM, kep, istatus
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph) :: c0
  real(dp), dimension(0:ddx_data % lmax, ddx_data % nsph) :: coef_bessel
  real(dp), allocatable :: vplm(:), basloc(:), vcos(:), vsin(:)

  ! initialize
  allocate(vplm(ddx_data % nbasis),basloc(ddx_data % nbasis),vcos(ddx_data % lmax+1),vsin(ddx_data % lmax+1))
  allocate(SK_rijn(0:lmax0),DK_rijn(0:lmax0))
  ic = 0 ; f(:,:)=0.d0
  !
  ! Compute c0 Eq.(98) QSM19.SISC
  !
  do isph = 1, ddx_data % nsph
    do ig = 1, ddx_data % ngrid
      if ( ddx_data % ui(ig,isph).ne.zero ) then
        ic = ic + 1
        nderphi = dot_product( gradphi(:,ic),ddx_data % cgrid(:,ig) )
        c0(:, isph) = c0(:,isph) + ddx_data % wgrid(ig)*ddx_data % ui(ig,isph)*nderphi*ddx_data % vgrid(:,ig)
      end if
    end do
  end do

  allocate (coefY(ddx_data % ncav, nbasis0, ddx_data % nsph), & 
            & C_ik(0:ddx_data % lmax, ddx_data % nsph), stat = istatus)
  ! memuse = memuse + ncav*nbasis0*nsph
  ! memmax = max(memmax,memuse)
  if ( istatus .ne. 0 ) then
      write(*,*)'wghpot_f : [1] allocation failed!'
      stop
  end if
  !
  ! Compute coef_bessel
  ! Here (der_i_l0/i_l0 - der_k_l0/k_l0)^(-1) is computed in Eq.(97)
  ! Here we consider three cases. One is when kappa is too small,
  ! kappa is too big, and kappa is moderate. To find these values,
  ! we use the formulas given in the book, Mathematical Methods for 
  ! Physicsts, Arfken and Weber.
  !
  do jsph = 1, ddx_data % nsph
    do l0 = 0, lmax0
      ! Case: kappa > tol_inf
      if (max(DI_ri(l0,jsph), SI_ri(l0,jsph)).gt.tol_inf) then
        termi = ddx_data % kappa
      ! Case: kappa < tol_zero
      !       One can ignore the other terms (including the second
      !       term given below) as kappa is so small the
      !       contributions are negligible
      else if (min(DI_ri(l0,jsph), SI_ri(l0,jsph)).lt.tol_zero) then
        termi = l0/ddx_data % rsph(jsph) + &
            & (l0 + one)*(ddx_data % kappa**2*ddx_data % rsph(jsph))/((two*l0 + one) * &
            & (two*l0 + three))
      ! Case: kappa is of normal size.
      ! NOTE: We notice a factor of kappa. The reason being while computing
      !       DI_ri the argument is kappa*r. Hence a factor of kappa.
      else
        termi = DI_ri(l0,jsph)/SI_ri(l0,jsph)*ddx_data % kappa
      end if
      !write(*,*) SI_ri(l0,jsph), termi

      ! Similar calculation for SK_ri to SI_ri
      if (SK_ri(l0,jsph).gt.tol_inf) then
        termk = - (l0 + one)/ddx_data % rsph(jsph) - &
            & l0*(ddx_data % kappa**2*ddx_data % rsph(jsph))/((two*l0 - one)*(two*l0 + one))
      else if (SK_ri(l0,jsph).lt.tol_zero) then
        termk = -ddx_data % kappa
      else
        termk = DK_ri(l0,jsph)/SK_ri(l0,jsph)*ddx_data % kappa
      end if

      !write(*,*) SK_ri(l0,jsph), termk
      coef_bessel(l0,jsph) = one/(termi - termk)
      C_ik(l0, jsph) = one/(termi - termk)
      !write(*,*) DI_ri(l0,jsph), SI_ri(l0,jsph), coef_bessel(l0,jsph)
      !write(*,*) (min(-DK_ri(l0,jsph), SK_ri(l0,jsph)).lt.tol_zero), &
      !    & DK_ri(l0,jsph), termk
    end do
  end do

  ! Computation of F0 using above terms
  ! kep: External grid poitns
  kep = 0
  do isph = 1, ddx_data % nsph
    do ig = 1, ddx_data % ngrid
      if (ddx_data % ui(ig,isph).gt.zero) then
        kep = kep + 1
        sumSijn = zero
        ! Loop to compute Sijn
        do jsph = 1, ddx_data % nsph
          ! (?) Use of sumSijn_pre (?)
          sumSijn_pre = sumSijn
          vij  = ddx_data % csph(:,isph) + ddx_data % rsph(isph)*ddx_data % cgrid(:,ig) - ddx_data % csph(:,jsph)
          rijn = sqrt(dot_product(vij,vij))
          sijn = vij/rijn
          
          ! Compute Bessel function of 2nd kind for the coordinates
          ! (s_ijn, r_ijn) and compute the basis function for s_ijn
          call SPHK_bessel(lmax0,rijn*ddx_data % kappa,NM,SK_rijn,DK_rijn)
          call ylmbas(sijn , rho, ctheta, stheta, cphi, &
                      & sphi, ddx_data % lmax, ddx_data % vscales, &
                      & basloc, vplm, vcos, vsin)

          do l0 = 0,lmax0
            ! term: k_l0(r_ijn)/k_l0(r_i)
            ! Uses approximation of SK_ri in double factorial terms when argument
            ! is less than one
            if (SK_ri(l0,jsph).gt.tol_inf) then
            ! Uses approximation of SK_ri in double factorial terms when argument
            ! is greater than l0
              term = (ddx_data % rsph(jsph)/rijn)**(l0+1)
            else if (SK_ri(l0,jsph).lt.tol_zero) then
              term = (ddx_data % rsph(jsph)/rijn)*exp(-ddx_data % kappa*(rijn-ddx_data % rsph(jsph)))
            else
              term = SK_rijn(l0)/SK_ri(l0,jsph)
            end if
            ! coef_Ylm : (der_i_l0/i_l0 - der_k_l0/k_l0)^(-1)*k_l0(r_ijn)/k_l0(r_i)
            coef_Ylm =  coef_bessel(l0,jsph)*term
            do m0 = -l0, l0
              ind0 = l0**2 + l0 + m0 + 1
              sumSijn = sumSijn + c0(ind0,jsph)*coef_Ylm*basloc(ind0)
              ! coefY : Intermediate calculations in Q matrix
              coefY(kep,ind0,jsph) = coef_Ylm*basloc(ind0)
            end do
          end do
        end do
        !
        ! Here Intermediate value of F_0 is computed Eq. (99)
        ! Mutilplication with Y_lm and weights will happen afterwards
        !write(6,*) sumSijn, epsp, eps, ddx_data % ui(ig,isph)
        f(ig,isph) = -(epsp/ddx_data % eps)*ddx_data % ui(ig,isph) * sumSijn
      end if
    end do
  end do 

  deallocate( vplm, basloc, vcos, vsin, SK_rijn, DK_rijn  )
  return
  end subroutine wghpot_f

  !
  ! Subroutine used for the GMRES solver
  ! NOTE: It is refered as matABx in the GMRES solver.
  !       Fortran is not case sensitive
  ! @param[in]      n : Size of the matrix
  ! @param[in]      x : Input vector
  ! @param[in, out] y : y=A*x
  !
  subroutine matABx(ddx_data, n, x, y )
  implicit none 
  type(ddx_type), intent(in)  :: ddx_data
  integer, intent(in) :: n
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: x
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(inout) :: y
  integer :: isph, istatus
  real(dp), allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)
  integer :: i
  ! allocate workspaces
  allocate( pot(ddx_data % ngrid), vplm(ddx_data % nbasis), basloc(ddx_data % nbasis), &
            & vcos(ddx_data % lmax+1), vsin(ddx_data % lmax+1), stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'Bx: allocation failed !'
    stop
  endif
  
  if (ddx_data % iprint .ge. 5) then
      call prtsph('X', ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, &
          & x)
  end if

  y = zero
  do isph = 1, ddx_data % nsph
    call calcv2_lpb(ddx_data, isph, pot, x, basloc, vplm, vcos, vsin )
    ! intrhs comes from ddCOSMO
    call intrhs(ddx_data % iprint, ddx_data % ngrid, &
                ddx_data % lmax, ddx_data % vwgrid, ddx_data % vgrid_nbasis, &
                & isph, pot, y(:,isph) )
    ! Action of off-diagonal blocks
    y(:,isph) = - y(:,isph)
    ! Add action of diagonal block
    y(:,isph) = y(:,isph) + x(:,isph)
  end do
  
  if (ddx_data % iprint .ge. 5) then
      call prtsph('Bx (off-diagonal)', ddx_data % nbasis, ddx_data % lmax, &
          & ddx_data % nsph, 0, y)
  end if
  deallocate( pot, basloc, vplm, vcos, vsin , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'matABx: allocation failed !'
    stop
  endif
  end subroutine matABx

  !!
  !! Scale the ddCOSMO solution vector
  !! @param[in]      direction : Direction of the scaling
  !! @param[in, out] vector    : ddCOSMO solution vector
  !!
  subroutine convert_ddcosmo(ddx_data, direction, vector)
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  integer :: isph, l, m, ind
  integer, intent(in) :: direction
  real(dp), intent(inout) :: vector(ddx_data % nbasis, ddx_data % nsph)
  real(dp) :: fac
  
  do isph = 1, ddx_data % nsph
    do l = 0, ddx_data % lmax
      ind = l*l + l + 1
      fac = four*pi/(two*dble(l) + one) 
      if (direction.eq.-1) fac = one/fac
      do m = -l, l
        vector(ind + m,isph) = fac*vector(ind + m,isph)
      end do
    end do
  end do
  return
  end subroutine convert_ddcosmo
  
  !!
  !! Scale the Force solution vector
  !! @param[in]      direction : Direction of the scaling
  !! @param[in, out] vector    : Force vector
  !!
  subroutine convert_force(ddx_data, direction, vector)
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  integer :: isph, l, m
  integer, intent(in) :: direction
  real(dp), intent(inout) :: vector(3, ddx_data % nsph)
  real(dp) :: fac
  
  do isph = 1, ddx_data % nsph
    fac = four*pi/(two*dble(l) + one) 
    if (direction.eq.-1) fac = one/fac
    do m = 1, 3
      vector(m,isph) = fac*vector(m,isph)
    end do
  end do
  return
  end subroutine convert_force
  
  !!
  !! Scale the ddLPB diagonal vector
  !! @param[in, out] vector    : ddLPB vector
  !!
  subroutine convert_ddlpb(ddx_data, vector)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    ! Output
    real(dp), intent(inout) :: vector(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    integer :: l, ind
    ! Loop over harmonics
    do l = 0, ddx_data % lmax
        ind = l*l + l + 1
        vector(ind-l:ind+l, :) = vector(ind-l:ind+l, :) / (ddx_data % vscales(ind)**2)
    end do
  end subroutine convert_ddlpb


  !!
  !! Solve the HSP problem
  !! @param[in]      rhs : Right hand side for HSP
  !! @param[in, out] Xe  : Solution vector
  !!
  subroutine lpb_hsp(ddx_data, rhs, Xe)
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: rhs
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(inout) :: Xe
  integer :: isph, istatus, n_iter, info, c1, c2, cr
  real(dp) :: tol, r_norm
  real(dp), allocatable :: work(:,:)
  integer, parameter  :: gmm = 20, gmj = 25
  
  allocate(work(ddx_data % nsph*ddx_data % nbasis, 0:2*gmj + gmm + 2 - 1),stat=istatus)
  if (istatus.ne.0) then
    write(*,*) ' LPB-HSP: failed allocation for GMRES'
    stop
  endif
   
  work = zero
  Xe = rhs
  
  matAB = 1
  tol_gmres = 1.0d-8
  !!
  !! Call GMRES solver
  !! @param[in]      Residue_print : Set to false by default. Prints the
  !!                                 intermediate residue.
  !! @param[in]      nsph*nylm     : Size of matrix
  !! @param[in]      gmj           : Integer truncation parameter, Default: 25
  !! @param[in]      gmm           : Integer dimension of the GMRES,
  !!                                 Default: 20
  !! @param[in]      rhs           : Right hand side
  !! @param[in, out] Xe            : Initial guess of the problem
  !! @param[in]      work          : Work space of size
  !!                               : nsph*nylm X (2*gmj + gmm + 2)
  !! @param[in]      tol_gmres     : GMRES tolerance
  !! @param[in]      Stopping      : Stopping criteria, Default set to
  !!                                 'rel' for relative. Other option
  !!                                 'abs' for absolute
  !! @param[in]      n_iter_gmres  : Number of GMRES iteration
  !! @param[in]      r_norm        : Residual measure
  !! @param[in]      matABx        : Subroutine A*x. Named matabx in file
  !! @param[in, out] info          : Flag after solve. 0 means within tolerance
  !!                                 1 means max number of iteration
  call gmresr(ddx_data, .false., ddx_data % nsph*ddx_data % nbasis, gmj, gmm, & 
             & rhs, Xe, work, tol_gmres,'rel', n_iter_gmres, r_norm, matABx, info)

  deallocate(work)
  endsubroutine lpb_hsp

  !
  ! Intermediate computation of BX_e
  ! @param[in]      isph   : Number of the sphere
  ! @param[in, out] pot    : Array of size ngrid
  ! @param[in]      x      : Input vector (Usually X_e)
  ! @param[in, out] basloc : Used to compute spherical harmonic
  ! @param[in, out] vplm   : Used to compute spherical harmonic
  ! @param[in, out] vcos   : Used to compute spherical harmonic
  ! @param[in, out] vsin   : Used to compute spherical harmonic
  !
  subroutine calcv2_lpb (ddx_data, isph, pot, x, basloc, vplm, vcos, vsin )
  type(ddx_type), intent(in) :: ddx_data
  integer, intent(in) :: isph
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: x
  real(dp), dimension(ddx_data % ngrid), intent(inout) :: pot
  real(dp), dimension(ddx_data % nbasis), intent(inout) :: basloc
  real(dp), dimension(ddx_data % nbasis), intent(inout) :: vplm
  real(dp), dimension(ddx_data % lmax+1), intent(inout) :: vcos
  real(dp), dimension(ddx_data % lmax+1), intent(inout) :: vsin
  real(dp), dimension(ddx_data % nbasis) :: fac_cosmo, fac_hsp
  integer :: its, ij, jsph
  real(dp) :: rho, ctheta, stheta, cphi, sphi
  real(dp) :: vij(3), sij(3)
  real(dp) :: vvij, tij, xij, oij

  pot = zero
  do its = 1, ddx_data % ngrid
    if (ddx_data % ui(its,isph).lt.one) then
      do ij = ddx_data % inl(isph), ddx_data % inl(isph+1)-1
        jsph = ddx_data % nl(ij)

        ! compute geometrical variables
        vij  = ddx_data % csph(:,isph) + ddx_data % rsph(isph)*ddx_data % cgrid(:,its) - ddx_data % csph(:,jsph)
        vvij = sqrt(dot_product(vij,vij))
        tij  = vvij/ddx_data % rsph(jsph)
        if ( tij.lt.one ) then
          sij = vij/vvij
          call ylmbas(sij, rho, ctheta, stheta, cphi, &
                      & sphi, ddx_data % lmax, &
                      & ddx_data % vscales, basloc, &
                      & vplm, vcos, vsin)
          call inthsp(ddx_data, vvij, ddx_data % rsph(jsph), jsph, basloc, fac_hsp)
          xij = fsw(tij, ddx_data % se, ddx_data % eta)
          if (ddx_data % fi(its,isph).gt.one) then
            oij = xij/ddx_data % fi(its, isph)
          else
            oij = xij
          end if
          pot(its) = pot(its) + oij*dot_product(fac_hsp,x(:,jsph))
        end if
      end do
    end if
  end do
  endsubroutine calcv2_lpb


  !
  ! Intermediate calculation in calcv2_lpb subroutine
  ! @param[in]  rijn    : Radius of sphers x_ijn
  ! @param[in]  ri      : Radius of sphers x_i
  ! @param[in]  isph    : Index of sphere
  ! @param[in]  basloc  : Spherical Harmonic
  ! @param[out] fac_hsp : Return bessel function ratio multiplied by 
  !                       the spherical harmonic Y_l'm'. Array of size nylm
  !
  subroutine inthsp(ddx_data, rijn, ri, isph, basloc, fac_hsp)
  use bessel
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  integer, intent(in) :: isph
  real(dp), intent(in) :: rijn, ri
  real(dp), dimension(ddx_data % nbasis), intent(in) :: basloc
  real(dp), dimension(ddx_data % nbasis), intent(inout) :: fac_hsp
  real(dp), dimension(0:ddx_data % lmax) :: SI_rijn, DI_rijn
  integer :: l, m, ind, NM

  SI_rijn = 0
  DI_rijn = 0
  fac_hsp = 0

  ! Computation of modified spherical Bessel function values      
  call SPHI_bessel(ddx_data % lmax,rijn*ddx_data % kappa,NM,SI_rijn,DI_rijn)
  
  do l = 0, ddx_data % lmax
    do  m = -l, l
      ind = l*l + l + 1 + m
      if ((SI_rijn(l).lt.tol_zero)) then
        write(*,*) 'DEBUG : Case 1'
        fac_hsp(ind) = (rijn/ri)**l*basloc(ind)
      else if (SI_ri(l,isph).lt.tol_zero) then
        write(*,*) 'DEBUG : Case 2'
        fac_hsp(ind) = (rijn/ri)**l*basloc(ind)
      else if(SI_rijn(l)/SI_ri(l,isph).gt.(rijn/ri)**l) then
        write(*,*) 'DEBUG : Case 3'
        fac_hsp(ind) = (rijn/ri)**l*basloc(ind)
      else if ( SI_ri(l,isph).gt.tol_inf) then
        write(*,*) 'DEBUG : Case 4'
        fac_hsp(ind) = zero
      else
        fac_hsp(ind) = SI_rijn(l)/SI_ri(l,isph)*basloc(ind)
      end if
    end do
  end do
  endsubroutine inthsp


  subroutine intcosmo(ddx_data, tij, basloc, fac_cosmo)
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  real(dp),  intent(in) :: tij
  real(dp), dimension(ddx_data % nbasis), intent(in) :: basloc
  real(dp), dimension(ddx_data % nbasis), intent(inout) :: fac_cosmo
  integer :: l, m, ind
  do l = 0, ddx_data % lmax
    do  m = -l, l
        ind = l*l + l + 1 + m
        fac_cosmo(ind) = tij**l*basloc(ind)
    end do
  end do
  end subroutine intcosmo

  !
  ! Update the RHS in outer iteration
  ! @param[in] rhs_cosmo_init : G_0
  ! @param[in] rhs_hsp_init   : F_0
  ! @param[in, out] rhs_cosmo : -C_1*X_r^(k-1) - C_2*X_e^(k-1) + G_0 + F_0
  ! @param[in, out] rhs_hsp   : -C_1*X_r^(k-1) - C_2*X_e^(k-1) + F_0
  ! @param[in] Xr             : X_r^(k-1)
  ! @param[in] Xe             : X_e^(k-1)
  !
  subroutine update_rhs(ddx_data, rhs_cosmo_init, rhs_hsp_init, rhs_cosmo, & 
      & rhs_hsp, Xr, Xe)
  use bessel
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph), intent(in) :: rhs_cosmo_init, &
      & rhs_hsp_init
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph), intent(inout) :: rhs_cosmo, rhs_hsp
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph) :: rhs_plus
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph), intent(in) :: Xr, Xe
  integer :: isph, jsph, ig, kep, ind, l1,m1, ind1, ind0, count, istatus
  real(dp), dimension(3) :: vij
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph) :: diff_re
  real(dp), dimension(nbasis0,ddx_data % nsph) :: diff0
  real(dp), dimension(ddx_data % nbasis,ddx_data % nbasis,ddx_data % nsph) :: smat
  real(dp), dimension(ddx_data % ncav) :: diff_ep
  real(dp) :: Qval, rijn, val
  integer :: c0, cr, c_qmat, c_init, c_ep0, c_ep1 !, nbasis_appro
      
  if (first_out_iter) then
    allocate(coefvec(ddx_data % ngrid, ddx_data % nbasis, ddx_data % nsph), &
            Pchi(ddx_data % nbasis, nbasis0, ddx_data % nsph), stat = istatus)
    if (istatus.ne.0) then
      write(*,*)'update_rhs : [1] allocation failed!'
      stop
    end if
  end if
  
  ! Compute P_chi matrix, Eq.(87)
  ! TODO: probably has to be declared somewhere
  ! and i have to recover mkpmat 
  if (first_out_iter) then      
    do jsph = 1, ddx_data % nsph
      call mkpmat(ddx_data, jsph, Pchi(:,:,jsph))
    end do    
  end if 
      
  ! Precompute coefvec of Qmat, Cost: linear scaling
  ! TODO: remove all the precomputations
  ! or, if they are really necessary, do them in a separate subroutine
  if (first_out_iter) then 
    do isph = 1,ddx_data % nsph
      do ig = 1,ddx_data % ngrid
        if (ddx_data % ui(ig, isph) .gt. 0) then
          do ind  = 1, ddx_data % nbasis 
            coefvec(ig,ind,isph) = ddx_data % wgrid(ig)*ddx_data % ui(ig,isph)*ddx_data % vgrid(ind,ig)
          end do
        end if
      end do
    end do
  end if
      
  ! precompute termimat
  ! TODO: same as before
  if (first_out_iter) then
    do jsph = 1, ddx_data % nsph
      do l1 = 0, ddx_data % lmax
        if (max(DI_ri(l1,jsph),SI_ri(l1,jsph)).gt.tol_inf) then
          termimat(l1,jsph) = ddx_data % kappa
        else if (min(DI_ri(l1,jsph),SI_ri(l1,jsph)).lt.tol_zero) then
          termimat(l1,jsph) = l1/ddx_data % rsph(jsph) + &
              & (l1 + 1)*(ddx_data % kappa**2*ddx_data % rsph(jsph))/((two*l1 + one) * &
              & (two*l1 + three))
        else
          termimat(l1,jsph) = DI_ri(l1,jsph)/SI_ri(l1,jsph)*ddx_data % kappa
        end if
      end do
    end do
  end if

  if (first_out_iter) then
    if (ddx_data % iprint .gt. 0) then  
      write(*,999) dble(c_init-c0)/dble(cr)
 999  format('Time of initializing Pchi, coefvec, termi: ',f8.3, &
          & ' secs.')
    end if
  end if

  ! diff_re = epsp/eps*l1/ri*Xr - i'(ri)/i(ri)*Xe,
  diff_re = zero 
  do jsph = 1, ddx_data % nsph
    do l1 = 0, ddx_data % lmax
      do m1 = -l1,l1
        ind1 = l1**2 + l1 + m1 + 1
        diff_re(ind1,jsph) = ((epsp/ddx_data % eps)*(l1/ddx_data % rsph(jsph)) * &
            & Xr(ind1,jsph) - termimat(l1,jsph)*Xe(ind1,jsph))
      end do
    end do
  end do

  ! diff0 = Pchi * diff_er, linear scaling
  ! TODO: probably doing PX on the fly is better 
  diff0 = zero 
  do jsph = 1, ddx_data % nsph
    do ind0 = 1, nbasis0
      diff0(ind0, jsph) = dot_product(diff_re(:,jsph), &
          & Pchi(:,ind0, jsph))
    end do
  end do

  ! diff_ep = diff0 * coefY,    COST: M^2*nbasis*Nleb
  diff_ep = zero
  do kep = 1, ddx_data % ncav
    val = zero
    do jsph = 1, ddx_data % nsph 
      do ind0 = 1, nbasis0
        val = val + diff0(ind0,jsph)*coefY(kep,ind0,jsph)
      end do
    end do
    diff_ep(kep) = val 
  end do

  rhs_plus = zero
  kep = 0
  do isph = 1, ddx_data % nsph
    do ig = 1, ddx_data % ngrid
      if (ddx_data % ui(ig,isph).gt.zero) then 
        kep = kep + 1
        do ind = 1, ddx_data % nbasis
          rhs_plus(ind,isph) = rhs_plus(ind,isph) + &
              & coefvec(ig,ind,isph)*diff_ep(kep)
        end do
      end if
    end do
  end do

  rhs_cosmo = rhs_cosmo_init - rhs_plus
  rhs_hsp = rhs_hsp_init - rhs_plus

  return
  end subroutine update_rhs  

  !
  ! Computation of P_chi
  ! @param[in]  isph : Sphere number
  ! @param[out] pmat : Matrix of size nbasis X (lmax0+1)^2, Fixed lmax0
  !
  subroutine mkpmat(ddx_data, isph, pmat)
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  integer,  intent(in) :: isph
  real(dp), dimension(ddx_data % nbasis, (lmax0+1)**2), intent(inout) :: pmat
  integer :: l, m, ind, l0, m0, ind0, its, nbasis0
  real(dp)  :: f, f0
  pmat(:,:) = zero
  do its = 1, ddx_data % ngrid
    if (ddx_data % ui(its,isph).ne.0) then
      do l = 0, ddx_data % lmax
        ind = l*l + l + 1
        do m = -l,l
          f = ddx_data % wgrid(its) * ddx_data % vgrid(ind+m,its) * ddx_data % ui(its,isph)
          do l0 = 0, lmax0
            ind0 = l0*l0 + l0 + 1
            do m0 = -l0, l0
              f0 = ddx_data % vgrid(ind0+m0,its)
              pmat(ind+m,ind0+m0) = pmat(ind+m,ind0+m0) + f * f0
            end do
          end do
        end do
      end do
    end if
  end do
  endsubroutine mkpmat
  
  !
  ! Computation of Adjoint
  ! @param[in] ddx_data: Input data file
  ! @param[in] psi     : psi_r
  subroutine ddx_lpb_adjoint(ddx_data, psi, Xadj_r, Xadj_e)
  implicit none
  type(ddx_type), intent(in)           :: ddx_data
  real(dp), intent(in)                 :: psi(ddx_data % nbasis, ddx_data % nsph)
  real(dp), intent(inout), dimension(ddx_data % nbasis,ddx_data % nsph) :: Xadj_r(:,:), Xadj_e(:,:)
  ! Local Variables
  ! ok             : Logical expression used in Jacobi solver
  ! n_iter         : Maximum number of iteration
  ! istatus        : Status of allocation
  ! isph           : Index for the sphere
  ! ibasis         : Index for the basis
  ! iteration      : Number of outer loop iterations
  ! inc            : Check for convergence threshold
  ! relative_num   : Numerator of the relative error
  ! relative_denom : Denominator of the relative error
  ! converged      : Convergence check for outer loop
  logical   :: ok = .false.
  integer   :: n_iter, istatus, isph, ibasis,  iteration
  real(dp)  :: inc, relative_num, relative_denom
  logical   :: converged = .false.
  
  !!
  !! Xadj_r         : Adjoint solution of Laplace equation
  !! Xadj_e         : Adjoint solution of HSP equation
  !! rhs_cosmo      : RHS corresponding to Laplace equation
  !! rhs_hsp        : RHS corresponding to HSP equation
  !! rhs_cosmo_init : Initial RHS for Laplace, psi_r in literature
  !! rhs_hsp_init   : Initial RHS for HSP, psi_e = 0 in literature
  !! X_r_k_1          : Solution of previous iterative step, holds Xadj_r_k_1
  !! X_e_k_1          : Solution of previous iterative step, holds Xadj_e_k_1
  real(dp), allocatable :: rhs_cosmo(:,:), rhs_hsp(:,:), &
                        & rhs_cosmo_init(:,:), rhs_hsp_init(:,:), &
                        & X_r_k_1(:,:), X_e_k_1(:,:)
                        
  ! Allocation
  allocate(rhs_cosmo(ddx_data % nbasis, ddx_data % nsph), &
           & rhs_hsp(ddx_data % nbasis, ddx_data % nsph), &
           & rhs_cosmo_init(ddx_data % nbasis, ddx_data % nsph), &
           & rhs_hsp_init(ddx_data % nbasis, ddx_data % nsph),&
           & X_r_k_1(ddx_data % nbasis, ddx_data % nsph),&
           & X_e_k_1(ddx_data % nbasis, ddx_data % nsph),&
           & stat = istatus)
  if(istatus .ne. 0) then
    write(*,*) 'Allocation failed in Forces LPB!'
    stop
  end if
  write(*,*) 'Computation of Forces for ddLPB'
  ! We compute the adjoint solution first
  write(*,*) 'Solution of adjoint system'
  ! Set maximum number of iteration
  n_iter = ddx_data % maxiter
  ! Set local variables
  isph = one
  ibasis = one
  iteration = one
  inc = zero
  relative_num = zero
  relative_denom = zero
  X_r_k_1 = zero
  X_e_k_1 = zero
  ! Logical variable for the first outer iteration
  first_out_iter = .true.
  ! Initial RHS
  ! rhs_cosmo_init = psi_r
  ! rhs_hsp_init   = psi_e (=0)
  rhs_cosmo_init = psi
  rhs_hsp_init = zero
  ! Updated RHS
  rhs_cosmo = rhs_cosmo_init
  rhs_hsp = zero
  ! Initial Xadj_r and Xadj_e
  Xadj_r = zero
  Xadj_e = zero
  ! Solve the adjoint system  
  do while (.not.converged)
    ! Solve A*X_adj_r = psi_r
    call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, ddx_data % ndiis, &
                    & 4, ddx_data % tol, rhs_cosmo, Xadj_r, n_iter, ok, lstarx, &
                    & ldm1x, hnorm)
    ! Remove the factor of 4pi/(2l+1)
    call convert_ddcosmo(ddx_data, 1, Xadj_r)
    
    ! Solve the HSP equation
    ! B*X_adj_e = psi_e
    ! For first iteration the rhs is zero for HSP equation. Hence, Xadj_e = 0
    if(iteration.ne.1) then
      call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, ddx_data % ndiis, &
                      & 4, ddx_data % tol, rhs_hsp, Xadj_e, n_iter, ok, bstarx, &
                      & ldm1x, hnorm)
      ! Remove the factor of 4pi/(2l+1) on diagonal
      call convert_ddlpb(ddx_data, Xadj_e)
    endif
    
    ! Update the RHS
    ! |rhs_r| = |psi_r|-|C1* C1*||Xadj_r|
    ! |rhs_e| = |psi_e| |C2* C2*||Xadj_e|
    call update_rhs_adj(ddx_data, rhs_cosmo_init, rhs_hsp_init,&
                        & rhs_cosmo, rhs_hsp, Xadj_r, Xadj_e)
    ! Stopping Criteria.
    ! Checking the relative error of Xadj_r
    relative_num = zero
    relative_denom = zero
    do isph = 1, ddx_data % nsph
      do ibasis = 1, ddx_data % nbasis
        relative_num = relative_num + (Xadj_r(ibasis, isph)-X_r_k_1(ibasis, isph))**2 &
                       & + (Xadj_e(ibasis, isph)-X_e_k_1(ibasis, isph))**2
        relative_denom = relative_denom + X_r_k_1(ibasis, isph)**2+X_e_k_1(ibasis, isph)**2
      end do
    end do
    relative_num = sqrt(relative_num)
    relative_denom = sqrt(relative_denom)

    ! Check for convergence
    inc = relative_num/relative_denom
    ! Store the previous step solution
    X_r_k_1 = Xadj_r
    X_e_k_1 = Xadj_e
    if ((iteration .gt. 1) .and. (inc.lt.ddx_data % tol)) then
      write(6,*) 'Reach tolerance.'
      converged = .true.
    end if
    write(6,*) iteration, inc
    iteration = iteration + 1
    first_out_iter = .false.
  end do
  if (ddx_data % iprint .ge. 5) then
    call prtsph('Xadj_r', ddx_data % nbasis, ddx_data % lmax, &
          & ddx_data % nsph, 0, Xadj_r)
    call prtsph('Xadj_e', ddx_data % nbasis, ddx_data % lmax, &
          & ddx_data % nsph, 0, Xadj_e)
  end if
  ! Deallocation
  deallocate(rhs_cosmo, rhs_hsp, rhs_cosmo_init, &
           & rhs_hsp_init, X_r_k_1, X_e_k_1, stat = istatus)
  if(istatus .ne. 0) then
    write(*,*) 'Deallocation failed in Forces LPB!'
    stop
  end if
  end subroutine ddx_lpb_adjoint
  
  !! Computation of Adjoint B, i.e., B*
  !> Apply adjoint single layer operator to spherical harmonics
  !! implementation is similar to lstarx in ddCOSMO
  !! Diagonal blocks are not counted here.
  subroutine bstarx(ddx_data, x, y)
  ! Inputs
  type(ddx_type), intent(in) :: ddx_data
  real(dp), intent(in)       :: x(ddx_data % nbasis, ddx_data % nsph)
  ! Output
  real(dp), intent(out)      :: y(ddx_data % nbasis, ddx_data % nsph)
  ! Local variables
  integer                    :: isph, igrid, istatus, l, ind
  real(dp), allocatable      :: xi(:,:), vplm(:), basloc(:), vcos(:), vsin(:)
  ! Allocate workspaces
  allocate(xi(ddx_data % ngrid, ddx_data % nsph), vplm(ddx_data % nbasis), &
        & basloc(ddx_data % nbasis), vcos(ddx_data % lmax+1), &
        & vsin(ddx_data % lmax+1), stat=istatus)
  if (istatus .ne. 0) then
      write(*, *) 'bstarx: allocation failed!'
      stop
  endif
  if (ddx_data % iprint .ge. 5) then
      call prtsph('X', ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, &
          & x)
  end if
  ! Initalize
  y = zero
  !! Expand x over spherical harmonics
  ! Loop over spheres      
  do isph = 1, ddx_data % nsph
      ! Loop over grid points
      do igrid = 1, ddx_data % ngrid
          xi(igrid, isph) = dot_product(x(:, isph), &
              & ddx_data % vgrid(:ddx_data % nbasis, igrid))
      end do
  end do
  !! Compute action
  ! Loop over spheres
  do isph = 1, ddx_data % nsph
      ! Compute NEGATIVE action of off-digonal blocks
      call adjrhs_lpb(ddx_data, isph, xi, y(:, isph), basloc, vplm, vcos, vsin)
      y(:, isph) = - y(:, isph)
  end do
  if (ddx_data % iprint .ge. 5) then
      call prtsph('B*X (off-diagonal)', ddx_data % nbasis, ddx_data % lmax, &
          & ddx_data % nsph, 0, y)
  end if
  !! Diagonal preconditioning for bstarx
  ! NOTE: Activate this if one is using GMRES solver
  !       Jacobi solver uses ldm1x for diagonals. Declare l and ind above.
  !do l = 0, ddx_data % lmax
  !  ind = l*l + l + 1
  !  y(ind-l:ind+l, :) = x(ind-l:ind+l, :) * (ddx_data % vscales(ind)**2)
  !end do
  ! Deallocate workspaces
  deallocate( xi, basloc, vplm, vcos, vsin , stat=istatus )
  if ( istatus.ne.0 ) then
      write(*,*) 'bstarx: allocation failed !'
      stop
  endif
  end subroutine bstarx

  !
  ! Taken from ddx_core routine adjrhs
  ! Called from bstarx
  ! Compute the Adjoint matix B*x
  !
  subroutine adjrhs_lpb( ddx_data, isph, xi, vlm, basloc, vplm, vcos, vsin )
  implicit none
  type(ddx_type), intent(in) :: ddx_data
  integer,  intent(in)    :: isph
  real(dp), dimension(ddx_data % ngrid, ddx_data % nsph), intent(in) :: xi
  real(dp), dimension((ddx_data % lmax+1)**2), intent(inout) :: vlm
  real(dp), dimension((ddx_data % lmax+1)**2), intent(inout) :: basloc, vplm
  real(dp), dimension(ddx_data % lmax+1), intent(inout) :: vcos, vsin

  integer :: ij, jsph, ig, l, ind, m
  real(dp)  :: vji(3), vvji, tji, sji(3), xji, oji, fac
  real(dp) :: rho, ctheta, stheta, cphi, sphi
  real(dp), dimension(ddx_data % nbasis) :: fac_hsp

  !loop over neighbors of i-sphere
  do ij = ddx_data % inl(isph), ddx_data % inl(isph+1)-1
    !j-sphere is neighbor
    jsph = ddx_data % nl(ij)
    !loop over integration points
    do ig = 1, ddx_data % ngrid
      !compute t_n^ji = | r_j + \rho_j s_n - r_i | / \rho_i
      vji  = ddx_data % csph(:,jsph) + ddx_data % rsph(jsph)* &
            & ddx_data % cgrid(:,ig) - ddx_data % csph(:,isph)
      vvji = sqrt(dot_product(vji,vji))
      tji  = vvji/ddx_data % rsph(isph)
      !point is INSIDE i-sphere (+ transition layer)
      if ( tji.lt.( one + (ddx_data % se+one)/two*ddx_data % eta ) ) then
        !compute s_n^ji
        sji = vji/vvji
        call ylmbas(sji, rho, ctheta, stheta, cphi, &
                      & sphi, ddx_data % lmax, &
                      & ddx_data % vscales, basloc, &
                      & vplm, vcos, vsin)
        call inthsp_adj(ddx_data, vvji, ddx_data % rsph(isph), isph, basloc, fac_hsp)
        !compute \chi( t_n^ji )
        xji = fsw( tji, ddx_data % se, ddx_data % eta )
        !compute W_n^ji
        if ( ddx_data % fi(ig,jsph).gt.one ) then
          oji = xji/ ddx_data % fi(ig,jsph)
        else
          oji = xji
        endif
        !compute w_n * xi(n,j) * W_n^ji
        fac = ddx_data % wgrid(ig) * xi(ig,jsph) * oji
        !loop over l
        do l = 0, ddx_data % lmax
          ind  = l*l + l + 1
          !loop over m
            do m = -l,l
              vlm(ind+m) = vlm(ind+m) + fac*fac_hsp(ind+m)
            enddo
        enddo
      endif
    enddo
  enddo
  end subroutine adjrhs_lpb
  
  !
  ! Intermediate calculation in adjrhs_lpb subroutine
  ! @param[in]  rijn    : Radius of sphers x_ijn
  ! @param[in]  ri      : Radius of sphers x_i
  ! @param[in]  isph    : Index of sphere
  ! @param[in]  basloc  : Spherical Harmonic
  ! @param[out] fac_hsp : Return bessel function ratio multiplied by 
  !                       the spherical harmonic Y_l'm'. Array of size nylm
  !
  subroutine inthsp_adj(ddx_data, rjin, rj, jsph, basloc, fac_hsp)
  use bessel
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  integer, intent(in) :: jsph
  real(dp), intent(in) :: rjin, rj
  real(dp), dimension(ddx_data % nbasis), intent(in) :: basloc
  real(dp), dimension(ddx_data % nbasis), intent(inout) :: fac_hsp
  real(dp), dimension(0:ddx_data % lmax) :: SI_rjin, DI_rjin
  integer :: l, m, ind, NM

  SI_rjin = 0
  DI_rjin = 0
  fac_hsp = 0

  ! Computation of modified spherical Bessel function values      
  call SPHI_bessel(ddx_data % lmax,rjin*ddx_data % kappa,NM,SI_rjin,DI_rjin)
  
  do l = 0, ddx_data % lmax
    do  m = -l, l
      ind = l*l + l + 1 + m
      if ((SI_rjin(l).lt.zero) .or. (SI_ri(l,jsph).lt.tol_zero) &
          & .or. (SI_rjin(l)/SI_ri(l,jsph).gt.(rjin/rj)**l)) then
        fac_hsp(ind) = (rjin/rj)**l*basloc(ind)
      else if ( SI_ri(l,jsph).gt.tol_inf) then
        fac_hsp(ind) = zero
      else
        fac_hsp(ind) = SI_rjin(l)/SI_ri(l,jsph)*basloc(ind)
      end if
    end do
  end do
  endsubroutine inthsp_adj
  
  !
  ! Update the RHS in outer iteration for adjoint system
  ! @param[in] rhs_r_init : psi_r
  ! @param[in] rhs_e_init : psi_e = 0
  ! @param[in, out] rhs_r : -C_1*\times Xadj_r^(k-1) - C_1*\times Xadj_e^(k-1)
  !                         + psi_r
  ! @param[in, out] rhs_e : -C_2*\times Xadj_r^(k-1) - C_2*\times Xadj_e^(k-1)
  ! @param[in] Xadj_r     : Xadj_r^(k-1)
  ! @param[in] Xadj_e     : Xadj_e^(k-1)
  !
  subroutine update_rhs_adj(ddx_data, rhs_r_init, rhs_e_init, rhs_r, & 
      & rhs_e, Xadj_r, Xadj_e)
  use bessel
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph), intent(in) :: rhs_r_init, &
                                                                       & rhs_e_init
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph), intent(in) :: Xadj_r, Xadj_e
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph), intent(inout) :: rhs_r, rhs_e
  ! Local Variables
  ! rhs_r_adj : C_1*(Xadj_r+Xadj_e)
  ! rhs_e_adj : C_2*(Xadj_r+Xadj_e)
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph) :: rhs_r_adj, rhs_e_adj
  ! isph    : Index for the sphere
  ! igrid   : Index for the grid points
  ! iext    : Index for the external grid points
  ! ibasis  : Index for the basis
  ! ibasis0 : Index for the fixed basis (nbasis0)
  ! il      : Index for lmax
  ! im      : Index for m:-l,...,l
  ! ind1    : l^2+l+m+1
  integer :: isph, igrid, iext, ibasis, ibasis0, il, im, ind1
  ! epsilon_ratio : epsilon_1/epsilon_2
  real(dp) :: epsilon_ratio
  ! val_1, val_2  : Intermediate summations variable
  real(dp) :: val_1, val_2  

  ! NOTE: Here allocation of coefvec does not take place, as they
  !       are global variables and hence only need to be allocated once
  !       which they have been in update_rhs.
      
  ! Compute the term 
  !    M        N_g
  !   ___ ___   ___     i  i           l    in
  !   \   \     \   w  U (x )Y    (s )---[Q]    [X +X ]
  !   /__ /__   /__  n     n  l'm'  n r_j   jlm   r  e il'm'
  !   i=1 l'm'  n=1
  ! \sum_{i,\ell',m'}\left[epsilon_1/epsilon_2\left\lbrace \sum_{n=1}^{N_g}
  ! \omega_n\chi_i((x_i^n)Y_{\ell'm'}(s_n)\right\rbrace*
  ! \left(Xadj_r(i\ell'm') - Xadj_e(i\ell'm')\right)\right]
  do isph = 1,ddx_data % nsph
    do igrid = 1,ddx_data % ngrid
      if (ddx_data % ui(igrid, isph) .gt. 0) then
        do ibasis  = 1, ddx_data % nbasis 
          coefvec(igrid,ibasis,isph) = ddx_data % wgrid(igrid)*&
                                & ddx_data % ui(igrid,isph)*ddx_data%vgrid(ibasis,igrid)* &
                                & (Xadj_r(ibasis,isph) + Xadj_e(ibasis,isph))
        end do
      end if
    end do
  end do

  !             i'_(l)(r_j)
  ! Termimat = -------------
  !             i_(l)(r_j)
  ! We do not have to compute Pchi and termimat as they have already being computed

  ! NOTE: These remain constant through the outer iteration and hence needs to be computed
  !       once. 
  if(first_out_iter) then
    ! Compute
    ! diff_re_c1 = \ell/rj
    ! diff_re_c2 = -i'_(\ell)(r_j)/i_(\ell)(r_j)
    ! epsilon_1/epsilon_2
    epsilon_ratio = epsp/ddx_data % eps
    diff_re_c1 = zero
    diff_re_c2 = zero
    do isph = 1, ddx_data % nsph
      do il = 0, ddx_data % lmax
        do im = -il,il
          ind1 = il**2 + il + im + 1
          diff_re_c1(ind1,isph) = epsilon_ratio*(il/ddx_data % rsph(isph))
          diff_re_c2(ind1,isph) = -termimat(il, isph)
        end do
      end do
    end do

    ! Compute
    ! diff0_c1 = Pchi * (\ell/rj)
    ! diff0_c2 = Pchi * (i'_(\ell)(r_j)/i_(\ell)(r_j))
    diff0_c1 = zero
    diff0_c2 = zero
    do isph = 1, ddx_data % nsph
      do ibasis0 = 1, nbasis0
        diff0_c1(ibasis0, isph) = dot_product(diff_re_c1(:,isph), &
            & Pchi(:,ibasis0, isph))
        diff0_c2(ibasis0, isph) = dot_product(diff_re_c2(:,isph), &
            & Pchi(:,ibasis0, isph))
      end do
    end do

    ! Compute 
    ! diff_ep_c1 = diff0_c1 * coefY,    COST: M^2*nbasis*Nleb
    ! diff_ep_c2 = diff0_c2 * coefY,    COST: M^2*nbasis*Nleb
    ! coefY = (i'_{\ell_0}(r_j)/i_{\ell_0}(r_j)-l'_{\ell_0}(r_j)/k_{\ell_0}(r_j))^{-1}
    !         x (k_{\ell_0}(r_jin)/k_{\ell_0}(r_j))Y_{\ell_0m_0}(s_jin)
    diff_ep_c1 = zero
    diff_ep_c2 = zero
    do iext = 1, ddx_data % ncav
      val_1 = zero
      val_2 = zero
      do isph = 1, ddx_data % nsph 
        do ibasis0 = 1, nbasis0
          val_1 = val_1 + diff0_c1(ibasis0,isph)*coefY(iext,ibasis0,isph)
          val_2 = val_2 + diff0_c2(ibasis0,isph)*coefY(iext,ibasis0,isph)
        end do
      end do
      diff_ep_c1(iext) = val_1
      diff_ep_c2(iext) = val_2
    end do
  endif

  rhs_r_adj = zero
  rhs_e_adj = zero
  iext = 0
  do isph = 1, ddx_data % nsph
    do igrid = 1, ddx_data % ngrid
      if (ddx_data % ui(igrid,isph).gt.zero) then 
        iext = iext + 1
        do ibasis = 1, ddx_data % nbasis
          rhs_r_adj(ibasis,isph) = rhs_r_adj(ibasis,isph) + &
              & coefvec(igrid,ibasis,isph)*diff_ep_c1(iext)
          rhs_e_adj(ibasis,isph) = rhs_e_adj(ibasis,isph) + &
              & coefvec(igrid,ibasis,isph)*diff_ep_c2(iext)
        end do
      end if
    end do
  end do

  ! Update the RHS
  rhs_r = rhs_r_init - rhs_r_adj
  rhs_e = rhs_e_init - rhs_e_adj

  return
  end subroutine update_rhs_adj
  
  !
  ! Subroutine to compute K^A counterpart for the HSP equation. Similar to fdoka.
  ! @param[in]  ddx_data  : Data type
  ! @param[in]  isph      : Index of sphere
  ! @param[in]  Xe        : Solution vector Xe
  ! @param[in]  Xadj_e    : Adjoint solution on evaluated on grid points Xadj_e_sgrid
  ! @param[in]  basloc    : Spherical harmonics Y_lm
  ! @param[in]  dbasloc   : Derivative of spherical harmonics \nabla^i(Y_lm)
  ! @param[in]  vplm      : Argument to call ylmbas
  ! @param[in]  vcos      : Argument to call ylmbas
  ! @param[in]  vsin      : Argument to call ylmbas
  ! @param[out] force_e   : Force of adjoint part
  subroutine fdoka_b_xe(ddx_data, isph, Xe, Xadj_e, basloc, dbsloc, &
                       & vplm, vcos, vsin, force_e)
  use bessel
  implicit none
  type(ddx_type), intent(in) :: ddx_data
  integer,                         intent(in)    :: isph
  real(dp),  dimension(ddx_data % nbasis, ddx_data % nsph), intent(in)   :: Xe
  real(dp),  dimension(ddx_data % ngrid),       intent(in)    :: Xadj_e
  real(dp),  dimension(ddx_data % nbasis),      intent(inout) :: basloc, vplm
  real(dp),  dimension(3, ddx_data % nbasis),    intent(inout) :: dbsloc
  real(dp),  dimension(ddx_data % lmax+1),      intent(inout) :: vcos, vsin
  real(dp),  dimension(3),           intent(inout) :: force_e
  
  ! Local Variables
  ! igrid  : Index of grid point
  ! ineigh : Index over Row space of neighbors
  ! jsph   : Index of neighbor sphere
  ! l     : Index for l, l:0,...,lmax
  ! m     : Index for m, m:-l,...,l
  ! ind    : l*l+l+1
  ! NM     : Argument for calling Bessel functions
  integer :: igrid, ineigh, jsph, l, ind, m, NM
  ! SI_rijn : Besssel function of first kind for rijn
  ! DI_rijn : Derivative of Besssel function of first kind for rijn
  real(dp), dimension(0:ddx_data % lmax) :: SI_rijn, DI_rijn
  ! rijn   : r_j*r_1^j(x_i^n) = |x_i^n-x_j|
  ! tij    : r_1^j(x_i^n)
  ! beta   : Eq.(53) Stamm.etal.18
  ! tlow   : Lower bound for switch region
  ! thigh  : Upper bound for switch region
  ! xij    : chi_j(x_i^n)
  ! oij    : omega^\eta_ijn
  ! f1     : First factor in alpha computation
  ! f2     : Second factor in alpha computation
  real(dp)  :: rijn, tij, beta, tlow, thigh, xij, oij, f1, f2, f3
  ! vij   : x_i^n-x_j
  ! sij   : e^j(x_i^n)
  ! alpha : Eq.(52) Stamm.etal.18
  ! va    : Eq.(54) Stamm.etal.18
  real(dp)  :: vij(3), sij(3), alpha(3), va(3), rj
  real(dp), external :: dnrm2
  
  SI_rijn = 0
  DI_rijn = 0
  
  tlow  = one - pt5*(one - ddx_data % se)*ddx_data % eta
  thigh = one + pt5*(one + ddx_data % se)*ddx_data % eta
    
  ! Loop over grid points
  do igrid = 1, ddx_data % ngrid
    va = zero
    do ineigh = ddx_data % inl(isph), ddx_data % inl(isph+1) - 1
      jsph = ddx_data % nl(ineigh)
      vij  = ddx_data % csph(:,isph) + &
            & ddx_data % rsph(isph)*ddx_data % cgrid(:,igrid) - &
            & ddx_data % csph(:,jsph)
      rijn = dnrm2(3, vij, 1)
      tij  = rijn/ddx_data % rsph(jsph)
      rj = ddx_data % rsph(jsph)
     
      if (tij.ge.thigh) cycle
      ! Computation of modified spherical Bessel function values      
      call SPHI_bessel(ddx_data % lmax, rijn*ddx_data % kappa, NM, SI_rijn, DI_rijn)
    
      sij  = vij/rijn
      !call dbasis(sij,basloc,dbsloc,vplm,vcos,vsin)
      call dbasis(ddx_data, sij, basloc, dbsloc, vplm, vcos, vsin)
      alpha  = zero
      do l = 0, ddx_data % lmax
        ind = l*l + l + 1
        if ((SI_rijn(l).lt.zero) .or. (SI_ri(l,jsph).lt.tol_zero) &
          & .or. (SI_rijn(l)/SI_ri(l,jsph).gt.(rijn/rj)**l)) then
          f2 = (rijn/rj)**l
          f1 = (l*(rijn/rj)**(l-1))/rj
        else if ( SI_ri(l,jsph).gt.tol_inf) then
          f2 = zero
          f1 = zero
        else
          f2 = SI_rijn(l)/(rijn*SI_ri(l, jsph))
          f1 = (DI_rijn(l)*ddx_data % kappa)/SI_ri(l, jsph)
        end if
        do m = -l, l
          alpha(:) = alpha(:) + (f1*sij(:)*basloc(ind+m) + &
                    & f2*dbsloc(:,ind+m))*Xe(ind+m,jsph)
        end do
      end do
      beta = compute_beta(ddx_data, SI_rijn, rijn, jsph, Xe(:,jsph),basloc)
      xij = fsw(tij, ddx_data % se, ddx_data % eta)
      if (ddx_data % fi(igrid,isph).gt.one) then
        oij = xij/ddx_data % fi(igrid,isph)
        f2  = -oij/ddx_data % fi(igrid,isph)
      else
        oij = xij
        f2  = zero
      end if
      f1 = oij
      va(:) = va(:) + f1*alpha(:) + beta*f2*ddx_data % zi(:,igrid,isph)
      if (tij .gt. tlow) then
        f3 = beta*dfsw(tij,ddx_data % se,ddx_data % eta)/ddx_data % rsph(jsph)
        if (ddx_data % fi(igrid,isph).gt.one) f3 = f3/ddx_data % fi(igrid,isph)
        va(:) = va(:) + f3*sij(:)
      end if
    end do
  force_e = force_e - ddx_data % wgrid(igrid)*Xadj_e(igrid)*va(:)
  end do
  return
  end subroutine fdoka_b_xe
  
  !
  ! Subroutine to compute K^A+K^C counterpart for the HSP equation. Similar to fdokb.
  ! @param[in]  ddx_data  : Data type
  ! @param[in]  isph      : Index of sphere
  ! @param[in]  Xe        : Solution vector Xe
  ! @param[in]  Xadj_e    : Adjoint solution on evaluated on grid points Xadj_e_sgrid
  ! @param[in]  basloc    : Spherical harmonics Y_lm
  ! @param[in]  dbasloc   : Derivative of spherical harmonics \nabla^i(Y_lm)
  ! @param[in]  vplm      : Argument to call ylmbas
  ! @param[in]  vcos      : Argument to call ylmbas
  ! @param[in]  vsin      : Argument to call ylmbas
  ! @param[out] force_e   : Force of adjoint part
  subroutine fdokb_b_xe(ddx_data, isph, Xe, Xadj_e, basloc, dbsloc, &
                        & vplm, vcos, vsin, force_e)
  use bessel
  implicit none
  type(ddx_type), intent(in) :: ddx_data
  integer,                         intent(in)    :: isph
  real(dp),  dimension(ddx_data % nbasis, ddx_data % nsph), intent(in)    :: Xe
  real(dp),  dimension(ddx_data % ngrid, ddx_data % nsph),  intent(in)    :: Xadj_e
  real(dp),  dimension(ddx_data % nbasis),      intent(inout) :: basloc,  vplm
  real(dp),  dimension(3, ddx_data % nbasis),   intent(inout) :: dbsloc
  real(dp),  dimension(ddx_data % lmax+1),      intent(inout) :: vcos, vsin
  real(dp),  dimension(3),           intent(inout) :: force_e

  ! Local Variables
  ! igrid  : Index of grid points
  ! jsph   : Index of jth sphere
  ! ksph   : Index of kth sphere
  ! ineigh : Row pointer over ith row
  ! l      : Index for l, l:0,...,lmax
  ! m      : Index for m, m:-l0,...,l0
  ! ind    : l*l+l+1
  ! jk     : Row pointer over kth row
  ! NM     : Argument for calling Bessel functions
  integer :: igrid, jsph, ksph, ineigh, l, m, ind, jk,  NM
  ! SI_rjin : Besssel function of first kind for rijn
  ! DI_rjin : Derivative of Besssel function of first kind for rijn
  ! SI_rjkn : Besssel function of first kind for rjkn
  ! DI_rjkn : Derivative of Besssel function of first kind for rjkn
  real(dp), dimension(0:ddx_data % lmax) :: SI_rjin, DI_rjin, SI_rjkn, DI_rjkn, SI_rjin_d,&
                                          & DI_rjin_d, SI_rjin_dh, DI_rjin_dh
  
  logical :: proc
  ! rjin    : r_i*r_1^i(x_j^n) = |x_j^n-x_i|
  ! tji     : r_1^i(x_j^n)
  ! xji     : chi_i(x_j^n)
  ! oji     : omega^\eta_jin
  ! fac     : \delta_fj_n*\omega^\eta_ji
  ! f1      : First factor in alpha computation
  ! f2      : Second factor in alpha computation
  ! beta_ji : Eq.(57) Stamm.etal.18
  ! dj      : Before Eq.(10) Stamm.etal.18
  ! tlow    : Lower bound for switch region
  ! thigh   : Upper bound for switch region
  real(dp)  :: rjin, tji, xji, oji, fac, f1, f2, beta_ji, dj, tlow, thigh, f1_d, f1_dh,&
               & f2_d, f2_dh, arg_d, arg_dh, step_size, f2_taylor, f2_taylor_dh
  ! beta_jk : Eq.(58) Stamm.etal.18
  ! rjkn    : r_k*r_1^k(x_j^n) = |x_j^n-x_k|
  ! tjk     : r_1^k(x_j^n)
  ! xjk     : chi_k(x_j^n)
  real(dp)  :: b, beta_jk, g1, g2, rjkn, tjk, xjk
  ! vji   : x_j^n-x_i
  ! sji   : e^i(x_j^n)
  ! vjk   : x_j^n-x_k
  ! sjk   : e^k(x_j^n)
  ! alpha : Eq.(56) Stamm.etal.18
  ! vb    : Eq.(60) Stamm.etal.18
  ! vc    : Eq.(59) Stamm.etal.18
  real(dp)  :: vji(3), sji(3), vjk(3), sjk(3), alpha(3), vb(3), vc(3)
  ! rho    : Argument for ylmbas
  ! ctheta : Argument for ylmbas
  ! stheta : Argument for ylmbas
  ! cphi   : Argument for ylmbas
  ! sphi   : Argument for ylmbas
  real(dp) :: rho, ctheta, stheta, cphi, sphi, ri
  real(dp), external :: dnrm2
  
  SI_rjin = 0
  DI_rjin = 0
  SI_rjin_d = 0
  DI_rjin_d = 0
  SI_rjin_dh = 0
  DI_rjin_dh = 0
  SI_rjkn = 0
  DI_rjkn = 0
  
  step_size = 0.000001
  arg_d =  0.094112531233527419
  arg_dh = arg_d + step_size

  tlow  = one - pt5*(one - ddx_data % se)*ddx_data % eta
  thigh = one + pt5*(one + ddx_data % se)*ddx_data % eta

  do igrid = 1, ddx_data % ngrid
    vb = zero
    vc = zero
    do ineigh = ddx_data % inl(isph), ddx_data % inl(isph+1) - 1
      jsph = ddx_data % nl(ineigh)
      vji  = ddx_data % csph(:,jsph) + &
              & ddx_data % rsph(jsph)*ddx_data % cgrid(:,igrid) - &
              & ddx_data % csph(:,isph)
      rjin = dnrm2(3, vji, 1)
      tji  = rjin/ddx_data % rsph(isph)
      ri = ddx_data % rsph(isph)

      if (tji.gt.thigh) cycle
      ! Computation of modified spherical Bessel function values      
      call SPHI_bessel(ddx_data % lmax, arg_d, NM, SI_rjin_d, DI_rjin_d)
      call SPHI_bessel(ddx_data % lmax, arg_dh, NM, SI_rjin_dh, DI_rjin_dh)
      
      call SPHI_bessel(ddx_data % lmax, rjin*ddx_data % kappa, NM, SI_rjin, DI_rjin)
      
      sji  = vji/rjin
      !call dbasis(sji,basloc,dbsloc,vplm,vcos,vsin)
      call dbasis(ddx_data, sji, basloc, dbsloc, vplm, vcos, vsin)
      alpha = zero
      do l = 0, ddx_data % lmax
        ind = l*l + l + 1
        if ((SI_rjin(l).lt.tol_zero) .or. (SI_ri(l,isph).lt.tol_zero) &
          & .or. (SI_rjin(l)/SI_ri(l,isph).gt.(rjin/ri)**l)) then
          f1 = (l*(rjin/ri)**(l-1))/ri
          f2 = (rjin/ri)**l
        else if ( SI_ri(l,isph).gt.tol_inf) then
          f1 = zero
          f2 = zero
        else
          f1 = (DI_rjin(l)*ddx_data % kappa)/SI_ri(l,isph)
          f2 = SI_rjin(l)/SI_ri(l,isph)
        end if
          
        f1_d = DI_rjin_d(l)
        f2_d = SI_rjin_d(l)
          
        f1_dh = DI_rjin_dh(l)
        f2_dh = SI_rjin_dh(l)
        
        f2_taylor = l/arg_d + &
            & (l + one)*(ddx_data % kappa**2*arg_d)/((two*l + one) * &
            & (two*l + three))
        f2_taylor_dh = l/arg_dh + &
            & (l + one)*(ddx_data % kappa**2*arg_dh)/((two*l + one) * &
            & (two*l + three))
        do m = -l, l
          alpha = alpha + (f1*sji*basloc(ind+m) + (f2/rjin)*dbsloc(:,ind+m))*Xe(ind+m,isph)
        end do
        write(*,*) ' DEBUG rjin                                 :', rjin*ddx_data % kappa
        write(*,*) ' DEBUG ri                                   :', ri*ddx_data % kappa
        write(*,*) ' DEBUG, Finite Differece derivative, ', l,' : ', (f2_dh-f2_d)/step_size
        write(*,*) ' DEBUG, Taylor Differece derivative, ', l,' : ', (f2_taylor_dh-f2_taylor)/step_size
        write(*,*) ' DEBUG, Analytical Derivative      , ', l,' : ', f1_d
        write(*,*) ' '
      end do
      xji = fsw(tji,ddx_data % se,ddx_data % eta)
      if (ddx_data % fi(igrid,jsph).gt.one) then
        oji = xji/ddx_data % fi(igrid,jsph)
      else
        oji = xji
      end if
      f1 = oji
      vb = vb + f1*alpha*Xadj_e(igrid,jsph)
      if (tji .gt. tlow) then
        ! Compute beta_jin, i.e., Eq.(57) Stamm.etal.18
        beta_ji = compute_beta(ddx_data, SI_rjin, rjin, isph, Xe(:,isph), basloc)
        if (ddx_data % fi(igrid,jsph) .gt. one) then
          dj  = one/ddx_data % fi(igrid,jsph)
          fac = dj*xji
          proc = .false.
          b    = zero
          do jk = ddx_data % inl(jsph), ddx_data % inl(jsph+1) - 1
            ksph = ddx_data % nl(jk)
            vjk  = ddx_data % csph(:,jsph) + &
                 & ddx_data % rsph(jsph)*ddx_data % cgrid(:,igrid) - &
                 & ddx_data % csph(:,ksph)
            rjkn = dnrm2(3, vjk, 1)
            tjk  = rjkn/ddx_data % rsph(ksph)
            ! Computation of modified spherical Bessel function values      
            call SPHI_bessel(ddx_data % lmax,rjkn*ddx_data % kappa,NM,SI_rjkn,DI_rjkn)
            if (ksph.ne.isph) then
              if (tjk .le. thigh) then
              proc = .true.
              sjk  = vjk/rjkn
              !call ylmbas(sjk,basloc,vplm,vcos,vsin)
              call ylmbas(sjk, rho, ctheta, stheta, cphi, sphi, &
                  & ddx_data % lmax, ddx_data % vscales, basloc, vplm, &
                  & vcos, vsin)
              beta_jk  = compute_beta(ddx_data, SI_rjkn, rjkn, ksph, Xe(:,ksph), basloc)
              xjk = fsw(tjk, ddx_data % se, ddx_data % eta)
              b   = b + beta_jk*xjk
              end if
            end if
          end do
          if (proc) then
            g1 = dj*dj*dfsw(tji,ddx_data % se,ddx_data % eta)/ddx_data % rsph(isph)
            g2 = g1*Xadj_e(igrid,jsph)*b
            vc = vc + g2*sji
          end if
        else
          dj  = one
          fac = zero
        end if
        f2 = (one-fac)*dj*dfsw(tji,ddx_data % se,ddx_data % eta)/ddx_data % rsph(isph)
        vb = vb + f2*Xadj_e(igrid,jsph)*beta_ji*sji
      end if 
    end do
    force_e = force_e + ddx_data % wgrid(igrid)*(vb - vc)
    end do
  return
  end subroutine fdokb_b_xe
  
  real(dp) function compute_beta(ddx_data, SI_rijn, rijn, jsph, Xe, basloc)
  implicit none
  type(ddx_type) :: ddx_data
  real(dp), dimension(0:ddx_data % lmax), intent(in) :: SI_rijn
  real(dp), dimension(ddx_data % nbasis), intent(in) :: basloc
  real(dp), dimension(ddx_data % nbasis), intent(in)   :: Xe
  integer, intent(in) :: jsph
  real(dp), intent(in) :: rijn

  integer :: l, m, ind
  real(dp)  :: ss, fac, rj
  ss = zero
  rj = ddx_data % rsph(jsph)

  ! loop over l
  do l = 0, ddx_data % lmax
    ind = l*l + l + 1
    do m = -l, l
      if ((SI_rijn(l).lt.zero) .or. (SI_ri(l,jsph).lt.1.0D-15) &
          & .or. (SI_rijn(l)/SI_ri(l,jsph).gt.(rijn/rj)**l)) then
        fac = (rijn/rj)**l
      else if ( SI_ri(l,jsph).gt.1.0D30) then
        fac = zero
      else
        fac = SI_rijn(l)/SI_ri(l,jsph)
      end if
      ss = ss + fac*basloc(ind+m)*Xe(ind+m)
    end do
  end do
     
  compute_beta = ss
  end function compute_beta
  
  !
  ! Subroutine to compute the derivative of C_1 and C_2
  ! @param[in] ddx_data   : Data Type
  ! @param[in] ksph       : Derivative with respect to x_k
  ! @param[in] Xr         : Solution of the Laplace problem
  ! @param[in] Xe         : Solution of the HSP problem
  ! @param[in] Xadj_r     : Solution of the Adjoint Laplace problem
  ! @param[in] Xadj_e     : Solution of the Adjoint HSP problem
  ! @param[inout] force_r : Force corresponding to Laplace problem
  ! @param[inout] force_e : Force corresponding to HSP problem
  ! Work in progress
  subroutine fdoc(ddx_data, ksph, Xr, Xe, Xadj_r, Xadj_e, Xadj_r_sgrid, &
                 & Xadj_e_sgrid, force_r, force_e)
  implicit none
  type(ddx_type) :: ddx_data
  integer , intent(in) :: ksph
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: Xr, Xe
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: Xadj_e, Xadj_r
  real(dp),  dimension(ddx_data % ngrid), intent(in)    :: Xadj_r_sgrid, Xadj_e_sgrid
  real(dp), dimension(3), intent(inout) :: force_r, force_e
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph) :: diff_re
  real(dp), dimension(nbasis0, ddx_data % nsph) :: diff0
  
  ! The matrix C_1[Xr]+C_2[Xe] has two terms depending on x_i, namely U_i(x_in)
  ! and [Q]_jl'm'^in. We compute these contributions seperately.
  ! First: Compute \nabla^k(U_i(x_in)) [Q]_jl'm'^in 
  call derivative_U(ddx_data, ksph, Xr, Xe, Xadj_r_sgrid, Xadj_e_sgrid, force_r,&
                    & force_e, diff_re, diff0 )
  ! Second: Compute U_i(x_in) [Q^k]_jl'm'^in 
  call derivative_Q(ddx_data, ksph, Xr, Xe, Xadj_r, Xadj_e, force_r, force_e, diff_re, &
                    & diff0 )
  return
  end subroutine fdoc
  
  !
  ! Subroutine to compute the derivative of PChi
  ! @param[in]  ddx_data : Data Type
  ! @param[in]  ksph     : Derivative with respect to x_k
  ! @param[out] dpmat    : P^kChi_i
  subroutine derivative_mkpmat(ddx_data, ksph, dpmat)
  implicit none
  type(ddx_type) :: ddx_data
  integer, intent (in) :: ksph
  real(dp), dimension(3, ddx_data % nbasis, nbasis0,&
                      & ddx_data % nsph), intent(out) :: dpmat
  ! Local Variables
  ! vik     : x_i^n-x_k
  ! sik     : e^k(x_i^n)
  ! f       : w_n*Y_lm(s_n)*Der_U_i(x_i^n)
  ! der_u_i : \nabla U_i(x_i^n)
  real(dp), dimension(3) :: vik, sik, f, der_u_i
  ! rikn  : r_k*r_1^k(x_k^n) = |x_i^n-x_k|
  ! tik   : r_1^k(x_i^n)
  ! fac   : Der_Chi(r_1^k(x_i^n))*e^k(x_i^n)/r_k
  ! f0    : Y_l0m0(s_n)
  real(dp) :: rikn, tik, fac, f0, swthr
  ! isph  : Index for sphere
  ! ig    : Index for grid point
  ! l     : l=0,lmax
  ! m     : m=-l, l
  ! ind   : l*l+l+1
  ! l0    : l=0,lmax0
  ! m0    : m0=-l0, l0
  ! ind0  : l0*l0+l0+1
  integer  :: isph, ig, l, m, ind, l0, m0, ind0, neighbor, ji
  real(dp), external :: dnrm2
  
  swthr = one + (ddx_data % se + one)*ddx_data % eta/two
  neighbor = zero
  
  ! Loop over spheres, i.e., different matrix PChi_i
  do isph = 1, ddx_data % nsph
    neighbor = zero
    do ji = ddx_data % inl(isph), ddx_data % inl(isph+1) - 1
      if (ksph .eq. ddx_data % nl(ji)) then
        neighbor = one
        exit
      end if
    end do
    ! Loop over grid points
    do ig = 1, ddx_data % ngrid
      der_u_i = zero
      if(ddx_data % ui(ig, isph) .gt. zero .and. ddx_data % ui(ig, isph) .lt. one) then
        vik = ddx_data % csph(:, isph) + ddx_data % rsph(isph)*ddx_data % cgrid(:,ig) -&
              & ddx_data % csph(:, ksph)
        rikn = dnrm2(3, vik, 1)
        tik  = rikn/ddx_data % rsph(ksph)
        ! Neighbor of x_i
        if(neighbor .eq. one .and. tik.lt.swthr .and. tik.gt.swthr-ddx_data % eta) then
          sik = vik/rikn
          fac = dfsw(tik, ddx_data % se, ddx_data % eta)/ddx_data % rsph(ksph)
          der_u_i = fac*sik
        ! ksph=isph
        else if (ksph .eq. isph) then
          der_u_i = -ddx_data % zi(:,ig,isph)
        ! Not neighbor or equal sphere
        else
          der_u_i = zero
        end if
        do l = 0, ddx_data % lmax
          ind = l*l + l + 1
          do m = -l,l
            f = ddx_data % wgrid(ig)*ddx_data % vgrid(ind+m, ig)*der_u_i
            do l0 = 0, lmax0
              ind0 = l0*l0 + l0 + 1
              do m0 = -l0, l0
                f0 = ddx_data % vgrid(ind0+m0, ig)
                dpmat(:,ind+m,ind0+m0,isph) = dpmat(:,ind+m,ind0+m0,isph) + f*f0
              end do
            end do
          end do
        end do
      end if
    end do
  end do
  return
  end subroutine derivative_mkpmat
  
  !
  ! Subroutine to compute the derivative of U_i(x_in)
  ! @param[in] ddx_data     : Data Type
  ! @param[in] ksph         : Derivative with respect to x_k
  ! @param[in] Xr           : Solution of the Laplace problem
  ! @param[in] Xe           : Solution of the HSP problem
  ! @param[in] Xadj_r_sgrid : Solution of the Adjoint Laplace problem
  ! @param[in] Xadj_e_sgrid : Solution of the Adjoint HSP problem
  ! @param[inout] force_r   : Force corresponding to Laplace problem
  ! @param[inout] force_e   : Force corresponding to HSP problem
  ! @param[out] diff_re     : epsilon_1/epsilon_2 * l'/r_j[Xr]_jl'm' 
  !                         - (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
  ! @param[out] diff0       : dot_product([PU_j]_l0m0^l'm', l'/r_j[Xr]_jl'm' -
  !                        (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm')
  subroutine derivative_U(ddx_data, ksph, Xr, Xe, Xadj_r_sgrid, Xadj_e_sgrid, force_r, & 
                          &force_e, diff_re, diff0)
  implicit none
  type(ddx_type), intent(in) :: ddx_data
  integer, intent(in) :: ksph
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: Xr, Xe
  real(dp), dimension(ddx_data % ngrid, ddx_data % nsph), intent(in) :: Xadj_r_sgrid,&
                                                                        & Xadj_e_sgrid
  real(dp), dimension(3), intent(inout) :: force_r, force_e
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(out) :: diff_re
  real(dp), dimension(nbasis0, ddx_data % nsph), intent(out) :: diff0
  real(dp), external :: dnrm2
  ! Local variable
  ! isph  : Index for sphere i
  ! jsph  : Index for sphere j
  ! igrid : Index for grid point
  ! l     : l=0, lmax
  ! m     : -l,..,l
  ! ind   : l^2+l+1
  ! l0    : l0=0, lmax0
  ! m0    : -l0, l0
  ! ind0  : l0^2+l0+1
  ! kep   : Index for external grid point
  integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, kep
  ! vik     : x_i^n-x_k
  ! der_u_i : \nabla U_i(x_i^n)
  ! sik     : e^k(x_i^n)
  real(dp), dimension(3) :: vik, der_u_i, sik
  ! val   : Intermediate variable to compute diff_ep
  real(dp) :: val
  ! coefvec : \omega_n*\nabla (U_i^e(x_in))*Y_lm(s_n)
  real(dp), dimension(3, ddx_data % ngrid, ddx_data % nbasis, ddx_data % nsph) :: coefvec
  ! phi_in : sum_{j=1}^N diff0_j * coefY_j
  real(dp), dimension(ddx_data % ngrid, ddx_data % nsph) :: phi_in
  ! Debug purpose
  real(dp), dimension(ddx_data % ncav, nbasis0, ddx_data % nsph) :: coefY_unit
  real(dp), dimension(ddx_data % nbasis, nbasis0, ddx_data % nsph):: Pchi_d
  
  ! For checking the derivative of U_i^e, we make the matrix Q constant
  Pchi_d = one
  coefY_unit = one

  ! Compute l'/r_j[Xr]_jl'm' -(i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
  do jsph = 1, ddx_data % nsph
    do l = 0, ddx_data % lmax
      do m = -l,l
        ind = l**2 + l + m + 1
        diff_re(ind,jsph) = ((epsp/ddx_data % eps)*(l/ddx_data % rsph(jsph)) * &
              & Xr(ind,jsph) - termimat(l,jsph)*Xe(ind,jsph))
      end do
    end do
  end do
  
  ! diff0 = Pchi * diff_re, linear scaling
  diff0 = zero 
  do jsph = 1, ddx_data % nsph
    do ind0 = 1, nbasis0
      diff0(ind0, jsph) = dot_product(diff_re(:,jsph), &
          & Pchi_d(:,ind0, jsph))
    end do
  end do
  
  ! phi_in = diff0 * coefY,    COST: M^2*nbasis*Nleb
  ! Here, summation over j takes place
  phi_in = zero
  kep = 0
  do igrid = 1, ddx_data % ngrid
    do isph = 1, ddx_data % nsph
      if(ddx_data % ui(igrid, isph) .gt. zero) then
        kep = kep + 1
        val = zero
        do jsph = 1, ddx_data % nsph 
          do ind0 = 1, nbasis0
            val = val + diff0(ind0,jsph)*coefY_unit(kep,ind0,jsph)
          end do
        end do
      end if
    phi_in(igrid, isph) = val
    end do
  end do
  call fdoga(ddx_data, ksph, Xadj_r_sgrid, phi_in, force_r)
  call fdoga(ddx_data, ksph, Xadj_e_sgrid, phi_in, force_e)
  
  end subroutine derivative_U
  
  !
  ! Subroutine to compute the derivative of [Q]_jl'm'^in
  ! @param[in] ddx_data   : Data Type
  ! @param[in] ksph       : Derivative with respect to x_k
  ! @param[in] Xr         : Solution of the Laplace problem
  ! @param[in] Xe         : Solution of the HSP problem
  ! @param[in] Xadj_r     : Solution of the Adjoint Laplace problem
  ! @param[in] Xadj_e     : Solution of the Adjoint HSP problem
  ! @param[inout] force_r : Force corresponding to Laplace problem
  ! @param[inout] force_e : Force corresponding to HSP problem
  ! @param[in] diff_re    : l'/r_j[Xr]_jl'm' -(i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
  ! @param[in] diff0      : dot_product([PU_j]_l0m0^l'm', l'/r_j[Xr]_jl'm' -
  !                        (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm')
  subroutine derivative_Q(ddx_data, ksph, Xr, Xe, Xadj_r, Xadj_e, force_r, force_e, & 
                         & diff_re, diff0)
  use bessel
  implicit none
  type(ddx_type), intent(in) :: ddx_data
  integer, intent(in) :: ksph
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: Xr, Xe
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: Xadj_r, Xadj_e
  real(dp), dimension(3), intent(inout) :: force_r, force_e
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph), intent(in) :: diff_re
  real(dp), dimension(nbasis0, ddx_data % nsph), intent(in) :: diff0
  real(dp), external :: dnrm2
  ! Local variable
  ! isph  : Index for sphere i
  ! jsph  : Index for sphere j
  ! igrid : Index for grid point
  ! l     : l=0, lmax
  ! m     : -l,..,l
  ! ind   : l^2+l+1
  ! l0    : l0=0, lmax0
  ! m0    : -l0, l0
  ! ind0  : l0^2+l0+1
  ! kep   : Index for external grid point
  ! i     : Index for dimension
  ! NM     : Argument for calling Bessel functions
  integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, kep, i, NM
  ! vik      : x_i^n-x_k
  ! der_u_i  : \nabla U_i(x_i^n)
  ! sik      : e^k(x_i^n)
  ! val_dim3 : Intermediate variable to compute diff_ep_der_Pchi
  real(dp), dimension(3) :: vik, der_u_i, sik, val_dim3
  ! tik   : r_1^k(x_i^n)
  ! tij   : r_1^j(x_i^n)
  ! fac   : Der_Chi(r_1^k(x_i^n))*e^k(x_i^n)/r_k
  ! rikn  : r_k*r_1^k(x_k^n) = |x_i^n-x_k|
  ! rijn  : r_j*r_1^j(x_j^n) = |x_i^n-x_j|
  ! f1    : First factor in alpha computation
  ! f2    : Second factor in alpha computation
  real(dp) :: tik, tij, fac, rikn, rijn, f1, f2
  ! vij   : x_i^n-x_j
  ! sij   : e^j(x_i^n)
  real(dp)  :: vij(3), sij(3)
  ! basloc : Y_lm(s_n)
  ! vplm   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % nbasis):: basloc, vplm
  ! dbasloc : Derivative of Y_lm(s_n)
  real(dp),  dimension(3, ddx_data % nbasis):: dbasloc
  ! vcos   : Argument to call ylmbas
  ! vsin   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % lmax+1):: vcos, vsin
  ! der_c_one_c_two : Derivative of C_1 and C_2
  real(dp), dimension(3, ddx_data % nbasis, ddx_data % nsph) :: der_c_one_c_two
  ! coefvec : \omega_n*U_i^e(x_in)*Y_lm(s_n)
  real(dp), dimension(ddx_data % ngrid, ddx_data % nbasis, ddx_data % nsph) :: coefvec
  ! diff0_der : Derivative of PChi matrix scaled by diff_re
  real(dp), dimension(3, nbasis0, ddx_data % nsph) :: diff0_der
  ! dpmat   : Derivative of PChi matrix. Size: 3 x nbasis x nbasis0 x nsph
  real(dp), dimension(3, ddx_data % nbasis, nbasis0, &
                      & ddx_data % nsph) :: dpmat
  ! diff_ep_der_Pchi : First product rule term
  real(dp), dimension(3, ddx_data % ncav) :: diff_ep_der_Pchi
  ! diff_ep_der_ki_Y : Second and Third product rule term
  real(dp), dimension(3) :: diff_ep_der_ki_Y
  ! SK_rijn : Besssel function of first kind for rijn
  ! DK_rijn : Derivative of Besssel function of first kind for rijn
  real(dp), dimension(0:ddx_data % lmax) :: SK_rijn, DK_rijn
  
  SK_rijn = 0
  DK_rijn = 0
  
  ! (?) This can be taken from update_rhs
  do isph = 1, ddx_data % nsph
    ! Loop over grid points
    do igrid = 1, ddx_data % ngrid
      ! Check if the characteristic function is positive
      if (ddx_data % ui(igrid, isph) .gt. zero) then
        do ind  = 1, ddx_data % nbasis 
          coefvec(igrid,ind,isph) = ddx_data % wgrid(igrid) &
                                  &*ddx_data % ui(igrid,isph)*ddx_data % vgrid(ind,igrid)
        end do
      end if
    end do
  end do
  
  ! Derivative of PChi_i with respect to x_k
  call derivative_mkpmat(ddx_data, ksph, dpmat)
  
  ! diff0_der = [PChi^k] * diff_re
  diff0_der = zero
  do jsph = 1, ddx_data % nsph
    do ind0 = 1, nbasis0
      do i = 1, 3
        diff0_der(i, ind0, jsph) = dot_product(diff_re(:,jsph), &
          & dpmat(i, :,ind0, jsph))
      end do
    end do
  end do
  
  ! diff_ep_der_Pchi = diff0_der * coefY,    COST: M^2*nbasis*Nleb
  diff_ep_der_Pchi = zero
  do kep = 1, ddx_data % ncav
    val_dim3 = zero
    do jsph = 1, ddx_data % nsph 
      do ind0 = 1, nbasis0
        val_dim3 = val_dim3 + diff0_der(:,ind0,jsph)*coefY(kep,ind0,jsph)
      end do
    end do
    diff_ep_der_Pchi(:, kep) = val_dim3(:)
  end do
  
  ! diff_ep_der_ki_Y = \sum_j \sum_l0m0 (C_ik [PU_j]_l0m0^l'm' (\nabla^k(k_l0^j(x_ni)*
  !                     Y_l0m0^j(x_ni) +k_l0^j(x_ni)\nabla^k(Y_l0m0^j(x_ni)))))
  diff_ep_der_ki_Y = zero
  der_c_one_c_two = zero
  kep = 0
  ! Loop over i-th sphers
  do isph = 1, ddx_data % nsph
    ! Loop over gird points
    do igrid = 1, ddx_data % ngrid
      ! Loop over j-th sphere
      do jsph = 1, ddx_data % nsph
        vij = ddx_data % csph(:, isph) + ddx_data % rsph(isph)*ddx_data % cgrid (:, igrid)&
              & - ddx_data % csph(:, jsph)
        rijn = dnrm2(3, vij, 1)
        tij = rijn/ddx_data % rsph(jsph)
        sij = vij / rijn
        ! Compute Y_lm(x_ni) and \nabla^k (Y_lm(x_ni))
        call dbasis(ddx_data, sij, basloc, dbasloc, vplm, vcos, vsin)
        ! Loop over lmax0
        do l0 = 0, lmax0
          ! Compute k_l0^ĵ(x_ni) and \nabla^k(k_l0^ĵ(x_ni))
          call SPHK_bessel(lmax0, rijn*ddx_data % kappa, NM, SK_rijn, DK_rijn)
          ind0 = l0*l0 + l0 + 1
          ! Loop over m0
          do m0 = -l0, l0
            ! f1 = k_l0'(x_ni-x_j)/k_l0(r_j) Y_l0m0(x_ni)
            f1 = (DK_rijn(l0)/SK_ri(l0, jsph)*ddx_data % kappa)*basloc(ind0 + m0)
            ! f2 = k_l0(x_in)/r_jr_1^j(x_ni)
            f2 = (SK_rijn(l0)/rijn*SK_ri(l0, jsph))
            if(ksph .eq. jsph) then
              diff_ep_der_ki_Y = diff_ep_der_ki_Y - (f1*sij + f2*dbasloc(:, ind0+m0))&
                                 & * C_ik(l0, jsph)*diff0(ind0, jsph)
            else if (ksph .eq. isph) then
              diff_ep_der_ki_Y = diff_ep_der_ki_Y + (f1*sij + f2*dbasloc(:, ind0+m0))&
                                 & * C_ik(l0, jsph)*diff0(ind0, jsph)
            else
              diff_ep_der_ki_Y = zero
            end if
          end do
        end do
        if (ddx_data % ui(igrid, jsph) .gt. zero) then
          kep = kep + 1
          ! Loop over lm
          do ind = 1, ddx_data % nbasis
            der_c_one_c_two(:, ind, jsph) = der_c_one_c_two (:, ind, jsph) &
                                           & + (diff_ep_der_Pchi(:, kep) &
                                           & + diff_ep_der_ki_Y(:)) &
                                           & * coefvec(igrid, ind, jsph)
          end do
        end if
      end do
    end do
    do ind = 1, ddx_data % nbasis
      force_r = force_r + Xadj_r(ind, isph)*der_c_one_c_two(:, ind, isph)
      force_e = force_e + Xadj_e(ind, isph)*der_c_one_c_two(:, ind, isph)
    end do
  end do
  end subroutine derivative_Q
  
  !
  ! Debug derivative of C1 and C2 in code
  ! @param[in, out] rhs_cosmo : C_1*X_r + C_2*X_e
  ! @param[in, out] rhs_hsp   : C_1*X_r + C_2*X_e
  ! @param[in] Xr             : X_r
  ! @param[in] Xe             : X_e
  !
  subroutine C1_C2(ddx_data, rhs_cosmo, rhs_hsp, Xr, Xe)
  use bessel
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph), intent(inout) :: rhs_cosmo, rhs_hsp
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph) :: rhs_plus
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph), intent(in) :: Xr, Xe
  
  integer :: isph, jsph, ig, kep, ind, l1,m1, ind1, ind0, count, istatus
  real(dp), dimension(3) :: vij
  real(dp), dimension(ddx_data % nbasis,ddx_data % nsph) :: diff_re
  real(dp), dimension(nbasis0,ddx_data % nsph) :: diff0
  real(dp), dimension(ddx_data % nbasis, nbasis0, ddx_data % nsph):: Pchi_d
  real(dp), dimension(ddx_data % ngrid, ddx_data % nbasis, ddx_data % nsph) ::coefvec_d
  real(dp), dimension(0:ddx_data % lmax, ddx_data % nsph):: termimat_d
  real(dp), dimension(ddx_data % ncav) :: diff_ep
  real(dp) :: val
  real(dp), dimension(ddx_data % ncav, nbasis0, ddx_data % nsph) :: coefY_unit
      
  ! Compute P_chi matrix, Eq.(87)
  !do jsph = 1, ddx_data % nsph
  !  call mkpmat(ddx_data, jsph, Pchi_d(:,:,jsph))
  !end do

  ! For checking the derivative of U_i^e, we make the matrix Q constant
  Pchi_d = one
  coefY_unit = one
      
  ! Precompute coefvec of Qmat, Cost: linear scaling
  do isph = 1,ddx_data % nsph
    do ig = 1,ddx_data % ngrid
      if (ddx_data % ui(ig, isph) .gt. 0) then
        do ind  = 1, ddx_data % nbasis 
          coefvec_d(ig,ind,isph) = ddx_data % wgrid(ig)*ddx_data % ui(ig,isph)*ddx_data % vgrid(ind,ig)
        end do
      end if
    end do
  end do
      
  ! precompute termimat
  do jsph = 1, ddx_data % nsph
    do l1 = 0, ddx_data % lmax
      if (max(DI_ri(l1,jsph),SI_ri(l1,jsph)).gt.tol_inf) then
        termimat_d(l1,jsph) = ddx_data % kappa
      else if (min(DI_ri(l1,jsph),SI_ri(l1,jsph)).lt.tol_zero) then
        termimat_d(l1,jsph) = l1/ddx_data % rsph(jsph) + &
            & (l1 + 1)*(ddx_data % kappa**2*ddx_data % rsph(jsph))/((two*l1 + one) * &
            & (two*l1 + three))
      else
        termimat_d(l1,jsph) = DI_ri(l1,jsph)/SI_ri(l1,jsph)*ddx_data % kappa
      end if
    end do
  end do

  ! diff_re = epsp/eps*l1/ri*Xr - i'(ri)/i(ri)*Xe,
  diff_re = zero 
  do jsph = 1, ddx_data % nsph
    do l1 = 0, ddx_data % lmax
      do m1 = -l1,l1
        ind1 = l1**2 + l1 + m1 + 1
        diff_re(ind1,jsph) = ((epsp/ddx_data % eps)*(l1/ddx_data % rsph(jsph)) * &
            & Xr(ind1,jsph) - termimat_d(l1,jsph)*Xe(ind1,jsph))
      end do
    end do
  end do

  ! diff0 = Pchi * diff_er, linear scaling
  ! TODO: probably doing PX on the fly is better 
  diff0 = zero 
  do jsph = 1, ddx_data % nsph
    do ind0 = 1, nbasis0
      diff0(ind0, jsph) = dot_product(diff_re(:,jsph), &
          & Pchi_d(:,ind0, jsph))
    end do
  end do

  ! diff_ep = diff0 * coefY,    COST: M^2*nbasis*Nleb
  ! Changing coefY to coefY_unit
  diff_ep = zero
  do kep = 1, ddx_data % ncav
    val = zero
    do jsph = 1, ddx_data % nsph 
      do ind0 = 1, nbasis0
        val = val + diff0(ind0,jsph)*coefY_unit(kep,ind0,jsph)
      end do
    end do
    diff_ep(kep) = val 
  end do

  rhs_plus = zero
  kep = 0
  do isph = 1, ddx_data % nsph
    do ig = 1, ddx_data % ngrid
      if (ddx_data % ui(ig,isph).gt.zero) then 
        kep = kep + 1
        do ind = 1, ddx_data % nbasis
          rhs_plus(ind,isph) = rhs_plus(ind,isph) + &
              & coefvec_d(ig,ind,isph)*diff_ep(kep)
        end do
      end if
    end do
  end do

  rhs_cosmo = rhs_plus
  rhs_hsp =  rhs_plus
  return
  end subroutine C1_C2  

  
  !
  ! Debug routine to check the adjoint system. Compute <y^T,Lx>
  ! @param[in] ddx_data : Data type
  ! @param[in] lx       : Matrix multiplication L*x
  ! @param[in] lstarx   : Matrix multiplication Lstar*x
  subroutine debug_adjoint(ddx_data, lx, lstarx)
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  external  :: lx, lstarx
  ! Local variable
  ! vector_y      : Lx
  ! vector_x      : Initial vector
  ! vector_y_star : Lstarx
  ! unit_vector   : Unit vector
  real(dp), dimension(ddx_data % n):: vector_y, vector_x, vector_y_star
  real(dp), dimension(ddx_data % n):: unit_vector
  ! in          : Index for n
  ! sum_adjoint  : Sum for adjoint system
  ! sum_original : Sum for original system
  integer :: in
  real(dp) :: sum_adjoint, sum_original
  
  vector_x = one
  vector_y = one
  vector_y_star = one
  unit_vector = one
  ! Call original system
  if (ddx_data % model .eq. 3) then 
    call lx(ddx_data, ddx_data % n, vector_x, vector_y)
  else
    call lx(ddx_data, vector_x, vector_y)
  end if
  ! Call adjoint system
  call lstarx(ddx_data, vector_x, vector_y_star)
  write(*,*) vector_y
  write(*,*) vector_y_star
  
  sum_adjoint = 0
  sum_original = 0
  do in = 1, ddx_data % n
    sum_adjoint = sum_adjoint + unit_vector(in)*vector_y_star(in)
    sum_original = sum_original + unit_vector(in)*vector_y(in)
  end do
  if (ddx_data % iprint .ge. 5) then
    write(*,*) 'The original system <y^T,Lx>: ', sum_original
    write(*,*) 'The adjoint system <y^T,Lstarx>: ', sum_adjoint
  end if
  return
  end subroutine debug_adjoint
  
    !
  ! Debug routine to check the derivative of system. Compute <y^T,L^Kx>
  ! @param[in] ddx_data : Data type
  ! @param[in] lx       : Matrix multiplication L*x
  subroutine debug_derivative(ddx_data)
  implicit none
  type(ddx_type), intent(inout)  :: ddx_data
  ! Local variable
  ! vector_y_cosmo  : Ax
  ! vector_x        : Initial vector
  ! unit_vector_1   : Unit vector dim: 1 \times M(l_max+1)^2
  ! vector_yh_cosmo : A(x+h)
  real(dp), dimension(ddx_data % n):: vector_y_cosmo, vector_x,&
                                      &unit_vector_1, vector_yh_cosmo
  real(dp), dimension(ddx_data % n):: vector_y_lpb, vector_yh_lpb
  ! unit_vector   : Unit vector dim: nbasis \times nsph
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph):: unit_vector
  ! 
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph):: vector_r, vector_rh
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph):: vector_e, vector_eh
  ! unit_vector_grid : Unit vector evaluated on grid points
  real(dp), dimension(ddx_data % ngrid, ddx_data % nsph):: unit_vector_grid
  ! derivative_cosmo : A^k(x)
  ! derivative_lpb   : A^k(x) 
  real(dp), dimension(3, ddx_data % nsph) :: derivative_cosmo
  real(dp), dimension(3, ddx_data % nsph) :: derivative_lpb
  real(dp), dimension(3, ddx_data % nsph) :: derivative_u_r
  real(dp), dimension(3, ddx_data % nsph) :: derivative_u_e
  ! vsin, vcos, vplm : Values used in basloc
  real(dp), dimension(ddx_data % lmax +1) :: vsin, vcos
  ! basloc : Y_lm
  real(dp), dimension(ddx_data % nbasis)  :: vplm, basloc
  ! dbasloc : Derivatives of Y_lm
  real(dp), dimension(3, ddx_data % nbasis)  :: dbsloc
  
  real(dp), dimension(ddx_data % nbasis, ddx_data % nsph) :: diff_re
  real(dp), dimension(nbasis0, ddx_data % nsph) :: diff0
  ! in          : Index for n
  ! isph        : Index for spheres
  integer :: in, isph, i, jsph, ibasis
  ! sum  : <y^T,Lx>
  ! step : Finite difference step length
  real(dp) :: sum_x_cosmo, sum_xh_cosmo, step
  real(dp) :: sum_x_lpb, sum_xh_lpb
  real(dp) :: sum_x_c1, sum_xh_c1
  
  ! Initial values
  ! x    : one
  ! y=Lx : one
  vector_x = one
  vector_y_cosmo = one
  vector_yh_cosmo = one
  vector_y_lpb = one
  vector_yh_lpb = one
  
  vector_r = zero
  vector_rh = zero
  vector_e = zero
  vector_eh = zero
  
  derivative_lpb = zero
  derivative_cosmo = zero
  derivative_u_e = zero
  derivative_u_r = zero
  
  unit_vector = one
  unit_vector_1 = one
  step = 0.000001
  
  isph = one
  i = one
  
  ! Call Lx for x
  ! Calling for Lx corrsponding to ddCOSMO and ddLPB different
  call lx(ddx_data, vector_x, vector_y_cosmo)
  ! ddLPB
  call matABx(ddx_data , ddx_data % n, vector_x, vector_y_lpb)
  ! C_1 and C_2
  call C1_C2(ddx_data, vector_r, vector_e, unit_vector, unit_vector)
  
  ! Values at x+h, change coordinates for x_i
  ddx_data % csph(i, isph) = ddx_data % csph(i, isph) + step
!   ! Call Lx for x+h
  call lx(ddx_data, vector_x, vector_yh_cosmo)
  ! ddLPB
  call matABx(ddx_data , ddx_data % n, vector_x, vector_yh_lpb)
  ! C_1 and C_2
  
  ! Values at x
  ddx_data % csph(i, isph) = ddx_data % csph(i, isph) - step
  ! Compute solution on grid points
  call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % vgrid, ddx_data % vgrid_nbasis, &
            & unit_vector , ddx_data % nbasis, zero, unit_vector_grid, &
            & ddx_data % ngrid)
  
  ! Compute unit_vector^T*B^k*unit_vector for x_1
  call fdoka(ddx_data, isph, unit_vector, unit_vector_grid(:, isph), &
                  & basloc, dbsloc, vplm, vcos, vsin, derivative_cosmo(:,isph))
  call fdokb(ddx_data, isph, unit_vector, unit_vector_grid, &
                  & basloc, dbsloc, vplm, vcos, vsin, derivative_cosmo(:, isph))
  call fdoka_b_xe(ddx_data, isph, unit_vector, unit_vector_grid(:, isph), &
                  & basloc, dbsloc, vplm, vcos, vsin, derivative_lpb(:,isph))
  call fdokb_b_xe(ddx_data, isph, unit_vector, unit_vector_grid, &
                  & basloc, dbsloc, vplm, vcos, vsin, derivative_lpb(:, isph))
  call derivative_U(ddx_data, isph, unit_vector, unit_vector, unit_vector_grid,&
                  & unit_vector(:, isph), derivative_u_r(:, isph), & 
                  & derivative_u_e(:, isph), diff_re, diff0)
  !call derivative_Q(ddx_data, isph, unit_vector, unit_vector, unit_vector, unit_vector, &
  !                 & derivative_u_r(:, isph), derivative_u_e(:, isph), diff_re, &
  !                & diff0)
  
  sum_x_cosmo = 0
  sum_xh_cosmo = 0
  sum_x_lpb = 0
  sum_xh_lpb = 0
  sum_x_c1 = 0
  
  do in = 1, ddx_data % n
    sum_x_cosmo = sum_x_cosmo + unit_vector_1(in)*vector_y_cosmo(in)
    sum_xh_cosmo = sum_xh_cosmo + unit_vector_1(in)*vector_yh_cosmo(in)
    sum_x_lpb = sum_x_lpb + unit_vector_1(in)*vector_y_lpb(in)
    sum_xh_lpb = sum_xh_lpb + unit_vector_1(in)*vector_yh_lpb(in)
  end do
  
  do ibasis = 1, ddx_data % nbasis
    do jsph = 1, ddx_data % nsph
      sum_x_c1 = sum_x_c1 + unit_vector(ibasis, jsph)*vector_r(ibasis, jsph)  + &
               & unit_vector(ibasis, jsph)*vector_e(ibasis, jsph)
    end do
  end do
  if (ddx_data % iprint .ge. 0) then
    ! Print x-direction force
    write(*,*) '<y^T,A^Kx>           : ', derivative_cosmo(: ,isph)
    ! Print finite difference
    write(*,*) 'Finite Difference  A : ', (sum_xh_cosmo-sum_x_cosmo)/step
    ! Print x-direction force
    write(*,*) '<y^T,B^Kx>           : ', derivative_lpb(: ,isph)
    ! Print finite difference
    write(*,*) 'Finite Difference  B : ', (sum_xh_lpb - sum_x_lpb)/step
    ! Print x-direction force for C_1Xr+C_2Xe, U derivative
    write(*,*) '<y^T, C^Kx>          : ', derivative_u_r(:,isph) + derivative_u_e(:,isph)
    write(*,*) 'Finite Difference  C1: ', sum_x_c1
  end if
  return
  end subroutine debug_derivative

end module ddx_lpb
