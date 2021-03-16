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
integer, parameter :: lmax0 = 6
integer, parameter :: nbasis0 = 49
real(dp),  parameter :: epsp = 1.0d0
!!
!! Taken from Chaoyu's MATLAB code
!!
!! wij          : (?) Maybe Eq. (54) 
real(dp), allocatable :: wij(:,:)
real(dp), allocatable :: coefvec(:,:,:), Pchi(:,:,:), &
                              & Qmat(:,:,:), coefY(:,:,:)
!! SI_ri        : Bessel' function of first kind
!! DI_ri        : Derivative of Bessel' function of first kind
!! SK_ri        : Bessel' function of second kind
!! DK_ri        : Derivative of Bessel' function of second kind
!! NOTE: Tolerance and n_iter_gmres can be given by the user
!! tol_gmres    : Tolerance of GMRES iteration
!! n_iter_gmres : Maximum number of GMRES itertation

real(dp), allocatable :: SI_ri(:,:), DI_ri(:,:), SK_ri(:,:), &
                              & DK_ri(:,:), termimat(:,:)
real(dp)              :: tol_gmres, n_iter_gmres

contains
  !!
  !! ddLPB calculation happens here
  !! @param[in] ddx_data : dd Data 
  !! @param[in] phi     : Boundary conditions
  !! @param[in] psi     : Electrostatic potential vector.
  !!                      Use of psi unknown
  !! @param[in] gradphi : Gradient of phi
  !!
  !! @param[out] sigma  : Solution of ddLPB
  !! @param[out] esolv  : Electrostatic solvation energy
  !!
  subroutine ddlpb(ddx_data, phi, psi, gradphi, sigma, esolv, charge, ndiis, niter, iconv)
  ! main ddLPB
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  logical                         :: converged = .false.
  integer                         :: iteration = 1
  integer, intent(in)             :: ndiis
  integer, intent(in)             :: niter
  integer, intent(in)             :: iconv
  real(dp), intent(inout)    :: esolv
  real(dp)                   :: inc, old_esolv
  real(dp), intent(inout)    :: sigma(ddx_data % nbasis,ddx_data % nsph)
  real(dp), intent(in)       :: phi(ddx_data % ncav), &
                                        & gradphi(3,ddx_data % ncav)
  real(dp), intent(in)       :: psi(ddx_data % nbasis, ddx_data % nsph)
  real(dp), intent(in)       :: charge(ddx_data % nsph)
  !!
  !! Xr         : Reaction potential solution (Laplace equation)
  !! Xe         : Extended potential solution (HSP equation)
  !! rhs_r      : Right hand side corresponding to Laplace equation
  !! rhs_e      : Right hand side corresponding to HSP equation
  !! rhs_r_init : Initial right hand side corresponding to Laplace equation
  !! rhs_e_init : Initial right hand side corresponding to HSP equation
  !!
  real(dp), allocatable ::   Xr(:,:), Xe(:,:), rhs_r(:,:), rhs_e(:,:), &
                                  & rhs_r_init(:,:), rhs_e_init(:,:)
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
  !! tol    : Tolerance for Jacobi solver
  !! n_iter : Number of iterative steps
  real(dp), allocatable :: g(:,:), f(:,:), g0(:), f0(:), phi_grid(:, :)
  integer                    :: isph
  integer                    :: i
  logical                    :: ok = .false.
  real(dp)              :: tol
  integer                    :: n_iter
  integer                    :: its
  !
  ! Allocate Bessel's functions of the first kind and the second kind
  ! and their derivatives
  !
  call ddlpb_init(ddx_data)

  allocate(g(ddx_data % ngrid,ddx_data % nsph),&
           & f(ddx_data % ngrid, ddx_data % nsph), &
           & phi_grid(ddx_data % ngrid, ddx_data % nsph))
  allocate(g0(ddx_data % nbasis),f0(ddx_data % nbasis))
  allocate(rhs_r(ddx_data % nbasis, ddx_data % nsph),&
           & rhs_e(ddx_data % nbasis, ddx_data % nsph))
  allocate(rhs_r_init(ddx_data % nbasis, ddx_data % nsph),&
           & rhs_e_init(ddx_data % nbasis, ddx_data % nsph))
  allocate(Xr(ddx_data % nbasis, ddx_data % nsph),&
           & Xe(ddx_data % nbasis, ddx_data % nsph))

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
    !! intrhs is a subroutine in dd_operators
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

  tol = 10.0d0**(-iconv)

  first_out_iter = .true.

  do while (.not.converged)

    !! Solve the ddCOSMO step
    !! A X_r = RHS_r (= G_X+G_0) 
    !! NOTE: Number of iterative steps can be user provided
    n_iter = 200
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
    call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, ndiis, 4, tol, &
                     & rhs_r, Xr, n_iter, ok, lx, ldm1x, hnorm)
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
      esolv = esolv + pt5*charge(isph)*Xr(1,isph)*(one/(two*sqrt(pi)))
    end do

    !! Check for convergence
    inc = abs(esolv - old_esolv)/abs(esolv)
    old_esolv = esolv
    if ((iteration.gt.1) .and. (inc.lt.tol)) then
      write(6,*) 'Reach tolerance.'
      converged = .true.
    end if
    write(6,*) iteration, esolv, inc
    iteration = iteration + 1

    ! to be removed
    first_out_iter = .false.
  end do

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

  allocate (coefY(ddx_data % ncav, nbasis0, ddx_data % nsph), stat = istatus)
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
      !write(*,*) DI_ri(l0,jsph), SI_ri(l0,jsph), coef_bessel(l0,jsph)
      !write(*,*) (min(-DK_ri(l0,jsph), SK_ri(l0,jsph)).lt.tol_zero), &
      !    & DK_ri(l0,jsph), termk
    end do
  end do

  ! Computation of F0 using above terms
  ! (?) Use of kep not known (?)
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
              ! NOTE: (?) Don't know the use of coefY till now (?)
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
            & vcos(ddx_data % lmax+1), vsin(ddx_data % lmax+1) , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'Bx: allocation failed !'
    stop
  endif

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

  deallocate( pot, basloc, vplm, vcos, vsin , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'matABx: allocation failed !'
    stop
  endif
  end subroutine matABx

  !!
  !! Scale the ddCOSMO solution vector
  !! @param[in]      direction : Direction of the scaling
  !! (?) NOTE: We send direction = 1, then why the case of -1 (?)
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
  call gmresr(ddx_data, .false., ddx_data % nsph*ddx_data % nbasis, gmj, gmm, rhs, Xe, work, tol_gmres,'rel', &
      & n_iter_gmres, r_norm, matABx, info)

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
      if ((SI_rijn(l).lt.zero) .or. (SI_ri(l,isph).lt.tol_zero) &
          & .or. (SI_rijn(l)/SI_ri(l,isph).gt.(rijn/ri)**l)) then
        fac_hsp(ind) = (rijn/ri)**l*basloc(ind)
      else if ( SI_ri(l,isph).gt.tol_inf) then
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
  ! @param[in, out] rhs_cosmo : -C_1*X_r^(k-1) - -C_2*X_e^(k-1) + G_0 + F_0
  ! @param[in, out] rhs_hsp   : -C_1*X_r^(k-1) - -C_2*X_e^(k-1) + F_0
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
  ! @param[out] pmat : Matrix of size nylm X (lmax0+1)^2, Fixed lmax0
  !
  subroutine mkpmat(ddx_data, isph, pmat )
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

end module ddx_lpb
