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
!! Logical variables for iterations
!!
!! first_outer_iter : Logical variable to check if the first outer iteration has been
!!                    performed
!! epsp             : Dielectric permittivity of the solvent. 1.0d0 is for H20
logical :: first_out_iter
real(dp),  parameter :: epsp = 1.0d0
!!
!! Local variables and their definitions that will be used throughout this file
!! isph    : Index for the sphere i
!! jsph    : Index for the sphere j
!! ksph    : Index for the sphere k
!! igrid   : Index for the grid points
!! icav    : Index for the external grid points
!! ibasis  : Index for the basis
!! ibasis0 : Index for the fixed basis (nbasis0)
!! l       : Index for lmax
!! m       : Index for m:-l,...,l
!! ind     : l^2 + l + m + 1
!! l0      : Index for lmax0
!! m0      : Index for m0:-l0,...,l0
!! ind0    : l0^2 + l0 + m0 + 1
!! ineigh  : Index over Row space of neighbors
!! rijn    : r_j*r_1^j(x_i^n) = |x_i^n-x_j|
!! tij     : r_1^j(x_i^n)
!! xij     : chi_j(x_i^n)
!! oij     : omega^\eta_ijn
!! vij     : x_i^n-x_j
!! sij     : e^j(x_i^n)
!! SI_rijn : Besssel function of first kind for rijn
!! DI_rijn : Derivative of Bessel function of first kind for rijn
!! tlow    : Lower bound for switch region
!! thigh   : Upper bound for switch region
!! basloc  : Y_lm(s_n)
!! dbasloc : Derivative of Y_lm(s_n)
!! i       : Index for dimension x,y,z

contains
  !!
  !! ddLPB calculation happens here
  !! @param[in] ddx_data : dd Data 
  !! @param[in] phi      : Boundary conditions
  !! @param[in] psi      : Electrostatic potential vector.
  !! @param[in] gradphi  : Gradient of phi
  !! @param[in] hessianphi  : Hessian of phi
  !! @param[out] esolv   : Electrostatic solvation energy
  !!
  subroutine ddlpb(ddx_data, phi, gradphi, hessianphi, psi, tol, esolv, &
          & force, info)
  ! main ddLPB
  ! Inputs
  type(ddx_type), intent(inout) :: ddx_data
  real(dp), dimension(ddx_data % constants % ncav), intent(in)       :: phi
  real(dp), dimension(3, ddx_data % constants % ncav), intent(in)    :: gradphi
  real(dp), dimension(3,3, ddx_data % constants % ncav), intent(in)  :: hessianphi
  real(dp), dimension(ddx_data % constants % nbasis, &
      & ddx_data % params % nsph), intent(in) :: psi
  real(dp), intent(in) :: tol
  ! Outputs
  real(dp), intent(out)      :: esolv
  real(dp), dimension(3, ddx_data % params % nsph), intent(out) :: force
  integer, intent(out) :: info
  integer                    :: istatus
  !!
  !! Xr         : Reaction potential solution (Laplace equation)
  !! Xe         : Extended potential solution (HSP equation)
  !!
  real(dp), allocatable ::   Xr(:,:), Xe(:,:), Xadj_r(:,:), Xadj_e(:,:)

  !! phi_grid: Phi evaluated at grid points
  real(dp), allocatable :: g(:,:), f(:,:), phi_grid(:, :)
  
  !
  ! Allocate Bessel's functions of the first kind, the second kind, their derivatives,
  ! coefvec, Pchi, and termimat
  !
  call ddlpb_init(ddx_data)

  allocate(Xr(ddx_data % constants % nbasis, ddx_data % params % nsph),&
           & Xe(ddx_data % constants % nbasis, ddx_data % params % nsph), &
           & Xadj_r(ddx_data % constants % nbasis, ddx_data % params % nsph),&
           & Xadj_e(ddx_data % constants % nbasis, ddx_data % params % nsph), &
           & g(ddx_data % params % ngrid, ddx_data % params % nsph),&
           & f(ddx_data % params % ngrid, ddx_data % params % nsph), &
           & phi_grid(ddx_data % params % ngrid, ddx_data % params % nsph), &
           & stat = istatus)
  if (istatus.ne.0) write(6,*) 'ddlpb allocation failed'

  !! Setting initial values to zero
  g = zero; f = zero
  phi_grid = zero;

  !! wghpot: Weigh potential at cavity points. Comes from ddCOSMO
  !!         Intermediate computation of G_0 Eq.(77) QSM19.SISC
  !!
  !! @param[in]  phi      : Boundary conditions (This is psi_0 Eq.(20) QSM19.SISC)
  !! @param[out] phi_grid : Phi evaluated at grid points
  !! @param[out] g        : Boundary conditions on solute-solvent boundary gamma_j_e
  !!
  !call wghpot(ddx_data, phi, phi_grid, g)
  call wghpot(ddx_data % constants % ncav, phi, ddx_data % params % nsph, ddx_data % params % ngrid, &
      & ddx_data % constants % ui, phi_grid, g)
  !!
  !! wghpot_f : Intermediate computation of F_0 Eq.(75) from QSM19.SISC
  !!
  !! @param[in]  gradphi : Gradient of psi_0
  !! @param[out] f       : Boundary conditions scaled by characteristic function
  !!
  call wghpot_f(ddx_data, gradphi,f)

  ! Call the subroutine to solve for Esolv
  call ddx_lpb_solve(ddx_data, g, f, Xr, Xe, tol, esolv)
  ! Start the Force computation
  if(ddx_data % params % force .eq. 1) then
    write(*,*) 'Computation of Forces for ddLPB'
    ! Call the subroutine adjoint to solve the adjoint solution
    call ddx_lpb_adjoint(ddx_data, psi, tol, Xadj_r, Xadj_e)
    !Call the subroutine to evaluate derivatives
    call ddx_lpb_force(ddx_data, hessianphi, phi_grid, gradphi, &
                     & Xr, Xe, Xadj_r, Xadj_e, force)

  endif
  call ddlpb_free(ddx_data)
  deallocate(Xr, Xe,&
           & Xadj_r, Xadj_e, &
           & g, f, phi_grid, stat = istatus)
  if (istatus.ne.0) write(6,*) 'ddlpb deallocation failed'

  return
  end subroutine ddlpb
  !!
  !! Allocate Bessel's functions of the first kind, the second kind,
  !! coefY, C_ik, coefvec, Pchi, and termimat
  !! @param[in] ddx_data : Data Type
  !!
  subroutine ddlpb_init(ddx_data)
  implicit none
  type(ddx_type), intent(inout)  :: ddx_data
  integer                     :: istatus, isph, igrid, ind, icav, l, l0, m0, ind0, jsph
  real(dp)                    :: termi, termk, term, rho, ctheta, stheta, cphi, sphi, rijn
  real(dp), dimension(3)      :: sijn ,vij
  real(dp), dimension(ddx_data % constants % nbasis) :: basloc, vplm
  real(dp), dimension(ddx_data % params % lmax+1) :: vcos, vsin
  real(dp), dimension(0:ddx_data % constants % lmax0) :: SK_rijn, DK_rijn
  integer :: ibasis, ibasis0
  real(dp) :: val
  return
  end subroutine ddlpb_init


  !!
  !! Find intermediate F0 in the RHS of the ddLPB model given in Eq.(82)
  !! @param[in]  gradphi : Gradient of psi_0
  !! @param[out] f       : Intermediate calculation of F0
  !!
  subroutine wghpot_f(ddx_data, gradphi, f )
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  real(dp), dimension(3, ddx_data % constants % ncav),       intent(in)  :: gradphi
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph),    intent(out) :: f

  integer :: isph, ig, ic, ind, ind0, jg, l, m, jsph
  real(dp) :: nderphi, sumSijn, rijn, coef_Ylm, sumSijn_pre, termi, &
      & termk, term
  real(dp), dimension(3) :: sijn, vij
  real(dp) :: rho, ctheta, stheta, cphi, sphi
  real(dp), allocatable :: SK_rijn(:), DK_rijn(:)

  integer :: l0, m0, icav, istatus
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph) :: c0
  real(dp), allocatable :: vplm(:), basloc(:), vcos(:), vsin(:)

  allocate(vplm(ddx_data % constants % nbasis), &
      & basloc(ddx_data % constants % nbasis), &
      & vcos(ddx_data % params % lmax+1), &
      & vsin(ddx_data % params % lmax+1))
  allocate(SK_rijn(0:ddx_data % constants % lmax0),DK_rijn(0:ddx_data % constants % lmax0))
  ic = 0 ; f(:,:)=0.d0
  c0 = zero
  !
  ! Compute c0 Eq.(98) QSM19.SISC
  !
  do isph = 1, ddx_data % params % nsph
    do ig = 1, ddx_data % params % ngrid
      if ( ddx_data % constants % ui(ig,isph).ne.zero ) then
        ic = ic + 1
        nderphi = dot_product( gradphi(:,ic),ddx_data % constants % cgrid(:,ig) )
        c0(:, isph) = c0(:,isph) + &
                     & ddx_data % constants % wgrid(ig)*ddx_data % constants % ui(ig,isph)*&
                     & nderphi*ddx_data % constants % vgrid(:,ig)
      end if
    end do
  end do


  ! Computation of F0 using above terms
  ! icav: External grid poitns
  icav = 0
  do isph = 1, ddx_data % params % nsph
    do ig = 1, ddx_data % params % ngrid
      if (ddx_data % constants % ui(ig,isph).gt.zero) then
        icav = icav + 1
        sumSijn = zero
        ! Loop to compute Sijn
        do jsph = 1, ddx_data % params % nsph
          sumSijn_pre = sumSijn
          vij  = ddx_data % params % csph(:,isph) + &
              & ddx_data % params % rsph(isph)* &
              & ddx_data % constants % cgrid(:,ig) - &
              & ddx_data % params % csph(:,jsph)
          rijn = sqrt(dot_product(vij,vij))
          sijn = vij/rijn
          
          ! Compute Bessel function of 2nd kind for the coordinates
          ! (s_ijn, r_ijn) and compute the basis function for s_ijn
          call modified_spherical_bessel_second_kind( &
              & ddx_data % constants % lmax0, rijn*ddx_data % params % kappa,&
                                                   & SK_rijn, DK_rijn)
          call ylmbas(sijn , rho, ctheta, stheta, cphi, &
                      & sphi, ddx_data % params % lmax, ddx_data % constants % vscales, &
                      & basloc, vplm, vcos, vsin)

          do l0 = 0, ddx_data % constants % lmax0
            term = SK_rijn(l0)/ddx_data % constants % SK_ri(l0,jsph)
            ! coef_Ylm : (der_i_l0/i_l0 - der_k_l0/k_l0)^(-1)*k_l0(r_ijn)/k_l0(r_i)
            coef_Ylm =  ddx_data % constants % C_ik(l0,jsph)*term
            do m0 = -l0, l0
              ind0 = l0**2 + l0 + m0 + 1
              sumSijn = sumSijn + c0(ind0,jsph)*coef_Ylm*basloc(ind0)
              !coefY(icav,ind0,jsph) = coef_Ylm*basloc(ind0)
            end do
          end do
        end do
        !
        ! Here Intermediate value of F_0 is computed Eq. (99)
        ! Mutilplication with Y_lm and weights will happen afterwards
        !write(6,*) sumSijn, epsp, eps, ddx_data % constants % ui(ig,isph)
        f(ig,isph) = -(epsp/ddx_data % params % eps)*ddx_data % constants % ui(ig,isph) * sumSijn
      end if
    end do
  end do 

  deallocate( vplm, basloc, vcos, vsin, SK_rijn, DK_rijn  )
  return
  end subroutine wghpot_f

  !
  ! Subroutine used for the GMRES solver
  ! @param[in]      n : Size of the matrix
  ! @param[in]      x : Input vector
  ! @param[in, out] y : y=A*x
  !
  subroutine bx(params, constants, workspace, x, y)
  implicit none 
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
  real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
  real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
  integer :: isph, istatus
  real(dp), allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)
  integer :: i
  ! allocate workspaces
  allocate( pot(params % ngrid), vplm(constants % nbasis), basloc(constants % nbasis), &
            & vcos(params % lmax+1), vsin(params % lmax+1) , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'Bx: allocation failed !'
    stop
  endif
  
!  if (ddx_data % iprint .ge. 5) then
!      call prtsph('X', ddx_data % constants % nbasis, ddx_data % params % lmax, ddx_data % params % nsph, 0, &
!          & x)
!  end if

  y = zero
  do isph = 1, params % nsph
    call calcv2_lpb(params, constants, isph, pot, x, basloc, vplm, vcos, vsin )
    ! intrhs comes from ddCOSMO
    call intrhs(1, constants % nbasis, params % ngrid, &
                constants % vwgrid, constants % vgrid_nbasis, &
                & pot, y(:,isph) )
    ! Action of off-diagonal blocks
    y(:,isph) = - y(:,isph)
    ! Add action of diagonal block
    y(:,isph) = y(:,isph) + x(:,isph)
  end do
  
!  if (ddx_data % iprint .ge. 5) then
!      call prtsph('Bx (off-diagonal)', ddx_data % constants % nbasis, ddx_data % params % lmax, &
!          & ddx_data % params % nsph, 0, y)
!  end if
  deallocate( pot, basloc, vplm, vcos, vsin , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'bx: allocation failed !'
    stop
  endif
  end subroutine bx

  !!
  !! Scale the ddCOSMO solution vector
  !! @param[in]      direction : Direction of the scaling
  !! @param[in, out] vector    : ddCOSMO solution vector
  !!
  subroutine convert_ddcosmo(ddx_data, direction, vector)
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  integer, intent(in) :: direction
  real(dp), intent(inout) :: vector(ddx_data % constants % nbasis, ddx_data % params % nsph)
  integer :: isph, l, m, ind
  real(dp) :: fac
  
  do isph = 1, ddx_data % params % nsph
    do l = 0, ddx_data % params % lmax
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
  subroutine lpb_hsp(ddx_data, tol, rhs, Xe)
  implicit none
  type(ddx_type), intent(inout)  :: ddx_data
  real(dp), intent(in) :: tol
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in) :: rhs
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(inout) :: Xe
  integer :: isph, istatus, n_iter, info, c1, c2, cr
  real(dp) :: r_norm
  
  Xe = rhs
  
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
  !! @param[in]      Stopping      : Stopping criteria, Default set to
  !!                                 'rel' for relative. Other option
  !!                                 'abs' for absolute
  !! @param[in]      r_norm        : Residual measure
  !! @param[in]      bx        : Subroutine A*x. Named matabx in file
  !! @param[in, out] info          : Flag after solve. 0 means within tolerance
  !!                                 1 means max number of iteration
  n_iter = ddx_data % params % maxiter
  call gmresr(ddx_data % params, ddx_data % constants, ddx_data % workspace, &
      & tol, rhs, Xe, n_iter, r_norm, bx, info)
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
  subroutine calcv2_lpb (params, constants, isph, pot, x, basloc, vplm, vcos, vsin )
  type(ddx_params_type), intent(in)  :: params
  type(ddx_constants_type), intent(in)  :: constants
  integer, intent(in) :: isph
  real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
  real(dp), dimension(params % ngrid), intent(inout) :: pot
  real(dp), dimension(constants % nbasis), intent(inout) :: basloc
  real(dp), dimension(constants % nbasis), intent(inout) :: vplm
  real(dp), dimension(params % lmax+1), intent(inout) :: vcos
  real(dp), dimension(params % lmax+1), intent(inout) :: vsin
  real(dp), dimension(constants % nbasis) :: fac_cosmo, fac_hsp
  integer :: its, ij, jsph
  real(dp) :: rho, ctheta, stheta, cphi, sphi
  real(dp) :: vij(3), sij(3)
  real(dp) :: vvij, tij, xij, oij

  pot = zero
  do its = 1, params % ngrid
    if (constants % ui(its,isph).lt.one) then
      do ij = constants % inl(isph), constants % inl(isph+1)-1
        jsph = constants % nl(ij)

        ! compute geometrical variables
        vij  = params % csph(:,isph) + params % rsph(isph)*constants % cgrid(:,its) - params % csph(:,jsph)
        vvij = sqrt(dot_product(vij,vij))
        tij  = vvij/params % rsph(jsph) 

        if ( tij.lt.one ) then
          sij = vij/vvij
          call ylmbas(sij, rho, ctheta, stheta, cphi, &
                      & sphi, params % lmax, &
                      & constants % vscales, basloc, &
                      & vplm, vcos, vsin)
          call inthsp(params, constants, vvij, params % rsph(jsph), jsph, basloc, fac_hsp)
          xij = fsw(tij, params % se, params % eta)
          if (constants % fi(its,isph).gt.one) then
            oij = xij/constants % fi(its, isph)
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
  subroutine inthsp(params, constants, rijn, ri, isph, basloc, fac_hsp)
  implicit none
  type(ddx_params_type), intent(in)  :: params
  type(ddx_constants_type), intent(in)  :: constants
  integer, intent(in) :: isph
  real(dp), intent(in) :: rijn, ri
  real(dp), dimension(constants % nbasis), intent(in) :: basloc
  real(dp), dimension(constants % nbasis), intent(inout) :: fac_hsp
  real(dp), dimension(0:params % lmax) :: SI_rijn, DI_rijn
  integer :: l, m, ind

  SI_rijn = 0
  DI_rijn = 0
  fac_hsp = 0

  ! Computation of modified spherical Bessel function values      
  call modified_spherical_bessel_first_kind(params % lmax, rijn*params % kappa, SI_rijn, DI_rijn)
  
  do l = 0, params % lmax
    do  m = -l, l
      ind = l*l + l + 1 + m
      fac_hsp(ind) = SI_rijn(l)/constants % SI_ri(l,isph)*basloc(ind)
    end do
  end do
  endsubroutine inthsp

  !
  ! Update the RHS in outer iteration
  ! @param[in] rhs_cosmo_init : G_0
  ! @param[in] rhs_hsp_init   : F_0
  ! @param[in, out] rhs_cosmo : -C_1*X_r^(k-1) - C_2*X_e^(k-1) + G_0 + F_0
  ! @param[in, out] rhs_hsp   : -C_1*X_r^(k-1) - C_2*X_e^(k-1) + F_0
  ! @param[in] Xr             : X_r^(k-1)
  ! @param[in] Xe             : X_e^(k-1)
  !
  subroutine update_rhs(params, constants, rhs_cosmo_init, rhs_hsp_init, rhs_cosmo, & 
      & rhs_hsp, Xr, Xe)
  implicit none
  type(ddx_params_type), intent(in)  :: params
  type(ddx_constants_type), intent(in)  :: constants
  real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: rhs_cosmo_init, &
      & rhs_hsp_init
  real(dp), dimension(constants % nbasis, params % nsph), intent(inout) :: rhs_cosmo, rhs_hsp
  real(dp), dimension(constants % nbasis, params % nsph) :: rhs_plus
  real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xr, Xe
  integer :: isph, jsph, igrid, icav, ind, l, m, ind0, istatus
  real(dp), dimension(3) :: vij
  real(dp), dimension(constants % nbasis, params % nsph) :: diff_re_c1_c2
  real(dp), dimension(constants % nbasis0, params % nsph) :: diff0
  real(dp), dimension(constants % nbasis, constants % nbasis, params % nsph) :: smat
  real(dp), dimension(constants % ncav) :: diff_ep
  real(dp) :: Qval, rijn, val
  integer :: c0, cr, c_qmat, c_init, c_ep0, c_ep1 !, nbasis_appro
      
  ! diff_re = epsp/eps*l1/ri*Xr - i'(ri)/i(ri)*Xe,
  diff_re_c1_c2 = zero
  do jsph = 1, params % nsph
    do l = 0, params % lmax
      do m = -l,l
        ind = l**2 + l + m + 1
        diff_re_c1_c2(ind,jsph) = ((epsp/params % eps)*&
                                & (l/params % rsph(jsph)) * &
                                & Xr(ind,jsph) - constants % termimat(l,jsph)*Xe(ind,jsph))
      end do
    end do
  end do

  ! diff0 = Pchi * diff_er, linear scaling
  ! TODO: probably doing PX on the fly is better 
  diff0 = zero 
  do jsph = 1, params % nsph
    do ind0 = 1, constants % nbasis0
      diff0(ind0, jsph) = dot_product(diff_re_c1_c2(:,jsph), &
          & constants % Pchi(:,ind0, jsph))
    end do
  end do

  ! diff_ep = diff0 * coefY,    COST: M^2*nbasis*Nleb
  diff_ep = zero
  do icav = 1, constants % ncav
    val = zero
    do jsph = 1, params % nsph 
      do ind0 = 1, constants % nbasis0
        val = val + diff0(ind0,jsph)*constants % coefY(icav,ind0,jsph)
      end do
    end do
    diff_ep(icav) = val
  end do

  rhs_plus = zero
  icav = 0
  do isph = 1, params % nsph
    do igrid = 1, params % ngrid
      if (constants % ui(igrid,isph).gt.zero) then
        icav = icav + 1
        do ind = 1, constants % nbasis
          rhs_plus(ind,isph) = rhs_plus(ind,isph) + &
              & constants % coefvec(igrid,ind,isph)*diff_ep(icav)
        end do
      end if
    end do
  end do

  rhs_cosmo = rhs_cosmo_init - rhs_plus
  rhs_hsp = rhs_hsp_init - rhs_plus

  return
  end subroutine update_rhs  

  !
  ! Computation for Solvation energy
  ! @param[in]  ddx_data   : Input data file
  ! @param[in]  g          : Intermediate matrix for computation of g0
  ! @param[in]  f          : Intermediate matrix for computation of f0
  ! @param[out] Xr         : Solution corresponding to COSMO
  ! @param[out] Xe         : Solution corresponding to HSP
  ! @param[out] esolv      : Solvation energy
  subroutine ddx_lpb_solve(ddx_data, g, f, &
                          & Xr, Xe, tol, esolv)
  implicit none
  type(ddx_type), intent(inout) :: ddx_data
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph), intent(in) :: g,f
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(out) :: Xr, Xe
  real(dp), intent(in) :: tol
  real(dp), intent(out) :: esolv

  ! Local variables
  !! g0      : Vector associated to psi_0 Eq.(77) QSM19.SISC
  !! f0      : Vector associated to partial_n_psi_0 Eq.(99) QSM19.SISC
  real(dp), dimension(ddx_data % constants % nbasis) :: g0, f0
  ! rhs_r_init : Initial RHS for COSMO, G0
  ! rhs_e_init : Initial RHS for HSP, F0
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph):: rhs_r_init, rhs_e_init
  !! rhs_r      : Right hand side corresponding to Laplace equation
  !! rhs_e      : Right hand side corresponding to HSP equation
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph):: rhs_r, rhs_e
  integer  :: iteration, n_iter, isph, info
  real(dp) :: old_esolv, inc
  logical  :: converged, ok

  ! Setting of the local variables
  old_esolv = zero; inc = zero
  rhs_r_init = zero; rhs_e_init = zero
  g0 = zero; f0 = zero

  !!
  !! Integrate Right hand side
  !! rhs_r_init: g0+f0
  !! rhs_e_init: f0
  !!
!  do isph = 1, ddx_data % params % nsph
!    !! intrhs is a subroutine in ddx_operators
!    !! @param[in]  isph : Sphere number, used for output
!    !! @param[in]  g    : Intermediate right side g
!    !! @param[out] g0   : Integrated right side Eq.(77) in QSM19.SISC
!    call intrhs(ddx_data % iprint, ddx_data % params % ngrid, &
!                ddx_data % params % lmax, ddx_data % constants % vwgrid, ddx_data % constants % vgrid_nbasis, &
!                & isph, g(:,isph), g0)
!    call intrhs(ddx_data % iprint, ddx_data % params % ngrid, &
!                ddx_data % params % lmax, ddx_data % constants % vwgrid, ddx_data % constants % vgrid_nbasis, &
!                & isph,f(:,isph),f0)
!    !! rhs
!    rhs_r_init(:,isph) = g0 + f0
!    rhs_e_init(:,isph) = f0
!  end do
  call intrhs(ddx_data % params % nsph, ddx_data % constants % nbasis, &
      & ddx_data % params % ngrid, ddx_data % constants % vwgrid, &
      & ddx_data % constants % vgrid_nbasis, g, rhs_r_init)
  call intrhs(ddx_data % params % nsph, ddx_data % constants % nbasis, &
      & ddx_data % params % ngrid, ddx_data % constants % vwgrid, &
      & ddx_data % constants % vgrid_nbasis, f, rhs_e_init)
  rhs_r_init = rhs_r_init + rhs_e_init


  rhs_r = rhs_r_init
  rhs_e = rhs_e_init


  ! Set values to default values
  first_out_iter = .true.
  converged = .false.
  ok = .false.
  iteration = one

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
    !!                             matrix to vector, i.e., L^{-1}x_r, comes from
    !!                             matvec.f90
    !! @param[in]      hnorm     : Use defined norm, comes from matvec.f90
    Xr = zero
    n_iter = ddx_data % params % maxiter
    call jacobi_diis(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, tol, rhs_r, Xr, n_iter, ddx_data % xs_rel_diff, &
        & lx_nodiag, ldm1x, hnorm, info)
    ! Scale by the factor of (2l+1)/4Pi
    call convert_ddcosmo(ddx_data, 1, Xr)

    !! Solve ddLPB step
    !! B X_e = RHS_e (= F_0)
    call lpb_hsp(ddx_data, tol, rhs_e, Xe)

    !! Update the RHS
    !! / RHS_r \ = / g + f \ - / c1 c2 \ / X_r \
    !! \ RHS_e /   \ f     /   \ c1 c2 / \ X_e /
    call update_rhs(ddx_data % params, ddx_data % constants, &
        & rhs_r_init, rhs_e_init, rhs_r, rhs_e, Xr, Xe)

    !! Compute energy
    !! esolv = pt5*sprod(nsph*nylm,xr,psi)
    esolv = zero
    do isph = 1, ddx_data % params % nsph
      esolv = esolv + pt5*ddx_data % params % charge(isph)*Xr(1,isph)*(one/sqrt4pi)
    end do

    !! Check for convergence
    inc = zero
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
  end subroutine ddx_lpb_solve

  !
  ! Computation of Adjoint
  ! @param[in] ddx_data: Input data file
  ! @param[in] psi     : psi_r
  subroutine ddx_lpb_adjoint(ddx_data, psi, tol, Xadj_r, Xadj_e)
  implicit none
  type(ddx_type), intent(inout)           :: ddx_data
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in) :: psi
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(out) :: Xadj_r, Xadj_e
  real(dp), intent(in) :: tol
  real(dp), external :: dnrm2
  ! Local Variables
  ! ok             : Logical expression used in Jacobi solver
  ! istatus        : Status of allocation
  ! isph           : Index for the sphere
  ! ibasis         : Index for the basis
  ! iteration      : Number of outer loop iterations
  ! inc            : Check for convergence threshold
  ! relative_num   : Numerator of the relative error
  ! relative_denom : Denominator of the relative error
  ! converged      : Convergence check for outer loop
  logical   :: ok
  integer   :: istatus, isph, ibasis,  iteration
  real(dp)  :: inc, relative_num, relative_denom
  logical   :: converged

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

  ! GMRES parameters
  integer :: info, n_iter
  real(dp) r_norm

  ! Allocation
  allocate(rhs_cosmo(ddx_data % constants % nbasis, ddx_data % params % nsph), &
           & rhs_hsp(ddx_data % constants % nbasis, ddx_data % params % nsph), &
           & rhs_cosmo_init(ddx_data % constants % nbasis, ddx_data % params % nsph), &
           & rhs_hsp_init(ddx_data % constants % nbasis, ddx_data % params % nsph),&
           & X_r_k_1(ddx_data % constants % nbasis, ddx_data % params % nsph),&
           & X_e_k_1(ddx_data % constants % nbasis, ddx_data % params % nsph),&
           & stat = istatus)
  if(istatus .ne. 0) then
    write(*,*) 'Allocation failed in adjoint LPB!'
    stop
  end if
  ! We compute the adjoint solution first
  write(*,*) 'Solution of adjoint system'
  ! Set local variables
  iteration = one
  inc = zero
  relative_num = zero
  relative_denom = zero
  X_r_k_1 = zero
  X_e_k_1 = zero
  ! Initial RHS
  ! rhs_cosmo_init = psi_r
  ! rhs_hsp_init   = psi_e (=0)
  rhs_cosmo_init = psi
  rhs_hsp_init = zero
  ! Updated RHS
  rhs_cosmo = rhs_cosmo_init
  rhs_hsp = rhs_hsp_init
  ! Initial Xadj_r and Xadj_e
  Xadj_r = zero
  Xadj_e = zero
  ! Logical variable for the first outer iteration
  first_out_iter = .true.
  converged = .false.
  ok = .false.
  ! Solve the adjoint system
  do while (.not.converged)
    ! Solve A*X_adj_r = psi_r
    ! Set the RHS to correct form by factoring with 4Pi/2l+1
    Xadj_r = zero
    call convert_ddcosmo(ddx_data, 1, rhs_cosmo)
    n_iter = ddx_data % params % maxiter
    call jacobi_diis(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, tol, rhs_cosmo, Xadj_r, n_iter, &
                    & ddx_data % xs_rel_diff, lstarx_nodiag, ldm1x, hnorm, info)

    ! Solve the HSP equation
    ! B*X_adj_e = psi_e
    ! For first iteration the rhs is zero for HSP equation. Hence, Xadj_e = 0
    Xadj_e = rhs_hsp
    if(iteration.ne.1) then
      n_iter = ddx_data % params % maxiter
      call gmresr(ddx_data % params, ddx_data % constants, &
          & ddx_data % workspace, tol, rhs_hsp, Xadj_e, n_iter, &
          & r_norm, bstarx, info)
    endif

    ! Update the RHS
    ! |rhs_r| = |psi_r|-|C1* C1*||Xadj_r|
    ! |rhs_e| = |psi_e| |C2* C2*||Xadj_e|
    call update_rhs_adj(ddx_data, rhs_cosmo_init, rhs_hsp_init,&
                        & rhs_cosmo, rhs_hsp, Xadj_r, Xadj_e)
    ! Stopping Criteria.
    ! Checking the relative error of Xadj_r
    inc  = zero
    inc = dnrm2(ddx_data % constants % n, Xadj_r + Xadj_e - X_r_k_1 - X_e_k_1, 1)/ &
        & dnrm2(ddx_data % constants % n, Xadj_r + Xadj_e, 1)

    ! Store the previous step solution
    X_r_k_1 = Xadj_r
    X_e_k_1 = Xadj_e
    if ((iteration .gt. 1) .and. (inc.lt.tol)) then
      write(6,*) 'Reach tolerance.'
      converged = .true.
    end if
    write(6,*) 'Adjoint computation :', iteration, inc
    iteration = iteration + 1
    first_out_iter = .false.
  end do
!  if (ddx_data % iprint .ge. 5) then
!    call prtsph('Xadj_r', ddx_data % constants % nbasis, ddx_data % params % lmax, &
!          & ddx_data % params % nsph, 0, Xadj_r)
!    call prtsph('Xadj_e', ddx_data % constants % nbasis, ddx_data % params % lmax, &
!          & ddx_data % params % nsph, 0, Xadj_e)
!  end if
  ! Deallocation
  deallocate(rhs_cosmo, rhs_hsp, rhs_cosmo_init, &
           & rhs_hsp_init, X_r_k_1, X_e_k_1, stat = istatus)
  if(istatus .ne. 0) then
    write(*,*) 'Deallocation failed in adjoint LPB!'
    stop
  end if
  end subroutine ddx_lpb_adjoint

  !
  ! Computation for Solvation energy
  ! @param[in]  ddx_data   : Input data file
  ! @param[in]  hessian    : Hessian of Psi
  ! @param[in]  phi_grid   : Phi evaluated at the grid point
  ! @param[in]  gradphi  : Gradient of phi
  ! @param[in]  Xr         : Solution corresponding to COSMO
  ! @param[in]  Xe         : Solution corresponding to HSP
  ! @param[in]  Xadj_r     : Solution corresponding to adjoint of the COSMO
  ! @param[in]  Xadj_e     : Solution corresponding to adjoint of the HSP
  ! @param[out] force      : Force
  subroutine ddx_lpb_force(ddx_data, hessian, phi_grid, gradphi, &
                         & Xr, Xe, Xadj_r, Xadj_e, force)
  implicit none
  type(ddx_type), intent(inout) :: ddx_data
  real(dp), dimension(3,3, ddx_data % constants % ncav) :: hessian
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph), intent(in) :: phi_grid
  real(dp), dimension(3, ddx_data % constants % ncav), intent(in)    :: gradphi
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in) :: Xr, Xe
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in) :: Xadj_r, Xadj_e
  real(dp), dimension(3, ddx_data % params % nsph), intent(out) :: force

  ! Local variable
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph) :: Xadj_r_sgrid, Xadj_e_sgrid
  real(dp), dimension(3, ddx_data % constants % ncav) :: normal_hessian_cav
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph) :: diff_re
  real(dp), dimension(ddx_data % constants % nbasis) :: basloc, vplm
  real(dp), dimension(3, ddx_data % constants % nbasis) :: dbasloc
  ! Local variables, used in force computation
  ! ef : Electrostatic force using potential of spheres
  real(dp), dimension(3, ddx_data % params % nsph) :: ef
  real(dp), dimension(ddx_data % params % lmax + 1) :: vsin, vcos
  !! scaled_Xr  : Reaction potential scaled by 1/(4Pi/2l+1)
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph) :: scaled_Xr
  ! icav_gr : Index for global cavity points for Laplace
  ! icav_ge : Index for global cavity points for HSP
  ! ok      : Input argument for Jacobi solver
  ! n_iter  : Number of iterative steps
  integer                    :: isph, icav, icav_gr, icav_ge, igrid
  integer                    :: i
  real(dp), external :: ddot

  Xadj_r_sgrid = zero; Xadj_e_sgrid = zero
  diff_re = zero; normal_hessian_cav = zero
  vsin = zero; vcos = zero; vplm = zero
  basloc = zero; dbasloc = zero; ef = zero
  force = zero

  ! Compute the derivative of the normal derivative of psi_0
  icav = 0
  do isph = 1, ddx_data % params % nsph
    do igrid = 1, ddx_data % params % ngrid
      if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
        icav = icav + 1
        do i = 1, 3
          normal_hessian_cav(:, icav) = normal_hessian_cav(:,icav) +&
                                  & hessian(:,i,icav)*ddx_data % constants % cgrid(i,igrid)
        end do
      end if
    end do
  end do

  ! Call dgemm to integrate the adjoint solution on the grid points
  call dgemm('T', 'N', ddx_data % params % ngrid, ddx_data % params % nsph, &
          & ddx_data % constants % nbasis, one, ddx_data % constants % vgrid, &
          & ddx_data % constants % vgrid_nbasis, &
          & Xadj_r , ddx_data % constants % nbasis, zero, Xadj_r_sgrid, &
          & ddx_data % params % ngrid)
  call dgemm('T', 'N', ddx_data % params % ngrid, ddx_data % params % nsph, &
          & ddx_data % constants % nbasis, one, ddx_data % constants % vgrid, &
          & ddx_data % constants % vgrid_nbasis, &
          & Xadj_e , ddx_data % constants % nbasis, zero, Xadj_e_sgrid, &
          & ddx_data % params % ngrid)

  ! Scale by the factor of 1/(4Pi/(2l+1))
  scaled_Xr = Xr
  call convert_ddcosmo(ddx_data, -1, scaled_Xr)

  do isph = 1, ddx_data % params % nsph
    ! Compute A^k*Xadj_r, using Subroutine from ddCOSMO
    call fdoka(ddx_data % params, ddx_data % constants, isph, scaled_Xr, Xadj_r_sgrid(:, isph), &
              & basloc, dbasloc, vplm, vcos, vsin, force(:,isph))
    call fdokb(ddx_data % params, ddx_data % constants, isph, scaled_Xr, Xadj_r_sgrid, basloc, &
              & dbasloc, vplm, vcos, vsin, force(:, isph))
    ! Compute B^k*Xadj_e
    call fdoka_b_xe(ddx_data, isph, Xe, Xadj_e_sgrid(:, isph), &
              & basloc, dbasloc, vplm, vcos, vsin, force(:,isph))
    call fdokb_b_xe(ddx_data, isph, Xe, Xadj_e_sgrid, &
              & basloc, dbasloc, vplm, vcos, vsin, force(:, isph))
    ! Compute C1 and C2 contributions
    diff_re = zero
    call fdouky(ddx_data, isph, &
                & Xr, Xe, &
                & Xadj_r_sgrid, Xadj_e_sgrid, &
                & Xadj_r, Xadj_e, &
                & force(:, isph), &
                & diff_re)
    call derivative_P(ddx_data, isph,&
                & Xr, Xe, &
                & Xadj_r_sgrid, Xadj_e_sgrid, &
                & diff_re, &
                & force(:, isph))
    ! Computation of G0
    call fdoga(ddx_data % params, ddx_data % constants, isph, Xadj_r_sgrid, phi_grid, force(:, isph))
  end do
  ! Computation of G0 continued
  ! NOTE: fdoga returns a positive summation
  force = -force
  icav = 0
  do isph = 1, ddx_data % params % nsph
    do igrid = 1, ddx_data % params % ngrid
      if(ddx_data % constants % ui(igrid, isph) .ne. zero) then
        icav = icav + 1
        ddx_data % zeta(icav) = -ddx_data % constants % wgrid(igrid) * &
                      & ddx_data % constants % ui(igrid, isph) * ddot(ddx_data % constants % nbasis, &
                      & ddx_data % constants % vgrid(1, igrid), 1, &
                      & Xadj_r(1, isph), 1)
        force(:, isph) = force(:, isph) + &
                      & ddx_data % zeta(icav)*gradphi(:, icav)
      end if
    end do
  end do
  call efld(ddx_data % constants % ncav, ddx_data % zeta, ddx_data % constants % ccav, &
              & ddx_data % params % nsph, ddx_data % params % csph, ef)
  do isph = 1, ddx_data % params % nsph
    force(:, isph) = force(:, isph) - ef(:, isph)*ddx_data % params % charge(isph)
  end do

  icav_gr = zero
  icav_ge = zero
  do isph = 1, ddx_data % params % nsph
    ! Computation of F0
    call fdouky_f0(ddx_data, isph, Xadj_r, Xadj_r_sgrid, &
        & gradphi, force(:, isph))
    call fdoco(ddx_data, isph, Xadj_r_sgrid, gradphi, &
           & normal_hessian_cav, icav_gr, force(:, isph))
    ! Computation of F0
    call fdouky_f0(ddx_data, isph, Xadj_e, Xadj_e_sgrid, &
        & gradphi, force(:, isph))
    call fdoco(ddx_data, isph, Xadj_e_sgrid, gradphi, &
           & normal_hessian_cav, icav_ge, force(:, isph))
  end do
  force = pt5*force

  end subroutine ddx_lpb_force

  !! Computation of Adjoint B, i.e., B*
  !> Apply adjoint single layer operator to spherical harmonics
  !! implementation is similar to lstarx in ddCOSMO
  !! Diagonal blocks are not counted here.
  subroutine bstarx(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
  real(dp), intent(in)       :: x(constants % nbasis, params % nsph)
  ! Output
  real(dp), intent(out)      :: y(constants % nbasis, params % nsph)
  ! Local variables
  integer                    :: isph, igrid, istatus
  real(dp), allocatable      :: xi(:,:), vplm(:), basloc(:), vcos(:), vsin(:)
  ! Allocate workspaces
  allocate(xi(params % ngrid, params % nsph), vplm(constants % nbasis), &
        & basloc(constants % nbasis), vcos(params % lmax+1), &
        & vsin(params % lmax+1), stat=istatus)
  if (istatus .ne. 0) then
      write(*, *) 'bstarx: allocation failed!'
      stop
  endif
!  if (ddx_data % iprint .ge. 5) then
!      call prtsph('X', ddx_data % constants % nbasis, ddx_data % params % lmax, ddx_data % params % nsph, 0, &
!          & x)
!  end if
  ! Initalize
  y = zero
  !! Expand x over spherical harmonics
  ! Loop over spheres
  do isph = 1, params % nsph
      ! Loop over grid points
      do igrid = 1, params % ngrid
          xi(igrid, isph) = dot_product(x(:, isph), &
              & constants % vgrid(:constants % nbasis, igrid))
      end do
  end do
  !! Compute action
  ! Loop over spheres
  do isph = 1, params % nsph
      ! Compute NEGATIVE action of off-digonal blocks
      call adjrhs_lpb(params, constants, workspace, isph, xi, y(:, isph), &
          & basloc, vplm, vcos, vsin)
      y(:, isph) = - y(:, isph)
      y(:,isph)  = y(:,isph) + x(:,isph)
  end do
!  if (ddx_data % iprint .ge. 5) then
!      call prtsph('B*X (off-diagonal)', ddx_data % constants % nbasis, ddx_data % params % lmax, &
!          & ddx_data % params % nsph, 0, y)
!  end if
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
  subroutine adjrhs_lpb(params, constants, workspace, isph, xi, vlm, basloc, vplm, vcos, vsin )
  implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
  integer,  intent(in)    :: isph
  real(dp), dimension(params % ngrid, params % nsph), intent(in) :: xi
  real(dp), dimension((params % lmax+1)**2), intent(inout) :: vlm
  real(dp), dimension((params % lmax+1)**2), intent(inout) :: basloc, vplm
  real(dp), dimension(params % lmax+1), intent(inout) :: vcos, vsin

  integer :: ij, jsph, ig, l, ind, m
  real(dp)  :: vji(3), vvji, tji, sji(3), xji, oji, fac
  real(dp) :: rho, ctheta, stheta, cphi, sphi
  real(dp), dimension(constants % nbasis) :: fac_hsp

  !loop over neighbors of i-sphere
  do ij = constants % inl(isph), constants % inl(isph+1)-1
    !j-sphere is neighbor
    jsph = constants % nl(ij)
    !loop over integration points
    do ig = 1, params % ngrid
      !compute t_n^ji = | r_j + \rho_j s_n - r_i | / \rho_i
      vji  = params % csph(:,jsph) + params % rsph(jsph)* &
            & constants % cgrid(:,ig) - params % csph(:,isph)
      vvji = sqrt(dot_product(vji,vji))
      tji  = vvji/params % rsph(isph)
      !point is INSIDE i-sphere (+ transition layer)
      if ( tji.lt.( one + (params % se+one)/two*params % eta ) ) then
        !compute s_n^ji
        sji = vji/vvji
        call ylmbas(sji, rho, ctheta, stheta, cphi, &
                      & sphi, params % lmax, &
                      & constants % vscales, basloc, &
                      & vplm, vcos, vsin)
        call inthsp_adj(params, constants, workspace, vvji, params % rsph(isph), isph, basloc, fac_hsp)
        !compute \chi( t_n^ji )
        xji = fsw( tji, params % se, params % eta )
        !compute W_n^ji
        if ( constants % fi(ig,jsph).gt.one ) then
          oji = xji/ constants % fi(ig,jsph)
        else
          oji = xji
        endif
        !compute w_n * xi(n,j) * W_n^ji
        fac = constants % wgrid(ig) * xi(ig,jsph) * oji
        !loop over l
        do l = 0, params % lmax
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
  subroutine inthsp_adj(params, constants, workspace, rjin, rj, jsph, basloc, fac_hsp)
  implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
  integer, intent(in) :: jsph
  real(dp), intent(in) :: rjin, rj
  real(dp), dimension(constants % nbasis), intent(in) :: basloc
  real(dp), dimension(constants % nbasis), intent(inout) :: fac_hsp
  real(dp), dimension(0:params % lmax) :: SI_rjin, DI_rjin
  integer :: l, m, ind

  SI_rjin = 0
  DI_rjin = 0
  fac_hsp = 0

  ! Computation of modified spherical Bessel function values      
  call modified_spherical_bessel_first_kind(params % lmax, rjin*params % kappa, SI_rjin, DI_rjin)
  
  do l = 0, params % lmax
    do  m = -l, l
      ind = l*l + l + 1 + m
      fac_hsp(ind) = SI_rjin(l)/constants % SI_ri(l,jsph)*basloc(ind)
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
  implicit none
  type(ddx_type), intent(in)  :: ddx_data
  real(dp), dimension(ddx_data % constants % nbasis,ddx_data % params % nsph), intent(in) :: rhs_r_init, &
                                                                       & rhs_e_init
  real(dp), dimension(ddx_data % constants % nbasis,ddx_data % params % nsph), intent(in) :: Xadj_r, Xadj_e
  real(dp), dimension(ddx_data % constants % nbasis,ddx_data % params % nsph), intent(out) :: rhs_r, rhs_e
  ! Local Variables
  ! rhs_r_adj : C_1*(Xadj_r+Xadj_e)
  ! rhs_e_adj : C_2*(Xadj_r+Xadj_e)
  real(dp), dimension(ddx_data % constants % nbasis,ddx_data % params % nsph) :: rhs_adj
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph) :: Xadj_sgrid
  integer :: isph, igrid, icav, ibasis, ibasis0, l, m, ind, jsph
  ! epsilon_ratio : epsilon_1/epsilon_2
  real(dp) :: epsilon_ratio
  ! val_1, val_2  : Intermediate summations variable
  real(dp) :: val

  epsilon_ratio = epsp/ddx_data % params % eps


  ! Call dgemm to integrate the adjoint solution on the grid points
  ! Summation over l' and m'
  call dgemm('T', 'N', ddx_data % params % ngrid, ddx_data % params % nsph, &
            & ddx_data % constants % nbasis, one, ddx_data % constants % vgrid, ddx_data % constants % vgrid_nbasis, &
            & Xadj_r + Xadj_e , ddx_data % constants % nbasis, zero, Xadj_sgrid, &
            & ddx_data % params % ngrid)


  rhs_adj = zero
  do isph = 1, ddx_data % params % nsph
    do ibasis = 1, ddx_data % constants % nbasis
      icav = zero
      !Summation over j and n
      do jsph = 1, ddx_data % params % nsph
        do igrid = 1, ddx_data % params % ngrid
          if (ddx_data % constants % ui(igrid,jsph).gt.zero) then
            icav = icav + 1
            rhs_adj(ibasis, isph) = rhs_adj(ibasis, isph) + &
                        & ddx_data % constants % wgrid(igrid)*&
                        & ddx_data % constants % ui(igrid, jsph)*&
                        & Xadj_sgrid(igrid, jsph)*&
                        & ddx_data % constants % diff_ep_adj(icav, ibasis, isph)
          end if
        end do
      end do
    end do
  end do

  do isph = 1, ddx_data % params % nsph
    do l = 0, ddx_data % params % lmax
      do m = -l, l
        ind = l**2 + l + m + 1
        rhs_r(ind, isph) = rhs_r_init(ind, isph) - &
                        & (epsilon_ratio*l*rhs_adj(ind, isph))/ddx_data % params % rsph(isph)
        rhs_e(ind, isph) = rhs_e_init(ind, isph) + ddx_data % constants % termimat(l,isph)*rhs_adj(ind, isph)
      end do
    end do
  end do

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
  subroutine fdoka_b_xe(ddx_data, isph, Xe, Xadj_e, basloc, dbasloc, &
                       & vplm, vcos, vsin, force_e)
  implicit none
  type(ddx_type), intent(in) :: ddx_data
  integer,                         intent(in)    :: isph
  real(dp),  dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in)   :: Xe
  real(dp),  dimension(ddx_data % params % ngrid),       intent(in)    :: Xadj_e
  real(dp),  dimension(ddx_data % constants % nbasis),      intent(inout) :: basloc, vplm
  real(dp),  dimension(3, ddx_data % constants % nbasis),    intent(inout) :: dbasloc
  real(dp),  dimension(ddx_data % params % lmax+1),      intent(inout) :: vcos, vsin
  real(dp),  dimension(3),           intent(inout) :: force_e

  ! Local Variables
  integer :: igrid, ineigh, jsph, l, ind, m
  real(dp), dimension(0:ddx_data % params % lmax) :: SI_rijn
  real(dp), dimension(0:ddx_data % params % lmax) :: DI_rijn
  ! beta   : Eq.(53) Stamm.etal.18
  ! tlow   : Lower bound for switch region
  ! thigh  : Upper bound for switch region
  ! f1     : First factor in alpha computation
  ! f2     : Second factor in alpha computation
  real(dp)  :: rijn, tij, beta, tlow, thigh, xij, oij, f1, f2, f3
  ! alpha : Eq.(52) Stamm.etal.18
  ! va    : Eq.(54) Stamm.etal.18
  real(dp)  :: vij(3), sij(3), alpha(3), va(3), rj
  real(dp), external :: dnrm2
  
  SI_rijn = 0
  DI_rijn = 0
  
  tlow  = one - pt5*(one - ddx_data % params % se)*ddx_data % params % eta
  thigh = one + pt5*(one + ddx_data % params % se)*ddx_data % params % eta

  ! Loop over grid points
  do igrid = 1, ddx_data % params % ngrid
    va = zero
    do ineigh = ddx_data % constants % inl(isph), ddx_data % constants % inl(isph+1) - 1
      jsph = ddx_data % constants % nl(ineigh)
      vij  = ddx_data % params % csph(:,isph) + &
            & ddx_data % params % rsph(isph)*ddx_data % constants % cgrid(:,igrid) - &
            & ddx_data % params % csph(:,jsph)
      rijn = dnrm2(3, vij, 1)
      tij  = rijn/ddx_data % params % rsph(jsph)
      rj = ddx_data % params % rsph(jsph)

      if (tij.ge.thigh) cycle
      ! Computation of modified spherical Bessel function values      
      call modified_spherical_bessel_first_kind(ddx_data % params % lmax, rijn*ddx_data % params % kappa, SI_rijn, DI_rijn)

      sij  = vij/rijn
      call dbasis(ddx_data % params, ddx_data % constants, sij, basloc, dbasloc, vplm, vcos, vsin)
      alpha  = zero
      do l = 0, ddx_data % params % lmax
        ind = l*l + l + 1
        f1 = (DI_rijn(l)*ddx_data % params % kappa)/ddx_data % constants % SI_ri(l,jsph);
        f2 = SI_rijn(l)/ddx_data % constants % SI_ri(l, jsph)
        do m = -l, l
          alpha(:) = alpha(:) + (f1*sij(:)*basloc(ind+m) + &
                    & (f2/rijn)*dbasloc(:,ind+m))*Xe(ind+m,jsph)
        end do
      end do
      beta = compute_beta(ddx_data, SI_rijn, rijn, jsph, Xe(:,jsph),basloc)
      xij = fsw(tij, ddx_data % params % se, ddx_data % params % eta)
      if (ddx_data % constants % fi(igrid,isph).gt.one) then
        oij = xij/ddx_data % constants % fi(igrid,isph)
        f2  = -oij/ddx_data % constants % fi(igrid,isph)
      else
        oij = xij
        f2  = zero
      end if
      f1 = oij
      va(:) = va(:) + f1*alpha(:) + beta*f2*ddx_data % constants % zi(:,igrid,isph)
      if (tij .gt. tlow) then
        f3 = beta*dfsw(tij,ddx_data % params % se,ddx_data % params % eta)/ddx_data % params % rsph(jsph)
        if (ddx_data % constants % fi(igrid,isph).gt.one) f3 = f3/ddx_data % constants % fi(igrid,isph)
        va(:) = va(:) + f3*sij(:)
      end if
    end do
  force_e = force_e - ddx_data % constants % wgrid(igrid)*Xadj_e(igrid)*va(:)
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
  subroutine fdokb_b_xe(ddx_data, isph, Xe, Xadj_e, basloc, dbasloc, &
                        & vplm, vcos, vsin, force_e)
  implicit none
  type(ddx_type), intent(in) :: ddx_data
  integer,                         intent(in)    :: isph
  real(dp),  dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in)    :: Xe
  real(dp),  dimension(ddx_data % params % ngrid, ddx_data % params % nsph),  intent(in)    :: Xadj_e
  real(dp),  dimension(ddx_data % constants % nbasis),      intent(inout) :: basloc,  vplm
  real(dp),  dimension(3, ddx_data % constants % nbasis),   intent(inout) :: dbasloc
  real(dp),  dimension(ddx_data % params % lmax+1),      intent(inout) :: vcos, vsin
  real(dp),  dimension(3),           intent(inout) :: force_e

  ! Local Variables
  ! jk     : Row pointer over kth row
  integer :: igrid, jsph, ksph, ineigh, l, m, ind, jk
  real(dp), dimension(0:ddx_data % params % lmax) :: SI_rjin, SI_rjkn
  real(dp), dimension(0:ddx_data % params % lmax) :: DI_rjin, DI_rjkn

  logical :: proc
  ! fac     : \delta_fj_n*\omega^\eta_ji
  ! f1      : First factor in alpha computation
  ! f2      : Second factor in alpha computation
  ! beta_ji : Eq.(57) Stamm.etal.18
  ! dj      : Before Eq.(10) Stamm.etal.18
  real(dp)  :: rjin, tji, xji, oji, fac, f1, f2, beta_ji, dj, tlow, thigh
  real(dp)  :: b, beta_jk, g1, g2, rjkn, tjk, xjk
  ! alpha : Eq.(56) Stamm.etal.18
  ! vb    : Eq.(60) Stamm.etal.18
  ! vc    : Eq.(59) Stamm.etal.18
  real(dp)  :: vji(3), sji(3), vjk(3), sjk(3), alpha(3), vb(3), vc(3)
  ! rho    : Argument for ylmbas
  ! ctheta : Argument for ylmbas
  ! stheta : Argument for ylmbas
  ! cphi   : Argument for ylmbas
  ! sphi   : Argument for ylmbas
  real(dp) :: rho, ctheta, stheta, cphi, sphi, ri, arg_bessel

  real(dp), external :: dnrm2
  
  SI_rjin = 0
  DI_rjin = 0
  SI_rjkn = 0
  DI_rjkn = 0

  tlow  = one - pt5*(one - ddx_data % params % se)*ddx_data % params % eta
  thigh = one + pt5*(one + ddx_data % params % se)*ddx_data % params % eta

  do igrid = 1, ddx_data % params % ngrid
    vb = zero
    vc = zero
    do ineigh = ddx_data % constants % inl(isph), ddx_data % constants % inl(isph+1) - 1
      jsph = ddx_data % constants % nl(ineigh)
      vji  = ddx_data % params % csph(:,jsph) + &
              & ddx_data % params % rsph(jsph)*ddx_data % constants % cgrid(:,igrid) - &
              & ddx_data % params % csph(:,isph)
      rjin = dnrm2(3, vji, 1)
      ri = ddx_data % params % rsph(isph)
      tji  = rjin/ri

      if (tji.gt.thigh) cycle

      call modified_spherical_bessel_first_kind(ddx_data % params % lmax, rjin*ddx_data % params % kappa, SI_rjin, DI_rjin)
      
      sji  = vji/rjin
      call dbasis(ddx_data % params, ddx_data % constants, sji, basloc, dbasloc, vplm, vcos, vsin)
      alpha = zero
      do l = 0, ddx_data % params % lmax
        ind = l*l + l + 1
        f1 = (DI_rjin(l)*ddx_data % params % kappa)/ddx_data % constants % SI_ri(l,isph);
        f2 = SI_rjin(l)/ddx_data % constants % SI_ri(l,isph)

        do m = -l, l
          alpha = alpha + (f1*sji*basloc(ind+m) + &
                 & (f2/rjin)*dbasloc(:,ind+m))*Xe(ind+m,isph)
        end do
      end do
      xji = fsw(tji,ddx_data % params % se,ddx_data % params % eta)
      if (ddx_data % constants % fi(igrid,jsph).gt.one) then
        oji = xji/ddx_data % constants % fi(igrid,jsph)
      else
        oji = xji
      end if
      f1 = oji
      vb = vb + f1*alpha*Xadj_e(igrid,jsph)
      if (tji .gt. tlow) then
        ! Compute beta_jin, i.e., Eq.(57) Stamm.etal.18
        beta_ji = compute_beta(ddx_data, SI_rjin, rjin, isph, Xe(:,isph), basloc)
        if (ddx_data % constants % fi(igrid,jsph) .gt. one) then
          dj  = one/ddx_data % constants % fi(igrid,jsph)
          fac = dj*xji
          proc = .false.
          b    = zero
          do jk = ddx_data % constants % inl(jsph), ddx_data % constants % inl(jsph+1) - 1
            ksph = ddx_data % constants % nl(jk)
            vjk  = ddx_data % params % csph(:,jsph) + &
                 & ddx_data % params % rsph(jsph)*ddx_data % constants % cgrid(:,igrid) - &
                 & ddx_data % params % csph(:,ksph)
            rjkn = dnrm2(3, vjk, 1)
            tjk  = rjkn/ddx_data % params % rsph(ksph)
            ! Computation of modified spherical Bessel function values      
            call modified_spherical_bessel_first_kind(ddx_data % params % lmax, rjkn*ddx_data % params % kappa, SI_rjkn, DI_rjkn)

            if (ksph.ne.isph) then
              if (tjk .le. thigh) then
              proc = .true.
              sjk  = vjk/rjkn
              call ylmbas(sjk, rho, ctheta, stheta, cphi, sphi, &
                  & ddx_data % params % lmax, ddx_data % constants % vscales, basloc, vplm, &
                  & vcos, vsin)
              beta_jk  = compute_beta(ddx_data, SI_rjkn, rjkn, ksph, Xe(:,ksph), basloc)
              xjk = fsw(tjk, ddx_data % params % se, ddx_data % params % eta)
              b   = b + beta_jk*xjk
              end if
            end if
          end do
          if (proc) then
            g1 = dj*dj*dfsw(tji,ddx_data % params % se,ddx_data % params % eta)/ddx_data % params % rsph(isph)
            g2 = g1*Xadj_e(igrid,jsph)*b
            vc = vc + g2*sji
          end if
        else
          dj  = one
          fac = zero
        end if
        f2 = (one-fac)*dj*dfsw(tji,ddx_data % params % se,ddx_data % params % eta)/ddx_data % params % rsph(isph)
        vb = vb + f2*Xadj_e(igrid,jsph)*beta_ji*sji
      end if 
    end do
    force_e = force_e + ddx_data % constants % wgrid(igrid)*(vb - vc)
    end do
  return
  end subroutine fdokb_b_xe
  
  real(dp) function compute_beta(ddx_data, SI_rijn, rijn, jsph, Xe, basloc)
  implicit none
  type(ddx_type) :: ddx_data
  real(dp), dimension(0:ddx_data % params % lmax), intent(in) :: SI_rijn
  real(dp), dimension(ddx_data % constants % nbasis), intent(in) :: basloc
  real(dp), dimension(ddx_data % constants % nbasis), intent(in)   :: Xe
  integer, intent(in) :: jsph
  real(dp), intent(in) :: rijn

  integer :: l, m, ind
  real(dp)  :: ss, fac
  ss = zero

  ! loop over l
  do l = 0, ddx_data % params % lmax
    do m = -l, l
      ind = l*l + l + m + 1
      fac = SI_rijn(l)/ddx_data % constants % SI_ri(l,jsph)
      ss = ss + fac*basloc(ind)*Xe(ind)
    end do
  end do
     
  compute_beta = ss
  end function compute_beta
  
  !
  ! Subroutine to compute the derivative of U_i(x_in) and \bf{k}^j_l0(x_in)Y^j_l0m0(x_in)
  ! fdouky : Force Derivative of U_i^e(x_in), k_l0, and Y_l0m0
  ! @param[in] ddx_data     : Data Type
  ! @param[in] ksph         : Derivative with respect to x_k
  ! @param[in] Xr           : Solution of the Laplace problem
  ! @param[in] Xe           : Solution of the HSP problem
  ! @param[in] Xadj_r_sgrid : Solution of the Adjoint Laplace problem evaluated at the
  !                           grid
  ! @param[in] Xadj_e_sgrid : Solution of the Adjoint HSP problem evaluated at the grid
  ! @param[in] Xadj_r       : Adjoint solution of the Laplace problem
  ! @param[in] Xadj_e       : Adjoint solution of the HSP problem
  ! @param[inout] force     : Force
  ! @param[out] diff_re     : epsilon_1/epsilon_2 * l'/r_j[Xr]_jl'm' 
  !                         - (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
  subroutine fdouky(ddx_data, ksph, Xr, Xe, &
                          & Xadj_r_sgrid, Xadj_e_sgrid, &
                          & Xadj_r, Xadj_e, &
                          & force, &
                          & diff_re)
  implicit none
  type(ddx_type), intent(in) :: ddx_data
  integer, intent(in) :: ksph
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in) :: Xr, Xe
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph), intent(in) :: Xadj_r_sgrid,&
                                                                        & Xadj_e_sgrid
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in) :: Xadj_r, Xadj_e
  real(dp), dimension(3), intent(inout) :: force
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(out) :: diff_re
  real(dp), external :: dnrm2
  ! Local variable
  integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, icav
  ! val_dim3 : Intermediate value array of dimension 3
  real(dp), dimension(3) :: sij, vij, val_dim3
  ! val   : Intermediate variable to compute diff_ep
  ! f1    : Intermediate variable for derivative of coefY_der
  ! f2    : Intermediate variable for derivative of coefY_der
  real(dp) :: val, f1, f2
  ! phi_in : sum_{j=1}^N diff0_j * coefY_j
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph) :: phi_in
  ! diff_ep_dim3 : 3 dimensional couterpart of diff_ep
  real(dp), dimension(3, ddx_data % constants % ncav) :: diff_ep_dim3
  ! sum_dim3 : Storage of sum
  real(dp), dimension(3, ddx_data % constants % nbasis, ddx_data % params % nsph) :: sum_dim3
  ! coefY_der : Derivative of k_l0 and Y_l0m0
  real(dp), dimension(3, ddx_data % constants % ncav, &
      & ddx_data % constants % nbasis0, ddx_data % params % nsph) :: coefY_der
  ! Debug purpose
  ! These variables can be taken from the subroutine update_rhs
  ! diff0       : dot_product([PU_j]_l0m0^l'm', l'/r_j[Xr]_jl'm' -
  !                        (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm')

  real(dp), dimension(ddx_data % constants % nbasis0, ddx_data % params % nsph) :: diff0
  real(dp) :: termi, termk, rijn
  ! basloc : Y_lm(s_n)
  ! vplm   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % constants % nbasis):: basloc, vplm
  ! dbasloc : Derivative of Y_lm(s_n)
  real(dp),  dimension(3, ddx_data % constants % nbasis):: dbasloc
  ! vcos   : Argument to call ylmbas
  ! vsin   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % params % lmax+1):: vcos, vsin
  real(dp), dimension(0:ddx_data % params % lmax) :: SK_rijn, DK_rijn


  ! Setting initial values to zero
  SK_rijn = zero
  DK_rijn = zero
  coefY_der = zero

  !Compute coefY = C_ik*Y_lm(x_in)*\bar(k_l0^j(x_in))
  icav = 0
  do isph = 1, ddx_data % params % nsph
    do igrid = 1, ddx_data % params % ngrid
      if (ddx_data % constants % ui(igrid,isph).gt.zero) then
        icav = icav + 1
        ! Loop to compute Sijn
        do jsph = 1, ddx_data % params % nsph
          vij  = ddx_data % params % csph(:,isph) + &
                & ddx_data % params % rsph(isph)*ddx_data % constants % cgrid(:,igrid) - &
                & ddx_data % params % csph(:,jsph)
          rijn = sqrt(dot_product(vij,vij))
          sij = vij/rijn

          call modified_spherical_bessel_second_kind( &
              & ddx_data % constants % lmax0, &
                                                     & rijn*ddx_data % params % kappa, &
                                                     & SK_rijn, DK_rijn)
          call dbasis(ddx_data % params, ddx_data % constants, sij, basloc, dbasloc, vplm, vcos, vsin)

          do l0 = 0, ddx_data % constants % lmax0
            f1 = (DK_rijn(l0)*ddx_data % params % kappa)/ddx_data % constants % SK_ri(l0,jsph)
            f2 = SK_rijn(l0)/ddx_data % constants % SK_ri(l0,jsph)
            do m0 = -l0, l0
              ind0 = l0**2 + l0 + m0 + 1
              ! coefY_der : Derivative of Bessel function and spherical harmonic
              ! Non-Diagonal entries
              if ((ksph .eq. isph) .and. (isph .ne. jsph)) then
                coefY_der(:,icav,ind0,jsph) = ddx_data % constants % C_ik(l0,jsph)*(f1*sij*basloc(ind0) + &
                                             & (f2/rijn)*dbasloc(:,ind0))
              elseif ((ksph .eq. jsph) .and. (isph .ne. jsph)) then
                coefY_der(:,icav,ind0,jsph) = -ddx_data % constants % C_ik(l0,jsph)*(f1*sij*basloc(ind0)+ &
                                             & (f2/rijn)*dbasloc(:,ind0))
              else
                coefY_der(:,icav,ind0,jsph) = zero
              endif
            end do ! End of loop m0
          end do ! End of l0
        end do ! End of loop jsph
      end if
    end do ! End of loop igrid
  end do ! End of loop isph

  diff_re = zero
  ! Compute l'/r_j[Xr]_jl'm' -(i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
  do jsph = 1, ddx_data % params % nsph
    do l = 0, ddx_data % params % lmax
      do m = -l,l
        ind = l**2 + l + m + 1
        diff_re(ind,jsph) = (epsp/ddx_data % params % eps)*(l/ddx_data % params % rsph(jsph)) * &
              & Xr(ind,jsph) - ddx_data % constants % termimat(l,jsph)*Xe(ind,jsph)
      end do
    end do
  end do

  ! diff0 = Pchi * diff_re, linear scaling
  diff0 = zero
  do jsph = 1, ddx_data % params % nsph
    do ind0 = 1, ddx_data % constants % nbasis0
      diff0(ind0, jsph) = dot_product(diff_re(:,jsph), &
          & ddx_data % constants % Pchi(:,ind0, jsph))
    end do
  end do
  ! phi_in = diff0 * coefY
  ! Here, summation over j takes place
  phi_in = zero
  icav = zero
  do isph = 1, ddx_data % params % nsph
    do igrid = 1, ddx_data % params % ngrid
      if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
        ! Extrenal grid point
        icav = icav + 1
        val = zero
        val_dim3(:) = zero
        do jsph = 1, ddx_data % params % nsph 
          do ind0 = 1, ddx_data % constants % nbasis0
            val = val + diff0(ind0,jsph)*ddx_data % constants % coefY(icav,ind0,jsph)
            val_dim3(:) = val_dim3(:) + diff0(ind0,jsph)*coefY_der(:, icav, ind0, jsph)
          end do
        end do
        diff_ep_dim3(:, icav) = val_dim3(:)
      end if
    phi_in(igrid, isph) = val
    end do
  end do
  ! Computation of derivative of U_i^e(x_in)
  call fdoga(ddx_data % params, ddx_data % constants, ksph, Xadj_r_sgrid, phi_in, force)
  call fdoga(ddx_data % params, ddx_data % constants, ksph, Xadj_e_sgrid, phi_in, force)

  sum_dim3 = zero
  icav = zero
  do isph = 1, ddx_data % params % nsph
    do igrid =1, ddx_data % params % ngrid
      if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
        icav = icav + 1
        do ind = 1, ddx_data % constants % nbasis
          sum_dim3(:,ind,isph) = sum_dim3(:,ind,isph) + &
                                & ddx_data % constants % coefvec(igrid, ind, isph)*diff_ep_dim3(:,icav)
        end do
      end if
    end do
  end do
  ! Computation of derivative of \bf(k)_j^l0(x_in)\times Y^j_l0m0(x_in)
  do isph = 1, ddx_data % params % nsph
    do ind = 1, ddx_data % constants % nbasis
      force = force + sum_dim3(:, ind, isph)*(Xadj_r(ind, isph) + &
             & Xadj_e(ind, isph))
    end do
  end do

  end subroutine fdouky

  !
  ! Subroutine to compute the derivative of U_i(x_in) and \bf{k}^j_l0(x_in)Y^j_l0m0(x_in)
  ! in F0
  ! fdouky_f0 : Force Derivative of U_i^e(x_in), k_l0, and Y_l0m0 in F0
  ! @param[in] ddx_data     : Data Type
  ! @param[in] ksph         : Derivative with respect to x_k
  ! @param[in] sol_sgrid    : Solution of the Adjoint problem evaluated at the grid
  ! @param[in] sol_adj      : Adjoint solution
  ! @param[in] gradpsi      : Gradient of Psi_0
  ! @param[inout] force     : Force
  subroutine fdouky_f0(ddx_data, ksph,&
                          & sol_adj, sol_sgrid, gradpsi, &
                          & force)
  implicit none
  type(ddx_type), intent(in) :: ddx_data
  integer, intent(in) :: ksph
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in) :: sol_adj
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph), intent(in) :: sol_sgrid
  real(dp), dimension(3, ddx_data % constants % ncav), intent(in) :: gradpsi
  real(dp), dimension(3), intent(inout) :: force
  real(dp), external :: dnrm2
  ! Local variable
  integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, icav
  ! val_dim3 : Intermediate value array of dimension 3
  real(dp), dimension(3) :: sij, vij, val_dim3
  ! val     : Intermediate variable to compute diff_ep
  ! f1      : Intermediate variable for derivative of coefY_der
  ! f2      : Intermediate variable for derivative of coefY_der
  ! nderpsi : Derivative of psi on grid points
  real(dp) :: val, f1, f2, nderpsi, sum_int
  ! phi_in : sum_{j=1}^N diff0_j * coefY_j
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph) :: phi_in
  ! diff_ep_dim3 : 3 dimensional couterpart of diff_ep
  real(dp), dimension(3, ddx_data % constants % ncav) :: diff_ep_dim3
  ! sum_dim3 : Storage of sum
  real(dp), dimension(3, ddx_data % constants % nbasis, ddx_data % params % nsph) :: sum_dim3
  ! coefY_der : Derivative of k_l0 and Y_l0m0
  real(dp), dimension(3, ddx_data % constants % ncav, &
      & ddx_Data % constants % nbasis0, ddx_data % params % nsph) :: coefY_der
  ! Debug purpose
  ! These variables can be taken from the subroutine update_rhs
  ! diff0       : dot_product([PU_j]_l0m0^l'm', l'/r_j[Xr]_jl'm' -
  !                        (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm')

  real(dp), dimension(ddx_data % constants % nbasis0, ddx_data % params % nsph) :: diff0
  real(dp) :: termi, termk, rijn
  ! vplm   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % constants % nbasis):: basloc, vplm
  real(dp),  dimension(3, ddx_data % constants % nbasis):: dbasloc
  ! vcos   : Argument to call ylmbas
  ! vsin   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % params % lmax+1):: vcos, vsin
  real(dp), dimension(0:ddx_data % params % lmax) :: SK_rijn, DK_rijn
  ! sum_Sjin : \sum_j [S]_{jin} Eq.~(97) [QSM20.SISC]
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph) :: sum_Sjin
  ! c0 : \sum_{n=1}^N_g w_n U_j^{x_nj}\partial_n psi_0(x_nj)Y_{l0m0}(s_n)
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph) :: c0_d

  ! Setting initial values to zero
  SK_rijn = zero
  DK_rijn = zero
  coefY_der = zero
  c0_d = zero

  icav = zero
  do isph = 1, ddx_data % params % nsph
    do igrid= 1, ddx_data % params % ngrid
      if ( ddx_data % constants % ui(igrid,isph) .gt. zero ) then
        icav = icav + 1
        nderpsi = dot_product( gradpsi(:,icav),ddx_data % constants % cgrid(:,igrid) )
        c0_d(:, isph) = c0_d(:,isph) + &
                     & ddx_data % constants % wgrid(igrid)* &
                     & ddx_data % constants % ui(igrid,isph)*&
                     & nderpsi* &
                     & ddx_data % constants % vgrid(:,igrid)
      end if
    end do
  end do

  ! Compute [S]_{jin}
  icav = 0
  do isph = 1, ddx_data % params % nsph
    do igrid = 1, ddx_data % params % ngrid
      if (ddx_data % constants % ui(igrid,isph).gt.zero) then
        icav = icav + 1
        sum_int = zero
        ! Loop to compute Sijn
        do jsph = 1, ddx_data % params % nsph
          vij  = ddx_data % params % csph(:,isph) + &
                & ddx_data % params % rsph(isph)*ddx_data % constants % cgrid(:,igrid) - &
                & ddx_data % params % csph(:,jsph)
          rijn = sqrt(dot_product(vij,vij))
          sij = vij/rijn

          call modified_spherical_bessel_second_kind( &
              & ddx_data % constants % lmax0, &
                                                     & rijn*ddx_data % params % kappa, &
                                                     & SK_rijn, DK_rijn)
          call dbasis(ddx_data % params, ddx_data % constants, sij, basloc, dbasloc, vplm, vcos, vsin)

          do l0 = 0, ddx_data % constants % lmax0
            f1 = (DK_rijn(l0)*ddx_data % params % kappa)/ddx_data % constants % SK_ri(l0,jsph)
            f2 = SK_rijn(l0)/ddx_data % constants % SK_ri(l0,jsph)
            do m0 = -l0, l0
              ind0 = l0**2 + l0 + m0 + 1
              sum_int = sum_int + c0_d(ind0,jsph)*ddx_data % constants % coefY(icav, ind0, jsph)
              ! coefY_der : Derivative of Bessel function and spherical harmonic
              ! Non-Diagonal entries
              if ((ksph .eq. isph) .and. (isph .ne. jsph)) then
                coefY_der(:,icav,ind0,jsph) = ddx_data % constants % C_ik(l0,jsph)*(f1*sij*basloc(ind0) + &
                                             & (f2/rijn)*dbasloc(:,ind0))
              elseif ((ksph .eq. jsph) .and. (isph .ne. jsph)) then
                coefY_der(:,icav,ind0,jsph) = -ddx_data % constants % C_ik(l0,jsph)*(f1*sij*basloc(ind0)+ &
                                             & (f2/rijn)*dbasloc(:,ind0))
              else
                coefY_der(:,icav,ind0,jsph) = zero
              endif
            end do ! End of loop m0
          end do ! End of l0
        end do ! End of loop jsph
      end if
      sum_Sjin(igrid,isph) = -(epsp/ddx_data % params % eps)*sum_int
    end do ! End of loop igrid
  end do ! End of loop isph

  ! Computation of derivative of U_i^e(x_in)
  call fdoga(ddx_data % params, ddx_data % constants, ksph, sol_sgrid, sum_Sjin, force)

  ! Here, summation over j takes place
  icav = zero
  do isph = 1, ddx_data % params % nsph
    do igrid = 1, ddx_data % params % ngrid
      if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
        ! Extrenal grid point
        icav = icav + 1
        val_dim3(:) = zero
        do jsph = 1, ddx_data % params % nsph
          do ind0 = 1, ddx_data % constants % nbasis0
            val_dim3(:) = val_dim3(:) + c0_d(ind0,jsph)*coefY_der(:, icav, ind0, jsph)
          end do
        end do
        diff_ep_dim3(:, icav) = val_dim3(:)
      end if
    end do
  end do

  sum_dim3 = zero
  icav = zero
  do isph = 1, ddx_data % params % nsph
    do igrid =1, ddx_data % params % ngrid
      if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
        icav = icav + 1
        do ind = 1, ddx_data % constants % nbasis
          sum_dim3(:,ind,isph) = sum_dim3(:,ind,isph) + &
                                & -(epsp/ddx_data % params % eps)* &
                                & ddx_data % constants % coefvec(igrid, ind, isph)*diff_ep_dim3(:,icav)
        end do
      end if
    end do
  end do

  ! Computation of derivative of \bf(k)_j^l0(x_in)\times Y^j_l0m0(x_in)
  do isph = 1, ddx_data % params % nsph
    do ind = 1, ddx_data % constants % nbasis
      force = force + sum_dim3(:, ind, isph)*sol_adj(ind, isph)
    end do
  end do
  end subroutine fdouky_f0

  
  !
  ! Subroutine to calculate the third derivative term in C1_C2 matrix, namely the derivative of PU_i
  ! @param[in]  ddx_data     : Input data file
  ! @param[in]  ksph         : Derivative with respect to x_k
  ! @param[in]  Xr           : Solution of the Laplace problem
  ! @param[in]  Xe           : Solution of the HSP problem
  ! @param[in]  Xadj_r_sgrid : Adjoint Laplace solution evaluated at grid point
  ! @param[in]  Xadj_e_sgrid : Adjoint HSP solution evaluated at grid point
  ! @param[in]  diff_re      : l'/r_j[Xr]_jl'm' -(i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
  ! @param[out] force        : Force
  subroutine derivative_P(ddx_data, ksph, &
                          & Xr, Xe, &
                          & Xadj_r_sgrid, Xadj_e_sgrid, &
                          & diff_re, force)
  implicit none
  type(ddx_type) :: ddx_data
  integer, intent(in) :: ksph
  real(dp), dimension(ddx_data % constants % nbasis, ddx_data % params % nsph), intent(in) :: Xr, Xe, diff_re
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph), intent(in) :: Xadj_r_sgrid, Xadj_e_sgrid
  real(dp), dimension(3), intent(inout) :: force
  ! Local variable
  ! igrid0: Index for grid point n0
  integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, igrid0, icav
  ! term  : SK_rijn/SK_rj
  ! termi : DI_ri/SI_ri
  ! termk : DK_ri/SK_ri
  ! sum_int : Intermediate sum
  ! sum_r   : Intermediate sum for Laplace
  ! sum_e   : Intermediate sum for HSP
  real(dp) :: term, termi, termk, sum_int, sum_r, sum_e
  real(dp) :: rijn
  real(dp)  :: vij(3), sij(3)
  ! phi_n_r : Phi corresponding to Laplace problem
  ! phi_n_e : Phi corresponding to HSP problem
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph) :: phi_n_r, phi_n_e
  ! coefY_d : sum_{l0m0} C_ik*term*Y_l0m0^j(x_in)*Y_l0m0(s_n)
  real(dp), dimension(ddx_data % constants % ncav, ddx_data % params % ngrid, ddx_data % params % nsph) :: coefY_d
  ! diff_re_sgrid : diff_re evaluated at grid point
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph) :: diff_re_sgrid
  ! basloc : Y_lm(s_n)
  ! vplm   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % constants % nbasis):: basloc, vplm
  ! dbasloc : Derivative of Y_lm(s_n)
  real(dp),  dimension(3, ddx_data % constants % nbasis):: dbasloc
  ! vcos   : Argument to call ylmbas
  ! vsin   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % params % lmax+1):: vcos, vsin
  ! SK_rijn : Besssel function of first kind for rijn
  ! DK_rijn : Derivative of Besssel function of first kind for rijn
  real(dp), dimension(0:ddx_data % params % lmax) :: SK_rijn, DK_rijn

  ! Intial allocation of vectors
  sum_int = zero
  sum_r = zero
  sum_e = zero
  phi_n_r = zero
  phi_n_e = zero
  coefY_d = zero
  diff_re_sgrid = zero
  basloc = zero
  vplm = zero
  dbasloc = zero
  vcos = zero
  vsin = zero
  SK_rijn = zero
  DK_rijn = zero

  ! Compute  summation over l0, m0
  ! Loop over the sphers j
  do jsph = 1, ddx_data % params % nsph
    ! Loop over the grid points n0
    do igrid0 = 1, ddx_data % params % ngrid
      icav = zero
      ! Loop over spheres i
      do isph = 1, ddx_data % params % nsph
        ! Loop over grid points n
        do igrid = 1, ddx_data % params % ngrid
          ! Check for U_i^{eta}(x_in)
          if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
            icav = icav + 1
            vij  = ddx_data % params % csph(:,isph) + &
                   & ddx_data % params % rsph(isph)*ddx_data % constants % cgrid(:,igrid) - &
                   & ddx_data % params % csph(:,jsph)
            rijn = sqrt(dot_product(vij,vij))
            sij = vij/rijn

            call modified_spherical_bessel_second_kind( &
                & ddx_data % constants % lmax0, rijn*ddx_data % params % kappa,&
                                                     & SK_rijn, DK_rijn)
            call dbasis(ddx_data % params, ddx_data % constants, &
                & sij, basloc, dbasloc, vplm, vcos, vsin)
            sum_int = zero
            ! Loop over l0
            do l0 = 0, ddx_data % constants % lmax0
              term = SK_rijn(l0)/ddx_data % constants % SK_ri(l0,jsph)
              ! Loop over m0
              do m0 = -l0,l0
                ind0 = l0**2 + l0 + m0 + 1
                sum_int = sum_int + ddx_data % constants % C_ik(l0, jsph) *term*basloc(ind0)&
                           & *ddx_data % constants % vgrid(ind0,igrid0)
              end do ! End of loop m0
            end do! End of loop l0
            coefY_d(icav, igrid0, jsph) = sum_int
          end if
        end do ! End of loop igrid
      end do! End of loop isph
    end do ! End of loop igrid0
  end do ! End of loop jsph

  ! Compute phi_in
  ! Loop over spheres j
  do jsph = 1, ddx_data % params % nsph
    ! Loop over grid points n0
    do igrid0 = 1, ddx_data % params % ngrid
      icav = zero
      sum_r = zero
      sum_e = zero
      ! Loop over sphers i
      do isph = 1, ddx_data % params % nsph
        ! Loop over grid points n
        do igrid = 1, ddx_data % params % ngrid
          if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
            icav = icav + 1
            sum_r = sum_r + coefY_d(icav, igrid0, jsph)*Xadj_r_sgrid(igrid, isph) &
                    & * ddx_data % constants % wgrid(igrid)*ddx_data % constants % ui(igrid, isph)
            sum_e = sum_e + coefY_d(icav, igrid0, jsph)*Xadj_e_sgrid(igrid, isph) &
                    & * ddx_data % constants % wgrid(igrid)*ddx_data % constants % ui(igrid, isph)
          end if
        end do
      end do
      phi_n_r(igrid0, jsph) = sum_r
      phi_n_e(igrid0, jsph) = sum_e
    end do! End of loop j
  end do ! End of loop igrid0


  call dgemm('T', 'N', ddx_data % params % ngrid, ddx_data % params % nsph, &
            & ddx_data % constants % nbasis, one, ddx_data % constants % vgrid, ddx_data % constants % vgrid_nbasis, &
            & diff_re , ddx_data % constants % nbasis, zero, diff_re_sgrid, &
            & ddx_data % params % ngrid)
  call fdoga(ddx_data % params, ddx_data % constants, ksph, diff_re_sgrid, phi_n_r, force)
  call fdoga(ddx_data % params, ddx_data % constants, ksph, diff_re_sgrid, phi_n_e, force)
  end subroutine derivative_P

  !
  ! Subroutine to calculate the derivative of C_{0l0m0}^j
  ! @param[in]  ddx_data           : Input data file
  ! @param[in]  ksph               : Derivative with respect to x_k
  ! @param[in]  sol_sgrid          : Solution of the Adjoint problem evaluated at the grid
  ! @param[in]  gradpsi            : Gradient of Psi_0
  ! @param[in]  normal_hessian_cav : Normal of the Hessian evaluated at cavity points
  ! @param[in]  icav_g             : Index of outside cavity point
  ! @param[out] force              : Force corresponding to HSP problem
  subroutine fdoco(ddx_data, ksph, sol_sgrid, gradpsi, normal_hessian_cav, icav_g, force)
  implicit none
  type(ddx_type) :: ddx_data
  integer, intent(in) :: ksph
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph), intent(in) :: sol_sgrid
  real(dp), dimension(3, ddx_data % constants % ncav), intent(in) :: gradpsi
  real(dp), dimension(3, ddx_data % constants % ncav), intent(in) :: normal_hessian_cav
  integer, intent(inout) :: icav_g
  real(dp), dimension(3), intent(inout) :: force
  ! Local variable
  integer :: isph, jsph, igrid, l, m, ind, l0, m0, ind0, igrid0, icav
  ! term  : SK_rijn/SK_rj
  ! termi : DI_ri/SI_ri
  ! termk : DK_ri/SK_ri
  ! sum_int : Intermediate sum
  ! hessian_contribution :
  real(dp) :: term, termi, termk, sum_int, hessian_contribution(3), nderpsi
  real(dp) :: rijn
  real(dp)  :: vij(3), sij(3)
  ! phi_n : Phi corresponding to Laplace problem
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph) :: phi_n
  ! coefY_d : sum_{l0m0} C_ik*term*Y_l0m0^j(x_in)*Y_l0m0(s_n)
  real(dp), dimension(ddx_data % constants % ncav, ddx_data % params % ngrid, ddx_data % params % nsph) :: coefY_d
  ! gradpsi_grid : gradpsi evaluated at grid point
  real(dp), dimension(ddx_data % params % ngrid, ddx_data % params % nsph) :: gradpsi_grid
  ! vplm   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % constants % nbasis):: basloc, vplm
  real(dp),  dimension(3, ddx_data % constants % nbasis):: dbasloc
  ! vcos   : Argument to call ylmbas
  ! vsin   : Argument to call ylmbas
  real(dp),  dimension(ddx_data % params % lmax+1):: vcos, vsin
  real(dp), dimension(0:ddx_data % params % lmax) :: SK_rijn, DK_rijn

  ! Intial allocation of vectors
  sum_int = zero
  phi_n = zero
  coefY_d = zero
  gradpsi_grid = zero
  basloc = zero
  vplm = zero
  dbasloc = zero
  vcos = zero
  vsin = zero
  SK_rijn = zero
  DK_rijn = zero

  ! Compute  summation over l0, m0
  ! Loop over the sphers j
  do jsph = 1, ddx_data % params % nsph
    ! Loop over the grid points n0
    do igrid0 = 1, ddx_data % params % ngrid
      icav = zero
      ! Loop over spheres i
      do isph = 1, ddx_data % params % nsph
        ! Loop over grid points n
        do igrid = 1, ddx_data % params % ngrid
          ! Check for U_i^{eta}(x_in)
          if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
            icav = icav + 1
            vij  = ddx_data % params % csph(:,isph) + &
                   & ddx_data % params % rsph(isph)*ddx_data % constants % cgrid(:,igrid) - &
                   & ddx_data % params % csph(:,jsph)
            rijn = sqrt(dot_product(vij,vij))
            sij = vij/rijn

            call modified_spherical_bessel_second_kind( &
                & ddx_data % constants % lmax0, rijn*ddx_data % params % kappa,&
                                                     & SK_rijn, DK_rijn)
            call dbasis(ddx_data % params, ddx_data % constants, &
                & sij, basloc, dbasloc, vplm, vcos, vsin)
            sum_int = zero
            ! Loop over l0
            do l0 = 0, ddx_data % constants % lmax0
              term = SK_rijn(l0)/ddx_data % constants % SK_ri(l0,jsph)
              ! Loop over m0
              do m0 = -l0,l0
                ind0 = l0**2 + l0 + m0 + 1
                sum_int = sum_int + ddx_data % constants % C_ik(l0, jsph) &
                           & *term*basloc(ind0) &
                           & *ddx_data % constants % vgrid(ind0,igrid0)
              end do ! End of loop m0
            end do! End of loop l0
            coefY_d(icav, igrid0, jsph) = sum_int
          end if
        end do ! End of loop igrid
      end do! End of loop isph
    end do ! End of loop igrid0
  end do ! End of loop jsph

  ! Compute phi_in
  ! Loop over spheres j
  do jsph = 1, ddx_data % params % nsph
    ! Loop over grid points n0
    do igrid0 = 1, ddx_data % params % ngrid
      icav = zero
      sum_int = zero
      ! Loop over sphers i
      do isph = 1, ddx_data % params % nsph
        ! Loop over grid points n
        do igrid = 1, ddx_data % params % ngrid
          if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
            icav = icav + 1
            sum_int = sum_int + coefY_d(icav, igrid0, jsph)*sol_sgrid(igrid, isph) &
                    & * ddx_data % constants % wgrid(igrid)*ddx_data % constants % ui(igrid, isph)
          end if
        end do
      end do
      phi_n(igrid0, jsph) = -(epsp/ddx_data % params % eps)*sum_int
    end do! End of loop j
  end do ! End of loop igrid

  icav = zero
  do isph = 1, ddx_data % params % nsph
    do igrid = 1, ddx_data % params % ngrid
      if(ddx_data % constants % ui(igrid, isph) .gt. zero) then
        icav = icav + 1
        nderpsi = dot_product( gradpsi(:,icav),ddx_data % constants % cgrid(:,igrid) )
        gradpsi_grid(igrid, isph) = nderpsi
      end if
    end do ! End of loop igrid
  end do ! End of loop i

  call fdoga(ddx_data % params, ddx_data % constants, ksph, gradpsi_grid, phi_n, force)

  ! Compute the Hessian contributions
  do igrid = 1, ddx_data % params % ngrid
    if(ddx_data % constants % ui(igrid, ksph) .gt. zero) then
      icav_g = icav_g + 1
      force = force + ddx_data % constants % wgrid(igrid)*ddx_data % constants % ui(igrid, ksph)*&
                       & phi_n(igrid, ksph)*normal_hessian_cav(:, icav_g)
    end if
  end do

  call fdops(ddx_data, ksph, phi_n, force)

  end subroutine fdoco

  !
  ! fdops : Force derivative of potential at spheres
  ! @param[in]  ddx_data : Input data file
  ! @param[in]  phi_n    : phi_n^j
  ! @param[in]  ksph     : Derivative with respect to x_k
  ! @param[out] force    : Force
  subroutine fdops(ddx_data, ksph, phi_n, force)
  type(ddx_type), intent(in) :: ddx_data
  integer, intent(in) :: ksph
  real(dp),  dimension(ddx_data % params % ngrid, ddx_data % params % nsph), intent(in)    :: phi_n
  real(dp),  dimension(3), intent(inout) :: force
  !
  ! Local variables
  integer :: jsph, i,j, igrid
  ! sum_int            : Intermediate sum
  ! normal_hessian_cav : Normal derivative of Hessian of psi
  ! vij_vijT           : vij*vij^T
  real(dp)  :: sum_int(3), vij(3), rijn, normal_hessian_cav(3), vij_vijT(3,3)
  ! identity_matrix : Identity matrix of size 3x3
  ! hessianv        : Hessian of psi evaluated by centers
  real(dp) :: identity_matrix(3,3), hessianv(3,3)
  real(dp), external :: dnrm2
  ! Create Identity matrix
  identity_matrix = zero
  do i = 1, 3
    identity_matrix(i,i) = one
  end do

  sum_int = zero
  ! Loop over spheres
  do jsph = 1, ddx_data % params % nsph
    ! Loop over grid points
    do igrid = 1, ddx_data % params % ngrid
      if(ddx_data % constants % ui(igrid, jsph) .ne. zero) then
        vij_vijT = zero
        hessianv = zero
        normal_hessian_cav = zero
        vij = zero
        vij = ddx_data % params % csph(:, jsph) + &
           & ddx_data % params % rsph(jsph)*ddx_data % constants % cgrid(:, igrid) - &
           & ddx_data % params % csph(:, ksph)
        do i = 1,3
          do j = 1,3
            vij_vijT(i,j) = vij(i)*vij(j)
          end do
        end do
        rijn = dnrm2(3, vij, 1)
        hessianv = 3*vij_vijT/(rijn**5)- &
                 & identity_matrix/(rijn**3)
        do i = 1, 3
          normal_hessian_cav = normal_hessian_cav + &
                             & hessianv(:,i)*ddx_data % constants % cgrid(i,igrid)
        end do
        sum_int = sum_int + &
                 & ddx_data % constants % ui(igrid, jsph)* &
                 & ddx_data % constants % wgrid(igrid)* &
                 & phi_n(igrid, jsph)* &
                 & normal_hessian_cav
      end if
    end do
  end do
  force = force - ddx_data % params % charge(ksph)*sum_int
  end subroutine fdops

  subroutine ddlpb_free(ddx_data)
  implicit none
  type(ddx_type), intent(inout) :: ddx_data
  integer :: istatus
  if(allocated(ddx_data % constants % SI_ri)) then
    deallocate(ddx_data % constants % SI_ri, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % DI_ri)) then
    deallocate(ddx_data % constants % DI_ri, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % SK_ri)) then
    deallocate(ddx_data % constants % SK_ri, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % DK_ri)) then
    deallocate(ddx_data % constants % DK_ri, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if


  if(allocated(ddx_data % constants % diff_ep_adj)) then
    deallocate(ddx_data % constants % diff_ep_adj, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if


  if(allocated(ddx_data % constants % termimat)) then
    deallocate(ddx_data % constants % termimat, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % coefY)) then
    deallocate(ddx_data % constants % coefY, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % C_ik)) then
    deallocate(ddx_data % constants % C_ik, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % coefvec)) then
    deallocate(ddx_data % constants % coefvec, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  if(allocated(ddx_data % constants % Pchi)) then
    deallocate(ddx_data % constants % Pchi, stat=istatus)
    if(istatus .ne. zero) then
      write(*,*) 'ddlpb_free: [1] deallocation failed'
      stop 1
    end if
  end if

  end subroutine ddlpb_free

end module ddx_lpb
