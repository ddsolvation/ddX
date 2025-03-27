!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!!
!! @version 1.0.0
!! @author Abhinav Jha
!! @date 2022-04-01

!> Core routines and parameters specific to gradients
module ddx_gradients
! Get the core-routines
use ddx_core
!
contains

!> Compute the gradients of the ddCOSMO matrix
subroutine contract_grad_L(params, constants, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx)
type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
      integer,                         intent(in)    :: isph
      real(dp),  dimension(constants % nbasis, params % nsph), intent(in)    :: sigma
      real(dp),  dimension(params % ngrid, params % nsph),       intent(in)    :: xi
      real(dp),  dimension(constants % nbasis),      intent(inout) :: basloc, vplm
      real(dp),  dimension(3, constants % nbasis),    intent(inout) :: dbsloc
      real(dp),  dimension(params % lmax+1),      intent(inout) :: vcos, vsin
      real(dp),  dimension(3),           intent(inout) :: fx

      call contract_gradi_Lik(params, constants, isph, sigma, xi(:, isph), basloc, dbsloc, vplm, vcos, vsin, fx )
      call contract_gradi_Lji(params, constants, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx )

end subroutine contract_grad_L

!> Contribution to the gradients of the ddCOSMO matrix
subroutine contract_gradi_Lik(params, constants, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx )
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
      integer,                         intent(in)    :: isph
      real(dp),  dimension(constants % nbasis, params % nsph), intent(in)    :: sigma
      real(dp),  dimension(params % ngrid),       intent(in)    :: xi
      real(dp),  dimension(constants % nbasis),      intent(inout) :: basloc, vplm
      real(dp),  dimension(3, constants % nbasis),    intent(inout) :: dbsloc
      real(dp),  dimension(params % lmax+1),      intent(inout) :: vcos, vsin
      real(dp),  dimension(3),           intent(inout) :: fx
      integer :: ig, ij, jsph, l, ind, m
      real(dp)  :: vvij, tij, xij, oij, t, fac, fl, f1, f2, f3, beta, tlow, thigh
      real(dp)  :: vij(3), sij(3), alp(3), va(3)
      real(dp), external :: dnrm2
      tlow  = one - pt5*(one - params % se)*params % eta
      thigh = one + pt5*(one + params % se)*params % eta

      do ig = 1, params % ngrid
        va = zero
        do ij = constants % inl(isph), constants % inl(isph+1) - 1
          jsph = constants % nl(ij)
          vij  = params % csph(:,isph) + &
              & params % rsph(isph)*constants % cgrid(:,ig) - &
              & params % csph(:,jsph)
          vvij = dnrm2(3, vij, 1)
          tij  = vvij/params % rsph(jsph)
          if (tij.ge.thigh) cycle
          if (tij.ne.zero) then
              sij = vij/vvij
          else
              sij = one
          end if
          call dbasis(params, constants, sij, basloc, dbsloc, vplm, vcos, vsin)
          alp  = zero
          t    = one
          do l = 1, params % lmax
            ind = l*l + l + 1
            fl  = dble(l)
            fac = t/(constants % vscales(ind)**2)
            do m = -l, l
              f2 = fac*sigma(ind+m,jsph)
              f1 = f2*fl*basloc(ind+m)
              alp(:) = alp(:) + f1*sij(:) + f2*dbsloc(:,ind+m)
            end do
            t = t*tij
          end do
          beta = intmlp(params, constants, tij,sigma(:,jsph),basloc)
          xij = fsw(tij, params % se, params % eta)
          if (constants % fi(ig,isph).gt.one) then
            oij = xij/constants % fi(ig,isph)
            f2  = -oij/constants % fi(ig,isph)
          else
            oij = xij
            f2  = zero
          end if
          f1 = oij/params % rsph(jsph)
          va(:) = va(:) + f1*alp(:) + beta*f2*constants % zi(:,ig,isph)
          if (tij .gt. tlow) then
            f3 = beta*dfsw(tij,params % se,params % eta)/params % rsph(jsph)
            if (constants % fi(ig,isph).gt.one) f3 = f3/constants % fi(ig,isph)
            va(:) = va(:) + f3*sij(:)
          end if
        end do
        fx = fx - constants % wgrid(ig)*xi(ig)*va(:)
      end do
end subroutine contract_gradi_Lik

!> Contribution to the gradients of the ddCOSMO matrix
subroutine contract_gradi_Lji(params, constants, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx )
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
      integer,                         intent(in)    :: isph
      real(dp),  dimension(constants % nbasis, params % nsph), intent(in)    :: sigma
      real(dp),  dimension(params % ngrid, params % nsph),  intent(in)    :: xi
      real(dp),  dimension(constants % nbasis),      intent(inout) :: basloc, vplm
      real(dp),  dimension(3, constants % nbasis),    intent(inout) :: dbsloc
      real(dp),  dimension(params % lmax+1),      intent(inout) :: vcos, vsin
      real(dp),  dimension(3),           intent(inout) :: fx

      integer :: ig, ji, jsph, l, ind, m, jk, ksph
      logical :: proc
      real(dp)  :: vvji, tji, xji, oji, t, fac, fl, f1, f2, beta, di, tlow, thigh
      real(dp)  :: b, g1, g2, vvjk, tjk, xjk
      real(dp)  :: vji(3), sji(3), alp(3), vb(3), vjk(3), sjk(3), vc(3)
      real(dp) :: rho, ctheta, stheta, cphi, sphi
      real(dp), external :: dnrm2

      tlow  = one - pt5*(one - params % se)*params % eta
      thigh = one + pt5*(one + params % se)*params % eta

      do ig = 1, params % ngrid
        vb = zero
        vc = zero
        do ji = constants % inl(isph), constants % inl(isph+1) - 1
          jsph = constants % nl(ji)
          vji  = params % csph(:,jsph) + &
              & params % rsph(jsph)*constants % cgrid(:,ig) - &
              & params % csph(:,isph)
          vvji = dnrm2(3, vji, 1)
          tji  = vvji/params % rsph(isph)
          if (tji.gt.thigh) cycle
          if (tji.ne.zero) then
              sji = vji/vvji
          else
              sji = one
          end if
          call dbasis(params, constants, sji, basloc, dbsloc, vplm, vcos, vsin)
          alp = zero
          t   = one
          do l = 1, params % lmax
            ind = l*l + l + 1
            fl  = dble(l)
            fac = t/(constants % vscales(ind)**2)
            do m = -l, l
              f2 = fac*sigma(ind+m,isph)
              f1 = f2*fl*basloc(ind+m)
              alp = alp + f1*sji + f2*dbsloc(:,ind+m)
            end do
            t = t*tji
          end do
          xji = fsw(tji, params % se, params % eta)
          if (constants % fi(ig,jsph).gt.one) then
            oji = xji/constants % fi(ig,jsph)
          else
            oji = xji
          end if
          f1 = oji/params % rsph(isph)
          vb = vb + f1*alp*xi(ig,jsph)
          if (tji .gt. tlow) then
            beta = intmlp(params, constants, tji, sigma(:,isph), basloc)
            if (constants % fi(ig,jsph) .gt. one) then
              di  = one/constants % fi(ig,jsph)
              fac = di*xji
              proc = .false.
              b    = zero
              do jk = constants % inl(jsph), constants % inl(jsph+1) - 1
                ksph = constants % nl(jk)
                vjk  = params % csph(:,jsph) + &
                    & params % rsph(jsph)*constants % cgrid(:,ig) - &
                    & params % csph(:,ksph)
                !vvjk = sqrt(dot_product(vjk,vjk))
                vvjk = dnrm2(3, vjk, 1)
                tjk  = vvjk/params % rsph(ksph)
                if (ksph.ne.isph) then
                  if (tjk .le. thigh) then
                    proc = .true.
                    sjk  = vjk/vvjk
                    !call ylmbas(sjk,basloc,vplm,vcos,vsin)
                    call ylmbas(sjk, rho, ctheta, stheta, cphi, sphi, &
                        & params % lmax, constants % vscales, basloc, vplm, &
                        & vcos, vsin)
                    g1  = intmlp(params, constants, tjk, sigma(:,ksph), basloc)
                    xjk = fsw(tjk, params % se, params % eta)
                    b   = b + g1*xjk
                  end if
                end if
              end do
              if (proc) then
                g1 = di*di*dfsw(tji, params % se, params % eta)/params % rsph(isph)
                g2 = g1*xi(ig,jsph)*b
                vc = vc + g2*sji
              end if
            else
              di  = one
              fac = zero
            end if
            f2 = (one-fac)*di*dfsw(tji, params % se, params % eta)/params % rsph(isph)
            vb = vb + f2*xi(ig,jsph)*beta*sji
          end if
        end do
        fx = fx + constants % wgrid(ig)*(vb - vc)
      end do
end subroutine contract_gradi_Lji

!> Gradient of the characteristic function U
subroutine contract_grad_U(params, constants, isph, xi, phi, fx )
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
      integer,                        intent(in)    :: isph
      real(dp),  dimension(params % ngrid, params % nsph), intent(in)    :: xi, phi
      real(dp),  dimension(3),          intent(inout) :: fx
      integer :: ig, ji, jsph
      real(dp)  :: vvji, tji, fac, swthr
      real(dp)  :: alp(3), vji(3), sji(3)
      real(dp), external :: dnrm2
      do ig = 1, params % ngrid
        alp = zero
        if (constants % ui(ig,isph) .gt. zero .and. constants % ui(ig,isph).lt.one) then
          alp = alp + phi(ig,isph)*xi(ig,isph)*constants % zi(:,ig,isph)
        end if
        do ji = constants % inl(isph), constants % inl(isph+1) - 1
          jsph  = constants % nl(ji)
          vji   = params % csph(:,jsph) + &
              & params % rsph(jsph)*constants % cgrid(:,ig) - &
              & params % csph(:,isph)
          !vvji  = sqrt(dot_product(vji,vji))
          vvji = dnrm2(3, vji, 1)
          tji   = vvji/params % rsph(isph)
          swthr = one + (params % se + 1.d0)*params % eta / 2.d0
          if (tji.lt.swthr .and. tji.gt.swthr-params % eta .and. constants % ui(ig,jsph).gt.zero) then
            sji = vji/vvji
            fac = - dfsw(tji, params % se, params % eta)/params % rsph(isph)
            alp = alp + fac*phi(ig,jsph)*xi(ig,jsph)*sji
          end if
        end do
        fx = fx - constants % wgrid(ig)*alp
      end do
end subroutine contract_grad_U

!> Subroutine to compute contraction of B matrix
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  isph       : Index of sphere
!! @param[in]  Xe         : Solution vector Xe
!! @param[in]  Xadj_e     : Adjoint solution on evaluated on grid points Xadj_e_sgrid
!! @param[out] force      : Force of adjoint part

subroutine contract_grad_B(params, constants, isph, Xe, Xadj_e, force)
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: isph
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xe
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: Xadj_e
    real(dp), dimension(3), intent(inout) :: force

    call contract_gradi_Bik(params, constants, isph, Xe(:,:), &
        & Xadj_e(:, isph), force)
    call contract_gradi_Bji(params, constants, isph, Xe(:,:), Xadj_e, force)

end subroutine contract_grad_B

!> Subroutine to compute contraction of C matrix
!!
!! @param[in]  params      : Input parameter file
!! @param[in]  constants   : Input constants file
!! @param[in]  workspace   : Input workspace
!! @param[in] Xr           : Solution of the Laplace problem
!! @param[in] Xe           : Solution of the HSP problem
!! @param[in] Xadj_r_sgrid : Solution of the Adjoint Laplace problem evaluated at the
!!                           grid
!! @param[in] Xadj_e_sgrid : Solution of the Adjoint HSP problem evaluated at the grid
!! @param[in] Xadj_r       : Adjoint solution of the Laplace problem
!! @param[in] Xadj_e       : Adjoint solution of the HSP problem
!! @param[inout] force     : Force
!! @param[out] diff_re     : epsilon_1/epsilon_2 * l'/r_j[Xr]_jl'm'
!!                         - (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
!! @param[inout] ddx_error: ddX error
subroutine contract_grad_C(params, constants, workspace, Xr, Xe, Xadj_r_sgrid, &
    & Xadj_e_sgrid, Xadj_r, Xadj_e, force, diff_re, ddx_error)
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xr, Xe
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: Xadj_r_sgrid, Xadj_e_sgrid
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xadj_r, Xadj_e
    real(dp), dimension(3, params % nsph), intent(inout) :: force
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: diff_re
    type(ddx_error_type), intent(inout) :: ddx_error

    call contract_grad_C_worker2(params, constants, workspace, Xr, Xe, Xadj_r_sgrid, &
        & Xadj_e_sgrid, Xadj_r, Xadj_e, force, diff_re, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "contract_grad_C_worker2 returned an error, exiting")
        return
    end if
    call contract_grad_C_worker1(params, constants, workspace, Xadj_r_sgrid, &
        & Xadj_e_sgrid, diff_re, force, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "contract_grad_C_worker1 returned an error, exiting")
        return
    end if

end subroutine contract_grad_C


!> Subroutine to compute contraction of F0
!! @param[in]  params             : Input parameter file
!! @param[in]  constants          : Input constants file
!! @param[in]  workspace          : Input workspace
!! @param[in] sol_sgrid           : Solution of the Adjoint problem evaluated at the grid
!! @param[in] sol_adj             : Adjoint solution
!! @param[in]  normal_hessian_cav : Normal of the Hessian evaluated at cavity points
!! @param[in]  icav_g             : Index of outside cavity point
!! @param[out] force              : Force
!! @param[inout] ddx_error: ddX error
subroutine contract_grad_f(params, constants, workspace, sol_adj, sol_sgrid, &
    & gradpsi, normal_hessian_cav, force, state, ddx_error)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: sol_adj
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: sol_sgrid
    real(dp), dimension(3, constants % ncav), intent(in) :: gradpsi
    real(dp), dimension(3, constants % ncav), intent(in) :: normal_hessian_cav
    real(dp), dimension(3, params % nsph), intent(inout) :: force
    type(ddx_error_type), intent(inout) :: ddx_error

    call contract_grad_f_worker1(params, constants, workspace, sol_adj, sol_sgrid, &
        & gradpsi, force, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "contract_grad_f_worker1 returned an error, exiting")
        return
    end if

    call contract_grad_f_worker2(params, constants, gradpsi, &
        & normal_hessian_cav, force, state, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "contract_grad_f_worker2 returned an error, exiting")
        return
    end if

end subroutine contract_grad_f

!> Subroutine to compute K^A counterpart for the HSP equation. Similar to contract_gradi_Lik.
!!
!! @param[in]  params    : Input parameter file
!! @param[in]  constants : Input constants file
!! @param[in]  isph      : Index of sphere
!! @param[in]  Xe        : Solution vector Xe
!! @param[in]  Xadj_e    : Adjoint solution on evaluated on grid points Xadj_e_sgrid
!! @param[out] force     : Force
subroutine contract_gradi_Bik(params, constants, isph, Xe, Xadj_e, force)
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: isph
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xe
    real(dp), dimension(params % ngrid), intent(in) :: Xadj_e
    real(dp), dimension(3), intent(inout) :: force

    ! Local Variables
    integer :: igrid, ineigh, jsph
    real(dp), dimension(0:params % lmax) :: SI_rijn
    real(dp), dimension(0:params % lmax) :: DI_rijn
    ! beta   : Eq.(53) Stamm.etal.18
    ! tlow   : Lower bound for switch region
    ! thigh  : Upper bound for switch region
    ! f1     : First factor in alpha computation
    ! f2     : Second factor in alpha computation
    real(dp)  :: rijn, tij, beta, tlow, thigh, xij, oij, f1, f2, f3
    ! alpha : Eq.(52) Stamm.etal.18
    ! va    : Eq.(54) Stamm.etal.18
    real(dp)  :: vij(3), sij(3), alpha(3), va(3), rj, vtij(3)
    real(dp), external :: dnrm2
    real(dp) :: work(params % lmax+1)
    complex(dp) :: work_complex(params % lmax+1)

    SI_rijn = 0
    DI_rijn = 0

    tlow  = one - pt5*(one - params % se)*params % eta
    thigh = one + pt5*(one + params % se)*params % eta

    ! Loop over grid points
    do igrid = 1, params % ngrid
        va = zero
        do ineigh = constants % inl(isph), constants % inl(isph+1) - 1
            jsph = constants % nl(ineigh)
            vij  = params % csph(:,isph) + &
                & params % rsph(isph)*constants % cgrid(:,igrid) - &
                & params % csph(:,jsph)
            rijn = dnrm2(3, vij, 1)
            tij  = rijn/params % rsph(jsph)
            rj = params % rsph(jsph)
            if (tij.ge.thigh) cycle
            ! Computation of modified spherical Bessel function values
            if (tij.ne.zero) then
                sij = vij/rijn
            else
                sij = one
            end if
            vtij = vij*params % kappa
            call fmm_l2p_bessel_grad(vtij, params % rsph(jsph)*params % kappa, &
                & params % lmax, constants % vscales, params % kappa, &
                & Xe(:, jsph), zero, alpha)
            call fmm_l2p_bessel_work(vtij, params % lmax, constants % vscales, &
                & constants % SI_ri(:, jsph), one, Xe(:, jsph), zero, beta, &
                & work_complex, work)
            xij = fsw(tij, params % se, params % eta)
            if (constants % fi(igrid,isph).gt.one) then
                oij = xij/constants % fi(igrid,isph)
                f2  = -oij/constants % fi(igrid,isph)
            else
                oij = xij
                f2  = zero
            end if
            f1 = oij
            va(:) = va(:) + f1*alpha(:) + beta*f2*constants % zi(:,igrid,isph)
            if (tij .gt. tlow) then
                f3 = beta*dfsw(tij,params % se,params % eta)/params % rsph(jsph)
                if (constants % fi(igrid,isph).gt.one) &
                    & f3 = f3/constants % fi(igrid,isph)
                va(:) = va(:) + f3*sij(:)
            end if
        end do
        force = force - constants % wgrid(igrid)*Xadj_e(igrid)*va(:)
    end do
end subroutine contract_gradi_Bik


!> Subroutine to compute K^A+K^C counterpart for the HSP equation. Similar to contract_gradi_Lji.
!!
!! @param[in]  params     : Input parameter file
!! @param[in]  constants  : Input constants file
!! @param[in]  workspace  : Input workspace
!! @param[in]  isph       : Index of sphere
!! @param[in]  Xe         : Solution vector Xe
!! @param[in]  Xadj_e     : Adjoint solution on evaluated on grid points Xadj_e_sgrid
!! @param[out] force      : Force of adjoint part
subroutine contract_gradi_Bji(params, constants, isph, Xe, Xadj_e, force)
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in)    :: isph
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xe
    real(dp), dimension(params % ngrid, params % nsph),  intent(in) :: Xadj_e
    real(dp), dimension(3), intent(inout) :: force

    ! Local Variables
    ! jk : Row pointer over kth row
    integer :: igrid, jsph, ksph, ineigh, jk
    real(dp), dimension(0:params % lmax) :: SI_rjin, SI_rjkn
    real(dp), dimension(0:params % lmax) :: DI_rjin, DI_rjkn

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
    real(dp)  :: vji(3), sji(3), vjk(3), alpha(3), vb(3), vc(3), &
        & vtji(3), vtjk(3)
    ! rho    : Argument for ylmbas
    ! ctheta : Argument for ylmbas
    ! stheta : Argument for ylmbas
    ! cphi   : Argument for ylmbas
    ! sphi   : Argument for ylmbas
    real(dp) :: ri

    real(dp), external :: dnrm2
    real(dp) :: work(params % lmax+1)
    complex(dp) :: work_complex(params % lmax+1)

    SI_rjin = 0
    DI_rjin = 0
    SI_rjkn = 0
    DI_rjkn = 0

    tlow  = one - pt5*(one - params % se)*params % eta
    thigh = one + pt5*(one + params % se)*params % eta

    do igrid = 1, params % ngrid
        vb = zero
        vc = zero
        do ineigh = constants % inl(isph), constants % inl(isph+1) - 1
            jsph = constants % nl(ineigh)
            vji  = params % csph(:,jsph) + &
                    & params % rsph(jsph)*constants % cgrid(:,igrid) - &
                    & params % csph(:,isph)
            rjin = dnrm2(3, vji, 1)
            ri = params % rsph(isph)
            tji  = rjin/ri
            if (tji.gt.thigh) cycle
            if (tji.ne.zero) then
                sji = vji/rjin
            else
                sji = one
            end if
            vtji = vji*params % kappa
            call fmm_l2p_bessel_grad(vtji, params % rsph(isph)*params % kappa, &
                & params % lmax, constants % vscales, params % kappa, &
                & Xe(:, isph), zero, alpha)
            xji = fsw(tji,params % se,params % eta)
            if (constants % fi(igrid,jsph).gt.one) then
                oji = xji/constants % fi(igrid,jsph)
            else
                oji = xji
            end if
            f1 = oji
            vb = vb + f1*alpha*Xadj_e(igrid,jsph)
            if (tji .gt. tlow) then
                call fmm_l2p_bessel_work(vtji, params % lmax, &
                    & constants % vscales, constants % SI_ri(:, isph), one, &
                    & Xe(:, isph), zero, beta_ji, work_complex, work)
            if (constants % fi(igrid,jsph) .gt. one) then
                dj  = one/constants % fi(igrid,jsph)
                fac = dj*xji
                proc = .false.
                b    = zero
                do jk = constants % inl(jsph), constants % inl(jsph+1) - 1
                    ksph = constants % nl(jk)
                    vjk  = params % csph(:,jsph) + &
                        & params % rsph(jsph)*constants % cgrid(:,igrid) - &
                        & params % csph(:,ksph)
                    rjkn = dnrm2(3, vjk, 1)
                    tjk  = rjkn/params % rsph(ksph)
                    if (ksph.ne.isph) then
                        if (tjk .le. thigh) then
                            proc = .true.
                            vtjk = vjk*params % kappa
                            call fmm_l2p_bessel_work(vtjk, params % lmax, &
                                & constants % vscales, &
                                & constants % SI_ri(:, ksph), one, Xe(:, ksph), &
                                & zero, beta_jk, work_complex, work)
                            xjk = fsw(tjk, params % se, params % eta)
                            b = b + beta_jk*xjk
                        end if
                    end if
                end do
                if (proc) then
                    g1 = dj*dj*dfsw(tji,params % se,params % eta) &
                        & /params % rsph(isph)
                    g2 = g1*Xadj_e(igrid,jsph)*b
                    vc = vc + g2*sji
                end if
              else
                  dj  = one
                  fac = zero
              end if
              f2 = (one-fac)*dj*dfsw(tji,params % se,params % eta) &
                  & /params % rsph(isph)
              vb = vb + f2*Xadj_e(igrid,jsph)*beta_ji*sji
            end if
        end do
        force = force + constants % wgrid(igrid)*(vb - vc)
    end do
end subroutine contract_gradi_Bji


!> Subroutine to calculate the third derivative term in C1_C2 matrix,
!! namely the derivative of PU_i
!!
!! @param[in]  params       : Input parameter file
!! @param[in]  constants    : Input constants file
!! @param[in]  workspace    : Input workspace
!! @param[in]  Xadj_r_sgrid : Adjoint Laplace solution evaluated at grid point
!! @param[in]  Xadj_e_sgrid : Adjoint HSP solution evaluated at grid point
!! @param[in]  diff_re      : l'/r_j[Xr]_jl'm' -(i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
!! @param[out] force        : Force
!! @param[inout] ddx_error: ddX error
subroutine contract_grad_C_worker1(params, constants, workspace, Xadj_r_sgrid, &
    & Xadj_e_sgrid, diff_re, force, ddx_error)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: &
        & diff_re
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: &
        & Xadj_r_sgrid, Xadj_e_sgrid
    real(dp), dimension(3, params % nsph), intent(inout) :: force
    type(ddx_error_type), intent(inout) :: ddx_error
    ! Local variable
    ! igrid0: Index for grid point n0
    integer :: isph, jsph, igrid, l0, m0, ind0, igrid0, icav, &
        & indl, inode, istat
    ! term  : SK_rijn/SK_rj
    ! termi : DI_ri/SI_ri
    ! termk : DK_ri/SK_ri
    ! sum_int : Intermediate sum
    ! sum_r   : Intermediate sum for Laplace
    ! sum_e   : Intermediate sum for HSP
    real(dp) :: sum_int, sum_r, sum_e
    real(dp) :: rijn
    real(dp)  :: vij(3), sij(3), vtij(3)

    ! local allocatable
    ! phi_n_r : Phi corresponding to Laplace problem
    ! phi_n_e : Phi corresponding to HSP problem
    ! coefY_d : sum_{l0m0} C_ik*term*Y_l0m0^j(x_in)*Y_l0m0(s_n)
    ! diff_re_sgrid : diff_re evaluated at grid point
    real(dp), allocatable :: phi_n_r(:,:), phi_n_e(:,:), coefY_d(:,:,:), &
        & diff_re_sgrid(:,:)

    ! basloc : Y_lm(s_n)
    ! vplm   : Argument to call ylmbas
    real(dp),  dimension(constants % nbasis):: basloc, vplm
    ! dbasloc : Derivative of Y_lm(s_n)
    real(dp),  dimension(3, constants % nbasis):: dbasloc
    ! vcos   : Argument to call ylmbas
    ! vsin   : Argument to call ylmbas
    real(dp),  dimension(params % lmax+1):: vcos, vsin
    ! SK_rijn : Besssel function of first kind for rijn
    ! DK_rijn : Derivative of Besssel function of first kind for rijn
    real(dp), dimension(0:params % lmax) :: SK_rijn, DK_rijn
    real(dp) :: coef(constants % nbasis0), work(constants % lmax0+1)
    complex(dp) :: work_complex(constants % lmax0+1)

    allocate(phi_n_r(params % ngrid, params % nsph), &
        & phi_n_e(params % ngrid, params % nsph), &
        & diff_re_sgrid(params % ngrid, params % nsph), stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, "allocation ddx_error in ddx_contract_grad_C_worker1")
        return
    end if
!
    sum_int = zero
    sum_r = zero
    sum_e = zero
    phi_n_r = zero
    phi_n_e = zero
    diff_re_sgrid = zero
    basloc = zero
    vplm = zero
    dbasloc = zero
    vcos = zero
    vsin = zero
    SK_rijn = zero
    DK_rijn = zero

    if (params % fmm .eq. 0) then
        allocate(coefY_d(constants % ncav, params % ngrid, params % nsph), &
            & stat=istat)
        if (istat.ne.0) then
            call update_error(ddx_error, "allocation ddx_error in fmm ddx_contract_grad_C_worker1")
            return
        end if
        coefY_d = zero
        ! Compute  summation over l0, m0
        ! Loop over the sphers j
        do jsph = 1, params % nsph
          ! Loop over the grid points n0
          do igrid0 = 1, params % ngrid
            icav = zero
            ! Loop over spheres i
            do isph = 1, params % nsph
              ! Loop over grid points n
              do igrid = 1, params % ngrid
                ! Check for U_i^{eta}(x_in)
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  vij  = params % csph(:,isph) + &
                         & params % rsph(isph)*constants % cgrid(:,igrid) - &
                         & params % csph(:,jsph)
                  rijn = sqrt(dot_product(vij,vij))
                  sij = vij/rijn

                  do l0 = 0, constants % lmax0
                    do m0 = -l0,l0
                      ind0 = l0**2 + l0 + m0 + 1
                      coef(ind0) = constants % vgrid(ind0, igrid0) * &
                          & constants % C_ik(l0, jsph)
                    end do
                  end do
                  vtij = vij*params % kappa
                  call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                      & constants % vscales, constants % SK_ri(:, jsph), one, &
                      & coef, zero, coefY_d(icav, igrid0, jsph), work_complex, work)
                end if
              end do
            end do
          end do
        end do

        ! Compute phi_in
        ! Loop over spheres j
        do jsph = 1, params % nsph
          ! Loop over grid points n0
          do igrid0 = 1, params % ngrid
            icav = zero
            sum_r = zero
            sum_e = zero
            ! Loop over sphers i
            do isph = 1, params % nsph
              ! Loop over grid points n
              do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  sum_r = sum_r + coefY_d(icav, igrid0, jsph)*Xadj_r_sgrid(igrid, isph) &
                          & * constants % wgrid(igrid)*constants % ui(igrid, isph)
                  sum_e = sum_e + coefY_d(icav, igrid0, jsph)*Xadj_e_sgrid(igrid, isph) &
                          & * constants % wgrid(igrid)*constants % ui(igrid, isph)
                end if
              end do
            end do
            phi_n_r(igrid0, jsph) = sum_r
            phi_n_e(igrid0, jsph) = sum_e
          end do
        end do
    else
        ! Compute phi_n_r at first
        ! Adjoint integration from spherical harmonics to grid points is not needed
        ! here as ygrid already contains grid values, we just need to scale it by
        ! weights of grid points
        do isph = 1, params % nsph
            workspace % tmp_grid(:, isph) = Xadj_r_sgrid(:, isph) * &
                & constants % wgrid(:) * constants % ui(:, isph)
        end do
        ! Adjoint FMM
        call tree_m2p_bessel_adj(params, constants, constants % lmax0, one, &
            & workspace % tmp_grid, zero, params % lmax, workspace % tmp_sph)
        call tree_l2p_bessel_adj(params, constants, one, workspace % tmp_grid, &
            & zero, workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l, workspace % tmp_node_m)
        call tree_m2m_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_m)
        ! Properly load adjoint multipole harmonics into tmp_sph
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph(1:constants % nbasis0, isph) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph) + &
                    & workspace % tmp_node_m(1:constants % nbasis0, inode)
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph(1:indl, isph) = &
                    & workspace % tmp_sph(1:indl, isph) + &
                    & workspace % tmp_node_m(:, inode)
            end do
        end if
        ! Scale by C_ik
        do isph = 1, params % nsph
            do l0 = 0, constants % lmax0
                ind0 = l0*l0 + l0 + 1
                workspace % tmp_sph(ind0-l0:ind0+l0, isph) = &
                    & workspace % tmp_sph(ind0-l0:ind0+l0, isph) * &
                    & constants % C_ik(l0, isph)
            end do
        end do
        ! Multiply by vgrid
        call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis0, &
            & one, constants % vgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, phi_n_r, params % ngrid)
        ! Compute phi_n_e now
        ! Adjoint integration from spherical harmonics to grid points is not needed
        ! here as ygrid already contains grid values, we just need to scale it by
        ! weights of grid points
        do isph = 1, params % nsph
            workspace % tmp_grid(:, isph) = Xadj_e_sgrid(:, isph) * &
                & constants % wgrid(:) * constants % ui(:, isph)
        end do
        ! Adjoint FMM
        call tree_m2p_bessel_adj(params, constants, constants % lmax0, one, &
            & workspace % tmp_grid, zero, params % lmax, workspace % tmp_sph)
        call tree_l2p_bessel_adj(params, constants, one, workspace % tmp_grid, &
            & zero, workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l, workspace % tmp_node_m)
        call tree_m2m_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_m)
        ! Properly load adjoint multipole harmonics into tmp_sph
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph(1:constants % nbasis0, isph) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph) + &
                    & workspace % tmp_node_m(1:constants % nbasis0, inode)
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph(1:indl, isph) = &
                    & workspace % tmp_sph(1:indl, isph) + &
                    & workspace % tmp_node_m(:, inode)
            end do
        end if
        ! Scale by C_ik
        do isph = 1, params % nsph
            do l0 = 0, constants % lmax0
                ind0 = l0*l0 + l0 + 1
                workspace % tmp_sph(ind0-l0:ind0+l0, isph) = &
                    & workspace % tmp_sph(ind0-l0:ind0+l0, isph) * &
                    & constants % C_ik(l0, isph)
            end do
        end do
        ! Multiply by vgrid
        call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis0, &
            & one, constants % vgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, phi_n_e, params % ngrid)
    end if
    call dgemm('T', 'N', params % ngrid, params % nsph, &
              & constants % nbasis, one, constants % vgrid, constants % vgrid_nbasis, &
              & diff_re , constants % nbasis, zero, diff_re_sgrid, &
              & params % ngrid)
    do isph = 1, params % nsph
        call contract_grad_U(params, constants, isph, diff_re_sgrid, phi_n_r, force(:, isph))
        call contract_grad_U(params, constants, isph, diff_re_sgrid, phi_n_e, force(:, isph))
    end do

    deallocate(phi_n_r, phi_n_e, diff_re_sgrid, stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, "deallocation ddx_error in ddx_contract_grad_C_worker1")
        return
    end if
    if (allocated(coefY_d)) then
        deallocate(coefY_d, stat=istat)
        if (istat.ne.0) then
            call update_error(ddx_error, "deallocation ddx_error in ddx_contract_grad_C_worker1")
            return
        end if
    end if
end subroutine contract_grad_C_worker1


!> Subroutine to compute the derivative of U_i(x_in) and \bf{k}^j_l0(x_in)Y^j_l0m0(x_in)
!!
!! contract_grad_C_worker2 : Force Derivative of U_i^e(x_in), k_l0, and Y_l0m0
!! @param[in]  params      : Input parameter file
!! @param[in]  constants   : Input constants file
!! @param[in]  workspace   : Input workspace
!! @param[in] ksph         : Derivative with respect to x_k
!! @param[in] Xr           : Solution of the Laplace problem
!! @param[in] Xe           : Solution of the HSP problem
!! @param[in] Xadj_r_sgrid : Solution of the Adjoint Laplace problem evaluated at the
!!                           grid
!! @param[in] Xadj_e_sgrid : Solution of the Adjoint HSP problem evaluated at the grid
!! @param[in] Xadj_r       : Adjoint solution of the Laplace problem
!! @param[in] Xadj_e       : Adjoint solution of the HSP problem
!! @param[inout] force     : Force
!! @param[out] diff_re     : epsilon_1/epsilon_2 * l'/r_j[Xr]_jl'm'
!!                         - (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
!! @param[inout] ddx_error: ddX error
subroutine contract_grad_C_worker2(params, constants, workspace, Xr, Xe, &
        & Xadj_r_sgrid, Xadj_e_sgrid, Xadj_r, Xadj_e, force, diff_re, ddx_error)
    !! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: Xr, Xe
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: &
        & Xadj_r_sgrid, Xadj_e_sgrid
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: &
        & Xadj_r, Xadj_e
    real(dp), dimension(3, params % nsph), intent(inout) :: force
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: &
        & diff_re
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), external :: dnrm2
    ! Local variable
    integer :: isph, jsph, igrid, l, m, ind, l0, ind0, icav, indl, inode, &
        & ksph, knode, jnode, knear, jsph_node, istat
    ! val_dim3 : Intermediate value array of dimension 3
    real(dp), dimension(3) :: vij, vtij
    ! val   : Intermediate variable to compute diff_ep
    real(dp) :: val
    ! large local are allocatable
    ! phi_in : sum_{j=1}^N diff0_j * coefY_j
    ! diff_ep_dim3 : 3 dimensional couterpart of diff_ep
    ! sum_dim3 : Storage of sum
    ! diff0       : dot_product([PU_j]_l0m0^l'm', l'/r_j[Xr]_jl'm' -
    !                        (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm')
    real(dp), allocatable :: diff_ep_dim3(:,:), phi_in(:,:), sum_dim3(:,:,:), &
        & diff0(:,:), diff1(:,:), diff1_grad(:,:,:), l2l_grad(:,:,:)
    real(dp), dimension(0:params % lmax) :: SK_rijn, DK_rijn
    complex(dp) :: work_complex(constants % lmax0+2)
    real(dp) :: work(constants % lmax0+2)

    allocate(diff_ep_dim3(3, constants % ncav), &
        & phi_in(params % ngrid, params % nsph), &
        & sum_dim3(3, constants % nbasis, params % nsph), &
        & diff0(constants % nbasis0, params % nsph), &
        & diff1(constants % nbasis0, params % nsph), &
        & diff1_grad((constants % lmax0+2)**2, 3, params % nsph), &
        & l2l_grad((params % pl+2)**2, 3, params % nsph), stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, "allocation ddx_error in ddx_contract_grad_C_worker2")
        return
    end if

    ! Setting initial values to zero
    SK_rijn = zero
    DK_rijn = zero

    diff_re = zero
    ! Compute l'/r_j[Xr]_jl'm' -(i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm'
    !$omp parallel do default(none) shared(params,diff_re,constants,xr,xe) &
    !$omp private(jsph,l,m,ind) schedule(dynamic)
    do jsph = 1, params % nsph
      do l = 0, params % lmax
        do m = -l,l
          ind = l**2 + l + m + 1
          diff_re(ind,jsph) = (params % epsp/params % eps)*(l/params % rsph(jsph)) * &
                & Xr(ind,jsph) - constants % termimat(l,jsph)*Xe(ind,jsph)
        end do
      end do
    end do

    ! diff0 = Pchi * diff_re, linear scaling
    diff0 = zero
    !$omp parallel do default(none) shared(params,constants,diff_re,diff0, &
    !$omp diff1,diff1_grad) private(jsph,l0,ind0) schedule(dynamic)
    do jsph = 1, params % nsph
      do l0 = 0, constants % lmax0
        do ind0 = l0*l0+1, l0*l0+2*l0+1
          diff0(ind0, jsph) = dot_product(diff_re(:,jsph), &
              & constants % Pchi(:,ind0, jsph))
          diff1(ind0, jsph) = diff0(ind0, jsph) * constants % C_ik(l0, jsph)
        end do
      end do
      ! Prepare diff1_grad
      call fmm_m2m_bessel_grad(constants % lmax0, constants % SK_ri(:, jsph), &
          & constants % vscales, diff1(:, jsph), diff1_grad(:, :, jsph))
    end do

    if (params % fmm .eq. 0) then
        ! phi_in = diff0 * coefY
        ! Here, summation over j takes place
        phi_in = zero
        icav = 0
        do isph = 1, params % nsph
          do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .gt. zero) then
              ! Extrenal grid point
              icav = icav + 1
              val = zero
              do jsph = 1, params % nsph
                !do ind0 = 1, constants % nbasis0
                !!====== This place requirs coefY, that is not precomputed anymore
                !  val = val + diff0(ind0,jsph)*constants % coefY(icav,ind0,jsph)
                !end do
                vij = params % csph(:, isph) + &
                    & params % rsph(isph)*constants % cgrid(:, igrid) - &
                    & params % csph(:, jsph)
                vtij = vij * params % kappa
                call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                    & constants % vscales, constants % SK_ri(:, jsph), one, &
                    & diff1(:, jsph), one, val, work_complex, work)
              end do
              phi_in(igrid, isph) = val
            end if
          end do
        end do
    else
        ! phi_in
        ! Load input harmonics into tree data
        workspace % tmp_sph = zero
        workspace % tmp_sph(1:constants % nbasis0, :) = diff1(:, :)
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(1:constants % nbasis0, inode) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph)
                workspace % tmp_node_m(constants % nbasis0+1:, inode) = zero
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(:, inode) = workspace % tmp_sph(1:indl, isph)
            end do
        end if
        ! Do FMM operations
        call tree_m2m_bessel_rotation(params, constants, workspace % tmp_node_m)
        call tree_m2l_bessel_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation(params, constants, workspace % tmp_node_l)
        call tree_l2p_bessel(params, constants, one, workspace % tmp_node_l, zero, &
            & phi_in)
        call tree_m2p_bessel(params, constants, constants % lmax0, one, &
            & params % lmax, workspace % tmp_sph, one, &
            & phi_in)
        ! Make phi_in zero at internal grid points
        !$omp parallel do default(none) shared(params,constants,phi_in) &
        !$omp private(isph,igrid) schedule(dynamic)
        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, isph) .eq. zero) then
                    phi_in(igrid, isph) = zero
                end if
            end do
        end do
        ! Get gradients of the L2L
        !$omp parallel do default(none) shared(params,constants,workspace, &
        !$omp l2l_grad) private(isph,igrid,inode) schedule(dynamic)
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_sph_l(:, isph) = workspace % tmp_node_l(:, inode)
            call fmm_l2l_bessel_grad(params % pl, &
                & constants % SI_ri(:, isph), constants % vscales, &
                & workspace % tmp_node_l(:, inode), &
                & l2l_grad(:, :, isph))
        end do
        workspace % tmp_sph = Xadj_r + Xadj_e
        call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
            & one, constants % vwgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, &
            & workspace % tmp_grid, params % ngrid)
        workspace % tmp_grid = workspace % tmp_grid * constants % ui
        ! Adjoint FMM with output tmp_sph2(:, :) which stores coefficients of
        ! harmonics of degree up to lmax+1
        call tree_m2p_bessel_nodiag_adj(params, constants, constants % lmax0+1, one, &
            & workspace % tmp_grid, zero, params % lmax+1, workspace % tmp_sph2)
        call tree_l2p_bessel_adj(params, constants, one, workspace % tmp_grid, zero, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_node_m)
        call tree_m2m_bessel_rotation_adj(params, constants, workspace % tmp_node_m)
        ! Properly load adjoint multipole harmonics into tmp_sph2 that holds
        ! harmonics of a degree up to lmax+1
        if(constants % lmax0+1 .lt. params % pm) then
            !$omp parallel do default(none) shared(params,constants, &
            !$omp workspace) private(isph,inode) schedule(dynamic)
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph2(1:(constants % lmax0+2)**2, isph) = &
                    & workspace % tmp_sph2(1:(constants % lmax0+2)**2, isph) + &
                    & workspace % tmp_node_m(1:(constants % lmax0+2)**2, inode)
            end do
        else
            indl = (params % pm+1)**2
            !$omp parallel do default(none) shared(params,constants,indl, &
            !$omp workspace) private(isph,inode) schedule(dynamic)
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph2(1:indl, isph) = &
                    & workspace % tmp_sph2(1:indl, isph) + &
                    & workspace % tmp_node_m(:, inode)
            end do
        end if
    end if

    do ksph = 1, params % nsph
        ! Computation of derivative of U_i^e(x_in)
        call contract_grad_U(params, constants, ksph, Xadj_r_sgrid, phi_in, force(:, ksph))
        call contract_grad_U(params, constants, ksph, Xadj_e_sgrid, phi_in, force(:, ksph))

        ! Aleksandr: my loop for the diff_ep_dim3
        diff_ep_dim3 = zero
        ! At first isph=ksph, jsph!=ksph
        icav = constants % icav_ia(ksph) - 1
        if (params % fmm .eq. 0) then
            do igrid = 1, params % ngrid
              if (constants % ui(igrid, ksph) .eq. zero) cycle
              icav = icav + 1
              do jsph = 1, params % nsph
                if (jsph .eq. ksph) cycle
                  vij  = params % csph(:,ksph) + &
                      & params % rsph(ksph)*constants % cgrid(:,igrid) - &
                      & params % csph(:,jsph)
                  vtij = vij * params % kappa
                  !call fmm_m2p_bessel_grad(vtij, &
                  !    & params % rsph(jsph)*params % kappa, &
                  !    & constants % lmax0, &
                  !    & constants % vscales, params % kappa, diff1(:, jsph), one, &
                  !    & diff_ep_dim3(:, icav))
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, diff1_grad(:, 1, jsph), one, &
                      & diff_ep_dim3(1, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, diff1_grad(:, 2, jsph), one, &
                      & diff_ep_dim3(2, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, diff1_grad(:, 3, jsph), one, &
                      & diff_ep_dim3(3, icav), work_complex, work)
              end do
            end do
        else
            knode = constants % snode(ksph)
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, ksph) .eq. zero) cycle
                icav = icav + 1
                ! Far-field
                call dgemv('T', (params % pl+2)**2, 3, -params % kappa, &
                    & l2l_grad(1, 1, ksph), &
                    & (params % pl+2)**2, constants % vgrid(1, igrid), 1, &
                    & one, diff_ep_dim3(1, icav), 1)
                ! Near-field
                do knear = constants % snear(knode), constants % snear(knode+1)-1
                    jnode = constants % near(knear)
                    do jsph_node = constants % cluster(1, jnode), &
                        & constants % cluster(2, jnode)
                        jsph = constants % order(jsph_node)
                        if (jsph .eq. ksph) cycle
                        vij = params % csph(:, ksph) - params % csph(:, jsph) + &
                            & params % rsph(ksph)*constants % cgrid(:, igrid)
                        vtij = vij * params % kappa
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, diff1_grad(:, 1, jsph), one, &
                            & diff_ep_dim3(1, icav), work_complex, work)
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, diff1_grad(:, 2, jsph), one, &
                            & diff_ep_dim3(2, icav), work_complex, work)
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, diff1_grad(:, 3, jsph), one, &
                            & diff_ep_dim3(3, icav), work_complex, work)
                    end do
                end do
            end do
        end if
        sum_dim3 = zero
        icav = constants % icav_ia(ksph) - 1
        do igrid =1, params % ngrid
            if(constants % ui(igrid, ksph) .gt. zero) then
                icav = icav + 1
                do ind = 1, constants % nbasis
                    sum_dim3(:,ind,ksph) = sum_dim3(:,ind,ksph) &
                        & + diff_ep_dim3(:,icav)*constants % ui(igrid, ksph) &
                        & *constants % vwgrid(ind, igrid)
                end do
            end if
        end do
        do ind = 1, constants % nbasis
            force(:, ksph) = force(:, ksph) + &
                & sum_dim3(:, ind, ksph)*(Xadj_r(ind, ksph) + &
                & Xadj_e(ind, ksph))
        end do

        ! Now jsph=ksph and isph!=ksph
        if (params % fmm .eq. 0) then
            diff_ep_dim3 = zero
            sum_dim3 = zero
            do isph = 1, params % nsph
              if (isph .eq. ksph) cycle
              icav = constants % icav_ia(isph) - 1
              do igrid = 1, params % ngrid
                  if (constants % ui(igrid, isph) .eq. zero) cycle
                  icav = icav + 1
                  vij  = params % csph(:,isph) + &
                      & params % rsph(isph)*constants % cgrid(:,igrid) - &
                      & params % csph(:,ksph)
                  vtij = vij * params % kappa
                  !call fmm_m2p_bessel_grad(vij * params % kappa, &
                  !    & params % rsph(ksph)*params % kappa, &
                  !    & constants % lmax0, &
                  !    & constants % vscales, -params % kappa, diff1(:, ksph), one, &
                  !    & diff_ep_dim3(:, icav))
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, diff1_grad(:, 1, ksph), one, &
                      & diff_ep_dim3(1, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, diff1_grad(:, 2, ksph), one, &
                      & diff_ep_dim3(2, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, diff1_grad(:, 3, ksph), one, &
                      & diff_ep_dim3(3, icav), work_complex, work)
              end do
            end do
            icav = zero
            do isph = 1, params % nsph
              do igrid =1, params % ngrid
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  do ind = 1, constants % nbasis
                    sum_dim3(:,ind,isph) = sum_dim3(:,ind,isph) &
                        & + diff_ep_dim3(:,icav)*constants % ui(igrid, isph) &
                        & *constants % vwgrid(ind, igrid)
                  end do
                end if
              end do
            end do
            ! Computation of derivative of \bf(k)_j^l0(x_in)\times Y^j_l0m0(x_in)
            do isph = 1, params % nsph
              do ind = 1, constants % nbasis
                force(:, ksph) = force(:, ksph) &
                    & + sum_dim3(:, ind, isph)*(Xadj_r(ind, isph) &
                    & + Xadj_e(ind, isph))
              end do
            end do
        else
            call dgemv('T', (constants % lmax0+2)**2, 3, params % kappa, &
                & diff1_grad(1, 1, ksph), (constants % lmax0+2)**2, &
                & workspace % tmp_sph2(1, ksph), 1, one, force(1, ksph), 1)
        end if
    end do
    deallocate(diff_ep_dim3, phi_in, sum_dim3, diff1_grad, l2l_grad, &
        & diff0, diff1, stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, "deallocation ddx_error in ddx_contract_grad_C_worker2")
        return
    end if
end subroutine contract_grad_C_worker2

!> Subroutine to compute the derivative of U_i(x_in) and \bf{k}^j_l0(x_in)Y^j_l0m0(x_in)
!! in F0
!!
!! contract_grad_f_worker1 : Force Derivative of U_i^e(x_in), k_l0, and Y_l0m0 in F0
!! @param[in]  params      : Input parameter file
!! @param[in]  constants   : Input constants file
!! @param[in]  workspace   : Input workspace
!! @param[in] ksph         : Derivative with respect to x_k
!! @param[in] sol_sgrid    : Solution of the Adjoint problem evaluated at the grid
!! @param[in] sol_adj      : Adjoint solution
!! @param[in] gradpsi      : Gradient of Psi_0
!! @param[inout] force     : Force
!! @param[inout] ddx_error: ddX error
subroutine contract_grad_f_worker1(params, constants, workspace, sol_adj, sol_sgrid, &
    & gradpsi, force, ddx_error)
    ! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: &
        & sol_adj
    real(dp), dimension(params % ngrid, params % nsph), intent(in) :: sol_sgrid
    real(dp), dimension(3, constants % ncav), intent(in) :: gradpsi
    real(dp), dimension(3, params % nsph), intent(inout) :: force
    type(ddx_error_type), intent(inout) :: ddx_error

    ! local
    real(dp), external :: dnrm2
    integer :: isph, jsph, igrid, ind, l0, ind0, icav, ksph, &
        & knode, jnode, knear, jsph_node, indl, inode, istat
    ! val_dim3 : Intermediate value array of dimension 3
    real(dp), dimension(3) :: vij, vtij
    ! val     : Intermediate variable to compute diff_ep
    ! nderpsi : Derivative of psi on grid points
    real(dp) :: nderpsi, sum_int

    ! local allocatable
    ! phi_in : sum_{j=1}^N diff0_j * coefY_j
    ! diff_ep_dim3 : 3 dimensional couterpart of diff_ep
    ! sum_dim3 : Storage of sum
    ! Debug purpose
    ! These variables can be taken from the subroutine update_rhs
    ! diff0       : dot_product([PU_j]_l0m0^l'm', l'/r_j[Xr]_jl'm' -
    !                        (i'_l'(r_j)/i_l'(r_j))[Xe]_jl'm')
    ! sum_Sjin : \sum_j [S]_{jin} Eq.~(97) [QSM20.SISC]
    ! c0 : \sum_{n=1}^N_g w_n U_j^{x_nj}\partial_n psi_0(x_nj)Y_{l0m0}(s_n)
    real(dp), allocatable :: phi_in(:,:), diff_ep_dim3(:,:), &
        & sum_dim3(:,:,:), diff0(:,:), sum_sjin(:,:), &
        & c0_d(:,:), c0_d1(:,:), c0_d1_grad(:,:,:), l2l_grad(:,:,:)
    real(dp), dimension(0:params % lmax) :: SK_rijn, DK_rijn
    complex(dp) :: work_complex(constants % lmax0 + 2)
    real(dp) :: work(constants % lmax0 + 2)

    allocate(phi_in(params % ngrid, params % nsph), &
        & diff_ep_dim3(3, constants % ncav), &
        & sum_dim3(3, constants % nbasis, params % nsph), &
        & c0_d(constants % nbasis0, params % nsph), &
        & c0_d1(constants % nbasis0, params % nsph), &
        & c0_d1_grad((constants % lmax0+2)**2, 3, params % nsph), &
        & sum_Sjin(params % ngrid, params % nsph), &
        & l2l_grad((params % pl+2)**2, 3, params % nsph), &
        & diff0(constants % nbasis0, params % nsph), stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, "allocation ddx_error in ddx_grad_f_worker1")
        return
    end if


    ! Setting initial values to zero
    SK_rijn = zero
    DK_rijn = zero
    c0_d = zero
    c0_d1 = zero
    c0_d1_grad = zero

    icav = zero
    do isph = 1, params % nsph
      do igrid= 1, params % ngrid
        if ( constants % ui(igrid,isph) .gt. zero ) then
          icav = icav + 1
          nderpsi = dot_product( gradpsi(:,icav),constants % cgrid(:,igrid) )
          c0_d(:, isph) = c0_d(:,isph) + &
                       & constants % wgrid(igrid)* &
                       & constants % ui(igrid,isph)*&
                       & nderpsi* &
                       & constants % vgrid(1:constants % nbasis0,igrid)
          do l0 = 0, constants % lmax0
              ind0 = l0*l0 + l0 + 1
              c0_d1(ind0-l0:ind0+l0, isph) = c0_d(ind0-l0:ind0+l0, isph) * &
                  & constants % C_ik(l0, isph)
          end do
          ! Prepare c0_d1_grad
          call fmm_m2m_bessel_grad(constants % lmax0, constants % SK_ri(:, isph), &
              & constants % vscales, c0_d1(:, isph), c0_d1_grad(:, :, isph))
        end if
      end do
    end do

    if (params % fmm .eq. 0) then
        ! Compute [S]_{jin}
        icav = 0
        do isph = 1, params % nsph
          do igrid = 1, params % ngrid
            if (constants % ui(igrid,isph).gt.zero) then
              icav = icav + 1
              sum_int = zero
              ! Loop to compute Sijn
              do jsph = 1, params % nsph
                vij = params % csph(:, isph) + &
                    & params % rsph(isph)*constants % cgrid(:, igrid) - &
                    & params % csph(:, jsph)
                vtij = vij * params % kappa
                call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                    & constants % vscales, constants % SK_ri(:, jsph), one, &
                    & c0_d1(:, jsph), one, sum_int, work_complex, work)
              end do ! End of loop jsph
              sum_Sjin(igrid,isph) = -(params % epsp/params % eps)*sum_int
            end if
          end do ! End of loop igrid
        end do ! End of loop isph
    else
        ! Load input harmonics into tree data
        workspace % tmp_sph = zero
        workspace % tmp_sph(1:constants % nbasis0, :) = c0_d1(:, :)
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(1:constants % nbasis0, inode) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph)
                workspace % tmp_node_m(constants % nbasis0+1:, inode) = zero
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(:, inode) = workspace % tmp_sph(1:indl, isph)
            end do
        end if
        ! Do FMM operations
        call tree_m2m_bessel_rotation(params, constants, workspace % tmp_node_m)
        call tree_m2l_bessel_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation(params, constants, workspace % tmp_node_l)
        call tree_l2p_bessel(params, constants, -params % epsp/params % eps, workspace % tmp_node_l, zero, &
            & sum_Sjin)
        call tree_m2p_bessel(params, constants, constants % lmax0, -params % epsp/params % eps, &
            & params % lmax, workspace % tmp_sph, one, &
            & sum_Sjin)
        ! Make phi_in zero at internal grid points
        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, isph) .eq. zero) then
                    sum_Sjin(igrid, isph) = zero
                end if
            end do
        end do
        ! Get gradients of the L2L
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_sph_l(:, isph) = workspace % tmp_node_l(:, inode)
            call fmm_l2l_bessel_grad(params % pl, &
                & constants % SI_ri(:, isph), constants % vscales, &
                & workspace % tmp_node_l(:, inode), &
                & l2l_grad(:, :, isph))
        end do
        workspace % tmp_sph = sol_adj
        call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
            & one, constants % vwgrid, constants % vgrid_nbasis, &
            & workspace % tmp_sph, constants % nbasis, zero, &
            & workspace % tmp_grid, params % ngrid)
        workspace % tmp_grid = workspace % tmp_grid * constants % ui
        ! Adjoint FMM with output tmp_sph2(:, :) which stores coefficients of
        ! harmonics of degree up to lmax+1
        call tree_m2p_bessel_nodiag_adj(params, constants, constants % lmax0+1, one, &
            & workspace % tmp_grid, zero, params % lmax+1, workspace % tmp_sph2)
        call tree_l2p_bessel_adj(params, constants, one, workspace % tmp_grid, zero, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_node_m)
        call tree_m2m_bessel_rotation_adj(params, constants, workspace % tmp_node_m)
        ! Properly load adjoint multipole harmonics into tmp_sph2 that holds
        ! harmonics of a degree up to lmax+1
        if(constants % lmax0+1 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph2(1:(constants % lmax0+2)**2, isph) = &
                    & workspace % tmp_sph2(1:(constants % lmax0+2)**2, isph) + &
                    & workspace % tmp_node_m(1:(constants % lmax0+2)**2, inode)
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_sph2(1:indl, isph) = &
                    & workspace % tmp_sph2(1:indl, isph) + &
                    & workspace % tmp_node_m(:, inode)
            end do
        end if
    end if
    do ksph = 1, params % nsph
        ! Computation of derivative of U_i^e(x_in)
        call contract_grad_U(params, constants, ksph, sol_sgrid, sum_Sjin, force(:, ksph))
        ! Aleksandr: my loop for the diff_ep_dim3
        diff_ep_dim3 = zero
        ! At first isph=ksph, jsph!=ksph
        icav = constants % icav_ia(ksph) - 1
        if (params % fmm .eq. 0) then
            do igrid = 1, params % ngrid
              if (constants % ui(igrid, ksph) .eq. zero) cycle
              icav = icav + 1
              do jsph = 1, params % nsph
                if (jsph .eq. ksph) cycle
                  vij  = params % csph(:,ksph) + &
                      & params % rsph(ksph)*constants % cgrid(:,igrid) - &
                      & params % csph(:,jsph)
                  vtij = vij * params % kappa
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, c0_d1_grad(:, 1, jsph), one, &
                      & diff_ep_dim3(1, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, c0_d1_grad(:, 2, jsph), one, &
                      & diff_ep_dim3(2, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, jsph), &
                      & -params % kappa, c0_d1_grad(:, 3, jsph), one, &
                      & diff_ep_dim3(3, icav), work_complex, work)
              end do
            end do
        else
            knode = constants % snode(ksph)
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, ksph) .eq. zero) cycle
                icav = icav + 1
                ! Far-field
                call dgemv('T', (params % pl+2)**2, 3, -params % kappa, &
                    & l2l_grad(1, 1, ksph), &
                    & (params % pl+2)**2, constants % vgrid(1, igrid), 1, &
                    & one, diff_ep_dim3(1, icav), 1)
                ! Near-field
                do knear = constants % snear(knode), constants % snear(knode+1)-1
                    jnode = constants % near(knear)
                    do jsph_node = constants % cluster(1, jnode), &
                        & constants % cluster(2, jnode)
                        jsph = constants % order(jsph_node)
                        if (jsph .eq. ksph) cycle
                        vij = params % csph(:, ksph) - params % csph(:, jsph) + &
                            & params % rsph(ksph)*constants % cgrid(:, igrid)
                        vtij = vij * params % kappa
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, c0_d1_grad(:, 1, jsph), one, &
                            & diff_ep_dim3(1, icav), work_complex, work)
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, c0_d1_grad(:, 2, jsph), one, &
                            & diff_ep_dim3(2, icav), work_complex, work)
                        call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & -params % kappa, c0_d1_grad(:, 3, jsph), one, &
                            & diff_ep_dim3(3, icav), work_complex, work)
                    end do
                end do
            end do
        end if
        sum_dim3 = zero
        icav = constants % icav_ia(ksph) - 1
        do igrid =1, params % ngrid
            if(constants % ui(igrid, ksph) .gt. zero) then
                icav = icav + 1
                do ind = 1, constants % nbasis
                    sum_dim3(:,ind,ksph) = sum_dim3(:,ind,ksph) &
                        & - (params % epsp/params % eps) &
                        & *diff_ep_dim3(:,icav)*constants % ui(igrid, ksph) &
                        & *constants % vwgrid(ind, igrid)
                end do
            end if
        end do
        do ind = 1, constants % nbasis
            force(:, ksph) = force(:, ksph) + &
                & sum_dim3(:, ind, ksph)*sol_adj(ind, ksph)
        end do

        ! Now jsph=ksph and isph!=ksph
        if (params % fmm .eq. 0) then
            diff_ep_dim3 = zero
            sum_dim3 = zero
            do isph = 1, params % nsph
              if (isph .eq. ksph) cycle
              icav = constants % icav_ia(isph) - 1
              do igrid = 1, params % ngrid
                  if (constants % ui(igrid, isph) .eq. zero) cycle
                  icav = icav + 1
                  vij  = params % csph(:,isph) + &
                      & params % rsph(isph)*constants % cgrid(:,igrid) - &
                      & params % csph(:,ksph)
                  vtij = vij * params % kappa
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, c0_d1_grad(:, 1, ksph), one, &
                      & diff_ep_dim3(1, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, c0_d1_grad(:, 2, ksph), one, &
                      & diff_ep_dim3(2, icav), work_complex, work)
                  call fmm_m2p_bessel_work(vtij, constants % lmax0+1, &
                      & constants % vscales, constants % SK_ri(:, ksph), &
                      & params % kappa, c0_d1_grad(:, 3, ksph), one, &
                      & diff_ep_dim3(3, icav), work_complex, work)
              end do
            end do

            icav = zero
            do isph = 1, params % nsph
              do igrid =1, params % ngrid
                if(constants % ui(igrid, isph) .gt. zero) then
                  icav = icav + 1
                  do ind = 1, constants % nbasis
                    sum_dim3(:,ind,isph) = sum_dim3(:,ind,isph) &
                      & - (params % epsp/params % eps) &
                      & *diff_ep_dim3(:,icav)*constants % ui(igrid, isph) &
                      & *constants % vwgrid(ind, igrid)
                  end do
                end if
              end do
            end do

            ! Computation of derivative of \bf(k)_j^l0(x_in)\times Y^j_l0m0(x_in)
            do isph = 1, params % nsph
              do ind = 1, constants % nbasis
                force(:, ksph) = force(:, ksph) + sum_dim3(:, ind, isph)*sol_adj(ind, isph)
              end do
            end do
        else
            call dgemv('T', (constants % lmax0+2)**2, 3, -params % epsp/params % eps*params % kappa, &
                & c0_d1_grad(1, 1, ksph), (constants % lmax0+2)**2, &
                & workspace % tmp_sph2(1, ksph), 1, one, force(1, ksph), 1)
        end if
    end do

    deallocate(phi_in, diff_ep_dim3, sum_dim3, c0_d, c0_d1, &
          & c0_d1_grad, sum_Sjin, l2l_grad, diff0, stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, "deallocation ddx_error in ddx_grad_f_worker1")
        return
    end if

end subroutine contract_grad_f_worker1

!! @param[inout] ddx_error: ddX error
subroutine contract_grad_f_worker2(params, constants, &
        & gradpsi, normal_hessian_cav, force, state, ddx_error)
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: gradpsi(3, constants % ncav), &
        & normal_hessian_cav(3, constants % ncav)
    real(dp), intent(inout) :: force(3, params % nsph)
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error

    integer :: icav, isph, igrid, istat
    real(dp) :: nderpsi
    real(dp), allocatable :: gradpsi_grid(:, :)

    allocate(gradpsi_grid(params % ngrid, params % nsph), stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, "allocation ddx_error in ddx_grad_f_worker2")
        return
    end if

    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .gt. zero) then
                icav = icav + 1
                nderpsi = dot_product(gradpsi(:, icav), &
                    & constants % cgrid(:, igrid))
                gradpsi_grid(igrid, isph) = nderpsi
            end if
        end do
    end do

    icav = 0
    do isph = 1, params % nsph
        call contract_grad_U(params, constants, isph, gradpsi_grid, &
            & state % phi_n, force(:, isph))
        do igrid = 1, params % ngrid
          if(constants % ui(igrid, isph) .gt. zero) then
            icav = icav + 1
            force(:, isph) = force(:, isph) &
                & + constants % wgrid(igrid)*constants % ui(igrid, isph) &
                & *state % phi_n(igrid, isph)*normal_hessian_cav(:, icav)
          end if
        end do
    end do

    deallocate(gradpsi_grid, stat=istat)
    if (istat.ne.0) then
        call update_error(ddx_error, "deallocation ddx_error in ddx_grad_f_worker2")
        return
    end if

end subroutine contract_grad_f_worker2

!> Sphere contribution to the ddPCM matrix gradient using N^2 code
subroutine gradr_sph(params, constants, isph, vplm, vcos, vsin, basloc, &
        & dbsloc, g, ygrid, fx)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: isph
    real(dp), intent(in) :: g(constants % nbasis, params % nsph), &
        & ygrid(params % ngrid, params % nsph)
    real(dp), intent(inout) :: vplm(constants % nbasis), vcos(params % lmax+1), &
        & vsin(params % lmax+1), basloc(constants % nbasis), &
        & dbsloc(3, constants % nbasis), fx(3)
    ! various scratch arrays
    real(dp) vik(3), sik(3), vki(3), ski(3), vkj(3), skj(3), vji(3), &
        & sji(3), va(3), vb(3), a(3)
    ! jacobian matrix
    real(dp) sjac(3,3)
    ! indexes
    integer its, ik, ksph, l, m, ind, jsph, icomp, jcomp
    ! various scalar quantities
    real(dp) cx, cy, cz, vvki, tki, gg, fl, fac, vvkj, tkj
    real(dp) tt, fcl, fjj, gi, fii, vvji, tji, qji
    real(dp) b, vvik, tik, qik, tlow, thigh, duj
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    real(dp), external :: dnrm2

    tlow  = one - pt5*(one - params % se)*params % eta
    thigh = one + pt5*(one + params % se)*params % eta

    ! first set of contributions:
    ! diagonal block, kc and part of kb

    fx = zero
    do its = 1, params % ngrid
        ! sum over ksph in neighbors of isph
        do ik = constants % inl(isph), constants % inl(isph+1) - 1
            ksph = constants % nl(ik)
            ! build geometrical quantities
            cx = params % csph(1,ksph) + params % rsph(ksph)*constants % cgrid(1,its)
            cy = params % csph(2,ksph) + params % rsph(ksph)*constants % cgrid(2,its)
            cz = params % csph(3,ksph) + params % rsph(ksph)*constants % cgrid(3,its)
            vki(1) = cx - params % csph(1,isph)
            vki(2) = cy - params % csph(2,isph)
            vki(3) = cz - params % csph(3,isph)
            !vvki = sqrt(vki(1)*vki(1) + vki(2)*vki(2) + &
            !    & vki(3)*vki(3))
            vvki = dnrm2(3, vki, 1)
            tki  = vvki/params % rsph(isph)

            ! contributions involving grad i of uk come from the switching
            ! region.
            ! note: ui avoids contributions from points that are in the
            ! switching between isph and ksph but are buried in a third
            ! sphere.
            if ((tki.gt.tlow).and.(tki.lt.thigh) .and. &
                & constants % ui(its,ksph).gt.zero) then
                ! other geometrical quantities
                ski = vki/vvki

                ! diagonal block kk contribution, with k in n(i)
                gg = zero
                do l = 0, params % lmax
                    ind = l*l + l + 1
                    fl = dble(l)
                    fac = twopi/(two*fl + one)
                    do m = -l, l
                        !! DEBUG comment
                        gg = gg + fac*constants % vgrid(ind+m,its)*g(ind+m,ksph)
                    end do
                end do

                ! kc contribution
                do jsph = 1, params % nsph
                    if (jsph.ne.ksph .and. jsph.ne.isph) then 
                        vkj(1) = cx - params % csph(1,jsph)
                        vkj(2) = cy - params % csph(2,jsph)
                        vkj(3) = cz - params % csph(3,jsph)
                        vvkj = sqrt(vkj(1)*vkj(1) + vkj(2)*vkj(2) + &
                            & vkj(3)*vkj(3))
                        vvkj = dnrm2(3, vkj, 1)
                        tkj  = vvkj/params % rsph(jsph)
                        skj  = vkj/vvkj
                        call ylmbas(skj, rho, ctheta, stheta, cphi, sphi, &
                            & params % lmax, constants % vscales, basloc, &
                            & vplm, vcos, vsin)
                        tt = one/tkj
                        do l = 0, params % lmax
                            ind = l*l + l + 1
                            fcl = - fourpi*dble(l)/(two*dble(l)+one)*tt
                            do m = -l, l
                                !! DEBUG comment
                                gg = gg + fcl*g(ind+m,jsph)*basloc(ind+m)
                            end do
                            tt = tt/tkj
                        end do
                        !call fmm_m2p(vkj, params % rsph(jsph), &
                        !    & params % lmax, constants % vscales_rel, -one, &
                        !    & g(:, jsph), one, gg)
                    end if
                end do

                ! part of kb contribution
                call ylmbas(ski, rho, ctheta, stheta, cphi, sphi, &
                    & params % lmax, constants % vscales, basloc, &
                    & vplm, vcos, vsin)
                tt = one/tki
                do l = 0, params % lmax
                    ind = l*l + l + 1
                    fcl = - four*pi*dble(l)/(two*dble(l)+one)*tt
                    do m = -l, l
                        !! DEBUG comment
                        gg = gg + fcl*g(ind+m,isph)*basloc(ind+m)
                    end do
                    tt = tt/tki
                end do
                !call fmm_m2p(vki, params % rsph(isph), &
                !    & params % lmax, constants % vscales_rel, -one, &
                !    & g(:, isph), one, gg)

                ! common step, product with grad i uj
                duj = dfsw(tki,params % se, params % eta)/params % rsph(isph)
                fjj = duj*constants % wgrid(its)*gg*ygrid(its,ksph)
                fx(1) = fx(1) - fjj*ski(1)
                fx(2) = fx(2) - fjj*ski(2)
                fx(3) = fx(3) - fjj*ski(3)
            end if
        end do

        ! diagonal block ii contribution
        if (constants % ui(its,isph).gt.zero.and.constants % ui(its,isph).lt.one) then
            gi = zero
            do l = 0, params % lmax
                ind = l*l + l + 1
                fl = dble(l)
                fac = twopi/(two*fl + one)
                do m = -l, l 
                    !! DEBUG comment
                    gi = gi + fac*constants % vgrid(ind+m,its)*g(ind+m,isph)
                    !gi = gi + pt5*constants % vgrid2(ind+m,its)*g(ind+m,isph)
                end do
            end do
            !do l = 0, (params % lmax+1)**2
            !    gi = gi + constants % vgrid2(l, its)*g(l, isph)
            !end do
            !gi = pt5 * gi
            fii = constants % wgrid(its)*gi*ygrid(its,isph)
            fx(1) = fx(1) + fii*constants % zi(1,its,isph)
            fx(2) = fx(2) + fii*constants % zi(2,its,isph)
            fx(3) = fx(3) + fii*constants % zi(3,its,isph)
        end if
    end do

    ! second set of contributions:
    ! part of kb and ka
    do its = 1, params % ngrid

        ! run over all the spheres except isph 
        do jsph = 1, params % nsph
            if (constants % ui(its,jsph).gt.zero .and. jsph.ne.isph) then
                ! build geometrical quantities
                cx = params % csph(1,jsph) + params % rsph(jsph)*constants % cgrid(1,its)
                cy = params % csph(2,jsph) + params % rsph(jsph)*constants % cgrid(2,its)
                cz = params % csph(3,jsph) + params % rsph(jsph)*constants % cgrid(3,its)
                vji(1) = cx - params % csph(1,isph)
                vji(2) = cy - params % csph(2,isph)
                vji(3) = cz - params % csph(3,isph)
                !vvji = sqrt(vji(1)*vji(1) + vji(2)*vji(2) + &
                !    &  vji(3)*vji(3))
                vvji = dnrm2(3, vji, 1)
                tji = vvji/params % rsph(isph)
                qji = one/vvji
                sji = vji/vvji

                ! build the jacobian of sji
                sjac = zero
                sjac(1,1) = - one
                sjac(2,2) = - one
                sjac(3,3) = - one
                do icomp = 1, 3
                    do jcomp = 1, 3
                        sjac(icomp,jcomp) = qji*(sjac(icomp,jcomp) &
                            & + sji(icomp)*sji(jcomp))
                    end do
                end do

                ! assemble the local basis and its gradient
                !call dbasis(sji,basloc,dbsloc,vplm,vcos,vsin)
                call dbasis(params, constants, sji,basloc,dbsloc,vplm,vcos,vsin)

                ! assemble the contribution
                a = zero
                tt = one/(tji)
                do l = 0, params % lmax
                    ind = l*l + l + 1
                    fl = dble(l)
                    fcl = - tt*fourpi*fl/(two*fl + one)
                    do m = -l, l
                        fac = fcl*g(ind+m,isph)
                        b = (fl + one)*basloc(ind+m)/(params % rsph(isph)*tji)

                        ! apply the jacobian to grad y
                        va(1) = sjac(1,1)*dbsloc(1,ind+m) + &
                            & sjac(1,2)*dbsloc(2,ind+m) + sjac(1,3)*dbsloc(3,ind+m)
                        va(2) = sjac(2,1)*dbsloc(1,ind+m) + &
                            & sjac(2,2)*dbsloc(2,ind+m) + sjac(2,3)*dbsloc(3,ind+m)
                        va(3) = sjac(3,1)*dbsloc(1,ind+m) + &
                            & sjac(3,2)*dbsloc(2,ind+m) + sjac(3,3)*dbsloc(3,ind+m)
                        a(1) = a(1) + fac*(sji(1)*b + va(1))
                        a(2) = a(2) + fac*(sji(2)*b + va(2))
                        a(3) = a(3) + fac*(sji(3)*b + va(3))
                    end do
                    tt = tt/tji
                end do
                fac = constants % ui(its,jsph)*constants % wgrid(its)*ygrid(its,jsph)
                fx(1) = fx(1) - fac*a(1)
                fx(2) = fx(2) - fac*a(2)
                fx(3) = fx(3) - fac*a(3)
            end if
        end do
    end do

    ! ka contribution
    do its = 1, params % ngrid
        cx = params % csph(1,isph) + params % rsph(isph)*constants % cgrid(1,its)
        cy = params % csph(2,isph) + params % rsph(isph)*constants % cgrid(2,its)
        cz = params % csph(3,isph) + params % rsph(isph)*constants % cgrid(3,its)
        a = zero

        ! iterate on all the spheres except isph
        do ksph = 1, params % nsph
            if (constants % ui(its,isph).gt.zero .and. ksph.ne.isph) then
                ! geometrical stuff
                vik(1) = cx - params % csph(1,ksph)
                vik(2) = cy - params % csph(2,ksph)
                vik(3) = cz - params % csph(3,ksph)
                !vvik = sqrt(vik(1)*vik(1) + vik(2)*vik(2) + & 
                !    & vik(3)*vik(3))
                vvik = dnrm2(3, vik, 1)
                tik = vvik/params % rsph(ksph)
                qik = one/vvik
                sik = vik/vvik

                ! build the jacobian of sik
                sjac = zero
                sjac(1,1) = one
                sjac(2,2) = one
                sjac(3,3) = one
                do icomp = 1, 3
                    do jcomp = 1, 3
                    sjac(icomp,jcomp) = qik*(sjac(icomp,jcomp) &
                        & - sik(icomp)*sik(jcomp))
                    end do
                end do

                ! if we are in the switching region, recover grad_i u_i
                vb = zero
                if (constants % ui(its,isph).lt.one) then
                    vb(1) = constants % zi(1,its,isph)
                    vb(2) = constants % zi(2,its,isph)
                    vb(3) = constants % zi(3,its,isph)
                end if

                ! assemble the local basis and its gradient
                !call dbasis(sik,basloc,dbsloc,vplm,vcos,vsin)
                call dbasis(params, constants, sik,basloc,dbsloc,vplm,vcos,vsin)

                ! assemble the contribution
                tt = one/(tik)
                do l = 0, params % lmax
                    ind = l*l + l + 1
                    fl = dble(l)
                    fcl = - tt*fourpi*fl/(two*fl + one)
                        do m = -l, l
                        fac = fcl*g(ind+m,ksph)
                        fac = - fac*basloc(ind+m)
                        a(1) = a(1) + fac*vb(1)
                        a(2) = a(2) + fac*vb(2) 
                        a(3) = a(3) + fac*vb(3)

                        fac = constants % ui(its,isph)*fcl*g(ind+m,ksph)
                        b = - (fl + one)*basloc(ind+m)/(params % rsph(ksph)*tik)

                        ! apply the jacobian to grad y
                        va(1) = sjac(1,1)*dbsloc(1,ind+m) + &
                            & sjac(1,2)*dbsloc(2,ind+m) + sjac(1,3)*dbsloc(3,ind+m)
                        va(2) = sjac(2,1)*dbsloc(1,ind+m) + &
                            & sjac(2,2)*dbsloc(2,ind+m) + sjac(2,3)*dbsloc(3,ind+m)
                        va(3) = sjac(3,1)*dbsloc(1,ind+m) + &
                            & sjac(3,2)*dbsloc(2,ind+m) + sjac(3,3)*dbsloc(3,ind+m)
                        a(1) = a(1) + fac*(sik(1)*b + va(1))
                        a(2) = a(2) + fac*(sik(2)*b + va(2))
                        a(3) = a(3) + fac*(sik(3)*b + va(3))
                    end do
                    tt = tt/tik
                end do
            end if
        end do
        fac = constants % wgrid(its)*ygrid(its,isph)
        fx(1) = fx(1) - fac*a(1)
        fx(2) = fx(2) - fac*a(2)
        fx(3) = fx(3) - fac*a(3)
    end do
end subroutine gradr_sph

!> Compute the ddPCM matrix gradient using FMMS (2 matvecs)
subroutine gradr_fmm(params, constants, workspace, g, ygrid, fx)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: g(constants % nbasis, params % nsph), &
        & ygrid(params % ngrid, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: fx(3, params % nsph)
    ! Local variables
    integer :: indl, indl1, l, isph, igrid, ik, ksph, &
        & jsph, jsph_node
    integer :: inear, inode, jnode
    real(dp) :: gg, c(3), vki(3), vvki, tki, gg3(3), tmp_gg, tmp_c(3)
    real(dp) :: tlow, thigh
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    real(dp), external :: ddot, dnrm2
    real(dp) :: work(params % lmax+2)
    !real(dp) :: l2g(params % ngrid, params % nsph)
    zx_coord_transform = 0
    zx_coord_transform(3, 2) = 1
    zx_coord_transform(2, 3) = 1
    zx_coord_transform(1, 1) = 1
    zy_coord_transform = 0
    zy_coord_transform(1, 2) = 1
    zy_coord_transform(2, 1) = 1
    zy_coord_transform(3, 3) = 1
    tlow  = one - pt5*(one - params % se)*params % eta
    thigh = one + pt5*(one + params % se)*params % eta
    fx = zero
    !! Scale input harmonics at first
    workspace % tmp_sph(1, :) = zero
    indl = 2
    do l = 1, params % lmax
        indl1 = (l+1)**2
        workspace % tmp_sph(indl:indl1, :) = l * g(indl:indl1, :)
        indl = indl1 + 1
    end do
    !! Compute gradient of M2M of tmp_sph harmonics at the origin and store it
    !! in tmp_sph_grad. tmp_sph_grad(:, 1, :), tmp_sph_grad(:, 2, :) and
    !! tmp_sph_grad(:, 3, :) correspond to the OX, OY and OZ axes. Variable
    !! tmp_sph2 is a temporary workspace here.
    call tree_grad_m2m(params, constants, workspace % tmp_sph, &
        & workspace % tmp_sph_grad, workspace % tmp_sph2)
    !! Adjoint full FMM matvec to get output multipole expansions from input
    !! external grid points. It is used to compute R_i^B fast as a contraction
    !! of a gradient stored in tmp_sph_grad and a result of adjoint matvec.
    ! Adjoint integration from spherical harmonics to grid points is not needed
    ! here as ygrid already contains grid values, we just need to scale it by
    ! weights of grid points
    do isph = 1, params % nsph
        workspace % tmp_grid(:, isph) = ygrid(:, isph) * &
            & constants % wgrid(:) * constants % ui(:, isph)
    end do
    ! Adjoint FMM with output tmp_sph2(:, :) which stores coefficients of
    ! harmonics of degree up to lmax+1
    call tree_m2p_adj(params, constants, params % lmax+1, one, &
        & workspace % tmp_grid, zero, workspace % tmp_sph2)
    call tree_l2p_adj(params, constants, one, workspace % tmp_grid, zero, &
        & workspace % tmp_node_l, workspace % tmp_sph_l)
    call tree_l2l_rotation_adj(params, constants, workspace % tmp_node_l)
    call tree_m2l_rotation_adj(params, constants, workspace % tmp_node_l, &
        & workspace % tmp_node_m)
    call tree_m2m_rotation_adj(params, constants, workspace % tmp_node_m)
    ! Properly load adjoint multipole harmonics into tmp_sph2 that holds
    ! harmonics of a degree up to lmax+1
    if(params % lmax+1 .lt. params % pm) then
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_sph2(:, isph) = workspace % tmp_sph2(:, isph) + &
                & workspace % tmp_node_m(1:constants % grad_nbasis, inode)
        end do
    else
        indl = (params % pm+1)**2
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_sph2(1:indl, isph) = &
                & workspace % tmp_sph2(1:indl, isph) + &
                & workspace % tmp_node_m(:, inode)
        end do
    end if
    ! Compute second term of R_i^B as a contraction
    do isph = 1, params % nsph
        call dgemv('T', constants % grad_nbasis, 3, one, &
            & workspace % tmp_sph_grad(1, 1, isph), constants % grad_nbasis, &
            & workspace % tmp_sph2(1, isph), 1, zero, fx(1, isph), 1)
    end do
    !! Direct far-field FMM matvec to get output local expansions from input
    !! multipole expansions. It will be used in R_i^A.
    !! As of now I compute potential at all external grid points, improved
    !! version shall only compute it at external points in a switch region
    ! Load input harmonics into tree data
    if(params % lmax .lt. params % pm) then
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_m(:constants % nbasis, inode) = &
                & workspace % tmp_sph(:, isph)
            workspace % tmp_node_m(constants % nbasis+1:, inode) = zero
        end do
    else
        indl = (params % pm+1)**2
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_m(:, inode) = workspace % tmp_sph(1:indl, isph)
        end do
    end if
    ! Perform direct FMM matvec to all external grid points
    call tree_m2m_rotation(params, constants, workspace % tmp_node_m)
    call tree_m2l_rotation(params, constants, workspace % tmp_node_m, &
        & workspace % tmp_node_l)
    call tree_l2l_rotation(params, constants, workspace % tmp_node_l)
    call tree_l2p(params, constants, one, workspace % tmp_node_l, zero, &
        & workspace % tmp_grid, workspace % tmp_sph_l)
    call tree_m2p(params, constants, params % lmax, one, &
        & workspace % tmp_sph, one, workspace % tmp_grid)
    !! Compute gradients of L2L if pl > 0
    if (params % pl .gt. 0) then
        call tree_grad_l2l(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_sph_l_grad, workspace % tmp_sph_l)
    end if
    !! Diagonal update of computed grid values, that is needed for R^C, a part
    !! of R^A and a part of R^B
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & pt5, constants % vgrid2, constants % vgrid_nbasis, g, &
        & constants % nbasis, -one, workspace % tmp_grid, params % ngrid)
    !! Scale temporary grid points by corresponding Lebedev weights and ygrid
    do igrid = 1, params % ngrid
        do isph = 1, params % nsph
            workspace % tmp_grid(igrid, isph) = &
                & workspace % tmp_grid(igrid, isph) *  constants % wgrid(igrid) * &
                & ygrid(igrid, isph)
        end do
    end do
    !! Compute all terms of grad_i(R). The second term of R_i^B is already
    !! taken into account and the first term is computed together with R_i^C.
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            ! Loop over all neighbouring spheres
            do ik = constants % inl(isph), constants % inl(isph+1) - 1
                ksph = constants % nl(ik)
                ! Only consider external grid points
                if(constants % ui(igrid, ksph) .eq. zero) cycle
                ! build geometrical quantities
                c = params % csph(:, ksph) + &
                    & params % rsph(ksph)*constants % cgrid(:, igrid)
                vki = c - params % csph(:, isph)
                !vvki = sqrt(vki(1)*vki(1) + vki(2)*vki(2) + &
                !    & vki(3)*vki(3))
                vvki = dnrm2(3, vki, 1)
                tki = vvki / params % rsph(isph)
                ! Only consider such points where grad U is non-zero
                if((tki.le.tlow) .or. (tki.ge.thigh)) cycle
                ! This is entire R^C and the first R^B component (grad_i of U
                ! of a sum of R_kj for index inequality j!=k)
                ! Indexes k and j are flipped compared to the paper
                gg = workspace % tmp_grid(igrid, ksph)
                ! Compute grad_i component of forces using precomputed
                ! potential gg
                !fx(:, isph) = fx(:, isph) - &
                !    & dfsw(tki, params % se, params % eta)/ &
                !    & params % rsph(isph)*constants % wgrid(igrid)*gg* &
                !    & ygrid(igrid, ksph)*(vki/vvki)
                fx(:, isph) = fx(:, isph) - &
                    & dfsw(tki, params % se, params % eta)/ &
                    & params % rsph(isph)*gg*(vki/vvki)
            end do
            ! contribution from the sphere itself
            if((constants % ui(igrid,isph).gt.zero) .and. &
                & (constants % ui(igrid,isph).lt.one)) then
                ! R^A component (grad_i of U of a sum of R_ij for index
                ! inequality j!=i)
                ! Indexes k and j are flipped compared to the paper
                gg = workspace % tmp_grid(igrid, isph)
                ! Compute grad_i component of forces using precomputed
                ! potential gg
                !fx(:, isph) = fx(:, isph) + constants % wgrid(igrid)*gg* &
                !    & ygrid(igrid, isph)*constants % zi(:, igrid, isph)
                fx(:, isph) = fx(:, isph) + gg*constants % zi(:, igrid, isph)
            end if
            if (constants % ui(igrid, isph) .gt. zero) then
                ! Another R^A component (grad_i of potential of a sum of R_ij
                ! for index inequality j!=i)
                ! Indexes k and j are flipped compared to the paper
                ! In case pl=0 MKL does not make gg3 zero reusing old value of
                ! gg3, so we have to clear it manually
                gg3 = zero
                call dgemv('T', params % pl**2, 3, one, &
                    & workspace % tmp_sph_l_grad(1, 1, isph), &
                    & (params % pl+1)**2, constants % vgrid2(1, igrid), 1, &
                    & zero, gg3, 1)
                ! Gradient of the near-field potential is a gradient of
                ! multipole expansion
                inode = constants % snode(isph)
                do inear = constants % snear(inode), constants % snear(inode+1)-1
                    jnode = constants % near(inear)
                    do jsph_node = constants % cluster(1, jnode), &
                        & constants % cluster(2, jnode)
                        jsph = constants % order(jsph_node)
                        if (isph .eq. jsph) cycle
                        c = params % csph(:, isph) + &
                            & params % rsph(isph)*constants % cgrid(:, igrid)
                        tmp_c = c - params % csph(:, jsph)
                        call fmm_m2p_work(tmp_c, &
                            & params % rsph(jsph), params % lmax+1, &
                            & constants % vscales_rel, one, &
                            & workspace % tmp_sph_grad(:, 1, jsph), zero, &
                            & tmp_gg, work)
                        gg3(1) = gg3(1) + tmp_gg
                        call fmm_m2p_work(tmp_c, &
                            & params % rsph(jsph), params % lmax+1, &
                            & constants % vscales_rel, one, &
                            & workspace % tmp_sph_grad(:, 2, jsph), zero, &
                            & tmp_gg, work)
                        gg3(2) = gg3(2) + tmp_gg
                        call fmm_m2p_work(tmp_c, &
                            & params % rsph(jsph), params % lmax+1, &
                            & constants % vscales_rel, one, &
                            & workspace % tmp_sph_grad(:, 3, jsph), zero, &
                            & tmp_gg, work)
                        gg3(3) = gg3(3) + tmp_gg
                    end do
                end do
                ! Accumulate all computed forces
                fx(:, isph) = fx(:, isph) - constants % wgrid(igrid)*gg3* &
                    & ygrid(igrid, isph)*constants % ui(igrid, isph)
            end if
        end do
    end do
end subroutine gradr_fmm

!> Gradient of the ddPCM matrix
subroutine gradr(params, constants, workspace, g, ygrid, fx)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: g(constants % nbasis, params % nsph), &
        & ygrid(params % ngrid, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: fx(3, params % nsph)
    ! Check which gradr to execute
    if (params % fmm .eq. 1) then
        call gradr_fmm(params, constants, workspace, g, ygrid, fx)
    else
        call gradr_dense(params, constants, workspace, g, ygrid, fx)
    end if
end subroutine gradr

!> Gradient of the ddPCM matrix using N^2 code
subroutine gradr_dense(params, constants, workspace, g, ygrid, fx)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: g(constants % nbasis, params % nsph), &
        & ygrid(params % ngrid, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: fx(3, params % nsph)
    ! Local variables
    integer :: isph
    ! Simply cycle over all spheres
    do isph = 1, params % nsph
        call gradr_sph(params, constants, isph, workspace % tmp_vplm, &
            & workspace % tmp_vcos, workspace % tmp_vsin, &
            & workspace % tmp_vylm, workspace % tmp_vdylm, &
            & g, ygrid, fx(:, isph))
    end do
end subroutine gradr_dense

!> Force term: interaction of the external electric field with the zeta
!! intermediate. This routine is called by the gradients of ddCOSMO
!! ddPCM and ddLPB.
subroutine zeta_grad(params, constants, state, e_cav, forces)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(inout) :: forces(3, params % nsph)
    real(dp), intent(in) :: e_cav(3, constants % ncav)
    ! local variables
    integer :: icav, isph, igrid

    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if (constants % ui(igrid, isph) .eq. zero) cycle
            icav = icav + 1
            forces(:, isph) = forces(:, isph) + pt5 &
                & *state % zeta(icav)*e_cav(:, icav)
        end do
    end do
end subroutine zeta_grad

end module ddx_gradients
