!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_harmonics.f90
!! This file contains routines to evaluate spherical harmonics, the multipole
!! expansioni M2P, the local expansion L2P and sphere-to-sphere M2M, M2L and
!! L2L operations of the FMM along with routines to precompute special
!! harmonics-related constants.
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

!> Harmonics-related core routines
module ddx_harmonics

! Get all compile-time definitions
use ddx_definitions

! Disable implicit types
implicit none

contains

!> Compute scaling factors of real normalized spherical harmonics
!!
!! Output values of scaling factors of \f$ Y_\ell^m \f$ harmonics are filled
!! only for non-negative \f$ m \f$ since scaling factor of \f$ Y_\ell^{-m} \f$
!! is the same as scaling factor of \f$ Y_\ell^m \f$.
!!
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[out] vscales: Array of scaling factors. Dimension is `(p+1)**2`
!! @param[out] vscales: Array of values 4pi/(2l+1). Dimension is `p+1`
!! @param[out] vscales_rel: Array of relative scaling factors.
!!      Dimension is `(p+1)**2`.
subroutine ylmscale(p, vscales, v4pi2lp1, vscales_rel)
    ! Input
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: vscales((p+1)**2), v4pi2lp1(p+1), &
        & vscales_rel((p+1)**2)
    ! Local variables
    real(dp) :: tmp, twolp1
    integer :: l, ind, m
    twolp1 = one
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        tmp = fourpi / twolp1
        v4pi2lp1(l+1) = tmp
        tmp = sqrt(tmp)
        vscales_rel(ind) = tmp
        vscales(ind) = one / tmp
        twolp1 = twolp1 + two
        tmp = vscales(ind) * sqrt2
        ! m != 0
        do m = 1, l
            tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
            vscales(ind+m) = tmp
            vscales(ind-m) = tmp
            vscales_rel(ind+m) = tmp * v4pi2lp1(l+1)
            vscales_rel(ind-m) = vscales_rel(ind+m)
        end do
    end do
end subroutine ylmscale

!> Compute FMM-related constants
!!
!! @param[in] dmax: Maximal degree of spherical harmonics to be evaluated.
!!      `dmax` >= 0
!! @param[in] pm: Maximal degree of the multipole expansion. `pm` >= 0.
!! @param[in] pl: Maximal degree of the local expansion. `pl` >= 0.
!! @param[out] vcnk: Array of squre roots of combinatorial factors C_n^k.
!!      Dimension is `(2*dmax+1)*(dmax+1)`.
!! @param[out] m2l_ztranslate_coef: Constants for M2L translation over OZ axis.
!!      Dimension is `(pm+1, pl+1, pl+1)`.
!! @param[out] m2l_ztranslate_coef: Constants for adjoint M2L translation over
!!      OZ axis. Dimension is `(pl+1, pl+1, pm+1)`.
subroutine fmm_constants(dmax, pm, pl, vcnk, m2l_ztranslate_coef, &
        & m2l_ztranslate_adj_coef)
    ! Inputs
    integer, intent(in) :: dmax, pm, pl
    ! Outputs
    real(dp), intent(out) :: vcnk((2*dmax+1)*(dmax+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1)
    ! Local variables
    integer :: i, indi, j, k, n, indjn
    real(dp) :: tmp1
    ! Compute combinatorial numbers C_n^k for n=0..dmax
    ! C_0^0 = 1
    vcnk(1) = one
    do i = 2, 2*dmax+1
        ! Offset to the C_{i-2}^{i-2}, next item to be stored is C_{i-1}^0
        indi = (i-1) * i / 2
        ! C_{i-1}^0 = 1
        vcnk(indi+1) = one
        ! C_{i-1}^{i-1} = 1
        vcnk(indi+i) = one
        ! C_{i-1}^{j-1} = C_{i-2}^{j-1} + C_{i-2}^{j-2}
        ! Offset to C_{i-3}^{i-3} is indi-i+1
        do j = 2, i-1
            vcnk(indi+j) = vcnk(indi-i+j+1) + vcnk(indi-i+j)
        end do
    end do
    ! Get square roots of C_n^k. sqrt(one) is one, so no need to update C_n^0
    ! and C_n^n
    do i = 3, 2*dmax+1
        indi = (i-1) * i / 2
        do j = 2, i-1
            vcnk(indi+j) = sqrt(vcnk(indi+j))
        end do
    end do
    ! Fill in m2l_ztranslate_coef and m2l_ztranslate_adj_coef
    do j = 0, pl
        do k = 0, j
            tmp1 = one
            do n = k, pm
                indjn = (j+n)*(j+n+1)/2 + 1
                m2l_ztranslate_coef(n-k+1, k+1, j-k+1) = &
                    & tmp1 * vcnk(indjn+j-k) * vcnk(indjn+j+k)
                m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1) = &
                    & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)
                tmp1 = -tmp1
            end do
        end do
    end do
end subroutine fmm_constants

!> Convert input cartesian coordinate into spherical coordinate
!!
!! Output coordinate \f$ (\rho, \theta, \phi) \f$ is presented by \f$ (\rho,
!! \cos \theta, \sin \theta, \cos \phi, \sin\phi) \f$.
!!
!! @param[in] x: Cartesian coordinate
!! @param[out] rho: \f$ \rho \f$
!! @param[out] ctheta: \f$ \cos \theta \f$
!! @param[out] stheta: \f$ \sin \theta \f$
!! @param[out] cphi: \f$ \cos \phi \f$
!! @param[out] sphi: \f$ \sin \phi \f$
subroutine carttosph(x, rho, ctheta, stheta, cphi, sphi)
    ! Input
    real(dp), intent(in) :: x(3)
    ! Output
    real(dp), intent(out) :: rho, ctheta, stheta, cphi, sphi
    ! Local variables
    real(dp) :: max12, ssq12
    ! Check x(1:2) = 0
    if ((x(1) .eq. zero) .and. (x(2) .eq. zero)) then
        rho = abs(x(3))
        ctheta = sign(one, x(3))
        stheta = zero
        cphi = one
        sphi = zero
        return
    end if
    ! In other situations use sum-of-scaled-squares technique
    ! Get norm of x(1:2) and cphi with sphi outputs
    if (abs(x(2)) .gt. abs(x(1))) then
        max12 = abs(x(2))
        ssq12 = one + (x(1)/x(2))**2
    else
        max12 = abs(x(1))
        ssq12 = one + (x(2)/x(1))**2
    end if
    stheta = max12 * sqrt(ssq12)
    cphi = x(1) / stheta
    sphi = x(2) / stheta
    ! Then compute rho, ctheta and stheta outputs
    if (abs(x(3)) .gt. max12) then
        rho = one + ssq12*(max12/x(3))**2
        rho = abs(x(3)) * sqrt(rho)
        stheta = stheta / rho
        ctheta = x(3) / rho
    else
        rho = ssq12 + (x(3)/max12)**2
        rho = max12 * sqrt(rho)
        stheta = stheta / rho
        ctheta = x(3) / rho
    end if
end subroutine carttosph

!> Compute all spherical harmonics up to a given degree at a given point
!!
!! Spherical harmonics are computed for a point \f$ x / \|x\| \f$. Cartesian
!! coordinate of input `x` is translated into a spherical coordinate \f$ (\rho,
!! \theta, \phi) \f$ that is represented by \f$ \rho, \cos \theta, \sin \theta,
!! \cos \phi \f$ and \f$ \sin \phi \f$. If \f$ \rho=0 \f$ nothing is computed,
!! only zero \f$ \rho \f$ is returned without doing anything else. If \f$
!! \rho>0 \f$ values \f$ \cos \theta \f$ and \f$ \sin \theta \f$ are computed.
!! If \f$ \sin \theta \ne 0 \f$ then \f$ \cos \phi \f$ and \f$ \sin \phi \f$
!! are computed.
!! Auxiliary values of associated Legendre polynomials \f$ P_\ell^m(\theta) \f$
!! are computed along with \f$ \cos (m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!! @param[in] x: Target point
!! @param[out] rho: Euclidian length of `x`
!! @param[out] ctheta: \f$ -1 \leq \cos \theta \leq 1\f$
!! @param[out] stheta: \f$ 0 \leq \sin \theta \leq 1\f$
!! @param[out] cphi: \f$ -1 \leq \cos \phi \leq 1\f$
!! @param[out] sphi: \f$ -1 \leq \sin \phi \leq 1\f$
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[in] vscales: Scaling factors of real normalized spherical harmonics.
!!      Dimension is `(p+1)**2`
!! @param[out] vylm: Values of spherical harmonics \f$ Y_\ell^m(x) \f$.
!!      Dimension is `(p+1)**2`
!! @param[out] vplm: Values of associated Legendre polynomials \f$ P_\ell^m(
!!      \theta) \f$. Dimension is `(p+1)**2`
!! @param[out] vcos: Array of alues of \f$ \cos(m\phi) \f$ of a dimension
!!      `(p+1)`
!! @param[out] vsin: array of values of \f$ \sin(m\phi) \f$ of a dimension
!!      `(p+1)`
subroutine ylmbas(x, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: x(3)
    real(dp), intent(in) :: vscales((p+1)**2)
    ! Outputs
    real(dp), intent(out) :: rho, ctheta, stheta, cphi, sphi
    real(dp), intent(out) :: vylm((p+1)**2), vplm((p+1)**2)
    real(dp), intent(out) :: vcos(p+1), vsin(p+1)
    ! Local variables
    integer :: l, m, ind
    real(dp) :: max12, ssq12, tmp
    ! Get rho cos(theta), sin(theta), cos(phi) and sin(phi) from the cartesian
    ! coordinates of x. To support full range of inputs we do it via a scale
    ! and a sum of squares technique.
    ! At first we compute x(1)**2 + x(2)**2
    if (x(1) .eq. zero) then
        max12 = abs(x(2))
        ssq12 = one
    else if (abs(x(2)) .gt. abs(x(1))) then
        max12 = abs(x(2))
        ssq12 = one + (x(1)/x(2))**2
    else
        max12 = abs(x(1))
        ssq12 = one + (x(2)/x(1))**2
    end if
    ! Then we compute rho
    if (x(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(x(3)) .gt. max12) then
        rho = one + ssq12 *(max12/x(3))**2
        rho = abs(x(3)) * sqrt(rho)
    else
        rho = ssq12 + (x(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case x=0 just exit without setting any other variable
    if (rho .eq. zero) then
        return
    end if
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = x(1) / stheta
        sphi = x(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = x(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                vylm(1) = vscales(1)
                return
            case (1)
                vylm(1) = vscales(1)
                vylm(2) = - vscales(4) * stheta * sphi
                vylm(3) = vscales(3) * ctheta
                vylm(4) = - vscales(4) * stheta * cphi
                return
        end select
        ! Evaluate associated Legendre polynomials, size of temporary workspace
        ! here needs only (p+1) elements so we use vcos as a temporary
        call polleg_work(ctheta, stheta, p, vplm, vcos)
        ! Evaluate cos(m*phi) and sin(m*phi) arrays
        call trgev(cphi, sphi, p, vcos, vsin)
        ! Construct spherical harmonics
        ! l = 0
        vylm(1) = vscales(1)
        ! l = 1
        vylm(2) = - vscales(4) * stheta * sphi
        vylm(3) = vscales(3) * ctheta
        vylm(4) = - vscales(4) * stheta * cphi
        ind = 3
        do l = 2, p
            ! Offset of a Y_l^0 harmonic in vplm and vylm arrays
            ind = ind + 2*l
            !l**2 + l + 1
            ! m = 0 implicitly uses `vcos(1) = 1`
            vylm(ind) = vscales(ind) * vplm(ind)
            do m = 1, l
                ! only P_l^m for non-negative m is used/defined
                tmp = vplm(ind+m) * vscales(ind+m)
                ! m > 0
                vylm(ind+m) = tmp * vcos(m+1)
                ! m < 0
                vylm(ind-m) = tmp * vsin(m+1)
            end do
        end do
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! Set spherical coordinates
        cphi = one
        sphi = zero
        ctheta = sign(one, x(3))
        stheta = zero
        ! Set output arrays vcos and vsin
        vcos = one
        vsin = zero
        ! Evaluate spherical harmonics. P_l^m = 0 for m > 0. In the case m = 0
        ! it depends if l is odd or even. Additionally, vcos = one and vsin =
        ! zero for all elements
        vylm = zero
        vplm = zero
        do l = 0, p, 2
            ind = l**2 + l + 1
            ! only case m = 0
            vplm(ind) = one
            vylm(ind) = vscales(ind)
        end do
        do l = 1, p, 2
            ind = l**2 + l + 1
            ! only case m = 0
            vplm(ind) = ctheta
            vylm(ind) = ctheta * vscales(ind)
        end do
    end if
end subroutine ylmbas

!> Compute arrays of \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$
!!
!! All values are computed recurrently from input \f$ \cos(\phi) \f$ and \f$
!! \sin(\phi) \f$ without accessing arccos or arcsin functions.
!!
!! @param[in] cphi: \f$ \cos(\phi) \f$. -1 <= `cphi` <= 1
!! @param[in] sphi: \f$ \sin(\phi) \f$. -1 <= `sphi` <= 1
!! @param[in] p: Maximal value of \f$ m \f$, for which to compute \f$ \cos(m
!!      \phi) \f$ and \f$ \sin(m\phi) \f$. `p` >= 0
!! @param[out] vcos: Array of \f$ \cos(m\phi) \f$ for \f$ m=0..p \f$. Dimension
!!      is `(p+1)`
!! @param[out] vsin: Array of \f$ \sin(m\phi) \f$ for \f$ m=0..p \f$. Dimension
!!      is `(p+1)`
subroutine trgev(cphi, sphi, p, vcos, vsin)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: cphi, sphi
    ! Output
    real(dp), intent(out) :: vcos(p+1), vsin(p+1)
    ! Local variables
    integer :: m
    real(dp) :: c4phi, s4phi
    !! Treat values of p from 0 to 3 differently
    select case(p)
        ! Do nothing if p < 0
        case (:-1)
            return
        ! p = 0
        case (0)
            vcos(1) = one
            vsin(1) = zero
            return
        ! p = 1
        case (1)
            vcos(1) = one
            vsin(1) = zero
            vcos(2) = cphi
            vsin(2) = sphi
            return
        ! p = 2
        case (2)
            vcos(1) = one
            vsin(1) = zero
            vcos(2) = cphi
            vsin(2) = sphi
            vcos(3) = cphi**2 - sphi**2
            vsin(3) = 2 * cphi * sphi
            return
        ! p = 3
        case (3)
            vcos(1) = one
            vsin(1) = zero
            vcos(2) = cphi
            vsin(2) = sphi
            vcos(3) = cphi**2 - sphi**2
            vsin(3) = 2 * cphi * sphi
            vcos(4) = vcos(3)*cphi - vsin(3)*sphi
            vsin(4) = vcos(3)*sphi + vsin(3)*cphi
            return
        ! p >= 4
        case default
            vcos(1) = one
            vsin(1) = zero
            vcos(2) = cphi
            vsin(2) = sphi
            vcos(3) = cphi**2 - sphi**2
            vsin(3) = 2 * cphi * sphi
            vcos(4) = vcos(3)*cphi - vsin(3)*sphi
            vsin(4) = vcos(3)*sphi + vsin(3)*cphi
            vcos(5) = vcos(3)**2 - vsin(3)**2
            vsin(5) = 2 * vcos(3) * vsin(3)
            c4phi = vcos(5)
            s4phi = vsin(5)
    end select
    ! Define cos(m*phi) and sin(m*phi) recurrently 4 values at a time
    do m = 6, p-2, 4
        vcos(m:m+3) = vcos(m-4:m-1)*c4phi - vsin(m-4:m-1)*s4phi
        vsin(m:m+3) = vcos(m-4:m-1)*s4phi + vsin(m-4:m-1)*c4phi
    end do
    ! Work with leftover
    select case(m-p)
        case (-1)
            vcos(p-1) = vcos(p-2)*cphi - vsin(p-2)*sphi
            vsin(p-1) = vcos(p-2)*sphi + vsin(p-2)*cphi
            vcos(p) = vcos(p-2)*vcos(3) - vsin(p-2)*vsin(3)
            vsin(p) = vcos(p-2)*vsin(3) + vsin(p-2)*vcos(3)
            vcos(p+1) = vcos(p-2)*vcos(4) - vsin(p-2)*vsin(4)
            vsin(p+1) = vcos(p-2)*vsin(4) + vsin(p-2)*vcos(4)
        case (0)
            vcos(p) = vcos(p-1)*cphi - vsin(p-1)*sphi
            vsin(p) = vcos(p-1)*sphi + vsin(p-1)*cphi
            vcos(p+1) = vcos(p-1)*vcos(3) - vsin(p-1)*vsin(3)
            vsin(p+1) = vcos(p-1)*vsin(3) + vsin(p-1)*vcos(3)
        case (1)
            vcos(p+1) = vcos(p)*cphi - vsin(p)*sphi
            vsin(p+1) = vcos(p)*sphi + vsin(p)*cphi
    end select
end subroutine trgev

!> Compute associated Legendre polynomials
!!
!! Only polynomials \f$ P_\ell^m (\cos \theta) \f$ with non-negative parameter
!! \f$ m \f$ are computed. Implemented via following recurrent formulas:
!! \f{align}{
!!      &P_0^0(\cos \theta) = 1\\
!!      &P_{m+1}^{m+1}(\cos \theta) = -(2m+1) \sin \theta P_m^m(\cos \theta) \\
!!      &P_{m+1}^m(\cos \theta) = \cos \theta (2m+1) P_m^m(\cos \theta) \\
!!      &P_\ell^m(\cos \theta) = \frac{1}{\ell-m} \left( \cos \theta (2\ell-1)
!!      P_{\ell-1}^m(\cos \theta) - (\ell+m-1)P_{\ell-2}^m(\cos \theta)
!!      \right), \quad \forall \ell \geq m+2.
!! \f}
!!
!! @param[in] ctheta: \f$ \cos(\theta) \f$. -1 <= `ctheta` <= 1
!! @param[in] stheta: \f$ \sin(\theta) \f$. 0 <= `stheta` <= 1
!! @param[in] p: Maximal degree of polynomials to compute. `p` >= 0
!! @param[out] vplm: Values of associated Legendre polynomials. Dimension is
!!      `(p+1)**2`
subroutine polleg(ctheta, stheta, p, vplm)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: ctheta, stheta
    ! Outputs
    real(dp), intent(out) :: vplm((p+1)**2)
    ! Temporary workspace
    real(dp) :: work(p+1)
    ! Call corresponding work routine
    call polleg_work(ctheta, stheta, p, vplm, work)
end subroutine polleg

!> Compute associated Legendre polynomials
!!
!! Only polynomials \f$ P_\ell^m (\cos \theta) \f$ with non-negative parameter
!! \f$ m \f$ are computed. Implemented via following recurrent formulas:
!! \f{align}{
!!      &P_0^0(\cos \theta) = 1\\
!!      &P_{m+1}^{m+1}(\cos \theta) = -(2m+1) \sin \theta P_m^m(\cos \theta) \\
!!      &P_{m+1}^m(\cos \theta) = \cos \theta (2m+1) P_m^m(\cos \theta) \\
!!      &P_\ell^m(\cos \theta) = \frac{1}{\ell-m} \left( \cos \theta (2\ell-1)
!!      P_{\ell-1}^m(\cos \theta) - (\ell+m-1)P_{\ell-2}^m(\cos \theta)
!!      \right), \quad \forall \ell \geq m+2.
!! \f}
!!
!! @param[in] ctheta: \f$ \cos(\theta) \f$. -1 <= `ctheta` <= 1
!! @param[in] stheta: \f$ \sin(\theta) \f$. 0 <= `stheta` <= 1
!! @param[in] p: Maximal degree of polynomials to compute. `p` >= 0
!! @param[out] vplm: Values of associated Legendre polynomials. Dimension is
!!      `(p+1)**2`
subroutine polleg_work(ctheta, stheta, p, vplm, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: ctheta, stheta
    ! Outputs
    real(dp), intent(out) :: vplm((p+1)**2)
    ! Temporary workspace
    real(dp), intent(out) :: work(p+1)
    ! Local variables
    integer :: m, l, indlm1, indlm2, indl
    real(dp) :: tmp, tmp1, tmp2, tmp3, tmp4, tmp5, fl, pl2m, pl1m, plm, pmm
    ! Easy cases
    select case (p)
        case (0)
            vplm(1) = one
            return
        case (1)
            vplm(1) = one
            vplm(3) = ctheta
            vplm(4) = -stheta
            return
        case (2)
            vplm(1) = one
            vplm(3) = ctheta
            vplm(4) = -stheta
            tmp1 = three * stheta
            tmp2 = ctheta * ctheta
            vplm(7) = 1.5d0*tmp2 - pt5
            vplm(8) = -tmp1 * ctheta
            vplm(9) = tmp1 * stheta
            return
        case (3)
            vplm(1) = one
            vplm(3) = ctheta
            vplm(4) = -stheta
            tmp1 = three * stheta
            tmp2 = ctheta * ctheta
            vplm(7) = 1.5d0*tmp2 - pt5
            vplm(8) = -tmp1 * ctheta
            vplm(9) = tmp1 * stheta
            tmp3 = 2.5d0 * tmp2 - pt5
            vplm(13) = tmp3*ctheta - ctheta
            vplm(14) = -tmp1 * tmp3
            tmp4 = -5d0 * stheta
            vplm(15) = tmp4 * vplm(8)
            vplm(16) = tmp4 * vplm(9)
            return
    end select
    ! Now p >= 4
    ! Precompute temporary values
    do l = 1, p
        work(l) = dble(2*l-1) * ctheta
    end do
    ! Save P_m^m for the inital loop m=0 which is P_0^0 now
    pmm = one
    ! Case m < p
    do m = 0, p-1
        ! P_{l-2}^m which is P_m^m now
        pl2m = pmm
        ! P_m^m
        vplm((m+1)**2) = pmm
        ! Temporary to reduce number of operations
        tmp1 = dble(2*m+1) * pmm
        ! Update P_m^m for the next iteration
        pmm = -stheta * tmp1
        ! P_{l-1}^m which is P_{m+1}^m now
        pl1m = ctheta * tmp1
        ! P_{m+1}^m
        vplm((m+1)*(m+3)) = pl1m
        ! P_l^m for l>m+1
        do l = m+2, p
            plm = work(l)*pl1m - dble(l+m-1)*pl2m
            plm = plm / dble(l-m)
            vplm(l*l+l+1+m) = plm
            ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
            pl2m = pl1m
            pl1m = plm
        end do
    end do
    ! Case m=p requires only to store P_m^m
    vplm((p+1)**2) = pmm
end subroutine polleg_work

!
! Subroutine to compute the Modified Spherical Bessel function of the first kind
! @param[in]  lmax     : Data type
! @param[in]  argument : Argument of Bessel function
! @param[out] SI       : Modified Bessel function of the first kind
! @param[out] DI       : Derivative of modified Bessel function of the first kind
subroutine modified_spherical_bessel_first_kind(lmax, argument, SI, DI, work)
    use Complex_Bessel
    !! Inputs
    integer, intent(in) :: lmax
    real(dp), intent(in) :: argument
    !! Outputs
    real(dp), dimension(0:lmax), intent(out) :: SI, DI
    !! Temporary workspace
    complex(dp), dimension(0:max(1, lmax)), intent(out) :: work
    !! Local Variables
    ! Complex_SI     : Modified Bessel functions of the first kind
    ! argument       : Argument for Complex_SI
    ! scaling_factor : sqrt(pi/2x)
    ! fnu            : Starting argument for I_J(x)
    ! NZ             : Number of components set to zero due to underflow
    ! ierr           : Erro indicator for I_J(x)
    ! l              : l=0,...,lmax
    complex(dp) :: complex_argument
    real(dp) :: scaling_factor
    real(dp), parameter :: fnu = 1.5d0
    integer :: NZ, ierr, l
    ! If argument is zero then SI is e_1, need to find what is the value of DI
    if (argument .eq. zero) then
        SI = zero
        SI(0) = one
        DI = zero ! ???
        return
    end if
    ! Compute I_(0.5+J) for J:1,...,lmax
    ! NOTE: cbesi computes I_(FNU+J-1)
    !fnu = 1.5
    scaling_factor = sqrt(PI/(2*argument))
    ! NOTE: Complex argument is required to call I_J(x)
    complex_argument = argument
    ! Compute for l = 0
    call cbesi(complex_argument, pt5, 1, 1, work(0), NZ, ierr)
    if (ierr .ne. 0) stop 'Error in computing Bessel function of first kind'
    ! Compute for l = 1,...,lmax
    if (lmax .gt. 0) then
        call cbesi(complex_argument, fnu, 1, lmax, work(1:lmax), NZ, ierr)
        if (ierr .ne. 0) stop 'Error in computing Bessel function of first kind'
    ! l=1 is needed for DI(0)
    else
        call cbesi(complex_argument, fnu, 1, 1, work(1), NZ, ierr)
        if (ierr .ne. 0) stop 'Error in computing Bessel function of first kind'
    end if
    ! Store the real part of the complex Bessel functions
    SI = real(work(0:lmax))
    ! Converting Modified Bessel to Spherical Modified Bessel
    SI = scaling_factor * SI
    ! Computation of Derivatives of SI
    DI(0) = scaling_factor * real(work(1))
    do l = 1, lmax
        DI(l)= SI(l-1) - ((l+1.0D0)*SI(l))/argument
    end do
end subroutine modified_spherical_bessel_first_kind

!
! Subroutine to compute the Modified Spherical Bessel function of the second kind
! @param[in]  lmax     : Size
! @param[in]  argument : Argument of Bessel function
! @param[out] SK       : Modified Bessel function of the second kind
! @param[out] DK       : Derivative of modified Bessel function of the second kind
subroutine modified_spherical_bessel_second_kind(lmax, argument, SK, DK, work)
    use Complex_Bessel
    !! Inputs
    integer , intent(in) :: lmax
    real(dp), intent(in) :: argument
    !! Outputs
    real(dp), dimension(0:lmax), intent(out) :: SK, DK
    !! Temporary workspace
    complex(dp), dimension(0:max(1, lmax)), intent(out) :: work
    ! Local Variables
    ! Complex_SK     : Modified Bessel functions of the second kind
    ! argument       : Argument for Complex_SK
    ! scaling_factor : sqrt(pi/2x)
    ! fnu            : Starting argument for K_J(x)
    ! NZ             : Number of components set to zero due to underflow
    ! ierr           : Error indicator for K_J(x)
    ! l              : l=0,...,lmax
    complex(dp) :: complex_argument
    real(dp) :: scaling_factor
    real(dp), parameter :: fnu=1.5d0
    integer :: NZ, ierr, l
    ! Compute K_(0.5+J) for J:1,...,lmax
    ! NOTE: cbesk computes K_(FNU+J-1)
    scaling_factor = sqrt(2/(PI*argument))
    ! NOTE: Complex argument is required to call K_J(x)
    complex_argument = argument
    ! Compute for l = 0
    call cbesk(complex_argument, pt5, 1, 1, work(0), NZ, ierr)
    if (ierr .ne. 0) stop 'Error in computing Bessel function of second kind'
    ! Compute for l = 1,...,lmax
    if (lmax .gt. 0) then
        call cbesk(complex_argument, fnu, 1, lmax, work(1:lmax), NZ, ierr)
        if (ierr .ne. 0) stop 'Error in computing Bessel function of second kind'
    ! l=1 is needed for DK(0)
    else
        call cbesk(complex_argument, fnu, 1, 1, work(1), NZ, ierr)
        if (ierr .ne. 0) stop 'Error in computing Bessel function of second kind'
    end if
    ! Store the real part of the complex Bessel functions
    SK = real(work(0:lmax))
    ! Converting Modified Bessel to Spherical Modified Bessel
    SK = scaling_factor * SK
    ! Computation of Derivatives of SK
    DK(0) = -scaling_factor * real(work(1))
    do l = 1, lmax
        DK(l) = -SK(l-1) - ((l+1.0D0)*SK(l))/argument
    end do
end subroutine modified_spherical_bessel_second_kind

!> Accumulate a multipole expansion induced by a particle of a given charge
!!
!! Computes the following sums:
!! \f[
!!      \forall \ell=0, \ldots, p, \quad \forall m=-\ell, \ldots, \ell : \quad
!!      M_\ell^m = \beta M_\ell^m + \frac{q \|c\|^\ell}{r^{\ell+1}}
!!      Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of output harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \beta \f$ is a scaling factor, \f$ q \f$ and \f$
!! c \f$ are a charge and coordinates of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{\ell+1} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^{-\ell-1} \f$.
!!
!! @param[in] c: Radius-vector from the particle to the center of harmonics
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output multipole spherical harmonics
!! @param[in] p: Maximal degree of output multipole spherical harmonics
!! @param[in] vscales: Normalization constants for spherical harmonics
!! @param[in] beta: Scaling factor for `dst_m`
!! @param[inout] dst_m: Multipole coefficients
!!
!! @sa fmm_m2p
subroutine fmm_p2m(c, src_q, dst_r, p, vscales, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Workspace
    real(dp) :: work(2*(p+1)*(p+2))
    ! Call corresponding work routine
    call fmm_p2m_work(c, src_q, dst_r, p, vscales, beta, dst_m, work)
end subroutine fmm_p2m

!> Accumulate a multipole expansion induced by a particle of a given charge
!!
!! Computes the following sums:
!! \f[
!!      \forall \ell=0, \ldots, p, \quad \forall m=-\ell, \ldots, \ell : \quad
!!      M_\ell^m = \beta M_\ell^m + \frac{q \|c\|^\ell}{r^{\ell+1}}
!!      Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of output harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \beta \f$ is a scaling factor, \f$ q \f$ and \f$
!! c \f$ are a charge and coordinates of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{\ell+1} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^{-\ell-1} \f$.
!!
!! @param[in] c: Radius-vector from the particle to the center of harmonics
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output multipole spherical harmonics
!! @param[in] p: Maximal degree of output multipole spherical harmonics
!! @param[in] vscales: Normalization constants for spherical harmonics
!! @param[in] beta: Scaling factor for `dst_m`
!! @param[inout] dst_m: Multipole coefficients
!! @param[out] work: Temporary workspace of a size (2*(p+1)*(p+2))
!!
!! @sa fmm_m2p
subroutine fmm_p2m_work(c, src_q, dst_r, p, vscales, beta, dst_m, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Workspace
    real(dp), intent(out), target :: work(2*(p+1)*(p+2))
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, t, rcoef
    integer :: n, ind, vylm1, vylm2, vplm1, vplm2, vcos1, vcos2, vsin1, vsin2
    ! In case src_q is zero just scale output properly
    if (src_q .eq. zero) then
        ! Zero init output if beta is also zero
        if (beta .eq. zero) then
            dst_m = zero
        ! Scale output by beta otherwise
        else
            dst_m = beta * dst_m
        end if
        ! Exit subroutine
        return
    end if
    ! Now src_q is non-zero
    ! Mark first and last elements of all work subarrays
    vylm1 = 1
    vylm2 = (p+1)**2
    vplm1 = vylm2 + 1
    vplm2 = 2 * vylm2
    vcos1 = vplm2 + 1
    vcos2 = vcos1 + p
    vsin1 = vcos2 + 1
    vsin2 = vsin1 + p
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, &
        & work(vylm1:vylm2), work(vplm1:vplm2), work(vcos1:vcos2), &
        & work(vsin1:vsin2))
    ! Harmonics are available only if rho > 0
    if (rho .ne. zero) then
        rcoef = rho / dst_r
        t = src_q / dst_r
        ! Ignore input `dst_m` in case of a zero scaling factor beta
        if (beta .eq. zero) then
            do n = 0, p
                ind = n*n + n + 1
                ! Array vylm is in the beginning of work array
                dst_m(ind-n:ind+n) = t * work(ind-n:ind+n)
                t = t * rcoef
            end do
        ! Update `dst_m` otherwise
        else
            do n = 0, p
                ind = n*n + n + 1
                ! Array vylm is in the beginning of work array
                dst_m(ind-n:ind+n) = beta*dst_m(ind-n:ind+n) + &
                    & t*work(ind-n:ind+n)
                t = t * rcoef
            end do
        end if
    ! Naive case of rho = 0
    else
        ! Ignore input `dst_m` in case of a zero scaling factor beta
        if (beta .eq. zero) then
            dst_m(1) = src_q / dst_r / sqrt4pi
            dst_m(2:) = zero
        ! Update `m` otherwise
        else
            dst_m(1) = beta*dst_m(1) + src_q/dst_r/sqrt4pi
            dst_m(2:) = beta * dst_m(2:)
        end if
    end if
end subroutine fmm_p2m_work

!> Accumulate potential, induced by multipole spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `v`
!! @param[inout] v: Value of induced potential
subroutine fmm_m2p(c, src_r, p, vscales_rel, alpha, src_m, beta, dst_v)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales_rel((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Temporary workspace
    real(dp) :: work(p+1)
    ! Call corresponding work routine
    call fmm_m2p_work(c, src_r, p, vscales_rel, alpha, src_m, beta, dst_v, &
        & work)
end subroutine fmm_m2p

!> Accumulate potential, induced by multipole spherical harmonics
!!
!! This function relies on a user-provided temporary workspace
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `v`
!! @param[inout] v: Value of induced potential
!! @param[out] work: Temporary workspace of size (p+1)
subroutine fmm_m2p_work(c, src_r, p, vscales_rel, alpha, src_m, &
        & beta, dst_v, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales_rel((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Workspace
    real(dp), intent(out), target :: work(p+1)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, fl, max12, ssq12, pl2m, pl1m, plm, pmm, cmphi, smphi, ylm
    integer :: l, m, indl
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    if (rho .eq. zero) then
        return
    end if
    ! Compute the actual induced potential
    rcoef = src_r / rho
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                ! l = 0
                dst_v = dst_v + alpha*rcoef*vscales_rel(1)*src_m(1)
                return
            case (1)
                ! l = 0
                tmp = src_m(1) * vscales_rel(1)
                ! l = 1
                tmp2 = ctheta * src_m(3) * vscales_rel(3)
                tmp3 = vscales_rel(4) * stheta
                tmp2 = tmp2 - tmp3*sphi*src_m(2)
                tmp2 = tmp2 - tmp3*cphi*src_m(4)
                dst_v = dst_v + alpha*rcoef*(tmp+rcoef*tmp2)
                return
        end select
        ! Now p>1
        ! Precompute alpha*rcoef^{l+1}
        work(1) = alpha * rcoef
        do l = 1, p
            work(l+1) = rcoef * work(l)
        end do
        ! Case m = 0
        ! P_{l-2}^m which is P_0^0 now
        pl2m = one
        dst_v = dst_v + work(1)*src_m(1)*vscales_rel(1)
        ! Update P_m^m for the next iteration
        pmm = -stheta
        ! P_{l-1}^m which is P_{m+1}^m now
        pl1m = ctheta
        ylm = pl1m * vscales_rel(3)
        dst_v = dst_v + work(2)*src_m(3)*ylm
        ! P_l^m for l>m+1
        do l = 2, p
            plm = dble(2*l-1)*ctheta*pl1m - dble(l-1)*pl2m
            plm = plm / dble(l)
            ylm = plm * vscales_rel(l*l+l+1)
            dst_v = dst_v + work(l+1)*src_m(l*l+l+1)*ylm
            ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
            pl2m = pl1m
            pl1m = plm
        end do
        ! Prepare cos(m*phi) and sin(m*phi) for m=1
        cmphi = cphi
        smphi = sphi
        ! Case 0<m<p
        do m = 1, p-1
            ! P_{l-2}^m which is P_m^m now
            pl2m = pmm
            ylm = pmm * vscales_rel((m+1)**2)
            tmp1 = cmphi*src_m((m+1)**2) + smphi*src_m(m*m+1)
            dst_v = dst_v + work(m+1)*ylm*tmp1
            ! Temporary to reduce number of operations
            tmp1 = dble(2*m+1) * pmm
            ! Update P_m^m for the next iteration
            pmm = -stheta * tmp1
            ! P_{l-1}^m which is P_{m+1}^m now
            pl1m = ctheta * tmp1
            ylm = pl1m * vscales_rel((m+1)*(m+3))
            tmp1 = cmphi*src_m((m+1)*(m+3)) + smphi*src_m((m+1)*(m+2)+1-m)
            dst_v = dst_v + work(m+2)*ylm*tmp1
            ! P_l^m for l>m+1
            do l = m+2, p
                plm = dble(2*l-1)*ctheta*pl1m - dble(l+m-1)*pl2m
                plm = plm / dble(l-m)
                ylm = plm * vscales_rel(l*l+l+1+m)
                tmp1 = cmphi*src_m(l*l+l+1+m) + smphi*src_m(l*l+l+1-m)
                dst_v = dst_v + work(l+1)*ylm*tmp1
                ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                pl2m = pl1m
                pl1m = plm
            end do
            ! Update cos(m*phi) and sin(m*phi) for the next iteration
            tmp1 = cmphi
            cmphi = cmphi*cphi - smphi*sphi
            smphi = tmp1*sphi + smphi*cphi
        end do
        ! Case m=p requires only to use P_m^m
        ylm = pmm * vscales_rel((p+1)**2)
        tmp1 = cmphi*src_m((p+1)**2) + smphi*src_m(p*p+1)
        dst_v = dst_v + work(p+1)*ylm*tmp1
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^0 = ctheta^l in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef. But before that we
        ! initialize alpha/r factor for l=0 case without taking possible sign
        ! into account.
        t = alpha * rcoef
        if (c(3) .lt. zero) then
            rcoef = -rcoef
        end if
        ! Proceed with accumulation of a potential
        indl = 1
        do l = 0, p
            ! Index of Y_l^0
            indl = indl + 2*l
            ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
            dst_v = dst_v + t*src_m(indl)*vscales_rel(indl)
            ! Update t
            t = t * rcoef
        end do
    end if
end subroutine fmm_m2p_work

!> Accumulate potential, induced by multipole spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales: Normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `v`
!! @param[inout] v: Value of induced potential
subroutine fmm_m2p_bessel(c, src_r, p, vscales, alpha, src_m, beta, dst_v)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Temporary workspace
    complex(dp) :: work_complex(max(2, p+1))
    real(dp) :: work(p+1)
    ! Local variables
    real(dp) :: src_sk(p+1)
    ! Call corresponding work routine
    call modified_spherical_bessel_second_kind(p, src_r, src_sk, work, &
        & work_complex)
    call fmm_m2p_bessel_work(c, p, vscales, src_sk, alpha, src_m, beta, &
        & dst_v, work_complex, work)
end subroutine fmm_m2p_bessel

!> Accumulate Bessel potential, induced by multipole spherical harmonics
!!
!! This function relies on a user-provided temporary workspace
!!
!! Computes the following sum:
!! \f[
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales: Normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `v`
!! @param[inout] v: Value of induced potential
!! @param[out] work: Temporary workspace of size (p+1)
subroutine fmm_m2p_bessel_work(c, p, vscales, SK_ri, alpha, src_m, beta, &
        & dst_v, work_complex, work)
    use complex_bessel
    !! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), vscales((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta, SK_ri(p+1)
    !! Output
    real(dp), intent(inout) :: dst_v
    !! Workspace
    complex(dp), intent(out) :: work_complex(p+1)
    real(dp), intent(out) :: work(p+1)
    !! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, tmp, tmp1, tmp2, &
        & tmp3, fl, max12, ssq12, pl2m, pl1m, plm, pmm, cmphi, smphi, ylm, t
    real(dp) :: scaling_factor
    complex(dp) :: complex_argument
    integer :: l, m, indl, ierr, NZ
    real(dp), parameter :: fnu=1.5d0
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    if (rho .eq. zero) then
        return
    end if
    ! Compute Bessel function
    scaling_factor = alpha * sqrt(two / (pi*rho))
    ! NOTE: Complex argument is required to call I_J(x)
    complex_argument = rho
    ! Compute for l = 0
    call cbesk(complex_argument, pt5, 1, 1, work_complex(1), NZ, ierr)
    ! Compute for l = 1,...,p
    if (p .gt. 0) then
        call cbesk(complex_argument, fnu, 1, p, work_complex(2:p+1), NZ, ierr)
    end if
    work = scaling_factor * real(work_complex) / SK_ri
    ! Compute the actual induced potential
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                ! l = 0
                dst_v = dst_v + work(1)*vscales(1)*src_m(1)
                return
            case (1)
                ! l = 0
                tmp = work(1) * src_m(1) * vscales(1)
                ! l = 1
                tmp2 = ctheta * src_m(3) * vscales(3)
                tmp3 = stheta * vscales(4)
                tmp2 = tmp2 - tmp3*sphi*src_m(2)
                tmp2 = tmp2 - tmp3*cphi*src_m(4)
                dst_v = dst_v + tmp + work(2)*tmp2
                return
        end select
        ! Now p>1
        ! Case m = 0
        ! P_{l-2}^m which is P_0^0 now
        pl2m = one
        dst_v = dst_v + work(1)*src_m(1)*vscales(1)
        ! Update P_m^m for the next iteration
        pmm = -stheta
        ! P_{l-1}^m which is P_{m+1}^m now
        pl1m = ctheta
        ylm = pl1m * vscales(3)
        dst_v = dst_v + work(2)*src_m(3)*ylm
        ! P_l^m for l>m+1
        do l = 2, p
            plm = dble(2*l-1)*ctheta*pl1m - dble(l-1)*pl2m
            plm = plm / dble(l)
            ylm = plm * vscales(l*l+l+1)
            dst_v = dst_v + work(l+1)*src_m(l*l+l+1)*ylm
            ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
            pl2m = pl1m
            pl1m = plm
        end do
        ! Prepare cos(m*phi) and sin(m*phi) for m=1
        cmphi = cphi
        smphi = sphi
        ! Case 0<m<p
        do m = 1, p-1
            ! P_{l-2}^m which is P_m^m now
            pl2m = pmm
            ylm = pmm * vscales((m+1)**2)
            tmp1 = cmphi*src_m((m+1)**2) + smphi*src_m(m*m+1)
            dst_v = dst_v + work(m+1)*ylm*tmp1
            ! Temporary to reduce number of operations
            tmp1 = dble(2*m+1) * pmm
            ! Update P_m^m for the next iteration
            pmm = -stheta * tmp1
            ! P_{l-1}^m which is P_{m+1}^m now
            pl1m = ctheta * tmp1
            ylm = pl1m * vscales((m+1)*(m+3))
            tmp1 = cmphi*src_m((m+1)*(m+3)) + smphi*src_m((m+1)*(m+2)+1-m)
            dst_v = dst_v + work(m+2)*ylm*tmp1
            ! P_l^m for l>m+1
            do l = m+2, p
                plm = dble(2*l-1)*ctheta*pl1m - dble(l+m-1)*pl2m
                plm = plm / dble(l-m)
                ylm = plm * vscales(l*l+l+1+m)
                tmp1 = cmphi*src_m(l*l+l+1+m) + smphi*src_m(l*l+l+1-m)
                dst_v = dst_v + work(l+1)*ylm*tmp1
                ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                pl2m = pl1m
                pl1m = plm
            end do
            ! Update cos(m*phi) and sin(m*phi) for the next iteration
            tmp1 = cmphi
            cmphi = cmphi*cphi - smphi*sphi
            smphi = tmp1*sphi + smphi*cphi
        end do
        ! Case m=p requires only to use P_m^m
        ylm = pmm * vscales((p+1)**2)
        tmp1 = cmphi*src_m((p+1)**2) + smphi*src_m(p*p+1)
        dst_v = dst_v + work(p+1)*ylm*tmp1
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^0 = ctheta^l in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef. But before that we
        ! initialize alpha/r factor for l=0 case without taking possible sign
        ! into account.
        t = one
        if (c(3) .lt. zero) then
            ! Proceed with accumulation of a potential
            indl = 1
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Y_l^0 contribution
                dst_v = dst_v + t*work(l+1)*src_m(indl)*vscales(indl)
                t = -t
            end do
        else
            ! Proceed with accumulation of a potential
            indl = 1
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Y_l^0 contribution
                dst_v = dst_v + work(l+1)*src_m(indl)*vscales(indl)
            end do
        end if
    end if
end subroutine fmm_m2p_bessel_work

subroutine fmm_m2p_bessel_grad(c, src_r, p, vscales, alpha, src_m, beta, dst_g)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales((p+2)*(p+2)), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_g(3)
    ! Temporary workspace
    complex(dp) :: work_complex(p+2)
    real(dp) :: work(p+2), src_sk(p+2), vcnk(1), src_m_grad((p+2)**2, 3)
    ! Call corresponding work routine
    call modified_spherical_bessel_second_kind(p+1, src_r, src_sk, work, &
        & work_complex)
    call fmm_m2m_bessel_grad(p, src_sk, vscales, vcnk, src_m, src_m_grad)
    call fmm_m2p_bessel_work(c, p+1, vscales, src_sk, -alpha, src_m_grad(:, 1), &
        & beta, dst_g(1), work_complex, work)
    call fmm_m2p_bessel_work(c, p+1, vscales, src_sk, -alpha, src_m_grad(:, 2), &
        & beta, dst_g(2), work_complex, work)
    call fmm_m2p_bessel_work(c, p+1, vscales, src_sk, -alpha, src_m_grad(:, 3), &
        & beta, dst_g(3), work_complex, work)
end subroutine fmm_m2p_bessel_grad

!> Adjoint M2P operation
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output multipole spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole coefficients. Dimension is `(p+1)**2`
subroutine fmm_m2p_adj(c, src_q, dst_r, p, vscales_rel, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales_rel((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Workspace
    real(dp) :: work(p+1)
    ! Call corresponding work routine
    call fmm_m2p_adj_work(c, src_q, dst_r, p, vscales_rel, beta, dst_m, work)
end subroutine fmm_m2p_adj

!> Adjoint M2P operation
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output multipole spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[out] work: Temporary workspace of a size ((p+1)*(p+1)+3*p)
subroutine fmm_m2p_adj_work(c, src_q, dst_r, p, vscales_rel, beta, dst_m, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales_rel((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Workspace
    real(dp), intent(out), target :: work(p+1)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, fl, max12, ssq12, pl2m, pl1m, plm, pmm, cmphi, smphi, ylm
    integer :: l, m, indl
    ! In case src_q is zero just scale output properly
    if (src_q .eq. zero) then
        ! Zero init output if beta is also zero
        if (beta .eq. zero) then
            dst_m = zero
        ! Scale output by beta otherwise
        else
            dst_m = beta * dst_m
        end if
        ! Exit subroutine
        return
    end if
    ! Now src_q is non-zero
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    if (rho .eq. zero) then
        return
    end if
    ! Compute actual induced potentials
    rcoef = dst_r / rho
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        t = src_q * rcoef
        select case(p)
            case (0)
                if (beta .eq. zero) then
                    dst_m(1) = t * vscales_rel(1)
                else
                    dst_m(1) = beta*dst_m(1) + t*vscales_rel(1)
                end if
                return
            case (1)
                if (beta .eq. zero) then
                    ! l = 0
                    dst_m(1) = t * vscales_rel(1)
                    ! l = 1
                    t = t * rcoef
                    tmp = -vscales_rel(4) * stheta
                    tmp2 = t * tmp
                    dst_m(2) = tmp2 * sphi
                    dst_m(3) = t * vscales_rel(3) * ctheta
                    dst_m(4) = tmp2 * cphi
                else
                    ! l = 0
                    dst_m(1) = beta*dst_m(1) + t*vscales_rel(1)
                    ! l = 1
                    t = t * rcoef
                    tmp = -vscales_rel(4) * stheta
                    tmp2 = t * tmp
                    dst_m(2) = beta*dst_m(2) + tmp2*sphi
                    dst_m(3) = beta*dst_m(3) + t*vscales_rel(3)*ctheta
                    dst_m(4) = beta*dst_m(4) + tmp2*cphi
                end if
                return
        end select
        ! Now p>1
        ! Precompute src_q*rcoef^{l+1}
        work(1) = src_q * rcoef
        do l = 1, p
            work(l+1) = rcoef * work(l)
        end do
        ! Overwrite output
        if (beta .eq. zero) then
            ! Case m = 0
            ! P_{l-2}^m which is P_0^0 now
            pl2m = one
            dst_m(1) = work(1) * vscales_rel(1)
            ! Update P_m^m for the next iteration
            pmm = -stheta
            ! P_{l-1}^m which is P_{m+1}^m now
            pl1m = ctheta
            ylm = pl1m * vscales_rel(3)
            dst_m(3) = work(2) * ylm
            ! P_l^m for l>m+1
            do l = 2, p
                plm = dble(2*l-1)*ctheta*pl1m - dble(l-1)*pl2m
                plm = plm / dble(l)
                ylm = plm * vscales_rel(l*l+l+1)
                dst_m(l*l+l+1) = work(l+1) * ylm
                ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                pl2m = pl1m
                pl1m = plm
            end do
            ! Prepare cos(m*phi) and sin(m*phi) for m=1
            cmphi = cphi
            smphi = sphi
            ! Case 0<m<p
            do m = 1, p-1
                ! P_{l-2}^m which is P_m^m now
                pl2m = pmm
                ylm = pmm * vscales_rel((m+1)**2)
                ylm = work(m+1) * ylm
                dst_m((m+1)**2) = cmphi * ylm
                dst_m(m*m+1) = smphi * ylm
                ! Temporary to reduce number of operations
                tmp1 = dble(2*m+1) * pmm
                ! Update P_m^m for the next iteration
                pmm = -stheta * tmp1
                ! P_{l-1}^m which is P_{m+1}^m now
                pl1m = ctheta * tmp1
                ylm = pl1m * vscales_rel((m+1)*(m+3))
                ylm = work(m+2) * ylm
                dst_m((m+1)*(m+3)) = cmphi * ylm
                dst_m((m+1)*(m+2)+1-m) = smphi * ylm
                ! P_l^m for l>m+1
                do l = m+2, p
                    plm = dble(2*l-1)*ctheta*pl1m - dble(l+m-1)*pl2m
                    plm = plm / dble(l-m)
                    ylm = plm * vscales_rel(l*l+l+1+m)
                    ylm = work(l+1) * ylm
                    dst_m(l*l+l+1+m) = cmphi * ylm
                    dst_m(l*l+l+1-m) = smphi * ylm
                    ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                    pl2m = pl1m
                    pl1m = plm
                end do
                ! Update cos(m*phi) and sin(m*phi) for the next iteration
                tmp1 = cmphi
                cmphi = cmphi*cphi - smphi*sphi
                smphi = tmp1*sphi + smphi*cphi
            end do
            ! Case m=p requires only to use P_m^m
            ylm = pmm * vscales_rel((p+1)**2)
            ylm = work(p+1) * ylm
            dst_m((p+1)**2) = cmphi * ylm
            dst_m(p*p+1) = smphi * ylm
        ! Update output
        else
            ! Case m = 0
            ! P_{l-2}^m which is P_0^0 now
            pl2m = one
            dst_m(1) = beta*dst_m(1) + work(1)*vscales_rel(1)
            ! Update P_m^m for the next iteration
            pmm = -stheta
            ! P_{l-1}^m which is P_{m+1}^m now
            pl1m = ctheta
            ylm = pl1m * vscales_rel(3)
            dst_m(3) = beta*dst_m(3) + work(2)*ylm
            ! P_l^m for l>m+1
            do l = 2, p
                plm = dble(2*l-1)*ctheta*pl1m - dble(l-1)*pl2m
                plm = plm / dble(l)
                ylm = plm * vscales_rel(l*l+l+1)
                dst_m(l*l+l+1) = beta*dst_m(l*l+l+1) + work(l+1)*ylm
                ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                pl2m = pl1m
                pl1m = plm
            end do
            ! Prepare cos(m*phi) and sin(m*phi) for m=1
            cmphi = cphi
            smphi = sphi
            ! Case 0<m<p
            do m = 1, p-1
                ! P_{l-2}^m which is P_m^m now
                pl2m = pmm
                ylm = pmm * vscales_rel((m+1)**2)
                ylm = work(m+1) * ylm
                dst_m((m+1)**2) = beta*dst_m((m+1)**2) + cmphi*ylm
                dst_m(m*m+1) = beta*dst_m(m*m+1) + smphi*ylm
                ! Temporary to reduce number of operations
                tmp1 = dble(2*m+1) * pmm
                ! Update P_m^m for the next iteration
                pmm = -stheta * tmp1
                ! P_{l-1}^m which is P_{m+1}^m now
                pl1m = ctheta * tmp1
                ylm = pl1m * vscales_rel((m+1)*(m+3))
                ylm = work(m+2) * ylm
                dst_m((m+1)*(m+3)) = beta*dst_m((m+1)*(m+3)) + cmphi*ylm
                dst_m((m+1)*(m+2)+1-m) = beta*dst_m((m+1)*(m+2)+1-m) + &
                    & smphi*ylm
                ! P_l^m for l>m+1
                do l = m+2, p
                    plm = dble(2*l-1)*ctheta*pl1m - dble(l+m-1)*pl2m
                    plm = plm / dble(l-m)
                    ylm = plm * vscales_rel(l*l+l+1+m)
                    ylm = work(l+1) * ylm
                    dst_m(l*l+l+1+m) = beta*dst_m(l*l+l+1+m) + cmphi*ylm
                    dst_m(l*l+l+1-m) = beta*dst_m(l*l+l+1-m) + smphi*ylm
                    ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                    pl2m = pl1m
                    pl1m = plm
                end do
                ! Update cos(m*phi) and sin(m*phi) for the next iteration
                tmp1 = cmphi
                cmphi = cmphi*cphi - smphi*sphi
                smphi = tmp1*sphi + smphi*cphi
            end do
            ! Case m=p requires only to use P_m^m
            ylm = pmm * vscales_rel((p+1)**2)
            ylm = work(p+1) * ylm
            dst_m((p+1)**2) = beta*dst_m((p+1)**2) + cmphi*ylm
            dst_m(p*p+1) = beta*dst_m(p*p+1) + smphi*ylm
        end if
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^0 = ctheta^l in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef. But before that we
        ! initialize alpha/r factor for l=0 case without taking possible sign
        ! into account.
        t = src_q * rcoef
        if (c(3) .lt. zero) then
            rcoef = -rcoef
        end if
        ! Proceed with accumulation of a potential
        indl = 1
        if (beta .eq. zero) then
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_m(indl) = t * vscales_rel(indl)
                dst_m(indl-l:indl-1) = zero
                dst_m(indl+1:indl+l) = zero
                ! Update t
                t = t * rcoef
            end do
        else
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_m(indl) = beta*dst_m(indl) + t*vscales_rel(indl)
                dst_m(indl-l:indl-1) = beta * dst_m(indl-l:indl-1)
                dst_m(indl+1:indl+l) = beta * dst_m(indl+1:indl+l)
                ! Update t
                t = t * rcoef
            end do
        end if
    end if
end subroutine fmm_m2p_adj_work

!> Accumulate potentials, induced by each multipole spherical harmonic
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[inout] mat: Values of potentials induced by each spherical harmonic
subroutine fmm_m2p_mat(c, r, p, vscales, mat)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), r, vscales((p+1)**2)
    ! Output
    real(dp), intent(out) :: mat((p+1)**2)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), rcoef, t, tmp
    integer :: n, ind
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    if (rho .eq. zero) then
        return
    end if
    ! Compute actual induced potentials
    rcoef = r / rho
    t = one
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        mat(ind-n:ind+n) = t / vscales(ind)**2 * vylm(ind-n:ind+n)
    end do
end subroutine fmm_m2p_mat

!> Accumulate potential, induced by local spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Local coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_v`
!! @param[inout] dst_v: Value of induced potential
subroutine fmm_l2p(c, src_r, p, vscales, alpha, src_l, beta, dst_v)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales((p+1)*(p+1)), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Workspace
    real(dp) :: work(p+1)
    ! Call corresponding work routine
    call fmm_l2p_work(c, src_r, p, vscales, alpha, src_l, beta, dst_v, work)
end subroutine fmm_l2p

!> Accumulate potential, induced by local spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Local coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_v`
!! @param[inout] dst_v: Value of induced potential
!! @param[out] work: Temporary workspace of a size (p+1)
subroutine fmm_l2p_work(c, src_r, p, vscales_rel, alpha, src_l, beta, dst_v, &
        & work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales_rel((p+1)*(p+1)), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Workspace
    real(dp), intent(out), target :: work(p+1)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, fl, max12, ssq12, pl2m, pl1m, plm, pmm, cmphi, smphi, ylm
    integer :: l, m, indl
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential depends only on the
    ! leading 1-st spherical harmonic
    if (rho .eq. zero) then
        dst_v = dst_v + alpha*src_l(1)*vscales_rel(1)
        return
    end if
    ! Compute the actual induced potential
    rcoef = rho / src_r
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                ! l = 0
                dst_v = dst_v + alpha*vscales_rel(1)*src_l(1)
                return
            case (1)
                ! l = 0
                tmp = src_l(1) * vscales_rel(1)
                ! l = 1
                tmp2 = ctheta * src_l(3) * vscales_rel(3)
                tmp3 = vscales_rel(4) * stheta
                tmp2 = tmp2 - tmp3*sphi*src_l(2)
                tmp2 = tmp2 - tmp3*cphi*src_l(4)
                dst_v = dst_v + alpha*(tmp+rcoef*tmp2)
                return
        end select
        ! Now p>1
        ! Precompute alpha*rcoef^l
        work(1) = alpha
        do l = 1, p
            work(l+1) = rcoef * work(l)
        end do
        ! Case m = 0
        ! P_{l-2}^m which is P_0^0 now
        pl2m = one
        dst_v = dst_v + work(1)*src_l(1)*vscales_rel(1)
        ! Update P_m^m for the next iteration
        pmm = -stheta
        ! P_{l-1}^m which is P_{m+1}^m now
        pl1m = ctheta
        ylm = pl1m * vscales_rel(3)
        dst_v = dst_v + work(2)*src_l(3)*ylm
        ! P_l^m for l>m+1
        do l = 2, p
            plm = dble(2*l-1)*ctheta*pl1m - dble(l-1)*pl2m
            plm = plm / dble(l)
            ylm = plm * vscales_rel(l*l+l+1)
            dst_v = dst_v + work(l+1)*src_l(l*l+l+1)*ylm
            ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
            pl2m = pl1m
            pl1m = plm
        end do
        ! Prepare cos(m*phi) and sin(m*phi) for m=1
        cmphi = cphi
        smphi = sphi
        ! Case 0<m<p
        do m = 1, p-1
            ! P_{l-2}^m which is P_m^m now
            pl2m = pmm
            ylm = pmm * vscales_rel((m+1)**2)
            tmp1 = cmphi*src_l((m+1)**2) + smphi*src_l(m*m+1)
            dst_v = dst_v + work(m+1)*ylm*tmp1
            ! Temporary to reduce number of operations
            tmp1 = dble(2*m+1) * pmm
            ! Update P_m^m for the next iteration
            pmm = -stheta * tmp1
            ! P_{l-1}^m which is P_{m+1}^m now
            pl1m = ctheta * tmp1
            ylm = pl1m * vscales_rel((m+1)*(m+3))
            tmp1 = cmphi*src_l((m+1)*(m+3)) + smphi*src_l((m+1)*(m+2)+1-m)
            dst_v = dst_v + work(m+2)*ylm*tmp1
            ! P_l^m for l>m+1
            do l = m+2, p
                plm = dble(2*l-1)*ctheta*pl1m - dble(l+m-1)*pl2m
                plm = plm / dble(l-m)
                ylm = plm * vscales_rel(l*l+l+1+m)
                tmp1 = cmphi*src_l(l*l+l+1+m) + smphi*src_l(l*l+l+1-m)
                dst_v = dst_v + work(l+1)*ylm*tmp1
                ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                pl2m = pl1m
                pl1m = plm
            end do
            ! Update cos(m*phi) and sin(m*phi) for the next iteration
            tmp1 = cmphi
            cmphi = cmphi*cphi - smphi*sphi
            smphi = tmp1*sphi + smphi*cphi
        end do
        ! Case m=p requires only to use P_m^m
        ylm = pmm * vscales_rel((p+1)**2)
        tmp1 = cmphi*src_l((p+1)**2) + smphi*src_l(p*p+1)
        dst_v = dst_v + work(p+1)*ylm*tmp1
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^m = ctheta^m in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef.
        rcoef = sign(rcoef, c(3))
        ! Init t and proceed with accumulation of a potential
        t = alpha
        indl = 1
        do l = 0, p
            ! Index of Y_l^0
            indl = indl + 2*l
            ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
            dst_v = dst_v + t*src_l(indl)*vscales_rel(indl)
            ! Update t
            t = t * rcoef
        end do
    end if
end subroutine fmm_l2p_work

!> Accumulate potential, induced by local spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Local coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_v`
!! @param[inout] dst_v: Value of induced potential
subroutine fmm_l2p_bessel(c, src_r, p, vscales, alpha, src_l, beta, dst_v)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales((p+1)*(p+1)), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Temporary workspace
    complex(dp) :: work_complex(max(2, p+1))
    real(dp) :: work(p+1)
    ! Local variables
    real(dp) :: src_si(p+1)
    ! Call corresponding work routine
    call modified_spherical_bessel_first_kind(p, src_r, src_si, work, &
        & work_complex)
    call fmm_l2p_bessel_work(c, p, vscales, src_si, alpha, src_l, beta, &
        & dst_v, work_complex, work)
end subroutine fmm_l2p_bessel

!> Accumulate Bessel potential, induced by local spherical harmonics
!!
!! This function relies on a user-provided temporary workspace
!!
!! Computes the following sum:
!! \f[
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales: Normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Local coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `v`
!! @param[inout] v: Value of induced potential
!! @param[out] work: Temporary workspace of size (p+1)
subroutine fmm_l2p_bessel_work(c, p, vscales, SI_ri, alpha, src_l, beta, &
        & dst_v, work_complex, work)
    use complex_bessel
    !! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), vscales((p+1)*(p+1)), alpha, &
        & src_l((p+1)*(p+1)), beta, SI_ri(p+1)
    !! Output
    real(dp), intent(inout) :: dst_v
    !! Workspace
    complex(dp), intent(out) :: work_complex(p+1)
    real(dp), intent(out) :: work(p+1)
    !! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, fl, max12, ssq12, pl2m, pl1m, plm, pmm, cmphi, smphi, ylm
    real(dp) :: scaling_factor
    complex(dp) :: complex_argument
    integer :: l, m, indl, ierr, NZ
    real(dp), parameter :: fnu=1.5d0
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential depends only on the
    ! leading 1-st spherical harmonic and the value I_0(0)=1
    if (rho .eq. zero) then
        dst_v = dst_v + alpha*src_l(1)*vscales(1)/SI_ri(1)
        return
    end if
    ! Compute Bessel function
    scaling_factor = alpha * sqrt(pi / (2*rho))
    ! NOTE: Complex argument is required to call I_J(x)
    complex_argument = rho
    ! Compute for l = 0
    call cbesi(complex_argument, pt5, 1, 1, work_complex(1), NZ, ierr)
    ! Compute for l = 1,...,p
    if (p .gt. 0) then
        call cbesi(complex_argument, fnu, 1, p, work_complex(2:p+1), NZ, ierr)
    end if
    work = scaling_factor * real(work_complex) / SI_ri
    ! Compute the actual induced potential
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                ! l = 0
                dst_v = dst_v + work(1)*vscales(1)*src_l(1)
                return
            case (1)
                ! l = 0
                tmp = work(1) * src_l(1) * vscales(1)
                ! l = 1
                tmp2 = ctheta * src_l(3) * vscales(3)
                tmp3 = stheta * vscales(4)
                tmp2 = tmp2 - tmp3*sphi*src_l(2)
                tmp2 = tmp2 - tmp3*cphi*src_l(4)
                dst_v = dst_v + tmp + work(2)*tmp2
                return
        end select
        ! Now p>1
        ! Case m = 0
        ! P_{l-2}^m which is P_0^0 now
        pl2m = one
        dst_v = dst_v + work(1)*src_l(1)*vscales(1)
        ! Update P_m^m for the next iteration
        pmm = -stheta
        ! P_{l-1}^m which is P_{m+1}^m now
        pl1m = ctheta
        ylm = pl1m * vscales(3)
        dst_v = dst_v + work(2)*src_l(3)*ylm
        ! P_l^m for l>m+1
        do l = 2, p
            plm = dble(2*l-1)*ctheta*pl1m - dble(l-1)*pl2m
            plm = plm / dble(l)
            ylm = plm * vscales(l*l+l+1)
            dst_v = dst_v + work(l+1)*src_l(l*l+l+1)*ylm
            ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
            pl2m = pl1m
            pl1m = plm
        end do
        ! Prepare cos(m*phi) and sin(m*phi) for m=1
        cmphi = cphi
        smphi = sphi
        ! Case 0<m<p
        do m = 1, p-1
            ! P_{l-2}^m which is P_m^m now
            pl2m = pmm
            ylm = pmm * vscales((m+1)**2)
            tmp1 = cmphi*src_l((m+1)**2) + smphi*src_l(m*m+1)
            dst_v = dst_v + work(m+1)*ylm*tmp1
            ! Temporary to reduce number of operations
            tmp1 = dble(2*m+1) * pmm
            ! Update P_m^m for the next iteration
            pmm = -stheta * tmp1
            ! P_{l-1}^m which is P_{m+1}^m now
            pl1m = ctheta * tmp1
            ylm = pl1m * vscales((m+1)*(m+3))
            tmp1 = cmphi*src_l((m+1)*(m+3)) + smphi*src_l((m+1)*(m+2)+1-m)
            dst_v = dst_v + work(m+2)*ylm*tmp1
            ! P_l^m for l>m+1
            do l = m+2, p
                plm = dble(2*l-1)*ctheta*pl1m - dble(l+m-1)*pl2m
                plm = plm / dble(l-m)
                ylm = plm * vscales(l*l+l+1+m)
                tmp1 = cmphi*src_l(l*l+l+1+m) + smphi*src_l(l*l+l+1-m)
                dst_v = dst_v + work(l+1)*ylm*tmp1
                ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                pl2m = pl1m
                pl1m = plm
            end do
            ! Update cos(m*phi) and sin(m*phi) for the next iteration
            tmp1 = cmphi
            cmphi = cmphi*cphi - smphi*sphi
            smphi = tmp1*sphi + smphi*cphi
        end do
        ! Case m=p requires only to use P_m^m
        ylm = pmm * vscales((p+1)**2)
        tmp1 = cmphi*src_l((p+1)**2) + smphi*src_l(p*p+1)
        dst_v = dst_v + work(p+1)*ylm*tmp1
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^0 = ctheta^l in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef. But before that we
        ! initialize alpha/r factor for l=0 case without taking possible sign
        ! into account.
        t = one
        if (c(3) .lt. zero) then
            ! Proceed with accumulation of a potential
            indl = 1
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_v = dst_v + t*work(l+1)*src_l(indl)*vscales(indl)
                t = -t
            end do
        else
            ! Proceed with accumulation of a potential
            indl = 1
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_v = dst_v + work(l+1)*src_l(indl)*vscales(indl)
            end do
        end if
    end if
end subroutine fmm_l2p_bessel_work

subroutine fmm_l2p_bessel_grad(c, src_r, p, vscales, alpha, src_l, beta, dst_g)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, vscales((p+2)*(p+2)), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_g(3)
    ! Temporary workspace
    complex(dp) :: work_complex(p+2)
    real(dp) :: work(p+2), src_si(p+2), vcnk(1), src_l_grad((p+2)**2, 3)
    ! Call corresponding work routine
    call modified_spherical_bessel_first_kind(p+1, src_r, src_si, work, &
        & work_complex)
    call fmm_l2l_bessel_grad(p, src_si, vscales, vcnk, src_l, src_l_grad)
    call fmm_l2p_bessel_work(c, p+1, vscales, src_si, -alpha, src_l_grad(:, 1), &
        & beta, dst_g(1), work_complex, work)
    call fmm_l2p_bessel_work(c, p+1, vscales, src_si, -alpha, src_l_grad(:, 2), &
        & beta, dst_g(2), work_complex, work)
    call fmm_l2p_bessel_work(c, p+1, vscales, src_si, -alpha, src_l_grad(:, 3), &
        & beta, dst_g(3), work_complex, work)
end subroutine fmm_l2p_bessel_grad

!> Adjoint of L2P
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output local spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Local coefficients. Dimension is `(p+1)**2`
subroutine fmm_l2p_adj(c, src_q, dst_r, p, vscales_rel, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales_rel((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)**2)
    ! Workspace
    real(dp) :: work(p+1)
    ! Call corresponding work routine
    call fmm_l2p_adj_work(c, src_q, dst_r, p, vscales_rel, beta, dst_l, work)
end subroutine fmm_l2p_adj

!> Adjoint L2P
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output local spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Local coefficients. Dimension is `(p+1)**2`
!! @param[out] work: Temporary workspace of a size (p+1)
subroutine fmm_l2p_adj_work(c, src_q, dst_r, p, vscales_rel, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales_rel((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)**2)
    ! Workspace
    real(dp), intent(out), target :: work(p+1)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, fl, max12, ssq12, pl2m, pl1m, plm, pmm, cmphi, smphi, ylm
    integer :: l, m, indl
    ! In case src_q is zero just scale output properly
    if (src_q .eq. zero) then
        ! Zero init output if beta is also zero
        if (beta .eq. zero) then
            dst_l = zero
        ! Scale output by beta otherwise
        else
            dst_l = beta * dst_l
        end if
        ! Exit subroutine
        return
    end if
    ! Now src_q is non-zero
    ! Get spherical coordinates
    if (c(1) .eq. zero) then
        max12 = abs(c(2))
        ssq12 = one
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = one + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = one + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (c(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = one + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=zero) induced potential depends only on the
    ! leading 1-st spherical harmonic, so in adjoint way we update only the leading
    if (rho .eq. zero) then
        if (beta .eq. zero) then
            dst_l(1) = src_q * vscales_rel(1)
            dst_l(2:) = zero
        else
            dst_l(1) = beta*dst_l(1) + src_q*vscales_rel(1)
            dst_l(2:) = beta * dst_l(2:)
        end if
        return
    end if
    ! Compute the actual induced potential
    rcoef = rho / dst_r
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                ! l = 0
                if (beta .eq. zero) then
                    dst_l(1) = src_q * vscales_rel(1)
                else
                    dst_l(1) = beta*dst_l(1) + src_q*vscales_rel(1)
                end if
                return
            case (1)
                ! l = 0 and l = 1
                tmp = src_q * rcoef
                tmp2 = -tmp * vscales_rel(4) * stheta
                if (beta .eq. zero) then
                    dst_l(1) = src_q * vscales_rel(1)
                    dst_l(2) = tmp2 * sphi
                    dst_l(3) = tmp * vscales_rel(3) * ctheta
                    dst_l(4) = tmp2 * cphi
                else
                    dst_l(1) = beta*dst_l(1) + src_q*vscales_rel(1)
                    dst_l(2) = beta*dst_l(2) + tmp2*sphi
                    dst_l(3) = beta*dst_l(3) + tmp*vscales_rel(3)*ctheta
                    dst_l(4) = beta*dst_l(4) + tmp2*cphi
                end if
                return
        end select
        ! Now p>1
        ! Precompute src_q*rcoef^l
        work(1) = src_q
        do l = 1, p
            work(l+1) = rcoef * work(l)
        end do
        ! Overwrite output
        if (beta .eq. zero) then
            ! Case m = 0
            ! P_{l-2}^m which is P_0^0 now
            pl2m = one
            dst_l(1) = work(1) * vscales_rel(1)
            ! Update P_m^m for the next iteration
            pmm = -stheta
            ! P_{l-1}^m which is P_{m+1}^m now
            pl1m = ctheta
            ylm = pl1m * vscales_rel(3)
            dst_l(3) = work(2) * ylm
            ! P_l^m for l>m+1
            do l = 2, p
                plm = dble(2*l-1)*ctheta*pl1m - dble(l-1)*pl2m
                plm = plm / dble(l)
                ylm = plm * vscales_rel(l*l+l+1)
                dst_l(l*l+l+1) = work(l+1) * ylm
                ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                pl2m = pl1m
                pl1m = plm
            end do
            ! Prepare cos(m*phi) and sin(m*phi) for m=1
            cmphi = cphi
            smphi = sphi
            ! Case 0<m<p
            do m = 1, p-1
                ! P_{l-2}^m which is P_m^m now
                pl2m = pmm
                ylm = pmm * vscales_rel((m+1)**2)
                ylm = work(m+1) * ylm
                dst_l((m+1)**2) = cmphi * ylm
                dst_l(m*m+1) = smphi * ylm
                ! Temporary to reduce number of operations
                tmp1 = dble(2*m+1) * pmm
                ! Update P_m^m for the next iteration
                pmm = -stheta * tmp1
                ! P_{l-1}^m which is P_{m+1}^m now
                pl1m = ctheta * tmp1
                ylm = pl1m * vscales_rel((m+1)*(m+3))
                ylm = work(m+2) * ylm
                dst_l((m+1)*(m+3)) = cmphi * ylm
                dst_l((m+1)*(m+2)+1-m) = smphi * ylm
                ! P_l^m for l>m+1
                do l = m+2, p
                    plm = dble(2*l-1)*ctheta*pl1m - dble(l+m-1)*pl2m
                    plm = plm / dble(l-m)
                    ylm = plm * vscales_rel(l*l+l+1+m)
                    ylm = work(l+1) * ylm
                    dst_l(l*l+l+1+m) = cmphi * ylm
                    dst_l(l*l+l+1-m) = smphi * ylm
                    ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                    pl2m = pl1m
                    pl1m = plm
                end do
                ! Update cos(m*phi) and sin(m*phi) for the next iteration
                tmp1 = cmphi
                cmphi = cmphi*cphi - smphi*sphi
                smphi = tmp1*sphi + smphi*cphi
            end do
            ! Case m=p requires only to use P_m^m
            ylm = pmm * vscales_rel((p+1)**2)
            ylm = work(p+1) * ylm
            dst_l((p+1)**2) = cmphi * ylm
            dst_l(p*p+1) = smphi * ylm
        ! Update output
        else
            ! Case m = 0
            ! P_{l-2}^m which is P_0^0 now
            pl2m = one
            dst_l(1) = beta*dst_l(1) + work(1)*vscales_rel(1)
            ! Update P_m^m for the next iteration
            pmm = -stheta
            ! P_{l-1}^m which is P_{m+1}^m now
            pl1m = ctheta
            ylm = pl1m * vscales_rel(3)
            dst_l(3) = beta*dst_l(3) + work(2)*ylm
            ! P_l^m for l>m+1
            do l = 2, p
                plm = dble(2*l-1)*ctheta*pl1m - dble(l-1)*pl2m
                plm = plm / dble(l)
                ylm = plm * vscales_rel(l*l+l+1)
                dst_l(l*l+l+1) = beta*dst_l(l*l+l+1) + work(l+1)*ylm
                ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                pl2m = pl1m
                pl1m = plm
            end do
            ! Prepare cos(m*phi) and sin(m*phi) for m=1
            cmphi = cphi
            smphi = sphi
            ! Case 0<m<p
            do m = 1, p-1
                ! P_{l-2}^m which is P_m^m now
                pl2m = pmm
                ylm = pmm * vscales_rel((m+1)**2)
                ylm = work(m+1) * ylm
                dst_l((m+1)**2) = beta*dst_l((m+1)**2) + cmphi*ylm
                dst_l(m*m+1) = beta*dst_l(m*m+1) + smphi*ylm
                ! Temporary to reduce number of operations
                tmp1 = dble(2*m+1) * pmm
                ! Update P_m^m for the next iteration
                pmm = -stheta * tmp1
                ! P_{l-1}^m which is P_{m+1}^m now
                pl1m = ctheta * tmp1
                ylm = pl1m * vscales_rel((m+1)*(m+3))
                ylm = work(m+2) * ylm
                dst_l((m+1)*(m+3)) = beta*dst_l((m+1)*(m+3)) + cmphi*ylm
                dst_l((m+1)*(m+2)+1-m) = beta*dst_l((m+1)*(m+2)+1-m) + &
                    & smphi*ylm
                ! P_l^m for l>m+1
                do l = m+2, p
                    plm = dble(2*l-1)*ctheta*pl1m - dble(l+m-1)*pl2m
                    plm = plm / dble(l-m)
                    ylm = plm * vscales_rel(l*l+l+1+m)
                    ylm = work(l+1) * ylm
                    dst_l(l*l+l+1+m) = beta*dst_l(l*l+l+1+m) + cmphi*ylm
                    dst_l(l*l+l+1-m) = beta*dst_l(l*l+l+1-m) + smphi*ylm
                    ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                    pl2m = pl1m
                    pl1m = plm
                end do
                ! Update cos(m*phi) and sin(m*phi) for the next iteration
                tmp1 = cmphi
                cmphi = cmphi*cphi - smphi*sphi
                smphi = tmp1*sphi + smphi*cphi
            end do
            ! Case m=p requires only to use P_m^m
            ylm = pmm * vscales_rel((p+1)**2)
            ylm = work(p+1) * ylm
            dst_l((p+1)**2) = beta*dst_l((p+1)**2) + cmphi*ylm
            dst_l(p*p+1) = beta*dst_l(p*p+1) + smphi*ylm
        end if
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^m = ctheta^m in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef.
        rcoef = sign(rcoef, c(3))
        ! Init t and proceed with accumulation of a potential
        t = src_q
        indl = 1
        ! Overwrite output
        if (beta .eq. zero) then
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_l(indl) = t * vscales_rel(indl)
                dst_l(indl-l:indl-1) = zero
                dst_l(indl+1:indl+l) = zero
                ! Update t
                t = t * rcoef
            end do
        ! Update output
        else
            do l = 0, p
                ! Index of Y_l^0
                indl = indl + 2*l
                ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
                dst_l(indl) = beta*dst_l(indl) + t*vscales_rel(indl)
                dst_l(indl-l:indl-1) = beta * dst_l(indl-l:indl-1)
                dst_l(indl+1:indl+l) = beta * dst_l(indl+1:indl+l)
                ! Update t
                t = t * rcoef
            end do
        end if
    end if
end subroutine fmm_l2p_adj_work

!> Transform coefficients of spherical harmonics to a new cartesion system
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$ R \f$
!! is a '(p+1)**2'-by-'(p+1)**2' matrix that transforms coefficients of the
!! source spherical harmonics to the destination spherical harmonics while
!! preserving values of any function on a unit sphere defined as a linear
!! combination of the source spherical harmonics.
!!
!! Input 3-by-3 matrix `r1` must be an orthogonal matrix \f$ R_1 \f$ of a
!! transform of new cartesion coordinates \f$ (\widetilde{y}, \widetilde{z},
!! \widetilde{x}) \f$ into initial cartesian coordinates \f$ (y, z, x) \f$.
!! This is due to the following equalities:
!! \f{align}{
!!      Y_1^{-1} (\theta, \phi) &= \sqrt{\frac{3}{4\pi}} \sin \theta \sin \phi
!!      = \sqrt{\frac{3}{4\pi}} y, \\ Y_1^0 (\theta, \phi) &=
!!      \sqrt{\frac{3}{4\pi}} \cos \theta = \sqrt{\frac{3}{4\pi}} z, \\
!!      Y_1^1 (\theta, \phi) &= \sqrt{\frac{3}{4\pi}} \sin \theta \cos \phi =
!!      \sqrt{\frac{3}{4\pi}} x.
!! \f}
!! So, to find a column-vector \f$ \widetilde{c} \f$ of coefficients of
!! spherical harmonics \f$ Y_1^{-1}, Y_1^0 \f$ and \f$ Y_1^1 \f$ in a new
!! system of coordinates \f$ (\widetilde{y}, \widetilde{z}, \widetilde{x}) \f$
!! the following system needs to be solved:
!! \f[
!!      \widetilde{c}^\top \cdot \begin{bmatrix} Y_1^{-1}
!!      (\widetilde{\theta}, \widetilde{\phi}) \\ Y_1^0 (\widetilde{\theta},
!!      \widetilde{\phi}) \\ Y_1^1 (\widetilde{\theta}, \widetilde{\phi})
!!      \end{bmatrix} = c ^\top \cdot \begin{bmatrix} Y_1^{-1} (\theta, \phi)
!!      \\ Y_1^0 (\theta, \phi) \\ Y_1^1 (\theta, \phi) \end{bmatrix}.
!! \f]
!! The solution has the following obvious form:
!! \f[
!!      \widetilde{c} = R_1^\top c.
!! \f]
!!
!! Translation of spherical harmonics of all other degrees is computed
!! recursively as is described in the following source:
!!      @cite ir-realharms-1996
!!      @cite ir-realharms-1998
!!
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] r1: Transformation from new to old cartesian coordinates
!! @param[in] alpha: Scalar multipler for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of transformed spherical harmonics
subroutine fmm_sph_transform(p, r1, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1(-1:1, -1:1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(2*p+1)*(2*p+3))
    ! Call corresponding work routine
    call fmm_sph_transform_work(p, r1, alpha, src, beta, dst, work)
end subroutine fmm_sph_transform

!> Transform spherical harmonics to a new cartesion system of coordinates
!!
!! This function implements @ref fmm_sph_transform with predefined values of
!! parameters \p alpha=one and \p beta=zero.
!! 
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] r1: Transformation from new to old cartesian coordinates
!! @param[in] alpha: Scalar multipler for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of transformed spherical harmonics
!! @param[out] work: Temporary workspace of a size (2*(2*p+1)*(2*p+3))
subroutine fmm_sph_transform_work(p, r1, alpha, src, beta, dst, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1(-1:1, -1:1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(out) :: dst((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(2*p+1)*(2*p+3))
    ! Local variables
    real(dp) :: u, v, w
    integer :: l, m, n, ind
    ! Pointers for a workspace
    real(dp), pointer :: r_prev(:, :), r(:, :), scal_uvw_m(:), scal_u_n(:), &
        & scal_v_n(:), scal_w_n(:), r_swap(:, :)
    ! In case alpha is zero just scale output
    if (alpha .eq. zero) then
        ! Set output to zero if beta is also zero
        if (beta .eq. zero) then
            dst = zero
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-zero
    ! In case beta is zero output is just overwritten without being read
    if (beta .eq. zero) then
        ! Compute rotations/reflections
        ! l = 0
        dst(1) = alpha * src(1)
        if (p .eq. 0) then
            return
        end if
        ! l = 1
        dst(2) = alpha * dot_product(r1(-1:1,-1), src(2:4))
        dst(3) = alpha * dot_product(r1(-1:1,0), src(2:4))
        dst(4) = alpha * dot_product(r1(-1:1,1), src(2:4))
        if (p .eq. 1) then
            return
        end if
        ! Set pointers
        n = 2*p + 1
        l = n ** 2
        r_prev(-p:p, -p:p) => work(1:l)
        m = 2*l
        r(-p:p, -p:p) => work(l+1:m)
        l = m + n
        scal_uvw_m(-p:p) => work(m+1:l)
        m = l + n
        scal_u_n(-p:p) => work(l+1:m)
        l = m + n
        scal_v_n(-p:p) => work(m+1:l)
        m = l + n
        scal_w_n(-p:p) => work(l+1:m)
        ! l = 2
        r(2, 2) = (r1(1, 1)*r1(1, 1) - r1(1, -1)*r1(1, -1) - &
            & r1(-1, 1)*r1(-1, 1) + r1(-1, -1)*r1(-1, -1)) / two
        r(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
        r(2, 0) = sqrt3 / two * (r1(1, 0)*r1(1, 0) - r1(-1, 0)*r1(-1, 0))
        r(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
        r(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
        r(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
        r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
        r(1, 0) = sqrt3 * r1(1, 0) * r1(0, 0)
        r(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
        r(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
        r(0, 2) = sqrt3 / two * (r1(0, 1)*r1(0, 1) - r1(0, -1)*r1(0, -1))
        r(0, 1) = sqrt3 * r1(0, 1) * r1(0, 0)
        r(0, 0) = (three*r1(0, 0)*r1(0, 0)-one) / two
        r(0, -1) = sqrt3 * r1(0, -1) * r1(0, 0)
        r(0, -2) = sqrt3 * r1(0, 1) * r1(0, -1)
        r(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
        r(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
        r(-1, 0) = sqrt3 * r1(-1, 0) * r1(0, 0)
        r(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(0, -1)*r1(-1, 0)
        r(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
        r(-2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
        r(-2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
        r(-2, 0) = sqrt3 * r1(1, 0) * r1(-1, 0)
        r(-2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
        r(-2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
        do m = -2, 2
            dst(7+m) = alpha * dot_product(r(-2:2, m), src(5:9))
        end do
        ! l > 2
        do l = 3, p
            ! Swap previous and current rotation matrices
            r_swap => r_prev
            r_prev => r
            r => r_swap
            ! Prepare scalar factors
            scal_uvw_m(0) = dble(l)
            do m = 1, l-1
                u = sqrt(dble(l*l-m*m))
                scal_uvw_m(m) = u
                scal_uvw_m(-m) = u
            end do
            u = two * dble(l)
            u = sqrt(dble(u*(u-one)))
            scal_uvw_m(l) = u
            scal_uvw_m(-l) = u
            scal_u_n(0) = dble(l)
            scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
            scal_w_n(0) = zero
            do n = 1, l-2
                u = sqrt(dble(l*l-n*n))
                scal_u_n(n) = u
                scal_u_n(-n) = u
                v = dble(l+n)
                v = sqrt(v*(v-one)) / two
                scal_v_n(n) = v
                scal_v_n(-n) = v
                w = dble(l-n)
                w = -sqrt(w*(w-one)) / two
                scal_w_n(n) = w
                scal_w_n(-n) = w
            end do
            u = sqrt(dble(2*l-1))
            scal_u_n(l-1) = u
            scal_u_n(1-l) = u
            scal_u_n(l) = zero
            scal_u_n(-l) = zero
            v = sqrt(dble((2*l-1)*(2*l-2))) / two
            scal_v_n(l-1) = v
            scal_v_n(1-l) = v
            v = sqrt(dble(2*l*(2*l-1))) / two
            scal_v_n(l) = v
            scal_v_n(-l) = v
            scal_w_n(l-1) = zero
            scal_w_n(l) = zero
            scal_w_n(-l) = zero
            scal_w_n(1-l) = zero
            ind = l*l + l + 1
            ! m = l, n = l
            v = r1(1, 1)*r_prev(l-1, l-1) - r1(1, -1)*r_prev(l-1, 1-l) - &
                & r1(-1, 1)*r_prev(1-l, l-1) + r1(-1, -1)*r_prev(1-l, 1-l)
            r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind+l) = src(ind+l) * r(l, l)
            ! m = l, n = -l
            v = r1(1, 1)*r_prev(1-l, l-1) - r1(1, -1)*r_prev(1-l, 1-l) + &
                & r1(-1, 1)*r_prev(l-1, l-1) - r1(-1, -1)*r_prev(l-1, 1-l)
            r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
            ! m = l, n = l-1
            u = r1(0, 1)*r_prev(l-1, l-1) - r1(0, -1)*r_prev(l-1, 1-l)
            v = r1(1, 1)*r_prev(l-2, l-1) - r1(1, -1)*r_prev(l-2, 1-l) - &
                & r1(-1, 1)*r_prev(2-l, l-1) + r1(-1, -1)*r_prev(2-l, 1-l)
            r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
            ! m = l, n = 1-l
            u = r1(0, 1)*r_prev(1-l, l-1) - r1(0, -1)*r_prev(1-l, 1-l)
            v = r1(1, 1)*r_prev(2-l, l-1) - r1(1, -1)*r_prev(2-l, 1-l) + &
                & r1(-1, 1)*r_prev(l-2, l-1) - r1(-1, -1)*r_prev(l-2, 1-l)
            r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
            ! m = l, n = 1
            u = r1(0, 1)*r_prev(1, l-1) - r1(0, -1)*r_prev(1, 1-l)
            v = r1(1, 1)*r_prev(0, l-1) - r1(1, -1)*r_prev(0, 1-l)
            w = r1(1, 1)*r_prev(2, l-1) - r1(1, -1)*r_prev(2, 1-l) + &
                & r1(-1, 1)*r_prev(-2, l-1) - r1(-1, -1)*r_prev(-2, 1-l)
            r(1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, l) = r(1, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
            ! m = l, n = -1
            u = r1(0, 1)*r_prev(-1, l-1) - r1(0, -1)*r_prev(-1, 1-l)
            v = r1(-1, 1)*r_prev(0, l-1) - r1(-1, -1)*r_prev(0, 1-l)
            w = r1(1, 1)*r_prev(-2, l-1) - r1(1, -1)*r_prev(-2, 1-l) - &
                & r1(-1, 1)*r_prev(2, l-1) + r1(-1, -1)*r_prev(2, 1-l)
            r(-1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, l) = r(-1, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
            ! m = l, n = 0
            u = r1(0, 1)*r_prev(0, l-1) - r1(0, -1)*r_prev(0, 1-l)
            v = r1(1, 1)*r_prev(1, l-1) - r1(1, -1)*r_prev(1, 1-l) + &
                & r1(-1, 1)*r_prev(-1, l-1) - r1(-1, -1)*r_prev(-1, 1-l)
            r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
            ! m = l, n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 1)*r_prev(n, l-1) - r1(0, -1)*r_prev(n, 1-l)
                v = r1(1, 1)*r_prev(n-1, l-1) - r1(1, -1)*r_prev(n-1, 1-l) - &
                    & r1(-1, 1)*r_prev(1-n, l-1) + r1(-1, -1)*r_prev(1-n, 1-l)
                w = r1(1, 1)*r_prev(n+1, l-1) - r1(1, -1)*r_prev(n+1, 1-l) + &
                    & r1(-1, 1)*r_prev(-n-1, l-1) - r1(-1, -1)*r_prev(-n-1, 1-l)
                r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, l) = r(n, l) / scal_uvw_m(l)
                !dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
                u = r1(0, 1)*r_prev(-n, l-1) - r1(0, -1)*r_prev(-n, 1-l)
                v = r1(1, 1)*r_prev(1-n, l-1) - r1(1, -1)*r_prev(1-n, 1-l) + &
                    & r1(-1, 1)*r_prev(n-1, l-1) - r1(-1, -1)*r_prev(n-1, 1-l)
                w = r1(1, 1)*r_prev(-n-1, l-1) - r1(1, -1)*r_prev(-n-1, 1-l) - &
                    & r1(-1, 1)*r_prev(n+1, l-1) + r1(-1, -1)*r_prev(n+1, 1-l)
                r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, l) = r(-n, l) / scal_uvw_m(l)
                !dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
            end do
            dst(ind+l) = alpha * dot_product(src(ind-l:ind+l), r(-l:l, l))
            ! m = -l, n = l
            v = r1(1, 1)*r_prev(l-1, 1-l) + r1(1, -1)*r_prev(l-1, l-1) - &
                & r1(-1, 1)*r_prev(1-l, 1-l) - r1(-1, -1)*r_prev(1-l, l-1)
            r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind-l) = src(ind+l)*r(l, -l)
            ! m = -l, n = -l
            v = r1(1, 1)*r_prev(1-l, 1-l) + r1(1, -1)*r_prev(1-l, l-1) + &
                & r1(-1, 1)*r_prev(l-1, 1-l) + r1(-1, -1)*r_prev(l-1, l-1)
            r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
            ! m = -l, n = l-1
            u = r1(0, 1)*r_prev(l-1, 1-l) + r1(0, -1)*r_prev(l-1, l-1)
            v = r1(1, 1)*r_prev(l-2, 1-l) + r1(1, -1)*r_prev(l-2, l-1) - &
                & r1(-1, 1)*r_prev(2-l, 1-l) - r1(-1, -1)*r_prev(2-l, l-1)
            r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
            ! m = -l, n = 1-l
            u = r1(0, 1)*r_prev(1-l, 1-l) + r1(0, -1)*r_prev(1-l, l-1)
            v = r1(1, 1)*r_prev(2-l, 1-l) + r1(1, -1)*r_prev(2-l, l-1) + &
                & r1(-1, 1)*r_prev(l-2, 1-l) + r1(-1, -1)*r_prev(l-2, l-1)
            r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
            ! m = -l, n = 1
            u = r1(0, 1)*r_prev(1, 1-l) + r1(0, -1)*r_prev(1, l-1)
            v = r1(1, 1)*r_prev(0, 1-l) + r1(1, -1)*r_prev(0, l-1)
            w = r1(1, 1)*r_prev(2, 1-l) + r1(1, -1)*r_prev(2, l-1) + &
                & r1(-1, 1)*r_prev(-2, 1-l) + r1(-1, -1)*r_prev(-2, l-1)
            r(1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, -l) = r(1, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
            ! m = -l, n = -1
            u = r1(0, 1)*r_prev(-1, 1-l) + r1(0, -1)*r_prev(-1, l-1)
            v = r1(-1, 1)*r_prev(0, 1-l) + r1(-1, -1)*r_prev(0, l-1)
            w = r1(1, 1)*r_prev(-2, 1-l) + r1(1, -1)*r_prev(-2, l-1) - &
                & r1(-1, 1)*r_prev(2, 1-l) - r1(-1, -1)*r_prev(2, l-1)
            r(-1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
            ! m = -l, n = 0
            u = r1(0, 1)*r_prev(0, 1-l) + r1(0, -1)*r_prev(0, l-1)
            v = r1(1, 1)*r_prev(1, 1-l) + r1(1, -1)*r_prev(1, l-1) + &
                & r1(-1, 1)*r_prev(-1, 1-l) + r1(-1, -1)*r_prev(-1, l-1)
            r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
            ! m = -l, n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 1)*r_prev(n, 1-l) + r1(0, -1)*r_prev(n, l-1)
                v = r1(1, 1)*r_prev(n-1, 1-l) + r1(1, -1)*r_prev(n-1, l-1) - &
                    & r1(-1, 1)*r_prev(1-n, 1-l) - r1(-1, -1)*r_prev(1-n, l-1)
                w = r1(1, 1)*r_prev(n+1, 1-l) + r1(1, -1)*r_prev(n+1, l-1) + &
                    & r1(-1, 1)*r_prev(-n-1, 1-l) + r1(-1, -1)*r_prev(-n-1, l-1)
                r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, -l) = r(n, -l) / scal_uvw_m(l)
                !dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
                u = r1(0, 1)*r_prev(-n, 1-l) + r1(0, -1)*r_prev(-n, l-1)
                v = r1(1, 1)*r_prev(1-n, 1-l) + r1(1, -1)*r_prev(1-n, l-1) + &
                    & r1(-1, 1)*r_prev(n-1, 1-l) + r1(-1, -1)*r_prev(n-1, l-1)
                w = r1(1, 1)*r_prev(-n-1, 1-l) + r1(1, -1)*r_prev(-n-1, l-1) - &
                    & r1(-1, 1)*r_prev(n+1, 1-l) - r1(-1, -1)*r_prev(n+1, l-1)
                r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
                !dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
            end do
            dst(ind-l) = alpha * dot_product(src(ind-l:ind+l), r(-l:l, -l))
            ! Now deal with m=1-l..l-1
            do m = 1-l, l-1
                ! n = l
                v = r1(1, 0)*r_prev(l-1, m) - r1(-1, 0)*r_prev(1-l, m)
                r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
                !dst(ind+m) = src(ind+l) * r(l, m)
                ! n = -l
                v = r1(1, 0)*r_prev(1-l, m) + r1(-1, 0)*r_prev(l-1, m)
                r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
                ! n = l-1
                u = r1(0, 0) * r_prev(l-1, m)
                v = r1(1, 0)*r_prev(l-2, m) - r1(-1, 0)*r_prev(2-l, m)
                r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
                r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
                ! n = 1-l
                u = r1(0, 0) * r_prev(1-l, m)
                v = r1(1, 0)*r_prev(2-l, m) + r1(-1, 0)*r_prev(l-2, m)
                r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
                r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
                ! n = 0
                u = r1(0, 0) * r_prev(0, m)
                v = r1(1, 0)*r_prev(1, m) + r1(-1, 0)*r_prev(-1, m)
                r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
                r(0, m) = r(0, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
                ! n = 1
                u = r1(0, 0) * r_prev(1, m)
                v = r1(1, 0) * r_prev(0, m)
                w = r1(1, 0)*r_prev(2, m) + r1(-1, 0)*r_prev(-2, m)
                r(1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
                r(1, m) = r(1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
                ! n = -1
                u = r1(0, 0) * r_prev(-1, m)
                v = r1(-1, 0) * r_prev(0, m)
                w = r1(1, 0)*r_prev(-2, m) - r1(-1, 0)*r_prev(2, m)
                r(-1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
                r(-1, m) = r(-1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
                ! n = 2-l..l-2, n != -1,0,1
                do n = 2, l-2
                    u = r1(0, 0) * r_prev(n, m)
                    v = r1(1, 0)*r_prev(n-1, m) - r1(-1, 0)*r_prev(1-n, m)
                    w = r1(1, 0)*r_prev(n+1, m) + r1(-1, 0)*r_prev(-1-n, m)
                    r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                    r(n, m) = r(n, m) / scal_uvw_m(m)
                    !dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                    u = r1(0, 0) * r_prev(-n, m)
                    v = r1(1, 0)*r_prev(1-n, m) + r1(-1, 0)*r_prev(n-1, m)
                    w = r1(1, 0)*r_prev(-n-1, m) - r1(-1, 0)*r_prev(n+1, m)
                    r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                    r(-n, m) = r(-n, m) / scal_uvw_m(m)
                    !dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
                end do
                dst(ind+m) = alpha * dot_product(src(ind-l:ind+l), r(-l:l, m))
            end do
        end do
    else
        ! Compute rotations/reflections
        ! l = 0
        dst(1) = beta*dst(1) + alpha*src(1)
        if (p .eq. 0) then
            return
        end if
        ! l = 1
        dst(2) = beta*dst(2) + alpha*dot_product(r1(-1:1,-1), src(2:4))
        dst(3) = beta*dst(3) + alpha*dot_product(r1(-1:1,0), src(2:4))
        dst(4) = beta*dst(4) + alpha*dot_product(r1(-1:1,1), src(2:4))
        if (p .eq. 1) then
            return
        end if
        ! Set pointers
        n = 2*p + 1
        l = n ** 2
        r_prev(-p:p, -p:p) => work(1:l)
        m = 2*l
        r(-p:p, -p:p) => work(l+1:m)
        l = m + n
        scal_uvw_m(-p:p) => work(m+1:l)
        m = l + n
        scal_u_n(-p:p) => work(l+1:m)
        l = m + n
        scal_v_n(-p:p) => work(m+1:l)
        m = l + n
        scal_w_n(-p:p) => work(l+1:m)
        ! l = 2
        r(2, 2) = (r1(1, 1)*r1(1, 1) - r1(1, -1)*r1(1, -1) - &
            & r1(-1, 1)*r1(-1, 1) + r1(-1, -1)*r1(-1, -1)) / two
        r(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
        r(2, 0) = sqrt3 / two * (r1(1, 0)*r1(1, 0) - r1(-1, 0)*r1(-1, 0))
        r(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
        r(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
        r(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
        r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
        r(1, 0) = sqrt3 * r1(1, 0) * r1(0, 0)
        r(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
        r(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
        r(0, 2) = sqrt3 / two * (r1(0, 1)*r1(0, 1) - r1(0, -1)*r1(0, -1))
        r(0, 1) = sqrt3 * r1(0, 1) * r1(0, 0)
        r(0, 0) = (three*r1(0, 0)*r1(0, 0)-one) / two
        r(0, -1) = sqrt3 * r1(0, -1) * r1(0, 0)
        r(0, -2) = sqrt3 * r1(0, 1) * r1(0, -1)
        r(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
        r(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
        r(-1, 0) = sqrt3 * r1(-1, 0) * r1(0, 0)
        r(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(0, -1)*r1(-1, 0)
        r(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
        r(-2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
        r(-2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
        r(-2, 0) = sqrt3 * r1(1, 0) * r1(-1, 0)
        r(-2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
        r(-2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
        do m = -2, 2
            dst(7+m) = beta*dst(7+m) + alpha*dot_product(r(-2:2, m), src(5:9))
        end do
        ! l > 2
        do l = 3, p
            ! Swap previous and current rotation matrices
            r_swap => r_prev
            r_prev => r
            r => r_swap
            ! Prepare scalar factors
            scal_uvw_m(0) = dble(l)
            do m = 1, l-1
                u = sqrt(dble(l*l-m*m))
                scal_uvw_m(m) = u
                scal_uvw_m(-m) = u
            end do
            u = two * dble(l)
            u = sqrt(dble(u*(u-one)))
            scal_uvw_m(l) = u
            scal_uvw_m(-l) = u
            scal_u_n(0) = dble(l)
            scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
            scal_w_n(0) = zero
            do n = 1, l-2
                u = sqrt(dble(l*l-n*n))
                scal_u_n(n) = u
                scal_u_n(-n) = u
                v = dble(l+n)
                v = sqrt(v*(v-one)) / two
                scal_v_n(n) = v
                scal_v_n(-n) = v
                w = dble(l-n)
                w = -sqrt(w*(w-one)) / two
                scal_w_n(n) = w
                scal_w_n(-n) = w
            end do
            u = sqrt(dble(2*l-1))
            scal_u_n(l-1) = u
            scal_u_n(1-l) = u
            scal_u_n(l) = zero
            scal_u_n(-l) = zero
            v = sqrt(dble((2*l-1)*(2*l-2))) / two
            scal_v_n(l-1) = v
            scal_v_n(1-l) = v
            v = sqrt(dble(2*l*(2*l-1))) / two
            scal_v_n(l) = v
            scal_v_n(-l) = v
            scal_w_n(l-1) = zero
            scal_w_n(l) = zero
            scal_w_n(-l) = zero
            scal_w_n(1-l) = zero
            ind = l*l + l + 1
            ! m = l, n = l
            v = r1(1, 1)*r_prev(l-1, l-1) - r1(1, -1)*r_prev(l-1, 1-l) - &
                & r1(-1, 1)*r_prev(1-l, l-1) + r1(-1, -1)*r_prev(1-l, 1-l)
            r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind+l) = src(ind+l) * r(l, l)
            ! m = l, n = -l
            v = r1(1, 1)*r_prev(1-l, l-1) - r1(1, -1)*r_prev(1-l, 1-l) + &
                & r1(-1, 1)*r_prev(l-1, l-1) - r1(-1, -1)*r_prev(l-1, 1-l)
            r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
            ! m = l, n = l-1
            u = r1(0, 1)*r_prev(l-1, l-1) - r1(0, -1)*r_prev(l-1, 1-l)
            v = r1(1, 1)*r_prev(l-2, l-1) - r1(1, -1)*r_prev(l-2, 1-l) - &
                & r1(-1, 1)*r_prev(2-l, l-1) + r1(-1, -1)*r_prev(2-l, 1-l)
            r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
            ! m = l, n = 1-l
            u = r1(0, 1)*r_prev(1-l, l-1) - r1(0, -1)*r_prev(1-l, 1-l)
            v = r1(1, 1)*r_prev(2-l, l-1) - r1(1, -1)*r_prev(2-l, 1-l) + &
                & r1(-1, 1)*r_prev(l-2, l-1) - r1(-1, -1)*r_prev(l-2, 1-l)
            r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
            ! m = l, n = 1
            u = r1(0, 1)*r_prev(1, l-1) - r1(0, -1)*r_prev(1, 1-l)
            v = r1(1, 1)*r_prev(0, l-1) - r1(1, -1)*r_prev(0, 1-l)
            w = r1(1, 1)*r_prev(2, l-1) - r1(1, -1)*r_prev(2, 1-l) + &
                & r1(-1, 1)*r_prev(-2, l-1) - r1(-1, -1)*r_prev(-2, 1-l)
            r(1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, l) = r(1, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
            ! m = l, n = -1
            u = r1(0, 1)*r_prev(-1, l-1) - r1(0, -1)*r_prev(-1, 1-l)
            v = r1(-1, 1)*r_prev(0, l-1) - r1(-1, -1)*r_prev(0, 1-l)
            w = r1(1, 1)*r_prev(-2, l-1) - r1(1, -1)*r_prev(-2, 1-l) - &
                & r1(-1, 1)*r_prev(2, l-1) + r1(-1, -1)*r_prev(2, 1-l)
            r(-1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, l) = r(-1, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
            ! m = l, n = 0
            u = r1(0, 1)*r_prev(0, l-1) - r1(0, -1)*r_prev(0, 1-l)
            v = r1(1, 1)*r_prev(1, l-1) - r1(1, -1)*r_prev(1, 1-l) + &
                & r1(-1, 1)*r_prev(-1, l-1) - r1(-1, -1)*r_prev(-1, 1-l)
            r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
            ! m = l, n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 1)*r_prev(n, l-1) - r1(0, -1)*r_prev(n, 1-l)
                v = r1(1, 1)*r_prev(n-1, l-1) - r1(1, -1)*r_prev(n-1, 1-l) - &
                    & r1(-1, 1)*r_prev(1-n, l-1) + r1(-1, -1)*r_prev(1-n, 1-l)
                w = r1(1, 1)*r_prev(n+1, l-1) - r1(1, -1)*r_prev(n+1, 1-l) + &
                    & r1(-1, 1)*r_prev(-n-1, l-1) - r1(-1, -1)*r_prev(-n-1, 1-l)
                r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, l) = r(n, l) / scal_uvw_m(l)
                !dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
                u = r1(0, 1)*r_prev(-n, l-1) - r1(0, -1)*r_prev(-n, 1-l)
                v = r1(1, 1)*r_prev(1-n, l-1) - r1(1, -1)*r_prev(1-n, 1-l) + &
                    & r1(-1, 1)*r_prev(n-1, l-1) - r1(-1, -1)*r_prev(n-1, 1-l)
                w = r1(1, 1)*r_prev(-n-1, l-1) - r1(1, -1)*r_prev(-n-1, 1-l) - &
                    & r1(-1, 1)*r_prev(n+1, l-1) + r1(-1, -1)*r_prev(n+1, 1-l)
                r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, l) = r(-n, l) / scal_uvw_m(l)
                !dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
            end do
            dst(ind+l) = beta*dst(ind+l) + &
                & alpha*dot_product(src(ind-l:ind+l), r(-l:l, l))
            ! m = -l, n = l
            v = r1(1, 1)*r_prev(l-1, 1-l) + r1(1, -1)*r_prev(l-1, l-1) - &
                & r1(-1, 1)*r_prev(1-l, 1-l) - r1(-1, -1)*r_prev(1-l, l-1)
            r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind-l) = src(ind+l)*r(l, -l)
            ! m = -l, n = -l
            v = r1(1, 1)*r_prev(1-l, 1-l) + r1(1, -1)*r_prev(1-l, l-1) + &
                & r1(-1, 1)*r_prev(l-1, 1-l) + r1(-1, -1)*r_prev(l-1, l-1)
            r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
            ! m = -l, n = l-1
            u = r1(0, 1)*r_prev(l-1, 1-l) + r1(0, -1)*r_prev(l-1, l-1)
            v = r1(1, 1)*r_prev(l-2, 1-l) + r1(1, -1)*r_prev(l-2, l-1) - &
                & r1(-1, 1)*r_prev(2-l, 1-l) - r1(-1, -1)*r_prev(2-l, l-1)
            r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
            ! m = -l, n = 1-l
            u = r1(0, 1)*r_prev(1-l, 1-l) + r1(0, -1)*r_prev(1-l, l-1)
            v = r1(1, 1)*r_prev(2-l, 1-l) + r1(1, -1)*r_prev(2-l, l-1) + &
                & r1(-1, 1)*r_prev(l-2, 1-l) + r1(-1, -1)*r_prev(l-2, l-1)
            r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
            ! m = -l, n = 1
            u = r1(0, 1)*r_prev(1, 1-l) + r1(0, -1)*r_prev(1, l-1)
            v = r1(1, 1)*r_prev(0, 1-l) + r1(1, -1)*r_prev(0, l-1)
            w = r1(1, 1)*r_prev(2, 1-l) + r1(1, -1)*r_prev(2, l-1) + &
                & r1(-1, 1)*r_prev(-2, 1-l) + r1(-1, -1)*r_prev(-2, l-1)
            r(1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, -l) = r(1, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
            ! m = -l, n = -1
            u = r1(0, 1)*r_prev(-1, 1-l) + r1(0, -1)*r_prev(-1, l-1)
            v = r1(-1, 1)*r_prev(0, 1-l) + r1(-1, -1)*r_prev(0, l-1)
            w = r1(1, 1)*r_prev(-2, 1-l) + r1(1, -1)*r_prev(-2, l-1) - &
                & r1(-1, 1)*r_prev(2, 1-l) - r1(-1, -1)*r_prev(2, l-1)
            r(-1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
            ! m = -l, n = 0
            u = r1(0, 1)*r_prev(0, 1-l) + r1(0, -1)*r_prev(0, l-1)
            v = r1(1, 1)*r_prev(1, 1-l) + r1(1, -1)*r_prev(1, l-1) + &
                & r1(-1, 1)*r_prev(-1, 1-l) + r1(-1, -1)*r_prev(-1, l-1)
            r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
            ! m = -l, n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 1)*r_prev(n, 1-l) + r1(0, -1)*r_prev(n, l-1)
                v = r1(1, 1)*r_prev(n-1, 1-l) + r1(1, -1)*r_prev(n-1, l-1) - &
                    & r1(-1, 1)*r_prev(1-n, 1-l) - r1(-1, -1)*r_prev(1-n, l-1)
                w = r1(1, 1)*r_prev(n+1, 1-l) + r1(1, -1)*r_prev(n+1, l-1) + &
                    & r1(-1, 1)*r_prev(-n-1, 1-l) + r1(-1, -1)*r_prev(-n-1, l-1)
                r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, -l) = r(n, -l) / scal_uvw_m(l)
                !dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
                u = r1(0, 1)*r_prev(-n, 1-l) + r1(0, -1)*r_prev(-n, l-1)
                v = r1(1, 1)*r_prev(1-n, 1-l) + r1(1, -1)*r_prev(1-n, l-1) + &
                    & r1(-1, 1)*r_prev(n-1, 1-l) + r1(-1, -1)*r_prev(n-1, l-1)
                w = r1(1, 1)*r_prev(-n-1, 1-l) + r1(1, -1)*r_prev(-n-1, l-1) - &
                    & r1(-1, 1)*r_prev(n+1, 1-l) - r1(-1, -1)*r_prev(n+1, l-1)
                r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
                !dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
            end do
            dst(ind-l) = beta*dst(ind-l) + &
                & alpha*dot_product(src(ind-l:ind+l), r(-l:l, -l))
            ! Now deal with m=1-l..l-1
            do m = 1-l, l-1
                ! n = l
                v = r1(1, 0)*r_prev(l-1, m) - r1(-1, 0)*r_prev(1-l, m)
                r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
                !dst(ind+m) = src(ind+l) * r(l, m)
                ! n = -l
                v = r1(1, 0)*r_prev(1-l, m) + r1(-1, 0)*r_prev(l-1, m)
                r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
                ! n = l-1
                u = r1(0, 0) * r_prev(l-1, m)
                v = r1(1, 0)*r_prev(l-2, m) - r1(-1, 0)*r_prev(2-l, m)
                r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
                r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
                ! n = 1-l
                u = r1(0, 0) * r_prev(1-l, m)
                v = r1(1, 0)*r_prev(2-l, m) + r1(-1, 0)*r_prev(l-2, m)
                r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
                r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
                ! n = 0
                u = r1(0, 0) * r_prev(0, m)
                v = r1(1, 0)*r_prev(1, m) + r1(-1, 0)*r_prev(-1, m)
                r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
                r(0, m) = r(0, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
                ! n = 1
                u = r1(0, 0) * r_prev(1, m)
                v = r1(1, 0) * r_prev(0, m)
                w = r1(1, 0)*r_prev(2, m) + r1(-1, 0)*r_prev(-2, m)
                r(1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
                r(1, m) = r(1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
                ! n = -1
                u = r1(0, 0) * r_prev(-1, m)
                v = r1(-1, 0) * r_prev(0, m)
                w = r1(1, 0)*r_prev(-2, m) - r1(-1, 0)*r_prev(2, m)
                r(-1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
                r(-1, m) = r(-1, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
                ! n = 2-l..l-2, n != -1,0,1
                do n = 2, l-2
                    u = r1(0, 0) * r_prev(n, m)
                    v = r1(1, 0)*r_prev(n-1, m) - r1(-1, 0)*r_prev(1-n, m)
                    w = r1(1, 0)*r_prev(n+1, m) + r1(-1, 0)*r_prev(-1-n, m)
                    r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                    r(n, m) = r(n, m) / scal_uvw_m(m)
                    !dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                    u = r1(0, 0) * r_prev(-n, m)
                    v = r1(1, 0)*r_prev(1-n, m) + r1(-1, 0)*r_prev(n-1, m)
                    w = r1(1, 0)*r_prev(-n-1, m) - r1(-1, 0)*r_prev(n+1, m)
                    r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                    r(-n, m) = r(-n, m) / scal_uvw_m(m)
                    !dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
                end do
                dst(ind+m) = beta*dst(ind+m) + &
                    & alpha*dot_product(src(ind-l:ind+l), r(-l:l, m))
            end do
        end do
    end if
end subroutine fmm_sph_transform_work

!> Rotate spherical harmonics around OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Call corresponding work routine
    call fmm_sph_rotate_oz_work(p, vcos, vsin, alpha, src, beta, dst)
end subroutine fmm_sph_rotate_oz

!> Rotate spherical harmonics around OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz_work(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Local variables
    integer :: l, m, ind
    real(dp) :: v1, v2, v3, v4
    ! In case alpha is zero just scale output
    if (alpha .eq. zero) then
        ! Set output to zero if beta is also zero
        if (beta .eq. zero) then
            dst = zero
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-zero
    ! In case beta is zero output is just overwritten without being read
    if (beta .eq. zero) then
        ! l = 0
        dst(1) = alpha*src(1)
        ! l > 0
        !!GCC$ unroll 4
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = alpha*src(ind)
            ! m != 0
            !!GCC$ unroll 4
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = alpha * (v1*v3-v2*v4)
                ! m < 0
                dst(ind-m) = alpha * (v1*v4+v2*v3)
            end do
        end do
    else
        ! l = 0
        dst(1) = beta*dst(1) + alpha*src(1)
        ! l > 0
        !!GCC$ unroll 4
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = beta*dst(ind) + alpha*src(ind)
            ! m != 0
            !!GCC$ unroll 4
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = beta*dst(ind+m) + alpha*(v1*v3-v2*v4)
                ! m < 0
                dst(ind-m) = beta*dst(ind-m) + alpha*(v1*v4+v2*v3)
            end do
        end do
    end if
end subroutine fmm_sph_rotate_oz_work

!> Rotate spherical harmonics around OZ axis in an opposite direction
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz_adj(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Call corresponding work routine
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src, beta, dst)
end subroutine fmm_sph_rotate_oz_adj

!> Rotate spherical harmonics around OZ axis in an opposite direction
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Local variables
    integer :: l, m, ind
    real(dp) :: v1, v2, v3, v4
    ! In case alpha is zero just scale output
    if (alpha .eq. zero) then
        ! Set output to zero if beta is also zero
        if (beta .eq. zero) then
            dst = zero
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-zero
    ! In case beta is zero output is just overwritten without being read
    if (beta .eq. zero) then
        ! l = 0
        dst(1) = alpha*src(1)
        ! l > 0
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = alpha*src(ind)
            ! m != 0
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = alpha * (v1*v3+v2*v4)
                ! m < 0
                dst(ind-m) = alpha * (v2*v3-v1*v4)
            end do
        end do
    else
        ! l = 0
        dst(1) = beta*dst(1) + alpha*src(1)
        ! l > 0
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = beta*dst(ind) + alpha*src(ind)
            ! m != 0
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = beta*dst(ind+m) + alpha*(v1*v3+v2*v4)
                ! m < 0
                dst(ind-m) = beta*dst(ind-m) + alpha*(v2*v3-v1*v4)
            end do
        end do
    end if
end subroutine fmm_sph_rotate_oz_adj_work

!> Transform spherical harmonics in the OXZ plane
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of transformation of coordinates in the OXZ plane (y
!! coordinate remains the same) presented by 2-by-2 matrix `r1xz`.
!!
!! Based on @ref fmm_sph_transform
!! by assuming `r1(-1, 0) = r1(-1, 1) = r1(0, -1) = r1(1, -1) = 0` and
!! `r1(-1, -1) = 1`, which corresponds to the following transformation matrix:
!! \f[
!!      R_1 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & a & b \\ 0 & c & d
!!      \end{bmatrix},
!! \f]
!! where unkown elements represent input `r1xz` 2x2 array:
!! \f[
!!      R_1^{xz} = \begin{bmatrix} a & b \\ c & d \end{bmatrix}.
!! \f]
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] r1xz: Transformation from new to old coordinates in the OXZ plane
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of transformed spherical harmonics
subroutine fmm_sph_rotate_oxz(p, ctheta, stheta, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: ctheta, stheta, alpha, src((p+1)**2), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)**2)
    ! Temporary workspace
    real(dp) :: work(2*(2*p+1)*(2*p+3)+p)
    ! Call corresponding work routine
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, alpha, src, beta, dst, &
        & work)
end subroutine fmm_sph_rotate_oxz

!> Transform spherical harmonics in the OXZ plane
!!
!! This function implements @ref fmm_sph_rotate_oxz with predefined values
!! of parameters \p alpha=one and \p beta=zero.
!! 
!! @param[in] p: maximum order of spherical harmonics
!! @param[in] r1xz: 2D transformation matrix in the OXZ plane
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[out] dst: coefficients of rotated spherical harmonics
!! @param[out] work: Temporary workspace of a size (2*(2*p+1)*(2*p+3))
subroutine fmm_sph_rotate_oxz_work(p, ctheta, stheta, alpha, src, beta, dst, &
        & work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: ctheta, stheta, alpha, src((p+1)**2), beta
    ! Output
    real(dp), intent(out) :: dst((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(4*p*p+13*p+4)
    ! Local variables
    real(dp) :: u, v, w, fl, fl2, tmp1, tmp2, vu(2), vv(2), vw(2), &
        ctheta2, stheta2, cstheta
    integer :: l, m, n, ind, k
    ! Pointers for a workspace
    real(dp), pointer :: r(:, :, :), r_prev(:, :, :), scal_uvw_m(:), &
        & scal_u_n(:), scal_v_n(:), scal_w_n(:), r_swap(:, :, :), vsqr(:)
    !! Spherical harmonics Y_l^m with negative m transform into harmonics Y_l^m
    !! with the same l and negative m, while harmonics Y_l^m with non-negative
    !! m transform into harmonics Y_l^m with the same l and non-negative m.
    !! Transformations for non-negative m will be stored in r(1, :, :) and for
    !! negative m will be in r(2, :, :)
    ! In case alpha is zero just scale output
    if (alpha .eq. zero) then
        ! Set output to zero if beta is also zero
        if (beta .eq. zero) then
            dst = zero
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-zero
    ! In case beta is zero output is just overwritten without being read
    if (beta .eq. zero) then
        ! Compute rotations/reflections
        ! l = 0
        dst(1) = alpha * src(1)
        if (p .eq. 0) then
            return
        end if
        ! l = 1
        dst(2) = alpha * src(2)
        dst(3) = alpha * (src(3)*ctheta - src(4)*stheta)
        dst(4) = alpha * (src(3)*stheta + src(4)*ctheta)
        if (p .eq. 1) then
            return
        end if
        ! Set pointers
        l = 2 * (p+1) * (p+1)
        r(1:2, 0:p, 0:p) => work(1:l)
        m = 2 * l
        r_prev(1:2, 0:p, 0:p) => work(l+1:m)
        l = m + p + 1
        scal_uvw_m(0:p) => work(m+1:l)
        m = l + p
        scal_u_n(0:p-1) => work(l+1:m)
        l = m + p + 1
        scal_v_n(0:p) => work(m+1:l)
        m = l + p - 2
        scal_w_n(1:p-2) => work(l+1:m)
        l = m + p
        vsqr(1:p) => work(m+1:l)
        ! l = 2, m >= 0
        ctheta2 = ctheta * ctheta
        cstheta = ctheta * stheta
        stheta2 = stheta * stheta
        r(1, 2, 2) = (ctheta2 + one) / two
        r(1, 1, 2) = cstheta
        r(1, 0, 2) = sqrt3 / two * stheta2
        dst(9) = alpha * (src(9)*r(1, 2, 2) + src(8)*r(1, 1, 2) + &
            & src(7)*r(1, 0, 2))
        r(1, 2, 1) = -cstheta
        r(1, 1, 1) = ctheta2 - stheta2
        r(1, 0, 1) = sqrt3 * cstheta
        dst(8) = alpha * (src(9)*r(1, 2, 1) + src(8)*r(1, 1, 1) + &
            & src(7)*r(1, 0, 1))
        r(1, 2, 0) = sqrt3 / two * stheta2
        r(1, 1, 0) = -sqrt3 * cstheta
        r(1, 0, 0) = (three*ctheta2-one) / two
        dst(7) = alpha * (src(9)*r(1, 2, 0) + src(8)*r(1, 1, 0) + &
            & src(7)*r(1, 0, 0))
        ! l = 2,  m < 0
        r(2, 1, 1) = ctheta
        r(2, 2, 1) = -stheta
        dst(6) = alpha * (src(6)*r(2, 1, 1) + src(5)*r(2, 2, 1))
        r(2, 1, 2) = stheta
        r(2, 2, 2) = ctheta
        dst(5) = alpha * (src(6)*r(2, 1, 2) + src(5)*r(2, 2, 2))
        ! l > 2
        vsqr(1) = one
        vsqr(2) = four
        do l = 3, p
            ! Swap previous and current rotation matrices
            r_swap => r_prev
            r_prev => r
            r => r_swap
            ! Prepare scalar factors
            fl = dble(l)
            fl2 = fl * fl
            vsqr(l) = fl2
            scal_uvw_m(0) = one / fl
            !!GCC$ unroll 4
            do m = 1, l-1
                u = sqrt(fl2 - vsqr(m))
                scal_uvw_m(m) = one / u
            end do
            u = two * dble(l)
            u = sqrt(dble(u*(u-one)))
            scal_uvw_m(l) = one / u
            scal_u_n(0) = dble(l)
            scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
            !!GCC$ unroll 4
            do n = 1, l-2
                u = sqrt(fl2-vsqr(n))
                scal_u_n(n) = u
            end do
            !!GCC$ unroll 4
            do n = 1, l-2
                v = dble(l+n)
                v = sqrt(v*v-v) / two
                scal_v_n(n) = v
                w = dble(l-n)
                w = -sqrt(w*w-w) / two
                scal_w_n(n) = w
            end do
            u = sqrt(dble(2*l-1))
            scal_u_n(l-1) = u
            v = sqrt(dble((2*l-1)*(2*l-2))) / two
            scal_v_n(l-1) = v
            v = sqrt(dble(2*l*(2*l-1))) / two
            scal_v_n(l) = v
            ind = l*l + l + 1
            ! m = l, n = l and m = -l, n = - l
            vv = ctheta*r_prev(:, l-1, l-1) + r_prev(2:1:-1, l-1, l-1)
            r(:, l, l) = vv * scal_v_n(l) * scal_uvw_m(l)
            tmp1 = src(ind+l) * r(1, l, l)
            tmp2 = src(ind-l) * r(2, l, l)
            ! m = l, n = l-1 and m = -l, n = 1-l
            vu = stheta * r_prev(:, l-1, l-1)
            vv = ctheta*r_prev(:, l-2, l-1) + r_prev(2:1:-1, l-2, l-1)
            r(:, l-1, l) = (vu*scal_u_n(l-1)+vv*scal_v_n(l-1)) * scal_uvw_m(l)
            tmp1 = tmp1 + src(ind+l-1)*r(1, l-1, l)
            tmp2 = tmp2 + src(ind-l+1)*r(2, l-1, l)
            ! m = l, n = 1 and m = -l, n = -1
            vu = stheta * r_prev(:, 1, l-1)
            vv(1) = ctheta * r_prev(1, 0, l-1)
            vv(2) = r_prev(1, 0, l-1)
            vw = ctheta*r_prev(:, 2, l-1) - r_prev(2:1:-1, 2, l-1)
            r(:, 1, l) = vu*scal_u_n(1) + vw*scal_w_n(1) + sqrt2*scal_v_n(1)*vv
            r(:, 1, l) = r(:, 1, l) * scal_uvw_m(l)
            tmp1 = tmp1 + src(ind+1)*r(1, 1, l)
            tmp2 = tmp2 + src(ind-1)*r(2, 1, l)
            ! m = l, n = 0
            u = stheta * r_prev(1, 0, l-1)
            v = ctheta*r_prev(1, 1, l-1) - r_prev(2, 1, l-1)
            r(1, 0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) * scal_uvw_m(l)
            tmp1 = tmp1 + src(ind)*r(1, 0, l)
            ! m = l, n = 2..l-2 and m = -l, n = 2-l..-2
            !!GCC$ unroll 4
            do n = 2, l-2
                vu = stheta * r_prev(:, n, l-1)
                vv = ctheta*r_prev(:, n-1, l-1) + r_prev(2:1:-1, n-1, l-1)
                vw = ctheta*r_prev(:, n+1, l-1) - r_prev(2:1:-1, n+1, l-1)
                vu = vu*scal_u_n(n) + vv*scal_v_n(n) + vw*scal_w_n(n)
                r(:, n, l) = vu * scal_uvw_m(l)
                tmp1 = tmp1 + src(ind+n)*r(1, n, l)
                tmp2 = tmp2 + src(ind-n)*r(2, n, l)
            end do
            dst(ind+l) = alpha * tmp1
            dst(ind-l) = alpha * tmp2
            ! Now deal with m = 0
            ! n = l and n = -l
            v = -stheta * r_prev(1, l-1, 0)
            u = scal_v_n(l) * scal_uvw_m(0)
            r(1, l, 0) = v * u
            tmp1 = src(ind+l) * r(1, l, 0)
            ! n = l-1
            u = ctheta * r_prev(1, l-1, 0)
            v = -stheta * r_prev(1, l-2, 0)
            w = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(1, l-1, 0) = w * scal_uvw_m(0)
            tmp1 = tmp1 + src(ind+l-1)*r(1, l-1, 0)
            ! n = 0
            u = ctheta * r_prev(1, 0, 0)
            v = -stheta * r_prev(1, 1, 0)
            w = u*scal_u_n(0) + v*scal_v_n(0)
            r(1, 0, 0) = w * scal_uvw_m(0)
            tmp1 = tmp1 + src(ind)*r(1, 0, 0)
            ! n = 1
            v = sqrt2*scal_v_n(1)*r_prev(1, 0, 0) + &
                & scal_w_n(1)*r_prev(1, 2, 0)
            u = ctheta * r_prev(1, 1, 0)
            w = scal_u_n(1)*u - stheta*v
            r(1, 1, 0) = w * scal_uvw_m(0)
            tmp1 = tmp1 + src(ind+1)*r(1, 1, 0)
            ! n = 2..l-2
            !!GCC$ unroll 4
            do n = 2, l-2
                v = scal_v_n(n)*r_prev(1, n-1, 0) + &
                    & scal_w_n(n)*r_prev(1, n+1, 0)
                u = ctheta * r_prev(1, n, 0)
                w = scal_u_n(n)*u - stheta*v
                r(1, n, 0) = w * scal_uvw_m(0)
                tmp1 = tmp1 + src(ind+n)*r(1, n, 0)
            end do
            dst(ind) = alpha * tmp1
            ! Now deal with m=1..l-1 and m=1-l..-1
            !!GCC$ unroll 4
            do m = 1, l-1
                ! n = l and n = -l
                vv = -stheta * r_prev(:, l-1, m)
                u = scal_v_n(l) * scal_uvw_m(m)
                r(:, l, m) = vv * u
                tmp1 = src(ind+l) * r(1, l, m)
                tmp2 = src(ind-l) * r(2, l, m)
                ! n = l-1 and n = 1-l
                vu = ctheta * r_prev(:, l-1, m)
                vv = -stheta * r_prev(:, l-2, m)
                vw = vu*scal_u_n(l-1) + vv*scal_v_n(l-1)
                r(:, l-1, m) = vw * scal_uvw_m(m)
                tmp1 = tmp1 + src(ind+l-1)*r(1, l-1, m)
                tmp2 = tmp2 + src(ind-l+1)*r(2, l-1, m)
                ! n = 0
                u = ctheta * r_prev(1, 0, m)
                v = -stheta * r_prev(1, 1, m)
                w = u*scal_u_n(0) + v*scal_v_n(0)
                r(1, 0, m) = w * scal_uvw_m(m)
                tmp1 = tmp1 + src(ind)*r(1, 0, m)
                ! n = 1
                v = sqrt2*scal_v_n(1)*r_prev(1, 0, m) + &
                    & scal_w_n(1)*r_prev(1, 2, m)
                u = ctheta * r_prev(1, 1, m)
                w = scal_u_n(1)*u - stheta*v
                r(1, 1, m) = w * scal_uvw_m(m)
                tmp1 = tmp1 + src(ind+1)*r(1, 1, m)
                ! n = -1
                u = ctheta * r_prev(2, 1, m)
                w = -stheta * r_prev(2, 2, m)
                v = u*scal_u_n(1) + w*scal_w_n(1)
                r(2, 1, m) = v * scal_uvw_m(m)
                tmp2 = tmp2 + src(ind-1)*r(2, 1, m)
                ! n = 2..l-2 and n = 2-l..-2
                !!GCC$ unroll 4
                do n = 2, l-2
                    vv = scal_v_n(n)*r_prev(:, n-1, m) + &
                        & scal_w_n(n)*r_prev(:, n+1, m)
                    vu = ctheta * r_prev(:, n, m)
                    vw = scal_u_n(n)*vu - stheta*vv
                    r(:, n, m) = vw * scal_uvw_m(m)
                    tmp1 = tmp1 + src(ind+n)*r(1, n, m)
                    tmp2 = tmp2 + src(ind-n)*r(2, n, m)
                end do
                dst(ind+m) = alpha * tmp1
                dst(ind-m) = alpha * tmp2
            end do
        end do
    else
        stop "Not Implemented"
    end if
end subroutine fmm_sph_rotate_oxz_work

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_ztranslate(z, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(p+1))
    ! Call corresponding work routine
    call fmm_m2m_ztranslate_work(z, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
end subroutine fmm_m2m_ztranslate

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_m2m_ztranslate_work(z, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, tmp3, res1, res2, pow_r1
    integer :: j, k, n, indj, indjn, indjk1, indjk2
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r2(:)
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_m = zero
        else
            dst_m = beta * dst_m
        end if
        return
    end if
    ! Now alpha is non-zero
    ! If harmonics have different centers
    if (z .ne. 0) then
        ! Prepare pointers
        pow_r2(1:p+1) => work(1:p+1)
        ! Get ratios r1 and r2
        r1 = src_r / dst_r
        r2 = z / dst_r
        ! Get powers of ratio r2/r1
        r2 = r2 / r1
        pow_r1 = r1
        pow_r2(1) = one
        do j = 2, p+1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Do actual M2M
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            do j = 0, p
                ! Offset for dst_m
                indj = j*j + j + 1
                ! k = 0
                tmp1 = alpha * pow_r1 * vscales(indj)
                pow_r1 = pow_r1 * r1
                tmp2 = tmp1
                res1 = zero
                ! Offset for vcnk
                indjk1 = j*(j+1)/2 + 1
                do n = 0, j
                    ! Offset for src_m
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp3 = pow_r2(n+1) / &
                        & vscales(indjn) * vcnk(indjk1+n)**2
                    res1 = res1 + tmp3*src_m(indjn)
                end do
                dst_m(indj) = tmp2 * res1
                ! k != 0
                do k = 1, j
                    tmp2 = tmp1
                    res1 = zero
                    res2 = zero
                    ! Offsets for vcnk
                    indjk1 = (j-k)*(j-k+1)/2 + 1
                    indjk2 = (j+k)*(j+k+1)/2 + 1
                    do n = 0, j-k
                        ! Offset for src_m
                        indjn = (j-n)**2 + (j-n) + 1
                        tmp3 = pow_r2(n+1) / &
                            & vscales(indjn) * vcnk(indjk1+n) * &
                            & vcnk(indjk2+n)
                        res1 = res1 + tmp3*src_m(indjn+k)
                        res2 = res2 + tmp3*src_m(indjn-k)
                    end do
                    dst_m(indj+k) = tmp2 * res1
                    dst_m(indj-k) = tmp2 * res2
                end do
            end do
        ! Update output if beta is non-zero
        else
            do j = 0, p
                ! Offset for dst_m
                indj = j*j + j + 1
                ! k = 0
                tmp1 = alpha * pow_r1 * vscales(indj)
                pow_r1 = pow_r1 * r1
                tmp2 = tmp1
                res1 = zero
                ! Offset for vcnk
                indjk1 = j * (j+1) /2 + 1
                do n = 0, j
                    ! Offset for src_m
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp3 = pow_r2(n+1) / &
                        & vscales(indjn) * vcnk(indjk1+n)**2
                    res1 = res1 + tmp3*src_m(indjn)
                end do
                dst_m(indj) = beta*dst_m(indj) + tmp2*res1
                ! k != 0
                do k = 1, j
                    tmp2 = tmp1
                    res1 = zero
                    res2 = zero
                    ! Offsets for vcnk
                    indjk1 = (j-k)*(j-k+1)/2 + 1
                    indjk2 = (j+k)*(j+k+1)/2 + 1
                    do n = 0, j-k
                        ! Offset for src_m
                        indjn = (j-n)**2 + (j-n) + 1
                        tmp3 = pow_r2(n+1) / &
                            & vscales(indjn) * vcnk(indjk1+n) * &
                            & vcnk(indjk2+n)
                        res1 = res1 + tmp3*src_m(indjn+k)
                        res2 = res2 + tmp3*src_m(indjn-k)
                    end do
                    dst_m(indj+k) = beta*dst_m(indj+k) + tmp2*res1
                    dst_m(indj-k) = beta*dst_m(indj-k) + tmp2*res2
                end do
            end do
        end if
    ! If harmonics are located at the same point
    else
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_m(k) = tmp1 * src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        ! Update output if beta is non-zero
        else
            r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_m(k) = beta*dst_m(k) + tmp1*src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        end if
    end if
end subroutine fmm_m2m_ztranslate_work

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_bessel_ztranslate(z, src_sk, dst_sk, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_sk(p+1), dst_sk(p+1), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*p+1)
    complex(dp) :: work_complex(2*p+1)
    ! Call corresponding work routine
    call fmm_m2m_bessel_ztranslate_work(z, src_sk, dst_sk, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m, work, work_complex)
end subroutine fmm_m2m_bessel_ztranslate

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from old to new centers of harmonics. Standard
!!      FMM M2M operation requires z to be OZ coordinate from new to old.
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_m2m_bessel_ztranslate_work(z, src_sk, dst_sk, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m, work, work_complex)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_sk(p+1), dst_sk(p+1), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*p+1)
    complex(dp), intent(out) :: work_complex(2*p+1)
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, tmp3, res1, res2, pow_r1
    integer :: j, k, n, indj, indn, l, NZ, ierr
    ! Pointers for temporary values of powers
    real(dp) :: fact(2*p+1)
    complex(dp) :: complex_argument
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_m = zero
        else
            dst_m = beta * dst_m
        end if
        return
    end if
    ! Get factorials
    fact(1) = one
    do j = 1, 2*p
        fact(j+1) = dble(j) * fact(j)
    end do
    ! Now alpha is non-zero
    ! If harmonics have different centers
    if (z .ne. 0) then
        complex_argument = z
        call cbesi(complex_argument, pt5, 1, 1, work_complex(1), NZ, ierr)
        work(1) = sqrt(pi/(2*z)) * real(work_complex(1))
        if (p .gt. 0) then
            call cbesi(complex_argument, 1.5d0, 1, 2*p, &
                & work_complex(2:2*p+1), NZ, ierr)
            work(2:2*p+1) = sqrt(pi/(2*z)) * real(work_complex(2:2*p+1))
        end if
        ! Do actual M2M
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            dst_m = zero
            do j = 0, p
                ! Offset for dst_m
                indj = j*j + j + 1
                ! k = 0
                k = 0
                res1 = zero
                do n = 0, p
                    ! Offset for src_m
                    indn = n*n + n + 1
                    tmp1 = vscales(indn) * &
                        & dble(2*j+1) * fact(j+1) * fact(n+1) / &
                        & src_sk(n+1) / vscales(indj) * dst_sk(j+1)
                    tmp3 = zero
                    do l = 0, min(j, n)
                        tmp2 = (two**(-l)) / &
                            & fact(l+1)*fact(2*l+1)/ &
                            & fact(l+1)/fact(l+1)/fact(j-l+1)/fact(n-l+1)* &
                            & (work(j+n-l+1) / (z**l))
                        tmp3 = tmp3 + tmp2
                    end do
                    tmp3 = tmp3 * tmp1
                    res1 = res1 + tmp3*src_m(indn)
                end do
                dst_m(indj) = res1
                ! k != 0
                do k = 1, j
                    res1 = zero
                    res2 = zero
                    do n = k, p
                        ! Offset for src_m
                        indn = n*n + n + 1
                        tmp1 = vscales(indn+k) * &
                            & dble(2*j+1) * fact(j-k+1) * fact(n+k+1) / &
                            & src_sk(n+1) / vscales(indj+k) * dst_sk(j+1)
                        tmp3 = zero
                        do l = k, min(j, n)
                            tmp2 = (two**(-l)) / &
                                & fact(l+k+1)*fact(2*l+1)/ &
                                & fact(l+1)/fact(l-k+1)/fact(j-l+1)/fact(n-l+1)* &
                                & (work(j+n-l+1) / (z**l))
                            tmp3 = tmp3 + tmp2
                        end do
                        tmp3 = tmp3 * tmp1
                        res1 = res1 + tmp3*src_m(indn+k)
                        res2 = res2 + tmp3*src_m(indn-k)
                    end do
                    dst_m(indj+k) = res1
                    dst_m(indj-k) = res2
                end do
            end do
        ! Update output if beta is non-zero
        else
            stop 111
        end if
    ! If harmonics are located at the same point
    else
        stop 222
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            !r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_m(k) = tmp1 * src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        ! Update output if beta is non-zero
        else
            !r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_m(k) = beta*dst_m(k) + tmp1*src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        end if
    end if
end subroutine fmm_m2m_bessel_ztranslate_work

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from old to new centers of harmonics. Standard
!!      FMM M2M operation requires z to be OZ coordinate from new to old.
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_m2m_bessel_derivative_ztranslate_work(src_sk, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: src_sk(p+2), vscales((p+2)*(p+2)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+2)*(p+2))
    ! Temporary workspace
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, tmp3, res1, res2, pow_r1
    integer :: j, k, n, indj, indn, l, NZ, ierr
    ! Pointers for temporary values of powers
    real(dp) :: fact(2*p+1), fact2(p+2)
    complex(dp) :: complex_argument
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_m = zero
        else
            dst_m = beta * dst_m
        end if
        return
    end if
    ! Get factorials
    fact(1) = one
    do j = 1, 2*p
        fact(j+1) = dble(j) * fact(j)
    end do
    fact2(1) = one
    do j = 1, p+1
        fact2(j+1) = dble(2*j+1) * fact2(j)
    end do
    ! Now alpha is non-zero
    ! Do actual M2M
    ! Overwrite output if beta is zero
    if (beta .eq. zero) then
        dst_m = zero
        ! j=0, k=0, n=1, l=0
        j = 0
        k = 0
        n = 1
        l = 0
        indj = 1
        indn = 3
        dst_m(1) = vscales(3) / src_sk(2) / vscales(1) * src_sk(1) / fact2(2) * &
            & src_m(3) * alpha
        ! j=1..p-1
        do j = 1, p-1
            ! Offset for dst_m
            indj = j*j + j + 1
            ! k = 0
            k = 0
            res1 = zero
            ! n=j-1
            n = j-1
            ! Offset for src_m
            indn = n*n + n + 1
            tmp1 = vscales(indn) * &
                & dble(2*j+1) * fact(j+1) * fact(n+1) / &
                & src_sk(n+1) / vscales(indj) * src_sk(j+1) * alpha
            tmp3 = zero
            ! l=n
            l = n
            tmp2 = (two**(-l)) / &
                & fact(l+1)*fact(2*l+1)/ &
                & fact(l+1)/fact(l+1)/ &
                & fact2(j+1)
            tmp3 = tmp3 + tmp2
            tmp3 = tmp3 * tmp1
            res1 = res1 + tmp3*src_m(indn)
            ! n=j+1
            n = j+1
            ! Offset for src_m
            indn = n*n + n + 1
            tmp1 = vscales(indn) * &
                & dble(2*j+1) * fact(j+1) * fact(n+1) / &
                & src_sk(n+1) / vscales(indj) * src_sk(j+1) * alpha
            tmp3 = zero
            ! l=j
            l = j
            tmp2 = (two**(-l)) / &
                & fact(l+1)*fact(2*l+1)/ &
                & fact(l+1)/fact(l+1)/ &
                & fact2(j+2)
            tmp3 = tmp3 + tmp2
            tmp3 = tmp3 * tmp1
            res1 = res1 + tmp3*src_m(indn)
            dst_m(indj) = res1
            ! k=1..j-1
            do k = 1, j-1
                res1 = zero
                res2 = zero
                ! n=j-1
                n = j-1
                ! Offset for src_m
                indn = n*n + n + 1
                tmp1 = vscales(indn+k) * &
                    & dble(2*j+1) * fact(j-k+1) * fact(n+k+1) / &
                    & src_sk(n+1) / vscales(indj+k) * src_sk(j+1) * alpha
                tmp3 = zero
                l = n
                tmp2 = (two**(-l)) / &
                    & fact(l+k+1)*fact(2*l+1)/ &
                    & fact(l+1)/fact(l-k+1)/ &
                    & fact2(j+1)
                tmp3 = tmp3 + tmp2
                tmp3 = tmp3 * tmp1
                res1 = res1 + tmp3*src_m(indn+k)
                res2 = res2 + tmp3*src_m(indn-k)
                ! n=j+1
                n = j+1
                ! Offset for src_m
                indn = n*n + n + 1
                tmp1 = vscales(indn+k) * &
                    & dble(2*j+1) * fact(j-k+1) * fact(n+k+1) / &
                    & src_sk(n+1) / vscales(indj+k) * src_sk(j+1) * alpha
                tmp3 = zero
                l = j
                tmp2 = (two**(-l)) / &
                    & fact(l+k+1)*fact(2*l+1)/ &
                    & fact(l+1)/fact(l-k+1)/ &
                    & fact2(j+2)
                tmp3 = tmp3 + tmp2
                tmp3 = tmp3 * tmp1
                res1 = res1 + tmp3*src_m(indn+k)
                res2 = res2 + tmp3*src_m(indn-k)
                dst_m(indj+k) = res1
                dst_m(indj-k) = res2
            end do
            ! k=j
            k = j
            res1 = zero
            res2 = zero
            ! n=j+1
            n = j+1
            ! Offset for src_m
            indn = n*n + n + 1
            tmp1 = vscales(indn+k) * &
                & dble(2*j+1) * fact(n+k+1) / &
                & src_sk(n+1) / vscales(indj+k) * src_sk(j+1) * alpha
            tmp3 = zero
            l = j
            tmp2 = (two**(-l)) / &
                & fact(l+k+1)*fact(2*l+1)/ &
                & fact(l+1)/ &
                & fact2(j+2)
            tmp3 = tmp3 + tmp2
            tmp3 = tmp3 * tmp1
            res1 = res1 + tmp3*src_m(indn+k)
            res2 = res2 + tmp3*src_m(indn-k)
            dst_m(indj+k) = res1
            dst_m(indj-k) = res2
        end do
        ! j=p, n=p-1, l=p-1
        j = p
        n = p-1
        ! Offset for dst_m
        indj = j*j + j + 1
        ! k = 0
        k = 0
        res1 = zero
        ! Offset for src_m
        indn = n*n + n + 1
        tmp1 = vscales(indn) * &
            & dble(2*j+1) * fact(j+1) * fact(n+1) / &
            & src_sk(n+1) / vscales(indj) * src_sk(j+1) * alpha
        tmp3 = zero
        l = p-1
        tmp2 = (two**(-l)) / &
            & fact(l+1)*fact(2*l+1)/ &
            & fact(l+1)/fact(l+1)/ &
            & fact2(j+1)
        tmp3 = tmp3 + tmp2
        tmp3 = tmp3 * tmp1
        res1 = res1 + tmp3*src_m(indn)
        dst_m(indj) = res1
        ! k=1..p-1
        do k = 1, p-1
            res1 = zero
            res2 = zero
            ! n=j-1
            n = j-1
            ! Offset for src_m
            indn = n*n + n + 1
            tmp1 = vscales(indn+k) * &
                & dble(2*j+1) * fact(j-k+1) * fact(n+k+1) / &
                & src_sk(n+1) / vscales(indj+k) * src_sk(j+1) * alpha
            tmp3 = zero
            l = n
            tmp2 = (two**(-l)) / &
                & fact(l+k+1)*fact(2*l+1)/ &
                & fact(l+1)/fact(l-k+1)/fact(j-l+1)/fact(n-l+1)/ &
                & fact2(j+1)
            tmp3 = tmp3 + tmp2
            tmp3 = tmp3 * tmp1
            res1 = res1 + tmp3*src_m(indn+k)
            res2 = res2 + tmp3*src_m(indn-k)
            dst_m(indj+k) = res1
            dst_m(indj-k) = res2
        end do
        ! j=p,k=p is impossible because n=p-1 or n=p+1 is impossible
        ! j=p+1, n=p, l=p
        j = p+1
        n = p
        ! Offset for dst_m
        indj = j*j + j + 1
        ! k = 0
        k = 0
        res1 = zero
        ! Offset for src_m
        indn = n*n + n + 1
        tmp1 = vscales(indn) * &
            & dble(2*j+1) * fact(j+1) * fact(n+1) / &
            & src_sk(n+1) / vscales(indj) * src_sk(j+1) * alpha
        tmp3 = zero
        l = n
        tmp2 = (two**(-l)) / &
            & fact(l+1)*fact(2*l+1)/ &
            & fact(l+1)/fact(l+1)/fact(j-l+1)/fact(n-l+1)/ &
            & fact2(j+1)
        tmp3 = tmp3 + tmp2
        tmp3 = tmp3 * tmp1
        res1 = res1 + tmp3*src_m(indn)
        dst_m(indj) = res1
        ! k != 0
        do k = 1, p
            res1 = zero
            res2 = zero
            ! n=j-1
            n = p
            ! Offset for src_m
            indn = n*n + n + 1
            tmp1 = vscales(indn+k) * &
                & dble(2*j+1) * fact(j-k+1) * fact(n+k+1) / &
                & src_sk(n+1) / vscales(indj+k) * src_sk(j+1) * alpha
            tmp3 = zero
            l = n
            tmp2 = (two**(-l)) / &
                & fact(l+k+1)*fact(2*l+1)/ &
                & fact(l+1)/fact(l-k+1)/fact(j-l+1)/fact(n-l+1)/ &
                & fact2(j+1)
            tmp3 = tmp3 + tmp2
            tmp3 = tmp3 * tmp1
            res1 = res1 + tmp3*src_m(indn+k)
            res2 = res2 + tmp3*src_m(indn-k)
            dst_m(indj+k) = res1
            dst_m(indj-k) = res2
        end do
    ! Update output if beta is non-zero
    else
        stop 111
    end if
end subroutine fmm_m2m_bessel_derivative_ztranslate_work

!> Adjoint M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of
!! multipole-to-multipole translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_ztranslate_adj(z, src_r, dst_r, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(p+1))
    ! Call corresponding work routine
    call fmm_m2m_ztranslate_adj_work(z, src_r, dst_r, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m, work)
end subroutine fmm_m2m_ztranslate_adj

!> Adjoint M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of
!! multipole-to-multipole translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_m2m_ztranslate_adj_work(z, src_r, dst_r, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, res1, res2
    integer :: j, k, n, indj, indn, indjk1, indjk2
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r1(:), pow_r2(:)
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_m = zero
        else
            dst_m = beta * dst_m
        end if
        return
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        ! Prepare pointers
        pow_r1(1:p+1) => work(1:p+1)
        pow_r2(1:p+1) => work(p+2:2*(p+1))
        ! Get powers of r1 and r2
        r1 = dst_r / src_r
        r2 = -z / src_r
        pow_r1(1) = r1
        pow_r2(1) = one
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Do actual adjoint M2M
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            do n = 0, p
                ! Offset for dst_m
                indn = n*n + n + 1
                tmp1 = alpha * pow_r1(n+1) / vscales(indn)
                ! k = 0
                res1 = zero
                do j = n, p
                    ! Offset for src_m
                    indj = j*j + j + 1
                    ! Offsets for vcnk
                    indjk1 = j*(j+1)/2 + 1
                    tmp2 = tmp1 * vscales(indj) * pow_r2(j-n+1) * &
                        & vcnk(indjk1+j-n)**2
                    res1 = res1 + tmp2*src_m(indj)
                end do
                dst_m(indn) = res1
                ! k != 1
                do k = 1, n
                    res1 = zero
                    res2 = zero
                    do j = n, p
                        ! Offset for src_m
                        indj = j*j + j + 1
                        ! Offsets for vcnk
                        indjk1 = (j-k)*(j-k+1)/2 + 1
                        indjk2 = (j+k)*(j+k+1)/2 + 1
                        tmp2 = tmp1 * vscales(indj) * pow_r2(j-n+1) * &
                            & vcnk(indjk1+j-n) * vcnk(indjk2+j-n)
                        res1 = res1 + tmp2*src_m(indj+k)
                        res2 = res2 + tmp2*src_m(indj-k)
                    end do
                    dst_m(indn+k) = res1
                    dst_m(indn-k) = res2
                end do
            end do
        ! Update output if beta is non-zero
        else
            do n = 0, p
                ! Offset for dst_m
                indn = n*n + n + 1
                tmp1 = alpha * pow_r1(n+1) / vscales(indn)
                ! k = 0
                res1 = zero
                do j = n, p
                    ! Offset for src_m
                    indj = j*j + j + 1
                    ! Offsets for vcnk
                    indjk1 = j*(j+1)/2 + 1
                    tmp2 = tmp1 * vscales(indj) * pow_r2(j-n+1) * &
                        & vcnk(indjk1+j-n)**2
                    res1 = res1 + tmp2*src_m(indj)
                end do
                dst_m(indn) = beta*dst_m(indn) + res1
                ! k != 1
                do k = 1, n
                    res1 = zero
                    res2 = zero
                    do j = n, p
                        ! Offset for src_m
                        indj = j*j + j + 1
                        ! Offsets for vcnk
                        indjk1 = (j-k)*(j-k+1)/2 + 1
                        indjk2 = (j+k)*(j+k+1)/2 + 1
                        tmp2 = tmp1 * vscales(indj) * pow_r2(j-n+1) * &
                            & vcnk(indjk1+j-n) * vcnk(indjk2+j-n)
                        res1 = res1 + tmp2*src_m(indj+k)
                        res2 = res2 + tmp2*src_m(indj-k)
                    end do
                    dst_m(indn+k) = beta*dst_m(indn+k) + res1
                    dst_m(indn-k) = beta*dst_m(indn-k) + res2
                end do
            end do
        end if
    ! If harmonics are located at the same point
    else
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            r1 = dst_r / src_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                res1 = zero
                do k = indj-j, indj+j
                    dst_m(k) = tmp1 * src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        ! Update output if beta is non-zero
        else
            r1 = dst_r / src_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                res1 = zero
                do k = indj-j, indj+j
                    dst_m(k) = beta*dst_m(k) + tmp1*src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        end if
    end if
end subroutine fmm_m2m_ztranslate_adj_work

!> Scale M2M, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation when harmonics are centered in the same point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_scale(src_r, dst_r, p, alpha, src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: src_r, dst_r, alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    r1 = src_r / dst_r
    tmp1 = alpha * r1
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_m(k) = dst_m(k) + src_m(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_m2m_scale

!> Adjoint scale M2M, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of a
!! multipole-to-multipole translation when harmonics are centered in the same
!! point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_scale_adj(src_r, dst_r, p, alpha, src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: src_r, dst_r, alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    r1 = dst_r / src_r
    tmp1 = alpha * r1
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_m(k) = dst_m(k) + src_m(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_m2m_scale_adj

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_rotation(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p + 19*p + 8)
    ! Call corresponding work routine
    call fmm_m2m_rotation_work(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
end subroutine fmm_m2m_rotation

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_m2m_rotation_work(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_m(:), tmp_m2(:), vcos(:), vsin(:)
    ! Convert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        ! Workspace here is 2*(p+1)
        call fmm_m2m_ztranslate_work(c(3), src_r, dst_r, p, vscales, vcnk, &
            & alpha, src_m, beta, dst_m, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_m(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_m2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_m, zero, tmp_m)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_m, zero, &
        & tmp_m2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_m2m_ztranslate_work(rho, src_r, dst_r, p, vscales, vcnk, one, &
        & tmp_m2, zero, tmp_m, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_m, zero, tmp_m2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_rotation_work

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_bessel_rotation(c, src_r, dst_r, kappa, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta, kappa
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p + 19*p + 8), src_sk(p+1), dst_sk(p+1), s1, s2
    complex(dp) :: work_complex(2*p+1), z1, z2
    integer :: NZ, ierr
    ! Compute Bessel functions
    z1 = src_r * kappa
    z2 = dst_r * kappa
    s1 = sqrt(2/(pi*real(z1)))
    s2 = sqrt(2/(pi*real(z2)))
    call cbesk(z1, pt5, 1, 1, work_complex(1), NZ, ierr)
    src_sk(1) = s1 * real(work_complex(1))
    call cbesk(z2, pt5, 1, 1, work_complex(1), NZ, ierr)
    dst_sk(1) = s2 * real(work_complex(1))
    if (p .gt. 0) then
        call cbesk(z1, 1.5d0, 1, p, work_complex(2:p+1), NZ, ierr)
        src_sk(2:p+1) = s1 * real(work_complex(2:p+1))
        call cbesk(z2, 1.5d0, 1, p, work_complex(2:p+1), NZ, ierr)
        dst_sk(2:p+1) = s2 * real(work_complex(2:p+1))
    end if
    ! Call corresponding work routine
    call fmm_m2m_bessel_rotation_work(c*kappa, src_sk, dst_sk, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work, work_complex)
end subroutine fmm_m2m_bessel_rotation

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from old to new centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_m2m_bessel_rotation_work(c, src_sk, dst_sk, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work, work_complex)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_sk(p+1), dst_sk(p+1), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    complex(dp), intent(out) :: work_complex(2*p+1)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_m(:), tmp_m2(:), vcos(:), vsin(:)
    ! Convert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    !if (stheta .eq. zero) then
    !    ! Workspace here is 2*(p+1)
    !    call fmm_m2m_bessel_ztranslate_work(c(3), src_r, dst_r, p, vscales, vcnk, &
    !        & alpha, src_m, beta, dst_m, work)
    !    return
    !end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_m(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_m2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_m, zero, tmp_m)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_m, zero, &
        & tmp_m2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_m2m_bessel_ztranslate_work(rho, src_sk, dst_sk, p, vscales, vcnk, one, &
        & tmp_m2, zero, tmp_m, work, work_complex)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_m, zero, tmp_m2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_bessel_rotation_work

!> Adjoint M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is aa adjoint matrix of a
!! multipole-to-multipole translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] src_m: expansion in old harmonics
!! @param[inout] dst_m: expansion in new harmonics
subroutine fmm_m2m_rotation_adj(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p + 19*p + 8)
    ! Call corresponding work routine
    call fmm_m2m_rotation_adj_work(c, src_r, dst_r, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
end subroutine fmm_m2m_rotation_adj

!> Adjoint M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is aa adjoint matrix of a
!! multipole-to-multipole translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] src_m: expansion in old harmonics
!! @param[inout] dst_m: expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_m2m_rotation_adj_work(c, src_r, dst_r, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_m, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_m(:), tmp_m2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. 0) then
        ! Workspace here is 2*(p+1)
        call fmm_m2m_ztranslate_adj_work(c(3), src_r, dst_r, p, vscales, &
            & vcnk, alpha, src_m, beta, dst_m, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_m(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_m2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_m, zero, tmp_m)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_m, zero, &
        & tmp_m2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_m2m_ztranslate_adj_work(rho, src_r, dst_r, p, vscales, vcnk, &
        & one, tmp_m2, zero, tmp_m, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_m, zero, tmp_m2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_rotation_adj_work

subroutine fmm_m2m_bessel_grad(p, sph_sk, vscales, vcnk, sph_m, sph_m_grad)
    !! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: sph_sk(p+2), vscales((p+2)**2), vcnk((2*p+3)*(p+2)), &
        & sph_m((p+1)**2)
    !! Output
    real(dp), intent(out) :: sph_m_grad((p+2)**2, 3)
    !! Local variables
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    real(dp) :: tmp((p+1)**2), tmp2((p+2)**2)
    ! Set coordinate transformations
    zx_coord_transform = zero
    zx_coord_transform(3, 2) = one
    zx_coord_transform(2, 3) = one
    zx_coord_transform(1, 1) = one
    zy_coord_transform = zero
    zy_coord_transform(1, 2) = one
    zy_coord_transform(2, 1) = one
    zy_coord_transform(3, 3) = one
    ! Transform input harmonics for OX axis
    call fmm_sph_transform(p, zx_coord_transform, one, sph_m, zero, tmp)
    ! Differentiate M2M along OZ
    call fmm_m2m_bessel_derivative_ztranslate_work(sph_sk, p, vscales, vcnk, &
        & one, tmp, zero, tmp2)
    ! Transform into output harmonics for OX axis
    call fmm_sph_transform(p+1, zx_coord_transform, one, tmp2, zero, &
        & sph_m_grad(:, 1))
    ! Transform input harmonics for OY axis
    call fmm_sph_transform(p, zy_coord_transform, one, sph_m, zero, tmp)
    ! Differentiate M2M along OZ
    call fmm_m2m_bessel_derivative_ztranslate_work(sph_sk, p, vscales, vcnk, &
        & one, tmp, zero, tmp2)
    ! Transform into output harmonics for OY axis
    call fmm_sph_transform(p+1, zy_coord_transform, one, tmp2, zero, &
        & sph_m_grad(:, 2))
    ! Differentiate for OZ axis
    call fmm_m2m_bessel_derivative_ztranslate_work(sph_sk, p, vscales, vcnk, &
        & one, sph_m, zero, sph_m_grad(:, 3))
end subroutine fmm_m2m_bessel_grad

!> Direct L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of local-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_ztranslate(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(p+1))
    ! Call corresponding work routine
    call fmm_l2l_ztranslate_work(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
end subroutine fmm_l2l_ztranslate

!> Direct L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of local-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_l2l_ztranslate_work(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2
    integer :: j, k, n, indj, indn
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r1(:), pow_r2(:)
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        ! Prepare pointers
        pow_r1(1:p+1) => work(1:p+1)
        pow_r2(1:p+1) => work(p+2:2*(p+1))
        ! Get powers of r1 and r2
        r1 = z / src_r
        r2 = dst_r / z
        pow_r1(1) = 1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Do actual L2L
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = alpha * pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1) * &
                    & vscales(indj)
                do n = j, p
                    indn = n*n + n + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) * &
                        & vfact(n-k+1) * vfact(n+k+1) / vfact(n-j+1) / &
                        & vfact(n-j+1)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    if (k .eq. 0) then
                        dst_l(indj) = dst_l(indj) + tmp2*src_l(indn)
                    else
                        dst_l(indj+k) = dst_l(indj+k) + tmp2*src_l(indn+k)
                        dst_l(indj-k) = dst_l(indj-k) + tmp2*src_l(indn-k)
                    end if
                end do
            end do
        end do
    ! If harmonics are located at the same point
    else
        r1 = dst_r / src_r
        tmp1 = alpha
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_l2l_ztranslate_work

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_l2l_bessel_ztranslate(z, src_si, dst_si, p, vscales, vcnk, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_si(p+1), dst_si(p+1), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*p+1)
    complex(dp) :: work_complex(2*p+1)
    ! Call corresponding work routine
    call fmm_l2l_bessel_ztranslate_work(z, src_si, dst_si, p, vscales, vcnk, &
        & alpha, src_l, beta, dst_l, work, work_complex)
end subroutine fmm_l2l_bessel_ztranslate

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from old to new centers of harmonics. Standard
!!      FMM M2M operation requires z to be OZ coordinate from new to old.
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_l2l_bessel_ztranslate_work(z, src_si, dst_si, p, vscales, vcnk, &
        & alpha, src_l, beta, dst_l, work, work_complex)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_si(p+1), dst_si(p+1), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*p+1)
    complex(dp), intent(out) :: work_complex(2*p+1)
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, tmp3, res1, res2, pow_r1
    integer :: j, k, n, indj, indn, l, NZ, ierr
    ! Pointers for temporary values of powers
    real(dp) :: fact(2*p+1)
    complex(dp) :: complex_argument
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_l = zero
        else
            dst_l = beta * dst_l
        end if
        return
    end if
    ! Get factorials
    fact(1) = one
    do j = 1, 2*p
        fact(j+1) = dble(j) * fact(j)
    end do
    ! Now alpha is non-zero
    ! If harmonics have different centers
    if (z .ne. 0) then
        complex_argument = z
        call cbesi(complex_argument, pt5, 1, 1, work_complex(1), NZ, ierr)
        work(1) = sqrt(pi/(2*z)) * real(work_complex(1))
        if (p .gt. 0) then
            call cbesi(complex_argument, 1.5d0, 1, 2*p, &
                & work_complex(2:2*p+1), NZ, ierr)
            work(2:2*p+1) = sqrt(pi/(2*z)) * real(work_complex(2:2*p+1))
        end if
        ! Do actual M2M
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            dst_l = zero
            do j = 0, p
                ! Offset for dst_l
                indj = j*j + j + 1
                ! k = 0
                k = 0
                res1 = zero
                do n = 0, p
                    ! Offset for src_l
                    indn = n*n + n + 1
                    tmp1 = vscales(indn) * ((-one)**(j+n)) * &
                        & dble(2*j+1) * fact(j+1) * fact(n+1) / &
                        & src_si(n+1) / vscales(indj) * dst_si(j+1) * alpha
                    tmp3 = zero
                    do l = 0, min(j, n)
                        tmp2 = (two**(-l)) / &
                            & fact(l+1)*fact(2*l+1)/ &
                            & fact(l+1)/fact(l+1)/fact(j-l+1)/fact(n-l+1)* &
                            & (work(j+n-l+1) / (z**l))
                        tmp3 = tmp3 + tmp2
                    end do
                    tmp3 = tmp3 * tmp1
                    res1 = res1 + tmp3*src_l(indn)
                end do
                dst_l(indj) = res1
                ! k != 0
                do k = 1, j
                    res1 = zero
                    res2 = zero
                    do n = k, p
                        ! Offset for src_l
                        indn = n*n + n + 1
                        tmp1 = vscales(indn+k) * ((-one)**(j+n)) * &
                            & dble(2*j+1) * fact(j-k+1) * fact(n+k+1) / &
                            & src_si(n+1) / vscales(indj+k) * dst_si(j+1) * alpha
                        tmp3 = zero
                        do l = k, min(j, n)
                            tmp2 = (two**(-l)) / &
                                & fact(l+k+1)*fact(2*l+1)/ &
                                & fact(l+1)/fact(l-k+1)/fact(j-l+1)/fact(n-l+1)* &
                                & (work(j+n-l+1) / (z**l))
                            tmp3 = tmp3 + tmp2
                        end do
                        tmp3 = tmp3 * tmp1
                        res1 = res1 + tmp3*src_l(indn+k)
                        res2 = res2 + tmp3*src_l(indn-k)
                    end do
                    dst_l(indj+k) = res1
                    dst_l(indj-k) = res2
                end do
            end do
        ! Update output if beta is non-zero
        else
            stop 111
        end if
    ! If harmonics are located at the same point
    else
        stop 222
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            !r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_l(k) = tmp1 * src_l(k)
                end do
                tmp1 = tmp1 * r1
            end do
        ! Update output if beta is non-zero
        else
            !r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_l(k) = beta*dst_l(k) + tmp1*src_l(k)
                end do
                tmp1 = tmp1 * r1
            end do
        end if
    end if
end subroutine fmm_l2l_bessel_ztranslate_work

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from old to new centers of harmonics. Standard
!!      FMM M2M operation requires z to be OZ coordinate from new to old.
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_l2l_bessel_derivative_ztranslate_work(src_si, p, vscales, vcnk, &
        & alpha, src_l, beta, dst_l)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: src_si(p+2), vscales((p+2)*(p+2)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+2)*(p+2))
    ! Temporary workspace
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, tmp3, res1, res2, pow_r1
    integer :: j, k, n, indj, indn, l, NZ, ierr
    ! Pointers for temporary values of powers
    real(dp) :: fact(2*p+1), fact2(p+2)
    complex(dp) :: complex_argument
    ! Derivative of the L2L is just a negative derivative of the M2M
    call fmm_m2m_bessel_derivative_ztranslate_work(src_si, p, vscales, vcnk, &
        & -alpha, src_l, beta, dst_l)
end subroutine fmm_l2l_bessel_derivative_ztranslate_work

!> Adjoint L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of
!! local-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_ztranslate_adj(z, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*(p+1))
    ! Call corresponding work routine
    call fmm_l2l_ztranslate_adj_work(z, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l, work)
end subroutine fmm_l2l_ztranslate_adj

!> Adjoint L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of
!! local-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_l2l_ztranslate_adj_work(z, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2
    integer :: j, k, n, indj, indn
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r1(:), pow_r2(:)
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        ! Prepare pointers
        pow_r1(1:p+1) => work(1:p+1)
        pow_r2(1:p+1) => work(p+2:2*(p+1))
        ! Get powers of r1 and r2
        r1 = -z / dst_r
        r2 = -src_r / z
        pow_r1(1) = one
        pow_r2(1) = one
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = alpha * pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1) * &
                    & vscales(indj)
                do n = j, p
                    indn = n*n + n + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) * &
                        & vfact(n-k+1) * vfact(n+k+1) / vfact(n-j+1) / &
                        & vfact(n-j+1)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    if (k .eq. 0) then
                        dst_l(indn) = dst_l(indn) + tmp2*src_l(indj)
                    else
                        dst_l(indn+k) = dst_l(indn+k) + tmp2*src_l(indj+k)
                        dst_l(indn-k) = dst_l(indn-k) + tmp2*src_l(indj-k)
                    end if
                end do
            end do
        end do
    ! If harmonics are located at the same point
    else
        r1 = src_r / dst_r
        tmp1 = alpha
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_l2l_ztranslate_adj_work

!> Scale L2L, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of a local-to-local
!! translation when harmonics are centered in the same point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_l: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_l: Multipole expansion in new harmonics
subroutine fmm_l2l_scale(src_r, dst_r, p, alpha, src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: src_r, dst_r, alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    r1 = dst_r / src_r
    tmp1 = alpha
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_l(k) = dst_l(k) + src_l(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_l2l_scale

!> Adjoint scale L2L, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation when harmonics are centered in the same point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_l: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_l: Multipole expansion in new harmonics
subroutine fmm_l2l_scale_adj(src_r, dst_r, p, alpha, src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: src_r, dst_r, alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    r1 = src_r / dst_r
    tmp1 = alpha
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_l(k) = dst_l(k) + src_l(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_l2l_scale_adj

!> Direct L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of a local-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_rotation(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p+19*p+8)
    ! Call corresponding work routine
    call fmm_l2l_rotation_work(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
end subroutine fmm_l2l_rotation

!> Direct L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of a local-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_l2l_rotation_work(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_l(:), tmp_l2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        ! Workspace here is 2*(p+1)
        call fmm_l2l_ztranslate_work(c(3), src_r, dst_r, p, vscales, vfact, &
            & alpha, src_l, beta, dst_l, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_l(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_l2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_l, zero, tmp_l)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_l, zero, &
        & tmp_l2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_l2l_ztranslate_work(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_l2, zero, tmp_l, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_l, zero, tmp_l2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_rotation_work

!> Direct L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_l2l_bessel_rotation(c, src_r, dst_r, kappa, p, vscales, vcnk, alpha, &
        & src_l, beta, dst_l)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_l((p+1)*(p+1)), beta, kappa
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p + 19*p + 8), src_si(p+1), dst_si(p+1), s1, s2
    complex(dp) :: work_complex(2*p+1), z1, z2
    integer :: NZ, ierr
    ! Compute Bessel functions
    z1 = src_r * kappa
    z2 = dst_r * kappa
    s1 = sqrt(pi/(2*real(z1)))
    s2 = sqrt(pi/(2*real(z2)))
    call cbesi(z1, pt5, 1, 1, work_complex(1), NZ, ierr)
    src_si(1) = s1 * real(work_complex(1))
    call cbesi(z2, pt5, 1, 1, work_complex(1), NZ, ierr)
    dst_si(1) = s2 * real(work_complex(1))
    if (p .gt. 0) then
        call cbesi(z1, 1.5d0, 1, p, work_complex(2:p+1), NZ, ierr)
        src_si(2:p+1) = s1 * real(work_complex(2:p+1))
        call cbesi(z2, 1.5d0, 1, p, work_complex(2:p+1), NZ, ierr)
        dst_si(2:p+1) = s2 * real(work_complex(2:p+1))
    end if
    ! Call corresponding work routine
    call fmm_l2l_bessel_rotation_work(c*kappa, src_si, dst_si, p, vscales, vcnk, alpha, &
        & src_l, beta, dst_l, work, work_complex)
end subroutine fmm_l2l_bessel_rotation

!> Direct L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from old to new centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_l2l_bessel_rotation_work(c, src_si, dst_si, p, vscales, vcnk, alpha, &
        & src_l, beta, dst_l, work, work_complex)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_si(p+1), dst_si(p+1), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    complex(dp), intent(out) :: work_complex(2*p+1)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_l(:), tmp_l2(:), vcos(:), vsin(:)
    ! Convert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    !if (stheta .eq. zero) then
    !    ! Workspace here is 2*(p+1)
    !    call fmm_m2m_bessel_ztranslate_work(c(3), src_r, dst_r, p, vscales, vcnk, &
    !        & alpha, src_m, beta, dst_m, work)
    !    return
    !end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_l(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_l2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_l, zero, tmp_l)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_l, zero, &
        & tmp_l2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_l2l_bessel_ztranslate_work(rho, src_si, dst_si, p, vscales, vcnk, one, &
        & tmp_l2, zero, tmp_l, work, work_complex)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_l, zero, tmp_l2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_bessel_rotation_work

subroutine fmm_l2l_bessel_grad(p, sph_si, vscales, vcnk, sph_l, sph_l_grad)
    !! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: sph_si(p+2), vscales((p+2)**2), vcnk((2*p+3)*(p+2)), &
        & sph_l((p+1)**2)
    !! Output
    real(dp), intent(out) :: sph_l_grad((p+2)**2, 3)
    !! Local variables
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    real(dp) :: tmp((p+1)**2), tmp2((p+2)**2)
    ! Set coordinate transformations
    zx_coord_transform = zero
    zx_coord_transform(3, 2) = one
    zx_coord_transform(2, 3) = one
    zx_coord_transform(1, 1) = one
    zy_coord_transform = zero
    zy_coord_transform(1, 2) = one
    zy_coord_transform(2, 1) = one
    zy_coord_transform(3, 3) = one
    ! Transform input harmonics for OX axis
    call fmm_sph_transform(p, zx_coord_transform, one, sph_l, zero, tmp)
    ! Differentiate L2L along OZ
    call fmm_l2l_bessel_derivative_ztranslate_work(sph_si, p, vscales, vcnk, &
        & one, tmp, zero, tmp2)
    ! Transform into output harmonics for OX axis
    call fmm_sph_transform(p+1, zx_coord_transform, one, tmp2, zero, &
        & sph_l_grad(:, 1))
    ! Transform input harmonics for OY axis
    call fmm_sph_transform(p, zy_coord_transform, one, sph_l, zero, tmp)
    ! Differentiate M2M along OZ
    call fmm_l2l_bessel_derivative_ztranslate_work(sph_si, p, vscales, vcnk, &
        & one, tmp, zero, tmp2)
    ! Transform into output harmonics for OY axis
    call fmm_sph_transform(p+1, zy_coord_transform, one, tmp2, zero, &
        & sph_l_grad(:, 2))
    ! Differentiate for OZ axis
    call fmm_l2l_bessel_derivative_ztranslate_work(sph_si, p, vscales, vcnk, &
        & one, sph_l, zero, sph_l_grad(:, 3))
end subroutine fmm_l2l_bessel_grad

!> Adjoint L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_rotation_adj(c, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p + 19*p + 8)
    ! Call corresponding work routine
    call fmm_l2l_rotation_adj_work(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l, work)
end subroutine fmm_l2l_rotation_adj

!> Adjoint L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_l2l_rotation_adj_work(c, src_r, dst_r, p, vscales, vfact, &
        & alpha, src_l, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_l(:), tmp_l2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        ! Workspace here is 2*(p+1)
        call fmm_l2l_ztranslate_adj_work(c(3), src_r, dst_r, p, vscales, &
            & vfact, alpha, src_l, beta, dst_l, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_l(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_l2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_l, zero, tmp_l)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_l, zero, &
        & tmp_l2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_l2l_ztranslate_adj_work(rho, src_r, dst_r, p, vscales, vfact, &
        & one, tmp_l2, zero, tmp_l, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_l, zero, tmp_l2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_rotation_adj_work

!> Direct M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of multipole-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (multipole) harmonics
!! @param[in] dst_r: Radius of new (local) harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old (multipole) harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new (local) harmonics
subroutine fmm_m2l_ztranslate(z, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l)
    ! Inputs
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
        & beta
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Temporary workspace
    real(dp) :: work((pm+2)*(pm+1))
    ! Call corresponding work routine
    call fmm_m2l_ztranslate_work(z, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
end subroutine fmm_m2l_ztranslate

!> Direct M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of multipole-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (multipole) harmonics
!! @param[in] dst_r: Radius of new (local) harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old (multipole) harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new (local) harmonics
!! @param[out] work: Temporary workspace of a size (pm+2)*(pm+1)
subroutine fmm_m2l_ztranslate_work(z, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
        & beta
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work((pm+2)*(pm+1))
    ! Local variables
    real(dp) :: tmp1, tmp2, r1, r2, res1, res2, pow_r2
    integer :: j, k, n, indj, indjn, indk1, indk2
    ! Pointers for temporary values of powers
    real(dp), pointer :: src_m2(:), pow_r1(:)
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_l = zero
        else
            dst_l = beta * dst_l
        end if
        return
    end if
    ! Now alpha is non-zero
    ! z cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (z .eq. zero) then
        return
    end if
    ! Prepare pointers
    n = (pm+1) ** 2
    src_m2(1:n) => work(1:n)
    pow_r1(1:pm+1) => work(n+1:n+pm+1)
    ! Get powers of r1 and r2
    r1 = src_r / z
    r2 = dst_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    tmp1 = abs(r1)
    do j = 0, pm
        indj = j*j + j + 1
        pow_r1(j+1) = tmp1 / vscales(indj)
        tmp1 = tmp1 * r1
    end do
    pow_r2 = one
    ! Reorder source harmonics from (degree, order) to (order, degree)
    ! Zero order k=0 at first
    do j = 0, pm
        indj = j*j + j + 1
        src_m2(j+1) = pow_r1(j+1) * src_m(indj)
    end do
    ! Non-zero orders next, a positive k followed by a negative -k
    indk1 = pm + 2
    do k = 1, pm
        n = pm - k + 1
        indk2 = indk1 + n
        !!GCC$ unroll 4
        do j = k, pm
            indj = j*j + j + 1
            src_m2(indk1+j-k) = pow_r1(j+1) * src_m(indj+k)
            src_m2(indk2+j-k) = pow_r1(j+1) * src_m(indj-k)
        end do
        indk1 = indk2 + n
    end do
    ! Do actual M2L
    ! Overwrite output if beta is zero
    if (beta .eq. zero) then
        do j = 0, pl
            ! Offset for dst_l
            indj = j*j + j + 1
            ! k = 0
            tmp1 = alpha * vscales(indj) * pow_r2
            pow_r2 = pow_r2 * r2
            res1 = zero
            !!GCC$ unroll 4
            do n = 0, pm
                res1 = res1 + m2l_ztranslate_coef(n+1, 1, j+1)*src_m2(n+1)
            end do
            dst_l(indj) = tmp1 * res1
            ! k != 0
            !!GCC$ unroll 4
            do k = 1, j
                ! Offsets for src_m2
                indk1 = pm + 2 + (2*pm-k+2)*(k-1)
                indk2 = indk1 + pm - k + 1
                res1 = zero
                res2 = zero
                !!GCC$ unroll 4
                do n = k, pm
                    res1 = res1 + &
                        & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                        & src_m2(indk1+n-k)
                    res2 = res2 + &
                        & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                        & src_m2(indk2+n-k)
                end do
                dst_l(indj+k) = tmp1 * res1
                dst_l(indj-k) = tmp1 * res2
            end do
        end do
    else
        do j = 0, pl
            ! Offset for dst_l
            indj = j*j + j + 1
            ! k = 0
            tmp1 = alpha * vscales(indj) * pow_r2
            pow_r2 = pow_r2 * r2
            res1 = zero
            do n = 0, pm
                res1 = res1 + m2l_ztranslate_coef(n+1, 1, j+1)*src_m2(n+1)
            end do
            dst_l(indj) = beta*dst_l(indj) + tmp1*res1
            ! k != 0
            do k = 1, j
                ! Offsets for src_m2
                indk1 = pm + 2 + (2*pm-k+2)*(k-1)
                indk2 = indk1 + pm - k + 1
                res1 = zero
                res2 = zero
                do n = k, pm
                    res1 = res1 + &
                        & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                        & src_m2(indk1+n-k)
                    res2 = res2 + &
                        & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                        & src_m2(indk2+n-k)
                end do
                dst_l(indj+k) = beta*dst_l(indj+k) + tmp1*res1
                dst_l(indj-k) = beta*dst_l(indj-k) + tmp1*res2
            end do
        end do
    end if
end subroutine fmm_m2l_ztranslate_work

!> Adjoint M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of
!! multipole-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (local) harmonics
!! @param[in] dst_r: Radius of new (multipole) harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] m2l_ztranslate_adj_coef:
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old (local) harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new (multipole) harmonics
subroutine fmm_m2l_ztranslate_adj(z, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m)
    ! Inputs
    integer, intent(in) :: pl, pm
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), alpha, &
        & src_l((pl+1)*(pl+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Temporary workspace
    real(dp) :: work((pl+2)*(pl+1))
    ! Call corresponding work routine
    call fmm_m2l_ztranslate_adj_work(z, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
end subroutine fmm_m2l_ztranslate_adj

!> Adjoint M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of
!! multipole-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (local) harmonics
!! @param[in] dst_r: Radius of new (multipole) harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] m2l_ztranslate_adj_coef:
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old (local) harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new (multipole) harmonics
subroutine fmm_m2l_ztranslate_adj_work(z, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
    ! Inputs
    integer, intent(in) :: pl, pm
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), alpha, &
        & src_l((pl+1)*(pl+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work((pl+2)*(pl+1))
    ! Local variables
    real(dp) :: tmp1, tmp2, tmp3, r1, r2, pow_r1, res1, res2
    integer :: j, k, n, indj, indn, indk1, indk2
    ! Pointers for temporary values of powers
    real(dp), pointer :: pow_r2(:), src_l2(:)
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_m = zero
        else
            dst_m = beta * dst_m
        end if
        return
    end if
    ! Now alpha is non-zero
    ! z cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (z .eq. zero) then
        return
    end if
    ! Prepare pointers
    n = (pl+1) ** 2
    src_l2(1:n) => work(1:n)
    pow_r2(1:pl+1) => work(n+1:n+pl+1)
    ! Get powers of r1 and r2
    r1 = -dst_r / z
    r2 = -src_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    tmp1 = one
    do j = 0, pl
        indj = j*j + j + 1
        pow_r2(j+1) = tmp1 * vscales(indj)
        tmp1 = tmp1 * r2
    end do
    pow_r1 = abs(r1)
    ! Reorder source harmonics from (degree, order) to (order, degree)
    ! Zero order k=0 at first
    do j = 0, pl
        indj = j*j + j + 1
        src_l2(j+1) = pow_r2(j+1) * src_l(indj)
    end do
    ! Non-zero orders next, a positive k followed by a negative -k
    indk1 = pl + 2
    do k = 1, pl
        n = pl - k + 1
        indk2 = indk1 + n
        do j = k, pl
            indj = j*j + j + 1
            src_l2(indk1+j-k) = pow_r2(j+1) * src_l(indj+k)
            src_l2(indk2+j-k) = pow_r2(j+1) * src_l(indj-k)
        end do
        indk1 = indk2 + n
    end do
    ! Do actual adjoint M2L
    ! Overwrite output if beta is zero
    if (beta .eq. zero) then
        do n = 0, pm
            indn = n*n + n + 1
            ! k = 0
            tmp1 = alpha * pow_r1 / vscales(indn)
            pow_r1 = pow_r1 * r1
            res1 = zero
            do j = 0, pl
                indj = j*j + j + 1
                res1 = res1 + m2l_ztranslate_adj_coef(j+1, 1, n+1)*src_l2(j+1)
            end do
            dst_m(indn) = tmp1 * res1
            ! k != 0
            do k = 1, n
                ! Offsets for src_l2
                indk1 = pl + 2 + (2*pl-k+2)*(k-1)
                indk2 = indk1 + pl - k + 1
                res1 = zero
                res2 = zero
                do j = k, pl
                    indj = j*j + j + 1
                    res1 = res1 + &
                        & m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1)* &
                        & src_l2(indk1+j-k)
                    res2 = res2 + &
                        & m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1)* &
                        & src_l2(indk2+j-k)
                end do
                dst_m(indn+k) = tmp1 * res1
                dst_m(indn-k) = tmp1 * res2
            end do
        end do
    ! Update output if beta is non-zero
    else
        do n = 0, pm
            indn = n*n + n + 1
            ! k = 0
            tmp1 = alpha * pow_r1 / vscales(indn)
            pow_r1 = pow_r1 * r1
            res1 = zero
            do j = 0, pl
                indj = j*j + j + 1
                res1 = res1 + m2l_ztranslate_adj_coef(j+1, 1, n+1)*src_l2(j+1)
            end do
            dst_m(indn) = beta*dst_m(indn) + tmp1*res1
            ! k != 0
            do k = 1, n
                ! Offsets for src_l2
                indk1 = pl + 2 + (2*pl-k+2)*(k-1)
                indk2 = indk1 + pl - k + 1
                res1 = zero
                res2 = zero
                do j = k, pl
                    indj = j*j + j + 1
                    res1 = res1 + &
                        & m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1)* &
                        & src_l2(indk1+j-k)
                    res2 = res2 + &
                        & m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1)* &
                        & src_l2(indk2+j-k)
                end do
                dst_m(indn+k) = beta*dst_m(indn+k) + tmp1*res1
                dst_m(indn-k) = beta*dst_m(indn-k) + tmp1*res2
            end do
        end do
    end if
end subroutine fmm_m2l_ztranslate_adj_work

!> Direct M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of a multipole-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_m2l_rotation(c, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l)
    ! Inputs
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
        & beta
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Temporary workspace
    real(dp) :: work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8)
    ! Call corresponding work routine
    call fmm_m2l_rotation_work(c, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
end subroutine fmm_m2l_rotation

!> Direct M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of a multipole-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8 where p is a
!!      maximum of pm and pl
subroutine fmm_m2l_rotation_work(c, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
        & beta
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Temporary workspace
    real(dp), intent(out), target :: &
        & work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n, p
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_ml(:), tmp_ml2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        ! Workspace here is (pm+2)*(pm+1)
        call fmm_m2l_ztranslate_work(c(3), src_r, dst_r, pm, pl, vscales, &
            & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
        return
    end if
    ! Prepare pointers
    p = max(pm, pl)
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_ml(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_ml2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(pm, vcos, vsin, alpha, src_m, zero, tmp_ml)
    ! Perform rotation in the OXZ plane, work size is 4*pm*pm+13*pm+4
    call fmm_sph_rotate_oxz_work(pm, ctheta, -stheta, one, tmp_ml, zero, &
        & tmp_ml2, work)
    ! OZ translation, workspace here is (pm+2)*(pm+1)
    call fmm_m2l_ztranslate_work(rho, src_r, dst_r, pm, pl, vscales, &
        & m2l_ztranslate_coef, one, tmp_ml2, zero, tmp_ml, work)
    ! Backward rotation in the OXZ plane, work size is 4*pl*pl+13*pl+4
    call fmm_sph_rotate_oxz_work(pl, ctheta, stheta, one, tmp_ml, zero, &
        & tmp_ml2, work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(pl, vcos, vsin, one, tmp_ml2, beta, dst_l)
end subroutine fmm_m2l_rotation_work

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2l_bessel_rotation(c, src_r, dst_r, kappa, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_l)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta, kappa
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(6*p*p + 19*p + 8), src_sk(p+1), dst_si(p+1), s1, s2
    complex(dp) :: work_complex(2*p+1), z1, z2
    integer :: NZ, ierr
    ! Compute Bessel functions
    z1 = src_r * kappa
    z2 = dst_r * kappa
    s1 = sqrt(2/(pi*real(z1)))
    s2 = sqrt(pi/(2*real(z2)))
    call cbesk(z1, pt5, 1, 1, work_complex(1), NZ, ierr)
    src_sk(1) = s1 * real(work_complex(1))
    call cbesi(z2, pt5, 1, 1, work_complex(1), NZ, ierr)
    dst_si(1) = s2 * real(work_complex(1))
    if (p .gt. 0) then
        call cbesk(z1, 1.5d0, 1, p, work_complex(2:p+1), NZ, ierr)
        src_sk(2:p+1) = s1 * real(work_complex(2:p+1))
        call cbesi(z2, 1.5d0, 1, p, work_complex(2:p+1), NZ, ierr)
        dst_si(2:p+1) = s2 * real(work_complex(2:p+1))
    end if
    ! Call corresponding work routine
    call fmm_m2l_bessel_rotation_work(c*kappa, src_sk, dst_si, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_l, work, work_complex)
end subroutine fmm_m2l_bessel_rotation

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from old to new centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_m2l_bessel_rotation_work(c, src_sk, dst_si, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_l, work, work_complex)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: c(3), src_sk(p+1), dst_si(p+1), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(6*p*p + 19*p + 8)
    complex(dp), intent(out) :: work_complex(2*p+1)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_m(:), tmp_m2(:), vcos(:), vsin(:)
    ! Convert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    !if (stheta .eq. zero) then
    !    ! Workspace here is 2*(p+1)
    !    call fmm_m2m_bessel_ztranslate_work(c(3), src_r, dst_r, p, vscales, vcnk, &
    !        & alpha, src_m, beta, dst_m, work)
    !    return
    !end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_m(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_m2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_m, zero, tmp_m)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, one, tmp_m, zero, &
        & tmp_m2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_m2l_bessel_ztranslate_work(rho, src_sk, dst_si, p, vscales, vcnk, one, &
        & tmp_m2, zero, tmp_m, work, work_complex)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, one, tmp_m, zero, tmp_m2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, one, tmp_m2, beta, dst_l)
end subroutine fmm_m2l_bessel_rotation_work

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2l_bessel_ztranslate(z, src_sk, dst_si, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_sk(p+1), dst_si(p+1), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp) :: work(2*p+1)
    complex(dp) :: work_complex(2*p+1)
    ! Call corresponding work routine
    call fmm_m2l_bessel_ztranslate_work(z, src_sk, dst_si, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_l, work, work_complex)
end subroutine fmm_m2l_bessel_ztranslate

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from old to new centers of harmonics. Standard
!!      FMM M2M operation requires z to be OZ coordinate from new to old.
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_m2l_bessel_ztranslate_work(z, src_sk, dst_si, p, vscales, vcnk, &
        & alpha, src_m, beta, dst_l, work, work_complex)
    use complex_bessel
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: z, src_sk(p+1), dst_si(p+1), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(dp), intent(out), target :: work(2*p+1)
    complex(dp), intent(out) :: work_complex(2*p+1)
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, tmp3, res1, res2, pow_r1
    integer :: j, k, n, indj, indn, l, NZ, ierr
    ! Pointers for temporary values of powers
    real(dp) :: fact(2*p+1)
    complex(dp) :: complex_argument
    ! In case alpha is zero just do a proper scaling of output
    if (alpha .eq. zero) then
        if (beta .eq. zero) then
            dst_l = zero
        else
            dst_l = beta * dst_l
        end if
        return
    end if
    ! Get factorials
    fact(1) = one
    do j = 1, 2*p
        fact(j+1) = dble(j) * fact(j)
    end do
    ! Now alpha is non-zero
    ! If harmonics have different centers
    if (z .ne. 0) then
        complex_argument = z
        call cbesk(complex_argument, pt5, 1, 1, work_complex(1), NZ, ierr)
        work(1) = sqrt(2/(pi*z)) * real(work_complex(1))
        if (p .gt. 0) then
            call cbesk(complex_argument, 1.5d0, 1, 2*p, &
                & work_complex(2:2*p+1), NZ, ierr)
            work(2:2*p+1) = sqrt(2/(pi*z)) * real(work_complex(2:2*p+1))
        end if
        ! Do actual M2M
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            dst_l = zero
            do j = 0, p
                ! Offset for dst_m
                indj = j*j + j + 1
                ! k = 0
                k = 0
                res1 = zero
                do n = 0, p
                    ! Offset for src_m
                    indn = n*n + n + 1
                    tmp1 = vscales(indn) * ((-one)**(n)) * &
                        & dble(2*j+1) * fact(j+1) * fact(n+1) / &
                        & src_sk(n+1) / vscales(indj) * dst_si(j+1)
                    tmp3 = zero
                    do l = 0, min(j, n)
                        tmp2 = ((-two)**(-l)) / &
                            & fact(l+1)*fact(2*l+1)/ &
                            & fact(l+1)/fact(l+1)/fact(j-l+1)/fact(n-l+1)* &
                            & (work(j+n-l+1) / (z**l))
                        tmp3 = tmp3 + tmp2
                    end do
                    tmp3 = tmp3 * tmp1
                    res1 = res1 + tmp3*src_m(indn)
                end do
                dst_l(indj) = res1
                ! k != 0
                do k = 1, j
                    res1 = zero
                    res2 = zero
                    do n = k, p
                        ! Offset for src_m
                        indn = n*n + n + 1
                        tmp1 = vscales(indn+k) * ((-one)**(n)) * &
                            & dble(2*j+1) * fact(j-k+1) * fact(n+k+1) / &
                            & src_sk(n+1) / vscales(indj+k) * dst_si(j+1)
                        tmp3 = zero
                        do l = k, min(j, n)
                            tmp2 = ((-two)**(-l)) / &
                                & fact(l+k+1)*fact(2*l+1)/ &
                                & fact(l+1)/fact(l-k+1)/fact(j-l+1)/fact(n-l+1)* &
                                & (work(j+n-l+1) / (z**l))
                            tmp3 = tmp3 + tmp2
                        end do
                        tmp3 = tmp3 * tmp1
                        res1 = res1 + tmp3*src_m(indn+k)
                        res2 = res2 + tmp3*src_m(indn-k)
                    end do
                    dst_l(indj+k) = res1
                    dst_l(indj-k) = res2
                end do
            end do
        ! Update output if beta is non-zero
        else
            stop 111
        end if
    ! If harmonics are located at the same point
    else
        stop 222
        ! Overwrite output if beta is zero
        if (beta .eq. zero) then
            !r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_l(k) = tmp1 * src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        ! Update output if beta is non-zero
        else
            !r1 = src_r / dst_r
            tmp1 = alpha * r1
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_l(k) = beta*dst_l(k) + tmp1*src_m(k)
                end do
                tmp1 = tmp1 * r1
            end do
        end if
    end if
end subroutine fmm_m2l_bessel_ztranslate_work

!> Adjoint M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of a
!! multipole-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] m2l_ztranslate_adj_oef:
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: expansion in new harmonics
subroutine fmm_m2l_rotation_adj(c, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m)
    ! Inputs
    integer, intent(in) :: pl, pm
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), alpha, &
        & src_l((pl+1)*(pl+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Temporary workspace
    real(dp) :: work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8)
    ! Call corresponding work routine
    call fmm_m2l_rotation_adj_work(c, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
end subroutine fmm_m2l_rotation_adj

!> Adjoint M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of a
!! multipole-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] m2l_ztranslate_adj_coef:
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8 where p is a
!!      maximum of pm and pl
subroutine fmm_m2l_rotation_adj_work(c, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
    ! Inputs
    integer, intent(in) :: pl, pm
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & m2l_ztranslate_adj_coef(pl+1, pl+1, pm+1), alpha, &
        & src_l((pl+1)*(pl+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Temporary workspace
    real(dp), intent(out), target :: &
        & work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n, p
    ! Pointers for temporary values of harmonics
    real(dp), pointer :: tmp_ml(:), tmp_ml2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. 0) then
        ! Workspace here is (pl+2)*(pl+1)
        call fmm_m2l_ztranslate_adj_work(c(3), src_r, dst_r, pl, pm, vscales, &
            & m2l_ztranslate_adj_coef, alpha, src_l, beta, dst_m, work)
        return
    end if
    ! Prepare pointers
    p = max(pm, pl)
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_ml(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_ml2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(pl, vcos, vsin, alpha, src_l, zero, tmp_ml)
    ! Perform rotation in the OXZ plane, work size is 4*pl*pl+13*pl+4
    call fmm_sph_rotate_oxz_work(pl, ctheta, -stheta, one, tmp_ml, zero, &
        & tmp_ml2, work)
    ! OZ translation, workspace here is (pl+2)*(pl+1)
    call fmm_m2l_ztranslate_adj_work(rho, src_r, dst_r, pl, pm, vscales, &
        & m2l_ztranslate_adj_coef, one, tmp_ml2, zero, tmp_ml, work)
    ! Backward rotation in the OXZ plane, work size is 4*pm*pm+13*pm+4
    call fmm_sph_rotate_oxz_work(pm, ctheta, stheta, one, tmp_ml, zero, &
        & tmp_ml2, work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(pm, vcos, vsin, one, tmp_ml2, beta, dst_m)
end subroutine fmm_m2l_rotation_adj_work

end module ddx_harmonics

