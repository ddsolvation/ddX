module pcm_fmm
use ddcosmo, only: inl, nl, lmax, zero, one, two, pi, basis, se, eta, dfsw, pt5, zi
implicit none
real(kind=8) :: sqrt_2=0, sqrt_four_pi=0
real(kind=8) :: total_time_matvec=0
integer :: total_count_matvec=0
integer :: pmax
integer, parameter :: ifmm_onfly=0
contains

! Init global constants
subroutine init_globals(pm, pl, vscales, ngrid, w, grid, vgrid)
! Parameters:
!   pm: maximum degree of multipole polynomials to compute
!   pl: maximum degree of local polynomials to compute
!   vscales: values of scaling factors for spherical harmonics
!   ngrid: number of Lebedev grid points
!   w: weights of Lebedev grid points
!   grid: coordinates of Lebedev grid points on unit sphere
!   vgrid: values of spherical harmonics at grid points
    integer, intent(in) :: pm, pl, ngrid
    real(kind=8), intent(out) :: vscales((pm+pl+1)*(pm+pl+1)), w(ngrid)
    real(kind=8), intent(out) :: grid(3, ngrid), vgrid((pl+1)*(pl+1), ngrid)
    integer :: i, n, m, indn, indm
    real(kind=8) :: c(3), ctheta, stheta, cphi, sphi, vplm((pl+1)*(pl+1))
    real(kind=8) :: vcos(pl+1), vsin(pl+1), tmp
    sqrt_2 = sqrt(dble(2))
    sqrt_four_pi = 4 * sqrt(atan(dble(1)))
    call scales_real_normal(pl+pm, vscales)
    call llgrid(ngrid, w, grid)
    do i = 1, ngrid
        c = grid(:, i)
        ctheta = c(3)
        stheta = sqrt(c(1)*c(1) + c(2)*c(2))
        if (stheta .ne. 0) then
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            call trgev(cphi, sphi, pl, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg(ctheta, stheta, pl, vplm)
        do n = 0, pl
            indn = n*n + n + 1
            vgrid(indn, i) = vscales(indn) * vplm(indn)
            do m = 1, n
                indm = indn + m
                tmp = vscales(indm) * vplm(indm)
                vgrid(indn+m, i) = tmp * vcos(1+m)
                vgrid(indn-m, i) = tmp * vsin(1+m)
            end do
        end do
    end do
end subroutine init_globals

! Compute associated Legendre polynomials P(l,m)(x)
! Uses recurrence formula
!   (l-m)P(l,m) = x(2l-1)P(l-1,m) - (l+m-1)P(l-2,m)
! Input x must be in range [-1;1]
! This function is simply copied from ddCOSMO
subroutine polleg(x, y, p, vplm)
! Parameters:
!   x: coordinate in interval [-1;1]
!   y: sqrt(1-x*2)
!   p: maximum degree of polynomials to compute
!   vplm: values of associated legendre polynomials P(l,m)(x)
    real(kind=8), intent(in) :: x, y
    integer, intent(in) :: p
    real(kind=8), intent(out) :: vplm((p+1)*(p+1))
    integer :: m, ind, l, ind2
    real(kind=8) :: fact, pmm, pmm1, pmmo, pll, fm, fl
    fact = 1
    pmm = 1
    do m = 0, p
        ind = (m + 1)*(m + 1)
        vplm(ind) = pmm
        if (m .eq. p) then
            return
        end if
        fm = dble(m)
        pmm1 = x * (2*fm+1) * pmm
        ind2 = ind + 2*m + 2
        vplm(ind2) = pmm1
        pmmo = pmm
        do l = m+2, p
            fl = dble(l)
            pll = x*(2*fl-1)*pmm1 - (fl+fm-1)*pmm
            pll = pll / (fl-fm)
            ind = l*l + l + 1
            vplm(ind+m) = pll
            pmm = pmm1
            pmm1 = pll
        end do
        pmm = -pmmo * fact * y
        fact = fact + 2
    end do
end subroutine polleg

! Compute scaling factors for normalized real spherical harmonics
subroutine scales_real_normal(p, vscales)
! Parameters:
!   p: maximum degree of spherical harmonics
!   vscales: values of scaling factors
    integer, intent(in) :: p
    real(kind=8), intent(out) :: vscales((p+1)*(p+1))
    real(kind=8) :: tmp
    integer :: l, ind, m
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        vscales(ind) = sqrt(dble(2*l+1)) / sqrt_four_pi
        tmp = vscales(ind) * sqrt_2
        ! m != 0
        do m = 1, l
            tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
            vscales(ind+m) = tmp
        end do
    end do
end subroutine scales_real_normal

! Compute multipole coefficients for particle of unit charge
! Based on normalized real spherical harmonics Y_l^m, scaled by r^(l+1). It
! means corresponding coefficients are scaled by r^(-l-1).
! This function is not needed for pcm, but it is useful for testing purposes
subroutine fmm_p2m(c, r, p, vscales, m)
! Parameters:
!   c: coordinates of charged particle (relative to center of harmonics)
!   r: radius of spherical harmonics
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   m: multipole coefficients
    real(kind=8), intent(in) :: c(3), r, vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: m((p+1)*(p+1))
    real(kind=8) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), t, tmp, rcoef
    integer :: n, k, ind
    stheta = c(1)*c(1) + c(2)*c(2)
    rho = sqrt(c(3)*c(3) + stheta)
    if (rho .ne. 0) then
        ctheta = c(3) / rho
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / rho
            call trgev(cphi, sphi, p, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg(ctheta, stheta, p, vplm)
        ! Now build harmonics to fill multipole coefficients
        rcoef = rho / r
        t = 1 / r
        do n = 0, p
            ind = n*n + n + 1
            m(ind) = t * vscales(ind) * vplm(ind)
            do k = 1, n
                tmp = t * vscales(ind+k) * vplm(ind+k)
                m(ind+k) = tmp * vcos(k+1)
                m(ind-k) = tmp * vsin(k+1)
            end do
            t = t * rcoef
        end do
    else
        m(1) = 1 / r
        m(2:) = 0
    end if
end subroutine fmm_p2m

! Compute potential, induced by multipole spherical harmonics
! Based on normalized scaled real spherical harmonics r^(l+1) Y_l^m
subroutine fmm_m2p(c, r, p, vscales, m, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   r: radius of spherical harmonics
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   m: multipole expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3), r, vscales((p+1)*(p+1)), m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    real(kind=8) :: tmp, tmp1, tmp2
    real(kind=8) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), rcoef
    integer :: n, k, ind
    stheta = c(1)*c(1) + c(2)*c(2)
    rho = sqrt(c(3)*c(3) + stheta)
    v = 0
    if (rho .eq. 0) then
        return
    end if
    ! rho is always > 0
    ctheta = c(3) / rho
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / rho
        call trgev(cphi, sphi, p, vcos, vsin)
    else
        cphi = 1
        sphi = 0
        vcos = 1
        vsin = 0
    end if
    call polleg(ctheta, stheta, p, vplm)
    rcoef = r / rho
    t = 1
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        ! k = 0
        tmp = m(ind) * vplm(ind) * vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = vplm(ind+k) * vscales(ind+k)
            tmp2 = m(ind+k)*vcos(k+1) + m(ind-k)*vsin(k+1)
            tmp = tmp + tmp1*tmp2
        end do
        v = v + t*tmp/vscales(ind)**2
    end do
end subroutine fmm_m2p

! Compute matrix of M2P, induced by multipole spherical harmonics
! Based on normalized scaled real spherical harmonics r^(l+1) Y_l^m
subroutine fmm_m2p_mat(c, r, p, vscales, mat)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   r: radius of spherical harmonics
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   ui: scalar multiplier for matrix
!   mat: transfer matrix from multipole expansion to point potential
    real(kind=8), intent(in) :: c(3), r, vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: mat((p+1)*(p+1))
    real(kind=8) :: tmp1, tmp2
    real(kind=8) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), rcoef
    integer :: n, k, ind
    stheta = c(1)*c(1) + c(2)*c(2)
    rho = sqrt(c(3)*c(3) + stheta)
    if (rho .eq. 0) then
        return
    end if
    ! rho is always > 0
    ctheta = c(3) / rho
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / rho
        call trgev(cphi, sphi, p, vcos, vsin)
    else
        cphi = 1
        sphi = 0
        vcos = 1
        vsin = 0
    end if
    call polleg(ctheta, stheta, p, vplm)
    rcoef = r / rho
    t = 1
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        tmp1 = t / vscales(ind)**2
        ! k = 0
        mat(ind) = tmp1 * vplm(ind) * vscales(ind)
        ! k != 0
        do k = 1, n
            tmp2 = tmp1 * vplm(ind+k) * vscales(ind+k)
            mat(ind+k) = tmp2 * vcos(k+1)
            mat(ind-k) = tmp2 * vsin(k+1)
        end do
    end do
end subroutine fmm_m2p_mat

! M2M baseline translation (p^4 operations)
! Baseline in terms of operation count: p^4
subroutine fmm_m2m_baseline(c, src_r, dst_r, p, vscales, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), fact(2*p+1), tmpk1, tmpk2, tmpk3, tmp1
    real(kind=8) :: tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, m, indj, indm, indn, indjn
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg(ctheta, stheta, p, vplm)
        r1 = src_r / dst_r
        r2 = r / dst_r
        pow_r1(1) = r1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Fill square roots of factorials
        fact(1) = 1
        do j = 2, 2*p+1
            fact(j) = sqrt(dble(j-1)) * fact(j-1)
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = vscales(indj) * fact(j-k+1) * fact(j+k+1)
                if (k .ne. 0) then
                    tmp1 = tmp1 * sqrt_2
                end if
                do n = 0, j
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp2 = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                        & vscales(indjn) / vscales(indn)
                    do m = max(k+n-j, -n), min(k+j-n, n)
                        indm = indn + abs(m)
                        cphi = vcos(1+abs(m))
                        sphi = vsin(1+abs(m))
                        tmpk1 = tmp2 / fact(n-m+1) / fact(n+m+1) / &
                            & fact(j-n-k+m+1) / fact(j-n+k-m+1) * &
                            & vplm(indm) * vscales(indm)
                        if (mod(abs(k-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        tmpk2 = src_m(indjn+abs(k-m)) * cphi
                        if ((m .ge. 0) .and. (m .le. k)) then
                            sphi = -sphi
                        end if
                        tmpk3 = -src_m(indjn+abs(k-m)) * sphi
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                            tmpk2 = tmpk2 + src_m(indjn-abs(k-m))*sphi
                            tmpk3 = tmpk3 + src_m(indjn-abs(k-m))*cphi
                        end if
                        if (m .gt. k) then
                            tmpk3 = -tmpk3
                        end if
                        dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            dst_m(indj-k) = dst_m(indj-k) + tmpk1*tmpk3
                        end if
                    end do
                end do
            end do
        end do
    else
        r1 = src_r / dst_r
        tmpk1 = r1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_m(k) = dst_m(k) + src_m(k)*tmpk1
            end do
            tmpk1 = tmpk1 * r1
        end do
    end if
end subroutine fmm_m2m_baseline

! Compute potential, induced by local spherical harmonics
! Based on normalized scaled real spherical harmonics of given radius
subroutine fmm_l2p(c, r, p, vscales, l, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   r: radius of spherical harmonics
!   p: maximum degree of local basis functions
!   vscales: normalization constants for Y_lm
!   l: local expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3), r, vscales((p+1)*(p+1)), l((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    real(kind=8) :: tmp, tmp1, tmp2
    real(kind=8) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), rcoef
    integer :: n, k, ind
    t = 1
    stheta = c(1)*c(1) + c(2)*c(2)
    rho = sqrt(c(3)*c(3) + stheta)
    if (rho .eq. 0) then
        ! Only first term
        v = t * l(1) / vscales(1)
        return
    end if
    v = 0
    ctheta = c(3) / rho
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / rho
        call trgev(cphi, sphi, p, vcos, vsin)
    else
        cphi = 1
        sphi = 0
        vcos = 1
        vsin = 0
    end if
    call polleg(ctheta, stheta, p, vplm)
    rcoef = rho / r
    do n = 0, p
        ind = n*n + n + 1
        ! k = 0
        tmp = l(ind) * vplm(ind) * vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = vplm(ind+k) * vscales(ind+k)
            tmp2 = l(ind+k)*vcos(k+1) + l(ind-k)*vsin(k+1)
            tmp = tmp + tmp1*tmp2
        end do
        v = v + t*tmp/vscales(ind)**2
        t = t * rcoef
    end do
end subroutine fmm_l2p

! Compute local expansion by given multipole expansion
! Baseline in terms of operation count: p^4
subroutine fmm_m2l_baseline(c, src_r, dst_r, pm, pl, vscales, src_m, dst_l)
! Parameters:
!   c: radius-vector from new (local) to old (multipole) centers of harmonics
!   src_r: radius of old (multipole) harmonics
!   dst_r: radius of new (local) harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm (of degree up to pl+pm)
!   src_m: expansion in old (multipole) harmonics
!   dst_l: expansion in new (local) harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: src_m((pm+1)*(pm+1))
    integer, intent(in) :: pm, pl
    real(kind=8), intent(inout) :: dst_l((pl+1)*(pl+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi
    real(kind=8) :: vcos(pm+pl+1), vsin(pm+pl+1)
    real(kind=8) :: vplm((pm+pl+1)*(pm+pl+1)), fact(2*(pm+pl)+1), tmpk1, tmpk2
    real(kind=8) :: tmpk3, tmp1, tmp2, pow_r1(pm+1), pow_r2(pl+1)
    integer :: j, k, n, m, indj, indmk, indn, indjn
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    ! r cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (r .eq. 0) then
        return
    end if
    ctheta = c(3) / r
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / r
        call trgev(cphi, sphi, pm+pl, vcos, vsin)
    else
        cphi = 1
        sphi = 0
        vcos = 1
        vsin = 0
    end if
    call polleg(ctheta, stheta, pm+pl, vplm)
    r1 = src_r / r
    r2 = dst_r / r
    pow_r1(1) = r1
    pow_r2(1) = 1
    do j = 2, pm+1
        pow_r1(j) = pow_r1(j-1) * r1
    end do
    do j = 2, pl+1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    ! Fill square roots of factorials
    fact(1) = 1
    do j = 2, 2*(pm+pl)+1
        fact(j) = sqrt(dble(j-1)) * fact(j-1)
    end do
    do j = 0, pl
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * pow_r2(j+1) / fact(j-k+1) / fact(j+k+1)
            if (k .ne. 0) then
                tmp1 = tmp1 * sqrt_2
            end if
            do n = 0, pm
                indn = n*n + n + 1
                indjn = (j+n)**2 + (j+n) + 1
                tmp2 = tmp1 * pow_r1(n+1) / vscales(indjn) / vscales(indn)
                if (mod(n, 2) .eq. 1) then
                    tmp2 = -tmp2
                end if
                do m = -n, n
                    indmk = indjn + abs(m-k)
                    cphi = vcos(1+abs(m-k))
                    sphi = vsin(1+abs(m-k))
                    tmpk1 = tmp2 / fact(n-m+1) / fact(n+m+1) * &
                        & fact(j+n-m+k+1) * fact(j+n+m-k+1) * vplm(indmk) * &
                        & vscales(indmk)
                    if (mod(abs(k+abs(m)-abs(k-m)), 4) .eq. 2) then
                        tmpk1 = -tmpk1
                    end if
                    tmpk2 = src_m(indn+abs(m)) * cphi
                    if ((m .ge. 0) .and. (m .le. k)) then
                        sphi = -sphi
                    end if
                    tmpk3 = -src_m(indn+abs(m)) * sphi
                    if (m .ne. k) then
                        tmpk1 = tmpk1 / sqrt_2
                    end if
                    if (m .ne. 0) then
                        tmpk1 = tmpk1 / sqrt_2
                        tmpk2 = tmpk2 + src_m(indn-abs(m))*sphi
                        tmpk3 = tmpk3 + src_m(indn-abs(m))*cphi
                    end if
                    if (m .lt. 0) then
                        tmpk3 = -tmpk3
                    end if
                    dst_l(indj+k) = dst_l(indj+k) + tmpk1*tmpk2
                    if (k .ne. 0) then
                        dst_l(indj-k) = dst_l(indj-k) + tmpk1*tmpk3
                    end if
                end do
            end do
        end do
    end do
end subroutine fmm_m2l_baseline

! Translate local expansion to another sphere
! Baseline in terms of operation count: p^4
subroutine fmm_l2l_baseline(c, src_r, dst_r, p, vscales, src_l, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm (of degree up to p)
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi
    real(kind=8) :: vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), fact(2*p+1), tmpk1, tmpk2, tmpk3
    real(kind=8) :: tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, m, indj, indmk, indn, indjn
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg(ctheta, stheta, p, vplm)
        r1 = r / src_r
        r2 = dst_r / r
        pow_r1(1) = 1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Fill square roots of factorials
        fact(1) = 1
        do j = 2, 2*p+1
            fact(j) = sqrt(dble(j-1)) * fact(j-1)
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = pow_r2(j+1) / fact(j-k+1) / fact(j+k+1) * vscales(indj)
                if (k .ne. 0) then
                    tmp1 = tmp1 * sqrt_2
                end if
                do n = j, p
                    indn = n*n + n + 1
                    indjn = (n-j)**2 + (n-j) + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indjn) / vscales(indn)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    do m = k+j-n, k+n-j
                        indmk = indjn + abs(m-k)
                        cphi = vcos(1+abs(m-k))
                        sphi = vsin(1+abs(m-k))
                        tmpk1 = tmp2 * fact(n-m+1) * fact(n+m+1) / &
                            & fact(n-j-m+k+1) / fact(n-j+m-k+1) * &
                            & vplm(indmk) * vscales(indmk)
                        if (mod(abs(k+abs(m-k)-abs(m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        tmpk2 = src_l(indn+abs(m)) * cphi
                        if ((m .ge. 0) .and. (m .le. k)) then
                            sphi = -sphi
                        end if
                        tmpk3 = -src_l(indn+abs(m)) * sphi
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                            tmpk2 = tmpk2 + src_l(indn-abs(m))*sphi
                            tmpk3 = tmpk3 + src_l(indn-abs(m))*cphi
                        end if
                        if (m .lt. 0) then
                            tmpk3 = -tmpk3
                        end if
                        dst_l(indj+k) = dst_l(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            dst_l(indj-k) = dst_l(indj-k) + tmpk1*tmpk3
                        end if
                    end do
                end do
            end do
        end do
    else
        r1 = dst_r / src_r
        tmpk1 = 1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmpk1
            end do
            tmpk1 = tmpk1 * r1
        end do
    end if
end subroutine fmm_l2l_baseline

! Optimized version of M2M, translation over OZ axis only
subroutine fmm_m2m_ztranslate(z, src_r, dst_r, p, vscales, src_m, dst_m)
! Parameters:
!   z: radius-vector from new to old centers of harmonics (z coordinate only)
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: r1, r2, fact(2*p+1), tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn, indjn
    if (z .ne. 0) then
        r1 = src_r / dst_r
        r2 = z / dst_r
        pow_r1(1) = r1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Fill square roots of factorials
        fact(1) = 1
        do j = 2, 2*p+1
            fact(j) = sqrt(dble(j-1)) * fact(j-1)
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = vscales(indj) * fact(j-k+1) * fact(j+k+1)
                do n = 0, j-k
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp2 = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                        & vscales(indjn) / fact(n+1) / fact(n+1) / &
                        & fact(j-n-k+1) / fact(j-n+k+1)
                    if (k .eq. 0) then
                        dst_m(indj) = dst_m(indj) + tmp2*src_m(indjn)
                    else
                        dst_m(indj+k) = dst_m(indj+k) + tmp2*src_m(indjn+k)
                        dst_m(indj-k) = dst_m(indj-k) + tmp2*src_m(indjn-k)
                    end if
                end do
            end do
        end do
    else
        r1 = src_r / dst_r
        tmp1 = r1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_m(k) = dst_m(k) + src_m(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_m2m_ztranslate

! Scale M2M, when spherical harmonics are centered in the same point
subroutine fmm_m2m_scale(src_r, dst_r, p, src_m, dst_m)
! Parameters:
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: r1, tmp1
    integer :: j, k, indj
    r1 = src_r / dst_r
    tmp1 = r1
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_m(k) = dst_m(k) + src_m(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_m2m_scale

! Compute and store translation over OZ axis in a special sparse matrix
subroutine fmm_m2m_get_ztrans_mat(z, src_r, dst_r, p, vscales, mat)
! Parameters:
!   z: radius-vector from new to old centers of harmonics (z coordinate only).
!       Value of z must be non-zero. In case z=0 call ordinary fmm_m2m_scale,
!       since this case is simply scaling, no need to compute and save matrix.
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   mat: translation matrix for spherical harmonics
    real(kind=8), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8) :: r1, r2, fact(2*p+1), tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn, indjn, indmat
    if (z .eq. 0) then
        return
    end if
    r1 = src_r / dst_r
    r2 = z / dst_r
    pow_r1(1) = r1
    pow_r2(1) = 1
    do j = 2, p+1
        pow_r1(j) = pow_r1(j-1) * r1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    ! Fill square roots of factorials
    fact(1) = 1
    do j = 2, 2*p+1
        fact(j) = sqrt(dble(j-1)) * fact(j-1)
    end do
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * fact(j-k+1) * fact(j+k+1)
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                if(indmat .gt. (p+1)*(p+2)*(p+3)/6) then
                    write(*,*) "wrong max index in fmm_m2m_get_ztrans_mat"
                    exit
                end if
                mat(indmat) = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                    & vscales(indjn) / fact(n+1) / fact(n+1) / &
                    & fact(j-n-k+1) / fact(j-n+k+1)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2m_get_ztrans_mat

! Use precomputed translation matrices over OZ axis
! Check if centers of spherical harmonics are different, otherwise use function
! fmm_m2m_scale, since case of the same center requires only scaling
subroutine fmm_m2m_use_ztrans_mat(p, mat, src_m, dst_m)
! Parameters:
!   p: maximum degree of spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    real(kind=8), intent(out) :: dst_m((p+1)*(p+1))
    integer :: j, k, n, indj, indn, indjn, indmat
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        dst_m(indj) = 0
        ! k = 0
        do n = 0, j
            indn = n*n + n + 1
            indjn = (j-n)**2 + (j-n) + 1
            dst_m(indj) = dst_m(indj) + mat(indmat)*src_m(indjn)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            dst_m(indj+k) = 0
            dst_m(indj-k) = 0
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                dst_m(indj+k) = dst_m(indj+k) + mat(indmat)*src_m(indjn+k)
                dst_m(indj-k) = dst_m(indj-k) + mat(indmat)*src_m(indjn-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2m_use_ztrans_mat

! Use precomputed matrices to perform adjoint translation over OZ axis
! Check if centers of spherical harmonics are different, otherwise use function
! fmm_m2m_scale, since case of the same center requires only scaling
subroutine fmm_m2m_adj_use_ztrans_mat(p, mat, src_m, dst_m)
! Parameters:
!   p: maximum degree of spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    real(kind=8), intent(out) :: dst_m((p+1)*(p+1))
    integer :: j, k, n, indj, indn, indjn, indmat
    indmat = 1
    dst_m = 0
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = 0, j
            indn = n*n + n + 1
            indjn = (j-n)**2 + (j-n) + 1
            dst_m(indjn) = dst_m(indjn) + mat(indmat)*src_m(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                dst_m(indjn+k) = dst_m(indjn+k) + mat(indmat)*src_m(indj+k)
                dst_m(indjn-k) = dst_m(indjn-k) + mat(indmat)*src_m(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2m_adj_use_ztrans_mat

! Use precomputed translation matrices over OZ axis
! Check if centers of spherical harmonics are different, otherwise use function
! fmm_m2m_scale, since case of the same center requires only scaling.
! Only difference with fmm_m2m_use_ztrans_mat is that output dst_m is treated
! as inout, not just out. This is done to improve performance.
subroutine fmm_m2m_use_ztrans_mat2(p, mat, src_m, dst_m)
! Parameters:
!   p: maximum degree of spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    integer :: j, k, n, indj, indn, indjn, indmat
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = 0, j
            indn = n*n + n + 1
            indjn = (j-n)**2 + (j-n) + 1
            dst_m(indj) = dst_m(indj) + mat(indmat)*src_m(indjn)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                dst_m(indj+k) = dst_m(indj+k) + mat(indmat)*src_m(indjn+k)
                dst_m(indj-k) = dst_m(indj-k) + mat(indmat)*src_m(indjn-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2m_use_ztrans_mat2

! Use precomputed matrices to perform adjoint translation over OZ axis
! Check if centers of spherical harmonics are different, otherwise use function
! fmm_m2m_scale, since case of the same center requires only scaling.
! Only difference with fmm_m2m_use_ztrans_mat is that output dst_m is treated
! as inout, not just out. This is done to improve performance.
subroutine fmm_m2m_adj_use_ztrans_mat2(p, mat, src_m, dst_m)
! Parameters:
!   p: maximum degree of spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    integer :: j, k, n, indj, indn, indjn, indmat
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = 0, j
            indn = n*n + n + 1
            indjn = (j-n)**2 + (j-n) + 1
            dst_m(indjn) = dst_m(indjn) + mat(indmat)*src_m(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                dst_m(indjn+k) = dst_m(indjn+k) + mat(indmat)*src_m(indj+k)
                dst_m(indjn-k) = dst_m(indjn-k) + mat(indmat)*src_m(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2m_adj_use_ztrans_mat2

! Optimized version of M2L, translation over OZ axis only
subroutine fmm_m2l_ztranslate(z, src_r, dst_r, pm, pl, vscales, src_m, dst_l)
! Parameters:
!   z: radius-vector from new to old centers of harmonics (z coordinate only)
!   src_r: radius of old (multipole) harmonics
!   dst_r: radius of new (local) harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm (of degree up to pl+pm)
!   src_m: expansion in old (multipole) harmonics
!   dst_l: expansion in new (local) harmonics
    real(kind=8), intent(in) :: z, src_r, dst_r
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: src_m((pm+1)*(pm+1))
    integer, intent(in) :: pm, pl
    real(kind=8), intent(inout) :: dst_l((pl+1)*(pl+1))
    real(kind=8) :: fact(2*(pm+pl)+1), tmp1, tmp2, pow_r1(pm+1), pow_r2(pl+1)
    real(kind=8) :: r1, r2
    integer :: j, k, n, indj, indn
    ! z cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (z .eq. 0) then
        return
    end if
    r1 = src_r / z
    r2 = dst_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    pow_r1(1) = abs(r1)
    pow_r2(1) = 1
    do j = 2, pm+1
        pow_r1(j) = pow_r1(j-1) * r1
    end do
    do j = 2, pl+1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    ! Fill square roots of factorials
    fact(1) = 1
    do j = 2, 2*(pm+pl)+1
        fact(j) = sqrt(dble(j-1)) * fact(j-1)
    end do
    do j = 0, pl
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * pow_r2(j+1) / fact(j-k+1) / fact(j+k+1)
            do n = k, pm
                indn = n*n + n + 1
                tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) / fact(n-k+1) / &
                    & fact(n+k+1) * fact(j+n+1) * fact(j+n+1)
                if (mod(n+k, 2) .eq. 1) then
                    tmp2 = -tmp2
                end if
                if (k .eq. 0) then
                    dst_l(indj) = dst_l(indj) + tmp2*src_m(indn)
                else
                    dst_l(indj+k) = dst_l(indj+k) + tmp2*src_m(indn+k)
                    dst_l(indj-k) = dst_l(indj-k) + tmp2*src_m(indn-k)
                end if
            end do
        end do
    end do
end subroutine fmm_m2l_ztranslate

! Compute and store translation over OZ axis in a special sparse matrix
subroutine fmm_m2l_get_ztrans_mat(z, src_r, dst_r, pm, pl, vscales, mat)
! Parameters:
!   z: radius-vector from new to old centers of harmonics (z coordinate only).
!       Value of z must be non-zero, since spheres with multipole and local
!       harmonics must be separated in space.
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm
!   mat: translation matrix for spherical harmonics
    real(kind=8), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1))
    integer, intent(in) :: pm, pl
    real(kind=8), intent(out) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    real(kind=8) :: r1, r2, fact(2*(pm+pl)+1), tmp1, pow_r1(pm+1), pow_r2(pl+1)
    integer :: j, k, n, indj, indn, indjn, indmat
    if (z .eq. 0) then
        return
    end if
    r1 = src_r / z
    r2 = dst_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    pow_r1(1) = abs(r1)
    pow_r2(1) = 1
    do j = 2, pm+1
        pow_r1(j) = pow_r1(j-1) * r1
    end do
    do j = 2, pl+1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    ! Fill square roots of factorials
    fact(1) = 1
    do j = 2, 2*(pm+pl)+1
        fact(j) = sqrt(dble(j-1)) * fact(j-1)
    end do
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * pow_r2(j+1) / fact(j-k+1) / fact(j+k+1)
            do n = k, pm
                indn = n*n + n + 1
                if(indmat .gt. (min(pm,pl)+1)*(min(pm,pl)+2) &
                    & *(3*max(pm,pl)+3-min(pm,pl))/6) then
                    write(*,*) "wrong max index in fmm_m2l_get_ztrans_mat"
                    exit
                end if
                mat(indmat) = tmp1 * pow_r1(n+1) / vscales(indn) &
                    & / fact(n-k+1) / fact(n+k+1) * fact(j+n+1) * fact(j+n+1)
                if (mod(n+k, 2) .eq. 1) then
                    mat(indmat) = -mat(indmat)
                end if
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2l_get_ztrans_mat

! Use precomputed translation matrices over OZ axis
subroutine fmm_m2l_use_ztrans_mat(pm, pl, mat, src_m, dst_l)
! Parameters:
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_m: expansion in old (multipole) harmonics
!   dst_l: expansion in new (local) harmonics
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    real(kind=8), intent(in) :: src_m((pm+1)*(pm+1))
    real(kind=8), intent(out) :: dst_l((pl+1)*(pl+1))
    integer :: j, k, n, indj, indn, indjn, indmat
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        dst_l(indj) = 0
        ! k = 0
        do n = 0, pm
            indn = n*n + n + 1
            dst_l(indj) = dst_l(indj) + mat(indmat)*src_m(indn)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            dst_l(indj+k) = 0
            dst_l(indj-k) = 0
            do n = k, pm
                indn = n*n + n + 1
                dst_l(indj+k) = dst_l(indj+k) + mat(indmat)*src_m(indn+k)
                dst_l(indj-k) = dst_l(indj-k) + mat(indmat)*src_m(indn-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2l_use_ztrans_mat

! Use precomputed matrices to perform adjoint translation over OZ axis
subroutine fmm_m2l_adj_use_ztrans_mat(pm, pl, mat, src_l, dst_m)
! Parameters:
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_l: expansion in old (local) harmonics
!   dst_m: expansion in new (multipole) harmonics
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    real(kind=8), intent(in) :: src_l((pl+1)*(pl+1))
    real(kind=8), intent(out) :: dst_m((pm+1)*(pm+1))
    integer :: j, k, n, indj, indn, indjn, indmat
    indmat = 1
    dst_m = 0
    do j = 0, pl
        indj = j*j + j + 1
        ! k = 0
        do n = 0, pm
            indn = n*n + n + 1
            dst_m(indn) = dst_m(indn) + mat(indmat)*src_l(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = k, pm
                indn = n*n + n + 1
                dst_m(indn+k) = dst_m(indn+k) + mat(indmat)*src_l(indj+k)
                dst_m(indn-k) = dst_m(indn-k) + mat(indmat)*src_l(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2l_adj_use_ztrans_mat

! Use precomputed translation matrices over OZ axis
! Only difference with fmm_m2l_use_ztrans_mat is that output dst_l is treated
! as inout, not just out. This is done to improve performance.
subroutine fmm_m2l_use_ztrans_mat2(pm, pl, mat, src_m, dst_l)
! Parameters:
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_m: expansion in old (multipole) harmonics
!   dst_l: expansion in new (local) harmonics
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    real(kind=8), intent(in) :: src_m((pm+1)*(pm+1))
    real(kind=8), intent(inout) :: dst_l((pl+1)*(pl+1))
    integer :: j, k, n, indj, indn, indjn, indmat
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        ! k = 0
        do n = 0, pm
            indn = n*n + n + 1
            dst_l(indj) = dst_l(indj) + mat(indmat)*src_m(indn)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = k, pm
                indn = n*n + n + 1
                dst_l(indj+k) = dst_l(indj+k) + mat(indmat)*src_m(indn+k)
                dst_l(indj-k) = dst_l(indj-k) + mat(indmat)*src_m(indn-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2l_use_ztrans_mat2

! Use precomputed matrices to perform adjoint translation over OZ axis
! Only difference with fmm_m2l_use_ztrans_mat is that output dst_m is treated
! as inout, not just out. This is done to improve performance.
subroutine fmm_m2l_adj_use_ztrans_mat2(pm, pl, mat, src_l, dst_m)
! Parameters:
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_l: expansion in old (local) harmonics
!   dst_m: expansion in new (multipole) harmonics
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    real(kind=8), intent(in) :: src_l((pl+1)*(pl+1))
    real(kind=8), intent(inout) :: dst_m((pm+1)*(pm+1))
    integer :: j, k, n, indj, indn, indjn, indmat
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        ! k = 0
        do n = 0, pm
            indn = n*n + n + 1
            dst_m(indn) = dst_m(indn) + mat(indmat)*src_l(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = k, pm
                indn = n*n + n + 1
                dst_m(indn+k) = dst_m(indn+k) + mat(indmat)*src_l(indj+k)
                dst_m(indn-k) = dst_m(indn-k) + mat(indmat)*src_l(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2l_adj_use_ztrans_mat2

! Optimized version of L2L, translation over OZ axis only
subroutine fmm_l2l_ztranslate(z, src_r, dst_r, p, vscales, src_l, dst_l)
! Parameters:
!   z: radius-vector from new to old centers of harmonics (z coordinate only)
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm (of degree up to p)
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    real(kind=8), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    real(kind=8) :: r1, r2, fact(2*p+1)
    real(kind=8) :: tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn
    if (z .ne. 0) then
        r1 = z / src_r
        r2 = dst_r / z
        pow_r1(1) = 1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Fill square roots of factorials
        fact(1) = 1
        do j = 2, 2*p+1
            fact(j) = sqrt(dble(j-1)) * fact(j-1)
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = pow_r2(j+1) / fact(j-k+1) / fact(j+k+1) * vscales(indj)
                do n = j, p
                    indn = n*n + n + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) * fact(n-k+1) * &
                        & fact(n+k+1) / fact(n-j+1) / fact(n-j+1)
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
    else
        r1 = dst_r / src_r
        tmp1 = 1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_l2l_ztranslate

! Scale L2L, when spherical harmonics are centered in the same point
subroutine fmm_l2l_scale(src_r, dst_r, p, src_l, dst_l)
! Parameters:
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    real(kind=8), intent(in) :: src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    real(kind=8) :: r1, tmp1
    integer :: j, k, indj
    r1 = dst_r / src_r
    tmp1 = 1
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_l(k) = dst_l(k) + src_l(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_l2l_scale

! Compute and store translation over OZ axis in a special sparse matrix
subroutine fmm_l2l_get_ztrans_mat(z, src_r, dst_r, p, vscales, mat)
! Parameters:
!   z: radius-vector from new to old centers of harmonics (z coordinate only).
!       Value of z must be non-zero. In case z=0 call ordinary fmm_l2l_scale,
!       since this case is simply scaling, no need to compute and save matrix.
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   mat: translation matrix for spherical harmonics
    real(kind=8), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8) :: r1, r2, fact(2*p+1), tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn, indmat
    if (z .eq. 0) then
        return
    end if
    r1 = z / src_r
    r2 = dst_r / z
    pow_r1(1) = 1
    pow_r2(1) = 1
    do j = 2, p+1
        pow_r1(j) = pow_r1(j-1) * r1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    ! Fill square roots of factorials
    fact(1) = 1
    do j = 2, 2*p+1
        fact(j) = sqrt(dble(j-1)) * fact(j-1)
    end do
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = pow_r2(j+1) / fact(j-k+1) / fact(j+k+1) * vscales(indj)
            do n = j, p
                indn = n*n + n + 1
                if(indmat .gt. (p+1)*(p+2)*(p+3)/6) then
                    write(*,*) "wrong max index in fmm_l2l_get_ztrans_mat"
                    exit
                end if
                mat(indmat) = tmp1 * pow_r1(n+1) / vscales(indn) &
                    & * fact(n-k+1) * fact(n+k+1) / fact(n-j+1) / fact(n-j+1)
                if (mod(n+j, 2) .eq. 1) then
                    mat(indmat) = -mat(indmat)
                end if
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_l2l_get_ztrans_mat

! Use precomputed translation matrices over OZ axis
! Check if centers of spherical harmonics are different, otherwise use function
! fmm_l2l_scale, since case of the same center requires only scaling
subroutine fmm_l2l_use_ztrans_mat(p, mat, src_l, dst_l)
! Parameters:
!   p: maximum degree of spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    real(kind=8), intent(out) :: dst_l((p+1)*(p+1))
    integer :: j, k, n, indj, indn, indmat
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        dst_l(indj) = 0
        ! k = 0
        do n = j, p
            indn = n*n + n + 1
            dst_l(indj) = dst_l(indj) + mat(indmat)*src_l(indn)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            dst_l(indj+k) = 0
            dst_l(indj-k) = 0
            do n = j, p
                indn = n*n + n + 1
                dst_l(indj+k) = dst_l(indj+k) + mat(indmat)*src_l(indn+k)
                dst_l(indj-k) = dst_l(indj-k) + mat(indmat)*src_l(indn-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_l2l_use_ztrans_mat

! Use precomputed matrices to perform adjoint translation over OZ axis
! Check if centers of spherical harmonics are different, otherwise use function
! fmm_l2l_scale, since case of the same center requires only scaling
subroutine fmm_l2l_adj_use_ztrans_mat(p, mat, src_l, dst_l)
! Parameters:
!   p: maximum degree of spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    real(kind=8), intent(out) :: dst_l((p+1)*(p+1))
    integer :: j, k, n, indj, indn, indmat
    indmat = 1
    dst_l = 0
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = j, p
            indn = n*n + n + 1
            dst_l(indn) = dst_l(indn) + mat(indmat)*src_l(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = j, p
                indn = n*n + n + 1
                dst_l(indn+k) = dst_l(indn+k) + mat(indmat)*src_l(indj+k)
                dst_l(indn-k) = dst_l(indn-k) + mat(indmat)*src_l(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_l2l_adj_use_ztrans_mat

! Use precomputed translation matrices over OZ axis
! Check if centers of spherical harmonics are different, otherwise use function
! fmm_l2l_scale, since case of the same center requires only scaling
! Only difference with fmm_l2l_use_ztrans_mat is that output dst_l is treated
! as inout, not just out. This is done to improve performance.
subroutine fmm_l2l_use_ztrans_mat2(p, mat, src_l, dst_l)
! Parameters:
!   p: maximum degree of spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    integer :: j, k, n, indj, indn, indmat
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = j, p
            indn = n*n + n + 1
            dst_l(indj) = dst_l(indj) + mat(indmat)*src_l(indn)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = j, p
                indn = n*n + n + 1
                dst_l(indj+k) = dst_l(indj+k) + mat(indmat)*src_l(indn+k)
                dst_l(indj-k) = dst_l(indj-k) + mat(indmat)*src_l(indn-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_l2l_use_ztrans_mat2

! Use precomputed matrices to perform adjoint translation over OZ axis
! Check if centers of spherical harmonics are different, otherwise use function
! fmm_l2l_scale, since case of the same center requires only scaling
! Only difference with fmm_l2l_use_ztrans_mat is that output dst_l is treated
! as inout, not just out. This is done to improve performance.
subroutine fmm_l2l_adj_use_ztrans_mat2(p, mat, src_l, dst_l)
! Parameters:
!   p: maximum degree of spherical harmonics
!   mat: translation matrix for spherical harmonics
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    integer :: j, k, n, indj, indn, indmat
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = j, p
            indn = n*n + n + 1
            dst_l(indn) = dst_l(indn) + mat(indmat)*src_l(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = j, p
                indn = n*n + n + 1
                dst_l(indn+k) = dst_l(indn+k) + mat(indmat)*src_l(indj+k)
                dst_l(indn-k) = dst_l(indn-k) + mat(indmat)*src_l(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_l2l_adj_use_ztrans_mat2

! Auxiliary routine to transform (rotate, reflect) spherical harmonics
! Factor u_lnm (lower case) from the paper:
!   ``Rotation Matrices for Real Spherical Harmonics. Direct Determination by
!   Recursion'', J.Ivanic and K.Ruedenberg, J. Phys. Chem. 1996, 100, 6342-6347
! This routine serves only for debugging purposes, as there are much faster
! implementations
subroutine fmm_sph_transform_get_u(l, n, m, c)
! Parameters:
!   l: order of spherical harmonics
!   n: index of initial spherical harmonic
!   m: index of output (rotated or reflected) spherical harmonic
!   c: value of factor u_lnm (lower case)
    integer, intent(in) :: l, n, m
    real(kind=8), intent(out) :: c
    c = l*l - n*n
    if (abs(m) .eq. l) then
        c = c / (2*l) / (2*l-1)
    else
        c = c / (l*l - m*m)
    end if
    c = sqrt(c)
end subroutine fmm_sph_transform_get_u

! Auxiliary routine to transform (rotate, reflect) spherical harmonics
! Factor v_lnm (lower case) from the paper:
!   ``Rotation Matrices for Real Spherical Harmonics. Direct Determination by
!   Recursion'', J.Ivanic and K.Ruedenberg, J. Phys. Chem. 1996, 100, 6342-6347
! This routine serves only for debugging purposes, as there are much faster
! implementations
subroutine fmm_sph_transform_get_v(l, n, m, c)
! Parameters:
!   l: order of spherical harmonics
!   n: index of initial spherical harmonic
!   m: index of output (rotated or reflected) spherical harmonic
!   c: value of factor v_lnm (lower case)
    integer, intent(in) :: l, n, m
    real(kind=8), intent(out) :: c
    integer :: k
    c = 1
    if (abs(m) .eq. l) then
        c = c / (2*l) / (2*l-1)
    else
        c = c / (l*l - m*m)
    end if
    if (n .eq. 0) then
        c = c * 2 * l * (l-1)
        c = -sqrt(c) / 2
    else
        k = abs(n)
        c = c * (l+k-1) * (l+k)
        c = sqrt(c) / 2
    end if
end subroutine fmm_sph_transform_get_v

! Auxiliary routine to transform (rotate, reflect) spherical harmonics
! Factor w_lnm (lower case) from the paper:
!   ``Rotation Matrices for Real Spherical Harmonics. Direct Determination by
!   Recursion'', J.Ivanic and K.Ruedenberg, J. Phys. Chem. 1996, 100, 6342-6347
! This routine serves only for debugging purposes, as there are much faster
! implementations
subroutine fmm_sph_transform_get_w(l, n, m, c)
! Parameters:
!   l: order of spherical harmonics
!   n: index of initial spherical harmonic
!   m: index of output (rotated or reflected) spherical harmonic
!   c: value of factor w_lnm (lower case)
    integer, intent(in) :: l, n, m
    real(kind=8), intent(out) :: c
    integer :: k
    if (n .eq. 0) then
        c = 0
        return
    end if
    c = 1
    if (abs(m) .eq. l) then
        c = c / (2*l) / (2*l-1)
    else
        c = c / (l*l - m*m)
    end if
    k = abs(n)
    c = c * (l-k-1) * (l-k)
    c = -sqrt(c) / 2
end subroutine fmm_sph_transform_get_w

! Auxiliary routine to transform (rotate, reflect) spherical harmonics
! Factor P_linm from the paper:
!   ``Rotation Matrices for Real Spherical Harmonics. Direct Determination by
!   Recursion'', J.Ivanic and K.Ruedenberg, J. Phys. Chem. 1996, 100, 6342-6347
! This routine serves only for debugging purposes, as there are much faster
! implementations
subroutine fmm_sph_transform_get_p(l, i, n, m, r, r_prev, c)
! Parameters:
!   l: order of spherical harmonics
!   i: -1, 0 or 1
!   n: index of initial spherical harmonic
!   m: index of output (rotated or reflected) spherical harmonic
!   r: transformation matrix from initial axes to rotated axes
!   r_prev: transformation matrix for all spherical harmonics of order l-1
!   c: value of factor P_linm
    integer, intent(in) :: l, i, n, m
    real(kind=8), intent(in) :: r(-1:1, -1:1), r_prev(1-l:l-1, 1-l:l-1)
    real(kind=8), intent(out) :: c
    if (m .eq. l) then
        c = r(i, 1)*r_prev(n, l-1) - r(i, -1)*r_prev(n, 1-l)
    else if (m .eq. -l) then
        c = r(i, 1)*r_prev(n, 1-l) + r(i, -1)*r_prev(n, l-1)
    else
        c = r(i, 0) * r_prev(n, m)
    end if
end subroutine fmm_sph_transform_get_p

! Auxiliary routine to transform (rotate, reflect) spherical harmonics
! Factor u_lnm * U_lnm from the paper:
!   ``Rotation Matrices for Real Spherical Harmonics. Direct Determination by
!   Recursion'', J.Ivanic and K.Ruedenberg, J. Phys. Chem. 1996, 100, 6342-6347
! This routine serves only for debugging purposes, as there are much faster
! implementations
subroutine fmm_sph_transform_get_uu(l, n, m, r, r_prev, c)
! Parameters:
!   l: order of spherical harmonics
!   n: index of initial spherical harmonic
!   m: index of output (rotated or reflected) spherical harmonic
!   r: transformation matrix from initial axes to rotated axes
!   r_prev: transformation matrix for all spherical harmonics of order l-1
!   c: value of factor u_lnm * U_lnm
    integer, intent(in) :: l, n, m
    real(kind=8), intent(in) :: r(-1:1, -1:1), r_prev(1-l:l-1, 1-l:l-1)
    real(kind=8), intent(out) :: c
    real(kind=8) :: c1
    call fmm_sph_transform_get_u(l, n, m, c)
    if (c .ne. 0) then
        call fmm_sph_transform_get_p(l, 0, n, m, r, r_prev, c1)
        c = c * c1
    end if
end subroutine fmm_sph_transform_get_uu

! Auxiliary routine to transform (rotate, reflect) spherical harmonics
! Factor v_lnm * V_lnm from the paper:
!   ``Rotation Matrices for Real Spherical Harmonics. Direct Determination by
!   Recursion'', J.Ivanic and K.Ruedenberg, J. Phys. Chem. 1996, 100, 6342-6347
! This routine serves only for debugging purposes, as there are much faster
! implementations
subroutine fmm_sph_transform_get_vv(l, n, m, r, r_prev, c)
! Parameters:
!   l: order of spherical harmonics
!   n: index of initial spherical harmonic
!   m: index of output (rotated or reflected) spherical harmonic
!   r: transformation matrix from initial axes to rotated axes
!   r_prev: transformation matrix for all spherical harmonics of order l-1
!   c: value of factor v_lnm * V_lnm
    integer, intent(in) :: l, n, m
    real(kind=8), intent(in) :: r(-1:1, -1:1), r_prev(1-l:l-1, 1-l:l-1)
    real(kind=8), intent(out) :: c
    real(kind=8) :: c1, c2
    call fmm_sph_transform_get_v(l, n, m, c)
    if (c .eq. 0) then
        return
    end if
    if (n .gt. 1) then
        call fmm_sph_transform_get_p(l, 1, n-1, m, r, r_prev, c1)
        call fmm_sph_transform_get_p(l, -1, 1-n, m, r, r_prev, c2)
        c = c * (c1-c2)
    else if (n .eq. 1) then
        call fmm_sph_transform_get_p(l, 1, 0, m, r, r_prev, c1)
        c = sqrt_2 * c * c1
    else if (n .eq. 0) then
        call fmm_sph_transform_get_p(l, 1, 1, m, r, r_prev, c1)
        call fmm_sph_transform_get_p(l, -1, -1, m, r, r_prev, c2)
        c = c * (c1+c2)
    else if (n .eq. -1) then
        call fmm_sph_transform_get_p(l, -1, 0, m, r, r_prev, c1)
        c = sqrt_2 * c * c1
    else
        call fmm_sph_transform_get_p(l, 1, n+1, m, r, r_prev, c1)
        call fmm_sph_transform_get_p(l, -1, -1-n, m, r, r_prev, c2)
        c = c * (c1+c2)
    end if
end subroutine fmm_sph_transform_get_vv

! Auxiliary routine to transform (rotate, reflect) spherical harmonics
! Factor w_lnm * W_lnm from the paper:
!   ``Rotation Matrices for Real Spherical Harmonics. Direct Determination by
!   Recursion'', J.Ivanic and K.Ruedenberg, J. Phys. Chem. 1996, 100, 6342-6347
! This routine serves only for debugging purposes, as there are much faster
! implementations
subroutine fmm_sph_transform_get_ww(l, n, m, r, r_prev, c)
! Parameters:
!   l: order of spherical harmonics
!   n: index of initial spherical harmonic
!   m: index of output (rotated or reflected) spherical harmonic
!   r: transformation matrix from initial axes to rotated axes
!   r_prev: transformation matrix for all spherical harmonics of order l-1
!   c: value of factor w_lnm * W_lnm
    integer, intent(in) :: l, n, m
    real(kind=8), intent(in) :: r(-1:1, -1:1), r_prev(1-l:l-1, 1-l:l-1)
    real(kind=8), intent(out) :: c
    real(kind=8) :: c1, c2
    call fmm_sph_transform_get_w(l, n, m, c)
    if (c .eq. 0) then
        return
    end if
    if (n .lt. 0) then
        call fmm_sph_transform_get_p(l, 1, n-1, m, r, r_prev, c1)
        call fmm_sph_transform_get_p(l, -1, 1-n, m, r, r_prev, c2)
        c = c * (c1-c2)
    else
        call fmm_sph_transform_get_p(l, 1, n+1, m, r, r_prev, c1)
        call fmm_sph_transform_get_p(l, -1, -1-n, m, r, r_prev, c2)
        c = c * (c1+c2)
    end if
end subroutine fmm_sph_transform_get_ww

! Transform (rotate or reflect) spherical harmonics
! Baseline implementation (very slow)
subroutine fmm_sph_transform(p, r1, src, dst)
! Parameters:
!   p: maximum order of spherical harmonics
!   r1: transformation matrix from initial axes to output axes
!   src: coefficients of initial spherical harmonics
!   dst: coefficients of output (rotated) spherical harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: r1(-1:1, -1:1), src((p+1)*(p+1))
    real(kind=8), intent(out) :: dst((p+1)*(p+1))
    real(kind=8) :: u, v, w, r_prev(-p:p, -p:p), r(-p:p, -p:p)
    integer :: l, m, n, ind
    ! l = 0
    dst(1) = src(1)
    ! l = 1
    do m = -1, 1
        dst(3+m) = 0
        do n = -1, 1
            r_prev(n, m) = r1(n, m)
            dst(3+m) = dst(3+m) + r1(n, m)*src(3+n)
        end do
    end do
    ! l > 2
    do l = 2, p
        ind = l*l + l + 1
        do m = -l, l
            dst(ind+m) = 0
            do n = -l, l
                call fmm_sph_transform_get_uu(l, n, m, r1, &
                    & r_prev(1-l:l-1, 1-l:l-1), u)
                call fmm_sph_transform_get_vv(l, n, m, r1, &
                    & r_prev(1-l:l-1, 1-l:l-1), v)
                call fmm_sph_transform_get_ww(l, n, m, r1, &
                    & r_prev(1-l:l-1, 1-l:l-1), w)
                r(n, m) = u + v + w
                dst(ind+m) = dst(ind+m) + r(n, m)*src(ind+n)
            end do
        end do
        r_prev(-l:l, -l:l) = r(-l:l, -l:l)
    end do
end subroutine fmm_sph_transform

! Transform (rotate or reflect) spherical harmonics
! More or less optimized version, 1000 times faster than fmm_sph_transform, but
! it is slower, than fmm_sph_transform3
subroutine fmm_sph_transform2(p, r1, src, dst)
! Parameters:
!   p: maximum order of spherical harmonics
!   r1: transformation matrix from initial axes to output axes
!   src: coefficients of initial spherical harmonics
!   dst: coefficients of output (rotated) spherical harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: r1(-1:1, -1:1), src((p+1)*(p+1))
    real(kind=8), intent(out) :: dst((p+1)*(p+1))
    real(kind=8) :: u, v, w, r_prev(-p:p, -p:p), r(-p:p, -p:p)
    real(kind=8) :: tmpm1, tmpn1, tmpn2, tmpn3, tmpu1, uu, vv, ww
    integer :: l, m, n, ind
    ! l = 0
    dst(1) = src(1)
    ! l = 1
    do m = -1, 1
        dst(3+m) = 0
        do n = -1, 1
            r_prev(n, m) = r1(n, m)
            dst(3+m) = dst(3+m) + r1(n, m)*src(3+n)
        end do
    end do
    ! l > 1
    do l = 2, p
        ind = l*l + l + 1
        do m = 0, l
            dst(ind+m) = 0
            dst(ind-m) = 0
            if (m .ne. l) then
                tmpm1 = l*l - m*m
            else
                tmpm1 = 2*l
                tmpm1 = tmpm1 * (tmpm1-1)
            end if
            tmpm1 = sqrt(tmpm1)
            do n = 0, l
                tmpn1 = sqrt(dble(l*l-n*n))
                tmpn2 = l + n
                tmpn2 = sqrt(tmpn2*(tmpn2-1))
                if (n .eq. 0) then
                    tmpn2 = -tmpn2 / sqrt_2
                    tmpn3 = 0
                else
                    tmpn2 = tmpn2 / 2
                    tmpn3 = l - n
                    tmpn3 = -sqrt(tmpn3*(tmpn3-1)) / 2
                end if
                u = tmpn1 / tmpm1
                v = tmpn2 / tmpm1
                w = tmpn3 / tmpm1
                ! define V
                if (n .eq. 0) then
                    if (m .eq. 0) then
                        vv = v*(r1(1, 0)*r_prev(1, 0)+r1(-1, 0)*r_prev(-1, 0))
                        r(0, 0) = vv
                    else if (m .ne. l) then
                        vv = v*(r1(1, 0)*r_prev(1, m)+ &
                            & r1(-1, 0)*r_prev(-1, m))
                        r(0, m) = vv
                        vv = v*(r1(1, 0)*r_prev(1, -m)+ &
                            & r1(-1, 0)*r_prev(-1, -m))
                        r(0, -m) = vv
                    else
                        vv = v*(r1(1, 1)*r_prev(1, l-1)- &
                            & r1(1, -1)*r_prev(1, 1-l)+ &
                            & r1(-1, 1)*r_prev(-1, l-1)- &
                            & r1(-1, -1)*r_prev(-1, 1-l))
                        r(0, l) = vv
                        vv = v*(r1(1, 1)*r_prev(1, 1-l)+ &
                            & r1(1, -1)*r_prev(1, l-1)+ &
                            & r1(-1, 1)*r_prev(-1, 1-l)+ &
                            & r1(-1, -1)*r_prev(-1, l-1))
                        r(0, -l) = vv
                    end if
                else if (n .eq. 1) then
                    if (m .eq. 0) then
                        vv = sqrt_2*v*r1(1, 0)*r_prev(0, 0)
                        r(1, 0) = vv
                        vv = sqrt_2*v*r1(-1, 0)*r_prev(0, 0)
                        r(-1, 0) = vv
                    else if (m .ne. l) then
                        vv = sqrt_2*v*r1(1, 0)*r_prev(0, m)
                        r(1, m) = vv
                        vv = sqrt_2*v*r1(1, 0)*r_prev(0, -m)
                        r(1, -m) = vv
                        vv = sqrt_2*v*r1(-1, 0)*r_prev(0, m)
                        r(-1, m) = vv
                        vv = sqrt_2*v*r1(-1, 0)*r_prev(0, -m)
                        r(-1, -m) = vv
                    else
                        vv = sqrt_2*v*(r1(1, 1)* &
                            & r_prev(0, l-1)-r1(1, -1)*r_prev(0, 1-l))
                        r(1, l) = vv
                        vv = sqrt_2*v*(r1(1, 1)* &
                            & r_prev(0, 1-l)+r1(1, -1)*r_prev(0, l-1))
                        r(1, -l) = vv
                        vv = sqrt_2*v*(r1(-1, 1)* &
                            & r_prev(0, l-1)-r1(-1, -1)*r_prev(0, 1-l))
                        r(-1, l) = vv
                        vv = sqrt_2*v*(r1(-1, 1)* &
                            & r_prev(0, 1-l)+r1(-1, -1)*r_prev(0, l-1))
                        r(-1, -l) = vv
                    end if
                else
                    if (m .eq. 0) then
                        vv = v*(r1(1, 0)*r_prev(n-1, 0)- &
                            & r1(-1, 0)*r_prev(1-n, 0))
                        r(n, 0) = vv
                        vv = v*(r1(1, 0)*r_prev(1-n, 0)+ &
                            & r1(-1, 0)*r_prev(n-1, 0))
                        r(-n, 0) = vv
                    else if (m .ne. l) then
                        vv = v*(r1(1, 0)*r_prev(n-1, m)- &
                            & r1(-1, 0)*r_prev(1-n, m))
                        r(n, m) = vv
                        vv = v*(r1(1, 0)*r_prev(n-1, -m)- &
                            & r1(-1, 0)*r_prev(1-n, -m))
                        r(n, -m) = vv
                        vv = v*(r1(1, 0)*r_prev(1-n, m)+ &
                            & r1(-1, 0)*r_prev(n-1, m))
                        r(-n, m) = vv
                        vv = v*(r1(1, 0)*r_prev(1-n, -m)+ &
                            & r1(-1, 0)*r_prev(n-1, -m))
                        r(-n, -m) = vv
                    else
                        vv = v*(r1(1, 1)*r_prev(n-1, l-1)- &
                            & r1(1, -1)*r_prev(n-1, 1-l) - &
                            & r1(-1, 1)*r_prev(1-n, l-1) + &
                            & r1(-1, -1)*r_prev(1-n, 1-l))
                        r(n, l) = vv
                        vv = v*(r1(1, 1)*r_prev(n-1, 1-l)+ &
                            & r1(1, -1)*r_prev(n-1, l-1) - &
                            & r1(-1, 1)*r_prev(1-n, 1-l) - &
                            & r1(-1, -1)*r_prev(1-n, l-1))
                        r(n, -l) = vv
                        vv = v*(r1(1, 1)*r_prev(1-n, l-1)- &
                            & r1(1, -1)*r_prev(1-n, 1-l) + &
                            & r1(-1, 1)*r_prev(n-1, l-1) - &
                            & r1(-1, -1)*r_prev(n-1, 1-l))
                        r(-n, l) = vv
                        vv = v*(r1(1, 1)*r_prev(1-n, 1-l)+ &
                            & r1(1, -1)*r_prev(1-n, l-1) + &
                            & r1(-1, 1)*r_prev(n-1, 1-l) + &
                            & r1(-1, -1)*r_prev(n-1, l-1))
                        r(-n, -l) = vv
                    end if
                end if
                ! define U only if u is not zero to avoid out-of-bounds
                ! access on array r_prev, which happens in case abs(n)=l
                if (u .ne. 0) then
                    tmpu1 = u * r1(0, 0)
                    if (m .ne. l) then
                        if ((n .ne. 0) .and. (m .ne. 0)) then
                            uu = tmpu1 * r_prev(n, m)
                            r(n, m) = r(n, m) + uu
                            uu = tmpu1 * r_prev(-n, m)
                            r(-n, m) = r(-n, m) + uu
                            uu = tmpu1 * r_prev(n, -m)
                            r(n, -m) = r(n, -m) + uu
                            uu = tmpu1 * r_prev(-n, -m)
                            r(-n, -m) = r(-n, -m) + uu
                        else if (n .ne. 0) then
                            uu = tmpu1 * r_prev(n, 0)
                            r(n, 0) = r(n, 0) + uu
                            uu = tmpu1 * r_prev(-n, 0)
                            r(-n, 0) = r(-n, 0) + uu
                        else if (m .ne. 0) then
                            uu = tmpu1 * r_prev(0, m)
                            r(0, m) = r(0, m) + uu
                            uu = tmpu1 * r_prev(0, -m)
                            r(0, -m) = r(0, -m) + uu
                        else
                            uu = tmpu1 * r_prev(0, 0)
                            r(0, 0) = r(0, 0) + uu
                        end if
                    else
                        if (n .ne. 0) then
                            uu = u * (r1(0, 1)*r_prev(n, l-1)- &
                                & r1(0, -1)*r_prev(n, 1-l))
                            r(n, l) = r(n, l) + uu
                            uu = u * (r1(0, 1)*r_prev(-n, l-1)- &
                                & r1(0, -1)*r_prev(-n, 1-l))
                            r(-n, l) = r(-n, l) + uu
                            uu = u * (r1(0, 1)*r_prev(n, 1-l)+ &
                                & r1(0, -1)*r_prev(n, l-1))
                            r(n, -l) = r(n, -l) + uu
                            uu = u * (r1(0, 1)*r_prev(-n, 1-l)+ &
                                & r1(0, -1)*r_prev(-n, l-1))
                            r(-n, -l) = r(-n, -l) + uu
                        else
                            uu = u * (r1(0, 1)*r_prev(0, l-1)- &
                                & r1(0, -1)*r_prev(0, 1-l))
                            r(0, l) = r(0, l) + uu
                            uu = u * (r1(0, 1)*r_prev(0, 1-l)+ &
                                & r1(0, -1)*r_prev(0, l-1))
                            r(0, -l) = r(0, -l) + uu
                        end if
                    end if
                end if
                ! define W only if w is not zero (to avoid out-of-bounds
                ! access) on array r_prev, which happens if n=0 or abs(n)=l or
                ! abs(n)=l-1
                if (w .ne. 0) then
                    if (m .eq. 0) then
                        ww = w*(r1(1, 0)*r_prev(n+1, 0)+ &
                            & r1(-1, 0)*r_prev(-n-1, 0))
                        r(n, 0) = r(n, 0) + ww
                        ww = w*(r1(1, 0)*r_prev(-n-1, 0)- &
                            & r1(-1, 0)*r_prev(n+1, 0))
                        r(-n, 0) = r(-n, 0) + ww
                    else if (m .ne. l) then
                        ww = w*(r1(1, 0)*r_prev(n+1, m)+ &
                            & r1(-1, 0)*r_prev(-n-1, m))
                        r(n, m) = r(n, m) + ww
                        ww = w*(r1(1, 0)*r_prev(n+1, -m)+ &
                            & r1(-1, 0)*r_prev(-n-1, -m))
                        r(n, -m) = r(n, -m) + ww
                        ww = w*(r1(1, 0)*r_prev(-n-1, m)- &
                            & r1(-1, 0)*r_prev(n+1, m))
                        r(-n, m) = r(-n, m) + ww
                        ww = w*(r1(1, 0)*r_prev(-n-1, -m)- &
                            & r1(-1, 0)*r_prev(n+1, -m))
                        r(-n, -m) = r(-n, -m) + ww
                    else
                        ww = w*(r1(1, 1)*r_prev(n+1, l-1)- &
                            & r1(1, -1)*r_prev(n+1, 1-l)+ &
                            & r1(-1, 1)*r_prev(-n-1, l-1)- &
                            & r1(-1, -1)*r_prev(-n-1, 1-l))
                        r(n, l) = r(n, l) + ww
                        ww = w*(r1(1, 1)*r_prev(n+1, 1-l)+ &
                            & r1(1, -1)*r_prev(n+1, l-1)+ &
                            & r1(-1, 1)*r_prev(-n-1, 1-l)+ &
                            & r1(-1, -1)*r_prev(-n-1, l-1))
                        r(n, -l) = r(n, -l) + ww
                        ww = w*(r1(1, 1)*r_prev(-n-1, l-1)- &
                            & r1(1, -1)*r_prev(-n-1, 1-l)- &
                            & r1(-1, 1)*r_prev(n+1, l-1)+ &
                            & r1(-1, -1)*r_prev(n+1, 1-l))
                        r(-n, l) = r(-n, l) + ww
                        ww = w*(r1(1, 1)*r_prev(-n-1, 1-l)+ &
                            & r1(1, -1)*r_prev(-n-1, l-1)- &
                            & r1(-1, 1)*r_prev(n+1, 1-l)- &
                            & r1(-1, -1)*r_prev(n+1, l-1))
                        r(-n, -l) = r(-n, -l) + ww
                    end if
                end if
                ! Apply computed rotations
                if ((m .eq. 0) .and. (n .eq. 0)) then
                    dst(ind) = dst(ind) + r(0, 0)*src(ind)
                else if (m .eq. 0) then
                    dst(ind) = dst(ind) + r(n, 0)*src(ind+n) + &
                        & r(-n, 0)*src(ind-n)
                else if (n .eq. 0) then
                    dst(ind+m) = dst(ind+m) + r(0, m)*src(ind)
                    dst(ind-m) = dst(ind-m) + r(0, -m)*src(ind)
                else
                    dst(ind+m) = dst(ind+m) + r(n, m)*src(ind+n) + &
                        & r(-n, m)*src(ind-n)
                    dst(ind-m) = dst(ind-m) + r(n, -m)*src(ind+n) + &
                        & r(-n, -m)*src(ind-n)
                end if
            end do
        end do
        if (l .ne. p) then
            r_prev(-l:l, -l:l) = r(-l:l, -l:l)
        end if
    end do
end subroutine fmm_sph_transform2

! Transform (rotate or reflect) spherical harmonics
! Most optimized version (faster than fmm_sph_transform2)
subroutine fmm_sph_transform3(p, r1, src, dst)
! Parameters:
!   p: maximum order of spherical harmonics
!   r1: transformation matrix from initial axes to output axes
!   src: coefficients of initial spherical harmonics
!   dst: coefficients of output (rotated or reflected) spherical harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: r1(-1:1, -1:1), src((p+1)*(p+1))
    real(kind=8), intent(out) :: dst((p+1)*(p+1))
    real(kind=8) :: u, v, w, r_prev(-p:p, -p:p), r(-p:p, -p:p)
    real(kind=8) :: scal_uvw_m(-p:p), scal_u_n(-p:p), scal_v_n(-p:p)
    real(kind=8) :: scal_w_n(-p:p)
    integer :: l, m, n, ind
    ! Compute rotations/reflections
    ! l = 0
    dst(1) = src(1)
    if (p .eq. 0) then
        return
    end if
    ! l = 1
    do m = -1, 1
        dst(3+m) = 0
        do n = -1, 1
            dst(3+m) = dst(3+m) + r1(n, m)*src(3+n)
        end do
    end do
    if (p .eq. 1) then
        return
    end if
    ! l = 2
    r(2, 2) = (r1(1, 1)*r1(1, 1) - r1(1, -1)*r1(1, -1) - &
        & r1(-1, 1)*r1(-1, 1) + r1(-1, -1)*r1(-1, -1)) / 2
    r(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
    r(2, 0) = sqrt(dble(3)) / 2 * (r1(1, 0)*r1(1, 0) - r1(-1, 0)*r1(-1, 0))
    r(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
    r(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
    r(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
    r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
    r(1, 0) = sqrt(dble(3)) * r1(1, 0) * r1(0, 0)
    r(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
    r(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
    r(0, 2) = sqrt(dble(3)) / 2 * (r1(0, 1)*r1(0, 1) - r1(0, -1)*r1(0, -1))
    r(0, 1) = sqrt(dble(3)) * r1(0, 1) * r1(0, 0)
    r(0, 0) = (3*r1(0, 0)*r1(0, 0)-1) / 2
    r(0, -1) = sqrt(dble(3)) * r1(0, -1) * r1(0, 0)
    r(0, -2) = sqrt(dble(3)) * r1(0, 1) * r1(0, -1)
    r(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
    r(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
    r(-1, 0) = sqrt(dble(3)) * r1(-1, 0) * r1(0, 0)
    r(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(0, -1)*r1(-1, 0)
    r(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
    r(-2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
    r(-2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
    r(-2, 0) = sqrt(dble(3)) * r1(1, 0) * r1(-1, 0)
    r(-2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
    r(-2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
    do m = -2, 2
        dst(7+m) = 0
        do n = -2, 2
            dst(7+m) = dst(7+m) + src(7+n)*r(n, m)
        end do
    end do
    ! l > 2
    do l = 3, p
        ! Prepare previous rotation matrix
        r_prev(1-l:l-1, 1-l:l-1) = r(1-l:l-1, 1-l:l-1)
        ! Prepare scalar factors
        scal_uvw_m(0) = l
        do m = 1, l-1
        u = sqrt(dble(l*l-m*m))
        scal_uvw_m(m) = u
        scal_uvw_m(-m) = u
        end do
        u = 2 * l
        u = sqrt(dble(u*(u-1)))
        scal_uvw_m(l) = u
        scal_uvw_m(-l) = u
        scal_u_n(0) = l
        scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt_2
        scal_w_n(0) = 0
        do n = 1, l-2
        u = sqrt(dble(l*l-n*n))
        scal_u_n(n) = u
        scal_u_n(-n) = u
        v = l + n
        v = sqrt(v*(v-1)) / 2
        scal_v_n(n) = v
        scal_v_n(-n) = v
        w = l - n
        w = -sqrt(w*(w-1)) / 2
        scal_w_n(n) = w
        scal_w_n(-n) = w
        end do
        u = sqrt(dble(2*l-1))
        scal_u_n(l-1) = u
        scal_u_n(1-l) = u
        scal_u_n(l) = 0
        scal_u_n(-l) = 0
        v = sqrt(dble((2*l-1)*(2*l-2))) / 2
        scal_v_n(l-1) = v
        scal_v_n(1-l) = v
        v = sqrt(dble(2*l*(2*l-1))) / 2
        scal_v_n(l) = v
        scal_v_n(-l) = v
        scal_w_n(l-1) = 0
        scal_w_n(l) = 0
        scal_w_n(-l) = 0
        scal_w_n(1-l) = 0
        ind = l*l + l + 1
        ! m = l, n = l
        v = r1(1, 1)*r_prev(l-1, l-1) - r1(1, -1)*r_prev(l-1, 1-l) - &
            & r1(-1, 1)*r_prev(1-l, l-1) + r1(-1, -1)*r_prev(1-l, 1-l)
        r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind+l) = src(ind+l) * r(l, l)
        ! m = l, n = -l
        v = r1(1, 1)*r_prev(1-l, l-1) - r1(1, -1)*r_prev(1-l, 1-l) + &
            & r1(-1, 1)*r_prev(l-1, l-1) - r1(-1, -1)*r_prev(l-1, 1-l)
        r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
        ! m = l, n = l-1
        u = r1(0, 1)*r_prev(l-1, l-1) - r1(0, -1)*r_prev(l-1, 1-l)
        v = r1(1, 1)*r_prev(l-2, l-1) - r1(1, -1)*r_prev(l-2, 1-l) - &
            & r1(-1, 1)*r_prev(2-l, l-1) + r1(-1, -1)*r_prev(2-l, 1-l)
        r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
        ! m = l, n = 1-l
        u = r1(0, 1)*r_prev(1-l, l-1) - r1(0, -1)*r_prev(1-l, 1-l)
        v = r1(1, 1)*r_prev(2-l, l-1) - r1(1, -1)*r_prev(2-l, 1-l) + &
            & r1(-1, 1)*r_prev(l-2, l-1) - r1(-1, -1)*r_prev(l-2, 1-l)
        r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
        ! m = l, n = 1
        u = r1(0, 1)*r_prev(1, l-1) - r1(0, -1)*r_prev(1, 1-l)
        v = r1(1, 1)*r_prev(0, l-1) - r1(1, -1)*r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(2, l-1) - r1(1, -1)*r_prev(2, 1-l) + &
            & r1(-1, 1)*r_prev(-2, l-1) - r1(-1, -1)*r_prev(-2, 1-l)
        r(1, l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, l) = r(1, l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
        ! m = l, n = -1
        u = r1(0, 1)*r_prev(-1, l-1) - r1(0, -1)*r_prev(-1, 1-l)
        v = r1(-1, 1)*r_prev(0, l-1) - r1(-1, -1)*r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(-2, l-1) - r1(1, -1)*r_prev(-2, 1-l) - &
            & r1(-1, 1)*r_prev(2, l-1) + r1(-1, -1)*r_prev(2, 1-l)
        r(-1, l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, l) = r(-1, l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
        ! m = l, n = 0
        u = r1(0, 1)*r_prev(0, l-1) - r1(0, -1)*r_prev(0, 1-l)
        v = r1(1, 1)*r_prev(1, l-1) - r1(1, -1)*r_prev(1, 1-l) + &
            & r1(-1, 1)*r_prev(-1, l-1) - r1(-1, -1)*r_prev(-1, 1-l)
        r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
        ! m = l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1)*r_prev(n, l-1) - r1(0, -1)*r_prev(n, 1-l)
            v = r1(1, 1)*r_prev(n-1, l-1) - r1(1, -1)*r_prev(n-1, 1-l) - &
                & r1(-1, 1)*r_prev(1-n, l-1) + r1(-1, -1)*r_prev(1-n, 1-l)
            w = r1(1, 1)*r_prev(n+1, l-1) - r1(1, -1)*r_prev(n+1, 1-l) + &
                & r1(-1, 1)*r_prev(-n-1, l-1) - r1(-1, -1)*r_prev(-n-1, 1-l)
            r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, l) = r(n, l) / scal_uvw_m(l)
            dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
            u = r1(0, 1)*r_prev(-n, l-1) - r1(0, -1)*r_prev(-n, 1-l)
            v = r1(1, 1)*r_prev(1-n, l-1) - r1(1, -1)*r_prev(1-n, 1-l) + &
                & r1(-1, 1)*r_prev(n-1, l-1) - r1(-1, -1)*r_prev(n-1, 1-l)
            w = r1(1, 1)*r_prev(-n-1, l-1) - r1(1, -1)*r_prev(-n-1, 1-l) - &
                & r1(-1, 1)*r_prev(n+1, l-1) + r1(-1, -1)*r_prev(n+1, 1-l)
            r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, l) = r(-n, l) / scal_uvw_m(l)
            dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
        end do
        ! m = -l, n = l
        v = r1(1, 1)*r_prev(l-1, 1-l) + r1(1, -1)*r_prev(l-1, l-1) - &
            & r1(-1, 1)*r_prev(1-l, 1-l) - r1(-1, -1)*r_prev(1-l, l-1)
        r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind-l) = src(ind+l)*r(l, -l)
        ! m = -l, n = -l
        v = r1(1, 1)*r_prev(1-l, 1-l) + r1(1, -1)*r_prev(1-l, l-1) + &
            & r1(-1, 1)*r_prev(l-1, 1-l) + r1(-1, -1)*r_prev(l-1, l-1)
        r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
        ! m = -l, n = l-1
        u = r1(0, 1)*r_prev(l-1, 1-l) + r1(0, -1)*r_prev(l-1, l-1)
        v = r1(1, 1)*r_prev(l-2, 1-l) + r1(1, -1)*r_prev(l-2, l-1) - &
            & r1(-1, 1)*r_prev(2-l, 1-l) - r1(-1, -1)*r_prev(2-l, l-1)
        r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
        ! m = -l, n = 1-l
        u = r1(0, 1)*r_prev(1-l, 1-l) + r1(0, -1)*r_prev(1-l, l-1)
        v = r1(1, 1)*r_prev(2-l, 1-l) + r1(1, -1)*r_prev(2-l, l-1) + &
            & r1(-1, 1)*r_prev(l-2, 1-l) + r1(-1, -1)*r_prev(l-2, l-1)
        r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
        ! m = -l, n = 1
        u = r1(0, 1)*r_prev(1, 1-l) + r1(0, -1)*r_prev(1, l-1)
        v = r1(1, 1)*r_prev(0, 1-l) + r1(1, -1)*r_prev(0, l-1)
        w = r1(1, 1)*r_prev(2, 1-l) + r1(1, -1)*r_prev(2, l-1) + &
            & r1(-1, 1)*r_prev(-2, 1-l) + r1(-1, -1)*r_prev(-2, l-1)
        r(1, -l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, -l) = r(1, -l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
        ! m = -l, n = -1
        u = r1(0, 1)*r_prev(-1, 1-l) + r1(0, -1)*r_prev(-1, l-1)
        v = r1(-1, 1)*r_prev(0, 1-l) + r1(-1, -1)*r_prev(0, l-1)
        w = r1(1, 1)*r_prev(-2, 1-l) + r1(1, -1)*r_prev(-2, l-1) - &
            & r1(-1, 1)*r_prev(2, 1-l) - r1(-1, -1)*r_prev(2, l-1)
        r(-1, -l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
        ! m = -l, n = 0
        u = r1(0, 1)*r_prev(0, 1-l) + r1(0, -1)*r_prev(0, l-1)
        v = r1(1, 1)*r_prev(1, 1-l) + r1(1, -1)*r_prev(1, l-1) + &
            & r1(-1, 1)*r_prev(-1, 1-l) + r1(-1, -1)*r_prev(-1, l-1)
        r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
        ! m = -l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1)*r_prev(n, 1-l) + r1(0, -1)*r_prev(n, l-1)
            v = r1(1, 1)*r_prev(n-1, 1-l) + r1(1, -1)*r_prev(n-1, l-1) - &
                & r1(-1, 1)*r_prev(1-n, 1-l) - r1(-1, -1)*r_prev(1-n, l-1)
            w = r1(1, 1)*r_prev(n+1, 1-l) + r1(1, -1)*r_prev(n+1, l-1) + &
                & r1(-1, 1)*r_prev(-n-1, 1-l) + r1(-1, -1)*r_prev(-n-1, l-1)
            r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, -l) = r(n, -l) / scal_uvw_m(l)
            dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
            u = r1(0, 1)*r_prev(-n, 1-l) + r1(0, -1)*r_prev(-n, l-1)
            v = r1(1, 1)*r_prev(1-n, 1-l) + r1(1, -1)*r_prev(1-n, l-1) + &
                & r1(-1, 1)*r_prev(n-1, 1-l) + r1(-1, -1)*r_prev(n-1, l-1)
            w = r1(1, 1)*r_prev(-n-1, 1-l) + r1(1, -1)*r_prev(-n-1, l-1) - &
                & r1(-1, 1)*r_prev(n+1, 1-l) - r1(-1, -1)*r_prev(n+1, l-1)
            r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
            dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
        end do
        ! Now deal with m=1-l..l-1
        do m = 1-l, l-1
            ! n = l
            v = r1(1, 0)*r_prev(l-1, m) - r1(-1, 0)*r_prev(1-l, m)
            r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            dst(ind+m) = src(ind+l) * r(l, m)
            ! n = -l
            v = r1(1, 0)*r_prev(1-l, m) + r1(-1, 0)*r_prev(l-1, m)
            r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
            ! n = l-1
            u = r1(0, 0) * r_prev(l-1, m)
            v = r1(1, 0)*r_prev(l-2, m) - r1(-1, 0)*r_prev(2-l, m)
            r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
            ! n = 1-l
            u = r1(0, 0) * r_prev(1-l, m)
            v = r1(1, 0)*r_prev(2-l, m) + r1(-1, 0)*r_prev(l-2, m)
            r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
            ! n = 0
            u = r1(0, 0) * r_prev(0, m)
            v = r1(1, 0)*r_prev(1, m) + r1(-1, 0)*r_prev(-1, m)
            r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
            r(0, m) = r(0, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
            ! n = 1
            u = r1(0, 0) * r_prev(1, m)
            v = r1(1, 0) * r_prev(0, m)
            w = r1(1, 0)*r_prev(2, m) + r1(-1, 0)*r_prev(-2, m)
            r(1, m) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, m) = r(1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
            ! n = -1
            u = r1(0, 0) * r_prev(-1, m)
            v = r1(-1, 0) * r_prev(0, m)
            w = r1(1, 0)*r_prev(-2, m) - r1(-1, 0)*r_prev(2, m)
            r(-1, m) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, m) = r(-1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
            ! n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 0) * r_prev(n, m)
                v = r1(1, 0)*r_prev(n-1, m) - r1(-1, 0)*r_prev(1-n, m)
                w = r1(1, 0)*r_prev(n+1, m) + r1(-1, 0)*r_prev(-1-n, m)
                r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, m) = r(n, m) / scal_uvw_m(m)
                dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                u = r1(0, 0) * r_prev(-n, m)
                v = r1(1, 0)*r_prev(1-n, m) + r1(-1, 0)*r_prev(n-1, m)
                w = r1(1, 0)*r_prev(-n-1, m) - r1(-1, 0)*r_prev(n+1, m)
                r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, m) = r(-n, m) / scal_uvw_m(m)
                dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
            end do
        end do
    end do
end subroutine fmm_sph_transform3

! Save matrix of reflection of spherical harmonics
! Corresponds to fmm_sph_transform3. Since this implements reflection, inverse
! operation is the same, as direct operation (P*P=I).
subroutine fmm_sph_get_reflect_mat(p, r1, mat)
! Parameters:
!   p: maximum order of spherical harmonics
!   r1: transformation matrix from initial axes to output axes. Be aware, that
!       this must come in a special order: OY, OZ and OX. This is because
!       Y_1^{-1}, Y_1^{0} and Y_1^{1} are corresponding to OY, OZ and OX. So,
!       r1(-1,-1) is a transform from old OY to new OY, r1(-1,0) is a transform
!       from old OY to new OZ, r1(1,-1) is a transform from old OX to new OY.
!   mat: reflection matrices for all degrees of spherical harmonics up to `p`
    integer, intent(in) :: p
    real(kind=8), intent(in) :: r1(-1:1, -1:1)
    real(kind=8), intent(out) :: mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8) :: u, v, w, r_prev(-p:p, -p:p), r(-p:p, -p:p)
    real(kind=8) :: scal_uvw_m(-p:p), scal_u_n(-p:p), scal_v_n(-p:p)
    real(kind=8) :: scal_w_n(-p:p)
    integer :: l, m, n, ind
    ! Compute rotations/reflections
    ! l = 0
    !dst(1) = src(1)
    mat(1) = 1
    if (p .eq. 0) then
        return
    end if
    ! l = 1
    do m = -1, 1
        !dst(3+m) = 0
        ind = 3*m + 6
        do n = -1, 1
            !dst(3+m) = dst(3+m) + r1(n, m)*src(3+n)
            mat(ind+n) = r1(n, m)
        end do
    end do
    if (p .eq. 1) then
        return
    end if
    ! l = 2
    r(2, 2) = (r1(1, 1)*r1(1, 1) - r1(1, -1)*r1(1, -1) - &
        & r1(-1, 1)*r1(-1, 1) + r1(-1, -1)*r1(-1, -1)) / 2
    r(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
    r(2, 0) = sqrt(dble(3)) / 2 * (r1(1, 0)*r1(1, 0) - r1(-1, 0)*r1(-1, 0))
    r(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
    r(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
    r(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
    r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
    r(1, 0) = sqrt(dble(3)) * r1(1, 0) * r1(0, 0)
    r(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
    r(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
    r(0, 2) = sqrt(dble(3)) / 2 * (r1(0, 1)*r1(0, 1) - r1(0, -1)*r1(0, -1))
    r(0, 1) = sqrt(dble(3)) * r1(0, 1) * r1(0, 0)
    r(0, 0) = (3*r1(0, 0)*r1(0, 0)-1) / 2
    r(0, -1) = sqrt(dble(3)) * r1(0, -1) * r1(0, 0)
    r(0, -2) = sqrt(dble(3)) * r1(0, 1) * r1(0, -1)
    r(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
    r(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
    r(-1, 0) = sqrt(dble(3)) * r1(-1, 0) * r1(0, 0)
    r(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(0, -1)*r1(-1, 0)
    r(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
    r(-2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
    r(-2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
    r(-2, 0) = sqrt(dble(3)) * r1(1, 0) * r1(-1, 0)
    r(-2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
    r(-2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
    do m = -2, 2
        !dst(7+m) = 0
        ind = 5*m + 23
        do n = -2, 2
            !dst(7+m) = dst(7+m) + src(7+n)*r(n, m)
            mat(ind+n) = r(n, m)
        end do
    end do
    ! l > 2
    do l = 3, p
        ! Prepare previous rotation matrix
        r_prev(1-l:l-1, 1-l:l-1) = r(1-l:l-1, 1-l:l-1)
        ! Prepare scalar factors
        scal_uvw_m(0) = l
        do m = 1, l-1
        u = sqrt(dble(l*l-m*m))
        scal_uvw_m(m) = u
        scal_uvw_m(-m) = u
        end do
        u = 2 * l
        u = sqrt(dble(u*(u-1)))
        scal_uvw_m(l) = u
        scal_uvw_m(-l) = u
        scal_u_n(0) = l
        scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt_2
        scal_w_n(0) = 0
        do n = 1, l-2
        u = sqrt(dble(l*l-n*n))
        scal_u_n(n) = u
        scal_u_n(-n) = u
        v = l + n
        v = sqrt(v*(v-1)) / 2
        scal_v_n(n) = v
        scal_v_n(-n) = v
        w = l - n
        w = -sqrt(w*(w-1)) / 2
        scal_w_n(n) = w
        scal_w_n(-n) = w
        end do
        u = sqrt(dble(2*l-1))
        scal_u_n(l-1) = u
        scal_u_n(1-l) = u
        scal_u_n(l) = 0
        scal_u_n(-l) = 0
        v = sqrt(dble((2*l-1)*(2*l-2))) / 2
        scal_v_n(l-1) = v
        scal_v_n(1-l) = v
        v = sqrt(dble(2*l*(2*l-1))) / 2
        scal_v_n(l) = v
        scal_v_n(-l) = v
        scal_w_n(l-1) = 0
        scal_w_n(l) = 0
        scal_w_n(-l) = 0
        scal_w_n(1-l) = 0
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
        r(1, l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, l) = r(1, l) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
        ! m = l, n = -1
        u = r1(0, 1)*r_prev(-1, l-1) - r1(0, -1)*r_prev(-1, 1-l)
        v = r1(-1, 1)*r_prev(0, l-1) - r1(-1, -1)*r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(-2, l-1) - r1(1, -1)*r_prev(-2, 1-l) - &
            & r1(-1, 1)*r_prev(2, l-1) + r1(-1, -1)*r_prev(2, 1-l)
        r(-1, l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
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
        r(1, -l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, -l) = r(1, -l) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
        ! m = -l, n = -1
        u = r1(0, 1)*r_prev(-1, 1-l) + r1(0, -1)*r_prev(-1, l-1)
        v = r1(-1, 1)*r_prev(0, 1-l) + r1(-1, -1)*r_prev(0, l-1)
        w = r1(1, 1)*r_prev(-2, 1-l) + r1(1, -1)*r_prev(-2, l-1) - &
            & r1(-1, 1)*r_prev(2, 1-l) - r1(-1, -1)*r_prev(2, l-1)
        r(-1, -l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
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
            r(1, m) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, m) = r(1, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
            ! n = -1
            u = r1(0, 0) * r_prev(-1, m)
            v = r1(-1, 0) * r_prev(0, m)
            w = r1(1, 0)*r_prev(-2, m) - r1(-1, 0)*r_prev(2, m)
            r(-1, m) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
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
        end do
        do m = -l, l
            ind = 2*l*(l+1)*(2*l+1)/3 + l + 1
            do n = -l, l
                mat(ind+m*(2*l+1)+n) = r(n, m)
            end do
        end do
    end do
end subroutine fmm_sph_get_reflect_mat

! Apply reflection matrices, saved by fmm_sph_get_reflect_mat
! This version does matrix-vector product without BLAS
subroutine fmm_sph_use_reflect_mat_baseline(p, mat, src, dst)
! Parameters:
!   p: maximum order of spherical harmonics
!   mat: matrices of reflections for all degrees of spherical harmnics
!   src: coefficients of initial spherical harmonics
!   dst: coefficients of output (rotated) spherical harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8), intent(in) :: src((p+1)*(p+1))
    real(kind=8), intent(out) :: dst((p+1)*(p+1))
    integer :: l, m, n, ind, indl
    ! l = 0
    dst(1) = src(1)
    if(p .eq. 0) then
        return
    end if
    do l = 1, p
        ! magical value for the offset to the current reflection matrix
        ind = 2*l*(l+1)*(2*l+1)/3 + l + 1
        ! offset for current spherical harmonics
        indl = l*l + l + 1
        do m = -l, l
            dst(indl+m) = 0
            do n = -l, l
                dst(indl+m) = dst(indl+m) + mat(ind+m*(2*l+1)+n)*src(indl+n)
            end do
        end do
    end do
end subroutine fmm_sph_use_reflect_mat_baseline

! Apply reflection matrices, saved by fmm_sph_get_reflect_mat
! This version does matrix-vector product with BLAS
subroutine fmm_sph_use_reflect_mat(p, mat, src, dst)
! Parameters:
!   p: maximum order of spherical harmonics
!   mat: matrices of reflections for all degrees of spherical harmnics
!   src: coefficients of initial spherical harmonics
!   dst: coefficients of output (rotated) spherical harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8), intent(in) :: src((p+1)*(p+1))
    real(kind=8), intent(out) :: dst((p+1)*(p+1))
    integer :: l, m, n, ind, indl
    real(kind=8) :: one=1.0, zero=0.0
    external :: dgemv
    ! l = 0
    dst(1) = src(1)
    if(p .eq. 0) then
        return
    end if
    do l = 1, p
        ! for small l do matvec manually
        if (l .lt. 6) then
            ! magical value for the offset to the current reflection matrix
            ind = 2*l*(l+1)*(2*l+1)/3 + l + 1
            ! offset for current spherical harmonics
            indl = l*l + l + 1
            do m = -l, l
                dst(indl+m) = 0
                do n = -l, l
                    dst(indl+m) = dst(indl+m) + &
                        & mat(ind+m*(2*l+1)+n)*src(indl+n)
                end do
            end do
        ! for larger l do matvec by blas
        else
            ! magical value for the offset to the current reflection matrix
            ind = l*(2*l-1)*(2*l+1)/3 + 1
            ! offset for current spherical harmonics
            indl = l*l + 1
            m = 2*l+1
            call dgemv('T', m, m, one, mat(ind), m, src(indl), 1, zero, &
                & dst(indl), 1)
        end if
    end do
end subroutine fmm_sph_use_reflect_mat

! Apply reflection matrices, saved by fmm_sph_get_reflect_mat
! Only difference with fmm_sph_reflect_mat is that output dst is treated as
! inout, not just out. This is done to improve performance
subroutine fmm_sph_use_reflect_mat2(p, mat, src, dst)
! Parameters:
!   p: maximum order of spherical harmonics
!   mat: matrices of reflections for all degrees of spherical harmnics
!   src: coefficients of initial spherical harmonics
!   dst: coefficients of output (rotated) spherical harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8), intent(in) :: src((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst((p+1)*(p+1))
    integer :: l, m, n, ind, indl
    real(kind=8) :: one=1.0, zero=0.0
    external :: dgemv
    ! l = 0
    dst(1) = dst(1) + src(1)
    if(p .eq. 0) then
        return
    end if
    do l = 1, p
        ! for small l do matvec manually
        if (l .lt. 6) then
            ! magical value for the offset to the current reflection matrix
            ind = 2*l*(l+1)*(2*l+1)/3 + l + 1
            ! offset for current spherical harmonics
            indl = l*l + l + 1
            do m = -l, l
                do n = -l, l
                    dst(indl+m) = dst(indl+m) + &
                        & mat(ind+m*(2*l+1)+n)*src(indl+n)
                end do
            end do
        ! for larger l do matvec by blas
        else
            ! magical value for the offset to the current reflection matrix
            ind = l*(2*l-1)*(2*l+1)/3 + 1
            ! offset for current spherical harmonics
            indl = l*l + 1
            m = 2*l+1
            call dgemv('T', m, m, one, mat(ind), m, src(indl), 1, one, &
                & dst(indl), 1)
        end if
    end do
end subroutine fmm_sph_use_reflect_mat2

! Transform spherical harmonics in OXZ plane (y coordinate remains the same)
! Based on fmm_sph_transform3 by assuming r1(-1,i)=0 and r1(i,-1)=0 for i=0,1
! and r1(-1,-1)=1
subroutine fmm_sph_transform3_oxz(p, r1, src, dst)
! Parameters:
!   p: maximum order of spherical harmonics
!   r1: transformation matrix from initial axes to output axes
!   src: coefficients of initial spherical harmonics
!   dst: coefficients of output (rotated) spherical harmonics
    integer, intent(in) :: p
    real(kind=8), intent(in) :: r1(0:1, 0:1), src((p+1)*(p+1))
    real(kind=8), intent(out) :: dst((p+1)*(p+1))
    real(kind=8) :: u, v, w, r_prev(-p:p, -p:p), r(-p:p, -p:p)
    real(kind=8) :: scal_uvw_m(-p:p), scal_u_n(-p:p), scal_v_n(-p:p)
    real(kind=8) :: scal_w_n(-p:p)
    integer :: l, m, n, ind
    ! Compute rotations/reflections
    ! l = 0
    dst(1) = src(1)
    if (p .eq. 0) then
        return
    end if
    ! l = 1
    dst(2) = src(2)
    dst(3) = src(3)*r1(0, 0) + src(4)*r1(1, 0)
    dst(4) = src(3)*r1(0, 1) + src(4)*r1(1, 1)
    if (p .eq. 1) then
        return
    end if
    ! l = 2
    r(2, 2) = (r1(1, 1)*r1(1, 1) + 1) / 2
    r(2, 1) = r1(1, 1)*r1(1, 0)
    r(2, 0) = sqrt(dble(3)) / 2 * r1(1, 0) * r1(1, 0)
    r(2, -1) = 0
    r(2, -2) = 0
    r(1, 2) = r1(1, 1)*r1(0, 1)
    r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
    r(1, 0) = sqrt(dble(3)) * r1(1, 0) * r1(0, 0)
    r(1, -1) = 0
    r(1, -2) = 0
    r(0, 2) = sqrt(dble(3)) / 2 * r1(0, 1) * r1(0, 1)
    r(0, 1) = sqrt(dble(3)) * r1(0, 1) * r1(0, 0)
    r(0, 0) = (3*r1(0, 0)*r1(0, 0)-1) / 2
    r(0, -1) = 0
    r(0, -2) = 0
    r(-1, 2) = 0
    r(-1, 1) = 0
    r(-1, 0) = 0
    r(-1, -1) = r1(0, 0)
    r(-1, -2) = r1(0, 1)
    r(-2, 2) = 0
    r(-2, 1) = 0
    r(-2, 0) = 0
    r(-2, -1) = r1(1, 0)
    r(-2, -2) = r1(1, 1)
    dst(9) = src(9)*r(2, 2) + src(8)*r(1, 2) + src(7)*r(0, 2)
    dst(8) = src(9)*r(2, 1) + src(8)*r(1, 1) + src(7)*r(0, 1)
    dst(7) = src(9)*r(2, 0) + src(8)*r(1, 0) + src(7)*r(0, 0)
    dst(6) = src(6)*r(-1, -1) + src(5)*r(-2, -1)
    dst(5) = src(6)*r(-1, -2) + src(5)*r(-2, -2)
    ! l > 2
    do l = 3, p
        ! Prepare previous rotation matrix
        r_prev(1-l:l-1, 1-l:l-1) = r(1-l:l-1, 1-l:l-1)
        ! Prepare scalar factors
        scal_uvw_m(0) = l
        do m = 1, l-1
        u = sqrt(dble(l*l-m*m))
        scal_uvw_m(m) = u
        scal_uvw_m(-m) = u
        end do
        u = 2 * l
        u = sqrt(dble(u*(u-1)))
        scal_uvw_m(l) = u
        scal_uvw_m(-l) = u
        scal_u_n(0) = l
        scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt_2
        scal_w_n(0) = 0
        do n = 1, l-2
        u = sqrt(dble(l*l-n*n))
        scal_u_n(n) = u
        scal_u_n(-n) = u
        v = l + n
        v = sqrt(v*(v-1)) / 2
        scal_v_n(n) = v
        scal_v_n(-n) = v
        w = l - n
        w = -sqrt(w*(w-1)) / 2
        scal_w_n(n) = w
        scal_w_n(-n) = w
        end do
        u = sqrt(dble(2*l-1))
        scal_u_n(l-1) = u
        scal_u_n(1-l) = u
        scal_u_n(l) = 0
        scal_u_n(-l) = 0
        v = sqrt(dble((2*l-1)*(2*l-2))) / 2
        scal_v_n(l-1) = v
        scal_v_n(1-l) = v
        v = sqrt(dble(2*l*(2*l-1))) / 2
        scal_v_n(l) = v
        scal_v_n(-l) = v
        scal_w_n(l-1) = 0
        scal_w_n(l) = 0
        scal_w_n(-l) = 0
        scal_w_n(1-l) = 0
        ind = l*l + l + 1
        ! m = l, n = l
        v = r1(1, 1)*r_prev(l-1, l-1) + r_prev(1-l, 1-l)
        r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind+l) = src(ind+l) * r(l, l)
        ! m = l, n = -l
        v = r1(1, 1)*r_prev(1-l, l-1) - r_prev(l-1, 1-l)
        r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
        ! m = l, n = l-1
        u = r1(0, 1) * r_prev(l-1, l-1)
        v = r1(1, 1)*r_prev(l-2, l-1) + r_prev(2-l, 1-l)
        r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
        ! m = l, n = 1-l
        u = r1(0, 1) * r_prev(1-l, l-1)
        v = r1(1, 1)*r_prev(2-l, l-1) - r_prev(l-2, 1-l)
        r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
        ! m = l, n = 1
        u = r1(0, 1) * r_prev(1, l-1)
        v = r1(1, 1) * r_prev(0, l-1)
        w = r1(1, 1)*r_prev(2, l-1) - r_prev(-2, 1-l)
        r(1, l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, l) = r(1, l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
        ! m = l, n = -1
        u = r1(0, 1) * r_prev(-1, l-1)
        v = r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(-2, l-1) + r_prev(2, 1-l)
        r(-1, l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, l) = r(-1, l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
        ! m = l, n = 0
        u = r1(0, 1) * r_prev(0, l-1)
        v = r1(1, 1)*r_prev(1, l-1) - r_prev(-1, 1-l)
        r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
        ! m = l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1) * r_prev(n, l-1)
            v = r1(1, 1)*r_prev(n-1, l-1) + r_prev(1-n, 1-l)
            w = r1(1, 1)*r_prev(n+1, l-1) - r_prev(-n-1, 1-l)
            r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, l) = r(n, l) / scal_uvw_m(l)
            dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
            u = r1(0, 1) * r_prev(-n, l-1)
            v = r1(1, 1)*r_prev(1-n, l-1) - r_prev(n-1, 1-l)
            w = r1(1, 1)*r_prev(-n-1, l-1) + r_prev(n+1, 1-l)
            r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, l) = r(-n, l) / scal_uvw_m(l)
            dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
        end do
        ! m = -l, n = l
        v = r1(1, 1)*r_prev(l-1, 1-l) - r_prev(1-l, l-1)
        r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind-l) = src(ind+l)*r(l, -l)
        ! m = -l, n = -l
        v = r1(1, 1)*r_prev(1-l, 1-l) + r_prev(l-1, l-1)
        r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
        ! m = -l, n = l-1
        u = r1(0, 1) * r_prev(l-1, 1-l)
        v = r1(1, 1)*r_prev(l-2, 1-l) - r_prev(2-l, l-1)
        r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
        ! m = -l, n = 1-l
        u = r1(0, 1) * r_prev(1-l, 1-l)
        v = r1(1, 1)*r_prev(2-l, 1-l) + r_prev(l-2, l-1)
        r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
        ! m = -l, n = 1
        u = r1(0, 1) * r_prev(1, 1-l)
        v = r1(1, 1) * r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(2, 1-l) + r_prev(-2, l-1)
        r(1, -l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, -l) = r(1, -l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
        ! m = -l, n = -1
        u = r1(0, 1) * r_prev(-1, 1-l)
        v = r_prev(0, l-1)
        w = r1(1, 1)*r_prev(-2, 1-l) - r_prev(2, l-1)
        r(-1, -l) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
        ! m = -l, n = 0
        u = r1(0, 1) * r_prev(0, 1-l)
        v = r1(1, 1)*r_prev(1, 1-l) + r_prev(-1, l-1)
        r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
        ! m = -l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1) * r_prev(n, 1-l)
            v = r1(1, 1)*r_prev(n-1, 1-l) - r_prev(1-n, l-1)
            w = r1(1, 1)*r_prev(n+1, 1-l) + r_prev(-n-1, l-1)
            r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, -l) = r(n, -l) / scal_uvw_m(l)
            dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
            u = r1(0, 1) * r_prev(-n, 1-l)
            v = r1(1, 1)*r_prev(1-n, 1-l) + r_prev(n-1, l-1)
            w = r1(1, 1)*r_prev(-n-1, 1-l) - r_prev(n+1, l-1)
            r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
            dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
        end do
        ! Now deal with m=1-l..l-1
        do m = 1-l, l-1
            ! n = l
            v = r1(1, 0) * r_prev(l-1, m)
            r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            dst(ind+m) = src(ind+l) * r(l, m)
            ! n = -l
            v = r1(1, 0) * r_prev(1-l, m)
            r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
            ! n = l-1
            u = r1(0, 0) * r_prev(l-1, m)
            v = r1(1, 0) * r_prev(l-2, m)
            r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
            ! n = 1-l
            u = r1(0, 0) * r_prev(1-l, m)
            v = r1(1, 0) * r_prev(2-l, m)
            r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
            ! n = 0
            u = r1(0, 0) * r_prev(0, m)
            v = r1(1, 0) * r_prev(1, m)
            r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
            r(0, m) = r(0, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
            ! n = 1
            u = r1(0, 0) * r_prev(1, m)
            v = r1(1, 0) * r_prev(0, m)
            w = r1(1, 0) * r_prev(2, m)
            r(1, m) = u*scal_u_n(1) + sqrt_2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, m) = r(1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
            ! n = -1
            u = r1(0, 0) * r_prev(-1, m)
            w = r1(1, 0) * r_prev(-2, m)
            r(-1, m) = u*scal_u_n(1) + w*scal_w_n(1)
            r(-1, m) = r(-1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
            ! n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 0) * r_prev(n, m)
                v = r1(1, 0) * r_prev(n-1, m)
                w = r1(1, 0) * r_prev(n+1, m)
                r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, m) = r(n, m) / scal_uvw_m(m)
                dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                u = r1(0, 0) * r_prev(-n, m)
                v = r1(1, 0) * r_prev(1-n, m)
                w = r1(1, 0) * r_prev(-n-1, m)
                r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, m) = r(-n, m) / scal_uvw_m(m)
                dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
            end do
        end do
    end do
end subroutine fmm_sph_transform3_oxz

! M2M translation by reflection (p^3 operations)
! Uses single reflection (done by a fmm_sph_transform2). fmm_m2m_fast is faster
! than this function, although result is the same
subroutine fmm_m2m_reflection(c, src_r, dst_r, p, vscales, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: src_m((p+1)*(p+1)), vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: stheta, c1(3), r1(3, 3), norm, nsgn, tmp_m((p+1)*(p+1)), r
    real(kind=8) :: tmp_m2((p+1)*(p+1))
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        call fmm_m2m_ztranslate(c(3), src_r, dst_r, p, vscales, src_m, dst_m)
        return
    end if
    r = sqrt(stheta + c(3)*c(3))
    c1(1) = c(2) / r
    c1(2) = c(3) / r
    c1(3) = c(1) / r
    if (c1(2) .ge. 0) then
        nsgn = -1
    else
        nsgn = 1
    end if
    c1(2) = c1(2) - nsgn
    norm = sqrt(c1(1)*c1(1) + c1(2)*c1(2) + c1(3)*c1(3))
    c1 = c1 / norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - 2*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_transform2(p, r1, src_m, tmp_m)
    tmp_m2 = 0
    call fmm_m2m_ztranslate(nsgn*r, src_r, dst_r, p, vscales, tmp_m, tmp_m2)
    call fmm_sph_transform2(p, r1, tmp_m2, tmp_m)
    dst_m = dst_m + tmp_m
end subroutine fmm_m2m_reflection

! M2M translation by reflection (p^3 operations)
! Uses single reflection (done by a fmm_sph_transform3). fmm_m2m_fast is faster
! than this function, although result is the same
subroutine fmm_m2m_reflection3(c, src_r, dst_r, p, vscales, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: src_m((p+1)*(p+1)), vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: stheta, c1(3), r1(3, 3), norm, nsgn, tmp_m((p+1)*(p+1)), r
    real(kind=8) :: tmp_m2((p+1)*(p+1))
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        call fmm_m2m_ztranslate(c(3), src_r, dst_r, p, vscales, src_m, dst_m)
        return
    end if
    r = sqrt(stheta + c(3)*c(3))
    c1(1) = c(2) / r
    c1(2) = c(3) / r
    c1(3) = c(1) / r
    if (c1(2) .ge. 0) then
        nsgn = -1
    else
        nsgn = 1
    end if
    c1(2) = c1(2) - nsgn
    norm = sqrt(c1(1)*c1(1) + c1(2)*c1(2) + c1(3)*c1(3))
    c1 = c1 / norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - 2*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_transform3(p, r1, src_m, tmp_m)
    tmp_m2 = 0
    call fmm_m2m_ztranslate(nsgn*r, src_r, dst_r, p, vscales, tmp_m, tmp_m2)
    call fmm_sph_transform3(p, r1, tmp_m2, tmp_m)
    dst_m = dst_m + tmp_m
end subroutine fmm_m2m_reflection3

! Rotate spherical harmonics around OZ axis
! Rotate on angle phi, presented by cos(m*phi) and sin(m*phi)
subroutine fmm_sph_rotate_oz(p, vcos, vsin, src, dst)
! Parameters:
!   p: maximum order of spherical harmonics
!   vcos: values of cos(m*phi)
!   vsin: values of sin(m*phi)
!   src: coefficients of initial spherical harmonics
!   dst: coefficients of rotated spherical harmonics
    real(kind=8), intent(in) :: vcos(p+1), vsin(p+1), src((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: dst((p+1)*(p+1))
    integer :: l, m, ind
    dst(1) = src(1)
    do l = 1, p
        ind = l*l + l + 1
        dst(ind) = src(ind)
        do m = 1, l
            dst(ind+m) = src(ind+m)*vcos(1+m) - src(ind-m)*vsin(1+m)
            dst(ind-m) = src(ind+m)*vsin(1+m) + src(ind-m)*vcos(1+m)
        end do
    end do
end subroutine fmm_sph_rotate_oz

! M2M translation by rotation around OZ and OY (p^3 operations)
! This is the fastest M2M operation I have implemented
subroutine fmm_m2m_fast(c, src_r, dst_r, p, vscales, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: src_m((p+1)*(p+1)), vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: r1(2, 2), tmp_m((p+1)*(p+1)), r, ctheta, stheta
    real(kind=8) :: tmp_m2((p+1)*(p+1)), cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vmsin(p+1)
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. 0) then
        call fmm_m2m_ztranslate(c(3), src_r, dst_r, p, vscales, src_m, dst_m)
        return
    end if
    r = sqrt(stheta + c(3)*c(3))
    ctheta = c(3) / r
    stheta = sqrt(stheta)
    cphi = c(1) / stheta
    sphi = c(2) / stheta
    stheta = stheta / r
    call trgev(cphi, sphi, p, vcos, vsin)
    vmsin = -vsin
    call fmm_sph_rotate_oz(p, vcos, vmsin, src_m, tmp_m)
    r1(1, 1) = ctheta
    r1(1, 2) = -stheta
    r1(2, 1) = stheta
    r1(2, 2) = ctheta
    call fmm_sph_transform3_oxz(p, r1, tmp_m, tmp_m2)
    tmp_m = 0
    call fmm_m2m_ztranslate(r, src_r, dst_r, p, vscales, tmp_m2, tmp_m)
    r1(1, 2) = stheta
    r1(2, 1) = -stheta
    call fmm_sph_transform3_oxz(p, r1, tmp_m, tmp_m2)
    call fmm_sph_rotate_oz(p, vcos, vsin, tmp_m2, tmp_m)
    dst_m = dst_m + tmp_m
end subroutine fmm_m2m_fast

! M2M translation by reflection (p^3 operations)
! Computes and uses reflection and OZ-translation matrices. This procedure
! helps debugging fmm_sph_get_reflect_mat, fmm_sph_use_reflect_mat,
! fmm_m2m_get_ztrans_mat and fmm_m2m_use_ztrans_mat
subroutine fmm_m2m_test_mat(c, src_r, dst_r, p, vscales, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: src_m((p+1)*(p+1)), vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: stheta, c1(3), r1(3, 3), norm, nsgn, tmp_m((p+1)*(p+1)), r
    real(kind=8) :: tmp_m2((p+1)*(p+1))
    real(kind=8) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same just do scaling
        if(c(3) .eq. 0) then
            call fmm_m2m_scale(src_r, dst_r, p, src_m, dst_m)
        ! Otherwise compute ztrans matrix and use it
        else
            call fmm_m2m_get_ztrans_mat(c(3), src_r, dst_r, p, vscales, &
                & ztrans_mat)
            call fmm_m2m_use_ztrans_mat2(p, ztrans_mat, src_m, dst_m)
        end if
        return
    end if
    r = sqrt(stheta + c(3)*c(3))
    c1(1) = c(2) / r
    c1(2) = c(3) / r
    c1(3) = c(1) / r
    if (c(3) .ge. 0) then
        nsgn = -1
    else
        nsgn = 1
    end if
    c1(2) = c1(2) - nsgn
    norm = sqrt(c1(1)*c1(1) + c1(2)*c1(2) + c1(3)*c1(3))
    c1 = c1 / norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - 2*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_get_reflect_mat(p, r1, reflect_mat)
    call fmm_sph_use_reflect_mat(p, reflect_mat, src_m, tmp_m)
    call fmm_m2m_get_ztrans_mat(nsgn*r, src_r, dst_r, p, vscales, &
        & ztrans_mat)
    call fmm_m2m_use_ztrans_mat(p, ztrans_mat, tmp_m, tmp_m2)
    call fmm_sph_use_reflect_mat2(p, reflect_mat, tmp_m2, dst_m)
end subroutine fmm_m2m_test_mat

! Compute and store reflection and OZ-translation matrices.
! This function does nothing if c is a zero-vector, if c(1)=c(2)=0, then this
! function computes only ztrans_mat. Zero cases must be treated with simple
! scaling by fmm_m2m_scale and c(1)=c(2)=0 must be treated by applying only
! translation over OZ axis with ztrans_mat matrix.
subroutine fmm_m2m_get_mat(c, src_r, dst_r, p, vscales, reflect_mat, &
        & ztrans_mat)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   reflect_mat: reflection matrix to align vector c with OZ axis
!   ztrans_mat: translation over OZ axis after reflection
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(in) :: vscales((p+1)*(p+1))
    real(kind=8), intent(out) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8), intent(out) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8) :: stheta, c1(3), r1(3, 3), norm, nsgn, r
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same do nothing
        if(c(3) .eq. 0) then
        ! Otherwise compute ztrans matrix
        else
            call fmm_m2m_get_ztrans_mat(c(3), src_r, dst_r, p, vscales, &
                & ztrans_mat)
        end if
        return
    end if
    ! Compute reflection of old coordinates into new ones to generate
    ! reflection of spherical harmonics
    r = sqrt(stheta + c(3)*c(3))
    c1(1) = c(2) / r
    c1(2) = c(3) / r
    c1(3) = c(1) / r
    if (c(3) .ge. 0) then
        nsgn = -1
    else
        nsgn = 1
    end if
    c1(2) = c1(2) - nsgn
    norm = sqrt(c1(1)*c1(1) + c1(2)*c1(2) + c1(3)*c1(3))
    c1 = c1 / norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - 2*c1(n)*c1(m)
        end do
    end do
    ! Compute reflection and translation
    call fmm_sph_get_reflect_mat(p, r1, reflect_mat)
    call fmm_m2m_get_ztrans_mat(nsgn*r, src_r, dst_r, p, vscales, ztrans_mat)
end subroutine fmm_m2m_get_mat

! M2M translation by precomputed reflection and OZ tranlation (p^3 operations)
! This function applies precomputed reflection and translation to make FMM
! matvecs much faster.
subroutine fmm_m2m_use_mat(c, src_r, dst_r, p, reflect_mat, ztrans_mat, &
        & src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   reflect_mat: reflection matrix to align vector c with OZ axis
!   ztrans_mat: translation over OZ axis after reflection
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(in) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8), intent(in) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: stheta, nsgn, tmp_m((p+1)*(p+1))
    real(kind=8) :: tmp_m2((p+1)*(p+1))
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same just scale
        if(c(3) .eq. 0) then
            call fmm_m2m_scale(src_r, dst_r, p, src_m, dst_m)
        ! Otherwise compute ztrans matrix and use it
        else
            call fmm_m2m_use_ztrans_mat2(p, ztrans_mat, src_m, dst_m)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_use_reflect_mat(p, reflect_mat, src_m, tmp_m)
    call fmm_m2m_use_ztrans_mat(p, ztrans_mat, tmp_m, tmp_m2)
    call fmm_sph_use_reflect_mat2(p, reflect_mat, tmp_m2, dst_m)
end subroutine fmm_m2m_use_mat

! Adjoint M2M translation by precomputed reflection and OZ tranlation
! This function applies precomputed reflection and translation to make FMM
! matvecs much faster (in O(p^3) operations).
subroutine fmm_m2m_adj_use_mat(c, src_r, dst_r, p, reflect_mat, ztrans_mat, &
        & src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   reflect_mat: reflection matrix to align vector c with OZ axis
!   ztrans_mat: translation over OZ axis after reflection
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(in) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8), intent(in) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: stheta, nsgn, tmp_m((p+1)*(p+1))
    real(kind=8) :: tmp_m2((p+1)*(p+1))
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same just scale
        if(c(3) .eq. 0) then
            call fmm_m2m_scale(src_r, dst_r, p, src_m, dst_m)
        ! Otherwise compute ztrans matrix and use it
        else
            call fmm_m2m_adj_use_ztrans_mat2(p, ztrans_mat, src_m, dst_m)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_use_reflect_mat(p, reflect_mat, src_m, tmp_m)
    call fmm_m2m_adj_use_ztrans_mat(p, ztrans_mat, tmp_m, tmp_m2)
    call fmm_sph_use_reflect_mat2(p, reflect_mat, tmp_m2, dst_m)
end subroutine fmm_m2m_adj_use_mat

! M2L translation by rotation around OZ and OY (p^3 operations)
! Based on fmm_m2m_fast
subroutine fmm_m2l_fast(c, src_r, dst_r, pm, pl, vscales, src_m, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old (multipole) harmonics
!   dst_l: expansion in new (local) harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: src_m((pm+1)*(pm+1))
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    integer, intent(in) :: pm, pl
    real(kind=8), intent(inout) :: dst_l((pl+1)*(pl+1))
    real(kind=8) :: r1(2, 2), tmp_ml((max(pm,pl)+1)*(max(pm,pl)+1))
    real(kind=8) :: r, ctheta, stheta
    real(kind=8) :: tmp_ml2((max(pm,pl)+1)*(max(pm,pl)+1)), cphi, sphi
    real(kind=8) :: vcos(max(pm,pl)+1), vsin(max(pm,pl)+1)
    real(kind=8) :: vmsin(max(pm,pl)+1)
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. 0) then
        call fmm_m2l_ztranslate(c(3), src_r, dst_r, pm, pl, vscales, src_m, &
            & dst_l)
        return
    end if
    r = sqrt(stheta + c(3)*c(3))
    ctheta = c(3) / r
    stheta = sqrt(stheta)
    cphi = c(1) / stheta
    sphi = c(2) / stheta
    stheta = stheta / r
    call trgev(cphi, sphi, max(pm,pl), vcos, vsin)
    vmsin = -vsin
    call fmm_sph_rotate_oz(pm, vcos, vmsin, src_m, tmp_ml)
    r1(1, 1) = ctheta
    r1(1, 2) = -stheta
    r1(2, 1) = stheta
    r1(2, 2) = ctheta
    call fmm_sph_transform3_oxz(pm, r1, tmp_ml, tmp_ml2)
    tmp_ml = 0
    call fmm_m2l_ztranslate(r, src_r, dst_r, pm, pl, vscales, tmp_ml2, tmp_ml)
    r1(1, 2) = stheta
    r1(2, 1) = -stheta
    call fmm_sph_transform3_oxz(pl, r1, tmp_ml, tmp_ml2)
    call fmm_sph_rotate_oz(pl, vcos, vsin, tmp_ml2, tmp_ml)
    dst_l = dst_l + tmp_ml(1:(pl+1)*(pl+1))
end subroutine fmm_m2l_fast

! M2L translation by reflection (p^3 operations)
! Computes and uses reflection and OZ-translation matrices. This procedure
! helps debugging fmm_sph_get_reflect_mat, fmm_sph_use_reflect_mat,
! fmm_m2l_get_ztrans_mat and fmm_m2l_use_ztrans_mat
subroutine fmm_m2l_test_mat(c, src_r, dst_r, pm, pl, vscales, src_m, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old (multipole) harmonics
!   dst_l: expansion in new (local) harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: src_m((pm+1)*(pm+1))
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    integer, intent(in) :: pm, pl
    real(kind=8), intent(inout) :: dst_l((pl+1)*(pl+1))
    real(kind=8) :: stheta, c1(3), r1(3, 3), norm, nsgn, r
    real(kind=8) :: tmp_ml((max(pm,pl)+1)*(max(pm,pl)+1))
    real(kind=8) :: tmp_ml2((max(pm,pl)+1)*(max(pm,pl)+1))
    real(kind=8) :: reflect_mat((max(pm,pl)+1) * (2*max(pm,pl)+1) &
        & * (2*max(pm,pl)+3) / 3)
    real(kind=8) :: ztrans_mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same do nothing (but it must not happen)
        if(c(3) .eq. 0) then
        ! Otherwise compute ztrans matrix and use it
        else
            call fmm_m2l_get_ztrans_mat(c(3), src_r, dst_r, pm, pl, vscales, &
                & ztrans_mat)
            call fmm_m2l_use_ztrans_mat2(pm, pl, ztrans_mat, src_m, dst_l)
        end if
        return
    end if
    r = sqrt(stheta + c(3)*c(3))
    c1(1) = c(2) / r
    c1(2) = c(3) / r
    c1(3) = c(1) / r
    if (c(3) .ge. 0) then
        nsgn = -1
    else
        nsgn = 1
    end if
    c1(2) = c1(2) - nsgn
    norm = sqrt(c1(1)*c1(1) + c1(2)*c1(2) + c1(3)*c1(3))
    c1 = c1 / norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - 2*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_get_reflect_mat(max(pm,pl), r1, reflect_mat)
    call fmm_sph_use_reflect_mat(pm, reflect_mat, src_m, tmp_ml)
    call fmm_m2l_get_ztrans_mat(nsgn*r, src_r, dst_r, pm, pl, vscales, &
        & ztrans_mat)
    call fmm_m2l_use_ztrans_mat(pm, pl, ztrans_mat, tmp_ml, tmp_ml2)
    call fmm_sph_use_reflect_mat2(pl, reflect_mat, tmp_ml2, dst_l)
end subroutine fmm_m2l_test_mat

! Compute and store reflection and OZ-translation matrices.
subroutine fmm_m2l_get_mat(c, src_r, dst_r, pm, pl, vscales, reflect_mat, &
        & ztrans_mat)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm
!   reflect_mat: reflection matrix to align vector c with OZ axis
!   ztrans_mat: translation over OZ axis after reflection
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(out) :: reflect_mat((max(pm,pl)+1) &
        & * (2*max(pm,pl)+1) * (2*max(pm,pl)+3) / 3)
    real(kind=8), intent(out) :: ztrans_mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    real(kind=8) :: stheta, c1(3), r1(3, 3), norm, nsgn, r
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same do nothing (but this must not happen)
        if(c(3) .eq. 0) then
        ! Otherwise compute ztrans matrix
        else
            call fmm_m2l_get_ztrans_mat(c(3), src_r, dst_r, pm, pl, vscales, &
                & ztrans_mat)
        end if
        return
    end if
    ! Compute reflection of old coordinates into new ones to generate
    ! reflection of spherical harmonics
    r = sqrt(stheta + c(3)*c(3))
    c1(1) = c(2) / r
    c1(2) = c(3) / r
    c1(3) = c(1) / r
    if (c(3) .ge. 0) then
        nsgn = -1
    else
        nsgn = 1
    end if
    c1(2) = c1(2) - nsgn
    norm = sqrt(c1(1)*c1(1) + c1(2)*c1(2) + c1(3)*c1(3))
    c1 = c1 / norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - 2*c1(n)*c1(m)
        end do
    end do
    ! Compute reflection and translation
    call fmm_sph_get_reflect_mat(max(pm,pl), r1, reflect_mat)
    call fmm_m2l_get_ztrans_mat(nsgn*r, src_r, dst_r, pm, pl, vscales, &
        & ztrans_mat)
end subroutine fmm_m2l_get_mat

! M2L translation by precomputed reflection and OZ tranlation (p^3 operations)
! This function applies precomputed reflection and translation to make FMM
! matvecs much faster.
subroutine fmm_m2l_use_mat(c, src_r, dst_r, pm, pl, reflect_mat, ztrans_mat, &
        & src_m, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   reflect_mat: reflection matrix to align vector c with OZ axis
!   ztrans_mat: translation over OZ axis after reflection
!   src_m: expansion in old harmonics
!   dst_l: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: reflect_mat((max(pm,pl)+1) &
        & * (2*max(pm,pl)+1) * (2*max(pm,pl)+3) / 3)
    real(kind=8), intent(in) :: ztrans_mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    real(kind=8), intent(in) :: src_m((pm+1)*(pm+1))
    real(kind=8), intent(inout) :: dst_l((pl+1)*(pl+1))
    real(kind=8) :: stheta, nsgn, tmp_ml((max(pm,pl)+1)*(max(pm,pl)+1))
    real(kind=8) :: tmp_ml2((max(pm,pl)+1)*(max(pm,pl)+1))
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same do nothing (but this must not happen)
        if(c(3) .eq. 0) then
        ! Otherwise use ztrans matrix
        else
            call fmm_m2l_use_ztrans_mat2(pm, pl, ztrans_mat, src_m, dst_l)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_use_reflect_mat(pm, reflect_mat, src_m, tmp_ml)
    call fmm_m2l_use_ztrans_mat(pm, pl, ztrans_mat, tmp_ml, tmp_ml2)
    call fmm_sph_use_reflect_mat2(pl, reflect_mat, tmp_ml2, dst_l)
end subroutine fmm_m2l_use_mat

! Adjoint M2L translation by precomputed reflection and OZ tranlation
! This function applies precomputed reflection and translation to make FMM
! matvecs much faster (in O(p^3) operations).
subroutine fmm_m2l_adj_use_mat(c, src_r, dst_r, pm, pl, reflect_mat, ztrans_mat, &
        & src_l, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   reflect_mat: reflection matrix to align vector c with OZ axis
!   ztrans_mat: translation over OZ axis after reflection
!   src_l: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: reflect_mat((max(pm,pl)+1) &
        & * (2*max(pm,pl)+1) * (2*max(pm,pl)+3) / 3)
    real(kind=8), intent(in) :: ztrans_mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    real(kind=8), intent(inout) :: dst_m((pm+1)*(pm+1))
    real(kind=8), intent(in) :: src_l((pl+1)*(pl+1))
    real(kind=8) :: stheta, nsgn, tmp_l((pl+1)*(pl+1)), tmp_m((pm+1)*(pm+1))
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same do nothing (but this must not happen)
        if(c(3) .eq. 0) then
        ! Otherwise use adjoint ztrans matrix
        else
            call fmm_m2l_adj_use_ztrans_mat2(pm, pl, ztrans_mat, src_l, dst_m)
        end if
        return
    end if
    ! Apply reflection->adjoint translation->reflection
    call fmm_sph_use_reflect_mat(pl, reflect_mat, src_l, tmp_l)
    call fmm_m2l_adj_use_ztrans_mat(pm, pl, ztrans_mat, tmp_l, tmp_m)
    call fmm_sph_use_reflect_mat2(pm, reflect_mat, tmp_m, dst_m)
end subroutine fmm_m2l_adj_use_mat

! L2L translation by rotation around OZ and OY (p^3 operations)
! Based on fmm_m2m_fast
subroutine fmm_l2l_fast(c, src_r, dst_r, p, vscales, src_l, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: src_l((p+1)*(p+1)), vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    real(kind=8) :: r1(2, 2), tmp_l((p+1)*(p+1)), r, ctheta, stheta
    real(kind=8) :: tmp_l2((p+1)*(p+1)), cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vmsin(p+1)
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. 0) then
        call fmm_l2l_ztranslate(c(3), src_r, dst_r, p, vscales, src_l, dst_l)
        return
    end if
    r = sqrt(stheta + c(3)*c(3))
    ctheta = c(3) / r
    stheta = sqrt(stheta)
    cphi = c(1) / stheta
    sphi = c(2) / stheta
    stheta = stheta / r
    call trgev(cphi, sphi, p, vcos, vsin)
    vmsin = -vsin
    call fmm_sph_rotate_oz(p, vcos, vmsin, src_l, tmp_l)
    r1(1, 1) = ctheta
    r1(1, 2) = -stheta
    r1(2, 1) = stheta
    r1(2, 2) = ctheta
    call fmm_sph_transform3_oxz(p, r1, tmp_l, tmp_l2)
    tmp_l = 0
    call fmm_l2l_ztranslate(r, src_r, dst_r, p, vscales, tmp_l2, tmp_l)
    r1(1, 2) = stheta
    r1(2, 1) = -stheta
    call fmm_sph_transform3_oxz(p, r1, tmp_l, tmp_l2)
    call fmm_sph_rotate_oz(p, vcos, vsin, tmp_l2, tmp_l)
    dst_l = dst_l + tmp_l
end subroutine fmm_l2l_fast

! L2L translation by reflection (p^3 operations)
! Computes and uses reflection and OZ-translation matrices. This procedure
! helps debugging fmm_sph_get_reflect_mat, fmm_sph_use_reflect_mat,
! fmm_l2l_get_ztrans_mat and fmm_l2l_use_ztrans_mat
subroutine fmm_l2l_test_mat(c, src_r, dst_r, p, vscales, src_l, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: src_l((p+1)*(p+1)), vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    real(kind=8) :: stheta, c1(3), r1(3, 3), norm, nsgn, tmp_l((p+1)*(p+1)), r
    real(kind=8) :: tmp_l2((p+1)*(p+1))
    real(kind=8) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same just do scaling
        if(c(3) .eq. 0) then
            call fmm_l2l_scale(src_r, dst_r, p, src_l, dst_l)
        ! Otherwise compute ztrans matrix and use it
        else
            call fmm_l2l_get_ztrans_mat(c(3), src_r, dst_r, p, vscales, &
                & ztrans_mat)
            call fmm_l2l_use_ztrans_mat2(p, ztrans_mat, src_l, dst_l)
        end if
        return
    end if
    r = sqrt(stheta + c(3)*c(3))
    c1(1) = c(2) / r
    c1(2) = c(3) / r
    c1(3) = c(1) / r
    if (c(3) .ge. 0) then
        nsgn = -1
    else
        nsgn = 1
    end if
    c1(2) = c1(2) - nsgn
    norm = sqrt(c1(1)*c1(1) + c1(2)*c1(2) + c1(3)*c1(3))
    c1 = c1 / norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - 2*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_get_reflect_mat(p, r1, reflect_mat)
    call fmm_sph_use_reflect_mat(p, reflect_mat, src_l, tmp_l)
    call fmm_l2l_get_ztrans_mat(nsgn*r, src_r, dst_r, p, vscales, &
        & ztrans_mat)
    call fmm_l2l_use_ztrans_mat(p, ztrans_mat, tmp_l, tmp_l2)
    call fmm_sph_use_reflect_mat2(p, reflect_mat, tmp_l2, dst_l)
end subroutine fmm_l2l_test_mat

! Compute and store reflection and OZ-translation matrices.
! This function does nothing if c is a zero-vector, if c(1)=c(2)=0, then this
! function computes only ztrans_mat. Zero cases must be treated with simple
! scaling by fmm_l2l_scale and c(1)=c(2)=0 must be treated by applying only
! translation over OZ axis with ztrans_mat matrix.
subroutine fmm_l2l_get_mat(c, src_r, dst_r, p, vscales, reflect_mat, &
        & ztrans_mat)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   reflect_mat: reflection matrix to align vector c with OZ axis
!   ztrans_mat: translation over OZ axis after reflection
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(in) :: vscales((p+1)*(p+1))
    real(kind=8), intent(out) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8), intent(out) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8) :: stheta, c1(3), r1(3, 3), norm, nsgn, tmp_m((p+1)*(p+1)), r
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same do nothing
        if(c(3) .eq. 0) then
        ! Otherwise compute ztrans matrix
        else
            call fmm_l2l_get_ztrans_mat(c(3), src_r, dst_r, p, vscales, &
                & ztrans_mat)
        end if
        return
    end if
    ! Compute reflection of old coordinates into new ones to generate
    ! reflection of spherical harmonics
    r = sqrt(stheta + c(3)*c(3))
    c1(1) = c(2) / r
    c1(2) = c(3) / r
    c1(3) = c(1) / r
    if (c(3) .ge. 0) then
        nsgn = -1
    else
        nsgn = 1
    end if
    c1(2) = c1(2) - nsgn
    norm = sqrt(c1(1)*c1(1) + c1(2)*c1(2) + c1(3)*c1(3))
    c1 = c1 / norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - 2*c1(n)*c1(m)
        end do
    end do
    ! Compute reflection and translation
    call fmm_sph_get_reflect_mat(p, r1, reflect_mat)
    call fmm_l2l_get_ztrans_mat(nsgn*r, src_r, dst_r, p, vscales, ztrans_mat)
end subroutine fmm_l2l_get_mat

! L2L translation by precomputed reflection and OZ tranlation (p^3 operations)
! This function applies precomputed reflection and translation to make FMM
! matvecs much faster.
subroutine fmm_l2l_use_mat(c, src_r, dst_r, p, reflect_mat, ztrans_mat, &
        & src_l, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   reflect_mat: reflection matrix to align vector c with OZ axis
!   ztrans_mat: translation over OZ axis after reflection
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(in) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8), intent(in) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    real(kind=8) :: stheta, nsgn, tmp_l((p+1)*(p+1))
    real(kind=8) :: tmp_l2((p+1)*(p+1))
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same just scale
        if(c(3) .eq. 0) then
            call fmm_l2l_scale(src_r, dst_r, p, src_l, dst_l)
        ! Otherwise use ztrans matrix
        else
            call fmm_l2l_use_ztrans_mat2(p, ztrans_mat, src_l, dst_l)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_use_reflect_mat(p, reflect_mat, src_l, tmp_l)
    call fmm_l2l_use_ztrans_mat(p, ztrans_mat, tmp_l, tmp_l2)
    call fmm_sph_use_reflect_mat2(p, reflect_mat, tmp_l2, dst_l)
end subroutine fmm_l2l_use_mat

! Adjoint L2L translation by precomputed reflection and OZ tranlation
! This function applies precomputed reflection and translation to make FMM
! matvecs much faster (in O(p^3) operations).
subroutine fmm_l2l_adj_use_mat(c, src_r, dst_r, p, reflect_mat, ztrans_mat, &
        & src_l, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   reflect_mat: reflection matrix to align vector c with OZ axis
!   ztrans_mat: translation over OZ axis after reflection
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(in) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8), intent(in) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    real(kind=8) :: stheta, nsgn, tmp_l((p+1)*(p+1))
    real(kind=8) :: tmp_l2((p+1)*(p+1))
    integer :: n, m
    stheta = c(1)*c(1) + c(2)*c(2)
    ! If no need for transformation, just do translation along z
    if (stheta .eq. 0) then
        ! If centers are the same just scale
        if(c(3) .eq. 0) then
            call fmm_l2l_scale(src_r, dst_r, p, src_l, dst_l)
        ! Otherwise use adjoint ztrans matrix
        else
            call fmm_l2l_adj_use_ztrans_mat2(p, ztrans_mat, src_l, dst_l)
        end if
        return
    end if
    ! Apply reflection->adjoint translation->reflection
    call fmm_sph_use_reflect_mat(p, reflect_mat, src_l, tmp_l)
    call fmm_l2l_adj_use_ztrans_mat(p, ztrans_mat, tmp_l, tmp_l2)
    call fmm_sph_use_reflect_mat2(p, reflect_mat, tmp_l2, dst_l)
end subroutine fmm_l2l_adj_use_mat

! Integrate spherical harmonics (grid -> coefficients)
subroutine int_grid(p, ngrid, w, vgrid, x, xlm)
! Parameters:
!   p: maximum degree of spherical harmonics
!   ngrid: number of Lebedev grid points
!   w: weights of Lebedev grid
!   vgrid: values of spherical harmonics at Lebedev grid points
!   x: values at grid points
!   xlm: resulting weights of spherical harmonics
    integer, intent(in) :: p, ngrid
    real(kind=8), intent(in) :: w(ngrid), vgrid((p+1)*(p+1), ngrid), x(ngrid)
    real(kind=8), intent(out) :: xlm((p+1)*(p+1))
    integer :: i
    xlm = 0
    do i = 1, ngrid
        xlm = xlm + x(i)*w(i)*vgrid(:,i)
    end do
end subroutine int_grid

! Adjoint for int_grid subroutine
subroutine int_grid_adj(p, ngrid, w, vgrid, xlm, xgrid)
! Parameters:
!   p: maximal degree of spherical harmonics
!   ngrid: number of Lebedev quadrature points
!   w: weights of quadrature points
!   vgrid: values of spherical harmonics at grid points
!   xlm: coefficients of spherical harmonics
!   xgrid: values at grid points
    integer, intent(in) :: p, ngrid
    real(kind=8), intent(in) :: w(ngrid), vgrid((p+1)*(p+1), ngrid)
    real(kind=8), intent(in) :: xlm((p+1)*(p+1))
    real(kind=8), intent(out) :: xgrid(ngrid)
    integer :: i
    xgrid = 0
    do i = 1, (p+1)*(p+1)
        xgrid = xgrid + vgrid(i,:)*xlm(i)
    end do
    xgrid = xgrid * w
end subroutine  int_grid_adj

! Integrate all external grid points with characteristic function ui
subroutine int_ui_grid(nsph, p, ngrid, w, vgrid, ui, x, xlm)
! Parameters:
!   nsph: number of spheres
!   p: maximum degree of spherical harmonics
!   ngrid: number of Lebedev grid points
!   w: weights of Lebedev grid
!   vgrid: values of spherical harmonics at Lebedev grid points
!   ui: characteristic function
!   x: values at grid points
!   xlm: resulting weights of spherical harmonics
    integer, intent(in) :: nsph, p, ngrid
    real(kind=8), intent(in) :: w(ngrid), vgrid((p+1)*(p+1), ngrid)
    real(kind=8), intent(in) :: x(ngrid, nsph), ui(ngrid, nsph)
    real(kind=8), intent(out) :: xlm((p+1)*(p+1), nsph)
    integer :: isph, igrid
    real(kind=8) :: tmp
    ! Multiply by ui and integrate from values on spheres to coefficients
    xlm = 0
    do isph = 1, nsph
        do igrid = 1, ngrid
            if (ui(igrid, isph) .ne. 0) then
                tmp = x(igrid, isph) * ui(igrid, isph)
                xlm(:, isph) = xlm(:, isph) + tmp*w(igrid)*vgrid(:, igrid)
            end if
        end do
    end do
end subroutine int_ui_grid

! Divide given cluster of spheres into two subclusters by inertial bisection
subroutine btree_node_divide(nsph, csph, n, ind, div)
! Parameters:
!   nsph: Number of all spheres
!   csph: Centers of all spheres
!   n: Number of spheres in given cluster
!   ind: Indexes of spheres in given cluster (sorted on exit)
!   div: Break point between two clusters. ind(1:div) belong to first cluster
!       and ind(div+1:n) belongs to second cluster
    integer, intent(in) :: nsph, n
    real(kind=8), intent(in) :: csph(3, nsph)
    integer, intent(inout) :: ind(n)
    integer, intent(out) :: div
    real(kind=8) :: c(3), tmp_csph(3, n), a(3, 3), w(3), work(9), scal(n)
    real(kind=8) :: alpha=1, beta=0
    integer :: i, l, r, lwork=9, info, tmp_ind(n)
    c = 0
    do i = 1, n
        c = c + csph(:, ind(i))
    end do
    c = c / n
    do i = 1, n
        tmp_csph(:, i) = csph(:, ind(i)) - c
    end do
    call dgemm('N', 'T', 3, 3, n, alpha, tmp_csph, 3, tmp_csph, 3, beta, a, 3)
    call dsyev('V', 'L', 3, a, 3, w, work, lwork, info)
    call dgemv('T', 3, n, alpha, tmp_csph, 3, a(:, 3), 1, beta, scal, 1)
    l = 1
    r = n
    do i = 1, n
        if (scal(i) .ge. 0) then
            tmp_ind(l) = ind(i)
            l = l + 1
        else
            tmp_ind(r) = ind(i)
            r = r - 1
        end if
    end do
    div = r
    ind = tmp_ind
end subroutine btree_node_divide

! Auxiliary routine for sorting
subroutine sort_by_index_partition(n, low, high, val, ind, ipiv)
    integer, intent(in) :: n, low, high
    real(kind=8), intent(in) :: val(n)
    integer, intent(inout) :: ind(n)
    integer, intent(out) :: ipiv
    integer :: i, j
    real(kind=8) :: pivot
    i = (low+high) / 2
    pivot = val(ind(i))
    i = low - 1
    j = high + 1
    do
        do
            i = i + 1
            if ((val(ind(i)) .lt. pivot) .or. (i .eq. high)) then
                exit
            end if
        end do
        do
            j = j - 1
            if ((val(ind(j)) .gt. pivot) .or. (j .eq. low)) then
                exit
            end if
        end do
        if (i .ge. j) then
            ipiv = j
            return
        end if
        ipiv = ind(i)
        ind(i) = ind(j)
        ind(j) = ipiv
    end do
end subroutine sort_by_index_partition

! Auxiliary routine for sorting
recursive subroutine sort_by_index_qsort(n, low, high, val, ind)
    integer, intent(in) :: n, low, high
    real(kind=8), intent(in) :: val(n)
    integer, intent(inout) :: ind(n)
    integer :: ipiv
    real(kind=8) :: pivot
    if (low .lt. high) then
        call sort_by_index_partition(n, low, high, val, ind, ipiv)
        call sort_by_index_qsort(n, low, ipiv-1, val, ind)
        call sort_by_index_qsort(n, ipiv+1, high, val, ind)
    end if
end subroutine sort_by_index_qsort

! Prepare binary tree (divide and compute bounding spheres)
! Number of clusters (nodes) is always 2*nsph-1
subroutine btree_init(nsph, csph, rsph, ind, cluster, children, parent, &
        cnode, rnode, snode)
! Parameters:
!   nsph: Number of all spheres
!   csph: Centers of all spheres
!   rsph: Radiuses of all spheres
!   ind: Permutation of spheres (to localize them)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   parent: parent of each cluster. 0 means no parent
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   snode: which node is leaf and contains only given sphere
    integer, intent(in) :: nsph
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph)
    integer, intent(inout) :: ind(nsph)
    integer, intent(out) :: cluster(2, 2*nsph-1), children(2, 2*nsph-1)
    integer, intent(out) :: parent(2*nsph-1)
    real(kind=8), intent(out) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    integer, intent(out) :: snode(nsph)
    integer :: i, j, n, s, e, div
    real(kind=8) :: r, r1, r2, c(3), c1(3), c2(3), d
    cluster(1, 1) = 1
    cluster(2, 1) = nsph
    parent(1) = 0
    j = 2
    ! Divide tree until leaves with single spheres inside
    do i = 1, 2*nsph-1
        s = cluster(1, i)
        e = cluster(2, i)
        n = e - s + 1
        if (n .gt. 1) then
            call btree_node_divide(nsph, csph, n, ind(s:e), div)
            cluster(1, j) = s
            cluster(2, j) = s + div - 1
            cluster(1, j+1) = s + div
            cluster(2, j+1) = e
            children(1, i) = j
            children(2, i) = j + 1
            parent(j) = i
            parent(j+1) = i
            j = j + 2
        else
            children(:, i) = 0
            snode(ind(s)) = i
        end if
    end do
    ! Compute bounding spheres
    do i = 2*nsph-1, 1, -1
        if (children(1, i) .eq. 0) then
            j = cluster(1, i)
            cnode(:, i) = csph(:, ind(j))
            rnode(i) = rsph(ind(j))
        else
            j = children(1, i)
            c1 = cnode(:, j)
            r1 = rnode(j)
            j = children(2, i)
            c2 = cnode(:, j)
            r2 = rnode(j)
            c = c1 - c2
            d = sqrt(c(1)*c(1) + c(2)*c(2) + c(3)*c(3))
            if ((r1-r2) .ge. d) then
                c = c1
                r = r1
            else if((r2-r1) .ge. d) then
                c = c2
                r = r2
            else
                r = (r1+r2+d) / 2
                c = c2 + c/d*(r-r2)
            end if
            cnode(:, i) = c
            rnode(i) = r
        end if
    end do
end subroutine btree_init

! Compute number of levels in a cluster tree
subroutine tree_get_height(n, p, h)
! Parameters:
!   n: Number of all nodes of a tree
!   p: parent of each node. 0 means no children
!   h: height of a given tree
    integer, intent(in) :: n, p(n)
    integer, intent(out) :: h
    integer :: level(n), i, j
    level(1) = 1
    h = 1
    do i = 2, n
        j = level(p(i)) + 1
        level(i) = j
        if (j .gt. h) then
            h = j
        end if
    end do
end subroutine tree_get_height

! Find near and far admissible pairs of tree nodes and store it in work array
subroutine tree_get_farnear_work(n, children, cnode, rnode, lwork, iwork, &
        & jwork, work, nnfar, nfar, nnnear, nnear)
! Parameters:
!   n: number of nodes
!   children: first and last children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   lwork: size of work array in dimension 2
!   iwork: index of current pair of nodes that needs to be checked for
!       admissibility. must be 0 for first call of this subroutine. if on exit
!       iwork is less or equal to jwork, that means lwork was too small, please
!       reallocate work array and copy all the values into new array and then
!       run procedure again.
!   jwork: amount of stored possible admissible pairs of nodes.
!       Please read iwork comments.
!   work: all the far and near pairs will be stored here
!   nnfar: total amount of far admissible pairs. valid only if iwork is
!       greater than jwork on exit.
!   nfar: amount of far admissible pairs for each node. valid only if iwork is
!       greater than jwork on exit.
!   nnnear: total amount of near admissible pairs. valid only if iwork is
!       greater than jwork on exit
!   nnear: amount of near admissible pairs for each node. valid only if iwork
!       is greater than jwork on exit
    integer, intent(in) :: n, children(2, n), lwork
    real(kind=8), intent(in) :: cnode(3, n), rnode(n)
    integer, intent(inout) :: iwork, jwork, work(3, lwork)
    integer, intent(out) :: nnfar, nfar(n), nnnear, nnear(n)
    integer :: j(2), npairs, k1, k2
    real(kind=8) :: c(3), r, d
    ! iwork is current temporary item in work array to process
    if (iwork .eq. 0) then
        work(1, 1) = 1
        work(2, 1) = 1
        iwork = 1
        jwork = 1
    end if
    ! jwork is total amount of temporary items in work array
    do while (iwork .le. jwork)
        j = work(1:2, iwork)
        c = cnode(:, j(1)) - cnode(:, j(2))
        r = rnode(j(1)) + rnode(j(2)) + max(rnode(j(1)), rnode(j(2)))
        !r = rnode(j(1)) + rnode(j(2))
        d = sqrt(c(1)*c(1) + c(2)*c(2) + c(3)*c(3))
        ! If node has no children then assume itself for purpose of finding
        ! far-field and near-filed interactions with children nodes of another
        ! node
        npairs = max(1, children(2, j(1))-children(1, j(1))+1) * &
            & max(1, children(2, j(2))-children(1, j(2))+1)
        if (d .ge. r) then
            ! Mark as far admissible pair
            !write(*,*) "FAR:", j
            work(3, iwork) = 1
        else if (npairs .eq. 1) then
            ! Mark as near admissible pair if both nodes are leaves
            !write(*,*) "NEAR:", j
            work(3, iwork) = 2
        else if (jwork+npairs .gt. lwork) then
            ! Exit procedure, since work array was too small
            !write(*,*) "SMALL LWORK"
            return
        else
            ! Mark as non-admissible pair and check all pairs of children nodes
            ! or pairs of one node (if it is a leaf node) with children of
            ! another node
            work(3, iwork) = 0
            if (children(1, j(1)) .eq. 0) then
                k1 = j(1)
                do k2 = children(1, j(2)), children(2, j(2))
                    jwork = jwork + 1
                    work(1, jwork) = k1
                    work(2, jwork) = k2
                end do
            else if(children(1, j(2)) .eq. 0) then
                k2 = j(2)
                do k1 = children(1, j(1)), children(2, j(1))
                    jwork = jwork + 1
                    work(1, jwork) = k1
                    work(2, jwork) = k2
                end do
            else
                do k1 = children(1, j(1)), children(2, j(1))
                    do k2 = children(1, j(2)), children(2, j(2))
                        jwork = jwork + 1
                        work(1, jwork) = k1
                        work(2, jwork) = k2
                    end do
                end do
            end if
            !write(*,*) "NON:", j
        end if
        iwork = iwork + 1
    end do
    nfar = 0
    nnear = 0
    do iwork = 1, jwork
        if (work(3, iwork) .eq. 1) then
            nfar(work(1, iwork)) = nfar(work(1, iwork)) + 1
        else if (work(3, iwork) .eq. 2) then
            nnear(work(1, iwork)) = nnear(work(1, iwork)) + 1
        end if
    end do
    iwork = jwork + 1
    nnfar = sum(nfar)
    nnnear = sum(nnear)
end subroutine tree_get_farnear_work

! Get near and far admissible pairs from work array of tree_get_farnear_work
! Works only for binary tree
subroutine tree_get_farnear(jwork, lwork, work, n, nnfar, nfar, sfar, far, &
        & nnnear, nnear, snear, near)
! Parameters:
!   jwork: Total number of checked pairs in work array
!   lwork: Total length of work array
!   work: Work array itself
!   n: Number of nodes
!   nnfar: Total number of all far-field interactions
!   nfar: Number of far-field interactions of each node
!   sfar: Index in far array of first far-field node for each node
!   far: Indexes of far-field nodes
!   nnnear: Total number of all near-field interactions
!   nnear: Number of near-field interactions of each node
!   snear: Index in near array of first near-field node for each node
!   near: Indexes of near-field nodes
    integer, intent(in) :: jwork, lwork, work(3, lwork), n, nnfar, nnnear
    integer, intent(in) :: nfar(n), nnear(n)
    integer, intent(out) :: sfar(n+1), far(nnfar), snear(n+1), near(nnnear)
    integer :: i, j
    integer :: cfar(n+1), cnear(n+1)
    sfar(1) = 1
    snear(1) = 1
    do i = 2, n+1
        sfar(i) = sfar(i-1) + nfar(i-1)
        snear(i) = snear(i-1) + nnear(i-1)
    end do
    cfar = sfar
    cnear = snear
    do i = 1, jwork
        if (work(3, i) .eq. 1) then
            ! Far
            j = work(1, i)
            if ((j .gt. n) .or. (j .le. 0)) then
                write(*,*) "ALARM", j
            end if
            far(cfar(j)) = work(2, i)
            cfar(j) = cfar(j) + 1
        else if (work(3, i) .eq. 2) then
            ! Near
            j = work(1, i)
            if ((j .gt. n) .or. (j .le. 0)) then
                write(*,*) "ALARM", j
            end if
            near(cnear(j)) = work(2, i)
            cnear(j) = cnear(j) + 1
        end if
    end do
end subroutine tree_get_farnear

! Transfer multipole coefficients for each node of tree
subroutine tree_m2m_baseline(nsph, nclusters, p, vscales, coef_sph, ind, &
        & cluster, children, cnode, rnode, coef_node)
! Parameters:
!   nsph: Number of all spheres
!   nclusters: Number of nodes in a tree
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   coef_sph: multipole coefficients of input spheres
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   coef_node: multipole coefficients of bounding spheres of nodes
    integer, intent(in) :: nsph, nclusters, p, ind(nsph), cluster(2, nclusters)
    real(kind=8), intent(in) :: vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: coef_sph((p+1)*(p+1), nsph)
    integer, intent(in) :: children(2, nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(out) :: coef_node((p+1)*(p+1), nclusters)
    integer :: i, j(2)
    real(kind=8) :: c1(3), c2(3), c(3), r1, r2, r
    do i = nclusters, 1, -1
        j = children(:, i)
        if (j(1) .eq. 0) then
            ! Assume that each leaf node contains single sphere with harmonics
            coef_node(:, i) = coef_sph(:, ind(cluster(1, i)))
        else
            c = cnode(:, i)
            r = rnode(i)
            c1 = cnode(:, j(1))
            r1 = rnode(j(1))
            c2 = cnode(:, j(2))
            r2 = rnode(j(2))
            coef_node(:, i) = 0
            call fmm_m2m_baseline(c1-c, r1, r, p, vscales, &
                & coef_node(:, j(1)), coef_node(:, i))
            call fmm_m2m_baseline(c2-c, r2, r, p, vscales, &
                & coef_node(:, j(2)), coef_node(:, i))
        end if
    end do
end subroutine tree_m2m_baseline

! Transfer multipole coefficients for each node of tree
subroutine tree_m2m_fast(nsph, nclusters, p, vscales, coef_sph, ind, cluster, &
        & children, cnode, rnode, coef_node)
! Parameters:
!   nsph: Number of all spheres
!   nclusters: Number of nodes in a tree
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   coef_sph: multipole coefficients of input spheres
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   coef_node: multipole coefficients of bounding spheres of nodes
    integer, intent(in) :: nsph, nclusters, p, ind(nsph), cluster(2, nclusters)
    real(kind=8), intent(in) :: vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: coef_sph((p+1)*(p+1), nsph)
    integer, intent(in) :: children(2, nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(out) :: coef_node((p+1)*(p+1), nclusters)
    integer :: i, j(2), k
    real(kind=8) :: c1(3), c(3), r1, r
    do i = nclusters, 1, -1
        j = children(:, i)
        if (j(1) .eq. 0) then
            ! In case of leaf node load input weights
            ! Assume that each leaf node contains single sphere with harmonics
            coef_node(:, i) = coef_sph(:, ind(cluster(1, i)))
        else
            ! In case of non-leaf compute weights by M2M
            c = cnode(:, i)
            r = rnode(i)
            coef_node(:, i) = 0
            do k = j(1), j(2)
                c1 = cnode(:, k)
                r1 = rnode(k)
                call fmm_m2m_fast(c1-c, r1, r, p, vscales, &
                    & coef_node(:, k), coef_node(:, i))
            end do
        end if
    end do
end subroutine tree_m2m_fast

! Compute reflection matrices for fast M2M
subroutine tree_m2m_get_mat(nclusters, children, cnode, rnode, p, vscales, &
        & reflect_mat, ztrans_mat)
! Parameters:
!   nclusters: Number of nodes of the input tree
!   children: children of each node
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   p: maximum degree of multipole basis functions
!   reflect_mat: reflection to OZ matrices
!   ztrans_mat: OZ translation matrices
    integer, intent(in) :: nclusters, children(2, nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: p
    real(kind=8), intent(in) :: vscales((p+1)*(p+1))
    real(kind=8), intent(out) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3, &
        & nclusters-1)
    real(kind=8), intent(out) :: ztrans_mat((p+1)*(p+2)*(p+3)/6, nclusters-1)
    real(kind=8) :: c(3)
    integer :: i, j(2), k
    ! Cycle over all nodes, except root one
    do i = nclusters, 1, -1
        j = children(:, i)
        ! In case of leaf node do nothing
        if (j(1) .eq. 0) then
            cycle
        end if
        do k = j(1), j(2)
            c = cnode(:, k) - cnode(:, i)
            call fmm_m2m_get_mat(c, rnode(k), rnode(i), p, vscales, &
                & reflect_mat(:, k-1), ztrans_mat(:, k-1))
        end do
    end do
end subroutine tree_m2m_get_mat

! Transfer multipole coefficients using precomputed matrices
subroutine tree_m2m_use_mat(nsph, nclusters, p, coef_sph, ind, cluster, &
        & children, cnode, rnode, reflect_mat, ztrans_mat, coef_node)
! Parameters:
!   nsph: Number of all spheres
!   nclusters: Number of nodes in a tree
!   p: maximum degree of multipole basis functions
!   coef_sph: multipole coefficients of input spheres
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   reflect_mat: reflection matrices for all node-parent transformations
!   ztrans_mat: OZ translation matrices for all node-parent transformations
!   coef_node: multipole coefficients of bounding spheres of nodes
    integer, intent(in) :: nsph, nclusters, p, ind(nsph), cluster(2, nclusters)
    real(kind=8), intent(in) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: ztrans_mat((p+1)*(p+2)*(p+3)/6, nclusters-1)
    real(kind=8), intent(in) :: coef_sph((p+1)*(p+1), nsph)
    integer, intent(in) :: children(2, nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(out) :: coef_node((p+1)*(p+1), nclusters)
    integer :: i, j(2), k
    real(kind=8) :: c1(3), c(3), r1, r
    do i = nclusters, 1, -1
        j = children(:, i)
        if (j(1) .eq. 0) then
            ! In case of leaf node load input weights
            ! Assume that each leaf node contains single sphere with harmonics
            coef_node(:, i) = coef_sph(:, ind(cluster(1, i)))
        else
            ! In case of non-leaf compute weights by M2M
            c = cnode(:, i)
            r = rnode(i)
            coef_node(:, i) = 0
            do k = j(1), j(2)
                c1 = cnode(:, k)
                r1 = rnode(k)
                call fmm_m2m_use_mat(c1-c, r1, r, p, reflect_mat(:, k-1), &
                    & ztrans_mat(:, k-1), coef_node(:, k), coef_node(:, i))
            end do
        end if
    end do
end subroutine tree_m2m_use_mat

! Adjoint transfer multipole coefficients using precomputed matrices
subroutine tree_m2m_adj_use_mat(nsph, nclusters, p, coef_sph, ind, cluster, &
        & children, cnode, rnode, reflect_mat, ztrans_mat, coef_node)
! Parameters:
!   nsph: Number of all spheres
!   nclusters: Number of nodes in a tree
!   p: maximum degree of multipole basis functions
!   coef_sph: output multipole coefficients of spheres
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   reflect_mat: reflection matrices for all node-parent transformations
!   ztrans_mat: OZ translation matrices for all node-parent transformations
!   coef_node: input multipole coefficients of bounding spheres of nodes
    integer, intent(in) :: nsph, nclusters, p, ind(nsph), cluster(2, nclusters)
    real(kind=8), intent(in) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: ztrans_mat((p+1)*(p+2)*(p+3)/6, nclusters-1)
    real(kind=8), intent(inout) :: coef_sph((p+1)*(p+1), nsph)
    integer, intent(in) :: children(2, nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(inout) :: coef_node((p+1)*(p+1), nclusters)
    integer :: i, j(2), k
    real(kind=8) :: c1(3), c(3), r1, r
    do i = 1, nclusters
        j = children(:, i)
        if (j(1) .eq. 0) then
            ! Output weights in case of a leaf node
            ! Assume that each leaf node contains single sphere with harmonics
            coef_sph(:, ind(cluster(1, i))) = coef_sph(:, ind(cluster(1, i))) + &
                & coef_node(:, i)
        else
            ! In case of non-leaf compute weights by adjoint M2M
            c = cnode(:, i)
            r = rnode(i)
            do k = j(1), j(2)
                c1 = cnode(:, k)
                r1 = rnode(k)
                call fmm_m2m_adj_use_mat(c-c1, r, r1, p, reflect_mat(:, k-1), &
                    & ztrans_mat(:, k-1), coef_node(:, i), coef_node(:, k))
            end do
        end if
    end do
end subroutine tree_m2m_adj_use_mat

! Apply M2P from entire tree to a given point
subroutine tree_m2p_treecode(c, leaf, p, vscales, nclusters, children, cnode, &
        & rnode, coef_node, v)
! Parameters:
!   c: coordinate where to compute potential
!   leaf: which leaf node contains given sphere with grid point c
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   nclusters: number of nodes in a tree
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   coef_node: multipole coefficients of bounding spheres of nodes
!   v: value of potential
    real(kind=8), intent(in) :: c(3), vscales((p+1)*(p+1)), cnode(3, nclusters)
    integer, intent(in) :: p, leaf, nclusters, children(2, nclusters)
    real(kind=8), intent(in) :: rnode(nclusters)
    real(kind=8), intent(in) :: coef_node((p+1)*(p+1), nclusters)
    real(kind=8), intent(out) :: v
    integer :: i, far(nclusters), j(2)
    real(kind=8) :: d(3), r, tmp_v
    far = 0
    ! far(i): 0 if M2P must be ignored, 1 if M2P must be applied and 2 if need
    ! to decide if M2P is applicable
    far(1) = 2
    v = 0
    do i = 1, nclusters
        ! If M2P is not applicable, ignore node
        if (far(i) .eq. 0) then
            cycle
        end if
        j = children(:, i)
        d = c - cnode(:, i)
        ! If c belongs to the origin leaf, ignore M2P
        if (leaf .eq. i) then
            far(i) = 0
            cycle
        ! Apply M2P for other leaf nodes (c is always outside for them)
        else if (j(1) .eq. 0) then
            far(i) = 1
        ! Check if point is outside sphere then apply M2P otherwise check
        ! hierarchically
        else
            r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
            ! constant 3 guarantees good convergence and accuracy
            if (r .gt. 3*rnode(i)) then
                far(i) = 1
            else
                far(i) = 0
                far(j(1)) = 2
                far(j(2)) = 2
            end if
        end if
        ! If M2P is needed
        if (far(i) .eq. 1) then
            call fmm_m2p(d, rnode(i), p, vscales, coef_node(:, i), tmp_v)
            v = v + tmp_v
        end if
    end do
end subroutine tree_m2p_treecode

! Apply matvec for ddPCM spherical harmonics by tree-code
subroutine pcm_matvec_grid_treecode(nsph, csph, rsph, ngrid, grid, w, vgrid, &
        & ui, p, vscales, ind, cluster, children, cnode, rnode, snode, &
        & coef_sph, coef_out)
! Parameters:
!   nsph: number of all spheres
!   csph: centers of all spheres
!   rsph: radiuses of all spheres
!   ngrid: number of Lebedev grid points on each sphere
!   grid: coordinates of Lebedev grid points on a unit sphere
!   ui: "outside" factor of grid points (0 is inside, 1 is "fully" outside)
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   snode: which node is leaf and contains only given sphere
!   coef_sph: multipole coefficients of bounding spheres of nodes
!   coef_out: output multipole coefficients of spherical harmonics
    integer, intent(in) :: nsph, ngrid, p, ind(nsph), cluster(2, 2*nsph-1)
    integer, intent(in) :: children(2, 2*nsph-1), snode(nsph)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((p+1)*(p+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    real(kind=8), intent(in) :: coef_sph((p+1)*(p+1), nsph)
    real(kind=8), intent(out) :: coef_out((p+1)*(p+1), nsph)
    real(kind=8) :: coef_sph_scaled((p+1)*(p+1), nsph)
    real(kind=8) :: coef_node((p+1)*(p+1), 2*nsph-1), c(3), x(ngrid)
    integer :: i, j, leaf, indi, indj
    do i = 0, p
        indi = i*i + i + 1
        do j = -i, i
            indj = indi + j
            coef_sph_scaled(indj, :) = i * coef_sph(indj, :)
        end do
    end do
    call tree_m2m_baseline(nsph, 2*nsph-1, p, vscales, coef_sph_scaled, ind, &
        & cluster, children, cnode, rnode, coef_node)
    do i = 1, nsph
        leaf = snode(i)
        do j = 1, ngrid
            if (ui(j, i) .eq. 0) then
                x(j) = 0
            else
                c = csph(:, i) + rsph(i)*grid(:,j)
                call tree_m2p_treecode(c, leaf, p, vscales, 2*nsph-1, &
                    & children, cnode, rnode, coef_node, x(j))
                x(j) = ui(j, i) * x(j)
            end if
        end do
        call int_grid(p, ngrid, w, vgrid, x, coef_out(:, i))
    end do
end subroutine pcm_matvec_grid_treecode

! Obtain local coefficients from multipole for each node of tree (M2L)
subroutine tree_m2l_baseline(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
        & pl, vscales, coef_m, coef_l)
! Parameters:
!   nclusters: Number of all clusters
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   nnfar: total number of admissible far-field pairs
!   sfar: offset for each node to get list of its far-field pairs
!   far: array that stores lists of far-field pairs
!   pm: maximum degree of multipole basis functions
!   pl: maximum degree of local basis functions
!   vscales: normalization constants for Y_lm
!   coef_m: input values of multipole expansions
!   coef_l: output values of local expansions
    integer, intent(in) :: nclusters, nnfar, sfar(nclusters+1), far(nnfar)
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: coef_m((pm+1)*(pm+1), nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(out) :: coef_l((pl+1)*(pl+1), nclusters)
    integer :: i, j, k
    real(kind=8) :: c(3), r1, r2
    do i = 1, nclusters
        coef_l(:, i) = 0
        do j = sfar(i), sfar(i+1)-1
            k = far(j)
            c = cnode(:, k) - cnode(:, i)
            r1 = rnode(k)
            r2 = rnode(i)
            call fmm_m2l_baseline(c, r1, r2, pm, pl, vscales, coef_m(:, k), &
                & coef_l(:, i))
        end do
    end do
end subroutine tree_m2l_baseline

! Obtain local coefficients from multipole for each node of tree (M2L)
subroutine tree_m2l_fast(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
        & pl, vscales, coef_m, coef_l)
! Parameters:
!   nclusters: Number of all clusters
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   nnfar: total number of admissible far-field pairs
!   sfar: offset for each node to get list of its far-field pairs
!   far: array that stores lists of far-field pairs
!   pm: maximum degree of multipole basis functions
!   pl: maximum degree of local basis functions
!   vscales: normalization constants for Y_lm
!   coef_m: input values of multipole expansions
!   coef_l: output values of local expansions
    integer, intent(in) :: nclusters, nnfar, sfar(nclusters+1), far(nnfar)
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: coef_m((pm+1)*(pm+1), nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(out) :: coef_l((pl+1)*(pl+1), nclusters)
    integer :: i, j, k
    real(kind=8) :: c(3), r1, r2
    do i = 1, nclusters
        coef_l(:, i) = 0
        do j = sfar(i), sfar(i+1)-1
            k = far(j)
            c = cnode(:, k) - cnode(:, i)
            r1 = rnode(k)
            r2 = rnode(i)
            call fmm_m2l_fast(c, r1, r2, pm, pl, vscales, coef_m(:, k), &
                & coef_l(:, i))
        end do
    end do
end subroutine tree_m2l_fast

! Compute reflection matrices for fast M2L
subroutine tree_m2l_get_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
        & pl, vscales, reflect_mat, ztrans_mat)
! Parameters:
!   nclusters: Number of all clusters
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   nnfar: total number of admissible far-field pairs
!   sfar: offset for each node to get list of its far-field pairs
!   far: array that stores lists of far-field pairs
!   pm: maximum degree of multipole basis functions
!   pl: maximum degree of local basis functions
!   vscales: normalization constants for Y_lm
!   reflect_mat: reflection to OZ matrices
!   ztrans_mat: OZ translation matrices
    integer, intent(in) :: nclusters, nnfar, sfar(nclusters+1), far(nnfar)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(out) :: reflect_mat((max(pm,pl)+1)*(2*max(pm,pl)+1) &
        & *(2*max(pm,pl)+3)/3, nnfar)
    real(kind=8), intent(out) :: ztrans_mat((min(pm,pl)+1)*(min(pm,pl)+2) &
        & *(3*max(pm,pl)+3-min(pm,pl))/6, nnfar)
    real(kind=8) :: c(3)
    integer :: i, j, k
    ! Cycle over all nodes, except root one
    do i = 1, nclusters
        do j = sfar(i), sfar(i+1)-1
            k = far(j)
            c = cnode(:, k) - cnode(:, i)
            call fmm_m2l_get_mat(c, rnode(k), rnode(i), pm, pl, vscales, &
                & reflect_mat(:, j), ztrans_mat(:, j))
        end do
    end do
end subroutine tree_m2l_get_mat

! Transfer multipole coefficients into local using precomputed matrices
subroutine tree_m2l_use_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
        & pl, reflect_mat, ztrans_mat, coef_m, coef_l)
! Parameters:
!   nclusters: Number of all clusters
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   nnfar: total number of admissible far-field pairs
!   sfar: offset for each node to get list of its far-field pairs
!   far: array that stores lists of far-field pairs
!   reflect_mat: reflection matrices for all node-parent transformations
!   ztrans_mat: OZ translation matrices for all node-parent transformations
!   coef_m: input values of multipole expansions
!   coef_l: output values of local expansions
    integer, intent(in) :: nclusters, nnfar, sfar(nclusters+1), far(nnfar)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: reflect_mat((max(pm,pl)+1)*(2*max(pm,pl)+1) &
        & *(2*max(pm,pl)+3)/3, nnfar)
    real(kind=8), intent(in) :: ztrans_mat((min(pm,pl)+1)*(min(pm,pl)+2) &
        & *(3*max(pm,pl)+3-min(pm,pl))/6, nnfar)
    real(kind=8), intent(in) :: coef_m((pm+1)*(pm+1), nclusters)
    real(kind=8), intent(out) :: coef_l((pl+1)*(pl+1), nclusters)
    integer :: i, j, k
    real(kind=8) :: c(3)
    do i = 1, nclusters
        coef_l(:, i) = 0
        do j = sfar(i), sfar(i+1)-1
            k = far(j)
            c = cnode(:, k) - cnode(:, i)
            call fmm_m2l_use_mat(c, rnode(k), rnode(i), pm, pl, &
                & reflect_mat(:, j), ztrans_mat(:, j), coef_m(:, k), &
                & coef_l(:, i))
        end do
    end do
end subroutine tree_m2l_use_mat

! Adjoint transfer multipole coefficients into local using precomputed matrices
subroutine tree_m2l_adj_use_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
        & pl, reflect_mat, ztrans_mat, coef_l, coef_m)
! Parameters:
!   nclusters: Number of all clusters
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   nnfar: total number of admissible far-field pairs
!   sfar: offset for each node to get list of its far-field pairs
!   far: array that stores lists of far-field pairs
!   reflect_mat: reflection matrices for all node-parent transformations
!   ztrans_mat: OZ translation matrices for all node-parent transformations
!   coef_l: input values of local expansions
!   coef_m: output values of multipole expansions
    integer, intent(in) :: nclusters, nnfar, sfar(nclusters+1), far(nnfar)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: pm, pl
    real(kind=8), intent(in) :: reflect_mat((max(pm,pl)+1)*(2*max(pm,pl)+1) &
        & *(2*max(pm,pl)+3)/3, nnfar)
    real(kind=8), intent(in) :: ztrans_mat((min(pm,pl)+1)*(min(pm,pl)+2) &
        & *(3*max(pm,pl)+3-min(pm,pl))/6, nnfar)
    real(kind=8), intent(in) :: coef_l((pl+1)*(pl+1), nclusters)
    real(kind=8), intent(out) :: coef_m((pm+1)*(pm+1), nclusters)
    integer :: i, j, k
    real(kind=8) :: c(3)
    coef_m = 0
    do i = 1, nclusters
        do j = sfar(i), sfar(i+1)-1
            k = far(j)
            c = cnode(:, k) - cnode(:, i)
            call fmm_m2l_adj_use_mat(c, rnode(k), rnode(i), pm, pl, &
                & reflect_mat(:, j), ztrans_mat(:, j), coef_l(:, i), &
                & coef_m(:, k))
        end do
    end do
end subroutine tree_m2l_adj_use_mat

! Transfer local coefficients for each node of tree
subroutine tree_l2l_baseline(nsph, nclusters, p, vscales, coef_node, ind, &
        & cluster, children, cnode, rnode, coef_sph)
! Parameters:
!   nsph: Number of all spheres
!   nclusters: Number of nodes in a tree
!   p: maximum degree of local basis functions
!   vscales: normalization constants for Y_lm
!   coef_node: local coefficients of bounding spherical harmonics of each node
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   coef_sph: local coefficients of output spherical hamonics
    integer, intent(in) :: nsph, p, nclusters, ind(nsph), cluster(2, nclusters)
    integer, intent(in) :: children(2, nclusters)
    real(kind=8), intent(in) :: vscales((p+1)*(p+1)), cnode(3, nclusters)
    real(kind=8), intent(inout) :: coef_node((p+1)*(p+1), nclusters)
    real(kind=8), intent(in) :: rnode(nclusters)
    real(kind=8), intent(out) :: coef_sph((p+1)*(p+1), nsph)
    integer :: i, j(2)
    real(kind=8) :: c1(3), c2(3), c(3), r1, r2, r
    do i = 1, nclusters
        j = children(:, i)
        if (j(1) .eq. 0) then
            ! Assume that each leaf node contains single sphere with harmonics
            coef_sph(:, ind(cluster(1, i))) = coef_node(:, i)
        else
            c = cnode(:, i)
            r = rnode(i)
            c1 = cnode(:, j(1))
            r1 = rnode(j(1))
            c2 = cnode(:, j(2))
            r2 = rnode(j(2))
            call fmm_l2l_baseline(c-c1, r, r1, p, vscales, &
                & coef_node(:, i), coef_node(:, j(1)))
            call fmm_l2l_baseline(c-c2, r, r2, p, vscales, &
                & coef_node(:, i), coef_node(:, j(2)))
        end if
    end do
end subroutine tree_l2l_baseline

! Transfer local coefficients for each node of tree
subroutine tree_l2l_fast(nsph, nclusters, p, vscales, coef_node, ind, &
        & cluster, children, cnode, rnode, coef_sph)
! Parameters:
!   nsph: Number of all spheres
!   nclusters: Number of all nodes in a tree
!   p: maximum degree of local basis functions
!   vscales: normalization constants for Y_lm
!   coef_node: local coefficients of bounding spherical harmonics of each node
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   coef_sph: local coefficients of output spherical hamonics
    integer, intent(in) :: nsph, nclusters, p, ind(nsph), cluster(2, nclusters)
    integer, intent(in) :: children(2, nclusters)
    real(kind=8), intent(in) :: vscales((p+1)*(p+1)), cnode(3, nclusters)
    real(kind=8), intent(inout) :: coef_node((p+1)*(p+1), nclusters)
    real(kind=8), intent(in) :: rnode(nclusters)
    real(kind=8), intent(out) :: coef_sph((p+1)*(p+1), nsph)
    integer :: i, j(2), k
    real(kind=8) :: c1(3), c(3), r1, r
    do i = 1, nclusters
        j = children(:, i)
        if (j(1) .eq. 0) then
            ! Output weights in case of a leaf node
            coef_sph(:, ind(cluster(1, i))) = coef_node(:, i)
        else
            ! Do L2L in case of non-leaf node
            c = cnode(:, i)
            r = rnode(i)
            do k = j(1), j(2)
                c1 = cnode(:, k)
                r1 = rnode(k)
                call fmm_l2l_fast(c-c1, r, r1, p, vscales, &
                    & coef_node(:, i), coef_node(:, k))
            end do
        end if
    end do
end subroutine tree_l2l_fast

! Compute reflection matrices for fast L2L
subroutine tree_l2l_get_mat(nclusters, children, cnode, rnode, p, vscales, &
        & reflect_mat, ztrans_mat)
! Parameters:
!   nclusters: Number of nodes of the input tree
!   children: children of each node
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   p: maximum degree of multipole basis functions
!   reflect_mat: reflection to OZ matrices
!   ztrans_mat: OZ translation matrices
    integer, intent(in) :: nclusters, children(2, nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: p
    real(kind=8), intent(in) :: vscales((p+1)*(p+1))
    real(kind=8), intent(out) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3, &
        & nclusters-1)
    real(kind=8), intent(out) :: ztrans_mat((p+1)*(p+2)*(p+3)/6, nclusters-1)
    real(kind=8) :: c(3)
    integer :: i, j(2), k
    ! Cycle over all nodes, except root one
    do i = 1, nclusters
        j = children(:, i)
        ! In case of leaf node do nothing
        if (j(1) .eq. 0) then
            cycle
        end if
        do k = j(1), j(2)
            c = cnode(:, i) - cnode(:, k)
            call fmm_l2l_get_mat(c, rnode(i), rnode(k), p, vscales, &
                & reflect_mat(:, k-1), ztrans_mat(:, k-1))
        end do
    end do
end subroutine tree_l2l_get_mat

! Transfer local coefficients using precomputed matrices
subroutine tree_l2l_use_mat(nsph, nclusters, p, coef_node, ind, cluster, &
        & children, cnode, rnode, reflect_mat, ztrans_mat, coef_sph)
! Parameters:
!   nsph: Number of all spheres
!   nclusters: Number of nodes in a tree
!   p: maximum degree of multipole basis functions
!   coef_node: local coefficients of bounding spheres of nodes
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   reflect_mat: reflection matrices for all node-parent transformations
!   ztrans_mat: OZ translation matrices for all node-parent transformations
!   coef_sph: local coefficients of output spheres
    integer, intent(in) :: nsph, nclusters, p, ind(nsph), cluster(2, nclusters)
    real(kind=8), intent(in) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: ztrans_mat((p+1)*(p+2)*(p+3)/6, nclusters-1)
    real(kind=8), intent(inout) :: coef_node((p+1)*(p+1), nclusters)
    integer, intent(in) :: children(2, nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(out) :: coef_sph((p+1)*(p+1), nsph)
    integer :: i, j(2), k
    real(kind=8) :: c1(3), c(3), r1, r
    do i = 1, nclusters
        j = children(:, i)
        if (j(1) .eq. 0) then
            ! Output weights in case of a leaf node
            ! Assume that each leaf node contains single sphere with harmonics
            coef_sph(:, ind(cluster(1, i))) = coef_node(:, i)
        else
            ! In case of non-leaf compute weights by L2L
            c = cnode(:, i)
            r = rnode(i)
            do k = j(1), j(2)
                c1 = cnode(:, k)
                r1 = rnode(k)
                call fmm_l2l_use_mat(c-c1, r, r1, p, reflect_mat(:, k-1), &
                    & ztrans_mat(:, k-1), coef_node(:, i), coef_node(:, k))
            end do
        end if
    end do
end subroutine tree_l2l_use_mat

! Adjoint transfer local coefficients using precomputed matrices
subroutine tree_l2l_adj_use_mat(nsph, nclusters, p, coef_node, ind, cluster, &
        & children, cnode, rnode, reflect_mat, ztrans_mat, coef_sph)
! Parameters:
!   nsph: Number of all spheres
!   nclusters: Number of nodes in a tree
!   p: maximum degree of multipole basis functions
!   coef_node: output local coefficients of bounding spheres of nodes
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   reflect_mat: reflection matrices for all node-parent transformations
!   ztrans_mat: OZ translation matrices for all node-parent transformations
!   coef_sph: local coefficients of input spheres
    integer, intent(in) :: nsph, nclusters, p, ind(nsph), cluster(2, nclusters)
    real(kind=8), intent(in) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: ztrans_mat((p+1)*(p+2)*(p+3)/6, nclusters-1)
    real(kind=8), intent(out) :: coef_node((p+1)*(p+1), nclusters)
    integer, intent(in) :: children(2, nclusters)
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(in) :: coef_sph((p+1)*(p+1), nsph)
    integer :: i, j(2), k
    real(kind=8) :: c1(3), c(3), r1, r
    do i = nclusters, 1, -1
        j = children(:, i)
        if (j(1) .eq. 0) then
            ! In case of leaf node load input weights
            ! Assume that each leaf node contains single sphere with harmonics
            coef_node(:, i) = coef_sph(:, ind(cluster(1, i)))
        else
            ! In case of non-leaf compute weights by adjoint L2L
            c = cnode(:, i)
            r = rnode(i)
            coef_node(:, i) = 0
            do k = j(1), j(2)
                c1 = cnode(:, k)
                r1 = rnode(k)
                call fmm_l2l_adj_use_mat(c1-c, r1, r, p, reflect_mat(:, k-1), &
                    & ztrans_mat(:, k-1), coef_node(:, k), coef_node(:, i))
            end do
        end if
    end do
end subroutine tree_l2l_adj_use_mat

! Apply L2P from local spherical harmonics to grid points and add near M2P
subroutine tree_l2p_m2p_fmm(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, coef_sph_m, coef_sph_l, nclusters, nnnear, &
        & snear, near, cluster, ind, ui, x)
    integer, intent(in) :: nsph, ngrid, pm, pl, nclusters, nnnear
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    integer, intent(in) :: cluster(2, nclusters), ind(nsph)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: vgrid((pl+1)*(pl+1), ngrid), ui(ngrid, nsph)
    real(kind=8), intent(in) :: coef_sph_m((pm+1)*(pm+1), nsph), w(ngrid)
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: coef_sph_l((pl+1)*(pl+1), nsph)
    real(kind=8), intent(out) :: x(ngrid, nsph)
    integer :: isph, i, j, inode, inear, jnode, isph_node, jsph_node, jsph
    integer :: igrid, indi
    real(kind=8) :: c(3), tmp_v, start, finish
    call cpu_time(start)
    x = 0
    ! Apply far-field L2P (from each sphere to its own grid points)
    do isph = 1, nsph
        x(:, isph) = 0
        do i = 0, pl
            indi = i*i + i + 1
            do j = indi-i, indi+i
                x(:, isph) = x(:, isph) + &
                    & coef_sph_l(j, isph)*vgrid(j, :)/vscales(indi)**2
            end do
        end do
    end do
    ! Apply near-field M2P (from input spheres to output grid points)
    do inode = 1, nclusters
        do inear = snear(inode), snear(inode+1)-1
            jnode = near(inear)
            do isph_node = cluster(1, inode), cluster(2, inode)
                isph = ind(isph_node)
                do igrid = 1, ngrid
                    if (ui(igrid, isph) .eq. 0) then
                        cycle
                    end if
                    do jsph_node = cluster(1, jnode), cluster(2, jnode)
                        jsph = ind(jsph_node)
                        if (isph .eq. jsph) then
                            cycle
                        end if
                        c = csph(:, isph) + rsph(isph)*grid(:, igrid)
                        call fmm_m2p(c-csph(:, jsph), rsph(jsph), pm, &
                            & vscales, coef_sph_m(:, jsph), tmp_v)
                        x(igrid, isph) = x(igrid, isph) + tmp_v
                    end do
                end do
            end do
        end do
    end do
    call cpu_time(finish)
end subroutine tree_l2p_m2p_fmm

! Find size of external grid for each sphere
subroutine get_ngrid_ext(nsph, ngrid, ui, ngrid_ext, ngrid_ext_sph)
    integer, intent(in) :: nsph, ngrid
    real(kind=8), intent(in) :: ui(ngrid, nsph)
    integer, intent(out) :: ngrid_ext, ngrid_ext_sph(nsph)
    integer :: i, j
    ngrid_ext = 0
    do i = 1, nsph
        ngrid_ext_sph(i) = 0
        do j = 1, ngrid
            if (ui(j, i) .ne. 0) then
                ngrid_ext_sph(i) = ngrid_ext_sph(i) + 1
                ngrid_ext = ngrid_ext + 1
            end if
        end do
    end do
end subroutine get_ngrid_ext

! Get number of all near-field sphere-to-external-grid-point pairs
subroutine get_ngrid_ext_near(nsph, ngrid_ext_sph, nclusters, nnear, snode, &
        & ngrid_ext_near)
    integer, intent(in) :: nsph, ngrid_ext_sph(nsph), nclusters, snode(nsph)
    integer, intent(in) :: nnear(nclusters)
    integer, intent(out) :: ngrid_ext_near
    integer :: isph, inode
    ngrid_ext_near = 0
    do isph = 1, nsph
        inode = snode(isph)
        ngrid_ext_near = ngrid_ext_near + ngrid_ext_sph(isph)*nnear(inode)
    end do
end subroutine get_ngrid_ext_near

! Find indexes of external grid points
! Output is stored in CSR (compressed sparse row) format, sometimes denotes by
! arrays A, IA and JA. grid_ext_ia with grid_ext_ja represent sparsity pattern
! as row and column indexes
subroutine get_grid_ext_ind(nsph, ngrid, ui, ngrid_ext, ngrid_ext_sph, &
        & grid_ext_ia, grid_ext_ja)
    integer, intent(in) :: nsph, ngrid, ngrid_ext, ngrid_ext_sph(nsph)
    real(kind=8), intent(in) :: ui(ngrid, nsph)
    integer, intent(out) :: grid_ext_ia(nsph+1), grid_ext_ja(ngrid_ext)
    integer :: i, j, k
    grid_ext_ia(1) = 1
    k = 1
    do i = 1, nsph
        grid_ext_ia(i+1) = grid_ext_ia(i) + ngrid_ext_sph(i)
        do j = 1, ngrid
            if (ui(j, i) .ne. 0) then
                grid_ext_ja(k) = j
                k = k + 1
            end if
        end do
    end do
end subroutine get_grid_ext_ind

! Save far-field L2P and near-field M2P, where P are external grid points
subroutine tree_l2p_m2p_get_mat(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, nclusters, nnnear, snear, near, cluster, snode, &
        & ind, ui, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
        & ngrid_ext_near, l2p_mat, m2p_mat)
!   snode: which node is leaf and contains only given sphere
    integer, intent(in) :: nsph, ngrid, pm, pl, nclusters, nnnear
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    integer, intent(in) :: cluster(2, nclusters), snode(nsph), ind(nsph)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: vgrid((pl+1)*(pl+1), ngrid), ui(ngrid, nsph)
    real(kind=8), intent(in) :: w(ngrid)
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    integer, intent(in) :: ngrid_ext, ngrid_ext_near, ngrid_ext_sph(nsph+1)
    integer, intent(in) :: grid_ext_ia(nsph+1), grid_ext_ja(ngrid_ext)
    real(kind=8), intent(out) :: l2p_mat((pl+1)*(pl+1), ngrid_ext)
    real(kind=8), intent(out) :: m2p_mat((pm+1)*(pm+1), ngrid_ext_near)
    integer :: isph, inode, inear, jnode, jsph, i, j, igrid_ext_sph, igrid_sph
    integer :: igrid_ext, igrid_ext_near, indi
    real(kind=8) :: x(ngrid, nsph), c(3), tmp_v
    ! Get far-field L2P matrices (from each sphere to its own grid points)
    do isph = 1, nsph
        do igrid_ext = grid_ext_ia(isph), grid_ext_ia(isph+1)-1
            igrid_sph = grid_ext_ja(igrid_ext)
            do i = 0, pl
                indi = i*i + i + 1
                l2p_mat(indi-i:indi+i, igrid_ext) = &
                    & vgrid(indi-i:indi+i, igrid_sph) / vscales(indi)**2
            end do
        end do
    end do
    ! Get near-field M2P matrices (to all external grid points of given sphere
    ! from all near-field spheres)
    igrid_ext_near = 1
    do isph = 1, nsph
        inode = snode(isph)
        do inear = snear(inode), snear(inode+1)-1
            jnode = near(inear)
            ! Assume near-field interaction is only possible between leaves
            jsph = ind(cluster(1, jnode))
            ! Self-interaction does not mean anything in this problem
            if (isph .eq. jsph) then
                i = igrid_ext_near
                j = i + ngrid_ext_sph(isph) - 1
                m2p_mat(:, i:j) = 0
                cycle
            end if
            ! M2P marices from spherical harmonics of jsph sphere to all
            ! external grid points of isph sphere
            do igrid_ext_sph = 1, ngrid_ext_sph(isph)
                igrid_sph = grid_ext_ja(grid_ext_ia(isph) + igrid_ext_sph - 1)
                c = csph(:, isph) + rsph(isph)*grid(:, igrid_sph)
                call fmm_m2p_mat(c-csph(:, jsph), rsph(jsph), pm, vscales, &
                    & m2p_mat(:, igrid_ext_sph+igrid_ext_near-1))
            end do
            igrid_ext_near = igrid_ext_near + ngrid_ext_sph(isph)
        end do
    end do
end subroutine tree_l2p_m2p_get_mat

! Use precomputed far-field L2P and near-field M2P
subroutine tree_l2p_m2p_use_mat(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, nclusters, nnnear, snear, near, cluster, snode, &
        & ind, ui, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
        & ngrid_ext_near, l2p_mat, m2p_mat, coef_sph_m, coef_sph_l, x)
!   snode: which node is leaf and contains only given sphere
    integer, intent(in) :: nsph, ngrid, pm, pl, nclusters, nnnear
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    integer, intent(in) :: cluster(2, nclusters), snode(nsph), ind(nsph)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: vgrid((pl+1)*(pl+1), ngrid), ui(ngrid, nsph)
    real(kind=8), intent(in) :: w(ngrid)
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    integer, intent(in) :: ngrid_ext, ngrid_ext_near, ngrid_ext_sph(nsph+1)
    integer, intent(in) :: grid_ext_ia(nsph+1), grid_ext_ja(ngrid_ext)
    real(kind=8), intent(in) :: l2p_mat((pl+1)*(pl+1), ngrid_ext)
    real(kind=8), intent(in) :: m2p_mat((pm+1)*(pm+1), ngrid_ext_near)
    real(kind=8), intent(in) :: coef_sph_m((pm+1)*(pm+1), nsph)
    real(kind=8), intent(in) :: coef_sph_l((pl+1)*(pl+1), nsph)
    real(kind=8), intent(out) :: x(ngrid, nsph)
    integer :: isph, inode, inear, jnode, jsph, i, j, igrid_ext_sph, igrid_sph
    integer :: igrid_ext, igrid_ext_near
    real(kind=8) :: c(3), tmp_v, y(ngrid)
    ! Get far-field L2P matrices (from each sphere to its own grid points)
    do isph = 1, nsph
        call dgemv('T', (pl+1)*(pl+1), ngrid_ext_sph(isph), 1.0d0, &
            & l2p_mat(1, grid_ext_ia(isph)), (pl+1)*(pl+1), &
            & coef_sph_l(1, isph), 1, 0.0d0, x(1, isph), 1)
    end do
    ! Get near-field M2P matrices (to all external grid points of given sphere
    ! from all near-field spheres)
    igrid_ext_near = 1
    do isph = 1, nsph
        inode = snode(isph)
        do inear = snear(inode), snear(inode+1)-1
            jnode = near(inear)
            ! Assume near-field interaction is only possible between leaves
            jsph = ind(cluster(1, jnode))
            ! Self-interaction does not mean anything in this problem
            if (isph .eq. jsph) then
                cycle
            end if
            ! Apply M2P from spherical harmonics of jsph sphere to all
            ! external grid points of isph sphere
            call dgemv('T', (pm+1)*(pm+1), ngrid_ext_sph(isph), 1.0d0, &
                & m2p_mat(1, igrid_ext_near), (pm+1)*(pm+1), &
                & coef_sph_m(1, jsph), 1, 1.0d0, x(1, isph), 1)
            igrid_ext_near = igrid_ext_near + ngrid_ext_sph(isph)
        end do
        y = 0
        do igrid_ext_sph = 1, ngrid_ext_sph(isph)
            igrid_sph = grid_ext_ja(grid_ext_ia(isph) + igrid_ext_sph - 1)
            y(igrid_sph) = x(igrid_ext_sph, isph)
        end do
        x(:, isph) = y
    end do
end subroutine tree_l2p_m2p_use_mat

! Adjoint for the L2P and M2P operations using stored matrices
subroutine tree_l2p_m2p_adj_use_mat(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, nclusters, nnnear, snear, near, cluster, snode, &
        & ind, ui, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
        & ngrid_ext_near, l2p_mat, m2p_mat, xgrid, coef_sph_m, coef_sph_l)
!   snode: which node is leaf and contains only given sphere
    integer, intent(in) :: nsph, ngrid, pm, pl, nclusters, nnnear
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    integer, intent(in) :: cluster(2, nclusters), snode(nsph), ind(nsph)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: vgrid((pl+1)*(pl+1), ngrid), ui(ngrid, nsph)
    real(kind=8), intent(in) :: w(ngrid)
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    integer, intent(in) :: ngrid_ext, ngrid_ext_near, ngrid_ext_sph(nsph+1)
    integer, intent(in) :: grid_ext_ia(nsph+1), grid_ext_ja(ngrid_ext)
    real(kind=8), intent(in) :: l2p_mat((pl+1)*(pl+1), ngrid_ext)
    real(kind=8), intent(in) :: m2p_mat((pm+1)*(pm+1), ngrid_ext_near)
    real(kind=8), intent(in) :: xgrid(ngrid, nsph)
    real(kind=8), intent(out) :: coef_sph_m((pm+1)*(pm+1), nsph)
    real(kind=8), intent(out) :: coef_sph_l((pl+1)*(pl+1), nsph)
    integer :: isph, inode, inear, jnode, jsph, i, j, igrid_ext_sph, igrid_sph
    integer :: igrid_ext, igrid_ext_near
    real(kind=8) :: x_ext(ngrid), c(3), tmp_v, y(ngrid)
    ! Use far-field L2P matrices and near-field M2P matrices (to all external grid
    ! points of given sphere)
    igrid_ext_near = 1
    coef_sph_m = 0
    coef_sph_l = 0
    x_ext = 0
    do isph = 1, nsph
        inode = snode(isph)
        ! Copy values from external grid points in a compact form
        do igrid_ext_sph = 1, ngrid_ext_sph(isph)
            igrid_sph = grid_ext_ja(grid_ext_ia(isph) + igrid_ext_sph - 1)
            x_ext(igrid_ext_sph) = xgrid(igrid_sph, isph) * ui(igrid_sph, isph)
        end do
        ! Use far-field L2P matrices (from each sphere to its own grid
        ! points)
        call dgemv('N', (pl+1)*(pl+1), ngrid_ext_sph(isph), 1.0d0, &
            & l2p_mat(1, grid_ext_ia(isph)), (pl+1)*(pl+1), &
            & x_ext, 1, 0.0d0, coef_sph_l(1, isph), 1)
        ! Use near-field M2P matrices
        do inear = snear(inode), snear(inode+1)-1
            jnode = near(inear)
            ! Assume near-field interaction is only possible between leaves
            jsph = ind(cluster(1, jnode))
            ! Self-interaction does not mean anything in this problem
            if (isph .eq. jsph) then
                cycle
            end if
            ! Actual adjoint operator from external grid points of isph sphere
            ! to spherical harmonics of jsph sphere
            call dgemv('N', (pm+1)*(pm+1), ngrid_ext_sph(isph), 1.0d0, &
                & m2p_mat(1, igrid_ext_near), (pm+1)*(pm+1), &
                & x_ext, 1, 1.0d0, coef_sph_m(1, jsph), 1)
            igrid_ext_near = igrid_ext_near + ngrid_ext_sph(isph)
        end do
    end do
end subroutine tree_l2p_m2p_adj_use_mat

! Apply matvec for ddPCM spherical harmonics by FMM
! Baseline in terms of p^4 operations
subroutine pcm_matvec_grid_fmm_baseline(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, coef_sph, &
        & coef_out)
    integer, intent(in) :: nsph, ngrid, pm, pl, ind(nsph), nclusters
    integer, intent(in) :: cluster(2, nclusters), children(2, nclusters), nnfar
    integer, intent(in) :: sfar(nclusters+1), far(nnfar), nnnear
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(in) :: coef_sph((pm+1)*(pm+1), nsph)
    real(kind=8), intent(out) :: coef_out((pl+1)*(pl+1), nsph)
    real(kind=8) :: coef_sph_scaled((pm+1)*(pm+1), nsph)
    real(kind=8) :: coef_node_m((pm+1)*(pm+1), nclusters)
    real(kind=8) :: coef_node_l((pl+1)*(pl+1), nclusters)
    real(kind=8) :: x(ngrid, nsph)
    integer :: i, j, indi, indj
    do i = 0, pm
        indi = i*i + i + 1
        do j = -i, i
            indj = indi + j
            coef_sph_scaled(indj, :) = i * coef_sph(indj, :)
        end do
    end do
    call tree_m2m_baseline(nsph, nclusters, pm, vscales, coef_sph_scaled, &
        & ind, cluster, children, cnode, rnode, coef_node_m)
    call tree_m2l_baseline(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
        & pl, vscales, coef_node_m, coef_node_l)
    call tree_l2l_baseline(nsph, nclusters, pl, vscales, coef_node_l, ind, &
        & cluster, children, cnode, rnode, coef_out)
    call tree_l2p_m2p_fmm(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, coef_sph_scaled, coef_out, nclusters, nnnear, &
        & snear, near, cluster, ind, ui, x)
    call int_ui_grid(nsph, pl, ngrid, w, vgrid, ui, x, coef_out)
end subroutine pcm_matvec_grid_fmm_baseline

! Apply matvec for ddPCM spherical harmonics by FMM
! Optimized in terms of p^3 operations
subroutine pcm_matvec_grid_fmm_fast(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, coef_sph, &
        & coef_out)
    integer, intent(in) :: nsph, ngrid, pm, pl, ind(nsph), nclusters
    integer, intent(in) :: cluster(2, nclusters), children(2, nclusters), nnfar
    integer, intent(in) :: sfar(nclusters+1), far(nnfar), nnnear
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    real(kind=8), intent(in) :: coef_sph((pm+1)*(pm+1), nsph)
    real(kind=8), intent(out) :: coef_out((pl+1)*(pl+1), nsph)
    real(kind=8) :: coef_sph_scaled((pm+1)*(pm+1), nsph)
    real(kind=8) :: coef_node_m((pm+1)*(pm+1), nclusters)
    real(kind=8) :: coef_node_l((pl+1)*(pl+1), nclusters)
    real(kind=8) :: x(ngrid, nsph)
    integer :: i, j, indi, indj
    do i = 0, pm
        indi = i*i + i + 1
        do j = -i, i
            indj = indi + j
            coef_sph_scaled(indj, :) = i * coef_sph(indj, :)
        end do
    end do
    call tree_m2m_fast(nsph, nclusters, pm, vscales, coef_sph_scaled, ind, &
        & cluster, children, cnode, rnode, coef_node_m)
    call tree_m2l_fast(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
        & pl, vscales, coef_node_m, coef_node_l)
    call tree_l2l_fast(nsph, nclusters, pl, vscales, coef_node_l, ind, &
        & cluster, children, cnode, rnode, coef_out)
    call tree_l2p_m2p_fmm(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, coef_sph_scaled, coef_out, nclusters, nnnear, &
        & snear, near, cluster, ind, ui, x)
    call int_ui_grid(nsph, pl, ngrid, w, vgrid, ui, x, coef_out)
end subroutine pcm_matvec_grid_fmm_fast

! Apply matvec for ddPCM spherical harmonics by FMM
! Computes and uses reflection and OZ translation matrices for M2M
subroutine pcm_matvec_grid_fmm_test_mat(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, snode, &
        & children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
        & ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, ngrid_ext_near, &
        & coef_sph, coef_out)
    integer, intent(in) :: nsph, ngrid, pm, pl, ind(nsph), nclusters
    integer, intent(in) :: cluster(2, nclusters), children(2, nclusters), nnfar
    integer, intent(in) :: sfar(nclusters+1), far(nnfar), nnnear, snode(nsph)
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: ngrid_ext, ngrid_ext_sph(nsph), grid_ext_ia(nsph+1)
    integer, intent(in) :: grid_ext_ja(ngrid_ext), ngrid_ext_near
    real(kind=8), intent(in) :: coef_sph((pm+1)*(pm+1), nsph)
    real(kind=8), intent(out) :: coef_out((pl+1)*(pl+1), nsph)
    real(kind=8) :: x(ngrid, nsph)
    real(kind=8) :: coef_sph_scaled((pm+1)*(pm+1), nsph)
    real(kind=8) :: coef_node_m((pm+1)*(pm+1), nclusters)
    real(kind=8) :: coef_node_l((pl+1)*(pl+1), nclusters)
    real(kind=8) :: m2m_l2l_reflect_mat((max(pm,pl)+1)*(2*max(pm,pl)+1) &
        & *(2*max(pm,pl)+3)/3, nclusters-1)
    real(kind=8) :: m2m_l2l_ztrans_mat((max(pm,pl)+1)*(max(pm,pl)+2) &
        & *(max(pm,pl)+3)/6, nclusters-1)
    real(kind=8) :: m2l_reflect_mat((max(pm,pl)+1)*(2*max(pm,pl)+1) &
        & *(2*max(pm,pl)+3)/3, nnfar)
    real(kind=8) :: m2l_ztrans_mat((min(pm,pl)+1)*(min(pm,pl)+2) &
        & *(3*max(pm,pl)+3-min(pm,pl))/6, nnfar)
    real(kind=8) :: l2p_mat((pl+1)*(pl+1), ngrid_ext)
    real(kind=8) :: m2p_mat((pm+1)*(pm+1), ngrid_ext_near)
    integer :: i, j, indi, indj
    do i = 0, pm
        indi = i*i + i + 1
        do j = -i, i
            indj = indi + j
            coef_sph_scaled(indj, :) = i * coef_sph(indj, :)
        end do
    end do
    call tree_m2m_get_mat(nclusters, children, cnode, rnode, pm, vscales, &
        & m2m_l2l_reflect_mat, m2m_l2l_ztrans_mat)
    call tree_m2m_use_mat(nsph, nclusters, pm, coef_sph_scaled, ind, cluster, &
        & children, cnode, rnode, m2m_l2l_reflect_mat, m2m_l2l_ztrans_mat, &
        & coef_node_m)
    call tree_m2l_get_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, pl, &
        & vscales, m2l_reflect_mat, m2l_ztrans_mat)
    call tree_m2l_use_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, pl, &
        & m2l_reflect_mat, m2l_ztrans_mat, coef_node_m, coef_node_l)
    call tree_l2l_get_mat(nclusters, children, cnode, rnode, pl, vscales, &
        & m2m_l2l_reflect_mat, m2m_l2l_ztrans_mat)
    call tree_l2l_use_mat(nsph, nclusters, pl, coef_node_l, ind, cluster, &
        & children, cnode, rnode, m2m_l2l_reflect_mat, m2m_l2l_ztrans_mat, &
        & coef_out)
    call tree_l2p_m2p_get_mat(nsph, csph, rsph, ngrid, grid, pm, pl, vscales, &
        & w, vgrid, nclusters, nnnear, snear, near, cluster, snode, ind, ui, &
        & ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, ngrid_ext_near, &
        & l2p_mat, m2p_mat)
    call tree_l2p_m2p_use_mat(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, nclusters, nnnear, snear, near, cluster, snode, &
        & ind, ui, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
        & ngrid_ext_near, l2p_mat, m2p_mat, coef_sph_scaled, coef_out, x)
    call int_ui_grid(nsph, pl, ngrid, w, vgrid, ui, x, coef_out)
end subroutine pcm_matvec_grid_fmm_test_mat

! Computes all M2M, M2L and L2L reflection and OZ translation matrices for FMM
subroutine pcm_matvec_grid_fmm_get_mat(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, snode, &
        & children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
        & m2m_reflect_mat, m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, &
        & l2l_reflect_mat, l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, &
        & grid_ext_ia, grid_ext_ja, l2p_mat, ngrid_ext_near, m2p_mat)
    integer, intent(in) :: nsph, ngrid, pm, pl, ind(nsph), nclusters
    integer, intent(in) :: cluster(2, nclusters), children(2, nclusters), nnfar
    integer, intent(in) :: sfar(nclusters+1), far(nnfar), nnnear, snode(nsph)
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: ngrid_ext, ngrid_ext_sph(nsph), grid_ext_ia(nsph+1)
    integer, intent(in) :: grid_ext_ja(ngrid_ext), ngrid_ext_near
    real(kind=8), intent(out) :: m2m_reflect_mat((pm+1)*(2*pm+1)*(2*pm+3)/3, &
        & nclusters-1)
    real(kind=8), intent(out) :: m2m_ztrans_mat((pm+1)*(pm+2)*(pm+3)/6, &
        & nclusters-1)
    real(kind=8), intent(out) :: l2l_reflect_mat((pl+1)*(2*pl+1)*(2*pl+3)/3, &
        & nclusters-1)
    real(kind=8), intent(out) :: l2l_ztrans_mat((pl+1)*(pl+2)*(pl+3)/6, &
        & nclusters-1)
    real(kind=8), intent(out) :: m2l_reflect_mat((max(pm,pl)+1) &
        & *(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3, nnfar)
    real(kind=8), intent(out) :: m2l_ztrans_mat((min(pm,pl)+1) &
        & *(min(pm,pl)+2)*(3*max(pm,pl)+3-min(pm,pl))/6, nnfar)
    real(kind=8), intent(out) :: l2p_mat((pl+1)*(pl+1), ngrid_ext)
    real(kind=8), intent(out) :: m2p_mat((pm+1)*(pm+1), ngrid_ext_near)
    call tree_m2m_get_mat(nclusters, children, cnode, rnode, pm, vscales, &
        & m2m_reflect_mat, m2m_ztrans_mat)
    call tree_m2l_get_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, pl, &
        & vscales, m2l_reflect_mat, m2l_ztrans_mat)
    call tree_l2l_get_mat(nclusters, children, cnode, rnode, pl, vscales, &
        & l2l_reflect_mat, l2l_ztrans_mat)
    call tree_l2p_m2p_get_mat(nsph, csph, rsph, ngrid, grid, pm, pl, vscales, &
        & w, vgrid, nclusters, nnnear, snear, near, cluster, snode, ind, ui, &
        & ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, ngrid_ext_near, &
        & l2p_mat, m2p_mat)
end subroutine pcm_matvec_grid_fmm_get_mat

! Matvec for ddPCM spherical harmonics by FMM by presaved matrices
subroutine pcm_matvec_grid_fmm_use_mat(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, snode, &
        & children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
        & m2m_reflect_mat, m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, &
        & l2l_reflect_mat, l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, &
        & grid_ext_ia, grid_ext_ja, l2p_mat, ngrid_ext_near, m2p_mat, &
        & coef_sph, coef_out)
    integer, intent(in) :: nsph, ngrid, pm, pl, ind(nsph), nclusters
    integer, intent(in) :: cluster(2, nclusters), children(2, nclusters), nnfar
    integer, intent(in) :: sfar(nclusters+1), far(nnfar), nnnear, snode(nsph)
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: ngrid_ext, ngrid_ext_sph(nsph), grid_ext_ia(nsph+1)
    integer, intent(in) :: grid_ext_ja(ngrid_ext), ngrid_ext_near
    real(kind=8), intent(in) :: coef_sph((pm+1)*(pm+1), nsph)
    real(kind=8), intent(out) :: coef_out((pl+1)*(pl+1), nsph)
    real(kind=8) :: x(ngrid, nsph)
    real(kind=8) :: coef_sph_scaled((pm+1)*(pm+1), nsph)
    real(kind=8) :: coef_node_m((pm+1)*(pm+1), nclusters)
    real(kind=8) :: coef_node_l((pl+1)*(pl+1), nclusters)
    real(kind=8), intent(in) :: m2m_reflect_mat((pm+1)*(2*pm+1)*(2*pm+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: m2m_ztrans_mat((pm+1)*(pm+2)*(pm+3)/6, &
        & nclusters-1)
    real(kind=8), intent(in) :: l2l_reflect_mat((pl+1)*(2*pl+1)*(2*pl+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: l2l_ztrans_mat((pl+1)*(pl+2)*(pl+3)/6, &
        & nclusters-1)
    real(kind=8), intent(in) :: m2l_reflect_mat((max(pm,pl)+1) &
        & *(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3, nnfar)
    real(kind=8), intent(in) :: m2l_ztrans_mat((min(pm,pl)+1) &
        & *(min(pm,pl)+2)*(3*max(pm,pl)+3-min(pm,pl))/6, nnfar)
    real(kind=8), intent(in) :: l2p_mat((pl+1)*(pl+1), ngrid_ext)
    real(kind=8), intent(in) :: m2p_mat((pm+1)*(pm+1), ngrid_ext_near)
    integer :: i, j, indi, indj
    integer :: counter=1
    real(kind=8) :: start, finish, time
    call cpu_time(start)
    do i = 0, pm
        indi = i*i + i + 1
        do j = -i, i
            indj = indi + j
            coef_sph_scaled(indj, :) = i * coef_sph(indj, :)
        end do
    end do
    call tree_m2m_use_mat(nsph, nclusters, pm, coef_sph_scaled, ind, cluster, &
        & children, cnode, rnode, m2m_reflect_mat, m2m_ztrans_mat, &
        & coef_node_m)
    call tree_m2l_use_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, pl, &
        & m2l_reflect_mat, m2l_ztrans_mat, coef_node_m, coef_node_l)
    call tree_l2l_use_mat(nsph, nclusters, pl, coef_node_l, ind, cluster, &
        & children, cnode, rnode, l2l_reflect_mat, l2l_ztrans_mat, &
        & coef_out)
    call tree_l2p_m2p_use_mat(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, nclusters, nnnear, snear, near, cluster, snode, &
        & ind, ui, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
        & ngrid_ext_near, l2p_mat, m2p_mat, coef_sph_scaled, coef_out, x)
    call int_ui_grid(nsph, pl, ngrid, w, vgrid, ui, x, coef_out)
    call cpu_time(finish)
    total_time_matvec = total_time_matvec + finish - start
    total_count_matvec = total_count_matvec + 1
end subroutine pcm_matvec_grid_fmm_use_mat

! Adjoint matvec for ddPCM spherical harmonics by FMM by presaved matrices
subroutine pcm_matvec_adj_grid_fmm_use_mat(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, snode, &
        & children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
        & m2m_reflect_mat, m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, &
        & l2l_reflect_mat, l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, &
        & grid_ext_ia, grid_ext_ja, l2p_mat, ngrid_ext_near, m2p_mat, &
        & coef_sph, coef_out)
    integer, intent(in) :: nsph, ngrid, pm, pl, ind(nsph), nclusters
    integer, intent(in) :: cluster(2, nclusters), children(2, nclusters), nnfar
    integer, intent(in) :: sfar(nclusters+1), far(nnfar), nnnear, snode(nsph)
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: ngrid_ext, ngrid_ext_sph(nsph), grid_ext_ia(nsph+1)
    integer, intent(in) :: grid_ext_ja(ngrid_ext), ngrid_ext_near
    real(kind=8), intent(in) :: coef_sph((pl+1)*(pl+1), nsph)
    real(kind=8), intent(out) :: coef_out((pm+1)*(pm+1), nsph)
    real(kind=8) :: xgrid(ngrid, nsph)
    real(kind=8) :: coef_sph_l((pl+1)*(pl+1), nsph)
    real(kind=8) :: coef_node_m((pm+1)*(pm+1), nclusters)
    real(kind=8) :: coef_node_l((pl+1)*(pl+1), nclusters)
    real(kind=8), intent(in) :: m2m_reflect_mat((pm+1)*(2*pm+1)*(2*pm+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: m2m_ztrans_mat((pm+1)*(pm+2)*(pm+3)/6, &
        & nclusters-1)
    real(kind=8), intent(in) :: l2l_reflect_mat((pl+1)*(2*pl+1)*(2*pl+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: l2l_ztrans_mat((pl+1)*(pl+2)*(pl+3)/6, &
        & nclusters-1)
    real(kind=8), intent(in) :: m2l_reflect_mat((max(pm,pl)+1) &
        & *(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3, nnfar)
    real(kind=8), intent(in) :: m2l_ztrans_mat((min(pm,pl)+1) &
        & *(min(pm,pl)+2)*(3*max(pm,pl)+3-min(pm,pl))/6, nnfar)
    real(kind=8), intent(in) :: l2p_mat((pl+1)*(pl+1), ngrid_ext)
    real(kind=8), intent(in) :: m2p_mat((pm+1)*(pm+1), ngrid_ext_near)
    integer :: i, j, indi, indj
    integer :: counter=1
    real(kind=8) :: start, finish, time
    call cpu_time(start)
    do i = 1, nsph
        call int_grid_adj(pl, ngrid, w, vgrid, coef_sph(:, i), xgrid(:, i))
    end do
    call tree_l2p_m2p_adj_use_mat(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, nclusters, nnnear, snear, near, cluster, snode, &
        & ind, ui, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
        & ngrid_ext_near, l2p_mat, m2p_mat, xgrid, coef_out, coef_sph_l)
    call tree_l2l_adj_use_mat(nsph, nclusters, pl, coef_node_l, ind, cluster, &
        & children, cnode, rnode, l2l_reflect_mat, l2l_ztrans_mat, &
        & coef_sph_l)
    call tree_m2l_adj_use_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, pl, &
        & m2l_reflect_mat, m2l_ztrans_mat, coef_node_l, coef_node_m)
    call tree_m2m_adj_use_mat(nsph, nclusters, pm, coef_out, ind, cluster, &
        & children, cnode, rnode, m2m_reflect_mat, m2m_ztrans_mat, &
        & coef_node_m)
    do i = 0, pm
        indi = i*i + i + 1
        do j = -i, i
            indj = indi + j
            coef_out(indj, :) = i * coef_out(indj, :)
        end do
    end do
    call cpu_time(finish)
    total_time_matvec = total_time_matvec + finish - start
    total_count_matvec = total_count_matvec + 1
end subroutine pcm_matvec_adj_grid_fmm_use_mat

! Compute PCM portion of forces (2 matvecs)
subroutine pcm_force_grid_fmm_use_mat(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, lmax, vscales, ind, nclusters, cluster, snode, &
        & children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
        & m2m_reflect_mat, m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, &
        & l2l_reflect_mat, l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, &
        & grid_ext_ia, grid_ext_ja, l2p_mat, ngrid_ext_near, m2p_mat, &
        & coef_sph, xgrid, force_out)
    integer, intent(in) :: nsph, ngrid, pm, pl, lmax, ind(nsph), nclusters
    integer, intent(in) :: cluster(2, nclusters), children(2, nclusters), nnfar
    integer, intent(in) :: sfar(nclusters+1), far(nnfar), nnnear, snode(nsph)
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: ngrid_ext, ngrid_ext_sph(nsph), grid_ext_ia(nsph+1)
    integer, intent(in) :: grid_ext_ja(ngrid_ext), ngrid_ext_near
    real(kind=8), intent(in) :: coef_sph((lmax+1)*(lmax+1), nsph)
    real(kind=8), intent(in) :: xgrid(ngrid, nsph)
    real(kind=8), intent(out) :: force_out(3, nsph)
    real(kind=8) :: coef_sph_scaled((pm+1)*(pm+1), nsph)
    real(kind=8) :: coef_sph_m_adj((pm+1)*(pm+1), nsph)
    real(kind=8) :: coef_sph_l((pl+1)*(pl+1), nsph)
    real(kind=8) :: coef_sph_l_adj((pl+1)*(pl+1), nsph)
    real(kind=8) :: coef_node_m((pm+1)*(pm+1), nclusters)
    real(kind=8) :: coef_node_l((pl+1)*(pl+1), nclusters)
    real(kind=8), intent(in) :: m2m_reflect_mat((pm+1)*(2*pm+1)*(2*pm+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: m2m_ztrans_mat((pm+1)*(pm+2)*(pm+3)/6, &
        & nclusters-1)
    real(kind=8), intent(in) :: l2l_reflect_mat((pl+1)*(2*pl+1)*(2*pl+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: l2l_ztrans_mat((pl+1)*(pl+2)*(pl+3)/6, &
        & nclusters-1)
    real(kind=8), intent(in) :: m2l_reflect_mat((max(pm,pl)+1) &
        & *(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3, nnfar)
    real(kind=8), intent(in) :: m2l_ztrans_mat((min(pm,pl)+1) &
        & *(min(pm,pl)+2)*(3*max(pm,pl)+3-min(pm,pl))/6, nnfar)
    real(kind=8), intent(in) :: l2p_mat((pl+1)*(pl+1), ngrid_ext)
    real(kind=8), intent(in) :: m2p_mat((pm+1)*(pm+1), ngrid_ext_near)
    integer :: i, j, indi, indj, l, m, isph, igrid, ik, ksph, jsph, jsph_node
    integer :: counter=1
    real(kind=8) :: start, finish, time
    integer :: inear, inode, jnode
    real(kind=8) :: fac, fl, gg, cx, cy, cz, c(3), vki(3), vvki, tki
    real(kind=8) :: tlow, thigh, xgrid2(ngrid, nsph), xgrid3(ngrid, nsph)
    real(kind=8) :: coef_sph_l_x((pl+1)*(pl+1)), coef_sph_l_y((pl+1)*(pl+1))
    real(kind=8) :: coef_sph_l_x2((pl+1)*(pl+1)), coef_sph_l_y2((pl+1)*(pl+1))
    real(kind=8) :: coef_sph_l_z((pl+1)*(pl+1))
    real(kind=8) :: coef_sph_m_x((lmax+2)*(lmax+2), nsph)
    real(kind=8) :: coef_sph_m_x2((lmax+2)*(lmax+2), nsph)
    real(kind=8) :: coef_sph_m_y((lmax+2)*(lmax+2), nsph)
    real(kind=8) :: coef_sph_m_y2((lmax+2)*(lmax+2), nsph)
    real(kind=8) :: coef_sph_m_z((lmax+2)*(lmax+2), nsph)
    real(kind=8) :: zx_coord_transform(3, 3), zy_coord_transform(3, 3)
    real(kind=8) :: fact(2*pl+1), tmp1, tmp2, tmp_gg, gg3(3)
    real(kind=8) :: zxz_mat((pl+1)*(2*pl+1)*(2*pl+3)/3)
    real(kind=8) :: zyz_mat((pl+1)*(2*pl+1)*(2*pl+3)/3)
    call cpu_time(start)
    zx_coord_transform = 0
    zx_coord_transform(3, 2) = 1
    zx_coord_transform(2, 3) = 1
    zx_coord_transform(1, 1) = 1
    zy_coord_transform = 0
    zy_coord_transform(1, 2) = 1
    zy_coord_transform(2, 1) = 1
    zy_coord_transform(3, 3) = 1
    force_out = 0
    tlow  = one - pt5*(one - se)*eta
    thigh = one + pt5*(one + se)*eta
    ! Fill square roots of factorials
    fact(1) = 1
    do j = 2, 2*pl+1
        fact(j) = sqrt(dble(j-1)) * fact(j-1)
    end do
    ! Adjoint full FMM matvec to get output multipole expansions from input
    ! external grid points. It will be used in R_i^B.
    do isph = 1, nsph
        xgrid3(:, isph) = xgrid(:, isph) * w(:)
    end do
    call tree_l2p_m2p_adj_use_mat(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, nclusters, nnnear, snear, near, cluster, snode, &
        & ind, ui, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
        & ngrid_ext_near, l2p_mat, m2p_mat, xgrid3, coef_sph_m_adj, &
        & coef_sph_l_adj)
    call tree_l2l_adj_use_mat(nsph, nclusters, pl, coef_node_l, ind, cluster, &
        & children, cnode, rnode, l2l_reflect_mat, l2l_ztrans_mat, &
        & coef_sph_l_adj)
    call tree_m2l_adj_use_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, pl, &
        & m2l_reflect_mat, m2l_ztrans_mat, coef_node_l, coef_node_m)
    call tree_m2m_adj_use_mat(nsph, nclusters, pm, coef_sph_m_adj, ind, cluster, &
        & children, cnode, rnode, m2m_reflect_mat, m2m_ztrans_mat, &
        & coef_node_m)
    ! Direct far-field FMM matvec to get output local expansions from input
    ! multipole expansions. It will be used in R_i^A.
    ! As of now I compute potential at all external grid points, improved
    ! version shall only compute it at external points in a switch region
    coef_sph_scaled = 0
    do i = 0, lmax
        indi = i*i + i + 1
        do j = -i, i
            indj = indi + j
            coef_sph_scaled(indj, :) = i * coef_sph(indj, :)
        end do
    end do
    call tree_m2m_use_mat(nsph, nclusters, pm, coef_sph_scaled, ind, cluster, &
        & children, cnode, rnode, m2m_reflect_mat, m2m_ztrans_mat, &
        & coef_node_m)
    call tree_m2l_use_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, pl, &
        & m2l_reflect_mat, m2l_ztrans_mat, coef_node_m, coef_node_l)
    call tree_l2l_use_mat(nsph, nclusters, pl, coef_node_l, ind, cluster, &
        & children, cnode, rnode, l2l_reflect_mat, l2l_ztrans_mat, &
        & coef_sph_l)
    xgrid2 = 0
    call tree_l2p_m2p_use_mat(nsph, csph, rsph, ngrid, grid, pm, pl, &
        & vscales, w, vgrid, nclusters, nnnear, snear, near, cluster, snode, &
        & ind, ui, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
        & ngrid_ext_near, l2p_mat, m2p_mat, coef_sph_scaled, coef_sph_l, xgrid2)
    !! Precompute M2M gradients for each sphere
    !! Derivative along OZ axis is a straightforward derivative of M2M along OZ
    !! axis with zero shift. Derivatives along OX and OY axes are computed with
    !! help of change of coordinates (OX/OY->OZ, OZ derivative, OZ->OX/OY).
    ! Get matrices of OX->OZ->OX and OY->OZ->OY
    call fmm_sph_get_reflect_mat(pl, zx_coord_transform, zxz_mat)
    call fmm_sph_get_reflect_mat(pl, zy_coord_transform, zyz_mat)
    ! At first init multipole coefficients for every axis
    do isph = 1, nsph
        coef_sph_m_z(:, isph) = coef_sph_scaled(:(lmax+1)*(lmax+1), isph)
        call fmm_sph_use_reflect_mat(lmax, zxz_mat, coef_sph_m_z(:, isph), &
            & coef_sph_m_x(:, isph))
        call fmm_sph_use_reflect_mat(lmax, zyz_mat, coef_sph_m_z(:, isph), &
            & coef_sph_m_y(:, isph))
    end do
    ! Fill zeros at lmax+1 degree since it is zero from the model
    coef_sph_m_x((lmax+1)*(lmax+1)+1:, :) = 0
    coef_sph_m_y((lmax+1)*(lmax+1)+1:, :) = 0
    coef_sph_m_z((lmax+1)*(lmax+1)+1:, :) = 0
    ! Get derivative of M2M from the origin to the origin along OZ axis
    do l = lmax+1, 1, -1
        indi = l*l + l + 1
        indj = (l-1)*(l-1) + (l-1) + 1
        do m = 1-l, l-1
            tmp1 = sqrt(dble(l*l-m*m)) * sqrt(dble(2*l+1)) / sqrt(dble(2*l-1))
            coef_sph_m_x(indi+m, :) = tmp1 * coef_sph_m_x(indj+m, :)
            coef_sph_m_y(indi+m, :) = tmp1 * coef_sph_m_y(indj+m, :)
            coef_sph_m_z(indi+m, :) = tmp1 * coef_sph_m_z(indj+m, :)
        end do
        coef_sph_m_x(indi-l, :) = 0
        coef_sph_m_x(indi+l, :) = 0
        coef_sph_m_y(indi-l, :) = 0
        coef_sph_m_y(indi+l, :) = 0
        coef_sph_m_z(indi-l, :) = 0
        coef_sph_m_z(indi+l, :) = 0
    end do
    coef_sph_m_x(1, :) = 0
    coef_sph_m_y(1, :) = 0
    coef_sph_m_z(1, :) = 0
    ! Derivatives involve 1/rsph factor
    do isph = 1, nsph
        coef_sph_m_x(:, isph) = coef_sph_m_x(:, isph) / rsph(isph)
        coef_sph_m_y(:, isph) = coef_sph_m_y(:, isph) / rsph(isph)
        coef_sph_m_z(:, isph) = coef_sph_m_z(:, isph) / rsph(isph)
    end do
    ! Now get coordinates OX and OY back
    do isph = 1, nsph
        call fmm_sph_use_reflect_mat(lmax+1, zxz_mat, coef_sph_m_x(:, isph), &
            & coef_sph_m_x2(:, isph))
        call fmm_sph_use_reflect_mat(lmax+1, zyz_mat, coef_sph_m_y(:, isph), &
            & coef_sph_m_y2(:, isph))
    end do
    ! Compute all terms of grad_i(R)
    do isph = 1, nsph
        do igrid = 1, ngrid
            ! Loop over all neighbouring spheres
            do ik = inl(isph), inl(isph+1) - 1
                ksph = nl(ik)
                ! Only consider outside grid points
                if(ui(igrid, ksph) .le. zero) cycle
                ! build geometrical quantities
                cx = csph(1,ksph) + rsph(ksph)*grid(1,igrid)
                cy = csph(2,ksph) + rsph(ksph)*grid(2,igrid)
                cz = csph(3,ksph) + rsph(ksph)*grid(3,igrid)
                vki(1) = cx - csph(1,isph)
                vki(2) = cy - csph(2,isph)
                vki(3) = cz - csph(3,isph)
                vvki = sqrt(vki(1)*vki(1) + vki(2)*vki(2) + &
                    & vki(3)*vki(3))
                tki = vvki/rsph(isph)
                ! Only consider such points where grad U is non-zero
                if((tki.le.tlow) .or. (tki.ge.thigh)) cycle
                gg = zero
                do l = 0, lmax
                    indi = l*l + l + 1
                    fl = dble(l)
                    fac = two*pi/(two*fl + one)
                    do m = -l, l 
                        ! This is R^D component (grad_i of U of a sum of R_kk
                        ! with index inequality k!=i for k being a neighbour
                        ! of i)
                        ! Indexes k and j are flipped compared to the paper
                        gg = gg + fac*basis(indi+m,igrid)*coef_sph(indi+m,ksph)
                    end do
                end do
                ! This is R^C and R^B component (grad_i of U of a sum of R_kj
                ! for index inequality j!=k)
                ! Indexes k and j are flipped compared to the paper
                gg = gg - xgrid2(igrid, ksph)
                ! Compute grad_i component of forces using precomputed
                ! potential gg
                force_out(:, isph) = force_out(:, isph) - &
                    & dfsw(tki,se,eta)/rsph(isph)*w(igrid)*gg* &
                    & xgrid(igrid, ksph)*(vki/vvki)
            end do
            ! contribution from the sphere itself
            if((ui(igrid,isph).gt.zero) .and. (ui(igrid,isph).lt.one)) then
                gg = zero
                do l = 0, lmax
                    indi = l*l + l + 1
                    fl = dble(l)
                    fac = two*pi/(two*fl + one)
                    do m = -l, l 
                        ! This is R^D component (grad_i of U of R_ii)
                        gg = gg + fac*basis(indi+m,igrid)*coef_sph(indi+m,isph)
                    end do
                end do
                ! R^A component (grad_i of U of a sum of R_ij for index
                ! inequality j!=i)
                ! Indexes k and j are flipped compared to the paper
                gg = gg - xgrid2(igrid, isph)
                ! Compute grad_i component of forces using precomputed
                ! potential gg
                force_out(:, isph) = force_out(:, isph) + &
                    & w(igrid)*gg*xgrid(igrid,isph)*zi(:, igrid, isph)
            end if
            if (ui(igrid, isph) .gt. zero) then
                ! Another R^A component (grad_i of potential of a sum of R_ij
                ! for index inequality j!=i)
                ! Indexes k and j are flipped compared to the paper
                gg3 = 0
                ! Get rotated local expansions
                call fmm_sph_use_reflect_mat(pl, zxz_mat, &
                    & coef_sph_l(:, isph), coef_sph_l_x)
                call fmm_sph_use_reflect_mat(pl, zyz_mat, &
                    & coef_sph_l(:, isph), coef_sph_l_y)
                coef_sph_l_z = coef_sph_l(:, isph)
                ! Gradient of the far-field potential is a gradient of local
                ! expansion
                do l = 0, pl-1
                    indi = l*l + l + 1
                    indj = (l+1)*(l+1) + (l+1) + 1
                    do m = -l, l
                        tmp1 = one / rsph(isph) / fact(l-m+1) / fact(l+m+1) / &
                            & vscales(indi) / vscales(indj) * fact(l-m+2) * &
                            & fact(l+m+2)
                        coef_sph_l_x(indi+m) = tmp1 * coef_sph_l_x(indj+m)
                        coef_sph_l_y(indi+m) = tmp1 * coef_sph_l_y(indj+m)
                        coef_sph_l_z(indi+m) = tmp1 * coef_sph_l_z(indj+m)
                    end do
                end do
                call fmm_sph_use_reflect_mat(pl, zxz_mat, &
                    & coef_sph_l_x, coef_sph_l_x2)
                call fmm_sph_use_reflect_mat(pl, zyz_mat, &
                    & coef_sph_l_y, coef_sph_l_y2)
                do l = 0, pl-1
                    indi = l*l + l + 1
                    do m = indi-l, indi+l
                        gg3(1) = gg3(1) - &
                            & coef_sph_l_x2(m)*vgrid(m, igrid)
                        gg3(2) = gg3(2) - &
                            & coef_sph_l_y2(m)*vgrid(m, igrid)
                        gg3(3) = gg3(3) - &
                            & coef_sph_l_z(m)*vgrid(m, igrid)
                    end do
                end do
                ! Gradient of the near-field potential is a gradient of
                ! multipole expansion
                inode = snode(isph)
                do inear = snear(inode), snear(inode+1)-1
                    jnode = near(inear)
                    do jsph_node = cluster(1, jnode), cluster(2, jnode)
                        jsph = ind(jsph_node)
                        if (isph .eq. jsph) cycle
                        c = csph(:, isph) + rsph(isph)*grid(:, igrid)
                        call fmm_m2p(c-csph(:, jsph), rsph(jsph), lmax+1, &
                            & vscales, coef_sph_m_x2(:, jsph), tmp_gg)
                        gg3(1) = gg3(1) + tmp_gg
                        call fmm_m2p(c-csph(:, jsph), rsph(jsph), lmax+1, &
                            & vscales, coef_sph_m_y2(:, jsph), tmp_gg)
                        gg3(2) = gg3(2) + tmp_gg
                        call fmm_m2p(c-csph(:, jsph), rsph(jsph), lmax+1, &
                            & vscales, coef_sph_m_z(:, jsph), tmp_gg)
                        gg3(3) = gg3(3) + tmp_gg
                    end do
                end do
                ! Accumulate all computed forces
                force_out(:, isph) = force_out(:, isph) - &
                    & w(igrid)*gg3*xgrid(igrid, isph)*ui(igrid, isph)
            end if
        end do
        ! Another R^B component (grad_i of potential of a sum of R_ji
        ! for index inequality j!=i)
        ! Indexes k and j are flipped compared to the paper
        gg3 = 0
        gg3(1) = dot_product(coef_sph_m_adj(:(lmax+2)*(lmax+2), isph), &
            & coef_sph_m_x2(:, isph))
        gg3(2) = dot_product(coef_sph_m_adj(:(lmax+2)*(lmax+2), isph), &
            & coef_sph_m_y2(:, isph))
        gg3(3) = dot_product(coef_sph_m_adj(:(lmax+2)*(lmax+2), isph), &
            & coef_sph_m_z(:, isph))
        force_out(:, isph) = force_out(:, isph) + gg3
    end do
    call cpu_time(finish)
    total_time_matvec = total_time_matvec + finish - start
    total_count_matvec = total_count_matvec + 1
end subroutine pcm_force_grid_fmm_use_mat

! Generate matrix of PCM equations (with FMM code)
subroutine tree_matrix_fmm(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, snode, &
        & children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
        & m2m_reflect_mat, m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, &
        & l2l_reflect_mat, l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, &
        & grid_ext_ia, grid_ext_ja, l2p_mat, ngrid_ext_near, m2p_mat, &
        & mat_out)
    integer, intent(in) :: nsph, ngrid, pm, pl, ind(nsph), nclusters
    integer, intent(in) :: cluster(2, nclusters), children(2, nclusters), nnfar
    integer, intent(in) :: sfar(nclusters+1), far(nnfar), nnnear, snode(nsph)
    integer, intent(in) :: snear(nclusters+1), near(nnnear)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: cnode(3, nclusters), rnode(nclusters)
    integer, intent(in) :: ngrid_ext, ngrid_ext_sph(nsph), grid_ext_ia(nsph+1)
    integer, intent(in) :: grid_ext_ja(ngrid_ext), ngrid_ext_near
    real(kind=8) :: coef_sph_m((pm+1)*(pm+1), nsph)
    real(kind=8), intent(out) :: mat_out((pl+1)*(pl+1), nsph, (pm+1)*(pm+1), nsph)
    real(kind=8) :: coef_sph_l((pl+1)*(pl+1), nsph)
    real(kind=8) :: coef_node_m((pm+1)*(pm+1), nclusters)
    real(kind=8) :: coef_node_l((pl+1)*(pl+1), nclusters)
    real(kind=8), intent(in) :: m2m_reflect_mat((pm+1)*(2*pm+1)*(2*pm+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: m2m_ztrans_mat((pm+1)*(pm+2)*(pm+3)/6, &
        & nclusters-1)
    real(kind=8), intent(in) :: l2l_reflect_mat((pl+1)*(2*pl+1)*(2*pl+3)/3, &
        & nclusters-1)
    real(kind=8), intent(in) :: l2l_ztrans_mat((pl+1)*(pl+2)*(pl+3)/6, &
        & nclusters-1)
    real(kind=8), intent(in) :: m2l_reflect_mat((max(pm,pl)+1) &
        & *(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3, nnfar)
    real(kind=8), intent(in) :: m2l_ztrans_mat((min(pm,pl)+1) &
        & *(min(pm,pl)+2)*(3*max(pm,pl)+3-min(pm,pl))/6, nnfar)
    real(kind=8), intent(in) :: l2p_mat((pl+1)*(pl+1), ngrid_ext)
    real(kind=8), intent(in) :: m2p_mat((pm+1)*(pm+1), ngrid_ext_near)
    integer :: i, j, indi, indj
    integer :: counter=1
    real(kind=8) :: start, finish, time
    do i = 1, (pm+1)*(pm+1)
        do j = 1, nsph
            coef_sph_m(:, :) = 0
            coef_sph_m(i, j) = 1
            coef_sph_l(:, :) = 0
            call pcm_matvec_grid_fmm_use_mat(nsph, csph, rsph, ngrid, grid, w, &
                    & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, snode, &
                    & children, &
                    & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
                    & m2m_reflect_mat, m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, &
                    & l2l_reflect_mat, l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, &
                    & grid_ext_ia, grid_ext_ja, l2p_mat, ngrid_ext_near, m2p_mat, &
                    & coef_sph_m, coef_sph_l)
             mat_out(:, :, i, j) = coef_sph_l
        end do
    end do
end subroutine tree_matrix_fmm

subroutine pcm_matvec_print_stats
    write(*, "(A,I0,A,ES9.3E2,A)") "FMM matvecs: ", total_count_matvec, &
        & " matvecs in ", total_time_matvec, " seconds"
end subroutine pcm_matvec_print_stats

! Apply matvec for ddPCM spherical harmonics by tree-code
! Per-sphere tree-code, which used not only M2P, but M2L and L2P also
subroutine pcm_matvec_grid_treecode2(nsph, csph, rsph, ngrid, grid, w, vgrid, &
        & ui, pm, pl, vscales, ind, cluster, children, cnode, rnode, &
        & coef_sph, coef_out)
! Parameters:
!   nsph: number of all spheres
!   csph: centers of all spheres
!   rsph: radiuses of all spheres
!   ngrid: number of Lebedev grid points on each sphere
!   grid: coordinates of Lebedev grid points on a unit sphere
!   ui: "outside" factor of grid points (0 is inside, 1 is "fully" outside)
!   pm: maximum degree of multipole basis functions
!   pl: maximum degree of local basis functions
!   vscales: normalization constants for Y_lm
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   coef_sph: multipole coefficients of bounding spheres of nodes
!   coef_out: output multipole coefficients of spherical harmonics
    integer, intent(in) :: nsph, ngrid, pm, pl, ind(nsph), cluster(2, 2*nsph-1)
    integer, intent(in) :: children(2, 2*nsph-1)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    real(kind=8), intent(in) :: coef_sph((pm+1)*(pm+1), nsph)
    real(kind=8), intent(out) :: coef_out((pl+1)*(pl+1), nsph)
    real(kind=8) :: coef_sph_scaled((pm+1)*(pm+1), nsph), coef_l((pl+1)*(pl+1))
    real(kind=8) :: coef_node((pm+1)*(pm+1), 2*nsph-1), d(3), x(ngrid)
    real(kind=8) :: y(ngrid), tmp_v, r
    integer :: i, j, k(2), indi, indj, far(2*nsph-1), igrid
    do i = 0, pm
        indi = i*i + i + 1
        do j = -i, i
            indj = indi + j
            coef_sph_scaled(indj, :) = i * coef_sph(indj, :)
        end do
    end do
    call tree_m2m_fast(nsph, 2*nsph-1, pm, vscales, coef_sph_scaled, ind, &
        & cluster, children, cnode, rnode, coef_node)
    do i = 1, nsph
        coef_l = 0
        x = 0
        far = 0
        ! far(i): 0 if M2P or M2L must be ignored, 1 if M2L must be applied,
        ! 2 if M2P must be applied and 3 if need to decide if M2P or M2L is
        ! applicable
        far(1) = 3
        do j = 1, 2*nsph-1
            ! If M2P or M2L is not applicable, ignore node
            if (far(j) .eq. 0) then
                cycle
            end if
            k = children(:, j)
            d = cnode(:, j) - csph(:, i)
            r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
            ! constant 2 guarantees good convergence and accuracy
            if (r .gt. 3*max(rsph(i), rnode(j))) then
                far(j) = 1
            else if(k(1) .eq. 0) then
                far(j) = 2
                ! If it contains the same sphere, then ignore it
                if (ind(cluster(1, j)) .eq. i) then
                    far(j) = 0
                end if
            else
                far(j) = 0
                far(k(1)) = 3
                far(k(2)) = 3
            end if
            ! If M2L is needed
            if (far(j) .eq. 1) then
                call fmm_m2l_fast(d, rnode(j), rsph(i), pm, pl, vscales, &
                    & coef_node(:, j), coef_l)
            ! If M2P is needed
            else if (far(j) .eq. 2) then
                do igrid = 1, ngrid
                    if (ui(igrid, i) .ne. 0) then
                        call fmm_m2p(rsph(i)*grid(:, igrid)-d, rnode(j), pm, &
                            vscales, coef_node(:, j), tmp_v)
                        x(igrid) = x(igrid) + tmp_v
                    end if
                end do
            end if
        end do
        ! Evaluate L2P and add it to M2P of near-field
        ! Apply normalization factor to coefficients
        do j = 0, pl
            coef_l(j*j+1:(j+1)*(j+1)) = coef_l(j*j+1:(j+1)*(j+1)) / &
                &(vscales(j*j+j+1)**2)
        end do
        ! Apply L2P
        y = 0
        do j = 1, (pl+1)*(pl+1)
            y = y + coef_l(j)*vgrid(j, :)
        end do
        ! Scale by ui
        x = (x+y) * ui(:, i)
        call int_grid(pl, ngrid, w, vgrid, x, coef_out(:, i))
    end do
end subroutine pcm_matvec_grid_treecode2

end module pcm_fmm

