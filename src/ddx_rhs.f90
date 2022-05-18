!> Routines to build rhs (phi and psi)
module ddx_rhs

use ddx_definitions
use ddx_core

implicit none

contains

!> Given a multipolar distribution, compute the potential at the target points
!> this is done with or without FMMs depending on the relevant flag
!> The multipoles must be centered on the ddx spheres.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points size (ncav)
!!
subroutine build_phi(params, constants, workspace, multipoles, &
        & mmax, phi_cav)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax
    real(dp), intent(inout) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav)
    if (params % fmm .eq. 0) then
        call build_phi_dense(params, constants, workspace, multipoles, &
            & params % csph, mmax, params % nsph, phi_cav, constants % ccav, &
            & constants % ncav)
    else if (params % fmm .eq. 1) then
        call build_phi_fmm(params, constants, workspace, multipoles, mmax, &
            & phi_cav)
    end if
end subroutine build_phi

!> Given a multipolar distribution, compute the potential at the target points
!> using a N^2 code. 
!> As this routine does not use the FMM machinery it is more flexible and
!> accepts arbitrary sources and targets.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2,nm)
!! @param[in] cm: centers of the multipolar distributions size (3,nm)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[in] nm: number of multipoles
!! @param[out] phi_cav: electric potential at the target points size (ncav)
!! @param[in] ccav: coordinates of the target points size (3,ncav)
!! @param[in] ncav: number of target points
!!
subroutine build_phi_dense(params, constants, workspace, multipoles, cm, &
        & mmax, nm, phi_cav, ccav, ncav)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax, nm, ncav
    real(dp), intent(inout) :: multipoles((mmax + 1)**2, nm)
    real(dp), intent(in) :: cm(3, nm)
    real(dp), intent(out) :: phi_cav(ncav)
    real(dp), intent(in) :: ccav(3, ncav)
    integer icav, im, l, m, i
    real(dp) :: v, c(3)
    real(dp) :: r

    ! m2p expects a strange normalization, hence the multipoles must be
    ! scaled by l, only for l>=1

    do im = 1, nm
        do l = 1, mmax
            i = l*l + l + 1
            do m = -l, l
                 multipoles(i+m,im) = multipoles(i+m,im)*dble(l)
            end do
        end do
    end do

    do icav = 1, ncav
        v = zero
        do im = 1, nm
            c(:) = ccav(:, icav) - cm(:, im)
            call fmm_m2p(c, one, mmax, constants % vscales_rel, &
                & one, multipoles(:, im), one, v)
        end do
        phi_cav(icav) = v
    end do

    ! restore the original multipoles

    do im = 1, nm
        do l = 1, mmax
            i = l*l + l + 1
            do m = -l, l
                 multipoles(i+m,im) = multipoles(i+m,im)/dble(l)
            end do
        end do
    end do

end subroutine build_phi_dense

!> Given a multipolar distribution, compute the potential at the target points
!> using FMMs.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2,nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points size (ncav)
!!
subroutine build_phi_fmm(params, constants, workspace, multipoles, mmax, &
        & phi_cav)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav)
end subroutine build_phi_fmm

!> Given a multipolar distribution, assemble the RHS psi.
!> The multipoles must be centered on the ddx spheres.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2,nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] psi: RHS for adjoint linear systems
!!     size ((lmax+1)**2,nsph), the internal lmax should be >= mmax
!!
subroutine build_psi(params, constants, workspace, multipoles, mmax, psi)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax+1)**2, params % nsph)
    real(dp), intent(out) :: psi((params % lmax+1)**2, params % nsph)
    integer :: isph, l, m, i
    real(dp) :: v

    psi = zero
    do isph = 1, params % nsph
        do l = 0, mmax
            v = fourpi/((two*dble(l) + one)*(params % rsph(isph)**(l)))
            i = l*l + l + 1
            do m = -l, l
                psi(i + m, isph) = v*multipoles(i + m, isph)
            end do
        end do
    end do
end subroutine build_psi

end module ddx_rhs
