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

    do icav = 1, ncav
        v = zero
        do im = 1, nm
            c(:) = ccav(:, icav) - cm(:, im)
            call fmm_m2p(c, one, mmax, constants % vscales_rel, &
                & one, multipoles(:, im), one, v)
        end do
        phi_cav(icav) = v
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
    ! local variables
    integer i, isph, igrid, icav, inode, tmp_max, l, m, ind

    ! P2M is not needed as we have already multipolar distributions
    ! we just have to copy the multipoles in the proper places
    do isph = 1, params % nsph
        inode = constants % snode(isph)
        workspace % tmp_sph(:, isph) = zero
        workspace % tmp_node_m(:, inode) = zero
        tmp_max = min(params % pm, mmax)
        do l = 0, tmp_max
            ind = l*l + l + 1
            do m = -l, l
                workspace % tmp_sph(ind+m, isph) = &
                    & multipoles(ind+m, isph)/(params%rsph(isph)**(l+1))
                workspace % tmp_node_m(ind+m, inode) = &
                    & multipoles(ind+m, isph)/(params%rsph(isph)**(l+1))
            end do
        end do
    end do

    call tree_m2m_rotation(params, constants, workspace % tmp_node_m)
    call tree_m2l_rotation(params, constants, workspace % tmp_node_m, &
        & workspace % tmp_node_l)
    call tree_l2l_rotation(params, constants, workspace % tmp_node_l)
    call tree_l2p(params, constants, one, workspace % tmp_node_l, zero, &
        & workspace % tmp_grid, workspace % tmp_sph_l)

    ! near field
    call tree_m2p(params, constants, params % lmax, one, &
        & workspace % tmp_sph, one, workspace % tmp_grid)

    ! Potential from each sphere to its own grid points
    call dgemm('T', 'N', params % ngrid, params % nsph, &
       & constants % nbasis, one, constants % vgrid2, &
       & constants % vgrid_nbasis, workspace % tmp_sph, &
       & constants % nbasis, one, workspace % tmp_grid, &
       & params % ngrid)

    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            ! Do not count internal grid points
            if(constants % ui(igrid, isph) .eq. zero) cycle
            icav = icav + 1
            phi_cav(icav) = workspace % tmp_grid(igrid, isph)
        end do
    end do

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
