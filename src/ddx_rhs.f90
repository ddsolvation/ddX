!> Routines to build rhs (phi and psi)
module ddx_rhs

use ddx_definitions
use ddx_core

implicit none

contains

!> Debug routine to test if the adjoint potential and the potential are
!> consistent.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!!
subroutine test(params, constants, workspace)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    ! local variables
    integer :: ncharges, nmultipoles, mmax, im, ic, l, m, ind, lm
    real(dp), allocatable :: charges(:), charges_coords(:, :)
    real(dp), allocatable :: multipoles(:, :), multipoles_coords(:, :)
    real(dp), allocatable :: phi(:), adj_phi(:, :)
    real(dp) :: ene, fac, c(3), box_size

    ncharges = 10
    nmultipoles = 1
    mmax = 10
    box_size = 10.0d0

    allocate(charges(ncharges), charges_coords(3, ncharges), &
        & multipoles((mmax + 1)**2, nmultipoles), &
        & multipoles_coords(3, nmultipoles), phi(ncharges), &
        & adj_phi((mmax + 1)**2, nmultipoles))

    ! setup with random numbers the coordinates and multipoles

    do im = 1, nmultipoles
        call random_number(c)
        multipoles_coords(:, im) = (c - 0.5d0)*box_size
        do lm = 1, (mmax + 1)**2
            call random_number(fac)
            multipoles(lm, im) = (fac - 0.5d0)*2.0d0
        end do
    end do

    do ic = 1, ncharges
        call random_number(c)
        charges_coords(:, im) = (c - 0.5d0)*box_size
        call random_number(fac)
        charges(ic) = (fac - 0.5d0)*2.0d0
    end do

    call build_phi_dense(params, constants, workspace, multipoles, &
        & multipoles_coords, mmax, nmultipoles, phi, charges_coords, ncharges)

    ene = 0.0d0
    do ic = 1, ncharges
        ene = ene + charges(ic)*phi(ic)
    end do
    write(6, *) 'energy multipoles -> charges:', ene 

    call build_adj_phi_dense(params, constants, workspace, charges, &
        & charges_coords, ncharges, multipoles_coords, mmax, nmultipoles, &
        & adj_phi)

    ene = 0.0d0
    do im = 1, nmultipoles
        do l = 0, mmax
            ind = l*l + l + 1
            fac = (-1.0d0)**l
            do m = -l, l
                ene = ene + fac*multipoles(ind + m, im)*adj_phi(ind + m, im)
            end do
        end do
    end do
    write(6, *) 'energy charges -> multipoles:', ene 

    deallocate(charges, charges_coords, multipoles, multipoles_coords, &
        & phi, adj_phi)
    stop
end subroutine test

!> Given a multipolar distribution, containing only charges, compute the
!> potential and its gradient at the target points using a simple N^2 code.
!> This is a test routine, used to debug the others.
!> The multipoles must be centered on the ddx spheres.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points size (ncav)
!! @param[out] e_cav: electric field at the target points
!!     size (3, ncav)
!!
subroutine test_field(params, constants, workspace, multipoles, &
        & mmax, phi_cav, e_cav)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax
    real(dp), intent(inout) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav)
    real(dp), intent(out) :: e_cav(3, constants % ncav)
    real(dp) :: c(3), v, ex, ey, ez, d2, d, d3, f
    integer :: icav, im

    f = sqrt4pi

    if (mmax .ne. 0) stop "error"
    do icav = 1, constants % ncav
        v = zero
        ex = zero
        ey = zero
        ez = zero
        do im = 1, params % nsph
            c(:) = constants % ccav(:, icav) - params % csph(:, im)

            d2 = c(1)*c(1) + c(2)*c(2) + c(3)*c(3)
            d = sqrt(d2)
            d3 = d*d2

            v = v + f*multipoles(1, im)/d
            ex = ex + f*multipoles(1, im)*c(1)/d3
            ey = ey + f*multipoles(1, im)*c(2)/d3
            ez = ez + f*multipoles(1, im)*c(3)/d3
        end do
        phi_cav(icav) = v
        e_cav(1, icav) = ex
        e_cav(2, icav) = ey
        e_cav(3, icav) = ez
    end do
end subroutine test_field

!> Given a multipolar distribution, compute the potential and its gradient
!> at the target points this is done with or without FMMs depending on the
!> relevant flag
!> The multipoles must be centered on the ddx spheres.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points size (ncav)
!! @param[out] e_cav: electric field at the target points
!!     size (3, ncav)
!!
subroutine build_e(params, constants, workspace, multipoles, &
        & mmax, phi_cav, e_cav)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax
    real(dp), intent(inout) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav)
    real(dp), intent(out) :: e_cav(3, constants % ncav)
    if (params % fmm .eq. 0) then
        call build_e_dense(params, constants, workspace, multipoles, &
            & params % csph, mmax, params % nsph, phi_cav, constants % ccav, &
            & constants % ncav, e_cav)
    else if (params % fmm .eq. 1) then
        call build_e_fmm(params, constants, workspace, multipoles, &
            & mmax, phi_cav, e_cav)
    end if
end subroutine build_e

!> Given a multipolar distribution, compute the potential and its gradient
!> at the target points using a N^2 code. 
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
!! @param[out] e_cav: electric field at the target points size (3, ncav)
!! @param[in] ccav: coordinates of the target points size (3,ncav)
!! @param[in] ncav: number of target points
!!
subroutine build_e_dense(params, constants, workspace, multipoles, cm, &
        & mmax, nm, phi_cav, ccav, ncav, e_cav)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax, nm, ncav
    real(dp), intent(inout) :: multipoles((mmax + 1)**2, nm)
    real(dp), intent(in) :: cm(3, nm)
    real(dp), intent(out) :: phi_cav(ncav)
    real(dp), intent(out) :: e_cav(3, ncav)
    real(dp), intent(in) :: ccav(3, ncav)
    real(dp), allocatable  :: tmp_m_grad(:, :, :)
    integer icav, im, l, m, i, info, indi, indj
    real(dp) :: v, ex, ey, ez, c(3), tmp1, tmp2

    ! allocate some space for the M2M gradients
    allocate(tmp_m_grad((mmax + 2)**2, 3, nm), stat=info)
    if (info .ne. 0) then
        stop "Allocation failed in build_e_dense!"
    end if

    ! call the helper routine for the M2M gradients
    call grad_m2m(params, constants, workspace, multipoles, mmax, nm, &
        & tmp_m_grad)

    ! loop over the targets and the sources and assemble the electric
    ! potential and field
    do icav = 1, ncav
        v = zero
        ex = zero
        ey = zero
        ez = zero
        do im = 1, nm
            c(:) = ccav(:, icav) - cm(:, im)
            call fmm_m2p(c, one, mmax, constants % vscales_rel, &
                & one, multipoles(:, im), one, v)

            call fmm_m2p(c, one, mmax + 1, constants % vscales_rel, &
                & one, tmp_m_grad(:, 1, im), one, ex)
            call fmm_m2p(c, one, mmax + 1, constants % vscales_rel, &
                & one, tmp_m_grad(:, 2, im), one, ey)
            call fmm_m2p(c, one, mmax + 1, constants % vscales_rel, &
                & one, tmp_m_grad(:, 3, im), one, ez)
        end do
        phi_cav(icav) = v
        e_cav(1, icav) = ex
        e_cav(2, icav) = ey
        e_cav(3, icav) = ez
    end do

    deallocate(tmp_m_grad, stat=info)
    if (info .ne. 0) then
        stop "Deallocation failed in build_e_dense!"
    end if
end subroutine build_e_dense

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

!> Given a multipolar distribution, compute the potential and its gradient
!> at the target points using FMMs.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2,nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points size (ncav)
!! @param[out] e_cav: electric field at the target points
!!     size(3, ncav)
!!
subroutine build_e_fmm(params, constants, workspace, multipoles, &
        & mmax, phi_cav, e_cav)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav)
    real(dp), intent(out) :: e_cav(3, constants % ncav)
    ! local variables
    integer :: info, isph, igrid, inode, jnode, jsph, jnear, icav
    real(dp) :: ex, ey, ez, c(3)
    real(dp), allocatable :: tmp_m_grad(:, :, :), tmp(:, :), &
        & grid_grad(:, :, :)
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform

    allocate(tmp_m_grad((mmax + 2)**2, 3, params % nsph), &
        & grid_grad(params % ngrid, 3, params % nsph), stat=info)
    if (info .ne. 0) then
        stop "Allocation failed in build_e_fmm!"
    end if

    ! compute the gradient of the m2m trasformation
    call grad_m2m(params, constants, workspace, multipoles, mmax, &
        & params % nsph, tmp_m_grad)

    ! copy the multipoles in the right places
    call load_m(params, constants, workspace, multipoles, mmax)

    ! perform the m2m, m2l and l2l steps
    call do_fmm(params, constants, workspace)

    ! near field potential (m2p)
    call tree_m2p(params, constants, params % lmax, one, &
        & workspace % tmp_sph, one, workspace % tmp_grid)

    ! far field potential, each sphere at its own points (l2p)
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & constants % nbasis, one, constants % vgrid2, &
        & constants % vgrid_nbasis, workspace % tmp_sph, &
        & constants % nbasis, one, workspace % tmp_grid, &
        & params % ngrid)

    ! near field gradients
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .eq. zero) cycle
            inode = constants % snode(isph)
            ex = zero
            ey = zero
            ez = zero
            do jnear = constants % snear(inode), constants % snear(inode+1) - 1
                jnode = constants % near(jnear)
                jsph = constants % order(constants % cluster(1, jnode))
                c = params % csph(:, isph) + constants % cgrid(:, igrid) &
                    & *params % rsph(isph) - params % csph(:, jsph)
                call fmm_m2p(c, one, mmax + 1, constants % vscales_rel, &
                    & one, tmp_m_grad(:, 1, jsph), one, ex)
                call fmm_m2p(c, one, mmax + 1, constants % vscales_rel, &
                    & one, tmp_m_grad(:, 2, jsph), one, ey)
                call fmm_m2p(c, one, mmax + 1, constants % vscales_rel, &
                    & one, tmp_m_grad(:, 3, jsph), one, ez)
            end do
            grid_grad(igrid, 1, isph) = ex
            grid_grad(igrid, 2, isph) = ey
            grid_grad(igrid, 3, isph) = ez
        end do
    end do

    ! far-field FMM gradients (only if pl > 0)
    if (params % pl .gt. 0) then
        call tree_grad_l2l(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_sph_l_grad, workspace % tmp_sph_l)
        call dgemm('T', 'N', params % ngrid, 3*params % nsph, &
            & (params % pl)**2, one, constants % vgrid2, &
            & constants % vgrid_nbasis, workspace % tmp_sph_l_grad, &
            & (params % pl+1)**2, one, grid_grad, params % ngrid)
    end if

    ! discard the internal points
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .eq. zero) cycle
            icav = icav + 1
            phi_cav(icav) = workspace % tmp_grid(igrid, isph)
            e_cav(:, icav) = grid_grad(igrid, :, isph)
        end do
    end do

    deallocate(grid_grad, stat=info)
    if (info .ne. 0) then
        stop "Deallocation failed in build_e_fmm!"
    end if

end subroutine build_e_fmm

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
    integer isph, igrid, icav

    ! copy the multipoles in the right places
    call load_m(params, constants, workspace, multipoles, mmax)

    ! perform the m2m, m2l and l2l steps
    call do_fmm(params, constants, workspace)

    ! near field
    call tree_m2p(params, constants, params % lmax, one, &
        & workspace % tmp_sph, one, workspace % tmp_grid)

    ! Potential from each sphere to its own grid points (l2p)
    call dgemm('T', 'N', params % ngrid, params % nsph, &
       & constants % nbasis, one, constants % vgrid2, &
       & constants % vgrid_nbasis, workspace % tmp_sph, &
       & constants % nbasis, one, workspace % tmp_grid, &
       & params % ngrid)

    ! discard the internal points
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .eq. zero) cycle
            icav = icav + 1
            phi_cav(icav) = workspace % tmp_grid(igrid, isph)
        end do
    end do

end subroutine build_phi_fmm

!> Given a multipolar distribution, assemble the RHS psi.
!> The multipoles must be centered on the ddx spheres.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
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

!> Given a multipolar distribution, load it into workspace % tmp_sph
!> and workspace % tmp_node_m to be used by the FMM
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2,nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!!
subroutine load_m(params, constants, workspace, multipoles, mmax)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    integer :: isph, l, ind, m, inode, tmp_max

    tmp_max = min(params % pm, mmax)

    ! P2M is not needed as we have already multipolar distributions
    ! we just have to copy the multipoles in the proper places
    do isph = 1, params % nsph
        inode = constants % snode(isph)
        workspace % tmp_sph(:, isph) = zero
        workspace % tmp_node_m(:, inode) = zero
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

end subroutine load_m

!> Given a multipolar distribution loaded in the workspace
!> perform the M2M, M2L, L2L and L2P steps.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!!
subroutine do_fmm(params, constants, workspace)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    call tree_m2m_rotation(params, constants, workspace % tmp_node_m)
    call tree_m2l_rotation(params, constants, workspace % tmp_node_m, &
        & workspace % tmp_node_l)
    call tree_l2l_rotation(params, constants, workspace % tmp_node_l)
    call tree_l2p(params, constants, one, workspace % tmp_node_l, zero, &
        & workspace % tmp_grid, workspace % tmp_sph_l)
end subroutine do_fmm

!> Given a multipolar distribution compute the action of dP on it 
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] mmax: maximum angular momentum of the multipolar distribution
!! @param[in] nm: number of multipoles
!! @param[out] tmp_m_grad: gradient of the M2M operator
!!     size ((mmax + 1)**2, 3, nm)
!!
subroutine grad_m2m(params, constants, workspace, multipoles, mmax, &
        & nm, tmp_m_grad)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax, nm
    real(dp), intent(in) :: multipoles((mmax + 1)**2, nm)
    real(dp), intent(out) :: tmp_m_grad((mmax + 2)**2, 3, nm)
    ! local variables
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    real(dp), allocatable :: tmp(:, :)
    integer :: info, im, l, indi, indj, m
    real(dp) :: tmp1, tmp2

    allocate(tmp((mmax + 2)**2, nm), stat=info)
    if (info .ne. 0) then
        stop "Allocation failed in grad_m2m!"
    end if

    zx_coord_transform = zero
    zx_coord_transform(3, 2) = one
    zx_coord_transform(2, 3) = one
    zx_coord_transform(1, 1) = one
    zy_coord_transform = zero
    zy_coord_transform(1, 2) = one
    zy_coord_transform(2, 1) = one
    zy_coord_transform(3, 3) = one

    tmp_m_grad(1:(mmax + 1)**2, 3, :) = multipoles
    do im = 1, nm
        call fmm_sph_transform(mmax, zx_coord_transform, one, &
            & multipoles(:, im), zero, tmp_m_grad(1:(mmax + 1)**2, 1, im))
        call fmm_sph_transform(mmax, zy_coord_transform, one, &
            & multipoles(:, im), zero, tmp_m_grad(1:(mmax + 1)**2, 2, im))
    end do

    do l = mmax+1, 1, -1
        indi = l*l + l + 1
        indj = indi - 2*l
        tmp1 = sqrt(dble(2*l+1)) / sqrt(dble(2*l-1))
        do m = 1-l, l-1
            tmp2 = sqrt(dble(l*l-m*m)) * tmp1
            tmp_m_grad(indi+m, :, :) = tmp2 * tmp_m_grad(indj+m, :, :)
        end do
        tmp_m_grad(indi+l, :, :) = zero
        tmp_m_grad(indi-l, :, :) = zero
    end do
    tmp_m_grad(1, :, :) = zero

    do im = 1, nm
        tmp_m_grad(:, 3, im) = tmp_m_grad(:, 3, im)
        tmp(:, im) = tmp_m_grad(:, 1, im)
        call fmm_sph_transform(mmax+1, zx_coord_transform, one, &
            & tmp(:, im), zero, tmp_m_grad(:, 1, im))
        tmp(:, im) = tmp_m_grad(:, 2, im)
        call fmm_sph_transform(mmax+1, zy_coord_transform, one, &
            & tmp(:, im), zero, tmp_m_grad(:, 2, im))
    end do

    deallocate(tmp, stat=info)
    if (info .ne. 0) then
        stop "Deallocation failed in grad_m2m!"
    end if

end subroutine grad_m2m

!> Given a multipolar distribution in real spherical harmonics and
!> centered on the spheres, compute the contributions stemming from
!> its potential to the ddX forces
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] mmax: maximum angular momentum of the multipolar distribution
!! @param[in] multipoles: 
!! @param[inout] forces: forces array size (3, nsph)
!!
subroutine grad_phi(params, constants, workspace, state, mmax, &
        & multipoles, forces, e_cav)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(inout) :: state
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(inout) :: forces(3, params % nsph)
    real(dp), intent(in) :: e_cav(3, constants % ncav)
    ! local variables
    integer :: isph, igrid, icav

    ! first contribution
    ! zeta * field
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            icav = icav + 1
            forces(:, isph) = forces(:, isph) - pt5 &
                & *state % zeta(icav)*e_cav(:, icav)
        end do
    end do

    ! second contribution

end subroutine grad_phi

!> Given a distribution of point charges at the cavity points, compute the
!> potential and its higher order derivatives up to pmax at the center of
!> the spheres. This is done with or without the FMMs depending on the
!> given flag.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points size (ncav)
!!
subroutine build_adj_phi(params, constants, workspace, charges, mmax, adj_phi)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax
    real(dp), intent(in) :: charges(constants % ncav)
    real(dp), intent(out) :: adj_phi((mmax + 1)**2, params % nsph)
    if (params % fmm .eq. 0) then
        call build_adj_phi_dense(params, constants, workspace, &
            & charges, constants % ccav, constants % ncav, params % csph, &
            & mmax, params % nsph, adj_phi)
    else if (params % fmm .eq. 1) then
        call build_adj_phi_fmm(params, constants, workspace, charges, mmax, &
            & adj_phi)
    end if
end subroutine build_adj_phi

!> Given a distribution of point charges at the cavity points, compute the
!> potential and its higher order derivatives up to pmax at the center of
!> the spheres. This is done using two loops
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!!
subroutine build_adj_phi_dense(params, constants, workspace, qcav, ccav, &
    ncav, cm, mmax, nm, adj_phi)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent (in) :: ncav, mmax, nm
    real(dp), intent(in) :: qcav(ncav), ccav(3, ncav), cm(3, nm)
    real(dp), intent(out) :: adj_phi((mmax + 1)**2, nm)
    ! local variables
    integer :: icav, im
    real(dp) :: c(3), work((mmax + 1)**2 + 3*mmax)

    do im = 1, nm
        adj_phi(:, im) = zero
        do icav = 1, ncav
            c(:) = cm(:, im) - ccav(:, icav)
            call fmm_m2p_adj_work(c, qcav(icav), one, mmax, &
                & constants % vscales_rel, one, adj_phi(:, im), work)
        end do
    end do

end subroutine build_adj_phi_dense

!> Given a distribution of point charges at the cavity points, compute the
!> potential and its higher order derivatives up to pmax at the center of
!> the spheres. This is done with FMMs.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points size (ncav)
!!
subroutine build_adj_phi_fmm(params, constants, workspace, charges, mmax, &
        & adj_phi)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: mmax
    real(dp), intent(in) :: charges(constants % ncav)
    real(dp), intent(out) :: adj_phi((mmax + 1)**2, params % nsph)
    ! local variables
    real(dp), allocatable :: tmp_grid(:, :), sph_m(:, :)
    integer :: info, icav, isph, igrid

    allocate(tmp_grid(params % ngrid, params % nsph), &
        & sph_m((mmax + 1)**2, params % nsph), stat=info)
    if (info .ne. 0) then
        stop "Allocation failed in build_adj_phi_fmm!"
    end if

    tmp_grid = zero

    ! expand the input
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if (constants % ui(igrid, isph) .eq. zero) cycle
            icav = icav + 1
            tmp_grid(igrid, isph) = charges(icav)
        end do
    end do
    adj_phi = zero

    call tree_m2p_adj(params, constants, mmax, one, tmp_grid, zero, sph_m)
    stop
    call tree_l2p_adj(params, constants, one, workspace % tmp_grid, zero, &
        & workspace % tmp_node_l, workspace % tmp_sph_l)
    call tree_l2l_rotation_adj(params, constants, workspace % tmp_node_l)
    call tree_m2l_rotation_adj(params, constants, workspace % tmp_node_l, &
        & workspace % tmp_node_m)
    call tree_m2m_rotation_adj(params, constants, workspace % tmp_node_m)

    deallocate(tmp_grid, stat=info)
    if (info .ne. 0) then
        stop "Deallocation failed in build_adj_phi_fmm"
    end if

end subroutine build_adj_phi_fmm

end module ddx_rhs
