!> Routines to build rhs (phi and psi)
module ddx_multipolar_solutes

use ddx_definitions
use ddx_core
use ddx_gradients

implicit none

contains

!> @defgroup Fortran_interface_multipolar Fortran interface: multipolar terms
!! Exposed multipolar modules in the Fortran API

!> Given a multipolar distribution, compute the required electrostatic
!! properties for the chosen model
!!
!> @ingroup Fortran_interface_multipolar
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] electrostatics, ddX electrostatic properties container
!! @param[inout] ddx_error: ddX error
!!
subroutine multipole_electrostatics(params, constants, workspace, multipoles, &
        & mmax, electrostatics, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    type(ddx_electrostatics_type), intent(out) :: electrostatics
    type(ddx_error_type), intent(inout) :: ddx_error

    call allocate_electrostatics(params, constants, electrostatics, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "allocate_electrostatics returned an error, exiting")
        return
    end if

    ! Compute the required electrostatic properties
    if (electrostatics % do_phi .and. electrostatics % do_e &
            & .and. electrostatics % do_g) then
        call multipole_electrostatics_2(params, constants, workspace, &
            & multipoles, 0, electrostatics % phi_cav, &
            & electrostatics % e_cav, electrostatics % g_cav, ddx_error)
    else if (electrostatics % do_phi .and. electrostatics % do_e) then
        call multipole_electrostatics_1(params, constants, workspace, &
            & multipoles, 0, electrostatics % phi_cav, &
            & electrostatics % e_cav, ddx_error)
    else
        call multipole_electrostatics_0(params, constants, workspace, &
            & multipoles, 0, electrostatics % phi_cav, ddx_error)
    end if

end subroutine multipole_electrostatics

!> Given a multipolar distribution, compute the potential, its gradient and
!> its hessian at the target points this is done with or without FMMs depending
!> on the relevant flag
!> The multipoles must be centered on the ddx spheres.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points, size (ncav)
!! @param[out] e_cav: electric field at the target points,
!!     size (3, ncav)
!! @param[out] g_cav: electric field gradient at the target points,
!!     size (3, 3, ncav)
!! @param[inout] ddx_error: ddX error
!!
subroutine multipole_electrostatics_2(params, constants, workspace, multipoles, &
        & mmax, phi_cav, e_cav, g_cav, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav)
    real(dp), intent(out) :: e_cav(3, constants % ncav)
    real(dp), intent(out) :: g_cav(3, 3, constants % ncav)

    if (params % fmm .eq. 0) then
        call build_g_dense(multipoles, params % csph, mmax, params % nsph, &
            & phi_cav, constants % ccav, constants % ncav, e_cav, g_cav, &
            & ddx_error)
    else if (params % fmm .eq. 1) then
        call build_g_fmm(params, constants, workspace, multipoles, &
            & mmax, phi_cav, e_cav, g_cav, ddx_error)
    end if
end subroutine multipole_electrostatics_2

!> Given a multipolar distribution, compute the potential, its gradient and
!> its hessian at the target points using a N^2 code.
!> As this routine does not use the FMM machinery it is more flexible and
!> accepts arbitrary sources and targets.
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2,nm)
!! @param[in] cm: centers of the multipolar distributions, size (3,nm)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[in] nm: number of multipoles
!! @param[out] phi_cav: electric potential at the target points, size (ncav)
!! @param[out] e_cav: electric field at the target points, size (3, ncav)
!! @param[out] g_cav: electric field at the target points, size (3, 3, ncav)
!! @param[in] ccav: coordinates of the target points, size (3,ncav)
!! @param[in] ncav: number of target points
!! @param[out] error_flag: 0 if everything went fine
!! @param[out] error_message: additional information in case of an error
!! @param[inout] ddx_error: ddX error
!!
subroutine build_g_dense(multipoles, cm, mmax, nm, phi_cav, ccav, ncav, &
        & e_cav, g_cav, ddx_error)
    implicit none
    integer, intent(in) :: mmax, nm, ncav
    real(dp), intent(in) :: multipoles((mmax + 1)**2, nm)
    real(dp), intent(in) :: cm(3, nm), ccav(3, ncav)
    real(dp), intent(out) :: phi_cav(ncav), e_cav(3, ncav), g_cav(3, 3, ncav)
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), allocatable  :: tmp_m_grad(:, :, :), vscales(:), &
        & vscales_rel(:), v4pi2lp1(:), tmp_m_hess(:, :, :, :), tmp(:, :)
    integer icav, im, info
    real(dp) :: v, ex, ey, ez, c(3), gxx, gxy, gxz, gyy, gyz, gzz

    ! allocate some space for the M2M gradients and precompute
    ! the quantities for the m2p
    allocate(tmp_m_grad((mmax + 2)**2, 3, nm), vscales((mmax + 3)**2), &
        & vscales_rel((mmax + 3)**2), v4pi2lp1(mmax + 3), &
        & tmp_m_hess((mmax + 3)**2, 3, nm, 3), tmp((mmax + 2)**2, 3), &
        & stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Allocation failed in build_g_dense!")
        return
    end if
    call ylmscale(mmax + 2, vscales, v4pi2lp1, vscales_rel)

    ! call the helper routine for the M2M gradient and hessian
    call grad_m2m(multipoles, mmax, nm, tmp_m_grad, ddx_error)
    tmp = tmp_m_grad(:, 1, :)
    call grad_m2m(tmp, mmax + 1, nm, tmp_m_hess(:, :, :, 1), ddx_error)
    tmp = tmp_m_grad(:, 2, :)
    call grad_m2m(tmp, mmax + 1, nm, tmp_m_hess(:, :, :, 2), ddx_error)
    tmp = tmp_m_grad(:, 3, :)
    call grad_m2m(tmp, mmax + 1, nm, tmp_m_hess(:, :, :, 3), ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "grad_m2m returned an error, exiting")
        return
    end if

    ! loop over the targets and the sources and assemble the electric
    ! potential and field
    do icav = 1, ncav
        v = zero
        ex = zero
        ey = zero
        ez = zero
        gxx = zero
        gxy = zero
        gxz = zero
        gyy = zero
        gyz = zero
        gzz = zero
        do im = 1, nm
            c(:) = ccav(:, icav) - cm(:, im)
            call fmm_m2p(c, one, mmax, vscales_rel, one, multipoles(:, im), &
                & one, v)

            call fmm_m2p(c, one, mmax + 1, vscales_rel, one, &
                & tmp_m_grad(:, 1, im), one, ex)
            call fmm_m2p(c, one, mmax + 1, vscales_rel, one, &
                & tmp_m_grad(:, 2, im), one, ey)
            call fmm_m2p(c, one, mmax + 1, vscales_rel, one, &
                & tmp_m_grad(:, 3, im), one, ez)

            call fmm_m2p(c, one, mmax + 2, vscales_rel, one, &
                & tmp_m_hess(:, 1, im, 1), one, gxx)
            call fmm_m2p(c, one, mmax + 2, vscales_rel, one, &
                & tmp_m_hess(:, 2, im, 1), one, gxy)
            call fmm_m2p(c, one, mmax + 2, vscales_rel, one, &
                & tmp_m_hess(:, 3, im, 1), one, gxz)

            call fmm_m2p(c, one, mmax + 2, vscales_rel, one, &
                & tmp_m_hess(:, 2, im, 2), one, gyy)
            call fmm_m2p(c, one, mmax + 2, vscales_rel, one, &
                & tmp_m_hess(:, 3, im, 2), one, gyz)

            call fmm_m2p(c, one, mmax + 2, vscales_rel, one, &
                & tmp_m_hess(:, 3, im, 3), one, gzz)
        end do
        phi_cav(icav) = v
        e_cav(1, icav) = ex
        e_cav(2, icav) = ey
        e_cav(3, icav) = ez
        g_cav(1, 1, icav) = gxx
        g_cav(1, 2, icav) = gxy
        g_cav(1, 3, icav) = gxz
        g_cav(2, 1, icav) = gxy
        g_cav(2, 2, icav) = gyy
        g_cav(2, 3, icav) = gyz
        g_cav(3, 1, icav) = gxz
        g_cav(3, 2, icav) = gyz
        g_cav(3, 3, icav) = gzz
    end do

    deallocate(tmp_m_grad, vscales, vscales_rel, v4pi2lp1, tmp_m_hess, &
        & tmp, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Deallocation failed in build_e_dense!")
        return
    end if
end subroutine build_g_dense

!> Given a multipolar distribution, compute the potential and its gradient
!> at the target points this is done with or without FMMs depending on the
!> relevant flag
!> The multipoles must be centered on the ddx spheres.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points, size (ncav)
!! @param[out] e_cav: electric field at the target points,
!!     size (3, ncav)
!! @param[inout] ddx_error: ddX error
!!
subroutine multipole_electrostatics_1(params, constants, workspace, multipoles, &
        & mmax, phi_cav, e_cav, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav)
    real(dp), intent(out) :: e_cav(3, constants % ncav)
    if (params % fmm .eq. 0) then
        call build_e_dense(multipoles, params % csph, mmax, params % nsph, &
            & phi_cav, constants % ccav, constants % ncav, e_cav, ddx_error)
    else if (params % fmm .eq. 1) then
        call build_e_fmm(params, constants, workspace, multipoles, &
            & mmax, phi_cav, e_cav, ddx_error)
    end if
end subroutine multipole_electrostatics_1

!> Given a multipolar distribution, compute the potential and its gradient
!> at the target points using a N^2 code. 
!> As this routine does not use the FMM machinery it is more flexible and
!> accepts arbitrary sources and targets.
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2,nm)
!! @param[in] cm: centers of the multipolar distributions, size (3,nm)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[in] nm: number of multipoles
!! @param[out] phi_cav: electric potential at the target points, size (ncav)
!! @param[out] e_cav: electric field at the target points, size (3, ncav)
!! @param[in] ccav: coordinates of the target points, size (3,ncav)
!! @param[in] ncav: number of target points
!! @param[inout] ddx_error: ddX error
!!
subroutine build_e_dense(multipoles, cm, mmax, nm, phi_cav, ccav, ncav, &
        & e_cav, ddx_error)
    implicit none
    integer, intent(in) :: mmax, nm, ncav
    real(dp), intent(in) :: cm(3, nm), ccav(3, ncav), &
        & multipoles((mmax + 1)**2, nm)
    real(dp), intent(out) :: phi_cav(ncav), e_cav(3, ncav)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! local variables
    real(dp), allocatable  :: tmp_m_grad(:, :, :), vscales(:), &
        & vscales_rel(:), v4pi2lp1(:)
    integer icav, im, info
    real(dp) :: v, ex, ey, ez, c(3)

    ! allocate some space for the M2M gradients and precompute
    ! the quantities for the m2p
    allocate(tmp_m_grad((mmax + 2)**2, 3, nm), vscales((mmax + 2)**2), &
        & vscales_rel((mmax + 2)**2), v4pi2lp1(mmax + 2), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Allocation failed in build_e_dense!")
        return
    end if
    call ylmscale(mmax + 1, vscales, v4pi2lp1, vscales_rel)

    ! call the helper routine for the M2M gradients
    call grad_m2m(multipoles, mmax, nm, tmp_m_grad, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "grad_m2m returned an error, exiting")
        return
    end if


    ! loop over the targets and the sources and assemble the electric
    ! potential and field
    do icav = 1, ncav
        v = zero
        ex = zero
        ey = zero
        ez = zero
        do im = 1, nm
            c(:) = ccav(:, icav) - cm(:, im)
            call fmm_m2p(c, one, mmax, vscales_rel, one, multipoles(:, im), &
                & one, v)

            call fmm_m2p(c, one, mmax + 1, vscales_rel, one, &
                & tmp_m_grad(:, 1, im), one, ex)
            call fmm_m2p(c, one, mmax + 1, vscales_rel, one, &
                & tmp_m_grad(:, 2, im), one, ey)
            call fmm_m2p(c, one, mmax + 1, vscales_rel, one, &
                & tmp_m_grad(:, 3, im), one, ez)
        end do
        phi_cav(icav) = v
        e_cav(1, icav) = ex
        e_cav(2, icav) = ey
        e_cav(3, icav) = ez
    end do

    deallocate(tmp_m_grad, vscales, vscales_rel, v4pi2lp1, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Deallocation failed in build_e_dense!")
        return
    end if
end subroutine build_e_dense

!> Given a multipolar distribution, compute the potential at the target points
!> this is done with or without FMMs depending on the relevant flag
!> The multipoles must be centered on the ddx spheres.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points, size (ncav)
!! @param[inout] ddx_error: ddX error
!!
subroutine multipole_electrostatics_0(params, constants, workspace, multipoles, &
        & mmax, phi_cav, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav)
    if (params % fmm .eq. 0) then
        call build_phi_dense(multipoles, params % csph, mmax, params % nsph, &
            & phi_cav, constants % ccav, constants % ncav, ddx_error)
    else if (params % fmm .eq. 1) then
        call build_phi_fmm(params, constants, workspace, multipoles, mmax, &
            & phi_cav)
    end if
end subroutine multipole_electrostatics_0

!> Given a multipolar distribution, compute the potential at the target points
!> using a N^2 code.
!> As this routine does not use the FMM machinery it is more flexible and
!> accepts arbitrary sources and targets.
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2,nm)
!! @param[in] cm: centers of the multipolar distributions, size (3,nm)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[in] nm: number of multipoles
!! @param[out] phi_cav: electric potential at the target points, size (ncav)
!! @param[in] ccav: coordinates of the target points, size (3,ncav)
!! @param[in] ncav: number of target points
!! @param[inout] ddx_error: ddX error
!!
subroutine build_phi_dense(multipoles, cm, mmax, nm, phi_cav, ccav, ncav, &
        & ddx_error)
    implicit none
    integer, intent(in) :: mmax, nm, ncav
    real(dp), intent(in) :: cm(3, nm), ccav(3, ncav), &
        & multipoles((mmax + 1)**2, nm)
    real(dp), intent(out) :: phi_cav(ncav)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! local variables
    integer icav, im, info
    real(dp) :: v, c(3)
    real(dp), allocatable :: vscales(:), vscales_rel(:), v4pi2lp1(:)

    ! precompute the quantities for the m2p
    allocate(vscales((mmax + 1)**2), vscales_rel((mmax + 1)**2), &
        & v4pi2lp1(mmax + 1), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, 'Allocation failed in build_phi_dense!')
        return
    end if
    call ylmscale(mmax, vscales, v4pi2lp1, vscales_rel)

    do icav = 1, ncav
        v = zero
        do im = 1, nm
            c(:) = ccav(:, icav) - cm(:, im)
            call fmm_m2p(c, one, mmax, vscales_rel, one, &
                & multipoles(:, im), one, v)
        end do
        phi_cav(icav) = v
    end do

    deallocate(vscales, vscales_rel, v4pi2lp1, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, 'Deallocation failed in build_phi_dense!')
        return
    end if

end subroutine build_phi_dense

!> Given a multipolar distribution, compute the potential and its gradient
!> at the target points using FMMs.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2,nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points, size (ncav)
!! @param[out] e_cav: electric field at the target points,
!!     size(3, ncav)
!! @param[inout] ddx_error: ddX error
!!
subroutine build_e_fmm(params, constants, workspace, multipoles, &
        & mmax, phi_cav, e_cav, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav), &
        & e_cav(3, constants % ncav)
    ! local variables
    integer :: info, isph, igrid, inode, jnode, jsph, jnear, icav
    real(dp) :: ex, ey, ez, c(3)
    real(dp), allocatable :: tmp_m_grad(:, :, :), grid_grad(:, :, :)

    allocate(tmp_m_grad((mmax + 2)**2, 3, params % nsph), &
        & grid_grad(params % ngrid, 3, params % nsph), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Allocation failed in build_e_fmm!")
        return
    end if

    ! compute the gradient of the m2m trasformation
    call grad_m2m(multipoles, mmax, params % nsph, tmp_m_grad, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "grad_m2m returned an error, exiting")
        return
    end if

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

    deallocate(tmp_m_grad, grid_grad, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Deallocation failed in build_e_fmm!")
        return
    end if

end subroutine build_e_fmm

!> Given a multipolar distribution, compute the potential, the field
!> and the field gradient at the target points using FMMs.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2,nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points, size (ncav)
!! @param[out] e_cav: electric field at the target points,
!!     size(3, ncav)
!! @param[out] g_cav: electric field gradient at the target points,
!!     size(3, 3, ncav)
!! @param[inout] ddx_error: ddX error
!!
subroutine build_g_fmm(params, constants, workspace, multipoles, &
        & mmax, phi_cav, e_cav, g_cav, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(out) :: phi_cav(constants % ncav), &
        & e_cav(3, constants % ncav), g_cav(3, 3, constants % ncav)
    ! local variables
    integer :: info, isph, igrid, inode, jnode, jsph, jnear, icav
    real(dp) :: ex, ey, ez, c(3), gxx, gxy, gxz, gyy, gyz, gzz
    real(dp), allocatable :: m_grad(:, :, :), grid_grad(:, :, :), &
        & m_hess_x(:, :, :), m_hess_y(:, :, :), m_hess_z(:, :, :), &
        & m_grad_component(:, :), grid_hess_x(:, :, :), &
        & grid_hess_y(:, :, :), grid_hess_z(:, :, :), grid_hessian2(:, :, :)

    allocate(m_hess_x((mmax+3)**2, 3, params % nsph), &
        & m_hess_y((mmax+3)**2, 3, params % nsph), &
        & m_hess_z((mmax+3)**2, 3, params % nsph), &
        & m_grad_component((mmax+2)**2, params % nsph), &
        & m_grad((mmax + 2)**2, 3, params % nsph), &
        & grid_grad(params % ngrid, 3, params % nsph), &
        & grid_hess_x(params % ngrid, 3, params % nsph), &
        & grid_hess_y(params % ngrid, 3, params % nsph), &
        & grid_hess_z(params % ngrid, 3, params % nsph), &
        & grid_hessian2(params % ngrid, 3, params % nsph), &
        & stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Allocation failed in build_g_fmm!")
        return
    end if

    ! compute the gradient of the m2m trasformation
    call grad_m2m(multipoles, mmax, params % nsph, m_grad, ddx_error)

    ! starting from the previously computed gradient, compute the
    ! hessian of the target multipoles
    m_grad_component = m_grad(:, 1, :)
    call grad_m2m(m_grad_component, mmax + 1, params % nsph, m_hess_x, ddx_error)
    m_grad_component = m_grad(:, 2, :)
    call grad_m2m(m_grad_component, mmax + 1, params % nsph, m_hess_y, ddx_error)
    m_grad_component = m_grad(:, 3, :)
    call grad_m2m(m_grad_component, mmax + 1, params % nsph, m_hess_z, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "grad_m2m returned an error, exiting")
        return
    end if


    ! copy the multipoles to the right places
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
                    & one, m_grad(:, 1, jsph), one, ex)
                call fmm_m2p(c, one, mmax + 1, constants % vscales_rel, &
                    & one, m_grad(:, 2, jsph), one, ey)
                call fmm_m2p(c, one, mmax + 1, constants % vscales_rel, &
                    & one, m_grad(:, 3, jsph), one, ez)
            end do
            grid_grad(igrid, 1, isph) = ex
            grid_grad(igrid, 2, isph) = ey
            grid_grad(igrid, 3, isph) = ez
        end do
    end do

    ! near field hessians
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .eq. zero) cycle
            inode = constants % snode(isph)
            gxx = zero
            gxy = zero
            gxz = zero
            gyy = zero
            gyz = zero
            gzz = zero
            do jnear = constants % snear(inode), constants % snear(inode+1) - 1
                jnode = constants % near(jnear)
                jsph = constants % order(constants % cluster(1, jnode))
                c = params % csph(:, isph) + constants % cgrid(:, igrid) &
                    & *params % rsph(isph) - params % csph(:, jsph)
                call fmm_m2p(c, one, mmax + 2, constants % vscales_rel, &
                    & one, m_hess_x(:, 1, jsph), one, gxx)
                call fmm_m2p(c, one, mmax + 2, constants % vscales_rel, &
                    & one, m_hess_x(:, 2, jsph), one, gxy)
                call fmm_m2p(c, one, mmax + 2, constants % vscales_rel, &
                    & one, m_hess_x(:, 3, jsph), one, gxz)
                call fmm_m2p(c, one, mmax + 2, constants % vscales_rel, &
                    & one, m_hess_y(:, 2, jsph), one, gyy)
                call fmm_m2p(c, one, mmax + 2, constants % vscales_rel, &
                    & one, m_hess_y(:, 3, jsph), one, gyz)
                call fmm_m2p(c, one, mmax + 2, constants % vscales_rel, &
                    & one, m_hess_z(:, 3, jsph), one, gzz)
            end do
            grid_hess_x(igrid, 1, isph) = gxx
            grid_hess_x(igrid, 2, isph) = gxy
            grid_hess_x(igrid, 3, isph) = gxz
            grid_hess_y(igrid, 1, isph) = gxy
            grid_hess_y(igrid, 2, isph) = gyy
            grid_hess_y(igrid, 3, isph) = gyz
            grid_hess_z(igrid, 1, isph) = gxz
            grid_hess_z(igrid, 2, isph) = gyz
            grid_hess_z(igrid, 3, isph) = gzz
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
        call tree_grad_l2l(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_sph_l_grad, workspace % tmp_sph_l)
    end if

    ! Far field FMM hessians
    if (params % pl .gt. 1) then
        ! Load previously computed gradient into leaves, since
        ! tree_grad_l2l currently takes local expansions of entire
        ! tree. In future it might be changed.
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_l(:, inode) = &
                & workspace % tmp_sph_l_grad(:, 1, isph)
        end do
        ! Get gradient of a gradient of L2L. Currently this uses input
        ! pl maximal degree of local harmonics but in reality we need
        ! only pl-1 maximal degree since coefficients of harmonics of a
        ! degree pl are zeros.
        call tree_grad_l2l(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_sph_l_grad2, workspace % tmp_sph_l)
        ! Apply L2P for every axis
        call dgemm('T', 'N', params % ngrid, 3*params % nsph, &
            & (params % pl-1)**2, one, constants % vgrid2, &
            & constants % vgrid_nbasis, workspace % tmp_sph_l_grad2, &
            & (params % pl+1)**2, zero, grid_hessian2, params % ngrid)
        ! Properly copy hessian
        grid_hess_x = grid_hess_x + grid_hessian2

        ! same for y
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_l(:, inode) = &
                & workspace % tmp_sph_l_grad(:, 2, isph)
        end do
        call tree_grad_l2l(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_sph_l_grad2, workspace % tmp_sph_l)
        call dgemm('T', 'N', params % ngrid, 3*params % nsph, &
            & (params % pl-1)**2, one, constants % vgrid2, &
            & constants % vgrid_nbasis, workspace % tmp_sph_l_grad2, &
            & (params % pl+1)**2, zero, grid_hessian2, params % ngrid)
        grid_hess_y = grid_hess_y + grid_hessian2

        ! same for z
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_l(:, inode) = &
                & workspace % tmp_sph_l_grad(:, 3, isph)
        end do
        call tree_grad_l2l(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_sph_l_grad2, workspace % tmp_sph_l)
        call dgemm('T', 'N', params % ngrid, 3*params % nsph, &
            & (params % pl-1)**2, one, constants % vgrid2, &
            & constants % vgrid_nbasis, workspace % tmp_sph_l_grad2, &
            & (params % pl+1)**2, zero, grid_hessian2, params % ngrid)
        grid_hess_z = grid_hess_z + grid_hessian2

    end if

    ! discard the internal points
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .eq. zero) cycle
            icav = icav + 1
            phi_cav(icav) = workspace % tmp_grid(igrid, isph)
            e_cav(:, icav) = grid_grad(igrid, :, isph)
            g_cav(1, :, icav) = grid_hess_x(igrid, :, isph)
            g_cav(2, :, icav) = grid_hess_y(igrid, :, isph)
            g_cav(3, :, icav) = grid_hess_z(igrid, :, isph)
        end do
    end do

    deallocate(m_hess_x, m_hess_y, m_hess_z, m_grad_component, m_grad, &
        & grid_grad, grid_hess_x, grid_hess_y, grid_hess_z, grid_hessian2, &
        & stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Deallocation failed in build_g_fmm!")
        return
    end if

end subroutine build_g_fmm

!> Given a multipolar distribution, compute the potential at the target points
!> using FMMs.
!! @param[in] params: ddx parameters
!! @param[in]  constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2,nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] phi_cav: electric potential at the target points, size (ncav)
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

!> @ingroup Fortran_interface_multipolar
!> Given a multipolar distribution, assemble the RHS psi.
!> The multipoles must be centered on the ddx spheres.
!! @param[in] params: ddx parameters
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2,nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[out] psi: RHS for adjoint linear systems,
!!     size ((lmax+1)**2,nsph), the internal lmax should be >= mmax
!!
subroutine multipole_psi(params, multipoles, mmax, psi)
    implicit none
    type(ddx_params_type), intent(in) :: params
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
end subroutine multipole_psi

!> Given a multipolar distribution, load it into workspace % tmp_sph
!> and workspace % tmp_node_m to be used by the FMM
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
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
!! @param[inout] ddx_error: ddX error
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

!> Given a multipolar distribution compute the action of dP on it, this
!> is required in the computation of the forces.
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipolar distribution
!! @param[in] nm: number of multipoles
!! @param[out] tmp_m_grad: gradient of the M2M operator,
!!     size ((mmax + 1)**2, 3, nm)
!! @param[inout] ddx_error: ddX error
!!
subroutine grad_m2m(multipoles, mmax, nm, tmp_m_grad, ddx_error)
    implicit none
    integer, intent(in) :: mmax, nm
    real(dp), intent(in) :: multipoles((mmax + 1)**2, nm)
    real(dp), intent(out) :: tmp_m_grad((mmax + 2)**2, 3, nm)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! local variables
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    real(dp), allocatable :: tmp(:, :)
    integer :: info, im, l, indi, indj, m
    real(dp) :: tmp1, tmp2

    allocate(tmp((mmax + 2)**2, nm), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Allocation failed in grad_m2m!")
        return
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
        call update_error(ddx_error, "Deallocation failed in grad_m2m!")
        return
    end if

end subroutine grad_m2m

!> Given a charge distribution centered on the spheres, compute the
!> contributions to the forces stemming from its electrostatic interactions.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] charges: charges, size (nsph)
!! @param[inout] forces: forces array, size (3, nsph)
!! @param[inout] ddx_error: ddX error
!!
subroutine grad_phi_for_charges(params, constants, workspace, state, &
        & charges, forces, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: charges(params % nsph)
    real(dp), intent(inout) :: forces(3, params % nsph)
    ! local variables
    integer :: info
    real(dp), allocatable :: multipoles(:, :)

    allocate(multipoles(1, params % nsph), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Allocation failed in grad_phi_for_charges")
        return
    end if
    ! convert the charges to multipoles
    multipoles(1, :) = charges/sqrt4pi
    call grad_phi(params, constants, workspace, state, 0, multipoles, forces, &
        & ddx_error)
    deallocate(multipoles, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Deallocation failed in grad_phi_for_charges")
        return
    end if

end subroutine grad_phi_for_charges

!> Given a multipolar distribution in real spherical harmonics and
!> centered on the spheres, compute the contributions to the forces
!> stemming from its electrostatic interactions.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[inout] state: ddx state
!! @param[in] mmax: maximum angular momentum of the multipolar distribution
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2, nsph)
!! @param[inout] forces: forces array, size (3, nsph)
!! @param[inout] ddx_error: ddX error
!!
subroutine grad_phi(params, constants, workspace, state, mmax, &
        & multipoles, forces, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(inout) :: forces(3, params % nsph)
    ! local variables
    integer :: info, im, lm, ind, l, m
    real(dp), allocatable :: adj_phi(:, :), m_grad(:, :, :)
    real(dp) :: fac

    ! get some space for the adjoint potential, note that we need it
    ! up to mmax + 1 as we are doing derivatives
    allocate(adj_phi((mmax + 2)**2, params % nsph), &
        & m_grad((mmax + 2)**2, 3, params % nsph), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Allocation failed in grad_phi!")
        return
    end if

    ! build the adjoint potential
    call build_adj_phi(params, constants, workspace, state % zeta, &
        & mmax + 1, adj_phi, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "build_adj_phi returned an error, exiting")
        return
    end if

    ! build the gradient of the M2M transformation
    call grad_m2m(multipoles, mmax, params % nsph, m_grad, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "grad_m2m returned an error, exiting")
        return
    end if

    ! contract the two ingredients to build the second contribution
    do im = 1, params % nsph
        do l = 1, mmax + 1
            ind = l*l + l + 1
            fac = (-one)**(l+1)
            do m = -l, l
                lm = ind + m
                forces(1, im) = forces(1, im) &
                    & + fac*pt5*m_grad(lm, 1, im)*adj_phi(lm, im)
                forces(2, im) = forces(2, im) &
                    & + fac*pt5*m_grad(lm, 2, im)*adj_phi(lm, im)
                forces(3, im) = forces(3, im) &
                    & + fac*pt5*m_grad(lm, 3, im)*adj_phi(lm, im)
            end do
        end do
    end do

    deallocate(adj_phi, m_grad, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Deallocation failed in grad_phi!")
        return
    end if

end subroutine grad_phi

!> Given a distribution of point charges at the cavity points, compute the
!> potential and its higher order derivatives up to pmax at the center of
!> the spheres. This is done with or without the FMMs depending on the
!> given flag. This is used in the computation of the forces.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] charges: charges (zeta) at cavity points, size(ncav)
!! @param[in] mmax: maximum angular momentum of the target multipoles
!! @param[out] adj_phi: electrostatic properties up to order mmax at
!!     the target multipoles, size ((mmax+1)**2, nsph)
!! @param[inout] ddx_error: ddX error
!!
subroutine build_adj_phi(params, constants, workspace, charges, mmax, &
        & adj_phi, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: charges(constants % ncav)
    real(dp), intent(out) :: adj_phi((mmax + 1)**2, params % nsph)
    if (params % fmm .eq. 0) then
        call build_adj_phi_dense(charges, constants % ccav, constants % ncav, &
            & params % csph, mmax, params % nsph, adj_phi, ddx_error)
    else if (params % fmm .eq. 1) then
        call build_adj_phi_fmm(params, constants, workspace, charges, mmax, &
            & adj_phi, ddx_error)
    end if
end subroutine build_adj_phi

!> Given a distribution of point charges at the cavity points, compute the
!> potential and its higher order derivatives up to pmax at the center of
!> the spheres. This is done using two loops and required in the computation
!> of the forces.
!! @param[in] qcav: charges, size(ncav)
!! @param[in] ccav: coordinates of the charges, size(3,ncav)
!! @param[in] ncav: number of charges
!! @param[in] cm: coordinates of the multipoles, size(3,nm)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[in] nm: number of multipoles
!! @param[out] adj_phi: adjoint potential, size((mmax+1)**2,nm)
!! @param[inout] ddx_error: ddX error
!!
subroutine build_adj_phi_dense(qcav, ccav, ncav, cm, mmax, nm, adj_phi, &
        & ddx_error)
    implicit none
    integer, intent (in) :: ncav, mmax, nm
    real(dp), intent(in) :: qcav(ncav), ccav(3, ncav), cm(3, nm)
    real(dp), intent(out) :: adj_phi((mmax + 1)**2, nm)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! local variables
    integer :: icav, im, info
    real(dp) :: c(3)
    real(dp), allocatable :: vscales(:), vscales_rel(:), v4pi2lp1(:), &
        & work(:)

    allocate(vscales((mmax + 1)**2), vscales_rel((mmax + 1)**2), &
        & v4pi2lp1(mmax + 1), work((mmax + 1)**2 + 3*mmax), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, 'Allocation failed in build_adj_phi_dense!')
        return
    end if
    call ylmscale(mmax, vscales, v4pi2lp1, vscales_rel)

    do im = 1, nm
        adj_phi(:, im) = zero
        do icav = 1, ncav
            c(:) = cm(:, im) - ccav(:, icav)
            call fmm_m2p_adj_work(c, qcav(icav), one, mmax, vscales_rel, &
                & one, adj_phi(:, im), work)
        end do
    end do

    deallocate(vscales, vscales_rel, v4pi2lp1, work, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, 'Deallocation failed in build_adj_phi_dense!')
        return
    end if

end subroutine build_adj_phi_dense

!> Given a distribution of point charges at the cavity points, compute the
!> potential and its higher order derivatives up to pmax at the center of
!> the spheres. This is done with FMMs and required in the computation of
!> the forces.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] charges: charges (zeta) at cavity points, size(ncav)
!! @param[in] mmax: maximum angular momentum of the target multipoles
!! @param[out] adj_phi: electrostatic properties up to order mmax at
!!     the target multipoles, size ((mmax+1)**2, nsph)
!! @param[inout] ddx_error: ddX error
!!
subroutine build_adj_phi_fmm(params, constants, workspace, charges, mmax, &
        & adj_phi, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: charges(constants % ncav)
    real(dp), intent(out) :: adj_phi((mmax + 1)**2, params % nsph)
    ! local variables
    real(dp), allocatable :: tmp_grid(:, :), sph_m(:, :), work(:)
    integer :: info, icav, isph, igrid, inode, indl, l, jnear, &
        & jnode, jsph, m, ind
    real(dp) :: c(3), fac

    allocate(tmp_grid(params % ngrid, params % nsph), &
        & sph_m((mmax + 1)**2, params % nsph), &
        & work((mmax + 1)**2 + 3*mmax), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Allocation failed in build_adj_phi_fmm!")
        return
    end if

    tmp_grid = zero

    ! expand the input from cavity points to all the Lebedev points
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if (constants % ui(igrid, isph) .eq. zero) cycle
            icav = icav + 1
            tmp_grid(igrid, isph) = charges(icav)
        end do
    end do
    adj_phi = zero

    ! near field. We do not call tree_m2p_adj as this routine skips the
    ! isph = jsph case which instead is needed.
    do isph = 1, params % nsph
        inode = constants % snode(isph)
        do jnear = constants % snear(inode), constants % snear(inode+1)-1
            jnode = constants % near(jnear)
            jsph = constants % order(constants % cluster(1, jnode))
            do igrid = 1, params % ngrid
                c = constants % cgrid(:, igrid)*params % rsph(isph) - &
                    & params % csph(:, jsph) + params % csph(:, isph)
                call fmm_m2p_adj_work(c, tmp_grid(igrid, isph), &
                    & params % rsph(jsph), mmax, constants % vscales_rel, one, &
                    & adj_phi(:, jsph), work)
            end do
        end do
    end do

    ! compute the far field and store it in tmp_node_m
    call tree_l2p_adj(params, constants, one, tmp_grid, zero, &
        & workspace % tmp_node_l, workspace % tmp_sph_l)
    call tree_l2l_rotation_adj(params, constants, workspace % tmp_node_l)
    call tree_m2l_rotation_adj(params, constants, workspace % tmp_node_l, &
        & workspace % tmp_node_m)
    call tree_m2m_rotation_adj(params, constants, workspace % tmp_node_m)

    ! move the far field from tmp_node_m to adj_phi
    if(mmax .lt. params % pm) then
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            adj_phi(:, isph) = adj_phi(:, isph) &
                & + workspace % tmp_node_m(1:(mmax + 1)**2, inode)
        end do
    else
        indl = (params % pm+1)**2
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            adj_phi(1:indl, isph) = adj_phi(1:indl, isph) &
                & + workspace % tmp_node_m(:, inode)
        end do
    end if

    ! scale by a factor to make it consistent with the non FMMized adjoint
    ! potential
    do isph = 1, params % nsph
        do l = 0, mmax
           ind = l*l + l + 1
           fac = - (-one/params % rsph(isph))**(l + 1)
           do m = -l, l
               adj_phi(ind + m, isph) = fac*adj_phi(ind + m, isph)
           end do
       end do
    end do

    deallocate(tmp_grid, work, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Deallocation failed in build_adj_phi_fmm")
        return
    end if

end subroutine build_adj_phi_fmm

!> Given a charge distribution centered on the spheres, compute the
!> contributions to the forces stemming from its electrostatic interactions.
!!
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[inout] state: ddx state
!! @param[in] charges: charges, size (nsph)
!! @param[inout] forces: forces array, size (3, nsph)
!! @param[inout] ddx_error: ddX error
!!
subroutine grad_e_for_charges(params, constants, workspace, state, &
        & charges, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: charges(params % nsph)
    real(dp), intent(inout) :: force(3, params % nsph)
    real(dp), allocatable :: multipoles(:, :)
    integer :: info

    allocate(multipoles(1, params % nsph), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Allocation failed in grad_e_for_charges")
        return
    end if
    ! convert the charges to multipoles
    multipoles(1, :) = charges/sqrt4pi
    call grad_e(params, constants, workspace, state, 0, multipoles, force, &
        & ddx_error)
    deallocate(multipoles, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, "Deallocation failed in grad_e_for_charges")
        return
    end if

end subroutine grad_e_for_charges

!> Given a multipolar distribution in real spherical harmonics and
!> centered on the spheres, compute the contributions to the forces
!> stemming from its electrostatic interactions in case of ddLPB F RHS.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[inout] state: ddx state
!! @param[in] mmax: maximum angular momentum of the multipolar distribution
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2, nsph)
!! @param[inout] forces: forces array, size (3, nsph)
!! @param[inout] ddx_error: ddX error
!!
subroutine grad_e(params, constants, workspace, state, mmax, &
        & multipoles, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax+1)**2, params % nsph)
    real(dp), intent(inout) :: force(3, params % nsph)

    real(dp), allocatable :: zeta_comp(:),  m_grad(:, :, :), &
        & m_hess_x(:, :, :), m_hess_y(:, :, :), m_hess_z(:, :, :), &
        & m_grad_component(:, :), adj_phi_x(:, :), adj_phi_y(:, :), &
        & adj_phi_z(:, :)
    integer :: info, isph, ind, l, m, lm
    real(dp) :: fac

    allocate(zeta_comp(constants % ncav), &
        & adj_phi_x((mmax+3)**2, params % nsph), &
        & adj_phi_y((mmax+3)**2, params % nsph), &
        & adj_phi_z((mmax+3)**2, params % nsph), &
        & m_hess_x((mmax+3)**2, 3, params % nsph), &
        & m_hess_y((mmax+3)**2, 3, params % nsph), &
        & m_hess_z((mmax+3)**2, 3, params % nsph), &
        & m_grad((mmax+2)**2, 3, params % nsph), &
        & m_grad_component((mmax+2)**2, params % nsph), &
        & stat=info)
    if (info.ne.0) then
        call update_error(ddx_error, "Allocation failed in grad_e")
        return
    end if

    ! compute the gradient of the target charges
    call grad_m2m(multipoles, mmax, params % nsph, m_grad, ddx_error)

    ! starting from the previously computed gradient, compute the
    ! hessian of the target multipoles
    m_grad_component = m_grad(:, 1, :)
    call grad_m2m(m_grad_component, mmax + 1, params % nsph, m_hess_x, ddx_error)
    m_grad_component = m_grad(:, 2, :)
    call grad_m2m(m_grad_component, mmax + 1, params % nsph, m_hess_y, ddx_error)
    m_grad_component = m_grad(:, 3, :)
    call grad_m2m(m_grad_component, mmax + 1, params % nsph, m_hess_z, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "grad_m2m returned an error, exiting")
        return
    end if


    ! now we compute the adjoint potential of the x, y and z components
    ! of the input dipoles separately
    zeta_comp = state % zeta_dip(1, :)
    call build_adj_phi(params, constants, workspace, zeta_comp, &
        & mmax + 2, adj_phi_x, ddx_error)
    zeta_comp = state % zeta_dip(2, :)
    call build_adj_phi(params, constants, workspace, zeta_comp, &
        & mmax + 2, adj_phi_y, ddx_error)
    zeta_comp = state % zeta_dip(3, :)
    call build_adj_phi(params, constants, workspace, zeta_comp, &
        & mmax + 2, adj_phi_z, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "build_adj_phi returned an error, exiting")
        return
    end if


    ! contract the hessian of the multipoles with the three adjoint
    ! potentials to get the forces contributions
    do isph = 1, params % nsph
        do l = 2, mmax + 2
            ind = l*l + l + 1
            fac = -(-one)**(l+1)/2.0d0
            do m = -l, l
                lm = ind + m
                force(1, isph) = force(1, isph) &
                    & + fac*m_hess_x(lm, 1, isph)*adj_phi_x(lm, isph) &
                    & + fac*m_hess_y(lm, 1, isph)*adj_phi_y(lm, isph) &
                    & + fac*m_hess_z(lm, 1, isph)*adj_phi_z(lm, isph)
                force(2, isph) = force(2, isph) &
                    & + fac*m_hess_x(lm, 2, isph)*adj_phi_x(lm, isph) &
                    & + fac*m_hess_y(lm, 2, isph)*adj_phi_y(lm, isph) &
                    & + fac*m_hess_z(lm, 2, isph)*adj_phi_z(lm, isph)
                force(3, isph) = force(3, isph) &
                    & + fac*m_hess_x(lm, 3, isph)*adj_phi_x(lm, isph) &
                    & + fac*m_hess_y(lm, 3, isph)*adj_phi_y(lm, isph) &
                    & + fac*m_hess_z(lm, 3, isph)*adj_phi_z(lm, isph)
            end do
        end do
    end do

    deallocate(zeta_comp, adj_phi_x, adj_phi_y, adj_phi_z, m_hess_x, &
        & m_hess_y, m_hess_z, m_grad, m_grad_component, stat=info)
    if (info.ne.0) then
        call update_error(ddx_error, "Allocation failed in grad_e")
        return
    end if

end subroutine grad_e

!> @ingroup Fortran_interface_multipolar
!> Given a multipolar distribution in real spherical harmonics and
!> centered on the spheres, compute the contributions to the forces
!> stemming from its electrostatic interactions.
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[inout] state: ddx state
!! @param[in] mmax: maximum angular momentum of the multipolar distribution
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!     size ((mmax+1)**2, nsph)
!! @param[inout] forces: forces array, size (3, nsph)
!! @param[inout] ddx_error: ddX error
!!
subroutine multipole_force_terms(params, constants, workspace, state, mmax, &
        & multipoles, forces, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax + 1)**2, params % nsph)
    real(dp), intent(inout) :: forces(3, params % nsph)

    call grad_phi(params, constants, workspace, state, mmax, multipoles, &
        & forces, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "multipole_force_terms: grad_phi returned an error, exiting")
        return
    end if

    if (params % model .eq. 3) then
        call grad_e(params, constants, workspace, state, mmax, multipoles, &
            & forces, ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "multipole_force_terms: grad_e returned an error, exiting")
            return
        end if
    end if
end subroutine multipole_force_terms

end module ddx_multipolar_solutes
