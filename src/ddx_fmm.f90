!> Module for integrating the FMM library with additional operation
!! to deal with the Lebedev points and various ddX requirements.
module ddx_fmm

use ddx_constants
use ddx_core
use fmmlib_interface
use mod_utils, only: ntot_sph_harm, fmm_error
use mod_harmonics, only: fmm_m2l

implicit none

contains

!> Given multipolar distributions at the center of the spheres,
!! compute some electrostatic properties at the cavity (exposed Lebedev
!! points).
!!
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!                        size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[inout] ddx_error: ddX error
!! @param[in] do_v: logical to enable the computation of the electric potential
!! @param[in] do_e: logical to enable the computation of the electric field
!! @param[in] do_v: logical to enable the computation of the electric field
!!                  gradient
!! @param[in] do_diag: logical to enable the computation of properties
!!                     of a multipolar distribution to its own Lebedev points
!! @param[out,optional] v_cav: electric potential at the cavity point,
!!                             size (ncav)
!! @param[out,optional] e_cav: electric field at the cavity point,
!!                             size (3, ncav)
!! @param[out,optional] g_cav: electric field gradient at the cavity point,
!!                             size (3, 3, ncav)
!!
subroutine sphere_to_cav_cart(params, constants, workspace, multipoles, mmax, &
        & ddx_error, do_v, do_e, do_g, do_diag, v_cav, e_cav, g_cav)
    implicit none
    type(ddx_params_type), intent(in):: params
    type(ddx_constants_type), intent(in):: constants
    type(ddx_workspace_type), intent(inout):: workspace
    type(ddx_error_type), intent(inout) :: ddx_error
    logical, intent(in) :: do_v, do_e, do_g, do_diag
    integer, intent(in):: mmax
    real(dp), intent(in):: multipoles((mmax+1)**2, params % nsph)
    real(dp), optional, intent(out):: v_cav(constants % ncav), &
        & e_cav(3, constants % ncav), g_cav(3, 3, constants % ncav)
    integer :: isph, icav, igrid, v_stat, e_stat, g_stat
    real(dp), allocatable :: v_grid(:, :), e_grid(:, :, :), g_grid(:, :, :, :)

    if (do_v .and. (.not.present(v_cav))) then
        call update_error(ddx_error, &
            & "In `sphere_to_cav_cart`: do_v requested but v_cav missing")
        return
    end if
    if (do_e .and. (.not.present(e_cav))) then
        call update_error(ddx_error, &
            & "In `sphere_to_cav_cart`: do_e requested but e_cav missing")
        return
    end if
    if (do_g .and. (.not.present(g_cav))) then
        call update_error(ddx_error, &
            & "In `sphere_to_cav_cart`: do_g requested but g_cav missing")
        return
    end if

    if (do_v) allocate(v_grid(params % ngrid, params % nsph), stat=v_stat)
    if (do_e) allocate(e_grid(3, params % ngrid, params % nsph), stat=e_stat)
    if (do_g) allocate(g_grid(3, 3, params % ngrid, params % nsph), stat=g_stat)
    if ((v_stat.ne.0).or.(e_stat.ne.0).or.(g_stat.ne.0)) then
        call update_error(ddx_error, &
            & "Allocation failed in `sphere_to_cav_cart`")
        return
    end if

    ! compute the properties at the grid
    call sphere_to_grid_cart(params, constants, workspace, multipoles, mmax, &
        & ddx_error, do_v, do_e, do_g, do_diag, v_grid=v_grid, e_grid=e_grid, &
        & g_grid=g_grid)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, "`sphere_to_cav_cart`:" // &
            & " `sphere_to_grid_cart` returned an error, exiting")
        return
    end if

    ! discard the internal points
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if (constants % ui(igrid, isph) .gt. zero) then
                icav = icav + 1
                if (do_v) v_cav(icav) = v_grid(igrid, isph)
                if (do_e) e_cav(:, icav) = e_grid(:, igrid, isph)
                if (do_g) g_cav(:, :, icav) = g_grid(:, :, igrid, isph)
            end if
        end do
    end do

    if (allocated(v_grid)) deallocate(v_grid, stat=v_stat)
    if (allocated(e_grid)) deallocate(e_grid, stat=e_stat)
    if (allocated(g_grid)) deallocate(g_grid, stat=g_stat)
    if ((v_stat.ne.0).or.(e_stat.ne.0).or.(g_stat.ne.0)) then
        call update_error(ddx_error, &
            & "Deallocation failed in `sphere_to_cav_cart`")
    end if
end subroutine sphere_to_cav_cart

!> Given multipolar distributions at the center of the spheres,
!! compute some electrostatic properties at the grid (all the Lebedev
!! points).
!!
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[inout] workspace: ddx workspace
!! @param[in] multipoles: multipoles as real spherical harmonics,
!!                        size ((mmax+1)**2, nsph)
!! @param[in] mmax: maximum angular momentum of the multipoles
!! @param[inout] ddx_error: ddX error
!! @param[in] do_v: logical to enable the computation of the electric potential
!! @param[in] do_e: logical to enable the computation of the electric field
!! @param[in] do_v: logical to enable the computation of the electric field
!!                  gradient
!! @param[in] do_diag: logical to enable the computation of properties
!!                     of a multipolar distribution to its own Lebedev points
!! @param[out,optional] v_grid: electric potential at the grid points,
!!                              size (ncav, nsph)
!! @param[out,optional] e_grid: electric field at the grid points,
!!                              size (3, ncav, nsph)
!! @param[out,optional] g_grid: electric field gradient at the grid points,
!!                              size (3, 3, ncav, nsph)
!!
subroutine sphere_to_grid_cart(params, constants, workspace, multipoles, mmax, &
        & ddx_error, do_v, do_e, do_g, do_diag, v_grid, e_grid, g_grid)
    implicit none
    type(ddx_params_type), intent(in):: params
    type(ddx_constants_type), intent(in):: constants
    type(ddx_workspace_type), intent(inout):: workspace
    type(ddx_error_type), intent(inout) :: ddx_error
    logical, intent(in) :: do_v, do_e, do_g, do_diag
    integer, intent(in) :: mmax
    real(dp), intent(in) :: multipoles((mmax+1)**2, params % nsph)
    real(dp), optional, intent(out) :: v_grid(params % ngrid, params % nsph), &
        & e_grid(3, params % ngrid, params % nsph), &
        & g_grid(3, 3, params % ngrid, params % nsph)
    real(dp), allocatable :: v(:), e(:, :), g(:, :, :)
    integer :: isph, v_stat, e_stat, g_stat

    if (do_v .and. (.not.present(v_grid))) then
        call update_error(ddx_error, &
            & "In `sphere_to_grid_cart`: do_v requested but v_grid missing")
        return
    end if
    if (do_e .and. (.not.present(e_grid))) then
        call update_error(ddx_error, &
            & "In `sphere_to_grid_cart`: do_e requested but e_grid missing")
        return
    end if
    if (do_g .and. (.not.present(g_grid))) then
        call update_error(ddx_error, &
            & "In `sphere_to_grid_cart`: do_g requested but g_grid missing")
        return
    end if

    if (do_v) v_grid = zero
    if (do_e) e_grid = zero
    if (do_g) g_grid = zero

    if (do_v) allocate(v(params % ngrid), stat=v_stat)
    if (do_e) allocate(e(3, params % ngrid), stat=e_stat)
    if (do_g) allocate(g(3, 3, params % ngrid), stat=g_stat)
    if ((v_stat.ne.0).or.(e_stat.ne.0).or.(g_stat.ne.0)) then
        call update_error(ddx_error, &
            & "Allocation failed in `sphere_to_grid_cart`")
        return
    end if

    call tree_p2m(workspace % fmm_obj, multipoles, mmax)
    call tree_m2m(workspace % fmm_obj)
    call tree_m2l(workspace % fmm_obj)
    call tree_l2l(workspace % fmm_obj)

    do isph = 1, params % nsph
        if (do_v) v = zero
        if (do_e) e = zero
        if (do_g) g = zero
        call cart_propfar_lebedev(workspace % fmm_obj, params, constants, &
            & isph, ddx_error, do_v, do_e, do_g, v=v, e=e, g=g)
        call cart_propnear_lebedev(workspace % fmm_obj, params, constants, &
            & isph, ddx_error, do_v, do_e, do_g, do_diag, v=v, e=e, g=g)
        if (do_v) v_grid(:, isph) = v
        if (do_e) e_grid(:, :, isph) = e
        if (do_g) g_grid(:, :, :, isph) = g
    end do

    if (allocated(v)) deallocate(v, stat=v_stat)
    if (allocated(e)) deallocate(e, stat=e_stat)
    if (allocated(g)) deallocate(g, stat=g_stat)
    if ((v_stat.ne.0).or.(e_stat.ne.0).or.(g_stat.ne.0)) then
        call update_error(ddx_error, &
            & "Deallocation failed in `sphere_to_grid_cart`")
    end if

end subroutine sphere_to_grid_cart

!> Compute the electrostatic properties due to the far field at
!! the Lebedev points of sphere isph.
!!
!! @param[in] fmm_obj: FMM object
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[in] isph: index of the sphere
!! @param[inout] ddx_error: ddX error
!! @param[in] do_v: logical to enable the computation of the electric potential
!! @param[in] do_e: logical to enable the computation of the electric field
!! @param[in] do_v: logical to enable the computation of the electric field
!!                  gradient
!! @param[out,optional] v: electric potential at the sphere Lebedev points,
!!                         size(ngrid)
!! @param[out,optional] e: electric field at the sphere Lebedev points,
!!                         size(3, ngrid)
!! @param[out,optional] g: electric field gradient at the sphere Lebedev points,
!!                         size(3, 3, ngrid)
!!
subroutine cart_propfar_lebedev(fmm_obj, params, constants, isph, &
        & ddx_error, do_v, do_e, do_g, v, e, g)
    implicit none
    type(fmm_type), intent(in) :: fmm_obj
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in) :: isph
    logical, intent(in) :: do_v, do_e, do_g
    real(dp), optional, intent(inout) :: v(params % ngrid), &
        & e(3, params % ngrid), g(3, 3, params % ngrid)
    real(dp) :: r_t, local_expansion((params % pl+1)**2)
    integer :: inode, l, ind, igrid, i, stat_1, stat_2, stat_3
    real(dp), allocatable :: local_expansion_grad(:, :), &
        & local_expansion_hess(:, :), grid_hess_component(:, :)
    real(dp) :: radii(1)

    if (do_v .and. (.not.present(v))) then
        call update_error(ddx_error, &
            & "In `cart_propfar_lebedev`: do_v requested but v missing")
        return
    end if
    if (do_e .and. (.not.present(e))) then
        call update_error(ddx_error, &
            & "In `cart_propfar_lebedev`: do_e requested but e missing")
        return
    end if
    if (do_g .and. (.not.present(g))) then
        call update_error(ddx_error, &
            & "In `cart_propfar_lebedev`: do_g requested but g missing")
        return
    end if

    if (do_e.or.do_g) allocate(local_expansion_grad((params % pl + 1)**2, 3), &
        & stat=stat_1)
    if (do_g) allocate(local_expansion_hess((params % pl + 1)**2, 3), &
        & grid_hess_component(3, params % ngrid), stat=stat_2)
    if ((stat_1.ne.0).or.(stat_2.ne.0)) then
        call update_error(ddx_error, &
            & "Allocation failed in `cart_propfar_lebedev`")
        return
    end if

    inode = fmm_obj % tree % particle_to_node(isph)
    r_t = fmm_obj % tree % node_dimension(inode)
    local_expansion = fmm_obj % local_expansion(:, inode)

    ! the next transformations assume spherical harmonics scaled by
    ! their radius, so if we were working without, we need to put it
    ! back on
    if (.not.fmm_obj%radii_scaling) then
        do l = 0, params % pl
            ind = l*l + l + 1
            local_expansion(ind-l:ind+l) = local_expansion(ind-l:ind+l)*r_t**l
        end do
    end if

    ! far field potential (L2P)
    call dgemv("T", (params % pl+1)**2, params % ngrid, one, &
        & constants % vgrid2, constants % vgrid_nbasis, &
        & local_expansion, 1, one, v, 1)

    ! if needed assemble the gradient of the L2L
     if ((do_e.or.do_g) .and. params % pl.gt.0) then
        radii(1) = params % rsph(isph)
        call grad_l2l(local_expansion, params % pl, 1, local_expansion_grad, &
            & radii, ddx_error)
     end if

    ! far field electric field
    if (do_e .and. params % pl.gt.0) then
        ! L2P
        call dgemm("T", "N", 3, params % ngrid, (params % pl)**2, &
            & one, local_expansion_grad, (params % pl+1)**2, constants % vgrid2, &
            & constants % vgrid_nbasis, one, e, 3)
    end if

    ! far field electric field gradient
    if (do_g .and. params % pl.gt.1) then
        radii(1) = params % rsph(isph)
        ! component by component
        do i = 1, 3
            ! compute the hessian of the L2L (just one component)
            call grad_l2l(local_expansion_grad(:, i), params % pl, 1, &
                & local_expansion_hess, radii, ddx_error)
            ! L2P
            call dgemm("T", "N", 3, params % ngrid, (params % pl)**2, &
                & one, local_expansion_hess, (params % pl+1)**2, &
                & constants % vgrid2, constants % vgrid_nbasis, zero, &
                & grid_hess_component, 3)
            ! transpose the result
            do igrid = 1, params % ngrid
                g(:, i, igrid) = g(:, i, igrid) + grid_hess_component(:, igrid)
            end do
        end do
    end if

    ! deallocate gradient and hessian of the L2L
    if (allocated(local_expansion_grad)) deallocate(local_expansion_grad, &
        & stat=stat_1)
    if (allocated(local_expansion_hess)) deallocate(local_expansion_hess, &
        & stat=stat_2)
    if (allocated(grid_hess_component)) deallocate(grid_hess_component, &
        & stat=stat_3)
    if ((stat_1.ne.0).or.(stat_2.ne.0).or.(stat_3.ne.0)) then
        call update_error(ddx_error, &
            & "Deallocation failed in `cart_propfar_lebedev`")
    end if

end subroutine cart_propfar_lebedev

!> Compute the electrostatic properties due to the near field at
!! the Lebedev points of sphere isph.
!!
!! @param[in] fmm_obj: FMM object
!! @param[in] params: ddx parameters
!! @param[in] constants: ddx constants
!! @param[in] isph: index of the sphere
!! @param[inout] ddx_error: ddX error
!! @param[in] do_v: logical to enable the computation of the electric potential
!! @param[in] do_e: logical to enable the computation of the electric field
!! @param[in] do_v: logical to enable the computation of the electric field
!!                  gradient
!! @param[in] do_v: logical to enable the computation of the properties of
!!                  a given sphere at its own Lebedev points
!! @param[out,optional] v: electric potential at the sphere Lebedev points,
!!                         size(ngrid)
!! @param[out,optional] e: electric field at the sphere Lebedev points,
!!                         size(3, ngrid)
!! @param[out,optional] g: electric field gradient at the sphere Lebedev points,
!!                         size(3, 3, ngrid)
!!
subroutine cart_propnear_lebedev(fmm_obj, params, constants, isph, &
        & ddx_error, do_v, do_e, do_g, do_diag, v, e, g)
    implicit none
    type(fmm_type), intent(in):: fmm_obj
    type(ddx_params_type), intent(in):: params
    type(ddx_constants_type), intent(in):: constants
    type(ddx_error_type), intent(inout) :: ddx_error
    integer, intent(in):: isph
    logical, intent(in):: do_v, do_e, do_g, do_diag
    real(dp), optional, intent(inout):: v(params % ngrid), &
        & e(3, params % ngrid), g(3, 3, params % ngrid)
    real(dp), allocatable:: local_tmp(:)
    type(fmm_tree_type), pointer:: t
    integer :: igrid, i, i_node, j, j_node, n_targets, stat
    integer, allocatable :: targets(:)
    real(dp):: r_s, r_t, c_st(3), x2_y2, z2

    if (do_v .and. (.not.present(v))) then
        call update_error(ddx_error, &
            & "In `cart_propnear_lebedev`: do_v requested but v missing")
        return
    end if
    if (do_e .and. (.not.present(e))) then
        call update_error(ddx_error, &
            & "In `cart_propnear_lebedev`: do_e requested but e missing")
        return
    end if
    if (do_g .and. (.not.present(g))) then
        call update_error(ddx_error, &
            & "In `cart_propnear_lebedev`: do_g requested but g missing")
        return
    end if

    t => fmm_obj%tree

    i_node = t%particle_to_node(isph)
    if(t%particle_list%ri(i_node+1) - t%particle_list%ri(i_node) > 1) then
        call update_error(ddx_error, &
            & "FMM error: a node has more than one particle!")
        return
    end if

    n_targets = t%near_nl%ri(i_node+1) - t%near_nl%ri(i_node)

    allocate(local_tmp(ntot_sph_harm(fmm_obj%pmax_le)), &
        & targets(n_targets+1), stat=stat)
    if (stat.ne.0) then
        call update_error(ddx_error, &
            & "Allocation failed in `cart_propnear_lebedev`")
        return
    end if

    ! Define the targets as the neighbors, and if requested the
    ! diagonal case.
    i = 1
    do j = t%near_nl%ri(i_node), t%near_nl%ri(i_node+1)-1
        j_node = t%near_nl%ci(j)
        targets(i) = j_node
        i = i + 1
    end do
    if (do_diag) then
        targets(i) = i_node
        n_targets = n_targets + 1
    end if

    do j = 1, n_targets
        j_node = targets(j)

        do igrid = 1, params % ngrid
            r_s = t%node_dimension(j_node)
            r_t = t%node_dimension(i_node)
            c_st = t%node_centroid(:,j_node) - t%node_centroid(:,i_node) &
                & - r_t*constants%cgrid(:,igrid)

            ! after using the radii to find the relative displacement,
            ! if they are disabled in the transformations, we set them to 1.0
            if (.not.fmm_obj%radii_scaling) then
                r_s = 1.0
                r_t = 1.0
            end if

            call fmm_m2l(c_st, r_s, r_t, fmm_obj%pmax_mm, fmm_obj%pmax_le, &
                & fmm_obj%multipoles(:,j_node), local_tmp)
            if(do_v) then
                v(igrid) = v(igrid) + sqrt(4.0*pi)*local_tmp(1)
            end if

            if(do_e) then
                e(3,igrid) = e(3,igrid) - sqrt(4.0/3.0*pi) * local_tmp(3)/r_t
                e(1,igrid) = e(1,igrid) - sqrt(4.0/3.0*pi) * local_tmp(4)/r_t
                e(2,igrid) = e(2,igrid) - sqrt(4.0/3.0*pi) * local_tmp(2)/r_t
            end if

            if(do_g) then
                x2_y2 = sqrt(16.0*pi/15.0) * local_tmp(9) * 3.0
                z2 = sqrt(16.0*pi/5.0) * local_tmp(7)
                g(3,3,igrid) = g(3,3,igrid) + z2/(r_t**2)
                g(1,1,igrid) = g(1,1,igrid) + (x2_y2-z2) / 2.0/(r_t**2)
                g(2,2,igrid) = g(2,2,igrid) - (x2_y2+z2) / 2.0/(r_t**2)
                g(2,1,igrid) = g(2,1,igrid) &
                    & + 3.0*sqrt(4.0*pi/15.0)*local_tmp(5)/(r_t**2) ! xy
                g(1,2,igrid) = g(1,2,igrid) &
                    & + 3.0*sqrt(4.0*pi/15.0)*local_tmp(5)/(r_t**2) ! xy
                g(3,1,igrid) = g(3,1,igrid) &
                    & + 3.0*sqrt(4.0*pi/15.0)*local_tmp(8)/(r_t**2) ! xz
                g(1,3,igrid) = g(1,3,igrid) &
                    & + 3.0*sqrt(4.0*pi/15.0)*local_tmp(8)/(r_t**2) ! xz
                g(3,2,igrid) = g(3,2,igrid) &
                    & + 3.0*sqrt(4.0*pi/15.0)*local_tmp(6)/(r_t**2) ! yz
                g(2,3,igrid) = g(2,3,igrid) &
                    & + 3.0*sqrt(4.0*pi/15.0) * local_tmp(6)/(r_t**2) ! yz
            end if
        end do
    end do

    deallocate(local_tmp, targets, stat=stat)
    if (stat.ne.0) then
        call update_error(ddx_error, &
            & "Deallocation failed in `cart_propnear_lebedev`")
    end if

end subroutine cart_propnear_lebedev

!> Given a multipolar distribution compute the action of dL on it.
!!
!! @param[in] local_expansion: local multipolar expansions,
!!                             size ((pl+1)**2, nl)
!! @param[in] pl: maximum angular momentum of the multipolar distribution
!! @param[in] nl: number of local expansions
!! @param[out] l_grad: action of dL on the local expansions,
!!                     size ((pl+1)**2, 3, nl)
!! @param[in] radii: radii of the multipolar expansions, size(nl)
!! @param[inout] ddx_error: ddX error
!!
subroutine grad_l2l(local_expansion, pl, nl, l_grad, radii, ddx_error)
    implicit none
    integer, intent(in) :: pl, nl
    real(dp), intent(in) :: local_expansion((pl + 1)**2, nl)
    real(dp), intent(out) :: l_grad((pl + 1)**2, 3, nl)
    real(dp), intent(in) :: radii(nl)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! local variables
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    real(dp), allocatable :: tmp(:, :)
    integer :: info, i, l, indi, indj, m
    real(dp) :: tmp1, tmp2

    if (pl .le. 0) return

    allocate(tmp((pl + 1)**2, nl), stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, &
            & "Allocation failed in `grad_l2l`")
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

    do i = 1, nl
        l_grad(:, 3, i) = local_expansion(:, i)
        call fmm_sph_transform(pl, zx_coord_transform, one, &
            & local_expansion(:, i), zero, l_grad(:, 1, i))
        call fmm_sph_transform(pl, zy_coord_transform, one, &
            & local_expansion(:, i), zero, l_grad(:, 2, i))
    end do

    do l = 1, pl
        indi = l*l + l + 1
        indj = indi - 2*l
        tmp1 = -sqrt(dble(2*l-1)) / sqrt(dble(2*l+1))
        do m = 1-l, l-1
            tmp2 = sqrt(dble(l*l-m*m)) * tmp1
            l_grad(indj+m, :, :) = tmp2 * l_grad(indi+m, :, :)
        end do
    end do

    ! Scale by 1/rsph(isph) and rotate harmonics of degree up to pl-1 back to
    ! the initial axis. Coefficient of pl-th degree is zero so we ignore it.
    do i = 1, nl
        l_grad(1:pl**2, 3, i) = l_grad(1:pl**2, 3, i)/radii(i)
        tmp(1:pl**2, i) = l_grad(1:pl**2, 1, i)/radii(i)
        call fmm_sph_transform(pl-1, zx_coord_transform, one, &
            & tmp(1:pl**2, i), zero, l_grad(1:pl**2, 1, i))
        tmp(1:pl**2, i) = l_grad(1:pl**2, 2, i)/radii(i)
        call fmm_sph_transform(pl-1, zy_coord_transform, one, &
            & tmp(1:pl**2, i), zero, l_grad(1:pl**2, 2, i))
    end do

    ! Set degree pl to zero to avoid problems if user actually uses it
    l = pl
    indi = l*l + l + 1
    l_grad(indi-l:indi+l, :, :) = zero
    deallocate(tmp, stat=info)
    if (info .ne. 0) then
        call update_error(ddx_error, &
            & "Deallocation failed in `grad_l2l`")
        return
    end if

end subroutine grad_l2l

end module ddx_fmm
