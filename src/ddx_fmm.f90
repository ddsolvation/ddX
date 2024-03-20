!> Module for additional FMM functionality
module ddx_fmm

use ddx_constants
use ddx_core
use fmmlib_interface
use mod_utils, only: ntot_sph_harm, fmm_error
use mod_harmonics, only: fmm_m2l

implicit none

contains

subroutine cart_propfar_lebedev(fmm_obj, params, constants, isph, &
        & do_v, v, do_e, e, do_g, g)
    implicit none

    type(fmm_type), intent(in) :: fmm_obj
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: isph
    logical, intent(in) :: do_v, do_e, do_g
    real(dp), intent(inout) :: v(params % ngrid), e(3, params % ngrid), &
        & g(3, 3, params % ngrid)

    real(dp) :: r_t, local_expansion((params % pl+1)**2)
    integer :: inode, l, ind
    real(dp), allocatable :: local_expansion_grad(:, :)
    real(dp) :: radii(1)

    if (do_e) then
        allocate(local_expansion_grad((params % pl + 1)**2, 3))
    end if

    inode = fmm_obj % tree % particle_to_node(isph)
    r_t = fmm_obj % tree % node_dimension(inode)
    local_expansion = fmm_obj % local_expansion(:, inode)

    if (.not.fmm_obj%radii_scaling) then
        do l = 0, params % pl
            ind = l*l + l + 1
            local_expansion(ind-l:ind+l) = local_expansion(ind-l:ind+l)*r_t**l
        end do
    end if

    call dgemv("T", (params % pl+1)**2, params % ngrid, one, &
        & constants % vgrid2, constants % vgrid_nbasis, &
        & local_expansion, 1, one, v, 1)

    if (do_e) then
        radii(1) = params % rsph(isph)
        call grad_l2l(local_expansion, params % pl, 1, local_expansion_grad, &
            radii)
        call dgemm("T", "N", 3, params % ngrid, (params % pl)**2, &
            & one, local_expansion_grad, (params % pl+1)**2, constants % vgrid2, &
            & constants % vgrid_nbasis, one, e, 3)
    end if

    if (allocated(local_expansion_grad)) deallocate(local_expansion_grad)

end subroutine

subroutine cart_propnear_lebedev(fmm_obj, params, constants, isph, &
        & do_v, v, do_e, e, do_g, g, do_diag)
    implicit none

    type(fmm_type), intent(in):: fmm_obj
    type(ddx_params_type), intent(in):: params
    type(ddx_constants_type), intent(in):: constants
    integer, intent(in):: isph
    logical, intent(in):: do_v, do_e, do_g, do_diag
    real(dp), intent(inout):: v(params % ngrid), e(3, params % ngrid), &
        & g(3, 3, params % ngrid)

    integer:: i_part, igrid

    real(dp), allocatable:: local_tmp(:)
    type(fmm_tree_type), pointer:: t
    integer:: i_node, j, j_node, j_particle
    real(dp):: r_s, r_t, c_st(3), x2_y2, z2, x2z_y2z, z3, xz2, yz2, x3_3xy2, y3_3x2y

    t => fmm_obj%tree

    i_node = t%particle_to_node(isph)
    if(t%particle_list%ri(i_node+1) - t%particle_list%ri(i_node) > 1) then
        call fmm_error("A node has more than one particle!")
    end if

    allocate(local_tmp(ntot_sph_harm(fmm_obj%pmax_le)))

    do j = t%near_nl%ri(i_node), t%near_nl%ri(i_node+1)-1
        j_node = t%near_nl%ci(j)
        j_particle = t%particle_list%ci(t%particle_list%ri(j_node))

        do igrid = 1, params % ngrid
            r_s = t%node_dimension(j_node)
            r_t = t%node_dimension(i_node)
            c_st = t%node_centroid(:,j_node) - t%node_centroid(:,i_node) &
                & - r_t*constants%cgrid(:,igrid)
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

            !if(do_grdE) then
            !    x2_y2 = sqrt(16.0*pi/15.0) * local_tmp(9) * 3.0
            !    z2 = sqrt(16.0*pi/5.0) * local_tmp(7)
            !    grdE(6) = grdE(6) + z2  ! zz
            !    grdE(1) = grdE(1) + (x2_y2-z2) / 2.0
            !    grdE(3) = grdE(3) - (x2_y2+z2) / 2.0
            !    grdE(2) = grdE(2) + 3.0*sqrt(4.0*pi/15.0) * local_tmp(5)  ! xy
            !    grdE(4) = grdE(4) + 3.0*sqrt(4.0*pi/15.0) * local_tmp(8)  ! xz
            !    grdE(5) = grdE(5) + 3.0*sqrt(4.0*pi/15.0) * local_tmp(6)  ! yz
            !end if
        end do
    end do

    ! add the diagonal
    if (do_diag) then
        do igrid = 1, params % ngrid
            r_s = t%node_dimension(i_node)
            r_t = t%node_dimension(i_node)
            c_st = - r_t*constants%cgrid(:,igrid)
            if (.not.fmm_obj%radii_scaling) then
                r_s = 1.0
                r_t = 1.0
            end if
            call fmm_m2l(c_st, r_s, r_t, fmm_obj%pmax_mm, fmm_obj%pmax_le, &
                         fmm_obj%multipoles(:,i_node), local_tmp)
            if(do_v) then
                v(igrid) = v(igrid) + sqrt(4.0*pi)*local_tmp(1)
            end if

            if(do_e) then
                e(3,igrid) = e(3,igrid) - sqrt(4.0/3.0*pi) * local_tmp(3)/r_t
                e(1,igrid) = e(1,igrid) - sqrt(4.0/3.0*pi) * local_tmp(4)/r_t
                e(2,igrid) = e(2,igrid) - sqrt(4.0/3.0*pi) * local_tmp(2)/r_t
            end if

            !if(do_grdE) then
            !    x2_y2 = sqrt(16.0*pi/15.0) * local_tmp(9) * 3.0
            !    z2 = sqrt(16.0*pi/5.0) * local_tmp(7)
            !    grdE(6) = grdE(6) + z2  ! zz
            !    grdE(1) = grdE(1) + (x2_y2-z2) / 2.0
            !    grdE(3) = grdE(3) - (x2_y2+z2) / 2.0
            !    grdE(2) = grdE(2) + 3.0*sqrt(4.0*pi/15.0) * local_tmp(5)  ! xy
            !    grdE(4) = grdE(4) + 3.0*sqrt(4.0*pi/15.0) * local_tmp(8)  ! xz
            !    grdE(5) = grdE(5) + 3.0*sqrt(4.0*pi/15.0) * local_tmp(6)  ! yz
            !end if
        end do
    end if

end subroutine

!> Given a multipolar distribution compute the action of dP on it, this
!> is required in the computation of the forces.
!! @param[in] local_expansion: local_expansion as real spherical harmonics,
!!     size ((pl+1)**2, nsph)
!! @param[in] pl: maximum angular momentum of the multipolar distribution
!! @param[in] nl: number of local_expansion
!! @param[out] l_grad: gradient of the M2M operator,
!!     size ((pl + 1)**2, 3, nl)
!!
subroutine grad_l2l(local_expansion, pl, nl, l_grad, radii)
    implicit none
    integer, intent(in) :: pl, nl
    real(dp), intent(in) :: local_expansion((pl + 1)**2, nl)
    real(dp), intent(out) :: l_grad((pl + 1)**2, 3, nl)
    real(dp), intent(in) :: radii(nl)
    ! local variables
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    real(dp), allocatable :: tmp(:, :)
    integer :: info, i, l, indi, indj, m
    real(dp) :: tmp1, tmp2

    if (pl .le. 0) return

    allocate(tmp((pl + 1)**2, nl), stat=info)
    if (info .ne. 0) then
        stop "Allocation failed"
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
        stop "Allocation failed"
    end if

end subroutine grad_l2l

end module ddx_fmm
