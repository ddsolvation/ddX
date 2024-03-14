!> Module for additional FMM functionality
module ddx_fmm

use ddx_constants
use fmmlib_interface
use mod_utils, only: ntot_sph_harm, fmm_error
use mod_harmonics, only: fmm_m2l

implicit none

contains

subroutine cart_propfar_lebedev(fmm_obj, params, constants, isph, &
        & do_v, v, do_e, e, do_g, g)
    implicit none

    type(fmm_type), intent(in):: fmm_obj
    type(ddx_params_type), intent(in):: params
    type(ddx_constants_type), intent(in):: constants
    integer, intent(in):: isph
    logical, intent(in):: do_v, do_e, do_g
    real(dp), intent(inout):: v(params % ngrid), e(3, params % ngrid), &
        & g(3, 3, params % ngrid)

    integer:: inode

    inode = fmm_obj % tree % particle_to_node(isph)

    call dgemv("T", (params % pl+1)**2, params % ngrid, one, &
        & constants % vgrid2, constants % vgrid_nbasis, &
        & fmm_obj % local_expansion(:, inode), 1, one, &
        & v, 1)

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

end module ddx_fmm
