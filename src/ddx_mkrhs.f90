!> Compute potential, its gradient and hessian in cavity points
!!
!! If ddx_data is FMM-ready then approximate output is computed by the FMM.
!!
!! @param[in] ddx_data: ddX object
!! @param[in] phi_flag: 1 if need to produce output potential
!! @param[out] phi_cav: Potential at cavity points. Referenced only if
!!      phi_flag=1
!! @param[in] grad_flag: 1 if need to produce gradient of potential
!! @param[out] gradphi_cav: Potential at cavity points. Referenced only if
!!      grad_flag=1
!! @param[in] hessian_flag: 1 if need to produce hessian of potential
!! @param[out] hessianphi_cav: Potential at cavity points. Referenced only if
!!      hessian_flag=1
subroutine mkrhs(ddx_data, phi_flag, phi_cav, grad_flag, gradphi_cav, &
        & hessian_flag, hessianphi_cav, psi)
    use ddx_core
    implicit none
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    integer, intent(in) :: phi_flag, grad_flag, hessian_flag
    ! Outputs
    real(dp), intent(out) :: phi_cav(ddx_data % constants % ncav)
    real(dp), intent(out) :: gradphi_cav(3, ddx_data % constants % ncav)
    real(dp), intent(out) :: hessianphi_cav(3, 3, ddx_data % constants % ncav)
    real(dp), intent(out) :: psi(ddx_data % constants % nbasis, &
        & ddx_data % params % nsph)
    ! Local variables
    integer :: isph, igrid, icav, inode, inear, jnear, jnode, jsph, i
    real(dp) :: d(3), v, tmpv, r, gradv(3), hessianv(3, 3), tmpd(3), epsp=one
    real(dp), allocatable :: grid_grad(:,:,:), grid_hessian(:,:,:,:), &
        & grid_hessian2(:,:,:)
    real(dp), external :: dnrm2
    real(dp) :: t

    if (grad_flag .eq. 1) allocate(grid_grad(ddx_data % params % ngrid, 3, &
        & ddx_data % params % nsph))
    if (hessian_flag .eq. 1) allocate(grid_hessian(ddx_data % params % ngrid, &
        & 3, 3, ddx_data % params % nsph), grid_hessian2(ddx_data % params % ngrid, &
        & 3, ddx_data % params % nsph))

    ! In case FMM is disabled compute phi and gradphi at cavity points by a
    ! naive double loop of a quadratic complexity
    if (ddx_data % params % fmm .eq. 0) then
        do icav = 1, ddx_data % constants % ncav
            v = zero
            gradv = zero
            hessianv = zero
            do isph = 1, ddx_data % params % nsph
                d = ddx_data % constants % ccav(:, icav) - &
                    & ddx_data % params % csph(:, isph)
                r = dnrm2(3, d, 1)
                d = d / r
                tmpv = ddx_data % params % charge(isph) / r
                v = v + tmpv
                tmpv = tmpv / r
                tmpd = tmpv * d
                tmpv = tmpv / r
                gradv = gradv - tmpd
                tmpd = three / r * tmpd
                hessianv(:, 1) = hessianv(:, 1) + tmpd*d(1)
                hessianv(:, 2) = hessianv(:, 2) + tmpd*d(2)
                hessianv(:, 3) = hessianv(:, 3) + tmpd*d(3)
                hessianv(1, 1) = hessianv(1, 1) - tmpv
                hessianv(2, 2) = hessianv(2, 2) - tmpv
                hessianv(3, 3) = hessianv(3, 3) - tmpv
            end do
            if (phi_flag .eq. 1) then
                phi_cav(icav) = v
            end if
            if (grad_flag .eq. 1) then
                gradphi_cav(:, icav) = gradv
            end if
            if (hessian_flag .eq. 1) then
                hessianphi_cav(:, :, icav) = hessianv
            end if
        end do
    ! Use the FMM otherwise
    else
        ! P2M step from centers of harmonics
        do isph = 1, ddx_data % params % nsph
            inode = ddx_data % constants % snode(isph)
            ddx_data % workspace % tmp_sph(1, isph) = &
                & ddx_data % params % charge(isph) &
                & / ddx_data % params % rsph(isph) / sqrt4pi
            ddx_data % workspace % tmp_sph(2:, isph) = zero
            ddx_data % workspace % tmp_node_m(1, inode) = &
                & ddx_data % workspace % tmp_sph(1, isph)
            ddx_data % workspace % tmp_node_m(2:, inode) = zero
        end do
        ! M2M, M2L and L2L translations

        call tree_m2m_rotation(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace % tmp_node_m)

        call tree_m2l_rotation(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace % tmp_node_m, ddx_data % workspace % tmp_node_l)

        call tree_l2l_rotation(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace % tmp_node_l)

        call tree_l2p(ddx_data % params, ddx_data % constants, one, &
            & ddx_data % workspace % tmp_node_l, zero, &
            & ddx_data % workspace % tmp_grid, ddx_data % workspace % tmp_sph_l)

        call tree_m2p(ddx_data % params, ddx_data % constants, &
            & ddx_data % params % lmax, one, ddx_data % workspace % tmp_sph, &
            & one, ddx_data % workspace % tmp_grid)

        ! Potential from each sphere to its own grid points
        call dgemm('T', 'N', ddx_data % params % ngrid, ddx_data % params % nsph, &
            & ddx_data % constants % nbasis, one, ddx_data % constants % vgrid2, &
            & ddx_data % constants % vgrid_nbasis, ddx_data % workspace % tmp_sph, &
            & ddx_data % constants % nbasis, one, ddx_data % workspace % tmp_grid, &
            & ddx_data % params % ngrid)

        ! Rearrange potential from all grid points to external only
        if (phi_flag.eq.1) then
            icav = 0
            do isph = 1, ddx_data % params % nsph
                do igrid = 1, ddx_data % params % ngrid
                    ! Do not count internal grid points
                    if(ddx_data % constants % ui(igrid, isph) .eq. zero) cycle
                    icav = icav + 1
                    phi_cav(icav) = ddx_data % workspace % tmp_grid(igrid, isph)
                end do
            end do
        end if

        ! Now compute near-field FMM gradients and hessians
        ! Cycle over all spheres
        if (grad_flag.eq.1 .or. hessian_flag.eq.1) then
            if (grad_flag.eq.1) grid_grad(:, :, :) = zero
            if (hessian_flag.eq.1) grid_hessian(:, :, :, :) = zero
            !$omp parallel do default(none) shared(ddx_data,grid_grad,grid_hessian, &
            !$omp grad_flag,hessian_flag) private(isph,igrid,inode,jnear,jnode, &
            !$omp jsph,d,r,tmpd,tmpv) schedule(dynamic)
            do isph = 1, ddx_data % params % nsph
                ! Cycle over all external grid points
                do igrid = 1, ddx_data % params % ngrid
                    if(ddx_data % constants % ui(igrid, isph) .eq. zero) cycle
                    ! Cycle over all near-field admissible pairs of spheres,
                    ! including pair (isph, isph) which is a self-interaction
                    inode = ddx_data % constants % snode(isph)
                    do jnear = ddx_data % constants % snear(inode), &
                        & ddx_data % constants % snear(inode+1) - 1
                        ! Near-field interactions are possible only between leaf
                        ! nodes, which must contain only a single input sphere
                        jnode = ddx_data % constants % near(jnear)
                        jsph = ddx_data % constants % order( &
                            & ddx_data % constants % cluster(1, jnode))
                        d = ddx_data % params % csph(:, isph) + &
                            & ddx_data % constants % cgrid(:, igrid) &
                            & *ddx_data % params % rsph(isph) &
                            & - ddx_data % params % csph(:, jsph)
!                       r = dnrm2(3, d, 1)
                        r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
                        d = d / r / r
                        tmpv = ddx_data % params % charge(jsph) / r
                        if (grad_flag.eq.1) grid_grad(igrid, :, isph) = &
                            & grid_grad(igrid, :, isph) - tmpv * d
                        tmpd = three * tmpv * d
                        tmpv = tmpv / r / r
                        if (hessian_flag.eq.1) then
                            grid_hessian(igrid, 1, :, isph) = &
                                & grid_hessian(igrid, 1, :, isph) + d(1)*tmpd
                            grid_hessian(igrid, 2, :, isph) = &
                                & grid_hessian(igrid, 2, :, isph) + d(2)*tmpd
                            grid_hessian(igrid, 3, :, isph) = &
                                & grid_hessian(igrid, 3, :, isph) + d(3)*tmpd
                            grid_hessian(igrid, 1, 1, isph) = &
                                & grid_hessian(igrid, 1, 1, isph) - tmpv
                            grid_hessian(igrid, 2, 2, isph) = &
                                & grid_hessian(igrid, 2, 2, isph) - tmpv
                            grid_hessian(igrid, 3, 3, isph) = &
                                & grid_hessian(igrid, 3, 3, isph) - tmpv
                        end if
                    end do
                end do
            end do
        end if

        ! Take into account far-field FMM gradients only if pl > 0
        if ((ddx_data % params % pl .gt. 0) .and. (grad_flag.eq.1)) then
            ! Get gradient of L2L
            call tree_grad_l2l(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace % tmp_node_l, &
                & ddx_data % workspace % tmp_sph_l_grad, &
                & ddx_data % workspace % tmp_sph_l)
            ! Apply L2P for every axis with -1 multiplier since grad over
            ! target point is equal to grad over source point
            call dgemm('T', 'N', ddx_data % params % ngrid, &
                & 3*ddx_data % params % nsph, (ddx_data % params % pl)**2, &
                & -one, ddx_data % constants % vgrid2, &
                & ddx_data % constants % vgrid_nbasis, &
                & ddx_data % workspace % tmp_sph_l_grad, &
                & (ddx_data % params % pl+1)**2, one, grid_grad, &
                & ddx_data % params % ngrid)
        end if
        ! Take into account far-field FMM hessians only if pl > 1
        if ((ddx_data % params % pl .gt. 1) .and. (hessian_flag.eq.1)) then
            do i = 1, 3
                ! Load previously computed gradient into leaves, since
                ! tree_grad_l2l currently takes local expansions of entire
                ! tree. In future it might be changed.
                do isph = 1, ddx_data % params % nsph
                    inode = ddx_data % constants % snode(isph)
                    ddx_data % workspace % tmp_node_l(:, inode) = ddx_data % &
                        & workspace % tmp_sph_l_grad(:, i, isph)
                end do
                ! Get gradient of a gradient of L2L. Currently this uses input
                ! pl maximal degree of local harmonics but in reality we need
                ! only pl-1 maximal degree since coefficients of harmonics of a
                ! degree pl are zeros.
                call tree_grad_l2l(ddx_data % params, ddx_data % constants, &
                    & ddx_data % workspace % tmp_node_l, &
                    & ddx_data % workspace % tmp_sph_l_grad2, &
                    & ddx_data % workspace % tmp_sph_l)
                ! Apply L2P for every axis
                call dgemm('T', 'N', ddx_data % params % ngrid, &
                    & 3*ddx_data % params % nsph, (ddx_data % params % pl-1)**2, &
                    & one, ddx_data % constants % vgrid2, &
                    & ddx_data % constants % vgrid_nbasis, &
                    & ddx_data % workspace % tmp_sph_l_grad2, &
                    & (ddx_data % params % pl+1)**2, zero, grid_hessian2, &
                    & ddx_data % params % ngrid)
                ! Properly copy hessian
                grid_hessian(:, i, :, :) = grid_hessian(:, i, :, :) + grid_hessian2
            end do
        end if

        ! Copy output for external grid points only
        icav = 0
        do isph = 1, ddx_data % params % nsph
            do igrid = 1, ddx_data % params % ngrid
                if(ddx_data % constants % ui(igrid, isph) .eq. zero) cycle
                icav = icav + 1
                if (grad_flag .eq. 1) gradphi_cav(:, icav) = &
                    & grid_grad(igrid, :, isph)
                if (hessian_flag .eq. 1) hessianphi_cav(:, :, icav) = &
                    & grid_hessian(igrid, :, :, isph)
            end do
        end do
    end if

    ! Vector psi
    psi(2:, :) = zero
    do isph = 1, ddx_data % params % nsph
        psi(1, isph) = sqrt4pi * ddx_data % params % charge(isph)
    end do

    ! deallocate temporary
    if (grad_flag .eq. 1) deallocate(grid_grad)
    if (hessian_flag .eq. 1) deallocate(grid_hessian, grid_hessian2)
end subroutine mkrhs
