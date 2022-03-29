!> Compute potential, its gradient and hessian in cavity points
!!
!! If ddx_data is FMM-ready then approximate output is computed by the FMM.
!!
!! @param[in] params: parameters data structure
!! @param[in] constants: constants data structure
!! @param[in] worskpace: temporary arrays data structure
!! @param[in] phi_flag: 1 if need to produce output potential
!! @param[out] phi_cav: Potential at cavity points. Referenced only if
!!      phi_flag=1
!! @param[in] grad_flag: 1 if need to produce gradient of potential
!! @param[out] gradphi_cav: Potential at cavity points. Referenced only if
!!      grad_flag=1
!! @param[in] hessian_flag: 1 if need to produce hessian of potential
!! @param[out] hessianphi_cav: Potential at cavity points. Referenced only if
!!      hessian_flag=1
subroutine mkrhs(params, constants, workspace, phi_flag, phi_cav, grad_flag, &
        & gradphi_cav, hessian_flag, hessianphi_cav, psi)
    use ddx_core
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: phi_flag, grad_flag, hessian_flag
    real(dp), intent(out) :: phi_cav(constants % ncav)
    real(dp), intent(out) :: gradphi_cav(3, constants % ncav)
    real(dp), intent(out) :: hessianphi_cav(3, 3, constants % ncav)
    real(dp), intent(out) :: psi(constants % nbasis, &
        & params % nsph)
    ! Local variables
    integer :: isph, igrid, icav, inode, inear, jnear, jnode, jsph, i
    real(dp) :: d(3), v, tmpv, r, gradv(3), hessianv(3, 3), tmpd(3), epsp=one
    real(dp), allocatable :: grid_grad(:,:,:), grid_hessian(:,:,:,:), &
        & grid_hessian2(:,:,:)
    real(dp), external :: dnrm2
    real(dp) :: t

    if (grad_flag .eq. 1) allocate(grid_grad(params % ngrid, 3, &
        & params % nsph))
    if (hessian_flag .eq. 1) allocate(grid_hessian(params % ngrid, &
        & 3, 3, params % nsph), grid_hessian2(params % ngrid, &
        & 3, params % nsph))

    ! In case FMM is disabled compute phi and gradphi at cavity points by a
    ! naive double loop of a quadratic complexity
    if (params % fmm .eq. 0) then
        do icav = 1, constants % ncav
            v = zero
            gradv = zero
            hessianv = zero
            do isph = 1, params % nsph
                d = constants % ccav(:, icav) - &
                    & params % csph(:, isph)
                r = dnrm2(3, d, 1)
                d = d / r
                tmpv = params % charge(isph) / r
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
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_sph(1, isph) = &
                & params % charge(isph) &
                & / params % rsph(isph) / sqrt4pi
            workspace % tmp_sph(2:, isph) = zero
            workspace % tmp_node_m(1, inode) = &
                & workspace % tmp_sph(1, isph)
            workspace % tmp_node_m(2:, inode) = zero
        end do
        ! M2M, M2L and L2L translations

        call tree_m2m_rotation(params, constants, &
            & workspace % tmp_node_m)

        call tree_m2l_rotation(params, constants, &
            & workspace % tmp_node_m, workspace % tmp_node_l)

        call tree_l2l_rotation(params, constants, &
            & workspace % tmp_node_l)

        call tree_l2p(params, constants, one, &
            & workspace % tmp_node_l, zero, &
            & workspace % tmp_grid, workspace % tmp_sph_l)

        call tree_m2p(params, constants, &
            & params % lmax, one, workspace % tmp_sph, &
            & one, workspace % tmp_grid)

        ! Potential from each sphere to its own grid points
        call dgemm('T', 'N', params % ngrid, params % nsph, &
            & constants % nbasis, one, constants % vgrid2, &
            & constants % vgrid_nbasis, workspace % tmp_sph, &
            & constants % nbasis, one, workspace % tmp_grid, &
            & params % ngrid)

        ! Rearrange potential from all grid points to external only
        if (phi_flag.eq.1) then
            icav = 0
            do isph = 1, params % nsph
                do igrid = 1, params % ngrid
                    ! Do not count internal grid points
                    if(constants % ui(igrid, isph) .eq. zero) cycle
                    icav = icav + 1
                    phi_cav(icav) = workspace % tmp_grid(igrid, isph)
                end do
            end do
        end if

        ! Now compute near-field FMM gradients and hessians
        ! Cycle over all spheres
        if (grad_flag.eq.1 .or. hessian_flag.eq.1) then
            if (grad_flag.eq.1) grid_grad(:, :, :) = zero
            if (hessian_flag.eq.1) grid_hessian(:, :, :, :) = zero
            !$omp parallel do default(none) shared(params,constants,workspace, &
            !$omp grid_grad,grid_hessian,grad_flag,hessian_flag) &
            !$omp private(isph,igrid,inode,jnear,jnode,jsph,d,r,tmpd,tmpv) &
            !$omp schedule(dynamic)
            do isph = 1, params % nsph
                ! Cycle over all external grid points
                do igrid = 1, params % ngrid
                    if(constants % ui(igrid, isph) .eq. zero) cycle
                    ! Cycle over all near-field admissible pairs of spheres,
                    ! including pair (isph, isph) which is a self-interaction
                    inode = constants % snode(isph)
                    do jnear = constants % snear(inode), &
                        & constants % snear(inode+1) - 1
                        ! Near-field interactions are possible only between leaf
                        ! nodes, which must contain only a single input sphere
                        jnode = constants % near(jnear)
                        jsph = constants % order( &
                            & constants % cluster(1, jnode))
                        d = params % csph(:, isph) + &
                            & constants % cgrid(:, igrid) &
                            & *params % rsph(isph) &
                            & - params % csph(:, jsph)
!                       r = dnrm2(3, d, 1)
                        r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
                        d = d / r / r
                        tmpv = params % charge(jsph) / r
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
        if ((params % pl .gt. 0) .and. (grad_flag.eq.1)) then
            ! Get gradient of L2L
            call tree_grad_l2l(params, constants, &
                & workspace % tmp_node_l, &
                & workspace % tmp_sph_l_grad, &
                & workspace % tmp_sph_l)
            ! Apply L2P for every axis with -1 multiplier since grad over
            ! target point is equal to grad over source point
            call dgemm('T', 'N', params % ngrid, &
                & 3*params % nsph, (params % pl)**2, &
                & -one, constants % vgrid2, &
                & constants % vgrid_nbasis, &
                & workspace % tmp_sph_l_grad, &
                & (params % pl+1)**2, one, grid_grad, &
                & params % ngrid)
        end if
        ! Take into account far-field FMM hessians only if pl > 1
        if ((params % pl .gt. 1) .and. (hessian_flag.eq.1)) then
            do i = 1, 3
                ! Load previously computed gradient into leaves, since
                ! tree_grad_l2l currently takes local expansions of entire
                ! tree. In future it might be changed.
                do isph = 1, params % nsph
                    inode = constants % snode(isph)
                    workspace % tmp_node_l(:, inode) = &
                        & workspace % tmp_sph_l_grad(:, i, isph)
                end do
                ! Get gradient of a gradient of L2L. Currently this uses input
                ! pl maximal degree of local harmonics but in reality we need
                ! only pl-1 maximal degree since coefficients of harmonics of a
                ! degree pl are zeros.
                call tree_grad_l2l(params, constants, &
                    & workspace % tmp_node_l, &
                    & workspace % tmp_sph_l_grad2, &
                    & workspace % tmp_sph_l)
                ! Apply L2P for every axis
                call dgemm('T', 'N', params % ngrid, &
                    & 3*params % nsph, (params % pl-1)**2, &
                    & one, constants % vgrid2, &
                    & constants % vgrid_nbasis, &
                    & workspace % tmp_sph_l_grad2, &
                    & (params % pl+1)**2, zero, grid_hessian2, &
                    & params % ngrid)
                ! Properly copy hessian
                grid_hessian(:, i, :, :) = grid_hessian(:, i, :, :) + grid_hessian2
            end do
        end if

        ! Copy output for external grid points only
        icav = 0
        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
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
    do isph = 1, params % nsph
        psi(1, isph) = sqrt4pi * params % charge(isph)
    end do

    ! deallocate temporary
    if (grad_flag .eq. 1) deallocate(grid_grad)
    if (hessian_flag .eq. 1) deallocate(grid_hessian, grid_hessian2)
end subroutine mkrhs
