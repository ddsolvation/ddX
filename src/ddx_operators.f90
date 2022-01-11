!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_operators.f90
!! Operators of ddX software
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

!> Operators shared among ddX methods
module ddx_operators
! Use underlying core routines
use ddx_core
implicit none

contains

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
    real(dp), allocatable :: grid_grad(:,:,:), grid_hessian(:,:,:,:), grid_hessian2(:,:,:)
    real(dp), external :: dnrm2
    real(dp) :: t
    if (ddx_data % params % force .eq. 1) then
        allocate(grid_grad(ddx_data % params % ngrid, 3, &
        & ddx_data % params % nsph), grid_hessian(ddx_data % params % ngrid, &
        & 3, 3, ddx_data % params % nsph), &
        & grid_hessian2(ddx_data % params % ngrid, 3, ddx_data % params % nsph))
    end if
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
        t = omp_get_wtime()
        do isph = 1, ddx_data % params % nsph
            inode = ddx_data % constants % snode(isph)
            ddx_data % workspace % tmp_sph(1, isph) = ddx_data % params % charge(isph) / &
                & ddx_data % params % rsph(isph) / sqrt4pi
            ddx_data % workspace % tmp_sph(2:, isph) = zero
            ddx_data % workspace % tmp_node_m(1, inode) = ddx_data % workspace % tmp_sph(1, isph)
            ddx_data % workspace % tmp_node_m(2:, inode) = zero
        end do
        write(6,*) 'centers', omp_get_wtime() - t
        ! M2M, M2L and L2L translations

        t = omp_get_wtime()
        call tree_m2m_rotation(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace % tmp_node_m)
        write(6,*) 'm2m rotation', omp_get_wtime() - t

        t = omp_get_wtime()
        call tree_m2l_rotation(ddx_data % params, ddx_data % constants, & 
            & ddx_data % workspace % tmp_node_m, ddx_data % workspace % tmp_node_l)
        write(6,*) 'm2l rotation', omp_get_wtime() - t

        t = omp_get_wtime()
        call tree_l2l_rotation(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace % tmp_node_l)
        write(6,*) 'l2l rotation', omp_get_wtime() - t

        t = omp_get_wtime()
        call tree_l2p(ddx_data % params, ddx_data % constants, one, &
            & ddx_data % workspace % tmp_node_l, zero, ddx_data % workspace % tmp_grid)
        write(6,*) 'l2p', omp_get_wtime() - t

        t = omp_get_wtime()
        call tree_m2p(ddx_data % params, ddx_data % constants, & 
            & ddx_data % params % lmax, one, ddx_data % workspace % tmp_sph, &
            & one, ddx_data % workspace % tmp_grid)
        write(6,*) 'm2p', omp_get_wtime() - t

        ! Potential from each sphere to its own grid points
        t = omp_get_wtime()
        call dgemm('T', 'N', ddx_data % params % ngrid, ddx_data % params % nsph, &
            & ddx_data % constants % nbasis, one, ddx_data % constants % vgrid2, &
            & ddx_data % constants % vgrid_nbasis, ddx_data % workspace % tmp_sph, ddx_data % constants % nbasis, &
            & one, ddx_data % workspace % tmp_grid, ddx_data % params % ngrid)
        write(6,*) 'dgemm', omp_get_wtime() - t
        ! Rearrange potential from all grid points to external only
        icav = 0
        do isph = 1, ddx_data % params % nsph
            do igrid = 1, ddx_data % params % ngrid
                ! Do not count internal grid points
                if(ddx_data % constants % ui(igrid, isph) .eq. zero) cycle
                icav = icav + 1
                phi_cav(icav) = ddx_data % workspace % tmp_grid(igrid, isph)
            end do
        end do
        if (ddx_data % params % force .eq. 1) then
            ! Now compute near-field FMM gradients and hessians
            ! Cycle over all spheres
            !$omp parallel do default(none) shared(ddx_data,grid_grad,grid_hessian) &
            !$omp private(isph,igrid,inode,jnear,jnode,jsph,d,r,tmpd,tmpv) schedule(static,1)
            do isph = 1, ddx_data % params % nsph
                grid_grad(:, :, isph) = zero
                grid_hessian(:, :, :, isph) = zero
                ! Cycle over all external grid points
                do igrid = 1, ddx_data % params % ngrid
                    if(ddx_data % constants % ui(igrid, isph) .eq. zero) cycle
                    ! Cycle over all near-field admissible pairs of spheres,
                    ! including pair (isph, isph) which is a self-interaction
                    inode = ddx_data % constants % snode(isph)
                    do jnear = ddx_data % constants % snear(inode), ddx_data % constants % snear(inode+1)-1
                        ! Near-field interactions are possible only between leaf
                        ! nodes, which must contain only a single input sphere
                        jnode = ddx_data % constants % near(jnear)
                        jsph = ddx_data % constants % order(ddx_data % constants % cluster(1, jnode))
                        d = ddx_data % params % csph(:, isph) + &
                            & ddx_data % constants % cgrid(:, igrid)*ddx_data % params % rsph(isph) - &
                            & ddx_data % params % csph(:, jsph)
                        r = dnrm2(3, d, 1)
                        d = d / r / r
                        tmpv = ddx_data % params % charge(jsph) / r
                        grid_grad(igrid, :, isph) = grid_grad(igrid, :, isph) - &
                            & tmpv * d
                        tmpd = three * tmpv * d
                        tmpv = tmpv / r / r
                        grid_hessian(igrid, 1, :, isph) = grid_hessian(igrid, 1, :, isph) + &
                            & d(1)*tmpd
                        grid_hessian(igrid, 2, :, isph) = grid_hessian(igrid, 2, :, isph) + &
                            & d(2)*tmpd
                        grid_hessian(igrid, 3, :, isph) = grid_hessian(igrid, 3, :, isph) + &
                            & d(3)*tmpd
                        grid_hessian(igrid, 1, 1, isph) = grid_hessian(igrid, 1, 1, isph) - tmpv
                        grid_hessian(igrid, 2, 2, isph) = grid_hessian(igrid, 2, 2, isph) - tmpv
                        grid_hessian(igrid, 3, 3, isph) = grid_hessian(igrid, 3, 3, isph) - tmpv
                    end do
                end do
            end do
            ! Take into account far-field FMM gradients only if pl > 0
            if (ddx_data % params % pl .gt. 0) then
                ! Get gradient of L2L
                call tree_grad_l2l(ddx_data % params, ddx_data % constants, ddx_data % workspace % tmp_node_l, &
                    & ddx_data % workspace % tmp_sph_l_grad, ddx_data % workspace % tmp_sph_l)
                ! Apply L2P for every axis with -1 multiplier since grad over
                ! target point is equal to grad over source point
                call dgemm('T', 'N', ddx_data % params % ngrid, 3*ddx_data % params % nsph, &
                    & (ddx_data % params % pl)**2, -one, ddx_data % constants % vgrid2, &
                    & ddx_data % constants % vgrid_nbasis, ddx_data % workspace % tmp_sph_l_grad, &
                    & (ddx_data % params % pl+1)**2, one, grid_grad, &
                    & ddx_data % params % ngrid)
            end if
            ! Take into account far-field FMM hessians only if pl > 1
            if (ddx_data % params % pl .gt. 1) then
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
                    ! degree pl zre zeros.
                    call tree_grad_l2l(ddx_data % params, ddx_data % constants, ddx_data % workspace % tmp_node_l, &
                        & ddx_data % workspace % tmp_sph_l_grad2, ddx_data % workspace % tmp_sph_l)
                    ! Apply L2P for every axis
                    call dgemm('T', 'N', ddx_data % params % ngrid, 3*ddx_data % params % nsph, &
                        & (ddx_data % params % pl-1)**2, one, ddx_data % constants % vgrid2, &
                        & ddx_data % constants % vgrid_nbasis, ddx_data % workspace % tmp_sph_l_grad2, &
                        & (ddx_data % params % pl+1)**2, zero, grid_hessian2, &
                        & ddx_data % params % ngrid)
                    ! Properly copy hessian
                    grid_hessian(:, i, :, :) = grid_hessian(:, i, :, :) + grid_hessian2
                end do
            end if
        end if
        ! Copy output for external grid points only
        icav = 0
        do isph = 1, ddx_data % params % nsph
            do igrid = 1, ddx_data % params % ngrid
                if(ddx_data % constants % ui(igrid, isph) .eq. zero) cycle
                icav = icav + 1
                if (ddx_data % params % force .eq. 1) then
                    gradphi_cav(:, icav) = grid_grad(igrid, :, isph)
                    hessianphi_cav(:, :, icav) = grid_hessian(igrid, :, :, isph)
                end if
            end do
        end do
    end if
    ! Vector psi
    psi(2:, :) = zero
    do isph = 1, ddx_data % params % nsph
        psi(1, isph) = sqrt4pi * ddx_data % params % charge(isph)
    end do
    if (ddx_data % params % force .eq. 1) then
        deallocate(grid_grad,grid_hessian,grid_hessian2)
    end if
end subroutine mkrhs

!> Single layer operator matvec without diagonal blocks
subroutine lx_nodiag(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: isph, l, ind
    !! Local arrays
    real(dp) :: pot(params % ngrid), vplm(constants % nbasis), &
        & vylm(constants % nbasis), vcos(params % lmax + 1), &
        & vsin(params % lmax + 1)
    !! Initialize
    y = zero
    !$omp parallel do default(none) shared(params,constants,x,y) &
    !$omp private(isph,pot,vylm,vplm,vcos,vsin) schedule(static,1)
    do isph = 1, params % nsph
        ! Compute NEGATIVE action of off-digonal blocks
        call calcv(params, constants, .false., isph, pot, &
            & x, vylm, vplm, vcos, vsin)
        call intrhs(1, constants % nbasis, params % ngrid, &
            & constants % vwgrid, constants % vgrid_nbasis, &
            & pot, y(:, isph))
        ! Action of off-diagonal blocks
        y(:, isph) = -y(:, isph)
    end do
end subroutine lx_nodiag

!> Single layer operator matvec with diagonal blocks
subroutine lx(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: isph, l, ind
    !! Initialize
    y = zero
    do isph = 1, params % nsph
        ! Compute NEGATIVE action of off-digonal blocks
        call calcv(params, constants, .false., isph, workspace % tmp_grid, &
            & x, workspace % tmp_vylm, workspace % tmp_vplm, &
            & workspace % tmp_vcos, workspace % tmp_vsin)
        call intrhs(1, constants % nbasis, params % ngrid, &
            & constants % vwgrid, constants % vgrid_nbasis, &
            & workspace % tmp_grid, y(:, isph))
    end do
    ! Loop over harmonics
    do l = 0, params % lmax
        ind = l*l + l + 1
        y(ind-l:ind+l, :) = -y(ind-l:ind+l, :) + &
            & x(ind-l:ind+l, :) / (constants % vscales(ind)**2)
    end do
end subroutine lx

!> Adjoint single layer operator matvec without diagonal blocks
subroutine lstarx_nodiag(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: isph, igrid, l, ind
    !! Initilize
    y = zero
    !! Expand x over spherical harmonics
    ! Loop over spheres      
    do isph = 1, params % nsph
        ! Loop over grid points
        do igrid = 1, params % ngrid
            workspace % tmp_grid(igrid, isph) = dot_product(x(:, isph), &
                & constants % vgrid(:constants % nbasis, igrid))
        end do
    end do
    !! Compute action
    ! Loop over spheres
    do isph = 1, params % nsph
        ! Compute NEGATIVE action of off-digonal blocks
        call adjrhs(params, constants, isph, workspace % tmp_grid, &
            & y(:, isph), workspace % tmp_vylm, workspace % tmp_vplm, &
            & workspace % tmp_vcos, workspace % tmp_vsin)
        y(:, isph) = - y(:, isph)
    end do
end subroutine lstarx_nodiag

!> Adjoint single layer operator matvec with diagonal blocks
subroutine lstarx(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: isph, igrid, l, ind
    !! Initilize
    y = zero
    !! Expand x over spherical harmonics
    ! Loop over spheres      
    do isph = 1, params % nsph
        ! Loop over grid points
        do igrid = 1, params % ngrid
            workspace % tmp_grid(igrid, isph) = dot_product(x(:, isph), &
                & constants % vgrid(:constants % nbasis, igrid))
        end do
    end do
    !! Compute action
    ! Loop over spheres
    do isph = 1, params % nsph
        ! Compute NEGATIVE action of off-digonal blocks
        call adjrhs(params, constants, isph, workspace % tmp_grid, &
            & y(:, isph), workspace % tmp_vylm, workspace % tmp_vplm, &
            & workspace % tmp_vcos, workspace % tmp_vsin)
    end do
    ! Loop over harmonics
    do l = 0, params % lmax
        ind = l*l + l + 1
        y(ind-l:ind+l, :) = -y(ind-l:ind+l, :) + &
            & x(ind-l:ind+l, :) / (constants % vscales(ind)**2)
    end do
end subroutine lstarx

!> Diagonal preconditioning for Lx operator
!!
!! Applies inverse diagonal (block) of the L matrix
subroutine ldm1x(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: l, ind
    !! Loop over harmonics
    do l = 0, params % lmax
        ind = l*l + l + 1
        y(ind-l:ind+l, :) = x(ind-l:ind+l, :) * (constants % vscales(ind)**2)
    end do
end subroutine ldm1x

!> Double layer operator matvec without diagonal blocks
subroutine dx_nodiag(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local parameter
    integer, parameter :: do_diag=0
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dx_dense(params, constants, workspace, do_diag, x, y)
    else
        call dx_fmm(params, constants, workspace, do_diag, x, y)
    end if
end subroutine dx_nodiag

!> Double layer operator matvec with diagonal blocks
subroutine dx(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local parameter
    integer, parameter :: do_diag=1
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dx_dense(params, constants, workspace, do_diag, x, y)
    else
        call dx_fmm(params, constants, workspace, do_diag, x, y)
    end if
end subroutine dx

!> Baseline implementation of double layer operator
subroutine dx_dense(params, constants, workspace, do_diag, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    real(dp) :: c(3), vij(3), sij(3)
    real(dp) :: vvij, tij, tt, f, f1, rho, ctheta, stheta, cphi, sphi
    integer :: its, isph, jsph, l, m, ind, lm, istatus
    real(dp), external :: dnrm2
    y = zero
    do isph = 1, params % nsph
        ! compute the "potential" from the other spheres
        ! at the exposed lebedv points of the i-th sphere 
        workspace % tmp_grid(:, 1) = zero
        do its = 1, params % ngrid
            if (constants % ui(its,isph).gt.zero) then
                c = params % csph(:,isph) + params % rsph(isph)* &
                    & constants % cgrid(:,its)
                do jsph = 1, params % nsph
                    if (jsph.ne.isph) then
                        ! build the geometrical variables
                        vij = c - params % csph(:,jsph)
                        !vvij = sqrt(dot_product(vij,vij))
                        vvij = dnrm2(3, vij, 1)
                        tij = vvij / params % rsph(jsph)
                        sij = vij/vvij 
                        ! build the local basis
                        call ylmbas(sij, rho, ctheta, stheta, cphi, sphi, &
                            & params % lmax, constants % vscales, &
                            & workspace % tmp_vylm, workspace % tmp_vplm, &
                            & workspace % tmp_vcos, workspace % tmp_vsin)
                        ! with all the required stuff, finally compute
                        ! the "potential" at the point 
                        tt = one/tij 
                        do l = 0, params % lmax
                            ind = l*l + l + 1
                            f = fourpi*dble(l)/(two*dble(l) + one)*tt
                            do m = -l, l
                                workspace % tmp_grid(its, 1) = &
                                    & workspace % tmp_grid(its, 1) + &
                                    & f*x(ind + m,jsph) * &
                                    & workspace % tmp_vylm(ind + m, 1)
                            end do
                            tt = tt/tij
                        end do
                    else if (do_diag .eq. 1) then
                        ! add the diagonal contribution
                        do l = 0, params % lmax
                            ind = l*l + l + 1
                            f = (two*dble(l) + one)/fourpi
                            do m = -l, l
                                workspace % tmp_grid(its, 1) = &
                                    & workspace % tmp_grid(its, 1) - &
                                    & pt5*x(ind + m,isph) * &
                                    & constants % vgrid(ind + m,its)/f
                            end do
                        end do
                    end if 
                end do
                workspace % tmp_grid(its, 1) = constants % ui(its, isph) * &
                    & workspace % tmp_grid(its, 1)
            end if
        end do
        ! now integrate the potential to get its modal representation
        call intrhs(1, constants % nbasis, params % ngrid, &
            & constants % vwgrid, constants % vgrid_nbasis, &
            & workspace % tmp_grid, y(:,isph))
    end do
end subroutine dx_dense

!> FMM-accelerated implementation of double layer operator
subroutine dx_fmm(params, constants, workspace, do_diag, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: isph, inode, l, indl, indl1, m
    real(dp) :: finish_time, start_time
    !! Scale input harmonics at first
    workspace % tmp_sph(1, :) = zero
    indl = 2
    do l = 1, params % lmax
        indl1 = (l+1)**2
        workspace % tmp_sph(indl:indl1, :) = l * x(indl:indl1, :)
        indl = indl1 + 1
    end do
    ! Load input harmonics into tree data
    if(params % lmax .lt. params % pm) then
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_m(1:constants % nbasis, inode) = &
                & workspace % tmp_sph(:, isph)
            workspace % tmp_node_m(constants % nbasis+1:, inode) = zero
        end do
    else
        indl = (params % pm+1)**2
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_m(:, inode) = workspace % tmp_sph(:indl, isph)
        end do
    end if
    ! Do FMM operations
    call tree_m2m_rotation(params, constants, workspace % tmp_node_m)
    call tree_m2l_rotation(params, constants, workspace % tmp_node_m, &
        & workspace % tmp_node_l)
    call tree_l2l_rotation(params, constants, workspace % tmp_node_l)
    call tree_l2p(params, constants, one, workspace % tmp_node_l, zero, &
        & workspace % tmp_grid)
    call tree_m2p(params, constants, params % lmax, one, workspace % tmp_sph, one, &
        & workspace % tmp_grid)
    ! Apply diagonal contribution if needed
    if(do_diag .eq. 1) then
        call dgemm('T', 'N', params % ngrid, params % nsph, &
            & constants % nbasis, -pt5, constants % vgrid2, &
            & constants % vgrid_nbasis, x, constants % nbasis, one, &
            & workspace % tmp_grid, params % ngrid)
    end if
    ! Multiply by characteristic function
    workspace % tmp_grid = workspace % tmp_grid * constants % ui
    ! now integrate the potential to get its modal representation
    ! output y is overwritten here
    call dgemm('N', 'N', constants % nbasis, params % nsph, params % ngrid, &
        & one, constants % vwgrid, constants % vgrid_nbasis, workspace % tmp_grid, &
        & params % ngrid, zero, y, constants % nbasis)
end subroutine dx_fmm

!> Adjoint double layer operator matvec without diagonal blocks
subroutine dstarx_nodiag(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer, parameter :: do_diag=0
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dstarx_dense(params, constants, workspace, do_diag, x, y)
    else
        call dstarx_fmm(params, constants, workspace, do_diag, x, y)
    end if
end subroutine dstarx_nodiag

!> Adjoint double layer operator matvec with diagonal blocks
subroutine dstarx(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer, parameter :: do_diag=1
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dstarx_dense(params, constants, workspace, do_diag, x, y)
    else
        call dstarx_fmm(params, constants, workspace, do_diag, x, y)
    end if
end subroutine dstarx

!> Baseline implementation of adjoint double layer operator
subroutine dstarx_dense(params, constants, workspace, do_diag, x, y)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    real(dp) :: c(3), vji(3), sji(3)
    real(dp) :: vvji, tji, fourpi, tt, f, f1, rho, ctheta, stheta, cphi, sphi
    integer :: its, isph, jsph, l, m, ind, lm, istatus
    real(dp), external :: dnrm2
    y = zero
    ! this loop is easily parallelizable
    do isph = 1, params % nsph
        do jsph = 1, params % nsph
            if (jsph.ne.isph) then
                do its = 1, params % ngrid
                    if (constants % ui(its,jsph).gt.zero) then
                        ! build the geometrical variables
                        vji = params % csph(:,jsph) + params % rsph(jsph) * &
                            & constants % cgrid(:,its) - params % csph(:,isph)
                        !vvji = sqrt(dot_product(vji,vji))
                        vvji = dnrm2(3, vji, 1)
                        tji = vvji/params % rsph(isph)
                        sji = vji/vvji
                        ! build the local basis
                        call ylmbas(sji, rho, ctheta, stheta, cphi, sphi, &
                            & params % lmax, constants % vscales, &
                            & workspace % tmp_vylm, workspace % tmp_vplm, &
                            & workspace % tmp_vcos, workspace % tmp_vsin)
                        tt = constants % ui(its,jsph)*dot_product(constants % vwgrid(:,its),x(:,jsph))/tji
                        do l = 0, params % lmax
                            ind = l*l + l + 1
                            f = dble(l)*tt/ constants % vscales(ind)**2
                            do m = -l, l
                                y(ind+m,isph) = y(ind+m,isph) + &
                                    & f*workspace % tmp_vylm(ind+m, 1)
                            end do
                            tt = tt/tji
                        end do
                    end if
                end do
            else if (do_diag .eq. 1) then
                do its = 1, params % ngrid
                    f = pt5*constants % ui(its,jsph)*dot_product(constants % vwgrid(:,its),x(:,jsph))
                    do l = 0, params % lmax
                        ind = l*l + l + 1
                        y(ind-l:ind+l,isph) = y(ind-l:ind+l,isph) - &
                            & f*constants % vgrid(ind-l:ind+l,its)/ &
                            & constants % vscales(ind)**2
                    end do
                end do
            end if
        end do
    end do
end subroutine dstarx_dense

!> FMM-accelerated implementation of adjoint double layer operator
subroutine dstarx_fmm(params, constants, workspace, do_diag, x, y)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    integer :: isph, inode, l, indl, indl1, m
    real(dp) :: finish_time, start_time
    ! Adjoint integration
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & one, constants % vwgrid, constants % vgrid_nbasis, x, &
        & constants % nbasis, zero, workspace % tmp_grid, params % ngrid)
    ! Multiply by characteristic function
    workspace % tmp_grid = workspace % tmp_grid * constants % ui
    ! Do FMM operations adjointly
    call tree_m2p_adj(params, constants, params % lmax, one, workspace % tmp_grid, &
        & zero, y)
    call tree_l2p_adj(params, constants, one, workspace % tmp_grid, zero, &
        & workspace % tmp_node_l)
    call tree_l2l_rotation_adj(params, constants, workspace % tmp_node_l)
    call tree_m2l_rotation_adj(params, constants, workspace % tmp_node_l, &
        & workspace % tmp_node_m)
    call tree_m2m_rotation_adj(params, constants, workspace % tmp_node_m)
    ! Adjointly move tree multipole harmonics into output
    if(params % lmax .lt. params % pm) then
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            y(:, isph) = y(:, isph) + &
                & workspace % tmp_node_m(1:constants % nbasis, inode)
        end do
    else
        indl = (params % pm+1)**2
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            y(1:indl, isph) = y(1:indl, isph) + workspace % tmp_node_m(:, inode)
        end do
    end if
    ! Scale output harmonics at last
    y(1, :) = zero
    indl = 2
    do l = 1, params % lmax
        indl1 = (l+1)**2
        y(indl:indl1, :) = l * y(indl:indl1, :)
        indl = indl1 + 1
    end do
    ! Apply diagonal contribution if needed
    if(do_diag .eq. 1) then
        call dgemm('N', 'N', constants % nbasis, params % nsph, &
            & params % ngrid, -pt5, constants % vgrid2, &
            & constants % vgrid_nbasis, workspace % tmp_grid, params % ngrid, &
            & one, y, constants % nbasis)
    end if
end subroutine dstarx_fmm

!> Apply \f$ R \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R x = - D x \f$ with excluded diagonal influence (blocks
!! D_ii are assumed to be zero).
subroutine rx(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    real(dp) :: fac
    !! Output `y` is cleaned here
    call dx_nodiag(params, constants, workspace, x, y)
    y = -y
end subroutine rx

!> Apply \f$ R^* \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R^* x = - D^* x \f$ with excluded diagonal influence (blocks
!! D_ii are assumed to be zero).
subroutine rstarx(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    real(dp) :: fac
    !! Output `y` is cleaned here
    call dstarx_nodiag(params, constants, workspace, x, y)
    y = -y
end subroutine rstarx

!> Apply \f$ R_\varepsilon \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\varepsilon x = (2\pi(\varepsilon + 1) / (\varepsilon
!! - 1) - D) x \f$.
subroutine repsx(params, constants, workspace, x, y)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    real(dp) :: fac
    !! Output `y` is cleaned here
    call dx(params, constants, workspace, x, y)
    ! Apply diagonal
    fac = twopi * (params % eps + one) / (params % eps - one)
    y = fac*x - y
end subroutine repsx

!> Apply \f$ R_\varepsilon^* \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R^*_\varepsilon x = (2\pi(\varepsilon + 1) / (\varepsilon
!! - 1) - D) x \f$.
subroutine rstarepsx(params, constants, workspace, x, y)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dstarx(params, constants, workspace, x, y)
    ! Apply diagonal
    fac = twopi * (params % eps + one) / (params % eps - one)
    y = fac*x - y
end subroutine rstarepsx

!> Apply \f$ R_\infty \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\infty x = (2\pi - D) x \f$.
!!
!! @param[in] ddx_data:
!! @param[in] x:
!! @param[out] y:
subroutine rinfx(params, constants, workspace, x, y)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dx(params, constants, workspace, x, y)
    ! Apply diagonal
    y = twopi*x - y
end subroutine rinfx

!> Apply \f$ R_\infty^* \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\infty^* x = (2\pi - D^*) x \f$.
!!
!! @param[in] ddx_data:
!! @param[in] x:
!! @param[out] y:
subroutine rstarinfx(params, constants, workspace, x, y)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dstarx(params, constants, workspace, x, y)
    ! Apply diagonal
    y = twopi*x - y
end subroutine rstarinfx

!> Apply preconditioner for 
subroutine apply_repsx_prec(params, constants, workspace, x, y)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    integer :: isph
    ! simply do a matrix-vector product with the transposed preconditioner 
    ! !!$omp parallel do default(shared) schedule(dynamic) &
    ! !!$omp private(isph)
    do isph = 1, params % nsph
        call dgemm('N', 'N', constants % nbasis, 1, constants % nbasis, one, &
            & constants % rx_prc(:, :, isph), constants % nbasis, x(:, isph), &
            & constants % nbasis, zero, y(:, isph), constants % nbasis)
    end do
    !call prtsph("rx_prec x", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, x)
    !call prtsph("rx_prec y", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, y)
end subroutine apply_repsx_prec

!> Apply preconditioner for 
subroutine apply_rstarepsx_prec(params, constants, workspace, x, y)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    integer :: isph
    ! simply do a matrix-vector product with the transposed preconditioner 
    ! !!$omp parallel do default(shared) schedule(dynamic) &
    ! !!$omp private(isph)
    do isph = 1, params % nsph
        call dgemm('T', 'N', constants % nbasis, 1, constants % nbasis, one, &
            & constants % rx_prc(:, :, isph), constants % nbasis, x(:, isph), &
            & constants % nbasis, zero, y(:, isph), constants % nbasis)
    end do
    !call prtsph("rx_prec x", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, x)
    !call prtsph("rx_prec y", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, y)
end subroutine apply_rstarepsx_prec

subroutine gradr(params, constants, workspace, g, ygrid, fx)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: g(constants % nbasis, params % nsph), &
        & ygrid(params % ngrid, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: fx(3, params % nsph)
    ! Check which gradr to execute
    if (params % fmm .eq. 1) then
        call gradr_fmm(params, constants, workspace, g, ygrid, fx)
    else
        call gradr_dense(params, constants, workspace, g, ygrid, fx)
    end if
end subroutine gradr

subroutine gradr_dense(params, constants, workspace, g, ygrid, fx)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: g(constants % nbasis, params % nsph), &
        & ygrid(params % ngrid, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: fx(3, params % nsph)
    ! Local variables
    integer :: isph
    real(dp) :: vplm(constants % nbasis), vcos(params % lmax+1), &
        & vsin(params % lmax+1), basloc(constants % nbasis), &
        & dbsloc(3, constants % nbasis)
    ! Simply cycle over all spheres
    do isph = 1, params % nsph
        call gradr_sph(params, constants, isph, vplm, vcos, vsin, basloc, &
            & dbsloc, g, ygrid, fx(:, isph))
    end do
end subroutine gradr_dense

subroutine gradr_sph(params, constants, isph, vplm, vcos, vsin, basloc, &
        & dbsloc, g, ygrid, fx)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: isph
    real(dp), intent(in) :: g(constants % nbasis, params % nsph), &
        & ygrid(params % ngrid, params % nsph)
    real(dp), intent(inout) :: vplm(constants % nbasis), vcos(params % lmax+1), &
        & vsin(params % lmax+1), basloc(constants % nbasis), &
        & dbsloc(3, constants % nbasis), fx(3)
    ! various scratch arrays
    real(dp) vik(3), sik(3), vki(3), ski(3), vkj(3), skj(3), vji(3), &
        & sji(3), va(3), vb(3), a(3), d(3)
    ! jacobian matrix
    real(dp) sjac(3,3)
    ! indexes
    integer its, ik, ksph, l, m, ind, jsph, icomp, jcomp
    ! various scalar quantities
    real(dp) cx, cy, cz, vvki, tki, gg, fl, fac, vvkj, tkj
    real(dp) tt, fcl, dij, fjj, gi, fii, vvji, tji, qji
    real(dp) b, vvik, tik, qik, tlow, thigh, duj
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    real(dp), external :: dnrm2

    tlow  = one - pt5*(one - params % se)*params % eta
    thigh = one + pt5*(one + params % se)*params % eta

    ! first set of contributions:
    ! diagonal block, kc and part of kb

    fx = zero
    do its = 1, params % ngrid
        ! sum over ksph in neighbors of isph
        do ik = constants % inl(isph), constants % inl(isph+1) - 1
            ksph = constants % nl(ik)
            ! build geometrical quantities
            cx = params % csph(1,ksph) + params % rsph(ksph)*constants % cgrid(1,its)
            cy = params % csph(2,ksph) + params % rsph(ksph)*constants % cgrid(2,its)
            cz = params % csph(3,ksph) + params % rsph(ksph)*constants % cgrid(3,its)
            vki(1) = cx - params % csph(1,isph)
            vki(2) = cy - params % csph(2,isph)
            vki(3) = cz - params % csph(3,isph)
            !vvki = sqrt(vki(1)*vki(1) + vki(2)*vki(2) + &
            !    & vki(3)*vki(3))
            vvki = dnrm2(3, vki, 1)
            tki  = vvki/params % rsph(isph)

            ! contributions involving grad i of uk come from the switching
            ! region.
            ! note: ui avoids contributions from points that are in the
            ! switching between isph and ksph but are buried in a third
            ! sphere.
            if ((tki.gt.tlow).and.(tki.lt.thigh) .and. &
                & constants % ui(its,ksph).gt.zero) then
                ! other geometrical quantities
                ski = vki/vvki

                ! diagonal block kk contribution, with k in n(i)
                gg = zero
                do l = 0, params % lmax
                    ind = l*l + l + 1
                    fl = dble(l)
                    fac = twopi/(two*fl + one)
                    do m = -l, l
                        !! DEBUG comment
                        gg = gg + fac*constants % vgrid(ind+m,its)*g(ind+m,ksph)
                    end do
                end do

                ! kc contribution
                do jsph = 1, params % nsph
                    if (jsph.ne.ksph .and. jsph.ne.isph) then 
                        vkj(1) = cx - params % csph(1,jsph)
                        vkj(2) = cy - params % csph(2,jsph)
                        vkj(3) = cz - params % csph(3,jsph)
                        vvkj = sqrt(vkj(1)*vkj(1) + vkj(2)*vkj(2) + &
                            & vkj(3)*vkj(3))
                        vvkj = dnrm2(3, vkj, 1)
                        tkj  = vvkj/params % rsph(jsph)
                        skj  = vkj/vvkj
                        call ylmbas(skj, rho, ctheta, stheta, cphi, sphi, &
                            & params % lmax, constants % vscales, basloc, &
                            & vplm, vcos, vsin)
                        tt = one/tkj
                        do l = 0, params % lmax
                            ind = l*l + l + 1
                            fcl = - fourpi*dble(l)/(two*dble(l)+one)*tt
                            do m = -l, l
                                !! DEBUG comment
                                gg = gg + fcl*g(ind+m,jsph)*basloc(ind+m)
                            end do
                            tt = tt/tkj
                        end do
                        !call fmm_m2p(vkj, params % rsph(jsph), &
                        !    & params % lmax, constants % vscales_rel, -one, &
                        !    & g(:, jsph), one, gg)
                    end if
                end do

                ! part of kb contribution
                call ylmbas(ski, rho, ctheta, stheta, cphi, sphi, &
                    & params % lmax, constants % vscales, basloc, &
                    & vplm, vcos, vsin)
                tt = one/tki
                do l = 0, params % lmax
                    ind = l*l + l + 1
                    fcl = - four*pi*dble(l)/(two*dble(l)+one)*tt
                    do m = -l, l
                        !! DEBUG comment
                        gg = gg + fcl*g(ind+m,isph)*basloc(ind+m)
                    end do
                    tt = tt/tki
                end do
                !call fmm_m2p(vki, params % rsph(isph), &
                !    & params % lmax, constants % vscales_rel, -one, &
                !    & g(:, isph), one, gg)

                ! common step, product with grad i uj
                duj = dfsw(tki,params % se, params % eta)/params % rsph(isph)
                fjj = duj*constants % wgrid(its)*gg*ygrid(its,ksph)
                fx(1) = fx(1) - fjj*ski(1)
                fx(2) = fx(2) - fjj*ski(2)
                fx(3) = fx(3) - fjj*ski(3)
            end if
        end do

        ! diagonal block ii contribution
        if (constants % ui(its,isph).gt.zero.and.constants % ui(its,isph).lt.one) then
            gi = zero
            do l = 0, params % lmax
                ind = l*l + l + 1
                fl = dble(l)
                fac = twopi/(two*fl + one)
                do m = -l, l 
                    !! DEBUG comment
                    gi = gi + fac*constants % vgrid(ind+m,its)*g(ind+m,isph)
                    !gi = gi + pt5*constants % vgrid2(ind+m,its)*g(ind+m,isph)
                end do
            end do
            !do l = 0, (params % lmax+1)**2
            !    gi = gi + constants % vgrid2(l, its)*g(l, isph)
            !end do
            !gi = pt5 * gi
            fii = constants % wgrid(its)*gi*ygrid(its,isph)
            fx(1) = fx(1) + fii*constants % zi(1,its,isph)
            fx(2) = fx(2) + fii*constants % zi(2,its,isph)
            fx(3) = fx(3) + fii*constants % zi(3,its,isph)
        end if
    end do

    ! second set of contributions:
    ! part of kb and ka
    do its = 1, params % ngrid

        ! run over all the spheres except isph 
        do jsph = 1, params % nsph
            if (constants % ui(its,jsph).gt.zero .and. jsph.ne.isph) then
                ! build geometrical quantities
                cx = params % csph(1,jsph) + params % rsph(jsph)*constants % cgrid(1,its)
                cy = params % csph(2,jsph) + params % rsph(jsph)*constants % cgrid(2,its)
                cz = params % csph(3,jsph) + params % rsph(jsph)*constants % cgrid(3,its)
                vji(1) = cx - params % csph(1,isph)
                vji(2) = cy - params % csph(2,isph)
                vji(3) = cz - params % csph(3,isph)
                !vvji = sqrt(vji(1)*vji(1) + vji(2)*vji(2) + &
                !    &  vji(3)*vji(3))
                vvji = dnrm2(3, vji, 1)
                tji = vvji/params % rsph(isph)
                qji = one/vvji
                sji = vji/vvji

                ! build the jacobian of sji
                sjac = zero
                sjac(1,1) = - one
                sjac(2,2) = - one
                sjac(3,3) = - one
                do icomp = 1, 3
                    do jcomp = 1, 3
                        sjac(icomp,jcomp) = qji*(sjac(icomp,jcomp) &
                            & + sji(icomp)*sji(jcomp))
                    end do
                end do

                ! assemble the local basis and its gradient
                !call dbasis(sji,basloc,dbsloc,vplm,vcos,vsin)
                call dbasis(params, constants, sji,basloc,dbsloc,vplm,vcos,vsin)

                ! assemble the contribution
                a = zero
                tt = one/(tji)
                do l = 0, params % lmax
                    ind = l*l + l + 1
                    fl = dble(l)
                    fcl = - tt*fourpi*fl/(two*fl + one)
                    do m = -l, l
                        fac = fcl*g(ind+m,isph)
                        b = (fl + one)*basloc(ind+m)/(params % rsph(isph)*tji)

                        ! apply the jacobian to grad y
                        va(1) = sjac(1,1)*dbsloc(1,ind+m) + &
                            & sjac(1,2)*dbsloc(2,ind+m) + sjac(1,3)*dbsloc(3,ind+m)
                        va(2) = sjac(2,1)*dbsloc(1,ind+m) + &
                            & sjac(2,2)*dbsloc(2,ind+m) + sjac(2,3)*dbsloc(3,ind+m)
                        va(3) = sjac(3,1)*dbsloc(1,ind+m) + &
                            & sjac(3,2)*dbsloc(2,ind+m) + sjac(3,3)*dbsloc(3,ind+m)
                        a(1) = a(1) + fac*(sji(1)*b + va(1))
                        a(2) = a(2) + fac*(sji(2)*b + va(2))
                        a(3) = a(3) + fac*(sji(3)*b + va(3))
                    end do
                    tt = tt/tji
                end do
                fac = constants % ui(its,jsph)*constants % wgrid(its)*ygrid(its,jsph)
                fx(1) = fx(1) - fac*a(1)
                fx(2) = fx(2) - fac*a(2)
                fx(3) = fx(3) - fac*a(3)
            end if
        end do
    end do

    ! ka contribution
    do its = 1, params % ngrid
        cx = params % csph(1,isph) + params % rsph(isph)*constants % cgrid(1,its)
        cy = params % csph(2,isph) + params % rsph(isph)*constants % cgrid(2,its)
        cz = params % csph(3,isph) + params % rsph(isph)*constants % cgrid(3,its)
        a = zero

        ! iterate on all the spheres except isph
        do ksph = 1, params % nsph
            if (constants % ui(its,isph).gt.zero .and. ksph.ne.isph) then
                ! geometrical stuff
                vik(1) = cx - params % csph(1,ksph)
                vik(2) = cy - params % csph(2,ksph)
                vik(3) = cz - params % csph(3,ksph)
                !vvik = sqrt(vik(1)*vik(1) + vik(2)*vik(2) + & 
                !    & vik(3)*vik(3))
                vvik = dnrm2(3, vik, 1)
                tik = vvik/params % rsph(ksph)
                qik = one/vvik
                sik = vik/vvik

                ! build the jacobian of sik
                sjac = zero
                sjac(1,1) = one
                sjac(2,2) = one
                sjac(3,3) = one
                do icomp = 1, 3
                    do jcomp = 1, 3
                    sjac(icomp,jcomp) = qik*(sjac(icomp,jcomp) &
                        & - sik(icomp)*sik(jcomp))
                    end do
                end do

                ! if we are in the switching region, recover grad_i u_i
                vb = zero
                if (constants % ui(its,isph).lt.one) then
                    vb(1) = constants % zi(1,its,isph)
                    vb(2) = constants % zi(2,its,isph)
                    vb(3) = constants % zi(3,its,isph)
                end if

                ! assemble the local basis and its gradient
                !call dbasis(sik,basloc,dbsloc,vplm,vcos,vsin)
                call dbasis(params, constants, sik,basloc,dbsloc,vplm,vcos,vsin)

                ! assemble the contribution
                tt = one/(tik)
                do l = 0, params % lmax
                    ind = l*l + l + 1
                    fl = dble(l)
                    fcl = - tt*fourpi*fl/(two*fl + one)
                        do m = -l, l
                        fac = fcl*g(ind+m,ksph)
                        fac = - fac*basloc(ind+m)
                        a(1) = a(1) + fac*vb(1)
                        a(2) = a(2) + fac*vb(2) 
                        a(3) = a(3) + fac*vb(3)

                        fac = constants % ui(its,isph)*fcl*g(ind+m,ksph)
                        b = - (fl + one)*basloc(ind+m)/(params % rsph(ksph)*tik)

                        ! apply the jacobian to grad y
                        va(1) = sjac(1,1)*dbsloc(1,ind+m) + &
                            & sjac(1,2)*dbsloc(2,ind+m) + sjac(1,3)*dbsloc(3,ind+m)
                        va(2) = sjac(2,1)*dbsloc(1,ind+m) + &
                            & sjac(2,2)*dbsloc(2,ind+m) + sjac(2,3)*dbsloc(3,ind+m)
                        va(3) = sjac(3,1)*dbsloc(1,ind+m) + &
                            & sjac(3,2)*dbsloc(2,ind+m) + sjac(3,3)*dbsloc(3,ind+m)
                        a(1) = a(1) + fac*(sik(1)*b + va(1))
                        a(2) = a(2) + fac*(sik(2)*b + va(2))
                        a(3) = a(3) + fac*(sik(3)*b + va(3))
                    end do
                    tt = tt/tik
                end do
            end if
        end do
        fac = constants % wgrid(its)*ygrid(its,isph)
        fx(1) = fx(1) - fac*a(1)
        fx(2) = fx(2) - fac*a(2)
        fx(3) = fx(3) - fac*a(3)
    end do
end subroutine gradr_sph

! Compute PCM portion of forces (2 matvecs)
subroutine gradr_fmm(params, constants, workspace, g, ygrid, fx)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: g(constants % nbasis, params % nsph), &
        & ygrid(params % ngrid, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: fx(3, params % nsph)
    ! Local variables
    integer :: i, j, indi, indj, indl, indl1, l, m, isph, igrid, ik, ksph, &
        & jsph, jsph_node
    real(dp) :: start, finish, time
    integer :: inear, inode, jnode
    real(dp) :: fac, fl, gg, c(3), vki(3), vvki, tki, gg3(3), tmp1, tmp2, tmp_gg
    real(dp) :: tlow, thigh
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    real(dp), external :: ddot, dnrm2
    real(dp) :: work(params % lmax+2)
    !real(dp) :: l2g(params % ngrid, params % nsph)
    zx_coord_transform = 0
    zx_coord_transform(3, 2) = 1
    zx_coord_transform(2, 3) = 1
    zx_coord_transform(1, 1) = 1
    zy_coord_transform = 0
    zy_coord_transform(1, 2) = 1
    zy_coord_transform(2, 1) = 1
    zy_coord_transform(3, 3) = 1
    tlow  = one - pt5*(one - params % se)*params % eta
    thigh = one + pt5*(one + params % se)*params % eta
    fx = zero
    !! Scale input harmonics at first
    workspace % tmp_sph(1, :) = zero
    indl = 2
    do l = 1, params % lmax
        indl1 = (l+1)**2
        workspace % tmp_sph(indl:indl1, :) = l * g(indl:indl1, :)
        indl = indl1 + 1
    end do
    !! Compute gradient of M2M of tmp_sph harmonics at the origin and store it
    !! in tmp_sph_grad. tmp_sph_grad(:, 1, :), tmp_sph_grad(:, 2, :) and
    !! tmp_sph_grad(:, 3, :) correspond to the OX, OY and OZ axes. Variable
    !! tmp_sph2 is a temporary workspace here.
    call tree_grad_m2m(params, constants, workspace % tmp_sph, &
        & workspace % tmp_sph_grad, workspace % tmp_sph2)
    !! Adjoint full FMM matvec to get output multipole expansions from input
    !! external grid points. It is used to compute R_i^B fast as a contraction
    !! of a gradient stored in tmp_sph_grad and a result of adjoint matvec.
    ! Adjoint integration from spherical harmonics to grid points is not needed
    ! here as ygrid already contains grid values, we just need to scale it by
    ! weights of grid points
    do isph = 1, params % nsph
        workspace % tmp_grid(:, isph) = ygrid(:, isph) * &
            & constants % wgrid(:) * constants % ui(:, isph)
    end do
    ! Adjoint FMM with output tmp_sph2(:, :) which stores coefficients of
    ! harmonics of degree up to lmax+1
    call tree_m2p_adj(params, constants, params % lmax+1, one, &
        & workspace % tmp_grid, zero, workspace % tmp_sph2)
    call tree_l2p_adj(params, constants, one, workspace % tmp_grid, zero, &
        & workspace % tmp_node_l)
    call tree_l2l_rotation_adj(params, constants, workspace % tmp_node_l)
    call tree_m2l_rotation_adj(params, constants, workspace % tmp_node_l, &
        & workspace % tmp_node_m)
    call tree_m2m_rotation_adj(params, constants, workspace % tmp_node_m)
    ! Properly load adjoint multipole harmonics into tmp_sph2 that holds
    ! harmonics of a degree up to lmax+1
    if(params % lmax+1 .lt. params % pm) then
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_sph2(:, isph) = workspace % tmp_sph2(:, isph) + &
                & workspace % tmp_node_m(1:constants % grad_nbasis, inode)
        end do
    else
        indl = (params % pm+1)**2
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_sph2(1:indl, isph) = &
                & workspace % tmp_sph2(1:indl, isph) + &
                & workspace % tmp_node_m(:, inode)
        end do
    end if
    ! Compute second term of R_i^B as a contraction
    do isph = 1, params % nsph
        call dgemv('T', constants % grad_nbasis, 3, one, &
            & workspace % tmp_sph_grad(1, 1, isph), constants % grad_nbasis, &
            & workspace % tmp_sph2(1, isph), 1, zero, fx(1, isph), 1)
    end do
    !! Direct far-field FMM matvec to get output local expansions from input
    !! multipole expansions. It will be used in R_i^A.
    !! As of now I compute potential at all external grid points, improved
    !! version shall only compute it at external points in a switch region
    ! Load input harmonics into tree data
    if(params % lmax .lt. params % pm) then
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_m(:constants % nbasis, inode) = &
                & workspace % tmp_sph(:, isph)
            workspace % tmp_node_m(constants % nbasis+1:, inode) = zero
        end do
    else
        indl = (params % pm+1)**2
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_m(:, inode) = workspace % tmp_sph(1:indl, isph)
        end do
    end if
    ! Perform direct FMM matvec to all external grid points
    call tree_m2m_rotation(params, constants, workspace % tmp_node_m)
    call tree_m2l_rotation(params, constants, workspace % tmp_node_m, &
        & workspace % tmp_node_l)
    call tree_l2l_rotation(params, constants, workspace % tmp_node_l)
    call tree_l2p(params, constants, one, workspace % tmp_node_l, zero, &
        & workspace % tmp_grid)
    call tree_m2p(params, constants, params % lmax, one, &
        & workspace % tmp_sph, one, workspace % tmp_grid)
    !! Compute gradients of L2L if pl > 0
    if (params % pl .gt. 0) then
        call tree_grad_l2l(params, constants, workspace % tmp_node_l, &
            & workspace % tmp_sph_l_grad, workspace % tmp_sph_l)
    end if
    !! Diagonal update of computed grid values, that is needed for R^C, a part
    !! of R^A and a part of R^B
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & pt5, constants % vgrid2, constants % vgrid_nbasis, g, &
        & constants % nbasis, -one, workspace % tmp_grid, params % ngrid)
    !! Scale temporary grid points by corresponding Lebedev weights and ygrid
    do igrid = 1, params % ngrid
        do isph = 1, params % nsph
            workspace % tmp_grid(igrid, isph) = &
                & workspace % tmp_grid(igrid, isph) *  constants % wgrid(igrid) * &
                & ygrid(igrid, isph)
        end do
    end do
    !! Compute all terms of grad_i(R). The second term of R_i^B is already
    !! taken into account and the first term is computed together with R_i^C.
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            ! Loop over all neighbouring spheres
            do ik = constants % inl(isph), constants % inl(isph+1) - 1
                ksph = constants % nl(ik)
                ! Only consider external grid points
                if(constants % ui(igrid, ksph) .eq. zero) cycle
                ! build geometrical quantities
                c = params % csph(:, ksph) + &
                    & params % rsph(ksph)*constants % cgrid(:, igrid)
                vki = c - params % csph(:, isph)
                !vvki = sqrt(vki(1)*vki(1) + vki(2)*vki(2) + &
                !    & vki(3)*vki(3))
                vvki = dnrm2(3, vki, 1)
                tki = vvki / params % rsph(isph)
                ! Only consider such points where grad U is non-zero
                if((tki.le.tlow) .or. (tki.ge.thigh)) cycle
                ! This is entire R^C and the first R^B component (grad_i of U
                ! of a sum of R_kj for index inequality j!=k)
                ! Indexes k and j are flipped compared to the paper
                gg = workspace % tmp_grid(igrid, ksph)
                ! Compute grad_i component of forces using precomputed
                ! potential gg
                !fx(:, isph) = fx(:, isph) - &
                !    & dfsw(tki, params % se, params % eta)/ &
                !    & params % rsph(isph)*constants % wgrid(igrid)*gg* &
                !    & ygrid(igrid, ksph)*(vki/vvki)
                fx(:, isph) = fx(:, isph) - &
                    & dfsw(tki, params % se, params % eta)/ &
                    & params % rsph(isph)*gg*(vki/vvki)
            end do
            ! contribution from the sphere itself
            if((constants % ui(igrid,isph).gt.zero) .and. &
                & (constants % ui(igrid,isph).lt.one)) then
                ! R^A component (grad_i of U of a sum of R_ij for index
                ! inequality j!=i)
                ! Indexes k and j are flipped compared to the paper
                gg = workspace % tmp_grid(igrid, isph)
                ! Compute grad_i component of forces using precomputed
                ! potential gg
                !fx(:, isph) = fx(:, isph) + constants % wgrid(igrid)*gg* &
                !    & ygrid(igrid, isph)*constants % zi(:, igrid, isph)
                fx(:, isph) = fx(:, isph) + gg*constants % zi(:, igrid, isph)
            end if
            if (constants % ui(igrid, isph) .gt. zero) then
                ! Another R^A component (grad_i of potential of a sum of R_ij
                ! for index inequality j!=i)
                ! Indexes k and j are flipped compared to the paper
                ! In case pl=0 MKL does not make gg3 zero reusing old value of
                ! gg3, so we have to clear it manually
                gg3 = zero
                call dgemv('T', params % pl**2, 3, one, &
                    & workspace % tmp_sph_l_grad(1, 1, isph), &
                    & (params % pl+1)**2, constants % vgrid2(1, igrid), 1, &
                    & zero, gg3, 1)
                ! Gradient of the near-field potential is a gradient of
                ! multipole expansion
                inode = constants % snode(isph)
                do inear = constants % snear(inode), constants % snear(inode+1)-1
                    jnode = constants % near(inear)
                    do jsph_node = constants % cluster(1, jnode), &
                        & constants % cluster(2, jnode)
                        jsph = constants % order(jsph_node)
                        if (isph .eq. jsph) cycle
                        c = params % csph(:, isph) + &
                            & params % rsph(isph)*constants % cgrid(:, igrid)
                        call fmm_m2p_work(c-params % csph(:, jsph), &
                            & params % rsph(jsph), params % lmax+1, &
                            & constants % vscales_rel, one, &
                            & workspace % tmp_sph_grad(:, 1, jsph), zero, &
                            & tmp_gg, work)
                        gg3(1) = gg3(1) + tmp_gg
                        call fmm_m2p_work(c-params % csph(:, jsph), &
                            & params % rsph(jsph), params % lmax+1, &
                            & constants % vscales_rel, one, &
                            & workspace % tmp_sph_grad(:, 2, jsph), zero, &
                            & tmp_gg, work)
                        gg3(2) = gg3(2) + tmp_gg
                        call fmm_m2p_work(c-params % csph(:, jsph), &
                            & params % rsph(jsph), params % lmax+1, &
                            & constants % vscales_rel, one, &
                            & workspace % tmp_sph_grad(:, 3, jsph), zero, &
                            & tmp_gg, work)
                        gg3(3) = gg3(3) + tmp_gg
                    end do
                end do
                ! Accumulate all computed forces
                fx(:, isph) = fx(:, isph) - constants % wgrid(igrid)*gg3* &
                    & ygrid(igrid, isph)*constants % ui(igrid, isph)
            end if
        end do
    end do
end subroutine gradr_fmm

end module ddx_operators

