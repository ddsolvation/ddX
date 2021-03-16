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

!> Compute potential in cavity points
!!
!! TODO: make use of FMM here
!! AJ: Added gradphi for ddLPB
subroutine mkrhs(ddx_data, phi_cav, gradphi_cav, psi)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    ! Outputs
    real(dp), intent(out) :: phi_cav(ddx_data % ncav)
    real(dp), intent(out) :: gradphi_cav(3, ddx_data % ncav)
    real(dp), intent(out) :: psi(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    integer :: isph, igrid, icav, inode, inear, jnear, jnode, jsph
    real(dp) :: d(3), v, dnorm, gradv(3), epsp=one
    real(dp) :: grid_grad(ddx_data % ngrid, 3, ddx_data % nsph)
    real(dp), external :: dnrm2
    ! In case FMM is disabled compute phi and gradphi at cavity points by a
    ! naive double loop of a quadratic complexity
    if (ddx_data % fmm .eq. 0) then
        do icav = 1, ddx_data % ncav
            v = zero
            gradv = zero
            do isph = 1, ddx_data % nsph
                d = ddx_data % ccav(:, icav) - ddx_data % csph(:, isph)
                dnorm = dnrm2(3, d, 1)
                v = v + ddx_data % charge(isph)/dnorm
                gradv = gradv - ddx_data % charge(isph)*d/(dnorm**3)
            end do
            phi_cav(icav) = v
            gradphi_cav(:,icav) = gradv
        end do
    ! Use the FMM otherwise
    else
        ! P2M step from centers of harmonics
        do isph = 1, ddx_data % nsph
            inode = ddx_data % snode(isph)
            ddx_data % tmp_sph(1, isph) = ddx_data % charge(isph) / &
                & ddx_data % rsph(isph) / sqrt4pi
            ddx_data % tmp_sph(2:, isph) = zero
            ddx_data % tmp_node_m(1, inode) = ddx_data % tmp_sph(1, isph)
            ddx_data % tmp_node_m(2:, inode) = zero
        end do
        ! M2M, M2L and L2L translations
        if(ddx_data % fmm_precompute .eq. 1) then
            call tree_m2m_reflection_use_mat(ddx_data, ddx_data % tmp_node_m)
            call tree_m2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_m, &
                & ddx_data % tmp_node_l)
            call tree_l2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_l)
            call tree_l2p(ddx_data, one, ddx_data % tmp_node_l, zero, &
                & ddx_data % tmp_grid)
            call tree_m2p_use_mat(ddx_data, ddx_data % lmax, one, &
                & ddx_data % tmp_sph, one, ddx_data % tmp_grid)
        else
            call tree_m2m_rotation(ddx_data, ddx_data % tmp_node_m)
            call tree_m2l_rotation(ddx_data, ddx_data % tmp_node_m, &
                & ddx_data % tmp_node_l)
            call tree_l2l_rotation(ddx_data, ddx_data % tmp_node_l)
            call tree_l2p(ddx_data, one, ddx_data % tmp_node_l, zero, &
                & ddx_data % tmp_grid)
            call tree_m2p(ddx_data, ddx_data % lmax, one, ddx_data % tmp_sph, &
                & one, ddx_data % tmp_grid)
        end if
        ! Potential from each sphere to its own grid points
        call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % l2grid, &
            & ddx_data % vgrid_nbasis, ddx_data % tmp_sph, ddx_data % nbasis, &
            & one, ddx_data % tmp_grid, ddx_data % ngrid)
        ! Rearrange potential from all grid points to external only
        icav = 0
        do isph = 1, ddx_data % nsph
            do igrid = 1, ddx_data % ngrid
                ! Do not count internal grid points
                if(ddx_data % ui(igrid, isph) .eq. zero) cycle
                icav = icav + 1
                phi_cav(icav) = ddx_data % tmp_grid(igrid, isph)
            end do
        end do
        ! Now compute near-field FMM gradients
        ! Cycle over all spheres
        do isph = 1, ddx_data % nsph
            grid_grad(:, :, isph) = zero
            ! Cycle over all external grid points
            do igrid = 1, ddx_data % ngrid
                if(ddx_data % ui(igrid, isph) .eq. zero) cycle
                ! Cycle over all near-field admissible pairs of spheres,
                ! including pair (isph, isph) which is a self-interaction
                inode = ddx_data % snode(isph)
                do jnear = ddx_data % snear(inode), ddx_data % snear(inode+1)-1
                    ! Near-field interactions are possible only between leaf
                    ! nodes, which must contain only a single input sphere
                    jnode = ddx_data % near(jnear)
                    jsph = ddx_data % order(ddx_data % cluster(1, jnode))
                    d = ddx_data % csph(:, isph) + &
                        & ddx_data % cgrid(:, igrid)*ddx_data % rsph(isph) - &
                        & ddx_data % csph(:, jsph)
                    dnorm = dnrm2(3, d, 1)
                    grid_grad(igrid, :, isph) = grid_grad(igrid, :, isph) - &
                        & ddx_data % charge(jsph)*d/(dnorm**3)
                end do
            end do
        end do
        ! Take into account far-field FMM gradients only if pl > 0
        if (ddx_data % pl .gt. 0) then
            ! Get gradient of L2L
            call tree_grad_l2l(ddx_data, ddx_data % tmp_node_l, &
                & ddx_data % tmp_sph_l_grad, ddx_data % tmp_sph_l)
            ! Apply L2P for every axis with -1 multiplier since grad over
            ! target point is equal to grad over source point
            call dgemm('T', 'N', ddx_data % ngrid, 3*ddx_data % nsph, &
                & (ddx_data % pl)**2, -one, ddx_data % l2grid, &
                & ddx_data % vgrid_nbasis, ddx_data % tmp_sph_l_grad, &
                & (ddx_data % pl+1)**2, one, grid_grad, &
                & ddx_data % ngrid)
        end if
        ! Copy output for external grid points only
        icav = 0
        do isph = 1, ddx_data % nsph
            do igrid = 1, ddx_data % ngrid
                if(ddx_data % ui(igrid, isph) .eq. zero) cycle
                icav = icav + 1
                gradphi_cav(:, icav) = grid_grad(igrid, :, isph)
            end do
        end do
    end if
    ! Vector psi
    psi(2:, :) = zero
    do isph = 1, ddx_data % nsph
        psi(1, isph) = sqrt4pi * ddx_data % charge(isph)
    end do
end subroutine mkrhs

!> Apply single layer operator to spherical harmonics
!!
!! Diagonal blocks are not counted here.
subroutine lx(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    integer :: isph, istatus
    real(dp), allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)
    ! Allocate workspaces
    allocate(pot(ddx_data % ngrid), vplm(ddx_data % nbasis), &
        & basloc(ddx_data % nbasis), vcos(ddx_data % lmax+1), &
        & vsin(ddx_data % lmax+1) , stat=istatus )
    if ( istatus.ne.0 ) then
        write(*,*) 'lx: allocation failed !'
        stop
    end if
    ! Debug printing
    if (ddx_data % iprint .ge. 5) then
        call prtsph('X', ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, &
            & x)
    end if
    ! Initialize
    y = zero
    ! !!$omp parallel do default(shared) private(isph,pot,basloc,vplm,vcos,vsin) &
    ! !!$omp schedule(dynamic)
    do isph = 1, ddx_data % nsph
        ! Compute NEGATIVE action of off-digonal blocks
        call calcv(ddx_data, .false., isph, pot, x, basloc, vplm, vcos, vsin)
        call intrhs(ddx_data % iprint, ddx_data % ngrid, ddx_data % lmax, &
            & ddx_data % vwgrid, ddx_data % vgrid_nbasis, isph, pot, y(:,isph))
        ! Action of off-diagonal blocks
        y(:,isph) = - y(:,isph)
    end do
    ! Debug printing
    if (ddx_data % iprint .ge. 5) then
        call prtsph('LX (off diagonal)', ddx_data % nbasis, ddx_data % lmax, &
            & ddx_data % nsph, 0, y)
    end if
    ! Deallocate workspaces
    deallocate(pot, basloc, vplm, vcos, vsin , stat=istatus)
    if (istatus .ne. 0) then
        write(*,*) 'lx: allocation failed !'
        stop
    endif
end subroutine lx

!> Apply adjoint single layer operator to spherical harmonics
!!
!! Diagonal blocks are not counted here.
subroutine lstarx(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    integer :: isph, igrid, istatus
    real(dp), allocatable :: xi(:,:), vplm(:), basloc(:), vcos(:), vsin(:)
    ! Allocate workspaces
    allocate(xi(ddx_data % ngrid, ddx_data % nsph), vplm(ddx_data % nbasis), &
        & basloc(ddx_data % nbasis), vcos(ddx_data % lmax+1), &
        & vsin(ddx_data % lmax+1), stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) 'lstarx: allocation failed!'
        stop
    endif
    if (ddx_data % iprint .ge. 5) then
        call prtsph('X', ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, &
            & x)
    end if
    ! Initilize
    y = zero
    ! !!$omp parallel do default(shared) private(isph,ig)
    !! Expand x over spherical harmonics
    ! Loop over spheres      
    do isph = 1, ddx_data % nsph
        ! Loop over grid points
        do igrid = 1, ddx_data % ngrid
            xi(igrid, isph) = dot_product(x(:, isph), &
                & ddx_data % vgrid(:ddx_data % nbasis, igrid))
        end do
    end do
    !! Compute action
    ! Loop over spheres
    ! !!$omp parallel do default(shared) private(isph,basloc,vplm,vcos,vsin) &
    ! !!$omp schedule(dynamic)
    do isph = 1, ddx_data % nsph
        ! Compute NEGATIVE action of off-digonal blocks
        call adjrhs(ddx_data, isph, xi, y(:, isph), basloc, vplm, vcos, vsin)
        y(:, isph) = - y(:, isph)
    end do
    if (ddx_data % iprint .ge. 5) then
        call prtsph('L*X (off-diagonal)', ddx_data % nbasis, ddx_data % lmax, &
            & ddx_data % nsph, 0, y)
    end if
    ! Deallocate workspaces
    deallocate( xi, basloc, vplm, vcos, vsin , stat=istatus )
    if ( istatus.ne.0 ) then
        write(*,*) 'lstarx: allocation failed !'
        stop
    endif
end subroutine lstarx

!> Diagonal preconditioning for Lx operator
!!
!! Applies inverse diagonal (block) of the L matrix
subroutine ldm1x(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    integer :: l, ind
    ! Loop over harmonics
    do l = 0, ddx_data % lmax
        ind = l*l + l + 1
        y(ind-l:ind+l, :) = x(ind-l:ind+l, :) * (ddx_data % vscales(ind)**2)
    end do
end subroutine ldm1x

!> Apply double layer operator to spherical harmonics
!!
!! @param[in] ddx_data
subroutine dx(ddx_data, do_diag, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Select implementation
    if (ddx_data % fmm .eq. 0) then
        call dx_dense(ddx_data, do_diag, x, y)
    else
        call dx_fmm(ddx_data, do_diag, x, y)
    end if
end subroutine dx

!> Baseline implementation of double layer operator
!!
!! @param[in] ddx_data
subroutine dx_dense(ddx_data, do_diag, x, y)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    real(dp), allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
    real(dp) :: c(3), vij(3), sij(3)
    real(dp) :: vvij, tij, tt, f, f1, rho, ctheta, stheta, cphi, sphi
    integer :: its, isph, jsph, l, m, ind, lm, istatus
    real(dp), external :: dnrm2
    ! Allocate temporaries
    allocate(vts(ddx_data % ngrid), vplm(ddx_data % nbasis), &
        & basloc(ddx_data % nbasis),vcos(ddx_data % lmax+1), &
        & vsin(ddx_data % lmax+1), stat=istatus)
    if (istatus.ne.0) then
        write(6,*) 'dx: allocation failed !'
        stop
    end if
    y = zero
    ! this loop is easily parallelizable
    ! !!$omp parallel do default(none) schedule(dynamic) &
    ! !!$omp private(isph,its,jsph,basloc,vplm,vcos,vsin,vij, &
    ! !!$omp vvij,tij,sij,tt,l,ind,f,m,vts,c) &
    ! !!$omp shared(nsph,ngrid,ui,csph,rsph,grid, &
    ! !!$omp lmax,fourpi,dodiag,x,y,basis)
    do isph = 1, ddx_data % nsph
        ! compute the "potential" from the other spheres
        ! at the exposed lebedv points of the i-th sphere 
        vts = zero
        do its = 1, ddx_data % ngrid
            if (ddx_data % ui(its,isph).gt.zero) then
                c = ddx_data % csph(:,isph) + ddx_data % rsph(isph)* &
                    & ddx_data % cgrid(:,its)
                do jsph = 1, ddx_data % nsph
                    if (jsph.ne.isph) then
                        ! build the geometrical variables
                        vij = c - ddx_data % csph(:,jsph)
                        !vvij = sqrt(dot_product(vij,vij))
                        vvij = dnrm2(3, vij, 1)
                        tij = vvij / ddx_data % rsph(jsph)
                        sij = vij/vvij 
                        ! build the local basis
                        call ylmbas(sij, rho, ctheta, stheta, cphi, sphi, &
                            & ddx_data % lmax, ddx_data % vscales, basloc, &
                            & vplm, vcos, vsin)
                        ! with all the required stuff, finally compute
                        ! the "potential" at the point 
                        tt = one/tij 
                        do l = 0, ddx_data % lmax
                            ind = l*l + l + 1
                            f = fourpi*dble(l)/(two*dble(l) + one)*tt
                            do m = -l, l
                                vts(its) = vts(its) + f*x(ind + m,jsph) * &
                                    & basloc(ind + m)
                            end do
                            tt = tt/tij
                        end do
                    else if (do_diag .eq. 1) then
                        ! add the diagonal contribution
                        do l = 0, ddx_data % lmax
                            ind = l*l + l + 1
                            f = (two*dble(l) + one)/fourpi
                            do m = -l, l
                                vts(its) = vts(its) - pt5*x(ind + m,isph) * &
                                    & ddx_data % vgrid(ind + m,its)/f
                            end do
                        end do
                    end if 
                end do
                !if(isph .eq. 1)
                vts(its) = ddx_data % ui(its,isph)*vts(its) 
            end if
        end do
        ! now integrate the potential to get its modal representation
        call intrhs(ddx_data % iprint, ddx_data % ngrid, ddx_data % lmax, &
            & ddx_data % vwgrid, ddx_data % vgrid_nbasis, isph, vts, y(:,isph))
    end do
    ! Clean up temporary data
    deallocate(vts,vplm,basloc,vcos,vsin,stat=istatus)
    if (istatus.ne.0) then
        write(6,*) 'dx: deallocation failed !'
        stop
    end if
end subroutine dx_dense

!> FMM-accelerated implementation of double layer operator
!!
!! @param[in] ddx_data
subroutine dx_fmm(ddx_data, do_diag, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    integer :: isph, inode, l, indl, indl1, m
    real(dp) :: finish_time, start_time
    ! Scale input harmonics at first
    ddx_data % tmp_sph(1, :) = zero
    indl = 2
    do l = 1, ddx_data % lmax
        indl1 = (l+1)**2
        ddx_data % tmp_sph(indl:indl1, :) = l * x(indl:indl1, :)
        indl = indl1 + 1
    end do
    ! Load input harmonics into tree data
    if(ddx_data % lmax .lt. ddx_data % pm) then
        do isph = 1, ddx_data % nsph
            inode = ddx_data % snode(isph)
            ddx_data % tmp_node_m(1:ddx_data % nbasis, inode) = &
                & ddx_data % tmp_sph(:, isph)
            ddx_data % tmp_node_m(ddx_data % nbasis+1:, inode) = zero
        end do
    else
        indl = (ddx_data % pm+1)**2
        do isph = 1, ddx_data % nsph
            inode = ddx_data % snode(isph)
            ddx_data % tmp_node_m(:, inode) = ddx_data % tmp_sph(:indl, isph)
        end do
    end if
    ! Do FMM operations
    if(ddx_data % fmm_precompute .eq. 1) then
        call tree_m2m_reflection_use_mat(ddx_data, ddx_data % tmp_node_m)
        call tree_m2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_m, &
            & ddx_data % tmp_node_l)
        call tree_l2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_l)
        call tree_l2p(ddx_data, one, ddx_data % tmp_node_l, zero, &
            & ddx_data % tmp_grid)
        call tree_m2p_use_mat(ddx_data, ddx_data % lmax, one, &
            & ddx_data % tmp_sph, one, ddx_data % tmp_grid)
    else
        call tree_m2m_rotation(ddx_data, ddx_data % tmp_node_m)
        call tree_m2l_rotation(ddx_data, ddx_data % tmp_node_m, &
            & ddx_data % tmp_node_l)
        call tree_l2l_rotation(ddx_data, ddx_data % tmp_node_l)
        call tree_l2p(ddx_data, one, ddx_data % tmp_node_l, zero, &
            & ddx_data % tmp_grid)
        call tree_m2p(ddx_data, ddx_data % lmax, one, ddx_data % tmp_sph, one, &
            & ddx_data % tmp_grid)
    end if
    ! Apply diagonal contribution if needed
    if(do_diag .eq. 1) then
        call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, -pt5, ddx_data % l2grid, &
            & ddx_data % vgrid_nbasis, x, ddx_data % nbasis, one, &
            & ddx_data % tmp_grid, ddx_data % ngrid)
    end if
    ! Multiply by characteristic function
    ddx_data % tmp_grid = ddx_data % tmp_grid * ddx_data % ui
    ! now integrate the potential to get its modal representation
    ! output y is overwritten here
    call dgemm('N', 'N', ddx_data % nbasis, ddx_data % nsph, ddx_data % ngrid, &
        & one, ddx_data % vwgrid, ddx_data % vgrid_nbasis, ddx_data % tmp_grid, &
        & ddx_data % ngrid, zero, y, ddx_data % nbasis)
end subroutine dx_fmm

!> Apply adjoint double layer operator to spherical harmonics
!!
!! @param[in] ddx_data
subroutine dstarx(ddx_data, do_diag, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Select implementation
    if (ddx_data % fmm .eq. 0) then
        call dstarx_dense(ddx_data, do_diag, x, y)
    else
        call dstarx_fmm(ddx_data, do_diag, x, y)
    end if
end subroutine dstarx

!> Baseline implementation of adjoint double layer operator
!!
!!
subroutine dstarx_dense(ddx_data, do_diag, x, y)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    real(dp), allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
    real(dp) :: c(3), vji(3), sji(3)
    real(dp) :: vvji, tji, fourpi, tt, f, f1, rho, ctheta, stheta, cphi, sphi
    integer :: its, isph, jsph, l, m, ind, lm, istatus
    real(dp), external :: dnrm2
    ! Allocate temporaries
    allocate(vts(ddx_data % ngrid), vplm(ddx_data % nbasis), &
        & basloc(ddx_data % nbasis),vcos(ddx_data % lmax+1), &
        & vsin(ddx_data % lmax+1), stat=istatus)
    if (istatus.ne.0) then
        write(6,*) 'dx: allocation failed !'
        stop
    end if
    y = zero
    ! this loop is easily parallelizable
    ! !!$omp parallel do default(none) schedule(dynamic) &
    ! !!$omp private(isph,its,jsph,basloc,vplm,vcos,vsin,vji, &
    ! !!$omp vvji,tji,sji,tt,l,ind,f,m,vts,c) &
    ! !!$omp shared(nsph,ngrid,ui,csph,rsph,grid,facl, &
    ! !!$omp lmax,fourpi,dodiag,x,y,basis,w,nbasis)
    do isph = 1, ddx_data % nsph
        do jsph = 1, ddx_data % nsph
            if (jsph.ne.isph) then
                do its = 1, ddx_data % ngrid
                    if (ddx_data % ui(its,jsph).gt.zero) then
                        ! build the geometrical variables
                        vji = ddx_data % csph(:,jsph) + ddx_data % rsph(jsph) * &
                            & ddx_data % cgrid(:,its) - ddx_data % csph(:,isph)
                        !vvji = sqrt(dot_product(vji,vji))
                        vvji = dnrm2(3, vji, 1)
                        tji = vvji/ddx_data % rsph(isph)
                        sji = vji/vvji
                        ! build the local basis
                        call ylmbas(sji, rho, ctheta, stheta, cphi, sphi, &
                            & ddx_data % lmax, ddx_data % vscales, basloc, &
                            & vplm, vcos, vsin)
                        tt = ddx_data % ui(its,jsph)*dot_product(ddx_data % vwgrid(:,its),x(:,jsph))/tji
                        do l = 0, ddx_data % lmax
                            ind = l*l + l + 1
                            f = dble(l)*tt/ ddx_data % vscales(ind)**2
                            do m = -l, l
                                y(ind+m,isph) = y(ind+m,isph) + f*basloc(ind+m)
                            end do
                            tt = tt/tji
                        end do
                    end if
                end do
            else if (do_diag .eq. 1) then
                do its = 1, ddx_data % ngrid
                    f = pt5*ddx_data % ui(its,jsph)*dot_product(ddx_data % vwgrid(:,its),x(:,jsph))
                    do l = 0, ddx_data % lmax
                        ind = l*l + l + 1
                        y(ind-l:ind+l,isph) = y(ind-l:ind+l,isph) - &
                            & f*ddx_data % vgrid(ind-l:ind+l,its)/ &
                            & ddx_data % vscales(ind)**2
                    end do
                end do
            end if
        end do
    end do
    deallocate(vts,vplm,basloc,vcos,vsin,stat=istatus)
    if (istatus.ne.0) then
        write(6,*) 'dx: deallocation failed !'
        stop
    end if
end subroutine dstarx_dense

!> FMM-accelerated implementation of adjoint double layer operator
!!
!!
subroutine dstarx_fmm(ddx_data, do_diag, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    integer :: isph, inode, l, indl, indl1, m
    real(dp) :: finish_time, start_time
    ! Adjoint integration
    call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, ddx_data % nbasis, &
        & one, ddx_data % vwgrid, ddx_data % vgrid_nbasis, x, &
        & ddx_data % nbasis, zero, ddx_data % tmp_grid, ddx_data % ngrid)
    ! Multiply by characteristic function
    ddx_data % tmp_grid = ddx_data % tmp_grid * ddx_data % ui
    ! Do FMM operations adjointly
    if(ddx_data % fmm_precompute .eq. 1) then
        call tree_m2p_use_mat_adj(ddx_data, ddx_data % lmax, one, &
            & ddx_data % tmp_grid, zero, y)
        call tree_l2p_adj(ddx_data, one, ddx_data % tmp_grid, zero, &
            & ddx_data % tmp_node_l)
        call tree_l2l_reflection_use_mat_adj(ddx_data, ddx_data % tmp_node_l)
        call tree_m2l_reflection_use_mat_adj(ddx_data, ddx_data % tmp_node_l, &
            & ddx_data % tmp_node_m)
        call tree_m2m_reflection_use_mat_adj(ddx_data, ddx_data % tmp_node_m)
    else
        call tree_m2p_adj(ddx_data, ddx_data % lmax, one, ddx_data % tmp_grid, &
            & zero, y)
        call tree_l2p_adj(ddx_data, one, ddx_data % tmp_grid, zero, &
            & ddx_data % tmp_node_l)
        call tree_l2l_rotation_adj(ddx_data, ddx_data % tmp_node_l)
        call tree_m2l_rotation_adj(ddx_data, ddx_data % tmp_node_l, &
            & ddx_data % tmp_node_m)
        call tree_m2m_rotation_adj(ddx_data, ddx_data % tmp_node_m)
    end if
    ! Adjointly move tree multipole harmonics into output
    if(ddx_data % lmax .lt. ddx_data % pm) then
        do isph = 1, ddx_data % nsph
            inode = ddx_data % snode(isph)
            y(:, isph) = y(:, isph) + &
                & ddx_data % tmp_node_m(1:ddx_data % nbasis, inode)
        end do
    else
        indl = (ddx_data % pm+1)**2
        do isph = 1, ddx_data % nsph
            inode = ddx_data % snode(isph)
            y(1:indl, isph) = y(1:indl, isph) + ddx_data % tmp_node_m(:, inode)
        end do
    end if
    ! Scale output harmonics at last
    y(1, :) = zero
    indl = 2
    do l = 1, ddx_data % lmax
        indl1 = (l+1)**2
        y(indl:indl1, :) = l * y(indl:indl1, :)
        indl = indl1 + 1
    end do
    ! Apply diagonal contribution if needed
    if(do_diag .eq. 1) then
        call dgemm('N', 'N', ddx_data % nbasis, ddx_data % nsph, &
            & ddx_data % ngrid, -pt5, ddx_data % l2grid, &
            & ddx_data % vgrid_nbasis, ddx_data % tmp_grid, ddx_data % ngrid, &
            & one, y, ddx_data % nbasis)
    end if
end subroutine dstarx_fmm

!> Apply \f$ R \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R x = - D x \f$ with excluded diagonal influence (blocks
!! D_ii are assumed to be zero).
!!
!! @param[in] ddx_data:
!! @param[in] x:
!! @param[out] y:
subroutine rx(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dx(ddx_data, 0, x, y)
    y = -y
end subroutine rx

!> Apply \f$ R^* \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R^* x = - D^* x \f$ with excluded diagonal influence (blocks
!! D_ii are assumed to be zero).
!!
!! @param[in] ddx_data:
!! @param[in] x:
!! @param[out] y:
subroutine rstarx(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dstarx(ddx_data, 0, x, y)
    y = -y
end subroutine rstarx

!> Apply \f$ R_\varepsilon \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\varepsilon x = (2\pi(\varepsilon + 1) / (\varepsilon
!! - 1) - D) x \f$.
!!
!! @param[in] ddx_data:
!! @param[in] x:
!! @param[out] y:
subroutine repsx(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dx(ddx_data, 1, x, y)
    ! Apply diagonal
    fac = twopi * (ddx_data % eps + one) / (ddx_data % eps - one)
    y = fac*x - y
end subroutine repsx

!> Apply \f$ R_\varepsilon^* \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R^*_\varepsilon x = (2\pi(\varepsilon + 1) / (\varepsilon
!! - 1) - D) x \f$.
!!
!! @param[in] ddx_data:
!! @param[in] x:
!! @param[out] y:
subroutine rstarepsx(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dstarx(ddx_data, 1, x, y)
    ! Apply diagonal
    fac = twopi * (ddx_data % eps + one) / (ddx_data % eps - one)
    y = fac*x - y
end subroutine rstarepsx

!> Apply \f$ R_\infty \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\infty x = (2\pi - D) x \f$.
!!
!! @param[in] ddx_data:
!! @param[in] x:
!! @param[out] y:
subroutine rinfx(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dx(ddx_data, 1, x, y)
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
subroutine rstarinfx(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(out) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dstarx(ddx_data, 1, x, y)
    ! Apply diagonal
    y = twopi*x - y
end subroutine rstarinfx

!> Apply preconditioner for 
subroutine apply_repsx_prec(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    real(dp), intent(inout) :: y(ddx_data % nbasis, ddx_data % nsph)
    integer :: isph
    ! simply do a matrix-vector product with the transposed preconditioner 
    ! !!$omp parallel do default(shared) schedule(dynamic) &
    ! !!$omp private(isph)
    do isph = 1, ddx_data % nsph
        call dgemm('N', 'N', ddx_data % nbasis, 1, ddx_data % nbasis, one, &
            & ddx_data % rx_prc(:, :, isph), ddx_data % nbasis, x(:, isph), &
            & ddx_data % nbasis, zero, y(:, isph), ddx_data % nbasis)
    end do
    !call prtsph("rx_prec x", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, x)
    !call prtsph("rx_prec y", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, y)
end subroutine apply_repsx_prec

!> Apply preconditioner for 
subroutine apply_rstarepsx_prec(ddx_data, x, y)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    real(dp), intent(in) :: x(ddx_data % nbasis, ddx_data % nsph)
    ! Output
    real(dp), intent(inout) :: y(ddx_data % nbasis, ddx_data % nsph)
    ! Local variables
    integer :: isph
    ! simply do a matrix-vector product with the transposed preconditioner 
    ! !!$omp parallel do default(shared) schedule(dynamic) &
    ! !!$omp private(isph)
    do isph = 1, ddx_data % nsph
        call dgemm('T', 'N', ddx_data % nbasis, 1, ddx_data % nbasis, one, &
            & ddx_data % rx_prc(:, :, isph), ddx_data % nbasis, x(:, isph), &
            & ddx_data % nbasis, zero, y(:, isph), ddx_data % nbasis)
    end do
    !call prtsph("rx_prec x", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, x)
    !call prtsph("rx_prec y", ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, y)
end subroutine apply_rstarepsx_prec

subroutine gradr(ddx_data, fx)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    ! Output
    real(dp), intent(out) :: fx(3, ddx_data % nsph)
    ! Check which gradr to execute
    if (ddx_data % fmm .eq. 1) then
        call gradr_fmm(ddx_data, fx)
    else
        call gradr_dense(ddx_data, fx)
    end if
end subroutine gradr

subroutine gradr_dense(ddx_data, fx)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    ! Output
    real(dp), intent(out) :: fx(3, ddx_data % nsph)
    ! Local variables
    integer :: isph
    real(dp) :: vplm(ddx_data % nbasis), vcos(ddx_data % lmax+1), &
        & vsin(ddx_data % lmax+1), basloc(ddx_data % nbasis), &
        & dbsloc(3, ddx_data % nbasis)
    ! Simply cycle over all spheres
    do isph = 1, ddx_data % nsph
        call gradr_sph(ddx_data, isph, vplm, vcos, vsin, basloc, dbsloc, fx(:, isph))
    end do
end subroutine gradr_dense

subroutine gradr_sph(ddx_data, isph, vplm, vcos, vsin, basloc, dbsloc, fx)
    ! Inputs
    type(ddx_type), intent(in) :: ddx_data
    integer, intent(in) :: isph
    real(dp), intent(out) :: vplm(ddx_data % nbasis), vcos(ddx_data % lmax+1), &
        & vsin(ddx_data % lmax+1), basloc(ddx_data % nbasis), &
        & dbsloc(3, ddx_data % nbasis), fx(3)
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

    tlow  = one - pt5*(one - ddx_data % se)*ddx_data % eta
    thigh = one + pt5*(one + ddx_data % se)*ddx_data % eta

    ! first set of contributions:
    ! diagonal block, kc and part of kb

    fx = zero
    do its = 1, ddx_data % ngrid
        ! sum over ksph in neighbors of isph
        do ik = ddx_data % inl(isph), ddx_data % inl(isph+1) - 1
            ksph = ddx_data % nl(ik)
            ! build geometrical quantities
            cx = ddx_data % csph(1,ksph) + ddx_data % rsph(ksph)*ddx_data % cgrid(1,its)
            cy = ddx_data % csph(2,ksph) + ddx_data % rsph(ksph)*ddx_data % cgrid(2,its)
            cz = ddx_data % csph(3,ksph) + ddx_data % rsph(ksph)*ddx_data % cgrid(3,its)
            vki(1) = cx - ddx_data % csph(1,isph)
            vki(2) = cy - ddx_data % csph(2,isph)
            vki(3) = cz - ddx_data % csph(3,isph)
            !vvki = sqrt(vki(1)*vki(1) + vki(2)*vki(2) + &
            !    & vki(3)*vki(3))
            vvki = dnrm2(3, vki, 1)
            tki  = vvki/ddx_data % rsph(isph)

            ! contributions involving grad i of uk come from the switching
            ! region.
            ! note: ui avoids contributions from points that are in the
            ! switching between isph and ksph but are buried in a third
            ! sphere.
            if ((tki.gt.tlow).and.(tki.lt.thigh) .and. &
                & ddx_data % ui(its,ksph).gt.zero) then
                ! other geometrical quantities
                ski = vki/vvki

                ! diagonal block kk contribution, with k in n(i)
                gg = zero
                do l = 0, ddx_data % lmax
                    ind = l*l + l + 1
                    fl = dble(l)
                    fac = twopi/(two*fl + one)
                    do m = -l, l
                        !! DEBUG comment
                        gg = gg + fac*ddx_data % vgrid(ind+m,its)*ddx_data % g(ind+m,ksph)
                    end do
                end do

                ! kc contribution
                do jsph = 1, ddx_data % nsph
                    if (jsph.ne.ksph .and. jsph.ne.isph) then 
                        vkj(1) = cx - ddx_data % csph(1,jsph)
                        vkj(2) = cy - ddx_data % csph(2,jsph)
                        vkj(3) = cz - ddx_data % csph(3,jsph)
                        vvkj = sqrt(vkj(1)*vkj(1) + vkj(2)*vkj(2) + &
                            & vkj(3)*vkj(3))
                        vvkj = dnrm2(3, vkj, 1)
                        tkj  = vvkj/ddx_data % rsph(jsph)
                        skj  = vkj/vvkj
                        call ylmbas(skj, rho, ctheta, stheta, cphi, sphi, &
                            & ddx_data % lmax, ddx_data % vscales, basloc, &
                            & vplm, vcos, vsin)
                        tt = one/tkj
                        do l = 0, ddx_data % lmax
                            ind = l*l + l + 1
                            fcl = - fourpi*dble(l)/(two*dble(l)+one)*tt
                            do m = -l, l
                                !! DEBUG comment
                                gg = gg + fcl*ddx_data % g(ind+m,jsph)*basloc(ind+m)
                            end do
                            tt = tt/tkj
                        end do
                        !call fmm_m2p(vkj, ddx_data % rsph(jsph), &
                        !    & ddx_data % lmax, ddx_data % vscales_rel, -one, &
                        !    & ddx_data % g(:, jsph), one, gg)
                    end if
                end do

                ! part of kb contribution
                call ylmbas(ski, rho, ctheta, stheta, cphi, sphi, &
                    & ddx_data % lmax, ddx_data % vscales, basloc, &
                    & vplm, vcos, vsin)
                tt = one/tki
                do l = 0, ddx_data % lmax
                    ind = l*l + l + 1
                    fcl = - four*pi*dble(l)/(two*dble(l)+one)*tt
                    do m = -l, l
                        !! DEBUG comment
                        gg = gg + fcl*ddx_data % g(ind+m,isph)*basloc(ind+m)
                    end do
                    tt = tt/tki
                end do
                !call fmm_m2p(vki, ddx_data % rsph(isph), &
                !    & ddx_data % lmax, ddx_data % vscales_rel, -one, &
                !    & ddx_data % g(:, isph), one, gg)

                ! common step, product with grad i uj
                duj = dfsw(tki,ddx_data % se, ddx_data % eta)/ddx_data % rsph(isph)
                fjj = duj*ddx_data % wgrid(its)*gg*ddx_data % ygrid(its,ksph)
                fx(1) = fx(1) - fjj*ski(1)
                fx(2) = fx(2) - fjj*ski(2)
                fx(3) = fx(3) - fjj*ski(3)
            end if
        end do

        ! diagonal block ii contribution
        if (ddx_data % ui(its,isph).gt.zero.and.ddx_data % ui(its,isph).lt.one) then
            gi = zero
            do l = 0, ddx_data % lmax
                ind = l*l + l + 1
                fl = dble(l)
                fac = twopi/(two*fl + one)
                do m = -l, l 
                    !! DEBUG comment
                    gi = gi + fac*ddx_data % vgrid(ind+m,its)*ddx_data % g(ind+m,isph)
                    !gi = gi + pt5*ddx_data % l2grid(ind+m,its)*ddx_data % g(ind+m,isph)
                end do
            end do
            !do l = 0, (ddx_data % lmax+1)**2
            !    gi = gi + ddx_data % l2grid(l, its)*ddx_data % g(l, isph)
            !end do
            !gi = pt5 * gi
            fii = ddx_data % wgrid(its)*gi*ddx_data % ygrid(its,isph)
            fx(1) = fx(1) + fii*ddx_data % zi(1,its,isph)
            fx(2) = fx(2) + fii*ddx_data % zi(2,its,isph)
            fx(3) = fx(3) + fii*ddx_data % zi(3,its,isph)
        end if
    end do

    ! second set of contributions:
    ! part of kb and ka
    do its = 1, ddx_data % ngrid

        ! run over all the spheres except isph 
        do jsph = 1, ddx_data % nsph
            if (ddx_data % ui(its,jsph).gt.zero .and. jsph.ne.isph) then
                ! build geometrical quantities
                cx = ddx_data % csph(1,jsph) + ddx_data % rsph(jsph)*ddx_data % cgrid(1,its)
                cy = ddx_data % csph(2,jsph) + ddx_data % rsph(jsph)*ddx_data % cgrid(2,its)
                cz = ddx_data % csph(3,jsph) + ddx_data % rsph(jsph)*ddx_data % cgrid(3,its)
                vji(1) = cx - ddx_data % csph(1,isph)
                vji(2) = cy - ddx_data % csph(2,isph)
                vji(3) = cz - ddx_data % csph(3,isph)
                !vvji = sqrt(vji(1)*vji(1) + vji(2)*vji(2) + &
                !    &  vji(3)*vji(3))
                vvji = dnrm2(3, vji, 1)
                tji = vvji/ddx_data % rsph(isph)
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
                call dbasis(ddx_data, sji,basloc,dbsloc,vplm,vcos,vsin)

                ! assemble the contribution
                a = zero
                tt = one/(tji)
                do l = 0, ddx_data % lmax
                    ind = l*l + l + 1
                    fl = dble(l)
                    fcl = - tt*fourpi*fl/(two*fl + one)
                    do m = -l, l
                        fac = fcl*ddx_data % g(ind+m,isph)
                        b = (fl + one)*basloc(ind+m)/(ddx_data % rsph(isph)*tji)

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
                fac = ddx_data % ui(its,jsph)*ddx_data % wgrid(its)* ddx_data % ygrid(its,jsph)
                fx(1) = fx(1) - fac*a(1)
                fx(2) = fx(2) - fac*a(2)
                fx(3) = fx(3) - fac*a(3)
            end if
        end do
    end do

    ! ka contribution
    do its = 1, ddx_data % ngrid
        cx = ddx_data % csph(1,isph) + ddx_data % rsph(isph)*ddx_data % cgrid(1,its)
        cy = ddx_data % csph(2,isph) + ddx_data % rsph(isph)*ddx_data % cgrid(2,its)
        cz = ddx_data % csph(3,isph) + ddx_data % rsph(isph)*ddx_data % cgrid(3,its)
        a = zero

        ! iterate on all the spheres except isph
        do ksph = 1, ddx_data % nsph
            if (ddx_data % ui(its,isph).gt.zero .and. ksph.ne.isph) then
                ! geometrical stuff
                vik(1) = cx - ddx_data % csph(1,ksph)
                vik(2) = cy - ddx_data % csph(2,ksph)
                vik(3) = cz - ddx_data % csph(3,ksph)
                !vvik = sqrt(vik(1)*vik(1) + vik(2)*vik(2) + & 
                !    & vik(3)*vik(3))
                vvik = dnrm2(3, vik, 1)
                tik = vvik/ddx_data % rsph(ksph)
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
                if (ddx_data % ui(its,isph).lt.one) then
                    vb(1) = ddx_data % zi(1,its,isph)
                    vb(2) = ddx_data % zi(2,its,isph)
                    vb(3) = ddx_data % zi(3,its,isph)
                end if

                ! assemble the local basis and its gradient
                !call dbasis(sik,basloc,dbsloc,vplm,vcos,vsin)
                call dbasis(ddx_data, sik,basloc,dbsloc,vplm,vcos,vsin)

                ! assemble the contribution
                tt = one/(tik)
                do l = 0, ddx_data % lmax
                    ind = l*l + l + 1
                    fl = dble(l)
                    fcl = - tt*fourpi*fl/(two*fl + one)
                        do m = -l, l
                        fac = fcl*ddx_data % g(ind+m,ksph)
                        fac = - fac*basloc(ind+m)
                        a(1) = a(1) + fac*vb(1)
                        a(2) = a(2) + fac*vb(2) 
                        a(3) = a(3) + fac*vb(3)

                        fac = ddx_data % ui(its,isph)*fcl*ddx_data % g(ind+m,ksph)
                        b = - (fl + one)*basloc(ind+m)/(ddx_data % rsph(ksph)*tik)

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
        fac = ddx_data % wgrid(its)*ddx_data % ygrid(its,isph)
        fx(1) = fx(1) - fac*a(1)
        fx(2) = fx(2) - fac*a(2)
        fx(3) = fx(3) - fac*a(3)
    end do
end subroutine gradr_sph

! Compute PCM portion of forces (2 matvecs)
subroutine gradr_fmm(ddx_data, fx)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    ! Output
    real(dp), intent(out) :: fx(3, ddx_data % nsph)
    ! Local variables
    integer :: i, j, indi, indj, indl, indl1, l, m, isph, igrid, ik, ksph, &
        & jsph, jsph_node
    real(dp) :: start, finish, time
    integer :: inear, inode, jnode
    real(dp) :: fac, fl, gg, c(3), vki(3), vvki, tki, gg3(3), tmp1, tmp2, tmp_gg
    real(dp) :: tlow, thigh
    real(dp), dimension(3, 3) :: zx_coord_transform, zy_coord_transform
    real(dp), external :: ddot, dnrm2
    real(dp) :: work((ddx_data % lmax+2)**2+3*(ddx_data % lmax+1))
    !real(dp) :: l2g(ddx_data % ngrid, ddx_data % nsph)
    zx_coord_transform = 0
    zx_coord_transform(3, 2) = 1
    zx_coord_transform(2, 3) = 1
    zx_coord_transform(1, 1) = 1
    zy_coord_transform = 0
    zy_coord_transform(1, 2) = 1
    zy_coord_transform(2, 1) = 1
    zy_coord_transform(3, 3) = 1
    tlow  = one - pt5*(one - ddx_data % se)*ddx_data % eta
    thigh = one + pt5*(one + ddx_data % se)*ddx_data % eta
    fx = zero
    !! Scale input harmonics at first
    ddx_data % tmp_sph(1, :) = zero
    indl = 2
    do l = 1, ddx_data % lmax
        indl1 = (l+1)**2
        ddx_data % tmp_sph(indl:indl1, :) = l * ddx_data % g(indl:indl1, :)
        indl = indl1 + 1
    end do
    !! Compute gradient of M2M of tmp_sph harmonics at the origin and store it
    !! in tmp_sph_grad. tmp_sph_grad(:, 1, :), tmp_sph_grad(:, 2, :) and
    !! tmp_sph_grad(:, 3, :) correspond to the OX, OY and OZ axes. Variable
    !! tmp_sph2 is a temporary workspace here.
    call tree_grad_m2m(ddx_data, ddx_data % tmp_sph, ddx_data % tmp_sph_grad, &
        & ddx_data % tmp_sph2)
    !! Adjoint full FMM matvec to get output multipole expansions from input
    !! external grid points. It is used to compute R_i^B fast as a contraction
    !! of a gradient stored in tmp_sph_grad and a result of adjoint matvec.
    ! Adjoint integration from spherical harmonics to grid points is not needed
    ! here as ygrid already contains grid values, we just need to scale it by
    ! weights of grid points
    do isph = 1, ddx_data % nsph
        ddx_data % tmp_grid(:, isph) = ddx_data % ygrid(:, isph) * &
            & ddx_data % wgrid(:) * ddx_data % ui(:, isph)
    end do
    ! Adjoint FMM with output tmp_sph2(:, :) which stores coefficients of
    ! harmonics of degree up to lmax+1
    if (ddx_data % fmm_precompute .eq. 1) then
        call tree_m2p_use_mat_adj(ddx_data, ddx_data % lmax+1, one, &
            & ddx_data % tmp_grid, zero, ddx_data % tmp_sph2)
        call tree_l2p_adj(ddx_data, one, ddx_data % tmp_grid, zero, &
            & ddx_data % tmp_node_l)
        call tree_l2l_reflection_use_mat_adj(ddx_data, ddx_data % tmp_node_l)
        call tree_m2l_reflection_use_mat_adj(ddx_data, ddx_data % tmp_node_l, &
            & ddx_data % tmp_node_m)
        call tree_m2m_reflection_use_mat_adj(ddx_data, ddx_data % tmp_node_m)
    else
        call tree_m2p_adj(ddx_data, ddx_data % lmax+1, one, ddx_data % tmp_grid, &
            & zero, ddx_data % tmp_sph2)
        call tree_l2p_adj(ddx_data, one, ddx_data % tmp_grid, zero, &
            & ddx_data % tmp_node_l)
        call tree_l2l_rotation_adj(ddx_data, ddx_data % tmp_node_l)
        call tree_m2l_rotation_adj(ddx_data, ddx_data % tmp_node_l, &
            & ddx_data % tmp_node_m)
        call tree_m2m_rotation_adj(ddx_data, ddx_data % tmp_node_m)
    end if
    ! Properly load adjoint multipole harmonics into tmp_sph2 that holds
    ! harmonics of a degree up to lmax+1
    if(ddx_data % lmax+1 .lt. ddx_data % pm) then
        do isph = 1, ddx_data % nsph
            inode = ddx_data % snode(isph)
            ddx_data % tmp_sph2(:, isph) = ddx_data % tmp_sph2(:, isph) + &
                & ddx_data % tmp_node_m(1:ddx_data % grad_nbasis, inode)
        end do
    else
        indl = (ddx_data % pm+1)**2
        do isph = 1, ddx_data % nsph
            inode = ddx_data % snode(isph)
            ddx_data % tmp_sph2(1:indl, isph) = &
                & ddx_data % tmp_sph2(1:indl, isph) + &
                & ddx_data % tmp_node_m(:, inode)
        end do
    end if
    ! Compute second term of R_i^B as a contraction
    do isph = 1, ddx_data % nsph
        call dgemv('T', ddx_data % grad_nbasis, 3, one, &
            & ddx_data % tmp_sph_grad(1, 1, isph), ddx_data % grad_nbasis, &
            & ddx_data % tmp_sph2(1, isph), 1, zero, fx(1, isph), 1)
    end do
    !! Direct far-field FMM matvec to get output local expansions from input
    !! multipole expansions. It will be used in R_i^A.
    !! As of now I compute potential at all external grid points, improved
    !! version shall only compute it at external points in a switch region
    ! Load input harmonics into tree data
    if(ddx_data % lmax .lt. ddx_data % pm) then
        do isph = 1, ddx_data % nsph
            inode = ddx_data % snode(isph)
            ddx_data % tmp_node_m(:ddx_data % nbasis, inode) = &
                & ddx_data % tmp_sph(:, isph)
            ddx_data % tmp_node_m(ddx_data % nbasis+1:, inode) = zero
        end do
    else
        indl = (ddx_data % pm+1)**2
        do isph = 1, ddx_data % nsph
            inode = ddx_data % snode(isph)
            ddx_data % tmp_node_m(:, inode) = ddx_data % tmp_sph(1:indl, isph)
        end do
    end if
    ! Perform direct FMM matvec to all external grid points
    if (ddx_data % fmm_precompute .eq. 1) then
        call tree_m2m_reflection_use_mat(ddx_data, ddx_data % tmp_node_m)
        call tree_m2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_m, &
            & ddx_data % tmp_node_l)
        call tree_l2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_l)
        call tree_l2p(ddx_data, one, ddx_data % tmp_node_l, zero, &
            & ddx_data % tmp_grid)
        call tree_m2p_use_mat(ddx_data, ddx_data % lmax, one, &
            & ddx_data % tmp_sph, one, ddx_data % tmp_grid)
    else
        call tree_m2m_rotation(ddx_data, ddx_data % tmp_node_m)
        call tree_m2l_rotation(ddx_data, ddx_data % tmp_node_m, &
            & ddx_data % tmp_node_l)
        call tree_l2l_rotation(ddx_data, ddx_data % tmp_node_l)
        call tree_l2p(ddx_data, one, ddx_data % tmp_node_l, zero, &
            & ddx_data % tmp_grid)
        call tree_m2p(ddx_data, ddx_data % lmax, one, ddx_data % tmp_sph, one, &
            & ddx_data % tmp_grid)
    end if
    !! Compute gradients of L2L if pl > 0
    if (ddx_data % pl .gt. 0) then
        call tree_grad_l2l(ddx_data, ddx_data % tmp_node_l, &
            & ddx_data % tmp_sph_l_grad, ddx_data % tmp_sph_l)
    end if
    !! Diagonal update of computed grid values, that is needed for R^C, a part
    !! of R^A and a part of R^B
    call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, ddx_data % nbasis, &
        & pt5, ddx_data % l2grid, ddx_data % vgrid_nbasis, ddx_data % g, &
        & ddx_data % nbasis, -one, ddx_data % tmp_grid, ddx_data % ngrid)
    !! Scale temporary grid points by corresponding Lebedev weights and ygrid
    do igrid = 1, ddx_data % ngrid
        do isph = 1, ddx_data % nsph
            ddx_data % tmp_grid(igrid, isph) = &
                & ddx_data % tmp_grid(igrid, isph) *  ddx_data % wgrid(igrid) * &
                & ddx_data % ygrid(igrid, isph)
        end do
    end do
    !! Compute all terms of grad_i(R). The second term of R_i^B is already
    !! taken into account and the first term is computed together with R_i^C.
    do isph = 1, ddx_data % nsph
        do igrid = 1, ddx_data % ngrid
            ! Loop over all neighbouring spheres
            do ik = ddx_data % inl(isph), ddx_data % inl(isph+1) - 1
                ksph = ddx_data % nl(ik)
                ! Only consider external grid points
                if(ddx_data % ui(igrid, ksph) .eq. zero) cycle
                ! build geometrical quantities
                c = ddx_data % csph(:, ksph) + &
                    & ddx_data % rsph(ksph)*ddx_data % cgrid(:, igrid)
                vki = c - ddx_data % csph(:, isph)
                !vvki = sqrt(vki(1)*vki(1) + vki(2)*vki(2) + &
                !    & vki(3)*vki(3))
                vvki = dnrm2(3, vki, 1)
                tki = vvki / ddx_data % rsph(isph)
                ! Only consider such points where grad U is non-zero
                if((tki.le.tlow) .or. (tki.ge.thigh)) cycle
                ! This is entire R^C and the first R^B component (grad_i of U
                ! of a sum of R_kj for index inequality j!=k)
                ! Indexes k and j are flipped compared to the paper
                gg = ddx_data % tmp_grid(igrid, ksph)
                ! Compute grad_i component of forces using precomputed
                ! potential gg
                !fx(:, isph) = fx(:, isph) - &
                !    & dfsw(tki, ddx_data % se, ddx_data % eta)/ &
                !    & ddx_data % rsph(isph)*ddx_data % wgrid(igrid)*gg* &
                !    & ddx_data % ygrid(igrid, ksph)*(vki/vvki)
                fx(:, isph) = fx(:, isph) - &
                    & dfsw(tki, ddx_data % se, ddx_data % eta)/ &
                    & ddx_data % rsph(isph)*gg*(vki/vvki)
            end do
            ! contribution from the sphere itself
            if((ddx_data % ui(igrid,isph).gt.zero) .and. &
                & (ddx_data % ui(igrid,isph).lt.one)) then
                ! R^A component (grad_i of U of a sum of R_ij for index
                ! inequality j!=i)
                ! Indexes k and j are flipped compared to the paper
                gg = ddx_data % tmp_grid(igrid, isph)
                ! Compute grad_i component of forces using precomputed
                ! potential gg
                !fx(:, isph) = fx(:, isph) + ddx_data % wgrid(igrid)*gg* &
                !    & ddx_data % ygrid(igrid, isph)*ddx_data % zi(:, igrid, isph)
                fx(:, isph) = fx(:, isph) + gg*ddx_data % zi(:, igrid, isph)
            end if
            if (ddx_data % ui(igrid, isph) .gt. zero) then
                ! Another R^A component (grad_i of potential of a sum of R_ij
                ! for index inequality j!=i)
                ! Indexes k and j are flipped compared to the paper
                call dgemv('T', ddx_data % pl**2, 3, one, &
                    & ddx_data % tmp_sph_l_grad(1, 1, isph), &
                    & (ddx_data % pl+1)**2, ddx_data % l2grid(1, igrid), 1, &
                    & zero, gg3, 1)
                ! Gradient of the near-field potential is a gradient of
                ! multipole expansion
                inode = ddx_data % snode(isph)
                do inear = ddx_data % snear(inode), ddx_data % snear(inode+1)-1
                    jnode = ddx_data % near(inear)
                    do jsph_node = ddx_data % cluster(1, jnode), &
                        & ddx_data % cluster(2, jnode)
                        jsph = ddx_data % order(jsph_node)
                        if (isph .eq. jsph) cycle
                        c = ddx_data % csph(:, isph) + &
                            & ddx_data % rsph(isph)*ddx_data % cgrid(:, igrid)
                        call fmm_m2p_work(c-ddx_data % csph(:, jsph), &
                            & ddx_data % rsph(jsph), ddx_data % lmax+1, &
                            & ddx_data % vscales_rel, one, &
                            & ddx_data % tmp_sph_grad(:, 1, jsph), zero, &
                            & tmp_gg, work)
                        gg3(1) = gg3(1) + tmp_gg
                        call fmm_m2p_work(c-ddx_data % csph(:, jsph), &
                            & ddx_data % rsph(jsph), ddx_data % lmax+1, &
                            & ddx_data % vscales_rel, one, &
                            & ddx_data % tmp_sph_grad(:, 2, jsph), zero, &
                            & tmp_gg, work)
                        gg3(2) = gg3(2) + tmp_gg
                        call fmm_m2p_work(c-ddx_data % csph(:, jsph), &
                            & ddx_data % rsph(jsph), ddx_data % lmax+1, &
                            & ddx_data % vscales_rel, one, &
                            & ddx_data % tmp_sph_grad(:, 3, jsph), zero, &
                            & tmp_gg, work)
                        gg3(3) = gg3(3) + tmp_gg
                    end do
                end do
                ! Accumulate all computed forces
                fx(:, isph) = fx(:, isph) - ddx_data % wgrid(igrid)*gg3* &
                    & ddx_data % ygrid(igrid, isph)*ddx_data % ui(igrid, isph)
            end if
        end do
    end do
end subroutine gradr_fmm

end module ddx_operators

