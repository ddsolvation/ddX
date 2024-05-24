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
! Get the solver-module
use ddx_solvers
! Get lpb core-module
use ddx_lpb_core
! Get gradient module
use ddx_gradients
implicit none

contains

!> Single layer operator matvec
subroutine lx(params, constants, workspace, x, y, ddx_error)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    !! Local variables
    integer :: isph, jsph, ij, l, ind, its, m, ijits
    real(dp) :: vij(3), tij, vvij, xij, oij, fac

    integer :: vgrid_nbasis, nbasis, ngrid, nsph, lmax
    real(dp) :: se, eta, work(params % lmax + 1), tmp_grid(params % ngrid)

    call time_push()

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue
    if (workspace % xs_time .eq. 0) continue

    nsph = params % nsph
    nbasis = constants % nbasis
    lmax = params % lmax

    if (params % matvecmem .eq. 1) then
        associate(inl => constants % inl, nl => constants % nl, &
                & l => constants % l)
            !$omp parallel do default(none) shared(inl,nl,l,x,y) &
            !$omp firstprivate(nsph,nbasis) &
            !$omp private(isph,ij,jsph) schedule(dynamic,1)
            do isph = 1, nsph
                y(:,isph) = 0.0d0
                do ij = inl(isph), inl(isph + 1) - 1
                    jsph = nl(ij)
                    call dgemv('n', nbasis, nbasis, 1.0d0, l(:,:,ij), &
                        & nbasis, x(:,jsph), 1, 1.0d0, y(:,isph), 1)
                end do
            end do
        end associate
    else
        se = params % se
        eta = params % eta
        ngrid = params % ngrid
        vgrid_nbasis = constants % vgrid_nbasis

        associate(inl => constants % inl, nl => constants % nl, &
                & csph => params % csph, rsph => params % rsph, &
                & cgrid => constants % cgrid, fi => constants % fi, &
                & vscales_rel => constants % vscales_rel, &
                & vwgrid => constants % vwgrid, &
                & ioverlap => constants % ioverlap, &
                & overlap => constants % overlap)

            !$omp parallel do default(none) schedule(dynamic,1) &
            !$omp firstprivate(nsph,se,eta,lmax,nbasis,ngrid, &
            !$omp vgrid_nbasis) shared(inl,nl,csph,rsph,cgrid,fi, &
            !$omp vscales_rel,x,y,vwgrid,ioverlap,overlap) private( &
            !$omp isph,tmp_grid,its,ij,jsph,vij,vvij,tij,xij,oij, &
            !$omp work,ijits)
            do isph = 1, nsph
                tmp_grid(:) = 0.0d0
                do ij = inl(isph), inl(isph+1)-1
                    jsph = nl(ij)
                    do ijits = ioverlap(ij), ioverlap(ij+1) - 1
                        its = overlap(ijits)
                        vij = csph(:,isph) + rsph(isph)*cgrid(:,its) &
                            & - csph(:,jsph)
                        vvij = sqrt(vij(1)*vij(1) + vij(2)*vij(2) &
                            & + vij(3)*vij(3))
                        tij  = vvij / rsph(jsph)
                        xij = fsw(tij, se, eta)
                        if (fi(its, isph).gt.1.0d0) then
                            oij = xij / fi(its, isph)
                        else
                            oij = xij
                        end if
                        call fmm_l2p_work(vij, rsph(jsph), lmax, &
                            & vscales_rel, oij, x(:, jsph), 1.0d0, &
                            & tmp_grid(its), work)
                    end do
                end do
                call ddintegrate(1, nbasis, ngrid, vwgrid, vgrid_nbasis, &
                    & tmp_grid, y(:, isph))
                y(:, isph) = - y(:, isph)
            end do
        end associate
    end if

    ! if required, add the diagonal.
    if (constants % dodiag) then
        associate(vscales => constants % vscales)
            !$omp parallel do collapse(2) default(none) &
            !$omp shared(vscales,x,y) firstprivate(nsph,lmax) &
            !$omp private(isph,l,ind,m,fac) schedule(static,10)
            do isph = 1, nsph
                do l = 0, lmax
                    ind = l*l + l + 1
                    fac = vscales(ind)**2
                    do m = -l, l
                        y(ind+m, isph) = y(ind+m, isph) + x(ind+m, isph)/fac
                    end do
                end do
            end do
        end associate
    end if
    call time_pull("lx")
end subroutine lx

!> Adjoint single layer operator matvec
subroutine lstarx(params, constants, workspace, x, y, ddx_error)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    !! Local variables
    integer :: isph, jsph, ij, indmat, igrid, l, ind, m, i, ijigrid
    real(dp) :: vji(3), vvji, tji, xji, oji, fac
    integer :: nsph, ngrid, nbasis, lmax
    real(dp) :: thigh, se, eta, work(params % lmax + 1), &
        & tmp(constants % nbasis)

    call time_push()

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    nsph = params % nsph
    nbasis = constants % nbasis

    if (params % matvecmem .eq. 1) then
        associate(inl => constants % inl, nl => constants % nl, &
                & itrnl => constants % itrnl, l => constants % l)
            !$omp parallel do default(none) schedule(dynamic,1) &
            !$omp firstprivate(nsph,nbasis) shared(x,y,inl,nl,itrnl,l) &
            !$omp private(isph,ij,jsph,indmat)
            do isph = 1, nsph
                y(:,isph) = 0.0d0
                do ij = inl(isph), inl(isph + 1) - 1
                    jsph = nl(ij)
                    indmat = itrnl(ij)
                    call dgemv('t', nbasis, nbasis, 1.0d0, &
                        & l(:,:,indmat), nbasis, x(:,jsph), 1, &
                        & 1.0d0, y(:,isph), 1)
                end do
            end do
        end associate
    else
        ngrid = params % ngrid
        thigh = one + (params % se+one)/two*params % eta
        se = params % se
        eta = params % eta
        lmax = params % lmax
        ! Expand x over spherical harmonics
        associate(vgrid => constants % vgrid, &
                & tmp_grid => workspace % tmp_grid)
            !$omp parallel do collapse(2) default(none) &
            !$omp firstprivate(nsph,ngrid,nbasis) shared(x,tmp_grid,vgrid) &
            !$omp private(isph,igrid) schedule(static,100)
            do isph = 1, nsph
                do igrid = 1, ngrid
                    tmp_grid(igrid,isph) = dot_product(x(:, isph), &
                        & vgrid(:nbasis, igrid))
                end do
            end do
        end associate
        associate(inl => constants % inl, nl => constants % nl, &
                & csph => params % csph, rsph => params % rsph, &
                & cgrid => constants % cgrid, fi => constants % fi, &
                & wgrid => constants % wgrid, &
                & tmp_grid => workspace % tmp_grid, &
                & vscales_rel => constants % vscales_rel, &
                & adj_ioverlap => constants % adj_ioverlap, &
                & adj_overlap => constants % adj_overlap)
            !$omp parallel do default(none) schedule(dynamic,1) &
            !$omp shared(inl,nl,csph,rsph,cgrid,fi,wgrid,tmp_grid,y, &
            !$omp adj_overlap,adj_ioverlap,vscales_rel) firstprivate( &
            !$omp nsph,ngrid,thigh,se,eta,lmax) private(isph,ij,jsph, &
            !$omp igrid,vji,vvji,tji,xji,oji,fac,work,i,ijigrid,tmp)
            do isph = 1, nsph
                tmp = 0.0d0
                do ij = inl(isph), inl(isph+1)-1
                    jsph = nl(ij)
                    do ijigrid = adj_ioverlap(ij), adj_ioverlap(ij+1) - 1
                        igrid = adj_overlap(ijigrid)

                        vji = csph(:,jsph) + rsph(jsph)*cgrid(:,igrid) &
                            & - csph(:,isph)
                        vvji = sqrt(vji(1)*vji(1) + vji(2)*vji(2) &
                            & + vji(3)*vji(3))
                        tji  = vvji/rsph(isph)
                        xji = fsw(tji, se, eta)
                        if (fi(igrid,jsph).gt.1.0d0) then
                            oji = xji/fi(igrid,jsph)
                        else
                            oji = xji
                        end if
                        fac = wgrid(igrid) * tmp_grid(igrid, jsph) * oji
                        call fmm_l2p_adj_work(vji, fac, rsph(isph), &
                            & lmax, vscales_rel, 1.0d0, tmp, work)
                    end do
                end do
                y(:, isph) = - tmp
            end do
        end associate
    end if
    if (constants % dodiag) then
        associate(vscales => constants % vscales)
            !$omp parallel do collapse(2) default(none) &
            !$omp shared(vscales,x,y) firstprivate(nsph,lmax) &
            !$omp private(isph,l,ind,m,fac) schedule(static,10)
            do isph = 1, nsph
                do l = 0, lmax
                    ind = l*l + l + 1
                    fac = vscales(ind)**2
                    do m = -l, l
                        y(ind+m, isph) = y(ind+m, isph) + x(ind+m, isph)/fac
                    end do
                end do
            end do
        end associate
    end if
    call time_pull("lstarx")
end subroutine lstarx

!> Diagonal preconditioning for Lx operator
!!
!! Applies inverse diagonal (block) of the L matrix
subroutine ldm1x(params, constants, workspace, x, y, ddx_error)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    !! Local variables
    integer :: isph, l, ind, nsph, lmax, m
    real(dp) :: fac

    ! dummy operation on unused interface arguments
    if ((ddx_error % flag .eq. 0) .or. &
        & (allocated(workspace % tmp_pot))) continue

    nsph = params % nsph
    lmax = params % lmax

    associate(vscales => constants % vscales)
        !$omp parallel do collapse(2) default(none) &
        !$omp shared(vscales,x,y) firstprivate(nsph,lmax) &
        !$omp private(isph,l,ind,m,fac) schedule(static,10)
        do isph = 1, nsph
            do l = 0, lmax
                ind = l*l + l + 1
                fac = vscales(ind)**2
                do m = -l, l
                    y(ind+m, isph) = x(ind+m, isph)*fac
                end do
            end do
        end do
    end associate
end subroutine ldm1x

!> Double layer operator matvec without diagonal blocks
subroutine dx(params, constants, workspace, x, y, ddx_error)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    !! Local parameter
    integer :: do_diag
    !! initialize do_diag
    do_diag = 0
    if (constants % dodiag) do_diag = 1
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dx_dense(params, constants, workspace, do_diag, x, y, ddx_error)
    else
        call dx_fmm(params, constants, workspace, do_diag, x, y, ddx_error)
    end if
end subroutine dx

!> Baseline implementation of double layer operator
subroutine dx_dense(params, constants, workspace, do_diag, x, y, ddx_error)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    real(dp) :: c(3), vij(3), sij(3)
    real(dp) :: vvij, tij, tt, f, rho, ctheta, stheta, cphi, sphi
    integer :: its, isph, jsph, l, m, ind
    real(dp), external :: dnrm2

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

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
        call ddintegrate(1, constants % nbasis, params % ngrid, &
            & constants % vwgrid, constants % vgrid_nbasis, &
            & workspace % tmp_grid, y(:,isph))
    end do
end subroutine dx_dense

!> FMM-accelerated implementation of double layer operator
subroutine dx_fmm(params, constants, workspace, do_diag, x, y, ddx_error)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: isph, inode, l, indl, indl1

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

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
        & workspace % tmp_grid, workspace % tmp_sph_l)
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

!> Adjoint double layer operator matvec 
subroutine dstarx(params, constants, workspace, x, y, ddx_error)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    !! Local variables
    integer :: do_diag
    !! initialize do_diag
    do_diag = 0
    if (constants % dodiag) do_diag = 1
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dstarx_dense(params, constants, workspace, do_diag, x, y, ddx_error)
    else
        call dstarx_fmm(params, constants, workspace, do_diag, x, y, ddx_error)
    end if
end subroutine dstarx

!> Baseline implementation of adjoint double layer operator
subroutine dstarx_dense(params, constants, workspace, do_diag, x, y, ddx_error)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_error_type), intent(inout) :: ddx_error
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    real(dp) :: vji(3), sji(3)
    real(dp) :: vvji, tji, tt, f, rho, ctheta, stheta, cphi, sphi
    integer :: its, isph, jsph, l, m, ind, iproc
    real(dp), external :: dnrm2

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0)continue

    y = zero
    !$omp parallel do default(none) shared(do_diag,params,constants, &
    !$omp workspace,x,y) private(isph,jsph,its,vji,vvji,tji,sji,rho, &
    !$omp ctheta,stheta,cphi,sphi,tt,l,ind,f,m,iproc) schedule(dynamic)
    do isph = 1, params % nsph
        iproc = omp_get_thread_num() + 1
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
                            & workspace % tmp_vylm(:(params % lmax + 1)**2, iproc), &
                            & workspace % tmp_vplm(:(params % lmax + 1)**2, iproc), &
                            & workspace % tmp_vcos(:(params % lmax + 1), iproc), &
                            & workspace % tmp_vsin(:(params % lmax + 1), iproc))
                        tt = constants % ui(its,jsph) &
                            & *dot_product(constants % vwgrid(:constants % nbasis, its), &
                            & x(:,jsph))/tji
                        do l = 0, params % lmax
                            ind = l*l + l + 1
                            f = dble(l)*tt/ constants % vscales(ind)**2
                            do m = -l, l
                                y(ind+m,isph) = y(ind+m,isph) + &
                                    & f*workspace % tmp_vylm(ind+m, iproc)
                            end do
                            tt = tt/tji
                        end do
                    end if
                end do
            else if (do_diag .eq. 1) then
                do its = 1, params % ngrid
                    f = pt5*constants % ui(its,jsph) &
                        & *dot_product(constants % vwgrid(:constants % nbasis, its), &
                        & x(:,jsph))
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
subroutine dstarx_fmm(params, constants, workspace, do_diag, x, y, ddx_error)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    integer, intent(in) :: do_diag
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_error_type), intent(inout) :: ddx_error
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    integer :: isph, inode, l, indl, indl1

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

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
        & workspace % tmp_node_l, workspace % tmp_sph_l)
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

!> Apply \f$ R_\varepsilon \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\varepsilon x = (2\pi(\varepsilon + 1) / (\varepsilon
!! - 1) - D) x \f$.
subroutine repsx(params, constants, workspace, x, y, ddx_error)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    !! Local variables
    real(dp) :: fac
    !! Output `y` is cleaned here
    call dx(params, constants, workspace, x, y, ddx_error)
    y = - y
    if (constants % dodiag) then
    ! Apply diagonal
        fac = twopi * (params % eps + one) / (params % eps - one)
        y = fac*x + y
    end if
end subroutine repsx

!> Apply \f$ R_\varepsilon^* \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R^*_\varepsilon x = (2\pi(\varepsilon + 1) / (\varepsilon
!! - 1) - D) x \f$.
!!
subroutine repsstarx(params, constants, workspace, x, y, ddx_error)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! Local variables
    real(dp) :: fac
    ! Output `y` is cleaned here
    call dstarx(params, constants, workspace, x, y, ddx_error)
    y = - y
    ! Apply diagonal
    if (constants % dodiag) then 
        fac = twopi * (params % eps + one) / (params % eps - one)
        y = fac*x + y
    end if
end subroutine repsstarx

!> Apply \f$ R_\infty \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\infty x = (2\pi - D) x \f$.
!!
!! @param[in] ddx_data:
!! @param[in] x:
!! @param[out] y:
!! @param[inout] ddx_error: ddX ddx_error
!!
subroutine rinfx(params, constants, workspace, x, y, ddx_error)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    !! note do_diag hardcoded to 1.
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dx_dense(params, constants, workspace, 1, x, y, ddx_error)
    else
        call dx_fmm(params, constants, workspace, 1, x, y, ddx_error)
    end if
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
!! @param[inout] ddx_error: ddX ddx_error
!!
subroutine rstarinfx(params, constants, workspace, x, y, ddx_error)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! Output `y` is cleaned here
    call dstarx(params, constants, workspace, x, y, ddx_error)
    ! Apply diagonal
    y = twopi*x - y
end subroutine rstarinfx

!> Apply preconditioner for ddPCM primal linear system
!! @param[inout] ddx_error: ddX ddx_error
subroutine prec_repsx(params, constants, workspace, x, y, ddx_error)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    integer :: isph

    ! dummy operation on unused interface arguments
    if ((ddx_error % flag .eq. 0) .or. &
   &    (allocated(workspace % tmp_pot))) continue

    ! simply do a matrix-vector product with the transposed preconditioner 
    !$omp parallel do default(shared) schedule(static,1) &
    !$omp private(isph)
    do isph = 1, params % nsph
        call dgemm('N', 'N', constants % nbasis, 1, constants % nbasis, one, &
            & constants % rx_prc(:, :, isph), constants % nbasis, x(:, isph), &
            & constants % nbasis, zero, y(:, isph), constants % nbasis)
    end do
end subroutine prec_repsx

!> Apply preconditioner for ddPCM adjoint linear system
subroutine prec_repsstarx(params, constants, workspace, x, y, ddx_error)
    implicit none
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    ! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    ! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! Local variables
    integer :: isph

    ! dummy operation on unused interface arguments
    if ((ddx_error % flag .eq. 0) .or. &
   &    (allocated(workspace % tmp_pot))) continue

    ! simply do a matrix-vector product with the transposed preconditioner 
    !$omp parallel do default(shared) schedule(static,1) &
    !$omp private(isph)
    do isph = 1, params % nsph
        call dgemm('T', 'N', constants % nbasis, 1, constants % nbasis, one, &
            & constants % rx_prc(:, :, isph), constants % nbasis, x(:, isph), &
            & constants % nbasis, zero, y(:, isph), constants % nbasis)
    end do
end subroutine prec_repsstarx

!> Adjoint HSP matrix vector product
subroutine bstarx(params, constants, workspace, x, y, ddx_error)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error
    ! Local variables
    integer :: isph, jsph, ij, indmat, ind, ijigrid, igrid, l, m, &
        & nsph, nbasis, lmax, ngrid
    real(dp) :: thigh, se, eta, kappa, vji(3), vvji, sji(3), tji, &
        & rho, ctheta, stheta, cphi, sphi, xji, oji, fac, &
        & basloc((params % lmax + 1)**2), vcos(params % lmax + 1), &
        & vsin(params % lmax + 1), vplm((params % lmax + 1)**2), &
        & tmp(constants % nbasis), fac_hsp(constants % nbasis), &
        & si_rjin(0:params % lmax), di_rjin(0:params % lmax) 
    complex(dp) :: work_complex(max(2, params % lmax+1))

    call time_push()

    nsph = params % nsph
    nbasis = constants % nbasis

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    if (params % matvecmem .eq. 1) then
        associate(inl => constants % inl, nl => constants % nl, &
                & itrnl => constants % itrnl, b => constants % b)
            !$omp parallel do default(none) schedule(dynamic,1) &
            !$omp firstprivate(nsph,nbasis) shared(x,y,inl,nl,itrnl,b) &
            !$omp private(isph,ij,jsph,indmat)
            do isph = 1, nsph
                y(:,isph) = zero
                do ij = inl(isph), inl(isph + 1) - 1
                    jsph = nl(ij)
                    indmat = itrnl(ij)
                    call dgemv('t', nbasis, nbasis, 1.0d0, &
                        & b(:,:,indmat), nbasis, x(:,jsph), 1, &
                        & 1.0d0, y(:,isph), 1)
                end do
            end do
        end associate
    else
        ngrid = params % ngrid
        thigh = one + (params % se+one)/two*params % eta
        se = params % se
        eta = params % eta
        lmax = params % lmax
        kappa = params % kappa

        associate(vgrid => constants % vgrid, &
                & tmp_grid => workspace % tmp_grid)
            !$omp parallel do collapse(2) default(none) &
            !$omp firstprivate(nsph,ngrid,nbasis) shared(x,tmp_grid,vgrid) &
            !$omp private(isph,igrid) schedule(static,100)
            do isph = 1, nsph
                do igrid = 1, ngrid
                    tmp_grid(igrid,isph) = dot_product(x(:, isph), &
                        & vgrid(:nbasis, igrid))
                end do
            end do
        end associate

        associate(inl => constants % inl, nl => constants % nl, &
                & csph => params % csph, rsph => params % rsph, &
                & cgrid => constants % cgrid, fi => constants % fi, &
                & wgrid => constants % wgrid, &
                & tmp_grid => workspace % tmp_grid, &
                & vscales => constants % vscales, &
                & adj_ioverlap => constants % adj_ioverlap, &
                & adj_overlap => constants % adj_overlap, &
                & si_ri => constants % si_ri)
            !$omp parallel do default(none) schedule(dynamic,1) &
            !$omp firstprivate(nsph,ngrid,thigh,se,eta,lmax,kappa) &
            !$omp shared(inl,nl,csph,rsph,cgrid,fi,wgrid,tmp_grid,y, &
            !$omp adj_overlap,adj_ioverlap,vscales,si_ri) &
            !$omp private(isph,ij,jsph,ijigrid,igrid,vji,vvji,tji,sji, &
            !$omp rho,ctheta,stheta,cphi,sphi,basloc,vplm,vcos,vsin,tmp, &
            !$omp si_rjin,di_rjin,work_complex,fac_hsp,xji,oji,fac,ind)
            do isph = 1, nsph
                tmp = 0.0d0
                do ij = inl(isph), inl(isph+1)-1
                    jsph = nl(ij)
                    do ijigrid = adj_ioverlap(ij), adj_ioverlap(ij+1) - 1
                        igrid = adj_overlap(ijigrid)

                        vji = csph(:,jsph) + rsph(jsph)*cgrid(:,igrid) &
                            & - csph(:,isph)
                        vvji = sqrt(vji(1)*vji(1) + vji(2)*vji(2) &
                            & + vji(3)*vji(3))
                        tji = vvji/rsph(isph)
                        sji = vji/vvji

                        call ylmbas(sji, rho, ctheta, stheta, cphi, sphi, &
                            & lmax, vscales, basloc, vplm, vcos, vsin)

                        si_rjin = 0.0d0
                        di_rjin = 0.0d0

                        call modified_spherical_bessel_first_kind(lmax, &
                            & vvji*kappa, si_rjin, di_rjin, work_complex)

                        do l = 0, lmax
                            ind = l*l + l + 1
                            do  m = -l, l
                                fac_hsp(ind+m) = SI_rjin(l)/SI_ri(l,isph)*basloc(ind+m)
                            end do
                        end do

                        xji = fsw(tji, se, eta)
                        if (fi(igrid,jsph).gt.1.0d0) then
                            oji = xji/fi(igrid,jsph)
                        else
                            oji = xji
                        end if

                        fac = wgrid(igrid) * tmp_grid(igrid,jsph) * oji
                        do l = 0, lmax
                            ind = l*l + l + 1
                            do m = -l, l
                                tmp(ind+m) = tmp(ind+m) + fac*fac_hsp(ind+m)
                            end do
                        end do
                    end do
                end do
                y(:, isph) = - tmp
            end do
        end associate
    end if

    ! add the diagonal if required
    if (constants % dodiag) y = y + x
    call time_pull("bstarx")

end subroutine bstarx

!> Primal HSP matrix vector product
subroutine bx(params, constants, workspace, x, y, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    type(ddx_error_type), intent(inout) :: ddx_error
    integer :: isph, jsph, ij, ijits, its

    integer :: nsph, nbasis, lmax, ngrid, vgrid_nbasis
    real(dp) :: se, eta, vij(3), vvij, tij, xij, oij, vtij(3), kappa
    complex(dp) :: work_complex(params % lmax+1)
    real(dp) :: work(params % lmax+1), tmp_grid(params % ngrid)

    call time_push()

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue
    if (workspace % xs_time .eq. 0) continue

    nsph = params % nsph
    nbasis = constants % nbasis
    lmax = params % lmax

    if (params % matvecmem .eq. 1) then
        !$omp parallel do default(none) shared(params,constants,x,y) &
        !$omp private(isph,ij,jsph) schedule(dynamic)
        do isph = 1, params % nsph
            y(:, isph) = 0.0d0
            do ij = constants % inl(isph), constants % inl(isph + 1) - 1
                jsph = constants % nl(ij)
                call dgemv('n', constants % nbasis, constants % nbasis, one, &
                    & constants % b(:,:,ij), constants % nbasis, x(:,jsph), 1, &
                    & one, y(:,isph), 1)
            end do
        end do
    else
        se = params % se
        eta = params % eta
        ngrid = params % ngrid
        kappa = params % kappa
        vgrid_nbasis = constants % vgrid_nbasis

        associate(inl => constants % inl, nl => constants % nl, &
                & csph => params % csph, rsph => params % rsph, &
                & cgrid => constants % cgrid, fi => constants % fi, &
                & vscales => constants % vscales, &
                & vwgrid => constants % vwgrid, &
                & ioverlap => constants % ioverlap, &
                & overlap => constants % overlap, &
                & si_ri => constants % si_ri)
            !$omp parallel do default(none) schedule(dynamic,1) &
            !$omp firstprivate(nsph,se,eta,lmax,nbasis,ngrid, &
            !$omp vgrid_nbasis,kappa) shared(inl,nl,csph,rsph,cgrid,fi, &
            !$omp vscales,x,y,vwgrid,ioverlap,overlap,si_ri) private( &
            !$omp isph,tmp_grid,its,ij,jsph,vij,vvij,tij,xij,oij,vtij, &
            !$omp work,work_complex,ijits)
            do isph = 1, nsph
                tmp_grid(:) = 0.0d0
                do ij = inl(isph), inl(isph+1)-1
                    jsph = nl(ij)
                    do ijits = ioverlap(ij), ioverlap(ij+1) - 1
                        its = overlap(ijits)
                        vij = csph(:,isph) + rsph(isph)*cgrid(:,its) &
                            & - csph(:,jsph)
                        vvij = sqrt(vij(1)*vij(1) + vij(2)*vij(2) &
                            & + vij(3)*vij(3))
                        tij  = vvij/rsph(jsph)
                        xij = fsw(tij, se, eta)
                        if (fi(its,isph).gt.1.0d0) then
                            oij = xij/fi(its, isph)
                        else
                            oij = xij
                        end if
                        vtij = vij*kappa
                        call fmm_l2p_bessel_work(vtij, lmax, vscales, &
                            & si_ri(:, jsph), oij, x(:, jsph), 1.0d0, &
                            & tmp_grid(its), work_complex, work)
                    end do
                end do
                call ddintegrate(1, nbasis, ngrid, vwgrid, vgrid_nbasis, &
                    & tmp_grid, y(:, isph))
                y(:, isph) = - y(:, isph)
            end do
        end associate
    end if
    if (constants % dodiag) y = y + x
    call time_pull("bx")
end subroutine bx

!> Adjoint ddLPB matrix-vector product
subroutine tstarx(params, constants, workspace, x, y, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y
    type(ddx_error_type), intent(inout) :: ddx_error

    real(dp), dimension(constants % nbasis, params % nsph, 2) :: temp_vector
    y = zero
    temp_vector = zero

    ! Compute AXr
    ! call LstarXr
    call lstarx(params, constants, workspace, x(:,:,1), temp_vector(:,:,1), &
        & ddx_error)
    ! Remove the scaling factor
    call convert_ddcosmo(params, constants, 1, temp_vector(:,:,1))
    y(:,:,1) = y(:,:,1) + temp_vector(:,:,1)
    ! Compute BXe
    call bstarx(params, constants, workspace, x(:,:,2), temp_vector(:,:,2), &
        & ddx_error)
    y(:,:,2) = y(:,:,2) + temp_vector(:,:,2)

    ! Reset temp_vector to zero
    temp_vector = zero
    ! Call CX
    call cstarx(params, constants, workspace, x, temp_vector, ddx_error)
    y = y + temp_vector
end subroutine tstarx

!> Apply the preconditioner to the primal HSP linear system
!! @param[inout] ddx_error: ddX ddx_error
subroutine bx_prec(params, constants, workspace, x, y, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    type(ddx_error_type), intent(inout) :: ddx_error

    ! dummy operation on unused interface arguments
    if ((ddx_error % flag .eq. 0) .or. &
   &    (allocated(workspace % tmp_pot))) continue

    y = x
end subroutine bx_prec

!> Apply the preconditioner to the ddLPB adjoint linear system
subroutine prec_tstarx(params, constants, workspace, x, y, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph, 2)
    real(dp), intent(inout) :: y(constants % nbasis, params % nsph, 2)
    type(ddx_error_type), intent(inout) :: ddx_error
    integer :: n_iter
    real(dp), dimension(params % maxiter) :: x_rel_diff
    real(dp) :: start_time

    start_time = omp_get_wtime()
    y(:,:,1) = x(:,:,1)
    call convert_ddcosmo(params, constants, 1, y(:,:,1))
    n_iter = params % maxiter
    ! the empirical 10^-2 factor reduces the number of macro iterations
    call jacobi_diis(params, constants, workspace, &
        & constants % inner_tol*1.0d-2, y(:,:,1), &
        & workspace % ddcosmo_guess, n_iter, x_rel_diff, lstarx, &
        & ldm1x, hnorm, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, 'prec_tstarx: ddCOSMO failed to ' // &
            & 'converge, exiting')
        return
    end if
    y(:,:,1) = workspace % ddcosmo_guess
    workspace % s_time = workspace % s_time + omp_get_wtime() - start_time

    start_time = omp_get_wtime()
    n_iter = params % maxiter
    call jacobi_diis(params, constants, workspace, &
        & constants % inner_tol*1.0d1, x(:,:,2), workspace % hsp_guess, &
        & n_iter, x_rel_diff, bstarx, bx_prec, hnorm, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, 'prec_tstarx: HSP failed to ' // &
            & 'converge, exiting')
        return
    end if
    y(:,:,2) = workspace % hsp_guess
    workspace % hsp_adj_time = workspace % hsp_adj_time &
        & + omp_get_wtime() - start_time

end subroutine prec_tstarx

!> Apply the preconditioner to the primal ddLPB linear system
!! |Yr| = |A^-1 0 |*|Xr|
!! |Ye|   |0 B^-1 | |Xe|
!! @param[in] ddx_data : dd Data
!! @param[in] x        : Input array
!! @param[out] y       : Linear system solution at current iteration
!! @param[inout] ddx_error: ddX ddx_error
subroutine prec_tx(params, constants, workspace, x, y, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph, 2)
    real(dp), intent(inout) :: y(constants % nbasis, params % nsph, 2)
    type(ddx_error_type), intent(inout) :: ddx_error
    integer :: n_iter
    real(dp), dimension(params % maxiter) :: x_rel_diff
    real(dp) :: start_time

    ! perform A^-1 * Yr
    start_time = omp_get_wtime()
    n_iter = params % maxiter
    call jacobi_diis(params, constants, workspace, &
        & constants % inner_tol, x(:,:,1), workspace % ddcosmo_guess, &
        & n_iter, x_rel_diff, lx, ldm1x, hnorm, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, 'prec_tx: ddCOSMO failed to ' // &
            & 'converge, exiting')
        return
    end if

    ! Scale by the factor of (2l+1)/4Pi
    y(:,:,1) = workspace % ddcosmo_guess
    call convert_ddcosmo(params, constants, 1, y(:,:,1))
    workspace % xs_time = workspace % xs_time + omp_get_wtime() - start_time

    ! perform B^-1 * Ye
    start_time = omp_get_wtime()
    n_iter = params % maxiter
    call jacobi_diis(params, constants, workspace, &
        & constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
        & n_iter, x_rel_diff, bx, bx_prec, hnorm, ddx_error)
    y(:,:,2) = workspace % hsp_guess
    workspace % hsp_time = workspace % hsp_time + omp_get_wtime() - start_time

    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, 'prec_tx: HSP failed to ' // &
            & 'converge, exiting')
        return
    end if

end subroutine prec_tx

!> ddLPB adjoint matrix-vector product
subroutine cstarx(params, constants, workspace, x, y, ddx_error)
    implicit none
    ! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y
    type(ddx_error_type), intent(inout) :: ddx_error
    ! local
    complex(dp) :: work_complex(constants % lmax0+1)
    real(dp) :: work(constants % lmax0+1)
    integer :: isph, igrid, jsph, ind, l0, ind0, indl, inode, l, m
    real(dp), dimension(3) :: vij, vtij
    real(dp) :: val, epsilon_ratio
    real(dp), allocatable :: scratch(:,:), scratch0(:,:)

    call time_push()

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    allocate(scratch(constants % nbasis, params % nsph), &
        & scratch0(constants % nbasis0, params % nsph))

    epsilon_ratio = params % epsp/params % eps

    ! TODO: maybe use ddeval_grid for code consistency
    scratch = - x(:,:,1) - x(:,:,2)
    call dgemm('T', 'N', params % ngrid, params % nsph, constants % nbasis, &
        & one, constants % vwgrid, constants % vgrid_nbasis, scratch, &
        & constants % nbasis, zero, workspace % tmp_grid, params % ngrid)

    scratch0 = zero
    if (params % fmm .eq. 0) then
        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                if (constants % ui(igrid, isph).gt.zero) then
                    val = workspace % tmp_grid(igrid,isph) &
                        & *constants % ui(igrid,isph)
                    ! quadratically scaling loop
                    do jsph = 1, params % nsph
                        vij  = params % csph(:,isph) + &
                            & params % rsph(isph)*constants % cgrid(:,igrid) - &
                            & params % csph(:,jsph)
                        vtij = vij * params % kappa
                        call fmm_m2p_bessel_adj_work(vtij, val, &
                            & constants % SK_ri(:, jsph), constants % lmax0, &
                            & constants % vscales, one, scratch0(:, jsph), &
                            & work_complex, work)
                    end do
                end if
            end do
        end do
    else
        ! Multiply by characteristic function
        workspace % tmp_grid = workspace % tmp_grid * constants % ui
        workspace % tmp_sph = zero
        ! Do FMM operations adjointly
        call tree_m2p_bessel_adj(params, constants, constants % lmax0, &
            & one, workspace % tmp_grid, zero, params % lmax, &
            & workspace % tmp_sph)
        call tree_l2p_bessel_adj(params, constants, one, &
            & workspace % tmp_grid, zero, workspace % tmp_node_l)
        call tree_l2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l)
        call tree_m2l_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_l, workspace % tmp_node_m)
        call tree_m2m_bessel_rotation_adj(params, constants, &
            & workspace % tmp_node_m)
        ! Adjointly move tree multipole harmonics into output
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                scratch0(:, isph) = workspace % tmp_sph(1:constants % nbasis0, isph) + &
                    & workspace % tmp_node_m(1:constants % nbasis0, inode)
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                scratch0(1:indl, isph) = &
                    & workspace % tmp_sph(1:indl, isph) + &
                    & workspace % tmp_node_m(:, inode)
                scratch0(indl+1:, isph) = zero
            end do
        end if
    end if
    ! Scale by C_ik
    do isph = 1, params % nsph
        do l0 = 0, constants % lmax0
            ind0 = l0*l0 + l0 + 1
            scratch0(ind0-l0:ind0+l0, isph) = scratch0(ind0-l0:ind0+l0, isph) * &
                & constants % C_ik(l0, isph)
        end do
    end do

    do jsph = 1, params % nsph
        call dgemv('n', constants % nbasis, constants % nbasis0, one, &
            & constants % pchi(1,1,jsph), constants % nbasis, &
            & scratch0(1,jsph), 1, zero, scratch(1,jsph), 1)
    end do

    do isph = 1, params % nsph
        do l = 0, params % lmax
            do m = -l, l
                ind = l**2 + l + m + 1
                y(ind,isph,1) = - (epsilon_ratio*l*scratch(ind,isph)) &
                    & /params % rsph(isph)
                y(ind,isph,2) = constants % termimat(l,isph)*scratch(ind,isph)
          end do
        end do
    end do

    call time_pull("cstarx")

end subroutine cstarx

!> ddLPB matrix-vector product
subroutine cx(params, constants, workspace, x, y, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y
    type(ddx_error_type), intent(inout) :: ddx_error

    integer :: isph, jsph, igrid, ind, l, m, indl, inode, info, &
        & nsph, lmax, nbasis, nbasis0, lmax0, ngrid
    real(dp) :: vij(3), vtij(3), val, epsp, eps, fac, &
        & work(constants % lmax0 + 1)
    complex(dp) :: work_complex(constants % lmax0 + 1)
    real(dp), allocatable :: diff_re(:,:), diff0(:,:)

    call time_push()
    allocate(diff_re(constants % nbasis, params % nsph), &
        & diff0(constants % nbasis0, params % nsph), stat=info)
    if (info.ne.0) then
        call update_error(ddx_error, "Allocation failed in cx")
        return
    end if

    nsph = params % nsph
    lmax = params % lmax
    lmax0 = constants % lmax0
    epsp = params % epsp
    eps = params % eps
    nbasis = constants % nbasis
    nbasis0 = constants % nbasis0
    ngrid = params % ngrid

    call time_push()
    associate(rsph => params % rsph, termimat => constants % termimat)
        !$omp parallel do collapse(2) default(none) schedule(static,100) &
        !$omp firstprivate(nsph,lmax,epsp,eps) private(jsph,l,ind,m) &
        !$omp shared(diff_re,rsph,x,termimat)
        do jsph = 1, nsph
            do l = 0, lmax
                ind = l**2 + l + 1
                do m = -l, l
                    diff_re(ind+m,jsph) = (epsp/eps)*(l/rsph(jsph)) &
                        & *x(ind+m,jsph,1) - termimat(l,jsph)*x(ind+m,jsph,2)
                end do
            end do
        end do
    end associate
    call time_pull("cx1")

    ! diff0 = Pchi * diff_er, linear scaling
    call time_push()
    associate(pchi => constants % pchi)
        !$omp parallel do default(none) schedule(static,10) &
        !$omp firstprivate(nsph,nbasis,nbasis0) private(jsph) &
        !$omp shared(pchi,diff_re,diff0)
        do jsph = 1, nsph
            call dgemv('t', nbasis, nbasis0, 1.0d0, pchi(1,1,jsph), nbasis, &
                & diff_re(1,jsph), 1, 0.0d0, diff0(1,jsph), 1)
        end do
    end associate
    call time_pull("cx2")

    ! Multiply diff0 by C_ik inplace
    call time_push()
    associate(c_ik => constants % c_ik)
        !$omp parallel do collapse(2) schedule(static,100) &
        !$omp firstprivate(nsph,lmax0) private(isph,l,ind,fac,m) &
        !$omp shared(diff0,c_ik)
        do isph = 1, nsph
            do l = 0, lmax0
                ind = l*l+l+1
                fac = c_ik(l, isph)
                do m = -l, l
                    diff0(ind+m, isph) = diff0(ind+m, isph) * fac
                end do
            end do
        end do
    end associate
    call time_pull("cx3")

    ! avoiding N^2 storage, this code does not use the cached coefY
    if (params % fmm .eq. 0) then
        do isph = 1, params % nsph
            y(:,isph,1) = zero
            do igrid = 1, params % ngrid
                if (constants % ui(igrid,isph).gt.zero) then
                    val = zero

                    ! quadratically scaling loop
                    do jsph = 1, params % nsph
                        vij  = params % csph(:,isph) + &
                            & params % rsph(isph)*constants % cgrid(:,igrid) - &
                            & params % csph(:,jsph)
                        vtij = vij * params % kappa
                        call fmm_m2p_bessel_work(vtij, constants % lmax0, &
                            & constants % vscales, constants % SK_ri(:, jsph), &
                            & one, diff0(:, jsph), one, val, work_complex, work)
                    end do
                    do ind = 1, constants % nbasis
                        y(ind,isph,1) = y(ind,isph,1) + val*&
                          & constants % ui(igrid,isph)*&
                          & constants % vwgrid(ind,igrid)
                    end do
                end if
            end do
        end do
    else
        ! Load input harmonics into tree data
        call time_push()
        workspace % tmp_node_m = zero
        workspace % tmp_node_l = zero
        workspace % tmp_sph = zero
        associate(tmp_sph => workspace % tmp_sph)
            !$omp parallel do collapse(2) schedule(static,100) &
            !$omp firstprivate(nsph,lmax0), shared(tmp_sph,diff0) &
            !$omp private(isph,l,ind)
            do isph = 1, nsph
                do l = 0, lmax0
                    ind = l*l+l+1
                    tmp_sph(ind-l:ind+l, isph) = diff0(ind-l:ind+l, isph)
                end do
            end do
        end associate
        if(constants % lmax0 .lt. params % pm) then
            associate(tmp_node_m => workspace % tmp_node_m, &
                    & snode => constants % snode)
                !$omp parallel do schedule(static,1) default(none) &
                !$omp firstprivate(nsph,nbasis0) &
                !$omp shared(snode,tmp_node_m,diff0) &
                !$omp private(isph,inode)
                do isph = 1, nsph
                    inode = snode(isph)
                    tmp_node_m(1:nbasis0, inode) = &
                        & diff0(1:nbasis0, isph)
                    tmp_node_m(nbasis0+1:, inode) = 0.0d0
                end do
            end associate
        else
            indl = (params % pm+1)**2
            associate(tmp_node_m => workspace % tmp_node_m, &
                    & snode => constants % snode)
                !$omp parallel do schedule(static,1) default(none) &
                !$omp firstprivate(nsph,nbasis0) &
                !$omp shared(snode,tmp_node_m,diff0) &
                !$omp private(isph,inode,indl)
                do isph = 1, nsph
                    inode = snode(isph)
                    tmp_node_m(:, inode) = diff0(1:indl, isph)
                end do
            end associate
        end if
        call time_pull("cx-fmmprep")
        call time_push()
        ! Do FMM operations
        call tree_m2m_bessel_rotation(params, constants, workspace % tmp_node_m)
        call time_pull("cx-m2m")
        call time_push()
        call tree_m2l_bessel_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call time_pull("cx-m2l")
        call time_push()
        call tree_l2l_bessel_rotation(params, constants, workspace % tmp_node_l)
        call time_pull("cx-l2l")
        call time_push()
        workspace % tmp_grid = zero
        call tree_l2p_bessel(params, constants, one, workspace % tmp_node_l, zero, &
            & workspace % tmp_grid)
        call time_pull("cx-l2p")
        call time_push()
        call tree_m2p_bessel(params, constants, constants % lmax0, one, &
            & params % lmax, workspace % tmp_sph, one, &
            & workspace % tmp_grid)
        call time_pull("cx-m2p")
        call time_push()

        associate(tmp_grid => workspace % tmp_grid, &
                 & vwgrid => constants % vwgrid, ui => constants % ui)
            !$omp parallel do collapse(2) schedule(static,100) &
            !$omp firstprivate(nsph,ngrid,nbasis), shared(y,tmp_grid, &
            !$omp vwgrid,ui) private(isph,igrid,ind,val)
            do isph = 1, nsph
                do ind = 1, nbasis
                    val = 0.0d0
                    do igrid = 1, ngrid
                         val = val + tmp_grid(igrid, isph)*&
                            & vwgrid(ind, igrid)*ui(igrid,isph)
                    end do
                    y(ind,isph,1) = val
                end do
            end do
        end associate
        call time_pull("cx-fmmend")
    end if

    y(:,:,2) = y(:,:,1)

    call time_pull("cx")
    deallocate(diff_re, diff0, stat=info)
    if (info.ne.0) then
        call update_error(ddx_error, "Deallocation failed in cx")
        return
    end if

end subroutine cx

end module ddx_operators
