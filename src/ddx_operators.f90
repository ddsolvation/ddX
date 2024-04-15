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
    integer :: isph, jsph, ij, l, ind, iproc

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    !! Initialize
    y = zero
!
!   incore code:
!
    if (params % matvecmem .eq. 1) then
        !$omp parallel do default(none) shared(params,constants,x,y) &
        !$omp private(isph,ij,jsph) schedule(dynamic)
        do isph = 1, params % nsph
            do ij = constants % inl(isph), constants % inl(isph + 1) - 1
                jsph = constants % nl(ij)
                call dgemv('n', constants % nbasis, constants % nbasis, one, &
                    & constants % l(:,:,ij), constants % nbasis, x(:,jsph), 1, &
                    & one, y(:,isph), 1)
            end do
        end do
    else
        call time_push()
        !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
        !$omp private(isph,iproc) schedule(dynamic)
        do isph = 1, params % nsph
            iproc = omp_get_thread_num() + 1
            ! Compute NEGATIVE action of off-digonal blocks
            call calcv(params, constants, isph, workspace % tmp_pot(:, iproc), x, &
                & workspace % tmp_work(:, iproc))
            call ddintegrate(1, constants % nbasis, params % ngrid, &
                & constants % vwgrid, constants % vgrid_nbasis, &
                & workspace % tmp_pot(:, iproc), y(:, isph))
            ! now, fix the sign.
            y(:, isph) = - y(:, isph)
        end do
        call time_pull("lx")
    end if
!
!   if required, add the diagonal.
!
    if (constants % dodiag) then 
        do isph = 1, params % nsph
            do l = 0, params % lmax
                ind = l*l + l + 1
                y(ind-l:ind+l, isph) = y(ind-l:ind+l, isph) + &
                     x(ind-l:ind+l, isph) / (constants % vscales(ind)**2)
            end do
        end do
    end if
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
    integer :: isph, jsph, ij, indmat, igrid, l, ind, iproc

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    y = zero
    if (params % matvecmem .eq. 1) then
        !$omp parallel do default(none) shared(params,constants,x,y) &
        !$omp private(isph,ij,jsph,indmat) schedule(dynamic)
        do isph = 1, params % nsph
            do ij = constants % inl(isph), constants % inl(isph + 1) - 1
                jsph = constants % nl(ij)
                indmat = constants % itrnl(ij)
                call dgemv('t', constants % nbasis, constants % nbasis, one, &
                    & constants % l(:,:,indmat), constants % nbasis, x(:,jsph), 1, &
                    & one, y(:,isph), 1)
            end do
        end do
    else
        ! Expand x over spherical harmonics
        !$omp parallel do default(none) shared(params,constants,workspace,x) &
        !$omp private(isph,igrid) schedule(dynamic)
        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                workspace % tmp_grid(igrid, isph) = dot_product(x(:, isph), &
                    & constants % vgrid(:constants % nbasis, igrid))
            end do
        end do
        ! Compute (negative) action
        !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
        !$omp private(isph,iproc) schedule(dynamic)
        do isph = 1, params % nsph
            iproc = omp_get_thread_num() + 1
            call adjrhs(params, constants, isph, workspace % tmp_grid, &
                & y(:, isph), workspace % tmp_work(:, iproc))
            ! fix the sign
            y(:, isph) = - y(:, isph)
        end do
    end if
    if (constants % dodiag) then
    ! Loop over harmonics
        do isph = 1, params % nsph
            do l = 0, params % lmax
                ind = l*l + l + 1
                y(ind-l:ind+l, isph) = y(ind-l:ind+l, isph) + &
                    & x(ind-l:ind+l, isph) / (constants % vscales(ind)**2)
            end do
        end do
    end if
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
    integer :: isph, l, ind

    ! dummy operation on unused interface arguments
    if ((ddx_error % flag .eq. 0) .or. &
   &    (allocated(workspace % tmp_pot))) continue

    !! Loop over harmonics
    !$omp parallel do default(none) shared(params,constants,x,y) &
    !$omp private(isph,l,ind) schedule(dynamic)
    do isph = 1, params % nsph
        do l = 0, params % lmax
            ind = l*l + l + 1
            y(ind-l:ind+l, isph) = x(ind-l:ind+l, isph) &
                & *(constants % vscales(ind)**2)
        end do
    end do
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
    integer :: isph, jsph, ij, indmat, iproc

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    y = zero
    if (params % matvecmem .eq. 1) then
        !$omp parallel do default(none) shared(params,constants,x,y) &
        !$omp private(isph,ij,jsph,indmat) schedule(dynamic)
        do isph = 1, params % nsph
            do ij = constants % inl(isph), constants % inl(isph + 1) - 1
                jsph = constants % nl(ij)
                indmat = constants % itrnl(ij)
                call dgemv('t', constants % nbasis, constants % nbasis, one, &
                    & constants % b(:,:,indmat), constants % nbasis, x(:,jsph), 1, &
                    & one, y(:,isph), 1)
            end do
        end do
    else
        !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
        !$omp private(isph) schedule(static,1)
        do isph = 1, params % nsph
            call dgemv('t', constants % nbasis, params % ngrid, one, constants % vgrid, &
                & constants % vgrid_nbasis, x(:, isph), 1, zero, &
                & workspace % tmp_grid(:, isph), 1)
        end do
        !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
        !$omp private(isph,iproc) schedule(dynamic)
        do isph = 1, params % nsph
            iproc = omp_get_thread_num() + 1
            call adjrhs_lpb(params, constants, isph, workspace % tmp_grid, &
                & y(:, isph), workspace % tmp_vylm(:, iproc), workspace % tmp_vplm(:, iproc), &
                & workspace % tmp_vcos(:, iproc), workspace % tmp_vsin(:, iproc), &
                & workspace % tmp_bessel(:, iproc))
            y(:,isph)  = - y(:,isph)
        end do
    end if

    ! add the diagonal if required
    if (constants % dodiag) y = y + x
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
    integer :: isph, jsph, ij, iproc

    ! dummy operation on unused interface arguments
    if (ddx_error % flag .eq. 0) continue

    y = zero
    if (params % matvecmem .eq. 1) then
        !$omp parallel do default(none) shared(params,constants,x,y) &
        !$omp private(isph,ij,jsph) schedule(dynamic)
        do isph = 1, params % nsph
            do ij = constants % inl(isph), constants % inl(isph + 1) - 1
                jsph = constants % nl(ij)
                call dgemv('n', constants % nbasis, constants % nbasis, one, &
                    & constants % b(:,:,ij), constants % nbasis, x(:,jsph), 1, &
                    & one, y(:,isph), 1)
            end do
        end do
    else
        call time_push()
        !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
        !$omp private(isph,iproc) schedule(dynamic)
        do isph = 1, params % nsph
          iproc = omp_get_thread_num() + 1
          call calcv2_lpb(params, constants, isph, workspace % tmp_pot(:, iproc), x)
          call ddintegrate(1, constants % nbasis, params % ngrid, constants % vwgrid, &
              & constants % vgrid_nbasis, workspace % tmp_pot(:, iproc), y(:,isph))
          y(:,isph) = - y(:,isph) 
        end do
        call time_pull("bx")
    end if
    if (constants % dodiag) y = y + x
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
    call jacobi_diis(params, constants, workspace, constants % inner_tol, &
        & y(:,:,1), workspace % ddcosmo_guess, n_iter, x_rel_diff, lstarx, &
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
    call jacobi_diis(params, constants, workspace, constants % inner_tol, &
        & x(:,:,2), workspace % hsp_guess, n_iter, x_rel_diff, bstarx, &
        & bx_prec, hnorm, ddx_error)
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
    call jacobi_diis(params, constants, workspace, constants % inner_tol, &
        & x(:,:,1), workspace % ddcosmo_guess, n_iter, x_rel_diff, lx, &
        & ldm1x, hnorm, ddx_error)
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
    call jacobi_diis(params, constants, workspace, constants % inner_tol, &
        & x(:,:,2), workspace % hsp_guess, n_iter, x_rel_diff, bx, &
        & bx_prec, hnorm, ddx_error)
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

    integer :: isph, jsph, igrid, ind, l, m, ind0
    real(dp), dimension(3) :: vij, vtij
    real(dp) :: val
    complex(dp) :: work_complex(constants % lmax0 + 1)
    real(dp) :: work(constants % lmax0 + 1)
    integer :: indl, inode, info

    real(dp), allocatable :: diff_re(:,:), diff0(:,:)

    allocate(diff_re(constants % nbasis, params % nsph), &
        & diff0(constants % nbasis0, params % nsph), stat=info)
    if (info.ne.0) then
        call update_error(ddx_error, "Allocation failed in cx")
        return
    end if

    ! diff_re = params % epsp/eps*l1/ri*Xr - i'(ri)/i(ri)*Xe,
    diff_re = zero
    !$omp parallel do default(none) shared(params,diff_re, &
    !$omp constants,x) private(jsph,l,m,ind)
    do jsph = 1, params % nsph
      do l = 0, params % lmax
        do m = -l,l
          ind = l**2 + l + m + 1
          diff_re(ind,jsph) = (params % epsp/params % eps)* &
              & (l/params % rsph(jsph))*x(ind,jsph,1) &
              & - constants % termimat(l,jsph)*x(ind,jsph,2)
        end do
      end do
    end do

    ! diff0 = Pchi * diff_er, linear scaling
    !$omp parallel do default(none) shared(constants,params, &
    !$omp diff_re,diff0) private(jsph)
    do jsph = 1, params % nsph
      call dgemv('t', constants % nbasis, constants % nbasis0, one, &
          & constants % pchi(1,1,jsph), constants % nbasis, &
          & diff_re(1,jsph), 1, zero, diff0(1,jsph), 1)
    end do

    ! Multiply diff0 by C_ik inplace
    do isph = 1, params % nsph
        do l = 0, constants % lmax0
            ind0 = l*l+l+1
            diff0(ind0-l:ind0+l, isph) = diff0(ind0-l:ind0+l, isph) * &
                & constants % C_ik(l, isph)
        end do
    end do
    ! avoiding N^2 storage, this code does not use the cached coefY
    y(:,:,1) = zero
    if (params % fmm .eq. 0) then
        !$omp parallel do default(none) shared(params,constants, &
        !$omp diff0,y) private(isph,igrid,val,vij,work, &
        !$omp ind0,ind,vtij,work_complex)
        do isph = 1, params % nsph
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
        workspace % tmp_node_m = zero
        workspace % tmp_node_l = zero
        workspace % tmp_sph = zero
        do isph = 1, params % nsph
            do l = 0, constants % lmax0
                ind0 = l*l+l+1
                workspace % tmp_sph(ind0-l:ind0+l, isph) = &
                    & diff0(ind0-l:ind0+l, isph)
            end do
        end do
        if(constants % lmax0 .lt. params % pm) then
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(1:constants % nbasis0, inode) = &
                    & workspace % tmp_sph(1:constants % nbasis0, isph)
                workspace % tmp_node_m(constants % nbasis0+1:, inode) = zero
            end do
        else
            indl = (params % pm+1)**2
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                workspace % tmp_node_m(:, inode) = workspace % tmp_sph(1:indl, isph)
            end do
        end if
        ! Do FMM operations
        call tree_m2m_bessel_rotation(params, constants, workspace % tmp_node_m)
        call tree_m2l_bessel_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call tree_l2l_bessel_rotation(params, constants, workspace % tmp_node_l)
        workspace % tmp_grid = zero
        call tree_l2p_bessel(params, constants, one, workspace % tmp_node_l, zero, &
            & workspace % tmp_grid)
        call tree_m2p_bessel(params, constants, constants % lmax0, one, &
            & params % lmax, workspace % tmp_sph, one, &
            & workspace % tmp_grid)

        do isph = 1, params % nsph
            do igrid = 1, params % ngrid
                do ind = 1, constants % nbasis
                    y(ind,isph,1) = y(ind,isph,1) + workspace % tmp_grid(igrid, isph)*&
                        & constants % vwgrid(ind, igrid)*&
                        & constants % ui(igrid,isph)
                end do
            end do
        end do
    end if

    y(:,:,2) = y(:,:,1)

    deallocate(diff_re, diff0, stat=info)
    if (info.ne.0) then
        call update_error(ddx_error, "Deallocation failed in cx")
        return
    end if

end subroutine cx

end module ddx_operators
