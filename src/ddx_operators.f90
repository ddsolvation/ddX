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
subroutine lx(params, constants, workspace, x, y)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: isph, jsph, ij, l, ind, iproc

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
subroutine lstarx(params, constants, workspace, x, y)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: isph, jsph, ij, indmat, igrid, l, ind, iproc
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
subroutine ldm1x(params, constants, workspace, x, y)
    implicit none
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
subroutine dx(params, constants, workspace, x, y)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local parameter
    integer :: do_diag
    !! initialize do_diag
    do_diag = 0
    if (constants % dodiag) do_diag = 1
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dx_dense(params, constants, workspace, do_diag, x, y)
    else
        call dx_fmm(params, constants, workspace, do_diag, x, y)
    end if
end subroutine dx

!> Baseline implementation of double layer operator
subroutine dx_dense(params, constants, workspace, do_diag, x, y)
    implicit none
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
        call ddintegrate(1, constants % nbasis, params % ngrid, &
            & constants % vwgrid, constants % vgrid_nbasis, &
            & workspace % tmp_grid, y(:,isph))
    end do
end subroutine dx_dense

!> FMM-accelerated implementation of double layer operator
subroutine dx_fmm(params, constants, workspace, do_diag, x, y)
    implicit none
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
subroutine dstarx(params, constants, workspace, x, y)
    implicit none
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    !! Temporaries
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Output
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    !! Local variables
    integer :: do_diag
    !! initialize do_diag
    do_diag = 0
    if (constants % dodiag) do_diag = 1
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dstarx_dense(params, constants, workspace, do_diag, x, y)
    else
        call dstarx_fmm(params, constants, workspace, do_diag, x, y)
    end if
end subroutine dstarx

!> Baseline implementation of adjoint double layer operator
subroutine dstarx_dense(params, constants, workspace, do_diag, x, y)
    implicit none
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
                        tt = constants % ui(its,jsph) &
                            & *dot_product(constants % vwgrid(:constants % nbasis, its), &
                            & x(:,jsph))/tji
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
subroutine dstarx_fmm(params, constants, workspace, do_diag, x, y)
    implicit none
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
subroutine repsx(params, constants, workspace, x, y)
    implicit none
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
subroutine repsstarx(params, constants, workspace, x, y)
    implicit none
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
subroutine rinfx(params, constants, workspace, x, y)
    implicit none
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
    !! note do_diag hardcoded to 1.
    !! Select implementation
    if (params % fmm .eq. 0) then
        call dx_dense(params, constants, workspace, 1, x, y)
    else
        call dx_fmm(params, constants, workspace, 1, x, y)
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
subroutine rstarinfx(params, constants, workspace, x, y)
    implicit none
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

!> Apply preconditioner for ddPCM primal linear system
subroutine prec_repsx(params, constants, workspace, x, y)
    implicit none
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
    !$omp parallel do default(shared) schedule(static,1) &
    !$omp private(isph)
    do isph = 1, params % nsph
        call dgemm('N', 'N', constants % nbasis, 1, constants % nbasis, one, &
            & constants % rx_prc(:, :, isph), constants % nbasis, x(:, isph), &
            & constants % nbasis, zero, y(:, isph), constants % nbasis)
    end do
end subroutine prec_repsx

!> Apply preconditioner for ddPCM adjoint linear system
subroutine prec_repsstarx(params, constants, workspace, x, y)
    implicit none
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
    !$omp parallel do default(shared) schedule(static,1) &
    !$omp private(isph)
    do isph = 1, params % nsph
        call dgemm('T', 'N', constants % nbasis, 1, constants % nbasis, one, &
            & constants % rx_prc(:, :, isph), constants % nbasis, x(:, isph), &
            & constants % nbasis, zero, y(:, isph), constants % nbasis)
    end do
end subroutine prec_repsstarx

!> Gradient of the ddPCM matrix
subroutine gradr(params, constants, workspace, g, ygrid, fx)
    implicit none
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

!> Gradient of the ddPCM matrix using N^2 code
subroutine gradr_dense(params, constants, workspace, g, ygrid, fx)
    implicit none
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

!> Sphere contribution to the ddPCM matrix gradient using N^2 code
subroutine gradr_sph(params, constants, isph, vplm, vcos, vsin, basloc, &
        & dbsloc, g, ygrid, fx)
    implicit none
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

!> Compute the ddPCM matrix gradient using FMMS (2 matvecs)
subroutine gradr_fmm(params, constants, workspace, g, ygrid, fx)
    implicit none
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
        & workspace % tmp_node_l, workspace % tmp_sph_l)
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
        & workspace % tmp_grid, workspace % tmp_sph_l)
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

!> Adjoint HSP matrix vector product
subroutine bstarx(params, constants, workspace, x, y)
    ! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph)
    real(dp), intent(out) :: y(constants % nbasis, params % nsph)
    ! Local variables
    integer :: isph, jsph, ij, indmat, igrid, iproc
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
            call adjrhs_lpb(params, constants, workspace, isph, workspace % tmp_grid, &
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
subroutine bx(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    integer :: isph, jsph, ij, iproc

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
        !$omp parallel do default(none) shared(params,constants,workspace,x,y) &
        !$omp private(isph,iproc) schedule(dynamic)
        do isph = 1, params % nsph
          iproc = omp_get_thread_num() + 1
          call calcv2_lpb(params, constants, isph, workspace % tmp_pot(:, iproc), x, &
              & workspace % tmp_vylm, workspace % tmp_vplm, workspace % tmp_vcos, &
              & workspace % tmp_vsin, workspace % tmp_bessel)
          call ddintegrate(1, constants % nbasis, params % ngrid, constants % vwgrid, &
              & constants % vgrid_nbasis, workspace % tmp_pot(:, iproc), y(:,isph))
          y(:,isph) = - y(:,isph) 
        end do
    end if
    if (constants % dodiag) y = y + x
end subroutine bx

!> Adjoint ddLPB matrix-vector product
subroutine tstarx(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y

    real(dp), dimension(constants % nbasis, params % nsph, 2) :: temp_vector
    y = zero
    temp_vector = zero

    ! Compute AXr
    ! call LstarXr
    call lstarx(params, constants, workspace, x(:,:,1), temp_vector(:,:,1))
    ! Remove the scaling factor
    call convert_ddcosmo(params, constants, 1, temp_vector(:,:,1))
    y(:,:,1) = y(:,:,1) + temp_vector(:,:,1)
    ! Compute BXe
    call bstarx(params, constants, workspace, x(:,:,2), temp_vector(:,:,2))
    y(:,:,2) = y(:,:,2) + temp_vector(:,:,2)

    ! Reset temp_vector to zero
    temp_vector = zero
    ! Call CX
    call cstarx(params, constants, workspace, x, temp_vector)
    y = y + temp_vector
end subroutine tstarx

!> Apply the preconditioner to the primal HSP linear system
subroutine bx_prec(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph), intent(out) :: y
    y = x
end subroutine bx_prec

!> Apply the preconditioner to the ddLPB adjoint linear system
subroutine prec_tstarx(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph, 2)
    real(dp), intent(inout) :: y(constants % nbasis, params % nsph, 2)
    integer :: n_iter
    real(dp) :: r_norm
    real(dp), dimension(params % maxiter) :: x_rel_diff

    y(:,:,1) = x(:,:,1)
    call convert_ddcosmo(params, constants, 1, y(:,:,1))
    n_iter = params % maxiter
    call jacobi_diis(params, constants, workspace, constants % inner_tol, y(:,:,1), &
        & workspace % ddcosmo_guess, n_iter, x_rel_diff, lstarx, ldm1x, hnorm)
    if (workspace % error_flag.ne.0) then
        workspace % error_message = 'prec_tstarx: ddCOSMO failed to converge'
        return
    end if
    y(:,:,1) = workspace % ddcosmo_guess

    n_iter = params % maxiter
    call jacobi_diis(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
        & n_iter, x_rel_diff, bstarx, bx_prec, hnorm)
    if (workspace % error_flag.ne.0) then
        workspace % error_message = 'prec_tstarx: HSP failed to converge'
        return
    end if
    y(:,:,2) = workspace % hsp_guess

end subroutine prec_tstarx

!> Apply the preconditioner to the primal ddLPB linear system
!! |Yr| = |A^-1 0 |*|Xr|
!! |Ye|   |0 B^-1 | |Xe|
!! @param[in] ddx_data : dd Data
!! @param[in] x        : Input array
!! @param[out] y       : Linear system solution at current iteration
subroutine prec_tx(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), intent(in) :: x(constants % nbasis, params % nsph, 2)
    real(dp), intent(inout) :: y(constants % nbasis, params % nsph, 2)
    integer :: n_iter
    real(dp) :: r_norm
    real(dp), dimension(params % maxiter) :: x_rel_diff

    ! perform A^-1 * Yr
    n_iter = params % maxiter
    call jacobi_diis(params, constants, workspace, constants % inner_tol, x(:,:,1), &
        & workspace % ddcosmo_guess, n_iter, x_rel_diff, lx, ldm1x, hnorm)
    if (workspace % error_flag.ne.0) then
        workspace % error_message = 'prec_tx: ddCOSMO failed to converge'
        return
    end if

    ! Scale by the factor of (2l+1)/4Pi
    y(:,:,1) = workspace % ddcosmo_guess
    call convert_ddcosmo(params, constants, 1, y(:,:,1))

    ! perform B^-1 * Ye
    n_iter = params % maxiter
    call jacobi_diis(params, constants, workspace, constants % inner_tol, x(:,:,2), workspace % hsp_guess, &
        & n_iter, x_rel_diff, bx, bx_prec, hnorm)
    y(:,:,2) = workspace % hsp_guess

    if (workspace % error_flag.ne.0) then
        workspace % error_message = 'prec_tx: HSP failed to converge'
        return
    end if
end subroutine prec_tx

!> ddLPB adjoint matrix-vector product
subroutine cstarx(params, constants, workspace, x, y)
    implicit none
    ! input/output
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y
    ! local
    real(dp), dimension(0:constants % lmax0) :: SK_rijn, DK_rijn
    real(dp), dimension(constants % nbasis) :: basloc, vplm
    real(dp), dimension(params % lmax + 1) :: vcos, vsin
    complex(dp) :: bessel_work(max(2, params % lmax+1))
    complex(dp) :: work_complex(constants % lmax0+1)
    real(dp) :: work(constants % lmax0+1)
    integer :: isph, igrid, jsph, ind, l0, m0, ind0, indl, inode, l, m
    real(dp), dimension(3) :: vij, vtij
    real(dp) :: val, epsilon_ratio
    real(dp), allocatable :: scratch(:,:), scratch0(:,:)

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
subroutine cx(params, constants, workspace, x, y)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(in) :: x
    real(dp), dimension(constants % nbasis, params % nsph, 2), intent(out) :: y

    integer :: isph, jsph, igrid, ind, l, m, ind0
    real(dp), dimension(3) :: sijn ,vij, vtij
    real(dp) :: term, rho, ctheta, stheta, cphi, sphi, rijn, val
    real(dp), dimension(constants % nbasis) :: basloc, vplm
    real(dp), dimension(params % lmax + 1) :: vcos, vsin
    real(dp), dimension(constants % lmax0 + 1) :: SK_rijn, DK_rijn
    complex(dp) :: work_complex(constants % lmax0 + 1)
    real(dp) :: work(constants % lmax0 + 1)
    integer :: indl, inode

    real(dp), allocatable :: diff_re(:,:), diff0(:,:)
    allocate(diff_re(constants % nbasis, params % nsph), &
        & diff0(constants % nbasis0, params % nsph))

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
        !$omp diff0,y) private(isph,igrid,val,vij,rijn,sijn,SK_rijn, &
        !$omp DK_rijn,work,rho,ctheta,stheta,cphi,sphi,basloc,vplm, &
        !$omp vcos,vsin,term,ind0,ind,vtij,work_complex)
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
    deallocate(diff_re, diff0)

end subroutine cx

end module ddx_operators

