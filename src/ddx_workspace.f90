!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_workspace.f90
!! Allocation of workspace for ddX
!!
!! TODO: deallocation function
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-06-22

!> Workspace for temporary buffers
module ddx_workspace
! Get ddx_params_type and all compile-time definitions
use ddx_parameters
! Get ddx_constants_type and all run-time constants
use ddx_constants

implicit none

type ddx_workspace_type
    !> Temporary workspace for scalar at the grid points. Dimension is
    !!      (ngrid, nproc)
    real(dp), allocatable :: tmp_pot(:,:)
    !> Temporary workspace for associated legendre polynomials. Dimension is
    !!      (vgrid_nbasis, nproc).
    real(dp), allocatable :: tmp_vplm(:, :)
    !> Temporary workspace for spherical harmomincs. Dimension is
    !!      (vgrid_nbasis, nproc).
    real(dp), allocatable :: tmp_vylm(:, :)
    !> Temporary workspace for derivatives of spherical harmomincs. Dimension is
    !!      (3, vgrid_nbasis, nproc).
    real(dp), allocatable :: tmp_vdylm(:, :, :)
    !> Temporary workspace of size
    !!      (vgrid_dmax + 1, nproc)
    real(dp), allocatable :: tmp_work(:, :)
    !> Temporary workspace for an array of cosinuses of a dimension
    !!      (vgrid_dmax+1, nproc).
    real(dp), allocatable :: tmp_vcos(:, :)
    !> Temporary workspace for an array of sinuses of a dimension
    !!      (vgrid_dmax+1, nproc).
    real(dp), allocatable :: tmp_vsin(:, :)
    !> Temporary workspace for Bessel functions. Dimension is
    !!      (max(2,lmax+1), nproc)
    complex(dp), allocatable :: tmp_bessel(:, :)
    !> Temporary workspace for multipole coefficients of a degree up to lmax
    !!      of each sphere. Dimension is (nbasis, nsph).
    real(dp), allocatable :: tmp_sph(:, :)
    !> Temporary workspace for multipole coefficients of a degree up to lmax+1
    !!      of each sphere. Dimension is (grad_nbasis, nsph). Allocated and
    !!      used only if fmm=1.
    real(dp), allocatable :: tmp_sph2(:, :)
    !> Temporary workspace for right hand side for solvers. Dimension is (nbasis, nsph).
    real(dp), allocatable :: tmp_rhs(:, :)
    !> Temporary workspace for a gradient of M2M of harmonics of a degree up to
    !!      lmax+1 of each sphere. Dimension is ((grad_nbasis, 3, nsph).
    real(dp), allocatable :: tmp_sph_grad(:, :, :)
    !> Temporary workspace for local coefficients of a degree up to pl
    !!      of each sphere. Dimension is ((pl+1)**2, nsph).
    real(dp), allocatable :: tmp_sph_l(:, :)
    !> Temporary workspace for a gradient of L2L of harmonics of a degree up to
    !!      pl of each sphere. Dimension is ((pl+1)**2, 3, nsph).
    real(dp), allocatable :: tmp_sph_l_grad(:, :, :)
    !> Temporary workspace for a gradient of L2L of harmonics of a degree up to
    !!      pl of each sphere. Dimension is ((pl+1)**2, 3, nsph).
    real(dp), allocatable :: tmp_sph_l_grad2(:, :, :)
    !> Temporary workspace for multipole coefficients of each node. Dimension
    !!      is ((pm+1)**2, nsph)
    real(dp), allocatable :: tmp_node_m(:, :)
    !> Temporary workspace for local coefficients of each node. Dimension is
    !!      ((pl+1)**2, nsph)
    real(dp), allocatable :: tmp_node_l(:, :)
    !> Temporary workspace for grid values of each sphere. Dimension is
    !!      (ngrid, nsph).
    real(dp), allocatable :: tmp_grid(:, :)
    !> Temporary workspace for grid values of each sphere. Dimension is
    !!      (ngrid, nsph).
    real(dp), allocatable :: tmp_grid2(:, :)
    !> Temporary workspace for values at cavity points. Dimension is
    !!      (ncav).
    real(dp), allocatable :: tmp_cav(:)
    !> Temporary electric field in cavity points. Dimension is (3, ncav).
    real(dp), allocatable :: tmp_efld(:, :)
    !> Jacobi solver temporary vector x. Dimension is (n).
    real(dp), allocatable :: tmp_x_new(:)
    !> Jacobi solver temporary result of matvec. Dimension is (n).
    real(dp), allocatable :: tmp_y(:)
    !> Jacobi solver temporary DIIS array. Dimension is (n, ndiis).
    real(dp), allocatable :: tmp_x_diis(:, :)
    !> Jacobi solver temporary DIIS array. Dimension is (n, ndiis).
    real(dp), allocatable :: tmp_e_diis(:, :)
    !> Jacobi solver temporary DIIS array. Dimension is (ndiis+1, ndiis+1).
    real(dp), allocatable :: tmp_bmat(:, :)
    !> GMRESR temporary workspace. Dimension is (n, 2*gmres_j+gmres_dim+2)
    real(dp), allocatable :: tmp_gmresr(:, :)
    !> ddLPB solutions for the microiterations
    real(dp), allocatable :: ddcosmo_guess(:,:), hsp_guess(:,:)
    !> Flag if there were an error
    integer :: error_flag = 2
    !> Last error message
    character(len=255) :: error_message
end type ddx_workspace_type

contains

subroutine workspace_init(params, constants, workspace, info)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    !! Outputs
    type(ddx_workspace_type), intent(out) :: workspace
    integer, intent(out) :: info
    !! Local variables
    character(len=255) :: string
    !! The code
    allocate(workspace % tmp_pot(params % ngrid, params % nproc), &
        & workspace % tmp_vplm(constants % vgrid_nbasis, params % nproc), &
        & workspace % tmp_vcos(constants % vgrid_dmax+1, params % nproc), &
        & workspace % tmp_vsin(constants % vgrid_dmax+1, params % nproc), &
        & workspace % tmp_work(constants % vgrid_dmax+1, params % nproc), &
        & stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        string = "workspace_init: `tmp_vplm`, `tmp_vcos` and `tmp_vsin` " &
            & // "allocations failed"
        workspace % error_message = string
        call params % print_func(string)
        info = 1
        return
    end if
    allocate(workspace % tmp_vylm(constants % vgrid_nbasis, params % nproc), &
        & workspace % tmp_vdylm(3, constants % vgrid_nbasis, params % nproc), &
        & stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_vylm` " &
            & // "and `tmp_vdylm` allocations failed"
        info = 1
        return
    end if
    allocate(workspace % tmp_sph(constants % nbasis, params % nsph), &
        & stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_sph` " // &
            & "allocation failed"
        info = 1
        return
    end if
    if (params % fmm .eq. 1) then
        allocate(workspace % tmp_sph2(constants % grad_nbasis, params % nsph), &
            & stat=info)
        if (info .ne. 0) then
            workspace % error_flag = 1
            workspace % error_message = "workspace_init: `tmp_sph2` " // &
                & "allocation failed"
            info = 1
            return
        end if
        allocate(workspace % tmp_sph_grad( &
            & constants % grad_nbasis, 3, params % nsph), &
            & stat=info)
        if (info .ne. 0) then
            workspace % error_flag = 1
            workspace % error_message = "workspace_init: `tmp_sph_grad` " // &
                & "allocation failed"
            info = 1
            return
        end if
        allocate(workspace % tmp_sph_l((params % pl+1)**2, params % nsph), &
            & stat=info)
        if (info .ne. 0) then
            workspace % error_flag = 1
            workspace % error_message = "workspace_init: `tmp_sph_l` " // &
                & "allocation failed"
            info = 1
            return
        end if
        allocate(workspace % tmp_sph_l_grad( &
            & (params % pl+1)**2, 3, params % nsph), &
            & stat=info)
        if (info .ne. 0) then
            workspace % error_flag = 1
            workspace % error_message = "workspace_init: `tmp_sph_l_grad` " // &
                & "allocation failed"
            info = 1
            return
        end if
        allocate(workspace % tmp_sph_l_grad2( &
            & (params % pl+1)**2, 3, params % nsph), &
            & stat=info)
        if (info .ne. 0) then
            workspace % error_flag = 1
            workspace % error_message = "workspace_init: `tmp_sph_l_grad2` " // &
                & "allocation failed"
            info = 1
            return
        end if
        allocate(workspace % tmp_node_m((params % pm+1)**2, &
            & constants % nclusters), stat=info)
        if (info .ne. 0) then
            workspace % error_flag = 1
            workspace % error_message = "workspace_init: `tmp_node_m` " // &
                & "allocation failed"
            info = 1
            return
        end if
        allocate(workspace % tmp_node_l((params % pl+1)**2, &
            & constants % nclusters), stat=info)
        if (info .ne. 0) then
            workspace % error_flag = 1
            workspace % error_message = "workspace_init: `tmp_node_l` " // &
                & "allocation failed"
            info = 1
            return
        end if
    end if
    allocate(workspace % tmp_grid(params % ngrid, params % nsph), &
        & stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_grid` " // &
            & "allocation failed"
        info = 1
        return
    end if
    allocate(workspace % tmp_grid2(params % ngrid, params % nsph), &
        & stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_grid2` " // &
            & "allocation failed"
        info = 1
        return
    end if
    allocate(workspace % tmp_cav(constants % ncav), stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_cav` " // &
            & "allocation failed"
        info = 1
        return
    end if
    allocate(workspace % tmp_efld(3, constants % ncav), stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_efld` " // &
            & "allocation failed"
        info = 1
        return
    end if
    allocate(workspace % tmp_x_new(constants % n), stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_x_new` " // &
            & "allocation failed"
        info = 1
        return
    end if
    allocate(workspace % tmp_y(constants % n), stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_y` " // &
            & "allocation failed"
        info = 1
        return
    end if
    allocate(workspace % tmp_x_diis(constants % n, 2*params % jacobi_ndiis), &
        & stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_x_diis` " // &
            & "allocation failed"
        info = 1
        return
    end if
    allocate(workspace % tmp_e_diis(constants % n, 2*params % jacobi_ndiis), &
        & stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_e_diis` " // &
            & "allocation failed"
        info = 1
        return
    end if
    allocate(workspace % tmp_bmat(2*params % jacobi_ndiis + 2, &
        & 2*params % jacobi_ndiis + 2), stat=info)
    if (info .ne. 0) then
        workspace % error_flag = 1
        workspace % error_message = "workspace_init: `tmp_bmat` " // &
            & "allocation failed"
        info = 1
        return
    end if
    ! In case of GMRESR iterative solver allocate corresponding temporary space
    if ((params % itersolver .eq. 2) .or. (params % model .eq. 3)) then
        allocate(workspace % tmp_gmresr(constants % n, &
            & 0:2*params % gmresr_j + params % gmresr_dim + 1), stat=info)
        if (info .ne. 0) then
            workspace % error_flag = 1
            workspace % error_message = "workspace_init: `tmp_gmresr` " // &
                & "allocation failed"
            info = 1
            return
        end if
    end if
    ! Allocations for LPB model
    if (params % model .eq. 3) then
        allocate(workspace % tmp_bessel(max(2, params % lmax+1), &
            & params % nproc), workspace % ddcosmo_guess(constants % nbasis, params % nsph), &
            & workspace % hsp_guess(constants % nbasis, params % nsph), stat=info)
        if (info .ne. 0) then
            workspace % error_flag = 1
            workspace % error_message = "workspace_init: `tmp_bessel` " // &
                & "allocation failed"
            info = 1
            return
        end if
    end if
    ! Clear error state
    info = 0
    workspace % error_flag = 0
    workspace % error_message = ""
end subroutine workspace_init

subroutine workspace_free(workspace, info)
    implicit none
    type(ddx_workspace_type), intent(out) :: workspace
    integer, intent(out) :: info
    integer :: istat

    istat = 0
    info = 0

    if (allocated(workspace % tmp_pot)) then
        deallocate(workspace % tmp_pot, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_pot` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_vplm)) then
        deallocate(workspace % tmp_vplm, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_vplm` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_vcos)) then
        deallocate(workspace % tmp_vcos, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_vcos` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_vsin)) then
        deallocate(workspace % tmp_vsin, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_vsin` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_work)) then
        deallocate(workspace % tmp_work, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_work` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_vylm)) then
        deallocate(workspace % tmp_vylm, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_vylm` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_vdylm)) then
        deallocate(workspace % tmp_vdylm, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_vdylm` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_sph)) then
        deallocate(workspace % tmp_sph, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_sph` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_sph2)) then
        deallocate(workspace % tmp_sph2, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_sph2` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_sph_grad)) then
        deallocate(workspace % tmp_sph_grad, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_sph_grad` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_sph_l)) then
        deallocate(workspace % tmp_sph_l, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_sph_l` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_sph_l_grad)) then
        deallocate(workspace % tmp_sph_l_grad, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_sph_l_grad` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_sph_l_grad2)) then
        deallocate(workspace % tmp_sph_l_grad2, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_sph_l_grad2` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_node_m)) then
        deallocate(workspace % tmp_node_m, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_node_m` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_node_l)) then
        deallocate(workspace % tmp_node_l, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_node_l` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_grid)) then
        deallocate(workspace % tmp_grid, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_grid` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_grid2)) then
        deallocate(workspace % tmp_grid2, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_grid2` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_cav)) then
        deallocate(workspace % tmp_cav, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_cav` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_efld)) then
        deallocate(workspace % tmp_efld, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_efld` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_x_new)) then
        deallocate(workspace % tmp_x_new, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_x_new` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_y)) then
        deallocate(workspace % tmp_y, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_y` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_x_diis)) then
        deallocate(workspace % tmp_x_diis, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_x_diis` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_e_diis)) then
        deallocate(workspace % tmp_e_diis, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_e_diis` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_bmat)) then
        deallocate(workspace % tmp_bmat, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_bmat` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_gmresr)) then
        deallocate(workspace % tmp_gmresr, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_gmresr` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_bessel)) then
        deallocate(workspace % tmp_bessel, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_bessel` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % ddcosmo_guess)) then
        deallocate(workspace % ddcosmo_guess, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`ddcosmo_guess` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % hsp_guess)) then
        deallocate(workspace % hsp_guess, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`hsp_guess` deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(workspace % tmp_rhs)) then
        deallocate(workspace % tmp_rhs, stat=istat)
        if (istat .ne. 0) then
            info = 1
            write(6, *) "`tmp_rhs` deallocation failed!"
            stop 1
        end if
    end if
end subroutine workspace_free

end module

