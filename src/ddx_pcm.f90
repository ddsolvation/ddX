!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_pcm.f90
!! PCM solver
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

module ddx_pcm
use ddx_core
use ddx_operators
use ddx_solvers
implicit none

contains

!> ddPCM solver
!!
!! Solves the problem within PCM model using a domain decomposition approach.
!!
!! @param[in] ddx_data: ddX object with all input information
!! @param[in] phi_cav: Potential at cavity points
!! @param[in] gradphi_cav: Gradient of a potential at cavity points
!! @param[in] psi: TODO
!! @param[out] esolv: Solvation energy
!! @param[out] force: Analytical forces
!! @param[out] info
subroutine ddpcm(ddx_data, phi_cav, gradphi_cav, psi, esolv, force, info)
    ! Inputs:
    type(ddx_type), intent(inout)  :: ddx_data
    real(dp), intent(in) :: phi_cav(ddx_data % ncav), &
        & gradphi_cav(3, ddx_data % ncav), psi(ddx_data % nbasis, ddx_data % nsph)
    ! Outputs
    real(dp), intent(out) :: esolv, force(3, ddx_data % nsph)
    integer, intent(out) :: info
    ! Local variables
    real(dp), allocatable :: vsin(:), vcos(:), vplm(:), basloc(:), &
        & dbsloc(:, :), fx(:, :), ef(:, :)
    integer :: istatus, isph, niter, igrid, icav, inear, inode, jnear, jnode, &
        & jsph
    real(dp) :: start_time, finish_time, tmp1, tmp2, d(3), dnorm
    logical :: ok
    double precision, external :: ddot, dnrm2
    external :: dgemm
    ! Zero initial guess on X (solution of the ddCOSMO system)
    ddx_data % xs = zero
    ! Get the energy
    call ddpcm_energy(ddx_data % params, ddx_data % constants, &
        & ddx_data % workspace, phi_cav, psi, ddx_data % xs, esolv, &
        & ddx_data % phi_grid, ddx_data % phi, ddx_data % phiinf, &
        & ddx_data % phieps, info)
    ! Get forces if needed
    if (ddx_data % force .eq. 1) then
        ! Zero initial guesses for adjoint systems
        ddx_data % s = zero
        ddx_data % y = zero
        ! Solve adjoint systems
        call ddpcm_adjoint(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, psi, ddx_data % s, ddx_data % y, info)
        ! Get forces, they are initialized with zeros in gradr
        call ddpcm_forces(ddx_data % params, ddx_data % constants, &
            & ddx_data % workspace, ddx_data % phi_grid, gradphi_cav, &
            & psi, ddx_data % phi, ddx_data % phieps, ddx_data % s, &
            & ddx_data % sgrid, ddx_data % y, ddx_data % ygrid, ddx_data % g, &
            & ddx_data % q, ddx_data % qgrid, ddx_data % xs, ddx_data % zeta, &
            & force, info)
    end if
end subroutine ddpcm

!> Solve primal ddPCM system to find solvation energy
!!
!! @param[in] params
!! @param[in] constants
!! @param[inout] workspace
!! @param[in] phi_cav
!! @param[in] psi
!! @param[inout] xs
!! @param[out] esolv
!! @param[in] info
subroutine ddpcm_energy(params, constants, workspace, phi_cav, psi, xs, &
        & esolv, phi_grid, phi, phiinf, phieps, info)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: phi_cav(constants % ncav), &
        & psi(constants % nbasis, params % nsph)
    real(dp), intent(inout) :: xs(constants % nbasis, params % nsph)
    !! Temporary buffers
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Outputs
    real(dp), intent(out) :: esolv, phi_grid(params % ngrid, params % nsph), &
        & phi(constants % nbasis, params % nsph), &
        & phiinf(constants % nbasis, params % nsph), &
        & phieps(constants % nbasis, params % nsph)
    integer, intent(out) :: info
    !! Local variables
    logical :: ok
    real(dp), external :: ddot
    external :: dgemm
    ! Unwrap sparsely stored potential at cavity points phi_cav into phi_grid
    ! and multiply it by characteristic function at cavity points ui
    call wghpot(constants % ncav, phi_cav, params % nsph, params % ngrid, &
        & constants % ui, phi_grid, workspace % tmp_grid)
    ! Integrate against spherical harmonics and Lebedev weights to get Phi
    call intrhs(params % nsph, constants % nbasis, params % ngrid, &
        & constants % vwgrid, constants % vgrid_nbasis, &
        & workspace % tmp_grid, phi)
    ! Compute Phi_infty
    call rinfx(params, constants, workspace, phi, phiinf)
    ! Set initial guess on Phi_epsilon as Phi
    phieps = phi
    ! Maximum number of iterations for an iterative solver
    info = params % maxiter
    ! Solve ddPCM system R_eps Phi_epsilon = Phi_infty
    call jacobi_diis(params, constants, workspace, constants % n, 0, &
        & params % ndiis, &
        & 4, params % tol, phiinf, phieps, info, ok, &
        & rx, apply_repsx_prec, hnorm)
    ! Solve ddCOSMO system L X = -Phi_epsilon with a zero initial guess
    info = params % maxiter
    call jacobi_diis(params, constants, workspace, constants % n, 0, &
        & params % ndiis, &
        & 4, params % tol, phieps, xs, info, ok, lx, &
        & ldm1x, hnorm)
    ! Solvation energy is computed
    esolv = pt5*ddot(constants % n, xs, 1, psi, 1)
end subroutine ddpcm_energy

subroutine ddpcm_adjoint(params, constants, workspace, psi, s, y, info)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: psi(constants % nbasis, params % nsph)
    !! Temporary buffers
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Outputs
    real(dp), intent(inout) :: s(constants % nbasis, params % nsph), &
        & y(constants % nbasis, params % nsph)
    integer, intent(out) :: info
    !! Local variables
    logical :: ok
    !! The code
    ! Solve the adjoint ddCOSMO system
    info = params % maxiter
    call jacobi_diis(params, constants, workspace, constants % n, 0, &
        & params % ndiis, 4, params % tol, psi, s, info, ok, &
        & lstarx, ldm1x, hnorm)
    ! Solve adjoint ddPCM system
    info = params % maxiter
    call jacobi_diis(params, constants, workspace, constants % n, 0, &
        & params % ndiis, 4, params % tol, s, y, info, ok, &
        & rstarx, apply_rstarepsx_prec, hnorm)
end subroutine ddpcm_adjoint

subroutine ddpcm_forces(params, constants, workspace, phi_grid, gradphi_cav, &
        & psi, phi, phieps, s, sgrid, y, ygrid, g, q, qgrid, xs, zeta, force, info)
    !! Inputs
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    real(dp), intent(in) :: phi_grid(params % ngrid, params % nsph), &
        & gradphi_cav(3, constants % ncav), &
        & psi(constants % nbasis, params % nsph), &
        & phi(constants % nbasis, params % nsph), &
        & phieps(constants % nbasis, params % nsph), &
        & s(constants % nbasis, params % nsph), &
        & y(constants % nbasis, params % nsph), &
        & xs(constants % nbasis, params % nsph)
    !! Temporary buffers
    type(ddx_workspace_type), intent(inout) :: workspace
    !! Outputs
    real(dp), intent(out) :: sgrid(params % ngrid, params % nsph), &
        & ygrid(params % ngrid, params % nsph), &
        & g(constants % nbasis, params % nsph), &
        & q(constants % nbasis, params % nsph), &
        & qgrid(params % ngrid, params % nsph), &
        & zeta(constants % ncav), &
        & force(3, params % nsph)
    integer, intent(out) :: info
    !! Local variables
    integer :: isph, icav, igrid, inode, jnode, jsph, jnear
    real(dp) :: tmp1, tmp2, d(3), dnorm
    real(dp), external :: ddot, dnrm2
    !! The code
    ! Get grid values of S and Y
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & constants % nbasis, one, constants % vgrid, constants % vgrid_nbasis, &
        & s, constants % nbasis, zero, sgrid, params % ngrid)
    call dgemm('T', 'N', params % ngrid, params % nsph, &
        & constants % nbasis, one, constants % vgrid, constants % vgrid_nbasis, &
        & y, constants % nbasis, zero, ygrid, params % ngrid)
    g = phi - phieps
    q = s - fourpi/(params % eps-one)*y
    qgrid = sgrid - fourpi/(params % eps-one)*ygrid
    ! gradr initializes forces with zeros
    call gradr(params, constants, workspace, g, ygrid, force)
    do isph = 1, params % nsph
        call fdoka(params, constants, isph, xs, sgrid(:, isph), &
            & workspace % tmp_vylm(:, 1), workspace % tmp_vdylm(:, :, 1), &
            & workspace % tmp_vplm(:, 1), workspace % tmp_vcos(:, 1), &
            & workspace % tmp_vsin(:, 1), force(:, isph)) 
        call fdokb(params, constants, isph, xs, sgrid, &
            & workspace % tmp_vylm(:, 1), workspace % tmp_vdylm(:, :, 1), &
            & workspace % tmp_vplm(:, 1), workspace % tmp_vcos(:, 1), &
            & workspace % tmp_vsin(:, 1), force(:, isph))
        call fdoga(params, constants, isph, qgrid, phi_grid, force(:, isph)) 
    end do
    force = -pt5 * force
    icav = 0
    do isph = 1, params % nsph
        do igrid = 1, params % ngrid
            if(constants % ui(igrid, isph) .ne. zero) then
                icav = icav + 1
                zeta(icav) = -pt5 * constants % wgrid(igrid) * &
                    & constants % ui(igrid, isph) * ddot(constants % nbasis, &
                    & constants % vgrid(1, igrid), 1, q(1, isph), 1)
                force(:, isph) = force(:, isph) + &
                    & zeta(icav)*gradphi_cav(:, icav)
            end if
        end do
    end do
    !! Last term where we compute gradients of potential at centers of atoms
    !! spawned by intermediate zeta.
    if(params % fmm .eq. 1) then
        !! This step can be substituted by a proper dgemm if zeta
        !! intermediate is converted from cavity points to all grid points
        !! with zero values at internal grid points
        ! P2M step
        icav = 0
        do isph = 1, params % nsph
            inode = constants % snode(isph)
            workspace % tmp_node_m(:, inode) = zero
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
                icav = icav + 1
                call fmm_p2m(constants % cgrid(:, igrid), zeta(icav), &
                    & one, params % pm, constants % vscales, one, &
                    & workspace % tmp_node_m(:, inode))
            end do
            workspace % tmp_node_m(:, inode) = workspace % tmp_node_m(:, inode) / &
                & params % rsph(isph)
        end do
        ! M2M, M2L and L2L translations
        call tree_m2m_rotation(params, constants, workspace % tmp_node_m)
        call tree_m2l_rotation(params, constants, workspace % tmp_node_m, &
            & workspace % tmp_node_l)
        call tree_l2l_rotation(params, constants, workspace % tmp_node_l)
        ! Now compute near-field FMM gradients
        ! Cycle over all spheres
        icav = 0
        workspace % tmp_efld = zero
        do isph = 1, params % nsph
            ! Cycle over all external grid points
            do igrid = 1, params % ngrid
                if(constants % ui(igrid, isph) .eq. zero) cycle
                icav = icav + 1
                ! Cycle over all near-field admissible pairs of spheres,
                ! including pair (isph, isph) which is a self-interaction
                inode = constants % snode(isph)
                do jnear = constants % snear(inode), constants % snear(inode+1)-1
                    ! Near-field interactions are possible only between leaf
                    ! nodes, which must contain only a single input sphere
                    jnode = constants % near(jnear)
                    jsph = constants % order(constants % cluster(1, jnode))
                    d = params % csph(:, isph) + &
                        & constants % cgrid(:, igrid)*params % rsph(isph) - &
                        & params % csph(:, jsph)
                    dnorm = dnrm2(3, d, 1)
                    workspace % tmp_efld(:, jsph) = workspace % tmp_efld(:, jsph) + &
                        & zeta(icav)*d/(dnorm**3)
                end do
            end do
        end do
        ! Take into account far-field FMM gradients only if pl > 0
        if (params % pl .gt. 0) then
            tmp1 = one / sqrt(3d0) / constants % vscales(1)
            do isph = 1, params % nsph
                inode = constants % snode(isph)
                tmp2 = tmp1 / params % rsph(isph)
                workspace % tmp_efld(3, isph) = workspace % tmp_efld(3, isph) + &
                    & tmp2*workspace % tmp_node_l(3, inode)
                workspace % tmp_efld(1, isph) = workspace % tmp_efld(1, isph) + &
                    & tmp2*workspace % tmp_node_l(4, inode)
                workspace % tmp_efld(2, isph) = workspace % tmp_efld(2, isph) + &
                    & tmp2*workspace % tmp_node_l(2, inode)
            end do
        end if
        do isph = 1, params % nsph
            force(:, isph) = force(:, isph) + &
                & workspace % tmp_efld(:, isph)*params % charge(isph)
        end do
    ! Naive quadratically scaling implementation
    else
        ! This routines actually computes -grad, not grad
        call efld(constants % ncav, zeta, constants % ccav, params % nsph, &
            & params % csph, workspace % tmp_efld)
        do isph = 1, params % nsph
            force(:, isph) = force(:, isph) - &
                & workspace % tmp_efld(:, isph)*params % charge(isph)
        end do
    end if
end subroutine ddpcm_forces

end module ddx_pcm
