!> Routines which are only used by the tests and should not be used by
!! external software.
module ddx_legacy
use ddx
implicit none

contains

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
        & gradphi_cav, hessian_flag, hessianphi_cav, psi, charges)
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
    real(dp), intent(in) :: charges(params % nsph)
    ! Local variables
    integer :: isph, igrid, icav, inode, jnear, jnode, jsph, i
    real(dp) :: d(3), v, tmpv, r, gradv(3), hessianv(3, 3), tmpd(3)
    real(dp), allocatable :: grid_grad(:,:,:), grid_hessian(:,:,:,:), &
        & grid_hessian2(:,:,:)
    real(dp), external :: dnrm2

    write(6, *) "Warning: subroutine mkrhs is deprecated"

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
                tmpv = charges(isph) / r
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
                & charges(isph) &
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
            !$omp grid_grad,grid_hessian,grad_flag,hessian_flag,charges) &
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
                        tmpv = charges(jsph) / r
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
        psi(1, isph) = sqrt4pi * charges(isph)
    end do

    ! deallocate temporary
    if (grad_flag .eq. 1) deallocate(grid_grad)
    if (hessian_flag .eq. 1) deallocate(grid_hessian, grid_hessian2)
end subroutine mkrhs

!> Main solver routine
!!
!! Solves the problem within COSMO model using a domain decomposition approach.
!!
!! @param[in] ddx_data: ddX object with all input information
!! @param[in] phi_cav: Potential at cavity points
!! @param[in] e_cav: Gradient of the eletric potential at cavity points
!! @param[in] hessianphi_cav: Hessian of the eletric potential at cavity points
!! @param[inout] state: ddx state (contains RHSs and solutions)
!! @param[in] psi: RHS of the adjoint problem
!! @param[in] tol
!! @param[out] esolv: Solvation energy
!! @param[out] force: Analytical forces
!! @param[inout] ddx_error: ddX error
!!
subroutine ddsolve_legacy(ddx_data, state, phi_cav, e_cav, hessianphi_cav, &
        & psi, tol, esolv, force, ddx_error)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(ddx_data % constants % ncav), &
        & e_cav(3, ddx_data % constants % ncav), &
        & hessianphi_cav(3, ddx_data % constants % ncav), &
        & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), tol
    ! Outputs
    real(dp), intent(out) :: esolv, force(3, ddx_data % params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error

    write(6, *) "Warning: subroutine ddsolve_legacy is deprecated"

    ! Find proper model
    select case(ddx_data % params % model)
        ! COSMO model
        case (1)
            call ddcosmo(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, phi_cav, psi, e_cav, &
                & tol, esolv, force, ddx_error)
        ! PCM model
        case (2)
            call ddpcm(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, phi_cav, psi, e_cav, &
                & tol, esolv, force, ddx_error)
        ! LPB model
        case (3)
            call ddlpb(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, phi_cav, e_cav, &
                & psi, tol, esolv, hessianphi_cav, force, ddx_error)
        ! Error case
        case default
            call update_error(ddx_error, "unsupported solvation " // &
                & " model in the dd solver.")
            return
    end select
end subroutine ddsolve_legacy

!> ddCOSMO solver
!!
!! Solves the problem within COSMO model using a domain decomposition approach.
!!
!> @ingroup Fortran_interface_ddcosmo
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] phi_cav: Potential at cavity points, size (ncav)
!! @param[in] psi: Representation of the solute potential in spherical
!!     harmonics, size (nbasis, nsph)
!! @param[in] tol: Tolerance for the linear system solver
!! @param[out] esolv: Solvation energy
!! @param[out] force: Solvation contribution to the forces
!! @param[inout] ddx_error: ddX error
!!
subroutine ddcosmo(params, constants, workspace, state, phi_cav, &
        & psi, e_cav, tol, esolv, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(constants % ncav), &
        & psi(constants % nbasis, params % nsph), tol
    real(dp), intent(out) :: esolv
    real(dp), intent(in) :: e_cav(3, constants % ncav)
    real(dp), intent(out), optional :: force(3, params % nsph)
    type(ddx_error_type), intent(inout) :: ddx_error

    write(6, *) "Warning: subroutine ddcosmo is deprecated"

    call cosmo_setup(params, constants, workspace, state, phi_cav, psi, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "ddlpb: ddcosmo_setup returned an error, exiting")
        return
    end if
    call cosmo_guess(params, constants, workspace, state, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "ddlpb: ddcosmo_guess returned an error, exiting")
        return
    end if
    call cosmo_solve(params, constants, workspace, state, tol, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "ddlpb: ddcosmo_solve returned an error, exiting")
        return
    end if

    call cosmo_energy(constants, state, esolv, ddx_error)

    ! Get forces if needed
    if (params % force .eq. 1) then
        ! solve the adjoint
        call cosmo_guess_adjoint(params, constants, workspace, state, ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddlpb: ddcosmo_guess_adjoint returned an error, exiting")
            return
        end if
        call cosmo_solve_adjoint(params, constants, workspace, state, tol, &
            & ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddlpb: ddcosmo_guess_adjoint returned an error, exiting")
            return
        end if

        ! evaluate the solvent unspecific contribution analytical derivatives
        force = zero
        call cosmo_solvation_force_terms(params, constants, workspace, &
            & state, e_cav, force, ddx_error)
    end if
end subroutine ddcosmo

!> ddPCM solver
!!
!! Solves the problem within PCM model using a domain decomposition approach.
!!
!> @ingroup Fortran_interface_ddpcm
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] phi_cav: Potential at cavity points, size (ncav)
!! @param[in] psi: Representation of the solute potential in spherical
!!     harmonics, size (nbasis, nsph)
!! @param[in] tol: Tolerance for the linear system solver
!! @param[out] esolv: Solvation energy
!! @param[out] force: Solvation contribution to the forces
!! @param[inout] ddx_error: ddX error
!!
subroutine ddpcm(params, constants, workspace, state, phi_cav, &
        & psi, e_cav, tol, esolv, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(in) :: phi_cav(constants % ncav), &
        & psi(constants % nbasis, params % nsph), tol
    real(dp), intent(out) :: esolv
    real(dp), intent(in) :: e_cav(3, constants % ncav)
    real(dp), intent(out), optional :: force(3, params % nsph)

    write(6, *) "Warning: subroutine ddpcm is deprecated"

    call pcm_setup(params, constants, workspace, state, phi_cav, psi, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "ddlpb: ddpcm_setup returned an error, exiting")
        return
    end if
    call pcm_guess(params, constants, workspace, state, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "ddlpb: ddpcm_guess returned an error, exiting")
        return
    end if
    call pcm_solve(params, constants, workspace, state, tol, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "ddlpb: ddpcm_solve returned an error, exiting")
        return
    end if

    call pcm_energy(constants, state, esolv, ddx_error)

    ! Get forces if needed
    if (params % force .eq. 1) then
        ! solve the adjoint
        call pcm_guess_adjoint(params, constants, workspace, state, ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddlpb: ddpcm_guess_adjoint returned an error, exiting")
            return
        end if
        call pcm_solve_adjoint(params, constants, workspace, state, &
            & tol, ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddlpb: ddpcm_guess_adjoint returned an error, exiting")
            return
        end if

        ! evaluate the solvent unspecific contribution analytical derivatives
        force = zero
        call pcm_solvation_force_terms(params, constants, workspace, &
            & state, e_cav, force, ddx_error)
    end if

end subroutine ddpcm

!> ddLPB solver
!!
!! Solves the LPB problem using a domain decomposition approach.
!!
!! @param[in] params: User specified parameters
!! @param[in] constants: Precomputed constants
!! @param[inout] workspace: Preallocated workspaces
!! @param[inout] state: ddx state (contains solutions and RHSs)
!! @param[in] phi_cav: Potential at cavity points, size (ncav)
!! @params[in] e_cav: Electric field at cavity points, size (3, ncav)
!! @param[in] psi: Representation of the solute potential in spherical
!!     harmonics, size (nbasis, nsph)
!! @param[in] tol: Tolerance for the linear system solver
!! @param[out] esolv: Solvation energy
!! @param[out] force: Solvation contribution to the forces
!! @param[inout] ddx_error: ddX error
!!
subroutine ddlpb(params, constants, workspace, state, phi_cav, e_cav, &
        & psi, tol, esolv, hessianphi_cav, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(inout) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(constants % ncav), &
        & e_cav(3, constants % ncav), &
        & psi( constants % nbasis,  params % nsph), tol
    real(dp), intent(out) :: esolv
    real(dp), intent(out), optional :: force(3, params % nsph)
    real(dp), intent(in), optional :: hessianphi_cav(3, 3, constants % ncav)
    type(ddx_error_type), intent(inout) :: ddx_error

    write(6, *) "Warning: subroutine ddlpb is deprecated"

    call lpb_setup(params, constants, workspace, state, phi_cav, &
        & e_cav, psi, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "ddlpb: ddlpb_setup returned an error, exiting")
        return
    end if
    call lpb_guess(params, constants, workspace, state, tol, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "ddlpb: ddlpb_guess returned an error, exiting")
        return
    end if
    call lpb_solve(params, constants, workspace, state, tol, ddx_error)
    if (ddx_error % flag .ne. 0) then
        call update_error(ddx_error, &
            & "ddlpb: ddlpb_solve returned an error, exiting")
        return
    end if

    ! Compute the solvation energy
    call lpb_energy(constants, state, esolv, ddx_error)

    ! Get forces if needed
    if(params % force .eq. 1) then
        call lpb_guess_adjoint(params, constants, workspace, state, tol, ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddlpb: ddlpb_guess_adjoint returned an error, exiting")
            return
        end if
        call lpb_solve_adjoint(params, constants, workspace, state, tol, ddx_error)
        if (ddx_error % flag .ne. 0) then
            call update_error(ddx_error, &
                & "ddlpb: ddlpb_solve_adjoint returned an error, exiting")
            return
        end if
        call lpb_solvation_force_terms(params, constants, workspace, &
            & state, hessianphi_cav, force, ddx_error)
    endif

end subroutine ddlpb


end module ddx_legacy
