!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx.f90
!! Main driver routine of the ddX with all the per-model solvers
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

module ddx
use ddx_core
use ddx_operators
use ddx_solvers
implicit none

contains

!> Main solver routine
!!
!! Solves the problem within COSMO model using a domain decomposition approach.
!!
!! @param[in] ddx_data: ddX object with all input information
!! @param[in] phi_cav: Potential at cavity points
!! @param[in] gradphi_cav: Gradient of a potential at cavity points
!! @param[in] psi: TODO
!! @param[out] esolv: Solvation energy
!! @param[out] force: Analytical forces
subroutine ddsolve(ddx_data, phi_cav, gradphi_cav, psi, esolv, force)
    ! Inputs:
    type(ddx_type), intent(inout)  :: ddx_data
    real(dp), intent(in) :: phi_cav(ddx_data % ncav), &
        & gradphi_cav(3, ddx_data % ncav), psi(ddx_data % nbasis, ddx_data % nsph)
    ! Outputs
    real(dp), intent(out) :: esolv, force(3, ddx_data % nsph)
    ! Find proper model
    select case(ddx_data % model)
        ! COSMO model
        case (1)
            call ddcosmo(ddx_data, phi_cav, gradphi_cav, psi, esolv, force)
        ! PCM model
        case (2)
            call ddpcm(ddx_data, phi_cav, gradphi_cav, psi, esolv, force)
        ! LPB model
        case (3)
            stop "LPB model is not yet fully supported"
        ! Error case
        case default
            stop "Non-supported model"
    end select
end subroutine ddsolve

!> ddCOSMO solver
!!
!! Solves the problem within COSMO model using a domain decomposition approach.
!!
!! @param[in] ddx_data: ddX object with all input information
!! @param[in] phi_cav: Potential at cavity points
!! @param[in] gradphi_cav: Gradient of a potential at cavity points
!! @param[in] psi: TODO
!! @param[out] esolv: Solvation energy
!! @param[out] force: Analytical forces
subroutine ddcosmo(ddx_data, phi_cav, gradphi_cav, psi, esolv, force)
    ! Inputs:
    type(ddx_type), intent(inout)  :: ddx_data
    real(dp), intent(in) :: phi_cav(ddx_data % ncav), &
        & gradphi_cav(3, ddx_data % ncav), psi(ddx_data % nbasis, ddx_data % nsph)
    ! Outputs
    real(dp), intent(out) :: esolv, force(3, ddx_data % nsph)
    ! Local variables
    real(dp), allocatable :: vsin(:), vcos(:), vplm(:), basloc(:), &
        & dbsloc(:, :), fx(:, :), ef(:, :)
    integer :: istatus, isph, niter, igrid, icav, inear, inode, jnear, jnode, &
        & jsph
    real(dp) :: start_time, finish_time, tmp1, tmp2, d(3), dnorm
    logical :: ok
    real(dp), external :: ddot, dnrm2
    external :: dgemm
    ! Unwrap sparsely stored potential at cavity points and multiply by ui
    call wghpot(ddx_data, phi_cav, ddx_data % phi_grid, ddx_data % tmp_grid)
    ! Integrate against spherical harmonics and Lebedev weights to get Phi
    call dgemm('N', 'N', ddx_data % nbasis, ddx_data % nsph, ddx_data % ngrid, &
        & one, ddx_data % vwgrid, ddx_data % vgrid_nbasis, ddx_data % tmp_grid, &
        & ddx_data % ngrid, zero, ddx_data % phi, ddx_data % nbasis)
    ! Solve ddCOSMO system L X = -Phi with a zero initial guess
    niter = ddx_data % maxiter
    ddx_data % xs = zero
    call cpu_time(start_time)
    call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, ddx_data % ndiis, &
        & 4, ddx_data % tol, ddx_data % phi, ddx_data % xs, niter, ok, lx, &
        & ldm1x, hnorm)
    call cpu_time(finish_time)
    if (ddx_data % iprint .ge. 1) then
        write(*, "(A,ES11.4E2,A)") " ddcosmo step time:", &
            & finish_time-start_time, " seconds"
        write(*, "(A,I0)") " ddcosmo step iterations: ", niter
    endif
    ! Solvation energy is computed
    esolv = pt5*ddot(ddx_data % n, ddx_data % xs, 1, psi, 1)
    ! Get forces if needed
    if (ddx_data % force .eq. 1) then
        ! Solve adjoint ddCOSMO system
        niter = ddx_data % maxiter
        ddx_data % s = zero
        call cpu_time(start_time)
        call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, &
            & ddx_data % ndiis, 4, ddx_data % tol, psi, ddx_data % s, niter, ok, &
            & lstarx, ldm1x, hnorm)
        call cpu_time(finish_time)
        if (ddx_data % iprint.ge.1) then
            write(*, "(A,ES11.4E2,A)") " adjoint ddcosmo step time:", &
                & finish_time-start_time, " seconds"
            write(*, "(A,I0)") " adjoint ddcosmo step iterations: ", niter
        endif
        call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % vgrid, ddx_data % vgrid_nbasis, &
            & ddx_data % s, ddx_data % nbasis, zero, ddx_data % sgrid, &
            & ddx_data % ngrid)
        ddx_data % q = ddx_data % s
        ddx_data % qgrid = ddx_data % sgrid
        allocate(vsin(ddx_data % lmax+1), vcos(ddx_data % lmax+1), &
            & vplm(ddx_data % nbasis), basloc(ddx_data % nbasis), &
            & dbsloc(3, ddx_data % nbasis), fx(3, ddx_data % nsph), stat=istatus)
        if (istatus.ne.0) write(6,*) 'ddpcm forces allocation failed'
        force = zero
        do isph = 1, ddx_data % nsph
            call fdoka(ddx_data, isph, ddx_data % xs, ddx_data % sgrid(:, isph), &
                & basloc, dbsloc, vplm, vcos, vsin, force(:,isph))
            call fdokb(ddx_data, isph, ddx_data % xs, ddx_data % sgrid, basloc, &
                & dbsloc, vplm, vcos, vsin, force(:, isph))
            call fdoga(ddx_data, isph, ddx_data % sgrid, ddx_data % phi_grid, &
                & force(:, isph))
        end do
        force = -pt5 * force
        icav = 0
        do isph = 1, ddx_data % nsph
            do igrid = 1, ddx_data % ngrid
                if(ddx_data % ui(igrid, isph) .ne. zero) then
                    icav = icav + 1
                    ddx_data % zeta(icav) = -pt5 * ddx_data % wgrid(igrid) * &
                        & ddx_data % ui(igrid, isph) * ddot(ddx_data % nbasis, &
                        & ddx_data % vgrid(1, igrid), 1, &
                        & ddx_data % s(1, isph), 1)
                    force(:, isph) = force(:, isph) + &
                        & ddx_data % zeta(icav)*gradphi_cav(:, icav)
                end if
            end do
        end do
        !! Last term where we compute gradients of potential at centers of atoms
        !! spawned by intermediate zeta.
        allocate(ef(3, ddx_data % nsph))
        if(ddx_data % fmm .eq. 1) then
            !! This step can be substituted by a proper dgemm if zeta
            !! intermediate is converted from cavity points to all grid points
            !! with zero values at internal grid points
            ! P2M step
            icav = 0
            do isph = 1, ddx_data % nsph
                inode = ddx_data % snode(isph)
                ddx_data % tmp_node_m(:, inode) = zero
                do igrid = 1, ddx_data % ngrid
                    if(ddx_data % ui(igrid, isph) .eq. zero) cycle
                    icav = icav + 1
                    call fmm_p2m(ddx_data % cgrid(:, igrid), ddx_data % zeta(icav), &
                        & one, ddx_data % pm, ddx_data % vscales, one, &
                        & ddx_data % tmp_node_m(:, inode))
                end do
                ddx_data % tmp_node_m(:, inode) = ddx_data % tmp_node_m(:, inode) / &
                    & ddx_data % rsph(isph)
            end do
            ! M2M, M2L and L2L translations
            if(ddx_data % fmm_precompute .eq. 1) then
                call tree_m2m_reflection_use_mat(ddx_data, ddx_data % tmp_node_m)
                call tree_m2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_m, &
                    & ddx_data % tmp_node_l)
                call tree_l2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_l)
            else
                call tree_m2m_rotation(ddx_data, ddx_data % tmp_node_m)
                call tree_m2l_rotation(ddx_data, ddx_data % tmp_node_m, &
                    & ddx_data % tmp_node_l)
                call tree_l2l_rotation(ddx_data, ddx_data % tmp_node_l)
            end if
            ! Now compute near-field FMM gradients
            ! Cycle over all spheres
            icav = 0
            ef = zero
            do isph = 1, ddx_data % nsph
                ! Cycle over all external grid points
                do igrid = 1, ddx_data % ngrid
                    if(ddx_data % ui(igrid, isph) .eq. zero) cycle
                    icav = icav + 1
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
                        ef(:, jsph) = ef(:, jsph) + &
                            & ddx_data % zeta(icav)*d/(dnorm**3)
                    end do
                end do
            end do
            ! Take into account far-field FMM gradients only if pl > 0
            if (ddx_data % pl .gt. 0) then
                tmp1 = one / sqrt(3d0) / ddx_data % vscales(1)
                do isph = 1, ddx_data % nsph
                    inode = ddx_data % snode(isph)
                    tmp2 = tmp1 / ddx_data % rsph(isph)
                    ef(3, isph) = ef(3, isph) + tmp2*ddx_data % tmp_node_l(3, inode)
                    ef(1, isph) = ef(1, isph) + tmp2*ddx_data % tmp_node_l(4, inode)
                    ef(2, isph) = ef(2, isph) + tmp2*ddx_data % tmp_node_l(2, inode)
                end do
            end if
            do isph = 1, ddx_data % nsph
                force(:, isph) = force(:, isph) + ef(:, isph)*ddx_data % charge(isph)
            end do
        ! Naive quadratically scaling implementation
        else
            ! This routines actually computes -grad, not grad
            call efld(ddx_data % ncav, ddx_data % zeta, ddx_data % ccav, &
                & ddx_data % nsph, ddx_data % csph, ef)
            do isph = 1, ddx_data % nsph
                force(:, isph) = force(:, isph) - ef(:, isph)*ddx_data % charge(isph)
            end do
        end if
        deallocate(ef)
    end if
end subroutine ddcosmo

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
subroutine ddpcm(ddx_data, phi_cav, gradphi_cav, psi, esolv, force)
    ! Inputs:
    type(ddx_type), intent(inout)  :: ddx_data
    real(dp), intent(in) :: phi_cav(ddx_data % ncav), &
        & gradphi_cav(3, ddx_data % ncav), psi(ddx_data % nbasis, ddx_data % nsph)
    ! Outputs
    real(dp), intent(out) :: esolv, force(3, ddx_data % nsph)
    ! Local variables
    real(dp), allocatable :: vsin(:), vcos(:), vplm(:), basloc(:), &
        & dbsloc(:, :), fx(:, :), ef(:, :)
    integer :: istatus, isph, niter, igrid, icav, inear, inode, jnear, jnode, &
        & jsph
    real(dp) :: start_time, finish_time, tmp1, tmp2, d(3), dnorm
    logical :: ok
    double precision, external :: ddot, dnrm2
    external :: dgemm
    ! Unwrap sparsely stored potential at cavity points and multiply by ui
    call wghpot(ddx_data, phi_cav, ddx_data % phi_grid, ddx_data % tmp_grid)
    ! Integrate against spherical harmonics and Lebedev weights to get Phi
    call dgemm('N', 'N', ddx_data % nbasis, ddx_data % nsph, ddx_data % ngrid, &
        & one, ddx_data % vwgrid, ddx_data % vgrid_nbasis, ddx_data % tmp_grid, &
        & ddx_data % ngrid, zero, ddx_data % phi, ddx_data % nbasis)
    ! Compute Phi_infty
    call rinfx(ddx_data, ddx_data % phi, ddx_data % phiinf)
    ! Set initial guess on Phi_epsilon as Phi
    ddx_data % phieps = ddx_data % phi
    ! Maximum number of iterations for an iterative solver
    niter = ddx_data % maxiter
    ! Solve ddPCM system R_eps Phi_epsilon = Phi_infty
    call cpu_time(start_time)
    call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, ddx_data % ndiis, &
        & 4, ddx_data % tol, ddx_data % phiinf, ddx_data % phieps, niter, ok, &
        & rx, apply_repsx_prec, hnorm)
    call cpu_time(finish_time)
    if (ddx_data % iprint .ge. 1) then
        write(*, "(A,ES11.4E2,A)") " ddpcm step time:", &
            & finish_time-start_time, " seconds"
        write(*, "(A,I0)") " ddpcm step iterations: ", niter
    endif
    ! Solve ddCOSMO system L X = -Phi_epsilon with a zero initial guess
    niter = ddx_data % maxiter
    ddx_data % xs = zero
    call cpu_time(start_time)
    call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, ddx_data % ndiis, &
        & 4, ddx_data % tol, ddx_data % phieps, ddx_data % xs, niter, ok, lx, &
        & ldm1x, hnorm)
    call cpu_time(finish_time)
    if (ddx_data % iprint .ge. 1) then
        write(*, "(A,ES11.4E2,A)") " ddcosmo step time:", &
            & finish_time-start_time, " seconds"
        write(*, "(A,I0)") " ddcosmo step iterations: ", niter
    endif
    ! Solvation energy is computed
    esolv = pt5*ddot(ddx_data % n, ddx_data % xs, 1, psi, 1)
    ! Get forces if needed
    if (ddx_data % force .eq. 1) then
        ! Solve adjoint ddCOSMO system
        niter = ddx_data % maxiter
        ddx_data % s = zero
        call cpu_time(start_time)
        call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, &
            & ddx_data % ndiis, 4, ddx_data % tol, psi, ddx_data % s, niter, ok, &
            & lstarx, ldm1x, hnorm)
        call cpu_time(finish_time)
        if (ddx_data % iprint.ge.1) then
            write(*, "(A,ES11.4E2,A)") " adjoint ddcosmo step time:", &
                & finish_time-start_time, " seconds"
            write(*, "(A,I0)") " adjoint ddcosmo step iterations: ", niter
        endif
        ! Solve adjoint ddPCM system
        niter = ddx_data % maxiter
        ddx_data % y = zero
        call cpu_time(start_time)
        call jacobi_diis(ddx_data, ddx_data % n, ddx_data % iprint, &
            & ddx_data % ndiis, 4, ddx_data % tol, ddx_data % s, ddx_data % y, &
            & niter, ok, rstarx, apply_rstarepsx_prec, hnorm)
        call cpu_time(finish_time)
        if (ddx_data % iprint .ge. 1) then
            write(*,"(A,ES11.4E2,A)") " adjoint ddpcm step time:", &
                & finish_time-start_time, " seconds"
            write(*,"(A,I0)") " adjoint ddpcm step iterations: ", &
                & niter
        end if
        call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % vgrid, ddx_data % vgrid_nbasis, &
            & ddx_data % s, ddx_data % nbasis, zero, ddx_data % sgrid, &
            & ddx_data % ngrid)
        call dgemm('T', 'N', ddx_data % ngrid, ddx_data % nsph, &
            & ddx_data % nbasis, one, ddx_data % vgrid, ddx_data % vgrid_nbasis, &
            & ddx_data % y, ddx_data % nbasis, zero, ddx_data % ygrid, &
            & ddx_data % ngrid)
        ddx_data % g = ddx_data % phi - ddx_data % phieps
        ddx_data % q = ddx_data % s - fourpi/(ddx_data % eps-one)*ddx_data % y
        ddx_data % qgrid = ddx_data % sgrid - &
            & fourpi/(ddx_data % eps-one)*ddx_data % ygrid
        allocate(vsin(ddx_data % lmax+1), vcos(ddx_data % lmax+1), &
            & vplm(ddx_data % nbasis), basloc(ddx_data % nbasis), &
            & dbsloc(3, ddx_data % nbasis), fx(3, ddx_data % nsph), stat=istatus)
        if (istatus.ne.0) write(6,*) 'ddpcm forces allocation failed'
        call gradr(ddx_data, force)
        do isph = 1, ddx_data % nsph
            call fdoka(ddx_data, isph, ddx_data % xs, ddx_data % sgrid(:, isph), &
                & basloc, dbsloc, vplm, vcos, vsin, force(:,isph)) 
            call fdokb(ddx_data, isph, ddx_data % xs, ddx_data % sgrid, basloc, &
                & dbsloc, vplm, vcos, vsin, force(:, isph))
            call fdoga(ddx_data, isph, ddx_data % qgrid, ddx_data % phi_grid, &
                & force(:, isph)) 
        end do
        force = -pt5 * force
        icav = 0
        do isph = 1, ddx_data % nsph
            do igrid = 1, ddx_data % ngrid
                if(ddx_data % ui(igrid, isph) .ne. zero) then
                    icav = icav + 1
                    ddx_data % zeta(icav) = -pt5 * ddx_data % wgrid(igrid) * &
                        & ddx_data % ui(igrid, isph) * ddot(ddx_data % nbasis, &
                        & ddx_data % vgrid(1, igrid), 1, &
                        & ddx_data % q(1, isph), 1)
                    force(:, isph) = force(:, isph) + &
                        & ddx_data % zeta(icav)*gradphi_cav(:, icav)
                end if
            end do
        end do
        !! Last term where we compute gradients of potential at centers of atoms
        !! spawned by intermediate zeta.
        allocate(ef(3, ddx_data % nsph))
        if(ddx_data % fmm .eq. 1) then
            !! This step can be substituted by a proper dgemm if zeta
            !! intermediate is converted from cavity points to all grid points
            !! with zero values at internal grid points
            ! P2M step
            icav = 0
            do isph = 1, ddx_data % nsph
                inode = ddx_data % snode(isph)
                ddx_data % tmp_node_m(:, inode) = zero
                do igrid = 1, ddx_data % ngrid
                    if(ddx_data % ui(igrid, isph) .eq. zero) cycle
                    icav = icav + 1
                    call fmm_p2m(ddx_data % cgrid(:, igrid), ddx_data % zeta(icav), &
                        & one, ddx_data % pm, ddx_data % vscales, one, &
                        & ddx_data % tmp_node_m(:, inode))
                end do
                ddx_data % tmp_node_m(:, inode) = ddx_data % tmp_node_m(:, inode) / &
                    & ddx_data % rsph(isph)
            end do
            ! M2M, M2L and L2L translations
            if(ddx_data % fmm_precompute .eq. 1) then
                call tree_m2m_reflection_use_mat(ddx_data, ddx_data % tmp_node_m)
                call tree_m2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_m, &
                    & ddx_data % tmp_node_l)
                call tree_l2l_reflection_use_mat(ddx_data, ddx_data % tmp_node_l)
            else
                call tree_m2m_rotation(ddx_data, ddx_data % tmp_node_m)
                call tree_m2l_rotation(ddx_data, ddx_data % tmp_node_m, &
                    & ddx_data % tmp_node_l)
                call tree_l2l_rotation(ddx_data, ddx_data % tmp_node_l)
            end if
            ! Now compute near-field FMM gradients
            ! Cycle over all spheres
            icav = 0
            ef = zero
            do isph = 1, ddx_data % nsph
                ! Cycle over all external grid points
                do igrid = 1, ddx_data % ngrid
                    if(ddx_data % ui(igrid, isph) .eq. zero) cycle
                    icav = icav + 1
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
                        ef(:, jsph) = ef(:, jsph) + &
                            & ddx_data % zeta(icav)*d/(dnorm**3)
                    end do
                end do
            end do
            ! Take into account far-field FMM gradients only if pl > 0
            if (ddx_data % pl .gt. 0) then
                tmp1 = one / sqrt(3d0) / ddx_data % vscales(1)
                do isph = 1, ddx_data % nsph
                    inode = ddx_data % snode(isph)
                    tmp2 = tmp1 / ddx_data % rsph(isph)
                    ef(3, isph) = ef(3, isph) + tmp2*ddx_data % tmp_node_l(3, inode)
                    ef(1, isph) = ef(1, isph) + tmp2*ddx_data % tmp_node_l(4, inode)
                    ef(2, isph) = ef(2, isph) + tmp2*ddx_data % tmp_node_l(2, inode)
                end do
            end if
            do isph = 1, ddx_data % nsph
                force(:, isph) = force(:, isph) + ef(:, isph)*ddx_data % charge(isph)
            end do
        ! Naive quadratically scaling implementation
        else
            ! This routines actually computes -grad, not grad
            call efld(ddx_data % ncav, ddx_data % zeta, ddx_data % ccav, &
                & ddx_data % nsph, ddx_data % csph, ef)
            do isph = 1, ddx_data % nsph
                force(:, isph) = force(:, isph) - ef(:, isph)*ddx_data % charge(isph)
            end do
        end if
        deallocate(ef)
    end if
end subroutine ddpcm

end module ddx
