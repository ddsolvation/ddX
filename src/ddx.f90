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

!> High-level module of the ddX software
module ddx
! Get ddcosmo-module
use ddx_cosmo
! Get ddpcm-module
use ddx_pcm
! Get ddlpb-module
use ddx_lpb
implicit none

contains

!> Main solver routine
!!
!! Solves the problem within COSMO model using a domain decomposition approach.
!!
!! @param[in] ddx_data: ddX object with all input information
!! @param[in] phi_cav: Potential at cavity points
!! @param[in] gradphi_cav: Gradient of a potential at cavity points
!! @param[in] psi: RHS of the adjoint problem
!! @param[in] tol
!! @param[out] esolv: Solvation energy
!! @param[out] force: Analytical forces
subroutine ddsolve(ddx_data, state, phi_cav, gradphi_cav, hessianphi_cav, &
        & psi, tol, esolv, force)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    type(ddx_state_type), intent(inout) :: state
    real(dp), intent(in) :: phi_cav(ddx_data % constants % ncav), &
        & gradphi_cav(3, ddx_data % constants % ncav), &
        & hessianphi_cav(3, ddx_data % constants % ncav), &
        & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), tol
    ! Outputs
    real(dp), intent(out) :: esolv, force(3, ddx_data % params % nsph)
    ! Find proper model
    select case(ddx_data % params % model)
        ! COSMO model
        case (1)
            call ddcosmo(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, phi_cav, psi, &
                & tol, esolv, force)
        ! PCM model
        case (2)
            call ddpcm(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, phi_cav, psi, &
                & tol, esolv, force)
        ! LPB model
        case (3)
            call ddlpb(ddx_data % params, ddx_data % constants, &
                & ddx_data % workspace, state, phi_cav, gradphi_cav, &
                & hessianphi_cav, psi, tol, esolv, force)
        ! Error case
        case default
            ddx_data % params % error_flag = 1
            ddx_data % params % error_message = "unsupported solvation " // &
                & " model in the dd solver."
            return
    end select
end subroutine ddsolve

end module ddx
