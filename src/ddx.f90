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
use ddx_cosmo
use ddx_pcm
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
!! @param[in] psi: TODO
!! @param[in] tol
!! @param[out] esolv: Solvation energy
!! @param[out] force: Analytical forces
subroutine ddsolve(ddx_data, phi_cav, gradphi_cav, hessianphi_cav, psi, tol, &
        & esolv, force, info)
    ! Inputs
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(in) :: phi_cav(ddx_data % constants % ncav), &
        & gradphi_cav(3, ddx_data % constants % ncav), &
        & hessianphi_cav(3, ddx_data % constants % ncav), &
        & psi(ddx_data % constants % nbasis, ddx_data % params % nsph), tol
    real(dp) :: psi_lpb(ddx_data % constants % nbasis, ddx_data % params % nsph)
    ! Outputs
    real(dp), intent(out) :: esolv, force(3, ddx_data % params % nsph)
    integer, intent(out) :: info
    ! Find proper model
    select case(ddx_data % params % model)
        ! COSMO model
        case (1)
            call ddcosmo(ddx_data, phi_cav, gradphi_cav, psi, tol, esolv, &
                & force, info)
        ! PCM model
        case (2)
            call ddpcm(ddx_data, phi_cav, gradphi_cav, psi, tol, esolv, &
                & force, info)
        ! LPB model
        case (3)
            ! Psi shall be divided by a factor 4pi for the LPB case
            ! It is intended to take into account this constant in the LPB
            psi_lpb = psi / fourpi
            call ddlpb(ddx_data, phi_cav, gradphi_cav, hessianphi_cav, psi_lpb, &
                & tol, esolv, force, info)
        ! Error case
        case default
            stop "Non-supported model"
    end select
end subroutine ddsolve

end module ddx
