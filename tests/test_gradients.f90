!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/ddx_driver.f90
!! Main ddx driver.
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

!> Standalone application of ddX: this program can compute the solvation
!! energy and forces for a solute made of point charges. All the
!! relevant steps are here outlined.
program test_gradients
use ddx
use ddx_multipolar_solutes
use omp_lib
implicit none

real(dp), parameter :: step = 1d-6
character(len=255) :: dummy_file_name = ''
character(len=255) :: fname
character(len=2047) :: banner
type(ddx_type) :: ddx_data
type(ddx_state_type) :: state
type(ddx_error_type) :: ddx_error
type(ddx_electrostatics_type) :: electrostatics
real(dp), allocatable :: psi(:, :), force(:, :), charges(:), &
    & multipoles(:, :), dr(:), numforce(:,:), numdr(:)
real(dp) :: tol, esolv, start_time, finish_time
integer :: i, isph, info, j
real(dp), allocatable :: x(:,:), s(:,:), tmp_lx(:,:)
real(dp) :: xpsi, slx, sphi, diff
real(dp), allocatable :: grad_xpsi(:,:), grad_slx(:,:), &
    & grad_sphi(:,:), grad_xpsi_num(:,:), grad_slx_num(:,:), &
    & grad_sphi_num(:,:), dr_xpsi(:), dr_slx(:), dr_sphi(:), &
    & dr_xpsi_num(:), dr_slx_num(:), dr_sphi_num(:), scratch_force(:, :)
real(dp), external :: ddot
real(dp), parameter :: threshold = 1e-8


call get_command_argument(1, fname)
call ddfromfile(fname, ddx_data, tol, charges, ddx_error)
call check_error(ddx_error)

call allocate_state(ddx_data % params, ddx_data % constants, state, &
    & ddx_error)
call check_error(ddx_error)

allocate( &
    & grad_xpsi(3,ddx_data%params%nsph), &
    & grad_slx(3,ddx_data%params%nsph), &
    & grad_sphi(3,ddx_data%params%nsph), &
    & grad_xpsi_num(3,ddx_data%params%nsph), &
    & grad_slx_num(3,ddx_data%params%nsph), &
    & grad_sphi_num(3,ddx_data%params%nsph), &
    & dr_xpsi(ddx_data%params%nsph), &
    & dr_slx(ddx_data%params%nsph), &
    & dr_sphi(ddx_data%params%nsph), &
    & dr_xpsi_num(ddx_data%params%nsph), &
    & dr_slx_num(ddx_data%params%nsph), &
    & dr_sphi_num(ddx_data%params%nsph), &
    & x(ddx_data%constants%nbasis,ddx_data%params%nsph), &
    & s(ddx_data%constants%nbasis,ddx_data%params%nsph), &
    & tmp_lx(ddx_data%constants%nbasis,ddx_data%params%nsph), &
    & force(3, ddx_data % params % nsph), &
    & scratch_force(3, ddx_data % params % nsph), &
    & dr(ddx_data % params % nsph), &
    & numforce(3, ddx_data % params % nsph), &
    & numdr(ddx_data % params % nsph), &
    & multipoles(1, ddx_data % params % nsph), stat=info)
if (info .ne. 0) then
    write(6, *) "Allocation failed in ddx_driver"
    stop 1
end if

multipoles(1, :) = charges/sqrt4pi
call multipole_electrostatics(ddx_data % params, ddx_data % constants, &
    & ddx_data % workspace, multipoles, 0, electrostatics, ddx_error)

allocate(psi(ddx_data % constants % nbasis, ddx_data % params % nsph), &
    & stat=info)
if (info .ne. 0) then
    write(6, *) "Allocation failed in ddx_driver"
    stop 1
end if
call multipole_psi(ddx_data % params, multipoles, 0, psi)


! Reference calculation
! call draco(ddx_data%params%csph, ddx_data%params%nsph, &
   !  & ddx_data%params%rsph)
call ddrun(ddx_data,state,electrostatics,psi,tol,esolv,ddx_error,force)


x(:,:) = state%xs(:,:)
s(:,:) = state%s(:,:)
! Compute the first Lagrangian term 0.5 X^T Psi
xpsi = pt5*ddot(ddx_data%constants%n,x,1,psi,1)

! Compute the second Lagrangian term -0.5 S^T L X
ddx_data%constants%dodiag = .true.
call lx(ddx_data%params,ddx_data%constants,ddx_data%workspace, &
    & x,tmp_lx,ddx_error)
slx = -pt5*ddot(ddx_data%constants%n,s,1,tmp_lx,1)

! Compute the third Lagrangian term - 0.5 S^T Phi
! Note that Phi is saved as -Phi, so no minus required
sphi = pt5*ddot(ddx_data%constants%n,s,1,state%phi,1)

write(6,*) xpsi, slx, sphi
! Check if the Lagrangian formulation returns the energy
diff = abs(esolv - (xpsi + slx + sphi))
if (diff.gt.threshold) then
    write(6, *) "Inconsistency"
    stop 1
end if


call sgradlx(ddx_data%params,ddx_data%constants,ddx_data%workspace, &
    & state,grad_slx,ddx_error)

call sgradphi(ddx_data%params,ddx_data%constants,ddx_data%workspace, &
    & state,grad_sphi,ddx_error,electrostatics%e_cav,multipoles)

scratch_force = zero
call sdrlx(ddx_data%params,ddx_data%constants,ddx_data%workspace, &
    & state,scratch_force,dr_slx,ddx_error)

scratch_force = zero
call sdrphi(ddx_data%params,ddx_data%constants,ddx_data%workspace, &
    & state,scratch_force,dr_sphi,ddx_error,electrostatics%e_cav,multipoles)

do isph = 1, ddx_data%params%nsph
    do j = 1, 3
        ddx_data%params%csph(j,isph) = ddx_data%params%csph(j,isph) - step
        call displaced_run(ddx_data,multipoles,tol,xpsi,slx,sphi,psi, &
             & esolv,force,x,s)

        grad_xpsi_num(j, isph) = - xpsi
        grad_slx_num(j, isph) = - slx
        grad_sphi_num(j, isph) = - sphi

        ddx_data%params%csph(j,isph) = ddx_data%params%csph(j,isph) + 2.0d0*step
        call displaced_run(ddx_data,multipoles,tol,xpsi,slx,sphi,psi, &
            & esolv,force,x,s)

        grad_xpsi_num(j, isph) = grad_xpsi_num(j, isph) + xpsi
        grad_slx_num(j, isph) = grad_slx_num(j, isph) + slx
        grad_sphi_num(j, isph) = grad_sphi_num(j, isph) + sphi
        grad_xpsi_num(j, isph) = grad_xpsi_num(j, isph)/(2.0d0*step)
        grad_slx_num(j, isph) = grad_slx_num(j, isph)/(2.0d0*step)
        grad_sphi_num(j, isph) = grad_sphi_num(j, isph)/(2.0d0*step)

        ddx_data%params%csph(j,isph) = ddx_data%params%csph(j,isph) - step
    end do
end do

diff = maxval(abs(grad_slx_num - grad_slx))
write(6,*) "Difference S^T grad L X", diff
if (diff.gt.threshold) then
    write(6, *) "Inconsistency"
    stop 1
end if

diff = maxval(abs(grad_sphi_num - grad_sphi))
write(6,*) "Difference S^T grad Phi", diff
if (diff.gt.threshold) then
    write(6, *) "Inconsistency"
    stop 1
end if


do isph = 1, ddx_data%params%nsph
        ddx_data%params%rsph(isph) = ddx_data%params%rsph(isph) - step
        call displaced_run(ddx_data,multipoles,tol,xpsi,slx,sphi,psi, &
             & esolv,force,x,s)

        dr_xpsi_num( isph) = - xpsi
        dr_slx_num(isph) = - slx
        dr_sphi_num( isph) = - sphi

        ddx_data%params%rsph(isph) = ddx_data%params%rsph(isph) + 2.0d0*step
        call displaced_run(ddx_data,multipoles,tol,xpsi,slx,sphi,psi, &
            & esolv,force,x,s)

        dr_xpsi_num( isph) = dr_xpsi_num( isph) + xpsi
        dr_slx_num(isph) = dr_slx_num(isph) + slx
        dr_sphi_num( isph) = dr_sphi_num(isph) + sphi
        dr_xpsi_num( isph) = dr_xpsi_num(isph)/(2.0d0*step)
        dr_slx_num(isph) = dr_slx_num(isph)/(2.0d0*step)
        dr_sphi_num( isph) = dr_sphi_num(isph)/(2.0d0*step)

        ddx_data%params%rsph(isph) = ddx_data%params%rsph(isph) - step
end do

write(6,*) "-----"
diff = maxval(abs(dr_slx_num - dr_slx))
write(6,*) "Difference S^T dr L X", diff
if (diff.gt.threshold) then
    write(6, *) "Inconsistency"
    stop 1
end if
diff = maxval(abs(dr_sphi_num - dr_sphi))
write(6,*) "Difference S^T dr Phi", diff
if (diff.gt.threshold) then
    write(6, *) "Inconsistency"
    stop 1
end if

deallocate(psi, multipoles, charges, force, scratch_force, &
    & numforce, dr, numdr, &
    & x, s, tmp_lx, &
    & grad_xpsi, grad_slx, grad_sphi, &
    & grad_xpsi_num, grad_slx_num, grad_sphi_num, &
    & dr_xpsi, dr_slx, dr_sphi, &
    & dr_xpsi_num, dr_slx_num, dr_sphi_num, &
    & stat=info)
if (info .ne. 0) then
    write(6, *) "Deallocation failed in ddx_driver"
    stop 1
end if
call deallocate_electrostatics(electrostatics, ddx_error)
call deallocate_state(state, ddx_error)
call deallocate_model(ddx_data, ddx_error)

contains

subroutine displaced_run(ddx_data,multipoles,tol,xpsi,slx,sphi, &
        & psi,esolv,force,x,s)
    type(ddx_type), intent(inout) :: ddx_data
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: xpsi, slx, sphi
    type(ddx_type) :: ddx_data2
    type(ddx_error_type) :: error2
    type(ddx_state_type) :: state2
    type(ddx_electrostatics_type) :: electrostatics2
    real(dp), intent(in) :: multipoles(1,ddx_data%params%nsph)
    real(dp), intent(out) :: esolv
    real(dp), intent(out) :: force(3,ddx_data%params%nsph)
    real(dp), intent(out) :: psi(ddx_data%constants%nbasis, &
        & ddx_data%params%nsph)
    real(dp), intent(in) :: &
        & x(ddx_data%constants%nbasis,ddx_data%params%nsph), &
        & s(ddx_data%constants%nbasis,ddx_data%params%nsph)

    ! Make a copy of the model
    call allocate_model(ddx_data%params%nsph,ddx_data%params%csph(1,:), &
        & ddx_data%params%csph(2,:),ddx_data%params%csph(3,:), &
        & ddx_data%params%rsph,ddx_data%params%model, &
        & ddx_data%params%lmax,ddx_data%params%ngrid,ddx_data%params%force, &
        & ddx_data%params%fmm,ddx_data%params%pm,ddx_data%params%pl, &
        & ddx_data%params%se,ddx_data%params%eta,ddx_data%params%eps, &
        & ddx_data%params%kappa,ddx_data%params%matvecmem, &
        & ddx_data%params%maxiter,ddx_data%params%jacobi_ndiis, &
        & ddx_data%params%nproc,dummy_file_name,ddx_data2,error2)
    call check_error(error2)

    call allocate_state(ddx_data2%params,ddx_data2%constants,state2, &
        & error2)
    call check_error(error2)


    call multipole_electrostatics(ddx_data2%params,ddx_data2%constants, &
        & ddx_data2%workspace,multipoles,0,electrostatics2,error2)
    call check_error(error2)

    call multipole_psi(ddx_data2%params,multipoles,0,psi)
    call ddrun(ddx_data2,state2,electrostatics2,psi,tol,esolv,ddx_error,force)

    ! Compute the first Lagrangian term 0.5 X^T Psi
    xpsi = pt5*ddot(ddx_data2%constants%n,x,1,psi,1)
    ! Compute the second Lagrangian term -0.5 S^T L X
    ddx_data2%constants%dodiag = .true.
    call lx(ddx_data2%params,ddx_data2%constants,ddx_data2%workspace, &
        & x,tmp_lx,error2)
    slx = -pt5*ddot(ddx_data2%constants%n,s,1,tmp_lx,1)
    ! Compute the third Lagrangian term - 0.5 S^T Phi
    ! Note that Phi is saved as -Phi, so no minus required
    sphi = pt5*ddot(ddx_data2%constants%n,s,1,state2%phi,1)

    call deallocate_electrostatics(electrostatics2,error2)
    call deallocate_state(state2,error2)
    call deallocate_model(ddx_data2,error2)
end subroutine displaced_run

subroutine sgradlx(params, constants, workspace, &
        & state, force, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(inout) :: force(3, params % nsph)
    integer :: isph
    force = zero
    do isph = 1, params % nsph
        call contract_grad_l(params,constants,isph,state%xs, &
            & state%sgrid,workspace%tmp_vylm(:,1), &
            & workspace%tmp_vdylm(:,:,1),workspace%tmp_vplm(:,1), &
            & workspace%tmp_vcos(:,1),workspace%tmp_vsin(:,1), &
            & force(:,isph))
    end do
    force = -pt5*force
end subroutine sgradlx

subroutine sdrlx(params, constants, workspace, &
        & state, force, dr, ddx_error)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(inout) :: force(3, params % nsph)
    real(dp), intent(inout) :: dr(params % nsph)
    integer :: isph
    force = zero
    dr = zero
    do isph = 1, params % nsph
        call contract_grad_l(params,constants,isph,state%xs, &
            & state%sgrid,workspace%tmp_vylm(:,1), &
            & workspace%tmp_vdylm(:,:,1),workspace%tmp_vplm(:,1), &
            & workspace%tmp_vcos(:,1),workspace%tmp_vsin(:,1), &
            & force(:,isph), dr=dr(isph))
    end do
    dr = -pt5*dr
    force = -pt5*force
end subroutine sdrlx

subroutine sgradphi(params, constants, workspace, &
        & state,force,ddx_error,e_cav,multipoles)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(inout) :: force(3, params % nsph)
    real(dp), intent(in) :: e_cav(3, constants % ncav)
    real(dp), intent(in) :: multipoles(1, ddx_data % params % nsph)
    integer :: isph
    force = zero
    do isph = 1, params % nsph
        call contract_grad_u(params, constants, isph, state % sgrid, &
            & state % phi_grid, force(:, isph))
    end do
    force = -pt5*force
    call zeta_grad(params, constants, state, e_cav, force)
    call multipole_force_terms(params,constants, &
        & workspace,state,0,multipoles,force,ddx_error)
end subroutine sgradphi


subroutine sdrphi(params, constants, workspace, &
        & state,force,dr,ddx_error,e_cav,multipoles)
    implicit none
    type(ddx_params_type), intent(in) :: params
    type(ddx_constants_type), intent(in) :: constants
    type(ddx_workspace_type), intent(inout) :: workspace
    type(ddx_state_type), intent(inout) :: state
    type(ddx_error_type), intent(inout) :: ddx_error
    real(dp), intent(inout) :: force(3, params % nsph)
    real(dp), intent(inout) :: dr(params % nsph)
    real(dp), intent(in) :: e_cav(3, constants % ncav)
    real(dp), intent(in) :: multipoles(1, ddx_data % params % nsph)
    integer :: isph
    dr = zero
    do isph = 1, params % nsph
        call contract_grad_u(params, constants, isph, state % sgrid, &
            & state % phi_grid, force(:, isph), dr=dr(isph))
    end do
    dr = -pt5*dr
    call zeta_grad_dr(params, constants, state, e_cav, dr)
end subroutine sdrphi

end program test_gradients
