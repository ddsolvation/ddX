!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/main_ddLPB.f90
!! Main ddx driver.
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @author Abhinav Jha
!! @date 2021-02-11

program main
use ddx_core
use ddx_operators
use ddx_solvers
use ddx_lpb
implicit none

character(len=255) :: fname
type(ddx_type) :: ddx_data
integer :: iprint, nproc, lmax, pmax, ngrid, iconv, igrad, n, force, fmm, model
integer :: niter, ndiis=25, fmm_precompute, itersolver, maxiter
logical :: ok
real(dp) :: eps, eta, tol, se, kappa
! esolv       : Electrostatic Solvation Energy
real(dp) :: esolv
!
! Quantities to be allocated by the user
! - Solute's parameters
!   x, y, z   : Coordinates of the spheres
!   rvdw      : Van Der Waal(VDW) radii
!   charge    : Point charges used to model the solute (or multipoles,
!               or qm density...)
!
real(dp), allocatable :: x(:), y(:), z(:), rvdw(:), charge(:)
!
! - Electrostatic potential 
!   phi         : Boundayry conditions (This is psi_0 Eq.(20) QSM19.SISC)
!   gradphi     : Gradient of the phi vector
!   psi         : Electrostatic potential vector of size nylm*n
!   xs          : (?)
!
real(dp), allocatable :: phi(:), gradphi(:,:), psi(:, :), xs(:, :)
real(dp), allocatable :: g(:, :), rhs(:, :)
!
! These constants are defined in ddX library already
! toang       : Conversion for Angstrom
! tokcal      : Conversion for Energy (?)
! tobohr      : Conversion from Angstrom to Bohr radius
!
!real(dp), parameter :: toang=0.52917721092d0, tokcal=627.509469d0
!real(dp), parameter :: tobohr=1d0/toang
!
! - ddLPB solution 
!   sigma (nylm,n)      : Solution of ddLPB
!
real*8, allocatable :: sigma(:,:)

integer :: i, j, info

!
! Here, we read all the ddCOSMO parameters from a file named Input.txt
! If one wants to use another input file, change file = '_', to the said file.
!
call getarg(1, fname)
write(*, *) "Reading input file ", fname
open (unit=100,file=fname,form='formatted',access='sequential')
!
! Scalar parameters. The variables are defined in the ddCOSMO module and are common to
! all the ddCOSMO routines (no need to declare them if ddcosmo.mod is loaded.)
! In ddLPB one more parameter is required, namely the Debye-Hückel constant (kappa)
! NOTE: Change the format so that epsilon_1, i.e., dielectric constant for solute
!       is taken.
!
read(100,*) iprint      ! printing flag
read(100,*) nproc       ! number of openmp threads
read(100,*) lmax        ! max angular momentum of spherical harmonics basis
read(100,*) pmax        ! max degree of harmonics for the FMM
read(100,*) ngrid       ! number of lebedev points
read(100,*) iconv       ! 10^(-iconv) is the convergence threshold for the iterative solver
read(100,*) igrad       ! whether to compute (1) or not (0) forces
read(100,*) eps         ! dielectric constant of the solvent
read(100,*) eta         ! regularization parameter
read(100,*) kappa       ! Debye-Hückel constant
!
read(100,*) n           ! number of atoms
!
allocate(x(n), y(n), z(n), rvdw(n), charge(n))
!
! We also read from the same file the charges, coordinates and vdw radii.
! In this example, the coordinates and radii are read in Angstrom and
! converted in Bohr before calling ddinit.
!
do i = 1, n
  read(100,*) charge(i), x(i), y(i), z(i), rvdw(i)
end do
!
! Conversion to Bohr radius
!
x    = x*tobohr
y    = y*tobohr
z    = z*tobohr
rvdw = rvdw*tobohr
!
close (100)
!
! model : 1 for COSMO, 2 for PCM, 3 for LPB
! force : 1 for forces, 0 for not
! fmm   : 1 for forces, 0 for not
!
model=3
force=0
fmm=0
se=-one
itersolver=1
tol=1d-1**iconv
maxiter=200
call ddinit(n, charge, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pmax, pmax, &
    & fmm_precompute, iprint, se, eta, eps, kappa, itersolver, tol, maxiter, &
    & ndiis, nproc, ddx_data, info)

allocate(phi(ddx_data % ncav), psi(ddx_data % nbasis,n), gradphi(3, ddx_data % ncav))

call mkrhs(ddx_data, phi, gradphi, psi)

niter = 200
! Now, call the ddLPB solver
!
allocate (sigma(ddx_data % nbasis ,n))
!
! @param[in] phi      : Boundary conditions
! @param[in] charge   : Charge of atoms
! @param[in] psi      : Electrostatic potential vector of size nylm*n
!                       Not sure where it is getting used (?)
! @param[in] gradphi  : Gradient of phi
!
! @param[out] sigma   : Solution of ddLPB
! @param[out] esolv   : Electrostatic solvation energy
!
call ddlpb(ddx_data, phi, psi, gradphi, sigma, esolv, charge, ndiis, niter, iconv)
!call cosmo(.false., .true., phi, xx, psi, sigma, esolv)
!
if (iprint.ge.3) call prtsph('Solution to the ddLPB equation',ddx_data % nbasis, ddx_data % lmax, ddx_data % nsph, 0, sigma)
!
write (6,'(1x,a,f14.6)') 'ddLPB Electrostatic Solvation Energy (kcal/mol):', esolv*tokcal
deallocate(phi, psi, gradphi)
deallocate(x, y, z, rvdw, charge)
call ddfree(ddx_data)

end program main

