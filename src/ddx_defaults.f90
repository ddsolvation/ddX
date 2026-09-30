!> Module defining default parameters, for consistency between the
!! various initialization routines.
module ddx_defaults
use ddx_definitions, only: dp
implicit none

integer, parameter :: default_force = 0
integer, parameter :: default_adjoint = 0
integer, parameter :: default_lmax = 6
integer, parameter :: default_ngrid = 302
integer, parameter :: default_incore = 0
integer, parameter :: default_maxiter = 100
integer, parameter :: default_jacobi_ndiis = 20
integer, parameter :: default_enable_fmm = 1
integer, parameter :: default_pm = 8
integer, parameter :: default_pl = 8
integer, parameter :: default_nproc = 1
integer, parameter :: default_switching = 0

real(dp), parameter :: default_kappa = 0.0d0
real(dp), parameter :: default_eta = 0.1d0
! the default value of the shift depends on the model for ddCOSMO
! we use an internal switching, for ddPCM and ddLPB a symmetric
! switching
real(dp), parameter :: default_se(3) = (/-1.0d0, 0.0d0, 0.0d0/)
real(dp), parameter :: default_eps_int = 1.0d0
real(dp), parameter :: default_tol = 1d-6

character(len=255), parameter :: default_logfile = ""

end module ddx_defaults
