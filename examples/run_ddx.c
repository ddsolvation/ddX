#include "ddx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void print_error(void* error) {
  char message[2000];
  ddx_get_error_message(error, message, 2000);
  printf("Model error: %s\n", message);
}

int main() {
  int n_spheres      = 12;
  double charges[12] = {-0.04192, -0.04192, -0.04198, -0.04192, -0.04192, -0.04198,
                        0.04193,  0.04193,  0.04197,  0.04193,  0.04193,  0.04197};

  // Radii in Aangstroem
  double radii[12] = {4.00253, 4.00253, 4.00253, 4.00253, 4.00253, 4.00253,
                      2.99956, 2.99956, 2.99956, 2.99956, 2.99956, 2.99956};
  for (int i = 0; i < 12; ++i) {
    radii[i] /= 0.5291772109;
  }

  // Coordinates in Aangstroem
  double centres[36];
  double x[12] = {0.00000, 0.00000, 0.00000, 0.00000,  0.00000,  0.00000,
                  0.00103, 0.00103, 0.00000, -0.00103, -0.00103, 0.00000};
  double y[12] = {2.29035, 2.29035, 0.00000, -2.29035, -2.29035, 0.00000,
                  4.05914, 4.05914, 0.00000, -4.05914, -4.05914, 0.00000};
  double z[12] = {1.32281, -1.32281, -2.64562, -1.32281, 1.32281, 2.64562,
                  2.34326, -2.34326, -4.68652, -2.34326, 2.34326, 4.68652};
  for (int i = 0; i < 12; ++i) {
    centres[3 * i + 0] = x[i] / 0.5291772109;
    centres[3 * i + 1] = y[i] / 0.5291772109;
    centres[3 * i + 2] = z[i] / 0.5291772109;
  }

  char banner[2048];
  ddx_get_banner(banner, 2048);
  printf("%s\n", banner);

  int pcm                = 2;
  int enable_forces      = 1;
  double epsilon         = 78.3553;
  double kappa           = 0.0;
  double eta             = 0.1;
  double shift           = 0.0;
  int lmax               = 8;
  int n_lebedev          = 302;
  int incore             = 0;
  int maxiter            = 100;
  int jacobi_n_diis      = 20;
  int enable_fmm         = 1;
  int fmm_multipole_lmax = 7;
  int fmm_local_lmax     = 6;
  int n_proc             = 1;
  int length_logfile     = 0;
  char logfile[1]        = {'\0'};

  void* error = ddx_allocate_error();

  //
  // Allocate model and print some info
  //
  void* model = ddx_allocate_model(pcm, enable_forces, epsilon, kappa, eta, shift, lmax,
                                   n_lebedev, incore, maxiter, jacobi_n_diis, enable_fmm,
                                   fmm_multipole_lmax, fmm_local_lmax, n_proc, n_spheres,
                                   centres, radii, length_logfile, logfile, error);
  if (ddx_get_error_flag(error) != 0) {
    print_error(error);
    return 1;
  }

  int nsph   = ddx_get_n_spheres(model);
  int nbasis = ddx_get_n_basis(model);
  int ncav   = ddx_get_n_cav(model);
  printf("nsph   = %4d\n", nsph);
  printf("nbasis = %4d\n", nbasis);
  printf("ncav   = %4d\n", ncav);


  // convert the charges to multipoles
  int nmultipoles = 1;
  double* solute_multipoles = (double*)malloc(sizeof(double) * nmultipoles * nsph);
  for (int i = 0; i < nmultipoles * nsph; ++i) {
    solute_multipoles[i] = charges[i] / sqrt(4.0 * M_PI);
  }


  // Compute the RHSs
  void* electrostatics = ddx_allocate_electrostatics(model, error);
  ddx_multipole_electrostatics(model, nsph, nmultipoles, solute_multipoles,
                               electrostatics, error);

  double* psi = (double*)malloc(sizeof(double) * nsph * nbasis);
  ddx_multipole_psi(model, nbasis, nsph, nmultipoles, solute_multipoles, psi, error);
  if (ddx_get_error_flag(error) != 0) {
    print_error(error);
    return 1;
  }

  //
  // Solve the PCM problem
  //
  void* state = ddx_allocate_state(model, error);
  if (ddx_get_error_flag(error) != 0) {
    print_error(state);
    return 1;
  }

  double tol = 1e-9;
  ddx_setup(model, state, electrostatics, nbasis, nsph, psi, error);
  ddx_fill_guess(model, state, tol, error);
  ddx_solve(model, state, tol, error);
  if (ddx_get_error_flag(error) != 0) {
    print_error(error);
    return 1;
  }
  printf("Forward system solved after %i iterations.\n", ddx_get_x_niter(state));

  ddx_fill_guess_adjoint(model, state, tol, error);
  ddx_solve_adjoint(model, state, tol, error);
  if (ddx_get_error_flag(error) != 0) {
    print_error(error);
    return 1;
  }
  printf("Adjoint system solved after %i iterations.\n", ddx_get_s_niter(state));

  //
  // Compute energy and forces
  //
  double energy = ddx_pcm_energy(model, state, error);

  double* forces = (double*)malloc(sizeof(double) * 3 * nsph);
  ddx_solvation_force_terms(model, state, electrostatics, nsph, forces, error);
  if (ddx_get_error_flag(error) != 0) {
    print_error(error);
    return 1;
  }

  ddx_multipole_force_terms(model, state, nsph, 0, nmultipoles, solute_multipoles,
                            forces, error);
  if (ddx_get_error_flag(error) != 0) {
    print_error(error);
    return 1;
  }

  //
  // Check results
  //
  int ret     = 0;
  double eref = -0.00017974013712832552;
  if (fabs(eref - energy) > 1e-6) {
    printf("Large deviation:   %15.9g\n", eref - energy);
    ret = 1;
  }

  // Finish up
  free(forces);
  ddx_deallocate_state(state, error);
  free(psi);
  free(solute_multipoles);
  ddx_deallocate_model(model, error);
  ddx_deallocate_electrostatics(electrostatics, error);
  return ret;
}
