#include "ddx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void print_model_error(void* model) {
  char message[256];
  ddx_get_error_message(model, message, 256);
  printf("Model error: %s\n", message);
}
void print_state_error(void* state) {
  char message[256];
  ddx_get_state_error_message(state, message, 256);
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
  int info               = 0;
  int length_logfile     = 0;
  char logfile[1]        = {'\0'};

  //
  // Allocate model and print some info
  //
  void* model = ddx_allocate_model(pcm, enable_forces, epsilon, kappa, eta, shift, lmax,
                                   n_lebedev, incore, maxiter, jacobi_n_diis, enable_fmm,
                                   fmm_multipole_lmax, fmm_local_lmax, n_proc, n_spheres,
                                   charges, centres, radii, length_logfile, logfile);
  if (ddx_get_error_flag(model) != 0) {
    print_model_error(model);
    return 1;
  }

  int nsph   = ddx_get_n_spheres(model);
  int nbasis = ddx_get_n_basis(model);
  int ncav   = ddx_get_n_cav(model);
  printf("nsph   = %4d\n", nsph);
  printf("nbasis = %4d\n", nbasis);
  printf("ncav   = %4d\n", ncav);

  //
  // Compute solute electrostatics
  //
  double* phi_cav = (double*)malloc(sizeof(double) * ncav);
  double* e_cav   = (double*)malloc(sizeof(double) * 3 * ncav);

  int nmultipoles           = 1;  // Just charges
  double* solute_multipoles = (double*)malloc(sizeof(double) * nmultipoles * nsph);
  for (int i = 0; i < nmultipoles * nsph; ++i) {
    solute_multipoles[i] = charges[i] / sqrt(4.0 * M_PI);
  }
  ddx_multipole_electrostatics_1(model, nsph, ncav, nmultipoles, solute_multipoles,
                                 phi_cav, e_cav);
  if (ddx_get_error_flag(model) != 0) {
    print_model_error(model);
    return 1;
  }

  double* psi = (double*)malloc(sizeof(double) * nsph * nbasis);
  ddx_multipole_psi(model, nbasis, nsph, nmultipoles, solute_multipoles, psi);
  if (ddx_get_error_flag(model) != 0) {
    print_model_error(model);
    return 1;
  }

  //
  // Solve the PCM problem
  //
  void* state = ddx_allocate_state(model);
  if (ddx_get_state_error_flag(state) != 0) {
    print_state_error(state);
    return 1;
  }

  double tol = 1e-9;
  ddx_pcm_setup(model, state, ncav, nbasis, nsph, phi_cav, psi);
  ddx_pcm_guess(model, state);
  ddx_pcm_solve(model, state, tol);
  if (ddx_get_error_flag(model) != 0) {
    print_model_error(model);
    return 1;
  }
  printf("Forward system solved after %i iterations.\n", ddx_get_x_niter(state));

  ddx_pcm_guess_adjoint(model, state);
  ddx_pcm_solve_adjoint(model, state, tol);
  if (ddx_get_error_flag(model) != 0) {
    print_model_error(model);
    return 1;
  }
  printf("Adjoint system solved after %i iterations.\n", ddx_get_s_niter(state));

  //
  // Compute energy and forces
  //
  double* solution_x = (double*)malloc(sizeof(double) * nsph * nbasis);
  ddx_get_x(state, nbasis, nsph, solution_x);
  double energy = 0.0;
  for (int i = 0; i < nsph * nbasis; ++i) {
    energy += 0.5 * solution_x[i] * psi[i];
  }

  double* forces = (double*)malloc(sizeof(double) * 3 * nsph);
  ddx_pcm_solvation_force_terms(model, state, nsph, forces);
  if (ddx_get_error_flag(model) != 0) {
    print_model_error(model);
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
  free(solution_x);
  ddx_deallocate_state(state);
  free(psi);
  free(solute_multipoles);
  free(e_cav);
  free(phi_cav);
  ddx_deallocate_model(model);
  return ret;
}
