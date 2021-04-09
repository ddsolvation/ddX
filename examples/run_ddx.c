#include "ddx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
  double charges[12] = {-0.04192, -0.04192, -0.04198, -0.04192, -0.04192, -0.04198,
                        0.04193,  0.04193,  0.04197,  0.04193,  0.04193,  0.04197};

  double radii[12] = {7.56368551, 7.56368551, 7.56368551, 7.56368551,
                      7.56368551, 7.56368551, 5.66834689, 5.66834689,
                      5.66834689, 5.66834689, 5.66834689, 5.66834689};

  double x[12] = {0.00000, 0.00000, 0.00000, 0.00000,  0.00000,  0.00000,
                  0.00103, 0.00103, 0.00000, -0.00103, -0.00103, 0.00000};
  double y[12] = {2.29035, 2.29035, 0.00000, -2.29035, -2.29035, 0.00000,
                  4.05914, 4.05914, 0.00000, -4.05914, -4.05914, 0.00000};
  double z[12] = {1.32281, -1.32281, -2.64562, -1.32281, 1.32281, 2.64562,
                  2.34326, -2.34326, -4.68652, -2.34326, 2.34326, 4.68652};
  for (int i = 0; i < 12; ++i) {
    x[i] *= 1 / 0.5291772109;
    y[i] *= 1 / 0.5291772109;
    z[i] *= 1 / 0.5291772109;
  }

  int pcm                = 2;
  int lmax               = 10;
  int n_lebedev          = 302;
  int compute_force      = 1;
  int fmm                = 1;
  int fmm_multipole_lmax = 20;
  int fmm_local_lmax     = 20;
  double se              = 0.0;
  double eta             = 0.1;
  double epsilon         = 78.3553;
  double kappa           = 0.0;
  int nproc              = 1;
  int info               = 0;
  void* model = ddx_init(12, charges, x, y, z, radii, pcm, lmax, n_lebedev, compute_force,
                         fmm, fmm_multipole_lmax, fmm_local_lmax, se, eta, epsilon, kappa,
                         nproc, &info);

  int nsph   = ddx_get_nsph(model);
  int nbasis = ddx_get_nbasis(model);
  int ncav   = ddx_get_ncav(model);
  printf("nsph   = %4d\n", nsph);
  printf("nbasis = %4d\n", nbasis);
  printf("ncav   = %4d\n", ncav);

  double* phi_cav     = (double*)malloc(sizeof(double) * ncav);
  double* gradphi_cav = (double*)malloc(sizeof(double) * 3 * ncav);
  double* psi         = (double*)malloc(sizeof(double) * nsph * nbasis);
  ddx_mkrhs(model, nsph, ncav, nbasis, phi_cav, gradphi_cav, psi);

  int itersolver = 1;
  double tol     = 1e-9;
  int maxiter    = 100;
  int ndiis      = 20;
  int iprint     = 1;
  double esolve  = 0;
  double* force  = (double*)malloc(sizeof(double) * 3 * nsph);
  ddx_solve(model, nsph, ncav, nbasis, phi_cav, gradphi_cav, psi, itersolver, tol,
            maxiter, ndiis, iprint, &esolve, force);

  int ret     = 0;
  double eref = -0.00017974013712832552;
  if (fabs(eref - esolve) > 1e-8) {
    printf("Large deviation:   %15.9g\n", eref - esolve);
    ret = 1;
  }

  // free(force);
  free(phi_cav);
  free(gradphi_cav);
  free(psi);
  ddx_finish(model);
  return ret;
}
