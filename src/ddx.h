#ifndef DDX_H
#define DDX_H

/** C header for interfacing with ddx */

#ifdef __cplusplus
extern "C" {
#endif

//
// Generic stuff
//
int ddx_supported_lebedev_grids(int n, int* grids);

// Get scaled ylm at a point and with respect to a cavity sphere
void ddx_scaled_ylm(const void* c_ddx, int lmax, const double* x, int sphere,
                    double* ylm);

void ddx_nuclear_contributions(const void* c_ddx, int nsph, int ncav, int nbasis,
                               double* phi, double* gradphi, double* psi);

//
// Setup object
//
void* ddx_allocate_model(int model, int enable_force, double solvent_epsilon,
                         double solvent_kappa, double eta, int lmax, int n_lebedev,
                         int incore, int maxiter, int jacobi_n_diis, int enable_fmm,
                         int fmm_multipole_lmax, int fmm_local_lmax, int n_proc,
                         int n_spheres, const double* sphere_charges,
                         const double* sphere_centres, const double* sphere_radii,
                         int* info);
void ddx_deallocate_model(void* ddx);

void ddx_get_error_message(const void* ddx, char* message, int maxlen);

// Generated block, see scripts/generate_cinterface.py
int ddx_get_enable_fmm(const void* ddx);
int ddx_get_enable_force(const void* ddx);
int ddx_get_jacobi_n_diis(const void* ddx);
int ddx_get_lmax(const void* ddx);
int ddx_get_incore(const void* ddx);
int ddx_get_maxiter(const void* ddx);
int ddx_get_model(const void* ddx);
int ddx_get_n_lebedev(const void* ddx);
int ddx_get_n_spheres(const void* ddx);
int ddx_get_n_proc(const void* ddx);
int ddx_get_fmm_local_lmax(const void* ddx);
int ddx_get_fmm_multipole_lmax(const void* ddx);
double ddx_get_solvent_epsilon(const void* ddx);
double ddx_get_eta(const void* ddx);
double ddx_get_solvent_kappa(const void* ddx);
void ddx_get_sphere_charges(const void* ddx, int nsph, double* c_charge);
void ddx_get_sphere_centres(const void* ddx, int nsph, double* c_csph);
void ddx_get_sphere_radii(const void* ddx, int nsph, double* c_rsph);
int ddx_get_n_basis(const void* ddx);
int ddx_get_n_cav(const void* ddx);
void ddx_get_cavity(const void* ddx, int ncav, double* c_ccav);
// end generated block

//
// State object
//
void* ddx_allocate_state(const void* ddx);
void ddx_deallocate_state(void* state);

void ddx_get_x(const void* state, int nbasis, int nsph, double* x);
int ddx_get_x_niter(const void* state);
void ddx_get_s(const void* state, int nbasis, int nsph, double* s);
int ddx_get_s_niter(const void* state);
// Hext one will change syntax later:
void ddx_get_xi(const void* state, const void* ddx, int ncav, double* xi);

// Cosmo
void ddx_cosmo_fill_guess(const void* ddx, void* state);
void ddx_cosmo_solve(const void* ddx, void* state, int ncav, const double* phi,
                     double tol);
void ddx_cosmo_adjoint(const void* ddx, void* state, int nbasis, int nsph,
                       const double* psi, double tol);
void ddx_cosmo_forces(const void* ddx, void* state, int nbasis, int nsph, int ncav,
                      const double* phi, const double* gradphi, const double* psi,
                      double* forces);

// PCM
void ddx_pcm_fill_guess(const void* ddx, void* state);
void ddx_pcm_solve(const void* ddx, void* state, int ncav, const double* phi, double tol);
void ddx_pcm_adjoint(const void* ddx, void* state, int nbasis, int nsph,
                     const double* psi, double tol);
void ddx_pcm_forces(const void* ddx, void* state, int nbasis, int nsph, int ncav,
                    const double* phi, const double* gradphi, const double* psi,
                    double* forces);

// LPB
// TODO

#ifdef __cplusplus
}
#endif

#endif /* DDX_H */
