#ifndef DDX_H
#define DDX_H

/** C header for interfacing with ddx */

#ifdef __cplusplus
extern "C" {
#endif

/** \name Generic quantities */
///@{
/** Fill the string banner with an informative header about ddx. At most maxlen characters
 * will be copied, so use this to specify the number of characters the char array can
 * hold. */
void ddx_get_banner(char* banner, int maxlen);

/** Fill the grids array with the number of grid points of the supported lebedev grids. At
 * most maxlen elements will be touched, so pass the size of the allocated integer array
 * as maxlen. Returns the number of entries in grids actually altered. E.g. if the
 * returned value is 5 the valid entries in the grid array after the call are grid[0],
 * ..., grid[4]. */
int ddx_supported_lebedev_grids(int maxlen, int* grids);

/** With reference to a atomic sphere `sphere` of radius `r` centred at `a` compute:
 * ```
 *       4π     |x - a|^l
 *     ------  ----------- Y_l^m(|x - a|)
 *     2l + 1       r^l
 * ```
 * where lmax should be identical to the value stored inside ddx or less.
 */
void ddx_scaled_ylm(const void* ddx, int lmax, const double* x, int sphere, double* ylm);
///@}

/** \name Allocate and manage the ddx model object */
///@{

/** Allocate a ddx model object.
 *
 * \param model            Integer describing the model to use.
 * \param enable_force     If 1 enable force calculations
 * \param solvent_epsilon  Relative dielectric permittivity
 * \param solvent_kappa    Debye-Hückel parameter (inverse screening length)
 * \param eta              Regularization parameter for the regularised characteristic
 *                         function chi_eta. Range [0,1]
 * \param shift            Shift for the regularized characteristic function chi_eta
 * \param lmax             Maximal degree of modelling spherical harmonics
 * \param n_lebedev        Number of Lebedev grid points to use
 * \param incore           If 1 store more large objects in memory
 * \param maxiter          Maximal number of iterations
 * \param jacobi_n_diis    Number of iterates stored in the DIIS space for acceleration
 * \param enable_fmm       If 1 enable the fast-multipole method for computations
 * \param fmm_multipole_lmax  Maximal degree of multipole spherical harmonics, ignored
 *                            in case `!enable_fmm`. Value `-1` means no far-field FFM
 *                            interactions are computed.
 * \param fmm_local_lmax   Maximal degree of local spherical harmonics, ignored in case
 *                         `use_fmm=false`. Value `-1` means no local FFM interactions
 *                         are computed.
 * \param n_proc           Number of processors to use
 * \param n_spheres        Number of cavity spheres
 * \param sphere_charges   The charges of the cavity spheres
 * \param sphere_centres   The centres of the cavity spheres as a *column-major*
 *                         (3, n_spheres) array.
 * \param sphere_radii     The radii of the spheres.
 * \param length_logfile   Length of the logfile name (0 if no log). Note that logs can
 *                         be quite verbose for larger runs and should only be used for
 *                         debugging.
 * \param logfile          The logfile name (empty if no log)
 */
void* ddx_allocate_model(int model, int enable_force, double solvent_epsilon,
                         double solvent_kappa, double eta, double shift, int lmax,
                         int n_lebedev, int incore, int maxiter, int jacobi_n_diis,
                         int enable_fmm, int fmm_multipole_lmax, int fmm_local_lmax,
                         int n_proc, int n_spheres, const double* sphere_charges,
                         const double* sphere_centres, const double* sphere_radii,
                         int length_logfile, const char* logfile);

/** Deallocate the model object */
void ddx_deallocate_model(void* ddx);

/** Get the error flag stored inside the model (0 is no error). It is the responsibility
 * of the user of DDX to regularly check for a non-zero error flag and take appropriate
 * actions if an error is present. */
int ddx_get_error_flag(const void* ddx);

/** Store the error message corresponding to the most recent error inside the passed char
 * array. at most maxlen characters are copied (so set this to the number of allocated
 * characters in the passed array to avoid a buffer overflow. */
void ddx_get_error_message(const void* ddx, char* message, int maxlen);
///@}

/** \name Getters for properties of the state object */
///@{

/** Copy the logfile used by the model to the passed logfile array. At most maxlen
 * characters are copied. */
void ddx_get_logfile(const void* ddx, char* logfile, int maxlen);

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
double ddx_get_shift(const void* ddx);
void ddx_get_sphere_charges(const void* ddx, int nsph, double* c_charge);
void ddx_get_sphere_centres(const void* ddx, int nsph, double* c_csph);
void ddx_get_sphere_radii(const void* ddx, int nsph, double* c_rsph);
int ddx_get_n_basis(const void* ddx);
int ddx_get_n_cav(const void* ddx);
void ddx_get_cavity(const void* ddx, int ncav, double* c_ccav);
// end generated block
///@}

/** \name Allocate and manage the state object */
///@{
/** Allocate an empty state from a model */
void* ddx_allocate_state(const void* ddx);

/** Deallocate a state object */
void ddx_deallocate_state(void* state);

/** Get the error flag stored inside the state (0 is no error). It is the responsibility
 * of the user of DDX to regularly check for a non-zero error flag and take appropriate
 * actions if an error is present. */
int ddx_get_state_error_flag(const void* state);

/** Store the error message corresponding to the most recent error inside the passed char
 * array. At most maxlen characters are copied. So best set maxlen to the number of
 * allocated characters in the passed array to avoid a buffer overflow. */
void ddx_get_state_error_message(const void* state, char* message, int maxlen);

/** Get the solution of the forward problem stored inside the state
 *  as a (nbasis, nsph) array in column-major ordering.
 *  \param state   DDX state
 *  \param nbasis  Number of basis functions used by DDX
 *  \param nsph    Number of cavity spheres
 *  \note This function just returns the data stored in the state,
 *  which only becomes valid once the appropriate ddx_cosmo_solve or
 *  ddx_pcm_solve have been called for the state.
 */
void ddx_get_x(const void* state, int nbasis, int nsph, double* x);

/** Get the number of iterations which where needed to obtain the solution
 *  returned by ddx_get_x */
int ddx_get_x_niter(const void* state);

/** Get the solution of the adjoint problem stored inside the state
 *  as a (nbasis, nsph) array in column-major ordering.
 *  \param state   DDX state
 *  \param nbasis  Number of basis functions used by DDX
 *  \param nsph    Number of cavity spheres
 *  \note This function just returns the data stored in the state,
 *  which only becomes valid once the appropriate ddx_cosmo_solve_adjoint or
 *  ddx_pcm_solve_adjoint have been called for the state.
 */
void ddx_get_s(const void* state, int nbasis, int nsph, double* s);

/** Get the number of iterations which where needed to obtain the adjoint solution
 *  returned by ddx_get_s */
int ddx_get_s_niter(const void* state);

/** Return the xi, the adjoint solution s projected onto the cavity points.
 *  \param state   DDX state
 *  \param ddx     DDX model
 *  \param ncav    Number of cavity points
 *  \param xi      Array to be filled with ncav
 */
void ddx_get_xi(const void* state, const void* ddx, int ncav, double* xi);
///@}

/** \name Problem setup and solution routines */
///@{

/** Setup a COSMO problem in the passed state.
 *  \param ddx     DDX model
 *  \param state   DDX state
 *  \param ncav    Number of cavity points
 *  \param nbasis  Number of basis functions used by DDX
 *  \param nsph    Number of cavity spheres
 *  \param phi_cav Phi array (ncav, )-shaped array
 *  \param psi     Psi array (nbasis, nsph)-shaped array (in column-major ordering)
 */
void ddx_cosmo_setup(const void* ddx, void* state, int ncav, int nbasis, int nsph,
                     const double* phi_cav, const double* psi);

/** In-place adjust the guess inside the state, getting ready to solve a COSMO problem.
 *  Avoid calling this step if you want to use the currently stored solution as an initial
 *  guess */
void ddx_cosmo_guess(const void* ddx, void* state);

/** In-place adjust the guess inside the state, getting ready to solve an adjoint COSMO
 *  problem. Avoid calling this step if you want to use the currently stored solution as
 *  an initial guess */
void ddx_cosmo_guess_adjoint(const void* ddx, void* state);

/** Solve the COSMO problem.
 *  \param tol Tolerance up to which the problem is solved */
void ddx_cosmo_solve(const void* ddx, void* state, double tol);

/** Solve the adjoint COSMO problem.
 *  \param tol Tolerance up to which the problem is solved */
void ddx_cosmo_solve_adjoint(const void* ddx, void* state, double tol);

/** Compute the COSMO force terms.
 *  \param ddx     DDX model
 *  \param state   DDX state
 *  \param nsph    Number of cavity spheres
 *  \param forces  Output force array (3, nsph) in column-major order
 */
void ddx_cosmo_solvation_force_terms(const void* ddx, void* state, int nsph,
                                     double* forces);

/** Setup a PCM problem in the passed state.
 *  \param ddx     DDX model
 *  \param state   DDX state
 *  \param ncav    Number of cavity points
 *  \param nbasis  Number of basis functions used by DDX
 *  \param nsph    Number of cavity spheres
 *  \param phi_cav Phi array (ncav, )-shaped array
 *  \param psi     Psi array (nbasis, nsph)-shaped array (in column-major ordering)
 */
void ddx_pcm_setup(const void* ddx, void* state, int ncav, int nbasis, int nsph,
                   const double* phi_cav, const double* psi);

/** In-place adjust the guess inside the state, getting ready to solve a PCM problem.
 *  Avoid calling this step if you want to use the currently stored solution as an initial
 *  guess */
void ddx_pcm_guess(const void* ddx, void* state);

/** In-place adjust the guess inside the state, getting ready to solve an adjoint PCM
 *  problem. Avoid calling this step if you want to use the currently stored solution as
 *  an initial guess */
void ddx_pcm_guess_adjoint(const void* ddx, void* state);

/** Solve the forward PCM problem.
 *  \param tol Tolerance up to which the problem is solved */
void ddx_pcm_solve(const void* ddx, void* state, double tol);

/** Solve the adjoint PCM problem.
 *  \param tol Tolerance up to which the problem is solved */
void ddx_pcm_solve_adjoint(const void* ddx, void* state, double tol);

/** Compute the PCM force terms.
 *  \param ddx     DDX model
 *  \param state   DDX state
 *  \param nsph    Number of cavity spheres
 *  \param forces  Output force array (3, nsph) in column-major order
 */
void ddx_pcm_solvation_force_terms(const void* ddx, void* state, int nsph,
                                   double* forces);

// TODO LPB

///@}

/** \name Multipolar solutes */
///@{

/** Build potential generated by a multipolar charge distribution.
 *  \param ddx      DDX model
 *  \param nsph     Number of cavity spheres
 *  \param ncav     Number of cavity points
 *  \param nmultipoles  Total number of multipoles. If multipoles up to mmax are used,
 *                      this is (mmax+1)^2. E.g. for charges, dipoles and quadrupoles
 *                      (mmax=2), this is 9.
 *  \param multipoles   Multipoles (nmultipoles, nsph)-shaped array in column-major order.
 *  \param phi_cav  Potential (phi) generated by multipoles: (ncav, )-shaped array
 */
void ddx_multipole_electrostatics_0(const void* ddx, int nsph, int ncav, int nmultipoles,
                                    const double* multipoles, double* phi_cav);

/** Build potential and electric field generated by a multipolar charge distribution.
 *  \param ddx      DDX model
 *  \param nsph     Number of cavity spheres
 *  \param ncav     Number of cavity points
 *  \param nmultipoles  Total number of multipoles. If multipoles up to mmax are used,
 *                      this is (mmax+1)^2. E.g. for charges, dipoles and quadrupoles
 *                      (mmax=2), this is 9.
 *  \param multipoles   Multipoles (nmultipoles, nsph)-shaped array in column-major order.
 *  \param phi_cav  Potential (phi) generated by multipoles: (ncav, )-shaped array
 *  \param e_cav    Electric field generated by multipoles: (3, ncav)-shaped array
 *                  (column major)
 */
void ddx_multipole_electrostatics_1(const void* ddx, int nsph, int ncav, int nmultipoles,
                                    const double* multipoles, double* phi_cav,
                                    double* e_cav);

/** Build potential, electric field and field gradient generated by a multipolar charge
 *  distribution.
 *  \param ddx      DDX model
 *  \param nsph     Number of cavity spheres
 *  \param ncav     Number of cavity points
 *  \param nmultipoles  Total number of multipoles. If multipoles up to mmax are used,
 *                      this is (mmax+1)^2. E.g. for charges, dipoles and quadrupoles
 *                      (mmax=2), this is 9.
 *  \param multipoles   Multipoles (nmultipoles, nsph)-shaped array in column-major order.
 *  \param phi_cav  Potential (phi) generated by multipoles: (ncav, )-shaped array
 *  \param e_cav    Electric field generated by multipoles: (3, ncav)-shaped array
 *                  (column major)
 *  \param g_cav    Electric field gradient gen. by multipoles: (3, 3, ncav)-shaped array
 *                  (column major)
 */
void ddx_multipole_electrostatics_2(const void* ddx, int nsph, int ncav, int nmultipoles,
                                    const double* multipoles, double* phi_cav,
                                    double* e_cav, double* g_cav);

/** Build the Psi generated by a multipolar charge distribution
 *  \param ddx     DDX model
 *  \param nbasis  Number of basis functions used by DDX
 *  \param nsph    Number of cavity spheres
 *  \param nmultipoles  Total number of multipoles. If multipoles up to mmax are used,
 *                      this is (mmax+1)^2. E.g. for charges, dipoles and quadrupoles
 *                      (mmax=2), this is 9.
 *  \param multipoles   Multipoles (nmultipoles, nsph)-shaped array in column-major order.
 *  \param psi     Psi array (nbasis, nsph)-shaped array (column major)
 */
void ddx_multipole_psi(const void* ddx, int nbasis, int nsph, int nmultipoles,
                       const double* multipoles, double* psi);

/** Compute the force terms generated by a multipolar charge distribution
 *  \param ddx     DDX model
 *  \param state   DDX state
 *  \param nsph    Number of cavity spheres
 *  \param ncav    Number of cavity points
 *  \param nmultipoles  Total number of multipoles. If multipoles up to mmax are used,
 *                      this is (mmax+1)^2. E.g. for charges, dipoles and quadrupoles
 *                      (mmax=2), this is 9.
 *  \param multipoles   Multipoles (nmultipoles, nsph)-shaped array in column-major order.
 *  \param e_cav   Electric field generated by multipoles: (3, ncav)-shaped array
 *                 (column major)
 *  \param forces  Output force array (3, nsph) in column-major order
 */
void ddx_multipole_forces(const void* ddx, void* state, int nsph, int ncav,
                          int nmultipoles, const double* multipoles, const double* e_cav,
                          double* forces);

///@}

#ifdef __cplusplus
}
#endif

#endif /* DDX_H */
