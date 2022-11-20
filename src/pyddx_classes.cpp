#include "ddx.h"
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <string>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;
using array_f_t = py::array_t<double, py::array::f_style | py::array::forcecast>;

class Model {
 public:
  Model(std::string model, array_f_t sphere_charges, array_f_t sphere_centres,
        array_f_t sphere_radii, double solvent_epsilon, double solvent_kappa, double eta,
        double shift, int lmax, int n_lebedev, bool incore, int maxiter,
        int jacobi_n_diis, bool enable_fmm, int fmm_multipole_lmax, int fmm_local_lmax,
        int n_proc, std::string logfile)
        : m_holder(nullptr), m_model(model) {
    int model_id = 0;
    if (model == "cosmo") {
      model_id = 1;
    } else if (model == "pcm") {
      model_id = 2;
    } else {
      throw py::value_error("Invalid model string: " + model);
    }
    if (model != "lpb" and solvent_kappa > 0.0) {
      throw py::value_error("Non-zero solvent_kappa only allowed for LPB");
    }

    // Check size of vdW and atomic data
    const size_t n_spheres = sphere_charges.size();
    if (sphere_charges.ndim() != 1) {
      throw py::value_error("Parameter sphere_charges is not a 1D array.");
    }
    if (sphere_radii.ndim() != 1) {
      throw py::value_error("Parameter sphere_radii is not a 1D array.");
    }
    if (n_spheres != sphere_radii.size()) {
      throw py::value_error("Length of 'sphere_charges' and 'sphere_radii' don't agree.");
    }
    if (sphere_centres.ndim() != 2) {
      throw py::value_error("sphere_centres is not a 2D array.");
    }
    if (3 != sphere_centres.shape(0) || n_spheres != sphere_centres.shape(1)) {
      throw py::value_error("sphere_centres should be a 3 x n_spheres array.");
    }
    if (maxiter < 0) {
      throw py::value_error("Maxiter should be non-zero.");
    }
    if (jacobi_n_diis < 0) {
      throw py::value_error("jacobi_n_diis should be positive.");
    }
    if (n_proc < 0) {
      throw py::value_error("n_proc should be positive.");
    }

    // Get supported Lebedev grids
    std::vector<int> supported_grids(100);
    const int n_supp_grids =
          ddx_supported_lebedev_grids(supported_grids.size(), supported_grids.data());
    supported_grids.resize(static_cast<size_t>(n_supp_grids));
    if (std::find(supported_grids.begin(), supported_grids.end(), n_lebedev) ==
        supported_grids.end()) {
      std::string msg = "Lebedev grid size '" + std::to_string(n_lebedev) +
                        "' not supported. Supported grid sizes are: ";
      for (size_t i = 0; i < supported_grids.size(); ++i) {
        const std::string separator = i < supported_grids.size() ? ", " : "";
        msg += std::to_string(supported_grids[i]) + separator;
      }
      throw py::value_error(msg);
    }

    if (shift == -100.0) {
      if (model == "cosmo") {
        shift = -1.0;
      } else if (model == "pcm" || model == "lpb") {
        shift = 0.0;
      } else {
        throw std::runtime_error("No default shift known for model: " + model);
      }
    }

    if (eta < 0 || eta > 1) {
      throw py::value_error(
            "Regularisation parameter 'eta' needs to be between 0 and 1.");
    }
    if (shift < -1 || shift > 1) {
      throw py::value_error(
            "Regularisation parameter 'shift' needs to be between -1 and 1. "
            "Alternatively use the value -100  to enable an automatic selection.");
    }
    if (solvent_epsilon <= 0) {
      throw py::value_error(
            "Dielectric permittivity 'solvent_epsilon' needs to be positive.");
    }
    if (lmax < 0) {
      throw py::value_error("Maximal spherical harmonics degree 'lmax' needs to be >= 0");
    }
    if (fmm_multipole_lmax < -1) {
      throw py::value_error(
            "Maximal spherical harmonics degree 'fmm_multipole_lmax' needs to >= -1 with "
            "-1 disabling far-field FMM contributions.");
    }
    if (fmm_local_lmax < -1) {
      throw py::value_error(
            "Maximal spherical harmonics degree 'fmm_local_lmax' needs to >= -1 with -1 "
            "disabling local FMM contributions.");
    }

    const int enable_force = 1;  // Always support force calculations.
    const int intfmm       = enable_fmm ? 1 : 0;
    const int intincore    = incore ? 1 : 0;

    m_holder = ddx_allocate_model(
          model_id, enable_force, solvent_epsilon, solvent_kappa, eta, shift, lmax,
          n_lebedev, intincore, maxiter, jacobi_n_diis, intfmm, fmm_multipole_lmax,
          fmm_local_lmax, n_proc, n_spheres, sphere_charges.data(), sphere_centres.data(),
          sphere_radii.data(), logfile.size(), logfile.c_str());
    throw_if_error();
  }

  ~Model() {
    if (m_holder) ddx_deallocate_model(m_holder);
    m_holder = nullptr;
  }
  Model(const Model&)            = delete;
  Model& operator=(const Model&) = delete;

  // Check for errors
  void throw_if_error() const {
    if (ddx_get_error_flag(holder()) != 0) {
      char message[256];
      ddx_get_error_message(holder(), message, 256);
      throw std::runtime_error(std::string(message));
    }
  }

  //
  // Accessors
  //

  // Return the holder pointer ... internal function. Use only if you know
  // what you are doing.
  void* holder() const { return m_holder; }

  bool has_fmm_enabled() const { return 1 == ddx_get_enable_fmm(m_holder); }
  bool has_force_enabled() const { return 1 == ddx_get_enable_force(m_holder); }
  int jacobi_n_diis() const { return ddx_get_jacobi_n_diis(m_holder); }
  int lmax() const { return ddx_get_lmax(m_holder); }
  int maxiter() const { return ddx_get_maxiter(m_holder); }
  bool incore() const { return 1 == ddx_get_incore(m_holder); }
  std::string model() const { return m_model; }
  int n_lebedev() const { return ddx_get_n_lebedev(m_holder); }
  int n_spheres() const { return ddx_get_n_spheres(m_holder); }
  int n_proc() const { return ddx_get_n_proc(m_holder); }
  std::string logfile() const {
    char filename[256];
    ddx_get_logfile(m_holder, filename, 256);
    return std::string(filename);
  }
  int fmm_local_lmax() const { return ddx_get_fmm_local_lmax(m_holder); }
  int fmm_multipole_lmax() const { return ddx_get_fmm_multipole_lmax(m_holder); }
  double solvent_epsilon() const { return ddx_get_solvent_epsilon(m_holder); }
  double eta() const { return ddx_get_eta(m_holder); }
  double shift() const { return ddx_get_shift(m_holder); }
  double solvent_kappa() const { return ddx_get_solvent_kappa(m_holder); }

  array_f_t sphere_charges() const {
    array_f_t result({n_spheres()});
    ddx_get_sphere_charges(m_holder, n_spheres(), result.mutable_data());
    return result;
  }
  array_f_t sphere_centres() const {
    array_f_t result({3, n_spheres()});
    ddx_get_sphere_centres(m_holder, n_spheres(), result.mutable_data());
    return result;
  }
  array_f_t sphere_radii() const {
    array_f_t result({n_spheres()});
    ddx_get_sphere_radii(m_holder, n_spheres(), result.mutable_data());
    return result;
  }

  // Get all input parameters as a dict
  py::dict input_parameters() const {
    return py::dict(
          "model"_a = model(), "sphere_charges"_a = sphere_charges(),
          "sphere_centres"_a = sphere_centres(), "sphere_radii"_a = sphere_radii(),
          "solvent_epsilon"_a = solvent_epsilon(), "solvent_kappa"_a = solvent_kappa(),
          "eta"_a = eta(), "shift"_a = shift(), "lmax"_a = lmax(),
          "n_lebedev"_a = n_lebedev(), "maxiter"_a = maxiter(), "incore"_a = incore(),
          "jacobi_n_diis"_a = jacobi_n_diis(), "enable_fmm"_a = has_fmm_enabled(),
          "fmm_multipole_lmax"_a = fmm_multipole_lmax(),
          "fmm_local_lmax"_a = fmm_local_lmax(), "n_proc"_a = n_proc(),
          "logfile"_a = logfile());
  }

  // Accessors to derived parameters
  int n_basis() const { return ddx_get_n_basis(m_holder); }
  int n_cav() const { return ddx_get_n_cav(m_holder); }

  array_f_t cavity() const {
    array_f_t coords({3, n_cav()});
    ddx_get_cavity(m_holder, n_cav(), coords.mutable_data());
    return coords;
  }

 private:
  void* m_holder;
  std::string m_model;
};

class State {
 public:
  State(std::shared_ptr<Model> model)
        : m_holder(ddx_allocate_state(model->holder())),
          m_model(model),
          m_solved(false),
          m_solved_adjoint(false) {}
  State(std::shared_ptr<Model> model, array_f_t phi, array_f_t psi) : State(model) {
    update_problem(phi, psi);
    fill_guess();
    fill_guess_adjoint();
  }
  ~State() {
    if (m_holder) ddx_deallocate_state(m_holder);
    m_holder = nullptr;
  }
  State(const State&)            = delete;
  State& operator=(const State&) = delete;

  void throw_if_error() const {
    if (ddx_get_state_error_flag(holder()) != 0) {
      char message[256];
      ddx_get_state_error_message(holder(), message, 256);
      throw std::runtime_error(std::string(message));
    }
    model()->throw_if_error();
  }

  //
  // Accessors
  //

  // Return the holder pointer ... internal function. Use only if you know
  // what you are doing.
  void* holder() const { return m_holder; }
  std::shared_ptr<Model> model() const { return m_model; }

  bool is_solved() const { return m_solved; }
  bool is_solved_adjoint() const { return m_solved_adjoint; }

  array_f_t x() const {
    check_solved();
    array_f_t x({m_model->n_basis(), m_model->n_spheres()});
    ddx_get_x(m_holder, m_model->n_basis(), m_model->n_spheres(), x.mutable_data());
    return x;
  }
  int x_n_iter() const {
    check_solved();
    return ddx_get_x_niter(m_holder);
  }
  array_f_t s() const {
    check_solved_adjoint();
    array_f_t s({m_model->n_basis(), m_model->n_spheres()});
    ddx_get_s(m_holder, m_model->n_basis(), m_model->n_spheres(), s.mutable_data());
    return s;
  }
  int s_n_iter() const {
    check_solved_adjoint();
    return ddx_get_s_niter(m_holder);
  }
  array_f_t xi() const {
    check_solved_adjoint();
    array_f_t xi({m_model->n_cav()});
    ddx_get_xi(m_holder, m_model->holder(), m_model->n_cav(), xi.mutable_data());
    return xi;
  }

  //
  // Solving COSMO / PCM
  //
  void update_problem(array_f_t phi, array_f_t psi) {
    if (psi.ndim() != 2 || psi.shape(0) != model()->n_basis() ||
        psi.shape(1) != model()->n_spheres()) {
      throw py::value_error("psi not of shape (n_basis, n_spheres) == (" +
                            std::to_string(model()->n_basis()) + ", " +
                            std::to_string(model()->n_spheres()) + ").");
    }
    if (phi.ndim() != 1 || phi.shape(0) != model()->n_cav()) {
      throw py::value_error("phi not of shape (n_cav, ) == (" +
                            std::to_string(model()->n_cav()) + ").");
    }

    if (model()->model() == "cosmo") {
      ddx_cosmo_setup(model()->holder(), holder(), model()->n_cav(), model()->n_basis(),
                      model()->n_spheres(), phi.data(), psi.data());
    } else if (model()->model() == "pcm") {
      ddx_pcm_setup(model()->holder(), holder(), model()->n_cav(), model()->n_basis(),
                    model()->n_spheres(), phi.data(), psi.data());
    } else {
      throw py::value_error("Model " + model()->model() + " not yet implemented.");
    }
    throw_if_error();
    m_solved         = false;
    m_solved_adjoint = false;
  }

  void fill_guess() {
    if (model()->model() == "cosmo") {
      ddx_cosmo_guess(model()->holder(), holder());
    } else if (model()->model() == "pcm") {
      ddx_pcm_guess(model()->holder(), holder());
    } else {
      throw py::value_error("Model " + model()->model() + " not yet implemented.");
    }
    throw_if_error();
    m_solved = false;
  }

  void fill_guess_adjoint() {
    if (model()->model() == "cosmo") {
      ddx_cosmo_guess_adjoint(model()->holder(), holder());
    } else if (model()->model() == "pcm") {
      ddx_pcm_guess_adjoint(model()->holder(), holder());
    } else {
      throw py::value_error("Model " + model()->model() + " not yet implemented.");
    }
    throw_if_error();
    m_solved_adjoint = false;
  }

  // Solve the forward COSMO / PCM System. The state is modified in-place.
  void solve(double tol) {
    if (model()->model() == "cosmo") {
      ddx_cosmo_solve(model()->holder(), holder(), tol);
    } else if (model()->model() == "pcm") {
      ddx_pcm_solve(model()->holder(), holder(), tol);
    } else {
      throw py::value_error("Model " + model()->model() + " not yet implemented.");
    }
    throw_if_error();
    m_solved = true;
  }

  // Solve the adjoint COSMO / PCM System. The state is modified in-place.
  void solve_adjoint(double tol) {
    if (model()->model() == "cosmo") {
      ddx_cosmo_solve_adjoint(model()->holder(), holder(), tol);
    } else if (model()->model() == "pcm") {
      ddx_pcm_solve_adjoint(model()->holder(), holder(), tol);
    } else {
      throw py::value_error("Model " + model()->model() + " not yet implemented.");
    }
    throw_if_error();
    m_solved_adjoint = true;
  }

  // Obtain the COSMO / PCM forces. The state is modified in-place
  array_f_t solvation_force_terms() {
    check_solved_adjoint();
    array_f_t forces({3, model()->n_spheres()});
    if (model()->model() == "cosmo") {
      ddx_cosmo_solvation_force_terms(model()->holder(), holder(), model()->n_spheres(),
                                      forces.mutable_data());
    } else if (model()->model() == "pcm") {
      ddx_pcm_solvation_force_terms(model()->holder(), holder(), model()->n_spheres(),
                                    forces.mutable_data());
    } else {
      throw py::value_error("Model " + model()->model() + " not yet implemented.");
    }
    throw_if_error();
    return forces;
  }

  //
  // Multipolar functions
  //
  array_f_t multipole_force_terms(array_f_t multipoles, array_f_t e) {
    if (e.ndim() != 2 || e.shape(0) != 3 || e.shape(1) != model()->n_cav()) {
      throw py::value_error("e not of shape (3, n_cav) == (3, " +
                            std::to_string(model()->n_cav()) + ").");
    }
    if (multipoles.ndim() != 2 || multipoles.shape(1) != model()->n_spheres()) {
      throw py::value_error(
            "'multipoles' should be an (nmultipoles, n_spheres) sized array");
    }
    int mmax = static_cast<int>(sqrt(multipoles.shape(0)) - 1.0);  // Multipole order
    if (mmax < 0 || (mmax + 1) * (mmax + 1) != multipoles.shape(0)) {
      throw py::value_error(
            "First axis of multipole array is not of form (mmax+1)^2 for an mmax >= 0");
    }

    array_f_t forces({3, model()->n_spheres()});
    ddx_multipole_forces(model()->holder(), holder(), model()->n_spheres(),
                         model()->n_cav(), multipoles.shape(0), multipoles.data(),
                         e.data(), forces.mutable_data());
    throw_if_error();
    return forces;
  }

 private:
  void check_solved() const {
    if (!m_solved)
      throw std::runtime_error(
            "State does currently not hold a solution to the forward problem. Call "
            "solve() first.");
  }
  void check_solved_adjoint() const {
    if (!m_solved_adjoint)
      throw std::runtime_error(
            "State does currently not hold a solution to the adjoint problem. Call "
            "solve_adjoint() first.");
  }

  void* m_holder;
  std::shared_ptr<Model> m_model;
  bool m_solved;          // Has forward problem been solved
  bool m_solved_adjoint;  // Has adjoint problem been solved
};

//
// General functions
//

py::array_t<double> scaled_ylm(std::shared_ptr<Model> model, array_f_t coord, int sphere,
                               py::array_t<double> out) {
  if (out.size() != model->n_basis()) {
    throw py::value_error("'out' should have the same size as `n_basis()`");
  }
  if (coord.size() != 3) {
    throw py::value_error("'coord' should have exactly three entries");
  }
  if (sphere >= model->n_spheres()) {
    throw py::value_error("'sphere' should be less than n_spheres()");
  }
  ddx_scaled_ylm(model->holder(), model->lmax(), coord.data(), sphere + 1,
                 out.mutable_data());
  model->throw_if_error();
  return out;
}
array_f_t scaled_ylm(std::shared_ptr<Model> model, array_f_t coord, int sphere) {
  return scaled_ylm(model, coord, sphere, array_f_t({model->n_basis()}));
}

//
// Multipolar electrostatics
//

py::dict multipole_electrostatics(std::shared_ptr<Model> model, array_f_t multipoles,
                                  int derivative_order) {
  if (multipoles.ndim() != 2 || multipoles.shape(1) != model->n_spheres()) {
    throw py::value_error(
          "'multipoles' should be an (nmultipoles, n_spheres) sized array");
  }
  int mmax = static_cast<int>(sqrt(multipoles.shape(0)) - 1.0);  // Multipole order
  if (mmax < 0 || (mmax + 1) * (mmax + 1) != multipoles.shape(0)) {
    throw py::value_error(
          "First axis of multipole array is not of form (mmax+1)^2 for an mmax >= 0");
  }

  array_f_t phi({model->n_cav()});
  if (derivative_order == 0) {
    ddx_multipole_electrostatics_0(model->holder(), model->n_spheres(), model->n_cav(),
                                   multipoles.shape(0), multipoles.data(),
                                   phi.mutable_data());
    model->throw_if_error();
    return py::dict("phi"_a = phi);
  } else if (derivative_order == 1) {
    array_f_t e({3, model->n_cav()});
    ddx_multipole_electrostatics_1(model->holder(), model->n_spheres(), model->n_cav(),
                                   multipoles.shape(0), multipoles.data(),
                                   phi.mutable_data(), e.mutable_data());
    model->throw_if_error();
    return py::dict("phi"_a = phi, "e"_a = e);
  } else if (derivative_order == 2) {
    array_f_t e({3, model->n_cav()});
    array_f_t g({3, 3, model->n_cav()});
    ddx_multipole_electrostatics_2(
          model->holder(), model->n_spheres(), model->n_cav(), multipoles.shape(0),
          multipoles.data(), phi.mutable_data(), e.mutable_data(), g.mutable_data());
    model->throw_if_error();
    return py::dict("phi"_a = phi, "e"_a = e, "g"_a = g);
  } else {
    throw py::value_error("'derivative_order' only implemented up to 2");
  }
}

array_f_t multipole_psi(std::shared_ptr<Model> model, array_f_t multipoles) {
  if (multipoles.ndim() != 2 || multipoles.shape(1) != model->n_spheres()) {
    throw py::value_error(
          "'multipoles' should be an (nmultipoles, n_spheres) sized array");
  }
  int mmax = static_cast<int>(sqrt(multipoles.shape(0)) - 1.0);  // Multipole order
  if (mmax < 0 || (mmax + 1) * (mmax + 1) != multipoles.shape(0)) {
    throw py::value_error(
          "First axis of multipole array is not of form (mmax+1)^2 for an mmax >= 0");
  }

  array_f_t psi({model->n_basis(), model->n_spheres()});
  ddx_multipole_psi(model->holder(), model->n_basis(), model->n_spheres(),
                    multipoles.shape(0), multipoles.data(), psi.mutable_data());
  model->throw_if_error();
  return psi;
}

//
// Python export
//

void export_pyddx_classes(py::module& m) {
  const char* init_docstring =
        "Setup solvation model in ddX.\n\n"
        "model:            'cosmo', 'pcm' or 'lpb'\n"
        "sphere_charges:   (n_spheres) array\n"
        "atomic_centers:   (n_spheres, 3) array\n"
        "sphere_radii:     (n_spheres) array\n"
        "solvent_epsilon:  Relative dielectric permittivity\n"
        "solvent_kappa:    Debye-Hückel parameter (inverse screening length)\n"
        "eta:              Regularization parameter for the regularised characteristic "
        "function chi_eta. Range [0,1]\n"
        "shift:            Shift for the regularized characteristic function chi_eta, "
        "Range [-1, 1]. The default (-100) selects a suitable value automatically.\n"
        "lmax:             Maximal degree of modelling spherical harmonics\n"
        "n_lebedev:        Number of Lebedev grid points to use\n"
        "incore:           Store more large objects in memory\n"
        "maxiter:          Maximal number of iterations\n"
        "jacobi_n_diis:    Number of iterates stored in the DIIS space for acceleration\n"
        "enable_fmm:       Use fast-multipole method (true) or not (false)\n"
        "fmm_multipole_lmax:  Maximal degree of multipole spherical harmonics, "
        "ignored "
        "in case `!enable_fmm`. Value `-1` means no far-field FFM interactions "
        "are computed. Using the same value as lmax is recommended. \n"
        "fmm_local_lmax:   Maximal degree of local spherical harmonics, ignored in "
        "case `use_fmm=false`. Value `-1` means no local FFM interactions are "
        "computed.\n"
        "n_proc:           Number of processors to use.\n"
        "logfile:          Logfile for debugging (very verbose output, use with care)\n";

  py::class_<Model, std::shared_ptr<Model>>(m, "Model",
                                            "Solvation model using ddX library.")
        .def(py::init<std::string, array_f_t, array_f_t, array_f_t, double, double,
                      double, double, int, int, int, int, int, bool, int, int, int,
                      std::string>(),
             init_docstring, "model"_a, "sphere_charges"_a, "sphere_centres"_a,
             "sphere_radii"_a, "solvent_epsilon"_a, "solvent_kappa"_a = 0.0,
             "eta"_a = 0.1, "shift"_a = -100, "lmax"_a = 7, "n_lebedev"_a = 302,
             "incore"_a = false, "maxiter"_a = 100, "jacobi_n_diis"_a = 20,
             "enable_fmm"_a = true, "fmm_multipole_lmax"_a = 7, "fmm_local_lmax"_a = 6,
             "n_proc"_a = 1, "logfile"_a = "")
        //
        .def_property_readonly("has_fmm_enabled", &Model::has_fmm_enabled)
        .def_property_readonly("has_force_enabled", &Model::has_force_enabled)
        .def_property_readonly("jacobi_n_diis", &Model::jacobi_n_diis)
        .def_property_readonly("lmax", &Model::lmax)
        .def_property_readonly("incore", &Model::incore)
        .def_property_readonly("maxiter", &Model::maxiter)
        .def_property_readonly("model", &Model::model)
        .def_property_readonly("n_lebedev", &Model::n_lebedev)
        .def_property_readonly("n_spheres", &Model::n_spheres)
        .def_property_readonly("n_proc", &Model::n_proc)
        .def_property_readonly("logfile", &Model::logfile)
        .def_property_readonly("fmm_local_lmax", &Model::fmm_local_lmax)
        .def_property_readonly("fmm_multipole_lmax", &Model::fmm_multipole_lmax)
        .def_property_readonly("solvent_epsilon", &Model::solvent_epsilon)
        .def_property_readonly("eta", &Model::eta)
        .def_property_readonly("shift", &Model::shift)
        .def_property_readonly("solvent_kappa", &Model::solvent_kappa)
        .def_property_readonly("sphere_charges", &Model::sphere_charges)
        .def_property_readonly("sphere_centres", &Model::sphere_centres)
        .def_property_readonly("sphere_radii", &Model::sphere_radii)
        .def_property_readonly("input_parameters", &Model::input_parameters)
        //
        .def_property_readonly("n_basis", &Model::n_basis)
        .def_property_readonly("n_cav", &Model::n_cav)
        .def_property_readonly("cavity", &Model::cavity)
        //
        .def(
              "scaled_ylm",
              [](std::shared_ptr<Model> model, array_f_t x, int sphere) {
                return scaled_ylm(model, x, sphere);
              },
              "x"_a, "sphere"_a,
              "With reference to a atomic sphere `sphere` of radius `r` centred at "
              "`a` "
              "compute 4π/(2l+1) * (|x-a|/r)^l * Y_l^m(|x-a|).")
        .def(
              "scaled_ylm",
              [](std::shared_ptr<Model> model, array_f_t x, int sphere,
                 py::array_t<double> out) { return scaled_ylm(model, x, sphere, out); },
              "coord"_a, "sphere"_a, "out"_a,
              "With reference to a atomic sphere `sphere` of radius `r` centred at "
              "`a` compute 4π/(2l+1) * (|x-a|/r)^l * Y_l^m(|x-a|).")
        //
        .def("multipole_electrostatics", &multipole_electrostatics,
             "Return the solute potential, electric field and field gradients for "
             "a solute represented by a distribution of multipoles on the cavity "
             "centres. `multipoles` is a (nmultipoles, n_spheres) array, where "
             "nmultipoles "
             "are the number of multipoles on each site (i.e. (mmax+x)^2 entries) where "
             "mmax is the maximal angular momentum of the multipoles. "
             "The order of potential derivatives is given by the "
             "'derivative_order' flag "
             "(0 for just the potential 'phi', 1 for 'phi' and field 'e', 2 for "
             "'phi', 'e' and field gradient 'g').",
             "multipoles"_a, "derivative_order"_a = 1)
        .def("multipole_psi", &multipole_psi,
             "Return the solute contribution to psi generated from a distribution "
             "of multipoles on the cavity centres.",
             "multipoles"_a)
        //
        ;

  py::class_<State, std::shared_ptr<State>>(
        m, "State", "Computational state and results of ddX models")
        .def(py::init<std::shared_ptr<Model>, array_f_t, array_f_t>(),
             "Construct a state from the model and a phi and psi to set up the problem.",
             "model"_a, "phi"_a, "psi"_a)
        .def_property_readonly("model", &State::model, "Model definition")
        .def_property_readonly("x", &State::x, "Solution of the forward problem.")
        .def_property_readonly(
              "x_n_iter", &State::x_n_iter,
              "Number of iterations required to solve the forward problem.")
        .def_property_readonly("s", &State::s, "Solution of the adjoint problem.")
        .def_property_readonly(
              "s_n_iter", &State::s_n_iter,
              "Number of iterations required to solve the adjoint problem.")
        .def_property_readonly("xi", &State::xi)
        .def_property_readonly("is_solved", &State::is_solved,
                               "Return whether the forward problem is solved.")
        .def_property_readonly("is_solved_adjoint", &State::is_solved_adjoint,
                               "Return whether the adjoint problem is solved.")
        //
        .def("update_problem", &State::update_problem,
             "Update the definition of of the forward and adjoint problem", "phi"_a,
             "psi"_a)
        .def("fill_guess", &State::fill_guess,
             "In-place construct an initial guess for the adjoint solver. Don't call "
             "this if you want to use the previous solution stored in this state as a "
             "guess")
        .def("fill_guess_adjoint", &State::fill_guess_adjoint,
             "In-place construct an initial guess for the forward solver. Don't call "
             "this if you want to use the previous solution stored in this state as a "
             "guess")
        .def("solve", &State::solve, "Solve the forward problem contained in the state.",
             "tol"_a = 1e-8)
        .def("solve_adjoint", &State::solve_adjoint,
             "Solve the adjoint problem contained in the state.", "tol"_a = 1e-8)
        .def("solvation_force_terms", &State::solvation_force_terms,
             "Compute and return the force terms of the solvation part of the solvation "
             "model.")
        .def("multipole_force_terms", &State::multipole_force_terms,
             "Obtain the solute force contributions from a solute represented using "
             "multipoles. ",
             "multipoles"_a, "e"_a);
}
