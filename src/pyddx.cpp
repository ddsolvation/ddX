#include "ddx.h"
#include <pybind11/eval.h>
#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void export_pyddx_data(py::module& m);     // defined in pyddx_data.cpp
void export_pyddx_classes(py::module& m);  // defined in pyddx_classes.cpp
void export_pyddx_kappa(py::module& m);    // defined below.

PYBIND11_MODULE(pyddx, m) {
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);

  m.def(
        "banner",
        []() {
          char message[2048];
          ddx_get_banner(message, 2048);
          return std::string(message);
        },
        "Return a nice banner describing ddx.");

  export_pyddx_kappa(m);
  export_pyddx_classes(m);
  export_pyddx_data(m);
}

void export_pyddx_kappa(py::module& m) {
  py::exec(R"(
def solvent_kappa(list_of_ions, temperature, solvent_epsilon):
    import numpy as np
    import scipy.constants
    from scipy.constants import epsilon_0  # vacuum permittivity (C/V/m)
    from scipy.constants import e          # electron charge (C)
    from scipy.constants import R          # Gas constant
    from scipy.constants import k          # Boltzmann constant (J/K)
    from scipy.constants import N_A        # Avogadro (molecules/mol)

    a0 = scipy.constants.value("atomic unit of length")

    # Compute ionic strength
    ionic_strength = 0.0
    for charge, concentration in list_of_ions:
        ionic_strength += concentration * charge**2

    ionic_strength_si = ionic_strength * 1000 * N_A  # Convert to molecules / m^3
    kappa_si = np.sqrt(e**2 * ionic_strength_si
                       / (epsilon_0 * solvent_epsilon * k * temperature))
    return kappa_si * a0  # Convert to Bohrs

solvent_kappa.__doc__ = (
    "Compute the Debye-Hückel parameter in inverse Bohr for a solvent " +
    "with dielectric constant solvent_epsilon and a list of dissolved " +
    "ions. The ion concentration is given in mol/l.\n\nExample:\n" +
    "kappa([(+1, 0.1), (-1, 0.1)], 78.54, 298.15)\ncomputes the parameter " +
    "for a 0.1 mol/l solution of sodium chloride in water at room temperature."
))",
           m.attr("__dict__"));
}
