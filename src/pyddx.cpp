#include "ddx.h"
#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void export_pyddx_data(py::module& m);     // defined in pyddx_data.cpp
void export_pyddx_classes(py::module& m);  // defined in pyddx_classes.cpp

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

  export_pyddx_classes(m);
  export_pyddx_data(m);
}
