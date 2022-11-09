# Define common source files
set(SRC
    ddx_definitions.f90
    ddx_parameters.f90
    ddx_constants.f90
    ddx_workspace.f90
    ddx_harmonics.f90
    ddx_core.f90
    ddx_lpb_core.f90
    ddx_operators.f90
    ddx_gradients.f90
    ddx_solvers.f90
    ddx_cosmo.f90
    ddx_pcm.f90
    ddx_lpb.f90
    ddx.f90
    ddx_mkrhs.f90
    llgnew.f
    cbessel.f90
    ddx_cinterface.f90
)

# Define source files to be build into the python version of ddx
set(SRC_PYDDX ${SRC} pyddx.cpp pyddx_classes.cpp pyddx_data.cpp)

# Define source files for the ddx.so shared library
set(SRC_DDX ${SRC} matrix_debug.f90)
