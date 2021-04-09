#!/usr/bin/env python3

REMAP = {
    "force": "enable_force",
    "eps": "solvent_epsilon",
    "kappa": "solvent_kappa",
    "pm": "fmm_multipole_lmax",
    "pl": "fmm_local_lmax",
    "fmm": "enable_fmm",
    "charge": "sphere_charges",
    "csph": "sphere_centres",
    "rsph": "sphere_radii",
    "ngrid": "n_lebedev",
    "ncav": "n_cav",
    "ccav": "cavity",
    "nsph": "n_spheres",
    "nbasis": "n_basis",
    "jacobi_ndiis": "jacobi_n_diis",
    "nproc": "n_proc",
    "matvecmem": "incore",
}

PARAMS_INT = ["model", "force", "lmax", "ngrid", "maxiter", "jacobi_ndiis",
              "fmm", "pm", "pl", "nsph", "nproc", "matvecmem"]
PARAMS_FLOAT = ["eta", "eps", "kappa"]
PARAMS_ARRAY = [("charge", "nsph", "nsph"),
                ("csph",   "nsph", "3, nsph"),
                ("rsph",   "nsph", "nsph")]
DATA_INT = ["nbasis", "ncav"]
DATA_ARRAY = [("ccav", "ncav", "3, ncav")]


def generate_type(name, type, container="params"):
    if type == "int":
        typestr = "integer(c_int)"
        ctypestr = "int"
    elif type == "float":
        typestr = "real(c_double)"
        ctypestr = "double"
    else:
        raise ValueError("error")
    outname = REMAP.get(name, name)

    fortran = f"""
function ddx_get_{outname}(c_ddx) result(c_{name}) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    type(ddx_setup), pointer :: ddx
    {typestr} :: c_{name}
    call c_f_pointer(c_ddx, ddx)
    c_{name} = ddx % {container} % {name}
end function
    """.strip() + "\n\n"

    header = f"{ctypestr} ddx_get_{outname}(const void* ddx);\n"

    return dict(header=header, fortran=fortran)


def generate_float_array(name, sizes, sizes_full, container="params"):
    outname = REMAP.get(name, name)
    n_dim = sizes_full.count(",") + 1
    mask = ", ".join(n_dim * ":")
    sizeargs = ", ".join(["int " + s.strip() for s in sizes.split(",")
                          if s.strip() != ""])

    fortran = f"""
subroutine ddx_get_{outname}(c_ddx, {sizes}, c_{name}) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: {sizes}
    real(c_double), intent(out) :: c_{name}({sizes_full})
    type(ddx_setup), pointer :: ddx
    call c_f_pointer(c_ddx, ddx)
    c_{name}({mask}) = ddx % {container} % {name}({mask})
end subroutine
    """.strip() + "\n\n"

    header = f"""
    void ddx_get_{outname}(const void* ddx, {sizeargs}, double* c_{name});
    """.strip() + "\n"

    return dict(header=header, fortran=fortran)


fortran_iface = "! Generated block, see scripts/generate_cinterface.py\n"
c_header = "// Generated block, see scripts/generate_cinterface.py\n"
for name in sorted(PARAMS_INT):
    res = generate_type(name, "int", "params")
    fortran_iface += res["fortran"]
    c_header += res["header"]
for name in sorted(PARAMS_FLOAT):
    res = generate_type(name, "float", "params")
    fortran_iface += res["fortran"]
    c_header += res["header"]
for data in sorted(PARAMS_ARRAY):
    res = generate_float_array(*data, "params")
    fortran_iface += res["fortran"]
    c_header += res["header"]

for name in sorted(DATA_INT):
    res = generate_type(name, "int", "constants")
    fortran_iface += res["fortran"]
    c_header += res["header"]
for data in sorted(DATA_ARRAY):
    res = generate_float_array(*data, "constants")
    fortran_iface += res["fortran"]
    c_header += res["header"]

fortran_iface += "! end generated block\n"
c_header += "// end generated block\n"

print(fortran_iface)
print()
print("---------------------------------")
print()
print(c_header)
