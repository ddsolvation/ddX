import pyddx
import numpy as np

tobohr = 1 / 0.52917721092
toang = 0.52917721092
tokcal = 627.509469

charges = np.array([
    -0.04192, -0.04192, -0.04198, -0.04192, -0.04192, -0.04198,
    0.04193, 0.04193,  0.04197,  0.04193,  0.04193,  0.04197
])
rvdw = tobohr * np.array([
    4.00253, 4.00253, 4.00253, 4.00253, 4.00253, 4.00253,
    2.99956, 2.99956, 2.99956, 2.99956, 2.99956, 2.99956
])
centres = tobohr * np.array([
    [ 0.00000,  2.29035,  1.32281],  # noqa: E201
    [ 0.00000,  2.29035, -1.32281],  # noqa: E201
    [ 0.00000,  0.00000, -2.64562],  # noqa: E201
    [ 0.00000, -2.29035, -1.32281],  # noqa: E201
    [ 0.00000, -2.29035,  1.32281],  # noqa: E201
    [ 0.00000,  0.00000,  2.64562],  # noqa: E201
    [ 0.00103,  4.05914,  2.34326],  # noqa: E201
    [ 0.00103,  4.05914, -2.34326],  # noqa: E201
    [ 0.00000,  0.00000, -4.68652],  # noqa: E201
    [-0.00103, -4.05914, -2.34326],  # noqa: E201
    [-0.00103, -4.05914,  2.34326],  # noqa: E201
    [ 0.00000,  0.00000,  4.68652],  # noqa: E201
]).T

def ddsolve(model_name, charges, rvdw, centres, do_forces=False):
    solvent_kappa = 0.0
    if model_name == "lpb":
        solvent_kappa = 0.1
    model = pyddx.Model(model_name, centres, rvdw, solvent_epsilon=1e10, lmax=10, n_lebedev=302,
        shift=0.0, jacobi_n_diis=25, enable_fmm=False, solvent_kappa=solvent_kappa)
    solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
    solute_field = model.multipole_electrostatics(solute_multipoles)
    solute_psi = model.multipole_psi(solute_multipoles)

    if model_name == "lpb":
        state = pyddx.State(model, solute_psi, solute_field["phi"], solute_field["e"])
    else:
        state = pyddx.State(model, solute_psi, solute_field["phi"])
    state.fill_guess()
    state.solve()
    energy = state.energy()

    forces = np.zeros(centres.shape)
    if do_forces:
        state.fill_guess_adjoint()
        state.solve_adjoint()
        forces = state.solvation_force_terms(solute_field)
        forces += state.multipole_force_terms(solute_multipoles);

    return energy, forces


def numerical_forces(model_name, charges, rvdw, centres, step=1e-3):
    forces = np.zeros(centres.shape)
    _, natoms = centres.shape
    for i in range(natoms):
        print(i, natoms)
        for j in range(3):
            centres[j, i] -= step
            e_minus, _ = ddsolve(model_name, charges, rvdw, centres, False)
            centres[j, i] += 2.0*step
            e_plus, _ = ddsolve(model_name, charges, rvdw, centres, False)
            forces[j, i] = (e_plus - e_minus)/(2*step)
            centres[j, i] -= step
    return forces

def test_ddcosmo():
    _, force = ddsolve("cosmo", charges, rvdw, centres, True)
    num_force = numerical_forces("cosmo", charges, rvdw, centres)
    assert np.max(np.abs(force - num_force)) < 1e-8

def test_ddpcm():
    _, force = ddsolve("pcm", charges, rvdw, centres, True)
    num_force = numerical_forces("pcm", charges, rvdw, centres)
    assert np.max(np.abs(force - num_force)) < 1e-8

def test_ddlpb():
    _, force = ddsolve("lpb", charges, rvdw, centres, True)
    num_force = numerical_forces("lpb", charges, rvdw, centres)
    assert np.max(np.abs(force - num_force)) < 1e-8
