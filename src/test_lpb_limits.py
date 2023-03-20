import pyddx
import numpy as np

tobohr = 1 / 0.52917721092

charges = np.array([0.08130, -0.20542, 0.06206, 0.06206])
rvdw = tobohr * np.array([4.00253, 3.63772, 2.99956, 2.99956])
centres = tobohr * np.array([
    [0.00000,  0.00000,  -1.16691],  # noqa: E201
    [0.00000,  0.00000,   1.42203],  # noqa: E201
    [0.00000,  1.76689,  -2.18736],  # noqa: E201
    [0.00000, -1.76689,  -2.18736],  # noqa: E201
]).T


def rundd(model, solvent_epsilon=78.3553, solvent_kappa=0.104):
    if model != "lpb":
        solvent_kappa = 0.0
    model = pyddx.Model(model, centres, rvdw, solvent_epsilon=solvent_epsilon,
                        solvent_kappa=solvent_kappa, lmax=10, n_lebedev=590)

    solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
    solute_field = model.multipole_electrostatics(solute_multipoles)
    solute_psi = model.multipole_psi(solute_multipoles)

    tol = 1e-12
    state = pyddx.State(model, solute_psi, solute_field["phi"], solute_field["e"])
    state.fill_guess(tol)
    state.solve(tol)
    return state.energy()


def test_lpb_limits():
    # TODO Also test solution vectors ?

    e_pcm = rundd("pcm")
    e_lpb0 = rundd("lpb", solvent_kappa=1e-6)
    e_lpbinf = rundd("lpb", solvent_kappa=18)
    e_cosmo = rundd("cosmo")

    assert abs(e_pcm - e_lpb0) < 1e-5
    assert abs(e_cosmo - e_lpbinf) < 1e-5
