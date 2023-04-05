import numpy as np

import pyddx

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
                        solvent_kappa=solvent_kappa, lmax=10, n_lebedev=590, shift=0.0)

    solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
    solute_field = model.multipole_electrostatics(solute_multipoles)
    solute_psi = model.multipole_psi(solute_multipoles)

    tol = 1e-12
    state = pyddx.State(model, solute_psi, solute_field["phi"], solute_field["e"])
    state.fill_guess(tol)
    state.fill_guess_adjoint(tol)
    state.solve(tol)
    state.solve_adjoint(tol)
    return state


def test_lpb_limits():
    s_pcm = rundd("pcm")
    s_lpb0 = rundd("lpb", solvent_kappa=1e-6)
    s_lpbinf = rundd("lpb", solvent_kappa=18)
    s_lpb_einf = rundd("lpb", solvent_epsilon=1e8)
    s_cosmo = rundd("cosmo")

    assert abs(s_lpb0.energy() - s_pcm.energy()) < 1e-5
    assert np.max(np.abs(s_lpb0.x - s_pcm.x)) < 5e-3
    assert np.max(np.abs(s_lpb0.s - s_pcm.s)) < 5e-3
    assert np.max(np.abs(s_lpb0.xi - s_pcm.xi)) < 5e-5

    assert abs(s_lpbinf.energy() - s_cosmo.energy()) < 1e-5
    assert np.max(np.abs(s_lpbinf.x - s_cosmo.x)) < 5e-3
    assert np.max(np.abs(s_lpbinf.s - s_cosmo.s)) < 5e-3
    assert np.max(np.abs(s_lpbinf.xi - s_cosmo.xi)) < 5e-5

    assert abs(s_lpb_einf.energy() - s_cosmo.energy()) < 1e-5
    assert np.max(np.abs(s_lpb_einf.x - s_cosmo.x)) < 5e-3
    assert np.max(np.abs(s_lpb_einf.s - s_cosmo.s)) < 5e-3
    assert np.max(np.abs(s_lpb_einf.xi - s_cosmo.xi)) < 5e-5
