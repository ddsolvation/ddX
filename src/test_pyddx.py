import pyddx
import numpy as np


def test_reference_pcm():
    tobohr = 1 / 0.52917721092

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

    ref = -0.00017974013712832552
    ref_force = np.array([
         [ 3.5183990823008846E-09,  1.1275703611421737E-05,  5.0328483106444900E-06],
         [ 3.5183990821928763E-09,  1.1275703611421761E-05, -5.0328483106445104E-06],
         [-1.1746398344638896E-19, -1.5253454073801846E-20, -1.2378811472579764E-05],
         [-3.5183990820914780E-09, -1.1275703611421757E-05, -5.0328483106445121E-06],
         [-3.5183990819739290E-09, -1.1275703611421764E-05,  5.0328483106444892E-06],
         [-2.8269124554634871E-19,  8.4593913117514543E-21,  1.2378811472579792E-05],
         [-4.9461915004768456E-09, -1.2439161462118743E-05, -7.1078808420645150E-06],
         [-4.9461915004767157E-09, -1.2439161462118744E-05,  7.1078808420645006E-06],
         [-6.9779879053998397E-21,  5.3074464425653005E-21,  1.4376717119412496E-05],
         [ 4.9461915004680485E-09,  1.2439161462118765E-05,  7.1078808420645040E-06],
         [ 4.9461915004679426E-09,  1.2439161462118765E-05, -7.1078808420645032E-06],
         [-2.0205380608822733E-21, -1.4299794681478375E-21, -1.4376717119412485E-05],
    ]).T

    model = pyddx.Model("pcm", centres, rvdw, solvent_epsilon=78.3553, lmax=10)
    solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
    solute_field = model.multipole_electrostatics(solute_multipoles)
    solute_psi = model.multipole_psi(solute_multipoles)

    state = pyddx.State(model, solute_psi, solute_field["phi"])
    state.fill_guess()
    state.fill_guess_adjoint()
    state.solve()
    state.solve_adjoint()

    force = state.solvation_force_terms()
    assert abs(state.energy() - ref) < 5e-9
    assert np.max(np.abs(force - ref_force)) < 1e-5


def test_reference_lpb():
    tobohr = 1 / 0.52917721092

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

    ref = -0.00018653364236802245  # XXX entirely unverified

    model = pyddx.Model("lpb", centres, rvdw, solvent_epsilon=78.3553, solvent_kappa=0.104, lmax=10)
    solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
    solute_field = model.multipole_electrostatics(solute_multipoles)
    solute_psi = model.multipole_psi(solute_multipoles)

    state = pyddx.State(model, solute_psi, solute_field["phi"], solute_field["e"])
    state.fill_guess()
    state.fill_guess_adjoint()
    state.solve()
    state.solve_adjoint()

    assert abs(state.energy() - ref) < 5e-9
