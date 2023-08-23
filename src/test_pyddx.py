import pyddx
import numpy as np

tobohr = 1 / 0.52917721092


def test_reference_pcm():
    # This is exactly the same setup as tests/data/ddpcm_force_fmm.in
    # and the reference values are taken from tests/data/ddpcm_force_fmm.out

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
         [ 3.5183990823008846E-09,  1.1275703611421737E-05,  5.0328483106444900E-06],  # noqa: E126,E201,E501
         [ 3.5183990821928763E-09,  1.1275703611421761E-05, -5.0328483106445104E-06],  # noqa: E126,E201,E501
         [-1.1746398344638896E-19, -1.5253454073801846E-20, -1.2378811472579764E-05],  # noqa: E126,E201,E501
         [-3.5183990820914780E-09, -1.1275703611421757E-05, -5.0328483106445121E-06],  # noqa: E126,E201,E501
         [-3.5183990819739290E-09, -1.1275703611421764E-05,  5.0328483106444892E-06],  # noqa: E126,E201,E501
         [-2.8269124554634871E-19,  8.4593913117514543E-21,  1.2378811472579792E-05],  # noqa: E126,E201,E501
         [-4.9461915004768456E-09, -1.2439161462118743E-05, -7.1078808420645150E-06],  # noqa: E126,E201,E501
         [-4.9461915004767157E-09, -1.2439161462118744E-05,  7.1078808420645006E-06],  # noqa: E126,E201,E501
         [-6.9779879053998397E-21,  5.3074464425653005E-21,  1.4376717119412496E-05],  # noqa: E126,E201,E501
         [ 4.9461915004680485E-09,  1.2439161462118765E-05,  7.1078808420645040E-06],  # noqa: E126,E201,E501
         [ 4.9461915004679426E-09,  1.2439161462118765E-05, -7.1078808420645032E-06],  # noqa: E126,E201,E501
         [-2.0205380608822733E-21, -1.4299794681478375E-21, -1.4376717119412485E-05],  # noqa: E126,E201,E501
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

    force = state.solvation_force_terms(solute_field)
    force += state.multipole_force_terms(solute_multipoles);
    assert abs(state.energy() - ref) < 5e-9
    assert np.max(np.abs(force - ref_force)) < 1e-8


def test_reference_lpb():
    # This is exactly the same setup as tests/data/ddlpb_force.txt
    # and reference values are taken from tests/ddlpb_esolv.f90

    charges = np.array([0.08130, -0.20542, 0.06206, 0.06206])
    rvdw = tobohr * np.array([4.00253, 3.63772, 2.99956, 2.99956])
    centres = tobohr * np.array([
        [0.00000,  0.00000,  -1.16691],  # noqa: E201
        [0.00000,  0.00000,   1.42203],  # noqa: E201
        [0.00000,  1.76689,  -2.18736],  # noqa: E201
        [0.00000, -1.76689,  -2.18736],  # noqa: E201
    ]).T

    def assert_lpb(ref, solvent_epsilon=78.3553, solvent_kappa=0.104):
        model = pyddx.Model("lpb", centres, rvdw, solvent_epsilon=solvent_epsilon,
                            solvent_kappa=solvent_kappa, lmax=7, n_lebedev=590,
                            fmm_multipole_lmax=20, fmm_local_lmax=20)

        solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
        solute_field = model.multipole_electrostatics(solute_multipoles)
        solute_psi = model.multipole_psi(solute_multipoles)

        tol = 1e-8
        state = pyddx.State(model, solute_psi, solute_field["phi"], solute_field["e"])
        state.fill_guess(tol)
        state.solve(tol)

        assert abs(state.energy() - ref) < 5e-8

    assert_lpb(-5.3280230267698165E-004, solvent_epsilon=2.0)
    assert_lpb(-1.0243211234919017E-003, solvent_epsilon=200.0)
    assert_lpb(-1.0233753433615445E-003, solvent_kappa=1 / 2)
    assert_lpb(-1.0162805489436332E-003, solvent_kappa=1 / 8)


def test_kappa():
    sodium_chloride = [(-1, 0.1), (+1, 0.1)]
    res = pyddx.solvent_kappa(sodium_chloride, solvent_epsilon=78.54, temperature=298.15)
    ref = 0.05499498832040873
    assert abs(res - ref) < 1e-8

    # CuSO4 0.3 mol/L
    copper_sulfate = [(-2, 0.3), (+2, 0.3)]
    res = pyddx.solvent_kappa(copper_sulfate, solvent_epsilon=78.54, temperature=298.15)
    ref = 0.1905082278652098
    assert abs(res - ref) < 1e-8

def test_ddrun():
    # This is exactly the same setup as tests/data/ddpcm_force_fmm.in
    # and the reference values are taken from tests/data/ddpcm_force_fmm.out

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
         [ 3.5183990823008846E-09,  1.1275703611421737E-05,  5.0328483106444900E-06],  # noqa: E126,E201,E501
         [ 3.5183990821928763E-09,  1.1275703611421761E-05, -5.0328483106445104E-06],  # noqa: E126,E201,E501
         [-1.1746398344638896E-19, -1.5253454073801846E-20, -1.2378811472579764E-05],  # noqa: E126,E201,E501
         [-3.5183990820914780E-09, -1.1275703611421757E-05, -5.0328483106445121E-06],  # noqa: E126,E201,E501
         [-3.5183990819739290E-09, -1.1275703611421764E-05,  5.0328483106444892E-06],  # noqa: E126,E201,E501
         [-2.8269124554634871E-19,  8.4593913117514543E-21,  1.2378811472579792E-05],  # noqa: E126,E201,E501
         [-4.9461915004768456E-09, -1.2439161462118743E-05, -7.1078808420645150E-06],  # noqa: E126,E201,E501
         [-4.9461915004767157E-09, -1.2439161462118744E-05,  7.1078808420645006E-06],  # noqa: E126,E201,E501
         [-6.9779879053998397E-21,  5.3074464425653005E-21,  1.4376717119412496E-05],  # noqa: E126,E201,E501
         [ 4.9461915004680485E-09,  1.2439161462118765E-05,  7.1078808420645040E-06],  # noqa: E126,E201,E501
         [ 4.9461915004679426E-09,  1.2439161462118765E-05, -7.1078808420645032E-06],  # noqa: E126,E201,E501
         [-2.0205380608822733E-21, -1.4299794681478375E-21, -1.4376717119412485E-05],  # noqa: E126,E201,E501
    ]).T

    model = pyddx.Model("pcm", centres, rvdw, solvent_epsilon=78.3553, lmax=10)
    solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
    solute_field = model.multipole_electrostatics(solute_multipoles)
    solute_psi = model.multipole_psi(solute_multipoles)

    state = pyddx.State(model, solute_psi, solute_field["phi"])
    energy, force = state.ddrun(solute_field)
    force += state.multipole_force_terms(solute_multipoles);

    assert abs(energy - ref) < 5e-9
    assert np.max(np.abs(force - ref_force)) < 1e-8
