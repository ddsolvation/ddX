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


def test_multipoles_functions():
    model = pyddx.Model("cosmo", centres, rvdw, solvent_epsilon=80, enable_fmm=False)

    solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
    solute_field_0 = model.multipole_electrostatics(solute_multipoles,
                                                    derivative_order=0)
    solute_field_1 = model.multipole_electrostatics(solute_multipoles,
                                                    derivative_order=1)
    solute_field_2 = model.multipole_electrostatics(solute_multipoles,
                                                    derivative_order=2)

    assert np.max(np.abs(solute_field_0["phi"] - solute_field_1["phi"])) < 1e-12
    assert np.max(np.abs(solute_field_0["phi"] - solute_field_2["phi"])) < 1e-12

    assert np.max(np.abs(solute_field_1["e"] - solute_field_1["e"])) < 1e-12

    model_fmm = pyddx.Model("cosmo", centres, rvdw, solvent_epsilon=80, enable_fmm=True)
    solute_field_fmm_0 = model_fmm.multipole_electrostatics(solute_multipoles,
                                                            derivative_order=0)
    solute_field_fmm_1 = model_fmm.multipole_electrostatics(solute_multipoles,
                                                            derivative_order=1)
    solute_field_fmm_2 = model_fmm.multipole_electrostatics(solute_multipoles,
                                                            derivative_order=2)

    assert np.max(np.abs(solute_field_0["phi"] - solute_field_fmm_0["phi"])) < 1e-12
    assert np.max(np.abs(solute_field_1["phi"] - solute_field_fmm_1["phi"])) < 1e-12
    assert np.max(np.abs(solute_field_1["e"] - solute_field_fmm_1["e"])) < 1e-12
    assert np.max(np.abs(solute_field_2["phi"] - solute_field_fmm_2["phi"])) < 1e-12
    assert np.max(np.abs(solute_field_2["e"] - solute_field_fmm_2["e"])) < 1e-12
    assert np.max(np.abs(solute_field_2["g"] - solute_field_fmm_2["g"])) < 1e-12
