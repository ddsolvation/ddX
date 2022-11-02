import numpy as np

from pytest import approx

import pyddx


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

    model = pyddx.Model("pcm", charges, centres, rvdw, solvent_epsilon=78.3553, lmax=10)
    nuclear = model.solute_nuclear_contribution()
    state = model.initial_guess()
    state = model.solve(state, nuclear["phi"], tol=1e-10)
    # state = model.adjoint_solve(state, nuclear["psi"], tol=1e-10)
    # TODO Test force

    print(state.get_reaction_potential_at_centers(range(1, 13)))
    print('ana')
    print(state.get_reaction_field_at_centers(range(1, 13)))
    print('num')
    print(state.get_reaction_field_at_centers(range(1, 13), True))

    energy = 0.5 * np.sum(state.x * nuclear["psi"])
    ref = -0.00017974013712832552
    assert abs(energy - ref) < 5e-9

# TODO Test COSMO
if __name__ == '__main__':
    test_reference_pcm()
