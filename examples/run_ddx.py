import pyddx
import numpy as np

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

print(pyddx.banner())

model = pyddx.Model("pcm", centres, rvdw, solvent_epsilon=78.3553)

# Compute solute contributions (here just charges)
solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
solute_field = model.multipole_electrostatics(solute_multipoles)
solute_psi = model.multipole_psi(solute_multipoles)

# Solve the problem
state = pyddx.State(model, solute_psi, solute_field["phi"])
state.fill_guess()
state.solve()
state.fill_guess_adjoint()
state.solve_adjoint()

# Show results
energy = 0.5 * np.sum(state.x * solute_psi)
force = state.solvation_force_terms()
print(energy)
print(force)
