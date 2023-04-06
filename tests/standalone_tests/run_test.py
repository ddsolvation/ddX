#!/usr/bin/python3
import sys
import os
import numpy as np

threshold = 1e-4

def read_log(path):
    energy = 0.0
    forces = []
    with open(path, "r", encoding="utf-8") as log:
        section = ""
        for line in log:
            if "Solvation energy (Hartree):" in line:
                tokens = line.split()
                energy = float(tokens[3]) 
            if section == 'forces':
                tokens = line.split()
                forces.append([float(tokens[1]), float(tokens[2]), \
                               float(tokens[3])])
            if 'Full forces (kcal/mol/A)' in line:
                section = 'forces'
    return energy, np.array(forces)

basename = sys.argv[1]
input_file = basename + ".txt"
output_file = basename + ".log"
ref_file = basename + ".ref"

os.system(f"./ddx_driver_testing {input_file} > {output_file}")

energy, forces = read_log(output_file)
ref_energy, ref_forces = read_log(ref_file)

print(f"Energy:         {energy:20.10f}")
print(f"Ref. energy:    {ref_energy:20.10f}")
assert (energy - ref_energy)/ref_energy < threshold

force_max_diff = np.max(np.abs(forces - ref_forces))
force_max_ref = np.max(np.abs(ref_forces))
print(f"Force max diff: {force_max_diff:20.10f}")
assert force_max_diff/force_max_ref < threshold


