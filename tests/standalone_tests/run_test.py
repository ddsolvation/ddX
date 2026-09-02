#!/usr/bin/python3
import sys
import os
import argparse
import subprocess

threshold = 1e-4

def list_diff(list1, list2):
    diff = []
    for i in range(len(list1)):
        tmp = []
        for j in range(len(list1[0])):
            tmp.append(list1[i][j] - list2[i][j])
        diff.append(tmp)
    return diff

def inf_norm(input_list):
    max_val = 0.0
    for i in range(len(input_list)):
        for j in range(len(input_list[0])):
            val = abs(input_list[i][j])
            if val > max_val:
                max_val = val
    return max_val

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
    return energy, forces


parser = argparse.ArgumentParser()
parser.add_argument("basename")
parser.add_argument("--exe", help="path to the ddx_driver executable to run")
parser.add_argument("--fpm", action="store_true",
                     help="invoke via 'fpm run' instead of a direct executable path")
args = parser.parse_args()

basename = args.basename
script_dir = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(script_dir, basename + ".txt")
output_file = os.path.join(script_dir, basename + ".log")
ref_file = os.path.join(script_dir, basename + ".ref")

if args.fpm:
    cmd = f"fpm run --target ddx_driver -- {input_file} > {output_file}"
elif args.exe:
    cmd = f'"{args.exe}" {input_file} > {output_file}'
else:
    parser.error("either --exe PATH or --fpm must be given")

status = subprocess.run(cmd, shell=True)
if status.returncode != 0:
    print(f"ERROR: command failed with exit code {status.returncode}")
    sys.exit(status.returncode)

energy, forces = read_log(output_file)
ref_energy, ref_forces = read_log(ref_file)

print(f"Energy:         {energy:20.10f}")
print(f"Ref. energy:    {ref_energy:20.10f}")
assert (energy - ref_energy)/ref_energy < threshold

force_diff = list_diff(forces, ref_forces)
force_max_diff = inf_norm(force_diff)
force_max_ref = inf_norm(ref_forces)
print(f"Force max diff: {force_max_diff:20.10f}")
assert force_max_diff/force_max_ref < threshold
