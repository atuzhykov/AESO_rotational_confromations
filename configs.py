import os

from ase import units

molecules_data = {
    "epoxy_group": {
        "mol_path": os.path.join("data", "mol_files", "epoxy_group.pdb"),
        "axis_bond": [1, 9],
        "atoms_to_rotate": [0],
        "atoms_to_fix": [3, 5, 4, 6, 13, 11, 12, 14],
        "cell": 15
    },
    "hydroxy_group": {
        "mol_path": os.path.join("data", "mol_files", "hydroxy_group.pdb"),
        "axis_bond": [0, 21],
        "atoms_to_rotate": [23],
        "atoms_to_fix": [3, 13, 5, 4, 6, 14, 15, 16],
        "cell": 15
    },
    "hydroxy_group_without_OH": {
        "mol_path": os.path.join("data", "mol_files", "hydroxy_group_without_OH.pdb"),
        "axis_bond": [1, 11],
        "atoms_to_rotate": [0, 21, 22],
        "cell": 15
    },
    "acrylate_group": {
        "mol_path": os.path.join("data", "mol_files", "acrylate_group.pdb"),
        "axis_bond": [25, 26],
        "atoms_to_rotate": [26, 27, 28],
        "atoms_to_fix": [17, 18, 19, 16, 7, 9, 8, 6],
        "cell": 15
    },
    "acrylate_group_1": {
        "mol_path": os.path.join("data", "mol_files", "acrylate_group.pdb"),
        "axis_bond": [25, 0],
        "atoms_to_rotate": [25, 27, 28, 29, 26],
        "atoms_to_fix": [17, 18, 19, 16, 7, 9, 8, 6],
        "cell": 25
    },
    "acrylate_group_1_unfix": {
        "mol_path": os.path.join("data", "mol_files", "acrylate_group.pdb"),
        "axis_bond": [25, 0],
        "atoms_to_rotate": [25, 27, 28, 29, 26],
        "atoms_to_fix": [],
        "cell": 25
    },
    "acrylate_group_2": {
        "mol_path": os.path.join("data", "mol_files", "acrylate_group.pdb"),
        "axis_bond": [3, 2],
        "atoms_to_rotate": [0, 1, 25, 27, 28, 29, 26],
        "atoms_to_fix": [17, 18, 19, 16, 7, 9, 8, 6],
        "cell": 25
    },
    "acrylate_group_3": {
        "mol_path": os.path.join("data", "mol_files", "acrylate_group.pdb"),
        "axis_bond": [3, 2],
        "atoms_to_rotate": [0, 1, 2, 25, 27, 28, 29, 26],
        "atoms_to_fix": [17, 18, 19, 16, 7, 9, 8, 6],
        "cell": 25
    },
    "acrylate_group_excited": {
        "mol_path": os.path.join("data", "mol_files", "acrylate_group_excited.pdb"),
        "axis_bond": [25, 26],
        "atoms_to_rotate": [26, 27, 28, 31],
        "atoms_to_fix": [17, 18, 19, 16, 7, 9, 8, 6],
        "cell": 25
    },
    "ethane": {
        "mol_path": os.path.join("data", "mol_files", "ethane.pdb"),
        "axis_bond": [1, 0],
        "atoms_to_rotate": [1, 2, 3, 4],
        "atoms_to_fix": [0, 5, 6, 7],
        "cell": 10
    },
    "ethylene": {
        "mol_path": os.path.join("data", "mol_files", "ethylene.pdb"),
        "axis_bond": [0, 1],
        "atoms_to_rotate": [2, 3, 1],
        "atoms_to_fix": [5, 4, 0],
        "cell": 10},

    "acrylate_group_based_photopolymer_3_C": {
        "mol_path": os.path.join("data", "mol_files",
                                 "acrylate_group_based_photopolymerization_shortened_3_carbons.pdb"),
        "axis_bond": [26, 27],
        "atoms_to_rotate": [27, 32, 33, 34],
        "atoms_to_fix": [6, 7, 8, 16, 17, 18, 9, 19],
        "cell": 15},

    "acrylate_group_based_photopolymer_3_C_exp_2": {
        "mol_path": os.path.join("data", "mol_files",
                                 "acrylate_group_based_photopolymerization_shortened_3_carbons.pdb"),
        "axis_bond": [25, 26],
        "atoms_to_rotate": [27, 32, 33, 34, 26, 28, 29],
        "atoms_to_fix": [6, 7, 8, 16, 17, 18, 9, 19],
        "cell": 25},
    "acrylate_group_based_photopolymer_3_C_exp_2_unfix": {
        "mol_path": os.path.join("data", "mol_files",
                                 "acrylate_group_based_photopolymerization_shortened_3_carbons.pdb"),
        "axis_bond": [25, 26],
        "atoms_to_rotate": [27, 32, 33, 34, 26, 28, 29],
        "atoms_to_fix": [],
        "cell": 25},
    "acrylate_group_based_photopolymer_3_C_exp_3": {
        "mol_path": os.path.join("data", "mol_files",
                                 "acrylate_group_based_photopolymerization_shortened_3_carbons.pdb"),
        "axis_bond": [25, 0],
        "atoms_to_rotate": [27, 32, 33, 34, 26, 28, 29, 25, 30, 31],
        "atoms_to_fix": [6, 7, 8, 16, 17, 18, 9, 19],
        "cell": 15},
    "acrylate_group_based_photopolymer_3_C_exp_4": {
        "mol_path": "data\\mol_files\\acrylate_group_based_photopolymerization_shortened_3_carbons.pdb",
        "axis_bond": [0, 2],
        "atoms_to_rotate": [27, 32, 33, 34, 26, 28, 29, 25, 30, 31, 0, 1],
        "atoms_to_fix": [6, 7, 8, 16, 17, 18, 9, 19],
        "cell": 15},
    "alpha relax 1":{
        "mol_path": os.path.join("data", "mol_files","alpha_relaxation_1.pdb"),
        "axis_bond": [0,9],
        "atoms_to_rotate": [8, 19, 20],
        "atoms_to_fix": [2, 3, 4, 5, 11, 12, 13, 14],
        "cell": 25},
    "alpha relax 2":{"mol_path": "data\\mol_files\\alpha_relaxation_2.pdb",
        "axis_bond": [0, 5],
        "atoms_to_rotate": [4, 18, 19],
        "atoms_to_fix": [15, 16, 17, 3, 26, 27, 28, 8],
        "cell": 30},

    "alpha relax 3":{"mol_path": "data\\mol_files\\alpha_relaxation_3.pdb",
        "axis_bond": [0, 6],
        "atoms_to_rotate": [5, 22, 23],
        "atoms_to_fix": [19, 20, 21, 4, 10, 32, 33, 34],
        "cell": 30},
    "alpha relax 4":{"mol_path": "data\\mol_files\\alpha_relaxation_4.pdb",
        "axis_bond": [0, 7],
        "atoms_to_rotate": [6, 26, 27],
        "atoms_to_fix": [23, 24, 25, 5, 12, 38, 39, 40],
        "cell": 35},
    "alpha relax 5":{
        "mol_path": os.path.join("data", "mol_files","alpha_relaxation_5.pdb"),
        "axis_bond": [19,21],
        "atoms_to_rotate": [20, 51, 52],
        "atoms_to_fix": [5, 6, 7, 63, 32, 33, 34, 64],
        "cell": 35},
    "alpha relax 6":{
        "mol_path": os.path.join("data", "mol_files","alpha_relaxation_6.pdb"),
        "axis_bond": [20,22],
        "atoms_to_rotate": [21, 37, 38],
        "atoms_to_fix": [8, 49, 50, 51, 34, 68, 69, 70],
        "cell": 60},
    "methyl_main_chain": {
        "mol_path": os.path.join("data", "mol_files","alpha_relaxation_4.pdb"),
        "axis_bond": [11, 12],
        "atoms_to_rotate": [12, 38, 39, 40],
        "atoms_to_fix": [23, 24, 25, 5],
        "cell": 20},
    "methyl_main_chain_shorter":{
        "mol_path": os.path.join("data", "mol_files","alpha_relaxation_1.pdb"),
        "axis_bond": [1,2],
        "atoms_to_rotate": [2, 3, 4, 5],
        "atoms_to_fix": [11, 12, 13, 14],
        "cell": 15},
    "methyl_main_chain_shorter_2":{
        "mol_path": os.path.join("data", "mol_files","alpha_relaxation_2.pdb"),
        "axis_bond": [3,2],
        "atoms_to_rotate": [3, 15, 16, 17],
        "atoms_to_fix": [9, 27, 28, 26],
        "cell": 15},
    }

TOTAL_ANGLE = 360
angle_step = 30

# Molecular Dynamic setting
# The desired temperature, in Kelvin.
temperature_K = 300
# A friction coefficient in inverse ASE time units.
friction = 5e-3
# The time step in ASE time units.
timestep = 2.0 * units.fs
# Number of molecular dynamics steps to be run.
steps = 50

# FIRE settings
# Convergence criterion of the forces on atoms.
fmax = 0.05
