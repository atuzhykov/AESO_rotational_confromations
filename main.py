import datetime
import json
import os
import shutil

import numpy as np
from ase.constraints import FixAtoms
from ase.io import read, write, Trajectory
from ase.md import Langevin
from ase.optimize import FIRE
from ase.visualize import view
from matplotlib import pyplot as plt

from aeso.configs import TOTAL_ANGLE, molecules_data, angle_step, temperature_K, friction, timestep, steps, fmax
from aeso.utils import rotate_submolecule, get_gauss_distribution_function, plot_gauss_distribution, \
    write_traj_xyz_and_plot_energy, record_energy

engine = "MACE" if os.name == 'nt' else "GPAW"

# Way of structure optimization. True if Langevin MD, False if FIRE optimizer
use_molecular_dynamics = False
# store relaxation checkpoints per angle step
store_checkpoints = False

if engine == "GPAW":
    from gpaw import GPAW, PW

    calc = GPAW(xc='PBE',
                mode=PW(300), symmetry='off', txt=None)
elif engine == "MACE":
    from mace.calculators import mace_off

    calc = mace_off(model="large", device='cuda')

for group in ["acrylate_group_based_photopolymer_3_C_exp_2_unfix"]:
    predictions = []
    mol_path = molecules_data[group]["mol_path"]
    axis_bond = molecules_data[group]["axis_bond"]
    atoms_to_rotate = molecules_data[group]["atoms_to_rotate"]
    atoms_to_fix = molecules_data[group]["atoms_to_fix"]
    cell = molecules_data[group]["cell"]

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    relaxation_method = "Langevin" if use_molecular_dynamics else "FIRE"
    experiment_dir = os.path.join("experiments",
                                  f"{group}_{engine}_{relaxation_method}_step_{angle_step}_degrees_lattice_{cell}_{timestamp}")
    os.makedirs(experiment_dir, exist_ok=True)
    json_path = os.path.join(experiment_dir, f"{group}_info.json")
    with open(json_path, 'w') as json_file:
        json.dump(molecules_data[group], json_file, indent=4)
    shutil.copy(mol_path, experiment_dir)

    ase_atoms = read(mol_path)
    positions = ase_atoms.get_positions()
    center_of_mass = np.mean(positions, axis=0)
    translation = np.array([cell / 2, cell / 2, cell / 2]) - center_of_mass
    positions = [list(np.array(pos) + translation) for pos in positions]
    ase_atoms.set_positions(positions)
    ase_atoms.set_cell([cell, cell, cell])
    ase_atoms.set_pbc(True)
    c = FixAtoms(indices=atoms_to_fix + atoms_to_rotate)
    ase_atoms.set_constraint(c)
    ase_atoms.calc = calc
    view(ase_atoms)
    write(os.path.join(experiment_dir, f'{group}.png'), ase_atoms)
    # initial structure optimization
    if use_molecular_dynamics:
        dyn = Langevin(ase_atoms, timestep, temperature_K=temperature_K, friction=friction)
        dyn.run(steps)
    else:
        dyn = FIRE(ase_atoms)
        dyn.run(fmax=fmax)

    rotated_mols = []
    rotated_mols.append(ase_atoms)
    predictions.append(ase_atoms.get_potential_energy())
    angle = angle_step
    while (angle <= TOTAL_ANGLE):
        energies = []
        print(f"{angle}° from {TOTAL_ANGLE}° with step: {angle_step}° was processed for {group}")
        # view(ase_atoms)
        ase_atoms = rotate_submolecule(ase_atoms, angle_step, axis_bond, atoms_to_rotate, cell)
        ase_atoms.calc = calc
        ase_atoms.set_constraint(c)

        # Structure optimization block: proceed and store energies minimization plot and
        # atoms positions trajectories for OVITO

        if store_checkpoints:
            traj_filename = os.path.join(experiment_dir, f"{engine}_{group}_step_{angle}_trajectory.traj")
            traj = Trajectory(traj_filename, 'w', ase_atoms)
        if use_molecular_dynamics:
            dyn = Langevin(ase_atoms, timestep, temperature_K=temperature_K, friction=friction)
            if store_checkpoints:
                dyn.attach(record_energy, interval=1, atoms=ase_atoms, energies=energies)
                dyn.attach(traj.write, interval=1)
            dyn.run(steps)
        else:
            dyn = FIRE(ase_atoms)
            if store_checkpoints:
                dyn.attach(record_energy, interval=1, atoms=ase_atoms, energies=energies)
                dyn.attach(traj.write, interval=1)
            dyn.run(fmax=fmax)

        if store_checkpoints:
            traj.write(ase_atoms)
            traj.close()
            write_traj_xyz_and_plot_energy(traj_filename,
                                           os.path.join(experiment_dir,
                                                        f"{engine}_{group}_step_{angle}_trajectory.xyz"),
                                           energies)

        rotated_mols.append(ase_atoms.copy())
        predictions.append(ase_atoms.get_potential_energy())

        angle += angle_step
    # for OVITO animation
    write(os.path.join(experiment_dir, f"{engine}_with_{os.path.splitext(os.path.basename(mol_path))[0]}_rotated.xyz"),
          rotated_mols, 'extxyz')
    min = np.min([predictions])
    diff_predictions = [item - min for item in predictions]

    np.savetxt(os.path.join(experiment_dir, f'{engine}_with_energy_barriers_predictions_{group}.txt'), diff_predictions)

    mu, sigma, gauss_func = get_gauss_distribution_function(diff_predictions)
    plot_gauss_distribution(diff_predictions, sigma, gauss_func,
                            os.path.join(experiment_dir, f"{engine}_with_gauss_{group}.png"))
    fig, ax = plt.subplots(figsize=(10, 5))
    angles = np.arange(0, TOTAL_ANGLE + 1, angle_step).tolist()
    ax.plot(angles, diff_predictions, '-o', label='Change in Predictions', color='orange')
    ax.set_xlabel('Angle (degrees)')
    ax.set_ylabel('Change in Predicted Energy (eV)')
    ax.set_title(f'Differential {engine} Predictions for {group}')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(experiment_dir, f"{engine}_{group}_energy_changes.png"))  # Save the figure to file
