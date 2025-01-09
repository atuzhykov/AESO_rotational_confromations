import os

import ase
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from scipy import stats
from scipy.spatial.transform import Rotation as R


def get_gauss_distribution_function(data):
    """
    Estimates the Gaussian distribution function from a list of data points.

    Args:
      data (list): A list of numerical data points.

    Returns:
      tuple: A tuple containing the mean (mu), standard deviation (sigma),
             and a function representing the Gaussian distribution (gauss_func).
    """
    mu = np.mean(data)

    sigma = np.std(data, ddof=1)  # Assuming sample standard deviation

    def gauss_func(x):
        return stats.norm.pdf(x, mu, sigma)

    return mu, sigma, gauss_func


def plot_gauss_distribution(data, sigma, gauss_func, path):
    """
    Plots the Gaussian distribution function based on the estimated parameters.

    Args:
      data (list): The original data list.
      mu (float): The estimated mean.
      sigma (float): The estimated standard deviation.
      gauss_func: The function representing the Gaussian distribution.
    """
    plt.clf()

    plt.figure()
    x_min = min(data) - 3 * sigma
    x_max = max(data) + 3 * sigma
    x_vals = np.linspace(x_min, x_max, 100)  # Create 100 points for smooth curve

    # Plot the Gaussian distribution function
    plt.plot(x_vals, gauss_func(x_vals), label="Gaussian Distribution")

    plt.hist(data, density=True, alpha=0.5, label="Data Histogram")

    plt.xlabel("x")
    plt.ylabel("Probability Density")
    plt.title("Gaussian Distribution for Data")
    plt.legend()
    plt.grid(True)
    plt.savefig(path)


def record_energy(atoms, energies):
    energy = atoms.get_potential_energy()
    energies.append(energy)


def write_traj_xyz_and_plot_energy(traj_path, output_file, energies):
    # Load the trajectory file
    atoms_list = read(traj_path, ':')
    with open(output_file, 'w') as f:

        for i, atoms in enumerate(atoms_list):
            f.write(f"{len(atoms)}\n")
            f.write(f'Lattice="{atoms.cell[0][0]} {atoms.cell[0][1]} {atoms.cell[0][2]} '
                    f'{atoms.cell[1][0]} {atoms.cell[1][1]} {atoms.cell[1][2]} '
                    f'{atoms.cell[2][0]} {atoms.cell[2][1]} {atoms.cell[2][2]}" '
                    f'Properties=species:S:1:pos:R:3 pbc="T T T"\n')
            for atom in atoms:
                f.write(f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n")

    # Plotting the energy profile
    min = np.min([energies])
    energies = [item - min for item in energies]
    plt.figure(figsize=(10, 5))
    plt.plot(energies, '-o', label='Energy')
    plt.xlabel('Step')
    plt.ylabel('Energy (eV)')
    plt.title('Energy Profile Over Trajectory')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file.replace('.xyz', '_energy.png'))


def rotate_submolecule(ase_atoms, angle_degrees, axis_atom_indices, atom_to_rotate_index, cell):
    positions = ase_atoms.positions
    pos1 = positions[axis_atom_indices[0]]
    pos2 = positions[axis_atom_indices[1]]
    midpoint = (pos1 + pos2) / 2

    rotation_axis = pos1 - pos2
    rotation_axis /= np.linalg.norm(rotation_axis)
    rotation_axis /= np.linalg.norm(rotation_axis)
    pos = positions[atom_to_rotate_index] - midpoint

    rotation = R.from_rotvec(np.radians(angle_degrees) * rotation_axis)
    rotated_pos_oxygen = rotation.apply(pos) + midpoint

    positions[atom_to_rotate_index] = rotated_pos_oxygen
    ase_atoms.set_positions(positions)
    return ase.Atoms(symbols=ase_atoms.symbols, positions=positions, pbc=True, cell=[cell, cell, cell])
