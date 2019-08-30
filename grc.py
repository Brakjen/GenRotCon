import math
import numpy as np
import chempy
from pyquaternion import Quaternion


def load_xyz(filename):
    """
    Load the xyz file.
    :param filename: str: path to file
    :return: tuple: ([atomic_labels], [coordinates])
    """
    with open(filename) as f:
        xyz = [line.split() for line in f.readlines()[2:]]

    atoms = []
    for i, line in enumerate(xyz):
        atoms.append(line[0])
        del(xyz[i][0])

    coords = np.asarray([list(map(float, atom)) for atom in xyz])
    return atoms, coords


def center_of_mass(atoms, coords):
    """
    Compute the center of mass.
    :param coords: numpy array
    :return: tuple: x, y, z coordinates of center of mass.
    """
    x = sum([chempy.Substance.from_formula(atom).mass * coord[0] for atom, coord in zip(atoms, coords)]) / sum(
        [chempy.Substance.from_formula(atom).mass for atom in atoms])
    y = sum([chempy.Substance.from_formula(atom).mass * coord[1] for atom, coord in zip(atoms, coords)]) / sum(
        [chempy.Substance.from_formula(atom).mass for atom in atoms])
    z = sum([chempy.Substance.from_formula(atom).mass * coord[2] for atom, coord in zip(atoms, coords)]) / sum(
        [chempy.Substance.from_formula(atom).mass for atom in atoms])

    return x, y, z


def center_molecule(coords, com):
    """
    Translate the  molecular coordinates such that the center of mass is in the origo.
    :param coords:
    :param com:
    :return: translated coordinates
    """
    coords_translated = [[] for i in range(len(coords))]
    for i, line in enumerate(coords_translated):
        for j in range(3):
            if j == 0:    # x
                coords_translated[i].append(coords[i][j] - com[0])
            elif j == 1:  # y
                coords_translated[i].append(coords[i][j] - com[1])
            elif j == 2:  # z
                coords_translated[i].append(coords[i][j] - com[2])

    return coords_translated


def quaternion_rotation(angle=math.pi/2, axis=(1.0, 0.0, 0.0), coordinates=[]):
    """
    Perform the quaternion rotation by "angle" radians around "axis"
    :param angle: float: by how many radians should the molecue be rotated
    :param axis: float: rotate around the given vector
    :return: tuple: rotated x, y, z coordinates
    """
    # Some shorthands
    gamma = math.cos(angle/2)
    theta = math.sin(angle/2)
    x = axis[0]  # Axis of rotation, x element
    y = axis[1]  # Axis of rotation, y element
    z = axis[2]  # Axis of rotation, z element

    # Define the quaternion
    q = (gamma, axis[0]*theta, axis[1]*theta, axis[2]*theta)

    # Construct quaternion representations for atomic coordinates
    #[atom.insert(0, 0.0) for atom in coordinates]

    # Test vector
    coords_rot = []
    for i, atom in enumerate(coordinates):
        px = atom[0]  # x element of vector to rotate
        py = atom[1]  # y element of vector to rotate
        pz = atom[2]  # z element of vector to rotate

        # q p q*
        px_rot = px*theta**2 * (x**2 - y**2 - z**2) + 2*x*theta**2 * (y*py + z*pz) + 2*gamma*theta * (y*pz -z*py)

        py_rot = py*theta**2 * (y**2 - z**2 - x**2) + 2*y*theta**2 * (x*px + z*pz) + 2*gamma*theta * (z*px - x*pz)

        pz_rot = pz*theta**2 * (z**2 - y**2 - x**2) + 2*z*theta**2 * (x*px + y*py) + 2*gamma*theta * (x*py - y*px)

        coords_rot.append([px_rot, py_rot, pz_rot])

    return coords_rot


atoms, coords = load_xyz("test.xyz")
com = center_of_mass(atoms, coords)
coords_centered = center_molecule(coords, com)

# Now rotate by 90 degrees around the x axis
coords_rotated = quaternion_rotation(coordinates=coords_centered)

with open("test_rot.xyz", "w") as f:
    f.write(f"{len(coords_rotated)}\n")
    f.write("\n")
    for atom, coord in zip(atoms, coords_rotated):
        f.write(f"{atom} {' '.join(list(map(str, coord)))}\n")