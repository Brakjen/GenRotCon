import math
import numpy as np
import chempy


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
    # Define the quaternion
    q = (math.cos(angle/2), axis[0]*math.sin(angle/2), axis[1]*math.sin(angle/2), axis[2]*math.sin(angle/2))

    # Construct quaternion representations for atomic coordinates
    [atom.insert(0, 0.0) for atom in coordinates]


    return coordinates


atoms, coords = load_xyz("test.xyz")
com = center_of_mass(atoms, coords)
coords_centered = center_molecule(coords, com)

print(quaternion_rotation(coordinates=coords_centered))

