import math
import chempy
import sys


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

    coords = [list(map(float, atom)) for atom in xyz]
    return atoms, coords


def quaternion_rotation(angle=math.pi/2, axis=(1.0, 0.0, 0.0), coordinates=[]):
    """
    First compute the center of mass (com).
    Then translate molecule such that the center of mass is in the origin.
    Finally perform the quaternion rotation by "angle" radians around "axis".

    :param angle: float: by how many radians should the molecue be rotated
    :param axis: float: rotate around the given vector
    :return: tuple: rotated x, y, z coordinates
    """
    # Compute the center of mass (com)
    com = [0, 0, 0]
    com[0] = sum([chempy.Substance.from_formula(atom).mass * coord[0] for atom, coord in zip(atoms, coordinates)]) / sum(
        [chempy.Substance.from_formula(atom).mass for atom in atoms])
    com[1] = sum([chempy.Substance.from_formula(atom).mass * coord[1] for atom, coord in zip(atoms, coordinates)]) / sum(
        [chempy.Substance.from_formula(atom).mass for atom in atoms])
    com[2] = sum([chempy.Substance.from_formula(atom).mass * coord[2] for atom, coord in zip(atoms, coordinates)]) / sum(
        [chempy.Substance.from_formula(atom).mass for atom in atoms])

    # Translate atomic coordinates so that COM is at the origo
    coords_translated = [[0, 0, 0] for i in range(len(coordinates))]
    for i, line in enumerate(coords_translated):
        for j in range(3):
            if j == 0:    # x
                coords_translated[i][0] = coords[i][j] - com[0]
            elif j == 1:  # y
                coords_translated[i][1] = coords[i][j] - com[1]
            elif j == 2:  # z
                coords_translated[i][2] = coords[i][j] - com[2]

    # Some shorthands for implementation of rotation formula
    gamma = math.cos(angle/2)
    theta = math.sin(angle/2)
    x = axis[0]  # Axis of rotation, x element
    y = axis[1]  # Axis of rotation, y element
    z = axis[2]  # Axis of rotation, z element

    # Perform the quaternion rotation
    coords_rot = [[] for i in range(len(coordinates))]
    for i, atom in enumerate(coordinates):
        px = atom[0]  # x element of vector to rotate
        py = atom[1]  # y element of vector to rotate
        pz = atom[2]  # z element of vector to rotate

        # Compute the x, y, and z elements of the rotated vector
        px_rot = px*theta**2 * (x**2 - y**2 - z**2) + 2*x*theta**2 * (y*py + z*pz) + 2*gamma*theta * (y*pz - z*py) \
                 + px*gamma**2

        py_rot = py*theta**2 * (y**2 - z**2 - x**2) + 2*y*theta**2 * (x*px + z*pz) + 2*gamma*theta * (z*px - x*pz) \
                 + py*gamma**2

        pz_rot = pz*theta**2 * (z**2 - y**2 - x**2) + 2*z*theta**2 * (x*px + y*py) + 2*gamma*theta * (x*py - y*px) \
                 + pz*gamma**2

        coords_rot[i] = [px_rot, py_rot, pz_rot]

    return coords_rot


help = ["--help", "-help", "-h", "-H", "h", "H", "help", "Help", "HELP", "--HELP"]
msg = """
This script generates 27 rotated conformers of the given molecule, 9 in each of the three dimensions.
It assumes the coordinates are given in a standard XYZ file.

Because the grid's rotational variance is periodic every 90 degrees, we limit the rotations to 90 degrees.
See Wheeler et al (2019, ChemRxiv) for more information: https://doi.org/10.26434/chemrxiv.8864204.v5

You can change the number of conformers to generate by editing the variable "INCREMEMENT" this script.
"""

if sys.argv[1] in help:
    sys.exit(msg)

# Load XYZ file
xyzfile = sys.argv[1]
assert xyzfile.endswith(".xyz"), "The molecular structure must be in an XYZ file: https://en.wikipedia.org/wiki/XYZ_file_format"
atoms, coords = load_xyz(xyzfile)

# Generate rotations in increments around x, y, and z axis
dims = {"x": [1, 0, 0],
        "y": [0, 1, 0],
        "z": [0, 0, 1]}

INCREMENT = 10  # This defines the number of rotations to perform in each dimension: n = 90 / INCREMENT
assert isinstance(INCREMENT, int), "The increment must be an integer!"

for dim in dims:
    for angle in range(0, 90, INCREMENT):
        rad = angle * math.pi / 180
        coords_rot = quaternion_rotation(angle=rad, coordinates=coords, axis=dims[dim])

        outputname = xyzfile.split(".")[0] + f"_{dim}_{angle}.xyz"
        with open(outputname, "w") as f:
            f.write(f"{len(coords_rot)}\n")
            f.write(f"Rotated by {angle} degrees around {dim} axis\n")
            for atom, coord in zip(atoms, coords_rot):
                f.write(f"{atom} {' '.join(list(map(str, coord)))}\n")

print(f"Number of rotational conformers generated: {3 * len(range(0, 90, INCREMENT))}")