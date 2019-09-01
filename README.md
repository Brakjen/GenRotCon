# GenRotCon
Generate rotational conformers easily for use in testing of integration grids' rotational variance.

Python3 version (3.6 or higher): grc.py
Python2 version (2.6 or higher): grc27.py

This script uses quaternions to perform spatial rotations of the submitted molecule.

## Requirements
The coordinates must be given in standard XYZ file format, and the atomic labels (first column) must be
the correct atomic symbol (i.e. "Ca" and not "ca" or "20").

## Usage
Run 
$python grc.py <coordinates.xyz>
to generate 24 rotational conformers, 8 in each dimension. The default is to rotate in range(10, 90, 10),
but the increment can be edited manually if desireable.

For a help message, run
$ python grc.py --help

To generate an XYZ file (animation.xyz) containing the full rotation in 1 degree increments around the x-axis, run
$ python grc.py <coordinates.xyz> --animation

## Algorithm
All atomic coordinates are translated such that the center of mass is at the origo.
Then the rotations are carried out, ensuring that the rotational center is the center of mass.
The axes of rotation are relative to the molecule's orientation as submitted by the user.

![](animation.gif)
