# GenRotCon
Generate rotational conformers easily for use in testing of integration grids' rotational variance.

![](animation.gif)

Python3 version (3.6 or higher): grc.py <br/>
Python2 version (2.6 or higher): grc27.py

This script uses quaternions to perform spatial rotations of the submitted molecule.

## Requirements
The coordinates must be given in standard XYZ file format, and the atomic labels (first column) must be
the correct atomic symbol (i.e. "Ca" and not "ca" or "20").

## Usage
Run <br/>
$python grc.py <coordinates.xyz> <br/>
to generate 24 rotational conformers, 8 in each dimension. The default is to rotate in range(10, 90, 10),
but the increment can be edited manually if desireable.

For a help message, run <br/>
$ python grc.py --help

To generate an XYZ file (animation.xyz) containing the full rotation in 1 degree increments around the x-axis, run <br/>
$ python grc.py <coordinates.xyz> --animation

## Algorithm
All atomic coordinates are translated such that the center of mass is at the origo.
Then the rotations are carried out, ensuring that the rotational center is the center of mass.
The axes of rotation are relative to the molecule's orientation as submitted by the user.

## Backstory
This script was motivated by the [recent and popular preprint in ChemRxiv](https://chemrxiv.org/articles/Popular_Integration_Grids_Can_Result_in_Large_Errors_in_DFT-Computed_Free_Energies/8864204/5), highlighting a (perhaps) lesser known error present in quantum chemistry calculations that use numerical integration grids. Energies, and especially Gibbs free energies of floppy molecules, are sensitive to the spatial orientation of the molecule, because the integration grid is not rotationally invariant. I wanted a quick and easy way to generate rotational conformers to use in testing, and so this little script was made. There are most likely several codes out there that do the same thing (in a better way), but then I would not have gotten this excursion into the wonderful world of quaternions...
