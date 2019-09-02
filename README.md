# GenRotCon
Generate rotational conformers easily for use in testing of integration grids' rotational variance.

![](animation.gif)

Usage summary: 
```bash
grc.py [-h] -x <path> [-n] [-s <N>] [-a] [-S <N>] [-A {x,y,z}] [-o <filename>]
```
##USAGE
To generate this help message, run 
```bash
$ python grc.py -h
```

To generate 8 rotational conformers in each dimension, giving
a total of 24 conformers, run
```bash
$ python grc.py --xyz <molecule.xyz> --step 10
```

Note that no rotations above 90 degrees are performed, as the grid's
rotational variance is periodic every 90 degrees.

To generate an animation XYZ file of 359 structures rotated around
the x axis, with output name <coolstuff.xyz>, without generating
the rotational conformers, run
```bash
$ python grc.py --norotation --xyz <molecule.xyz> --animation --animationstep 1 --animationaxis x --outputname coolstuff.xyz
```

or the shorthand version

```bash
$ python grc.py -na -x <molecule.xyz> -S 1 -A x -o coolstuff.xyz
```

##REQUIREMENTS PYTHON 3
Python version required: 3.6 or higher <br/>
Libraries needed: argparse version 1.1 <br/>
It assumes the coordinates are given in a standard XYZ file, and that the atomic labels <br/>
in the first column are correctly typed atomic symbols <br/>
(e.g. "H", not "h" or "1"; "Ca", not "ca" or "CA" or "20") 

##REQUIREMENTS PYTHON 2
Python version required: 2.6 or higher <br/>
Libraries needed: argparse version 1.1 <br/>
It assumes the coordinates are given in a standard XYZ file, and that the atomic labels <br/>
in the first column are correctly typed atomic symbols <br/>
(e.g. "H", not "h" or "1"; "Ca", not "ca" or "CA" or "20")

##BACKGROUND INFORMATION
See Wheeler et al (2019, ChemRxiv) for more information
on the integration grid's lack of rotational invariance:
https://doi.org/10.26434/chemrxiv.8864204.v5

##AUTHOR INFORMATION
This script was made by <br/>
Anders Brakestad <br/>
PhD Candidate in Computational Chemistry <br/>
UiT The Arctic University of Tromsø <br/>
anders.m.brakestad@uit.no
