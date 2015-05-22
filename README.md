distforce
=========

This repository contains a fix for LAMMPS that pseudo integrates and distributes forces on virtual COM particles onto physical atoms of the same molecule. 
This fix was designed for the Febuary 2014 distribution of LAMMPS. 

To utilize this code, add COM mass sites to your data file. Give them a small mass, but not nonzero. This will cause LAMMPS to break.
In the input file, set the charge and LJ parameters to zero. Also make sure not to give the COM vitrual sites no velocity.
Before the run command in the input write:

fix fix-id distforce group-id nmolecule

the group-id should be the virtual sites.

Todo:

Improve fix to allow for multiple virtual sites per molecule. 
Have mass of virtual sites be the total mass

New fix:

gatherforce
