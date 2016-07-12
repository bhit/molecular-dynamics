# molecular-dynamics
Some of my work on molecular dynamics. 

## Notebooks
[Simulating a water drop](http://nbviewer.jupyter.org/github/bhit/molecular-dynamics/blob/master/Simulating%20a%20water%20drop%20in%20molecular%20dynamics.ipynb)

## Tools
[Bond orientation](https://github.com/bhit/molecular-dynamics/blob/master/bond_angle.c): a C++ code which calculates the bond orientation of a triatomic molecule with respect to some given vector.

## Water models
This is a collection of files to generate water boxes using `moltemplate`, which generates input files for the Lammps molecular dynamics solver. 

* [TIP4P/2005](https://github.com/bhit/molecular-dynamics/blob/master/tip4p-2005.tar.gz): parameters from the Lammps manual.
* [TIP4P-Ew]: a reparameterisation with Ewald summation, parameters from Mario Orsi, Queen Mary University of London
