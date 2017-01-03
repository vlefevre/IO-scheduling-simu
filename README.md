# IO-scheduling-simu
Source code for heuristics of periodic IO-scheduling

## To compile
make

## To launch simulations
./simu [$OUTPUT="output.txt"]

The first optional argument is the name of the output file, which contains all the results (best system efficiency, corresponding dilation and best period)
for all sets of applications tested. If nothing is given as argument, the
simulator will write into the file output.txt.

!!!TO USE THE LINEAR OPTIMIZATION FUNCTION YOU NEED TO HAVE IBM ILOG CPLEX INSTALLED!!!

In main.cpp are defined 4 applications (PlasmaPhysics, AstroPhysics, Turbulence1, Turbulence2), and one set of applications.
You can add applications and change the set, then re-compile to test with different applications.
In main.cpp are also defined the 3 parameters Tmin, Tmax and epsilon.

The function Simulation::print_pattern() can be used to display information on its pattern.
