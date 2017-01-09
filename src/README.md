# IO-scheduling-simu
Source code for heuristics of periodic IO-scheduling

## To compile
make

## To launch simulations
./simu [OUTPUT="results.txt"]

The first optional argument is the name of the output file, which contains all the results (best system efficiency, corresponding dilation and best period)
for all sets of applications tested. If nothing is given as argument, the
simulator will write into the file output.txt.

The machine used by default is a machine with 640 cores (so all the sets of applications use 640 cores at total) with B=3GB and b=0.01GB. However to avoid rounding issues, b should always be set
to 1 and the #define BWRATIO should be set to 1/Real_b. For example to set up a machine with B=8GB and b=0.05GB, set the define BWRATIO to 20, and your machine should be created by Machine(8*BWRATIO,1).

110 sets of applications are provided:
	The first ten are the ones presented in the paper, they are the ones based on real scientific applications (Turbulence1,Turbulence2,AstroPhysics,PlasmaPhysics). To reproduce the results of the paper you need to run them independently to set the value of Tmin to the following ones:
		Set 1 : Tmin = 900
		Set 2 : Tmin = 16000
		Set 3 : Tmin = 16000
		Set 4 : Tmin = 16000
		Set 5 : Tmin = 495000
		Set 6 : Tmin = 16000
		Set 7 : Tmin = 4544
		Set 8 : Tmin = 495000
		Set 9 : Tmin = 16000
		Set 10 : Tmin = 16000
	The next 100 sets were generated randomly, using the provided script generate_sets.py. For these ones, Tmin was chosen exactly to be the longest execution time of an instance of any application without congestion.
To choose which set to use, just change the values of the main loop inside the main function.

In the main.cpp, you can also try to change Tmax and epsilon parameters. We always used Tmax = 10*Tmin and epsilon=1/100.

The function Simulation::print_pattern() can be used to display information on its pattern.

## Results

The application displays in the console the time spent in creating all the schedules for each set.

The application also creates two files, one of the name you chose in the command-line (by default "results.txt") and one file "temp.txt".
The first file contains one line for each set of applications with on it: the number of the set, then the best system efficiency, then the corresponding dilation and finally the corresponding period.
The second file can be used to have more information on what the algorithm did during the last set. You can see every period tried, the system efficiency and dilation as long as the number of instances of each application in the period.
The file "results_paper.txt" contains our results on the 110 sets.
