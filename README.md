# SelfingAdaptation
C++ program for simulating a population of partially self-fertilizing individuals forced to adapt from pre-existing standing genetic variation. This program is a multiallelic modified version of the script used in Abu Awad, Diala; Roze, Denis (2018), Effects of partial selfing on the equilibrium genetic variance, mutation load and inbreeding depression under stabilizing selection, Dryad, Dataset, https://doi.org/10.5061/dryad.jk4r6

To compile this program the header file MersenneTwister.h is included. It is based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus, written by Richard J. Wagner.

The file main.cpp initialises the life-cycle, recursion.cpp contains the loop representing the life-cycle (processes of selection, measuring fitness, mutation and reproduction), the rec.cpp file contains the recombination function, the header fisher.h, contains defintions for different structures used, ranbin.cpp contains codes for distributions used within the program (Gaussian distribution, Poisson distribution, etc.), fichiers.cpp contains the code to point to the parameter file and generate output files with summary informatio. An example parameter file parametres.txt is also included. 

To compile this program using GNU you need to type this command in the terminal, making sure you are in the right working directory (i.e. where all the necessary files are):
g++ -o sims *.cpp *.h -lm

"sims" can be replaced by any other name you wish to name the executable. To launch, the command in the termal is simply:
./sims

