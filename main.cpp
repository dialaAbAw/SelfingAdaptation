// main() function: reads parameter values in input file,
// broadcast to all processors, initiate sprng random number generator
// and runs the simulation.

#include "fisher.h"
#include "MersenneTwister.h"
#include <string>
#include <iostream>
using namespace std;

// input and output files:

MTRand rnd;
FILE * fichierE;
FILE * fichierS;

int main()
{	
	// definitions of variables:

	int Nt, n, m, nbS, NbGen, mini, pas, smp, nbr, nbrep, nbopt, optionmut, eq;
	double s, sig, a, Q, U, L, mutM, Nnewopt, diff;

	// opens input and output files:

	ouvrirFichierE();
	ouvrirFichierS();
    
    bool lire = true;

    do{

    lire = lireFichier(Nt, s, n, m, sig, a, Q, U, nbS, L, mutM, NbGen, mini, pas, smp, nbr, nbrep, nbopt, Nnewopt, optionmut, eq, diff);
    ecrireParametres(Nt, s, n, m, sig, a, Q, U, nbS, L, mutM, NbGen, mini, pas, smp, nbr, nbrep, nbopt, Nnewopt, optionmut, eq, diff);

    recursion(Nt, s, n, m, sig, a, Q, U, nbS, L, mutM, NbGen, mini, pas, smp, nbr, nbrep, nbopt, Nnewopt, optionmut, eq, diff);
        
    }while(lire);
	
	// closes files:

	fclose(fichierE);
	fclose(fichierS);
	
	return 0;
}
