// Functions to open input and output files,
// read parameter values from input file and 
// write them in output file.

#include "fisher.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * fichierE;
extern FILE * fichierS;

// opens input file:

void ouvrirFichierE()    
{						 
	fichierE = fopen(fichierLecture,"r");
}


// opens output file:

void ouvrirFichierS()   
{
	fichierS = fopen(fichierEcriture,"a");
}


// reads parameter values from input file,
// returns 1 if end of input file, else returns 0

bool lireFichier(int &Nr, double &sr, int &nr, int &mr, double &sigr, double &ar, double &Qr, 
		double &Ur, int &nbSr, double &Lr, double &mutMr, int &NbGenr, int &minir, int &pasr, int &smpr, int &nbr,
		int &nbrepr, int &nboptr, double &Nnewoptr, int &optionmutr, int &eqr, double &diffr)
{					 
	int x;
	bool term = true;
	do {x = fgetc(fichierE);} while (!((x == '*') || (x == EOF)));
		// each parameter set must start with *
    if (x == EOF){
	cout << "\nEnd of input file\n";
        term = false;
    }
	else
	{
		fscanf(fichierE,"%d ",&Nr);
		fscanf(fichierE,"%lf ",&sr);
		fscanf(fichierE,"%d ",&nr);
		fscanf(fichierE,"%d ",&mr);
		fscanf(fichierE,"%lf ",&sigr);
		fscanf(fichierE,"%lf ",&ar);
		fscanf(fichierE,"%lf ",&Qr);
		fscanf(fichierE,"%lf ",&Ur);
		fscanf(fichierE,"%d ",&nbSr);
		fscanf(fichierE,"%lf ",&Lr);
		fscanf(fichierE,"%lf ",&mutMr);
		fscanf(fichierE,"%d ",&NbGenr);
		fscanf(fichierE,"%d ",&minir);
		fscanf(fichierE,"%d ",&pasr);
		fscanf(fichierE,"%d ",&smpr);
        fscanf(fichierE,"%d ",&nbr);
		fscanf(fichierE,"%d ",&nbrepr);
		fscanf(fichierE,"%d ",&nboptr);
		fscanf(fichierE,"%lf ",&Nnewoptr);
		fscanf(fichierE,"%d ",&optionmutr);
		fscanf(fichierE,"%d ",&eqr);
		fscanf(fichierE,"%lf ",&diffr);
		
		term = true;
	} 
	return term;
}


// writes parameter values in output file:

void ecrireParametres(int Nv, double sv, int nv, int mv, double sigv, double av, double Qv, double Uv, int nbSv, double Lv, double mutMv, int NbGenv, int miniv, int pasv, int smpv, int nbrv,
						int nbrepv, int nboptv, double Nnewoptv, int optionmutv, int eqv, double diffv)
{
	fprintf(fichierS,"\n_________________________________________\n");
	fprintf(fichierS,"\nN = %d", Nv);
	fprintf(fichierS,", s = %g", sv);
	fprintf(fichierS,", n = %d", nv);
	fprintf(fichierS,", m = %d", mv);
	fprintf(fichierS,", sigma = %g", sigv);
	fprintf(fichierS,", alpha = %g", av);
	fprintf(fichierS,", Q = %g", Qv);
	fprintf(fichierS,", U = %g", Uv);
	fprintf(fichierS,", nbS = %d", nbSv);
	fprintf(fichierS,", L = %g", Lv);
	fprintf(fichierS,"\nmutM = %g", mutMv);
	fprintf(fichierS,", generations = %d", NbGenv);
	fprintf(fichierS,", mini = %d", miniv);
	fprintf(fichierS,", pas = %d", pasv);
	fprintf(fichierS,", smp = %d", smpv);
    fprintf(fichierS,", nbr = %d", nbrv);
    fprintf(fichierS,", nbrep = %d", nbrepv);
    fprintf(fichierS,", nbopt = %d", nboptv);
    fprintf(fichierS,", Nnewopt = %g", Nnewoptv);
    fprintf(fichierS,", optionmut = %d", optionmutv);
    fprintf(fichierS,", eq = %d", eqv);
    fprintf(fichierS,", diff = %g", diffv);
}
