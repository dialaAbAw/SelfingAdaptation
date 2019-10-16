// Header file: definitions of global variables, function prototypes

#ifndef FISHER_H
#define FISHER_H

#include <vector>
#include <iostream>
#include "MersenneTwister.h"
using namespace std;


// Global variables:

#define fichierLecture "parametres.txt"     // names of input
#define fichierEcriture "resultats.txt"		// and output files

// "chr": represents a chromosome

struct chr
{
	double mod; 	// neutral locus (infinite sites model, each value in the vector corresponds to one mutation)
	unsigned char * sel; // selected loci (for each locus, holds the number of the allele in the "mutations" table)
	int nbchr;   // number of copies of this genome segment in the population 
};

// "Nall" is used to count allele frequencies at the neutral loci:

struct Nall
{
	double all;
	double freq;
};

// Function prototypes:

void ouvrirFichierE();
void ouvrirFichierS();
void ecrireParametres(int Nv, double sv, int nv, int mv, double sigv, double av, double Qv, double Uv, int nbSv, double Lv, double mutMv, int NbGenv, int miniv, int pasv, int smpv, int nbrv,
						int nbrepv, int nboptv, double Nnewoptv, int optionmutv, int eqv, double diffv);
bool lireFichier(int &Nr, double &sr, int &nr, int &mr, double &sigr, double &ar, double &Qr, 
		double &Ur, int &nbSr, double &Lr, double &mutMr, int &NbGenr, int &minir, int &pasr, int &smpr, int &nbrr,
		int &nbrepr, int &nboptr, double &Nnewoptr, int &optionmutr, int &eqr, double &diffr);
void recursion(int Nt, double sv, int nv, int mv, double sigv, double av, double Qv, double Uv, int nbSv, double Lv, 
			   double mutMv, int NbGenv, int miniv, int pasv, int smpv, int nbrv,
			   int nbrepv, int nboptv, double Nnewoptv, int optionmutv, int eqv, double diffv);
double gammln(const double xx);
double poisdev(const double xm);
double gasdev();
double binldev(const double pp, const int n);
void rec(chr &res, chr &c1, chr &c2, int nbCo, int nS);

#endif
