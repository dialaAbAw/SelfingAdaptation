#include "fisher.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <climits>
#include <cmath>
#include <csignal>
#include <vector>
#include <algorithm>
using namespace std;

extern MTRand rnd;
extern FILE * fichierE;
extern FILE * fichierS;


/*----------------------------------------------------------
Function recursion: iterates the life cycle.
Parameters are:
Nv: population size (nb of diploid organisms)
sv: selfing rate
nv: number of phenotypic dimensions
mv: pleiotropy
sigv: size of mutational steps along each phenotypic axis
om2v: strength of selection
Qv: curvature of fitness function (cf Gros et al 2009 Genetics)
Uv: mutation rate per genome (on selected loci)
nbSv: number of selected loci per genome
Lv: genome map length (mean nb of cross-overs per meiosis)
mutMv: mutation rate at neutral locus
NbGenv: number of generations
pasv: time interval between measurements of neutral diversity and writing in result file
nbrepv : number of repetitions for parameter set
nbopt : number of optimum with new values after equilibrium
Nnewopt : new value of optimum after equilibrium
optionmut : keep or not mutations after changing some optimums (if = 1, no mutations anymore)
eqv: Number of generations for testing if population is at equilibrium
diffv: The differences between means fitness you allow between the last 1000 generations and the 1000 previous ones
-----------------------------------------------------------*/


void recursion(int Nv, double sv, int nv, int mv, double sigv, double om2v, double Qv, double Uv, int nbSv, double Lv,
		double mutMv, int NbGenv, int miniv, int pasv, int smpv, int nbr, int nbrepv, int nbopt, double Nnewoptv, int optionmut, int eqv,double diffv)
{
	 // Keeps track of number of repetitions
    int rep=0;

    do{
        rep++;
        
	// variables:

	vector<Nall> Fr;
	Nall allTemp;
    int i, j, k, nm, gen, mut, chr1, ind, nb, nb2, nb3, nb4, parent, prem1, par1, frs,
	    free, site, ns, part, cmpt, nbCo;
	double w, wbar, wmax, varw, rd, rd2, div, pp, hom, d, x, muM, sbar, lnsbar, sbmoy, lnwbar, lnvarw, le, divs, divsTot, Ut;
	double wout, lnwout, wself, lnwself,Wt, Wtb;

	// various fixed quantities:
	
	int twoN = 2*Nv;
	int fourN = 4*Nv;
	int twoN_1 = 2*Nv - 1;
	int nbS_1 = nbSv - 1;
	int nS = nbSv * nv;
	int twonN = 2 * Nv * nv;
	int twon = 2 * nv;
    int np1 = nv + 1;
	double hQ = Qv / 2.0;
    int twoSmp = 2 * smpv;
	int twonSmp = 2 * nv * smpv;
	vector<int> store;
    double VE = 1.0 / double(nv);
    //double VE = 0;
    double sqrn = sqrt(VE);
    double VS = om2v + VE;
    double twoVS = 2.0 * VS;
    double twoom2 = 2.0 * om2v;
    int uchmax1 = UCHAR_MAX + 1;
    double sigs = sqrt(sigv);
    
    // creates result file (with parameter values in file name):
    
    char nomFichier[256];
    stringstream nomF;
    nomF << "result_N" << Nv << "_s" << sv << "_n" << nv << "_m" << mv << "_sig" << sigv << "_om2_" << om2v << "_Q" << Qv
    << "_U" << Uv << "_nbS" << nbSv << "_L" << Lv << "_u" << mutMv << "_nb"<<nbr<<"_rep" << rep<<".txt";
    nomF >> nomFichier;
    ofstream fout;
    

     
     /*char nomFichier3[256];
     stringstream nomF3;
     nomF3 << "effets_N" << Nv << "_s" << sv << "_n" << nv << "_m" << mv << "_sig" << sigv << "_om2" << om2v << "_Q" << Qv
     << "_U" << Uv << "_nbS" << nbSv << "_L" << Lv << "_u" << mutMv << ".txt";
     nomF3 >> nomFichier3;
     ofstream fout3;*/
    
    fout.open(nomFichier);
    
    //fout3.open(nomFichier3);

	// Table "Chrm" holds all the types of chromosomes present in the population:
	
	chr * Chrm = new chr [fourN];
	
	// table "pop" contains 2N pointers to "chr" elements of the Chrm table
	// (one pointer for each chromosome in the population).
	// "temp" is used during the formation of the next generation:

	chr ** pop = new chr *[twoN];
	chr ** temp = new chr *[twoN];

    // "mutations" will hold the phenotypic effects of mutant alleles at each locus
	// on each phenotypic axis:
	
    vector<vector<float *> > mutations;
    vector<float *> locus;
    
    for (i = 0; i < nbSv; i++) // one element per locus
        mutations.push_back(locus);
    
    // initialization: only one allele exists at each locus, with effects 0 on all traits
    
    for (i = 0; i < nbSv; i++)
    {
        float * all0 = new float [np1];
        for (j = 0; j < nv; j++)   // phenotypic effects of the allele
            all0[j] = 0;
        all0[nv] = twoN;    // the last element is the number of chromosomes in the population carrying this allele
        mutations[i].push_back(all0);
    }
    
    // traits affected by the different loci:
    
    vector<unsigned char *> traits;
    for (i = 0; i < nbSv; i++)
    {
        unsigned char * pleio = new unsigned char [mv];
        store.clear();
        for (j = 0; j < mv; j++)
        {
            do
                k = int(nv * rnd.randExc());
            while (find(store.begin(), store.end(), k) != store.end());
            store.push_back(k);
            pleio[j] = k;
        }
   
        traits.push_back(pleio);
    }
    
    
    // "pheno" will store the phenotypic distance from the optimum
	// generated by mutations present on the chromosome, for each
	// of the 2N chromosomes and along each axis (2Nn values).
	// "phenoTemp" is used during the formation of the next generation:

	double * pheno = new double [twonN];
	double * phenoTemp = new double [twonN];
	double * phenoDep = new double [twonSmp];
    
    int * par = new int [twoN];
    int * parDep = new int [twoSmp];

	// "Wtot" will hold the fitness of each individual:

	double * Wtot = new double [Nv];
	
	//Saves optimum

    double * opt = new double [nv];
    
    // phenotypic moments:
	
	double * m = new double [nv];
	double * ms = new double [nv];
	double * v = new double [nv];
	double * vs = new double [nv];
	double * m3 = new double [nv];
	double * m4 = new double [nv];
	double * m5 = new double [nv];
	double * m6 = new double [nv];
	
	double * mout = new double [nv];
	double * mself = new double [nv];
	double * vout = new double [nv];
	double * vself = new double [nv];
	
	//double * Vg = new double [nv];
	//double * Vgintra = new double [nv];

	// for time length measure:

	time_t debut, fin;
	struct tm *ptr;
	debut = time(0);
	
	// initialization: there is a single type
	// of chromosome (containing only 0's), which is the first entry of the "Chrm"
	// table. The number of copies of this genome segment is 2N:

	for (i = 0; i < fourN; i++)
    {
        Chrm[i].sel = new unsigned char [nbSv];
        Chrm[i].nbchr = 0;
    }
    for (i = 0; i < nbSv; i++)
        Chrm[0].sel[i] = 0;
    Chrm[0].nbchr = twoN;
    
    chr chrTemp;
    chrTemp.sel = new unsigned char [nbSv];
    chrTemp.nbchr = 0;

	// all 2N pointers (one for each chromosome in the population)
	// point towards Chrm[0]:
	
	for (i = 0; i < twoN; i++)
		pop[i] = &Chrm[0];

	//Creation of the optimum value
        for (i = 0; i < nv; i++){
            opt[i] = 0;
        cout<<opt[i]<<" ";}
        cout<<endl;

	// the phenotypic distance from the optimum is 0 along each phenotypic axis:
	
	for (i = 0; i < twonN; i++)
		pheno[i] = 0;
    
    int pas_c = 1;
    
    //Storing fitness values
        Wt=0;
        Wtb=0;
 
int step = 0;
gen = 0;
Ut = Uv;

 do { gen++;
    // generations:
        
        // removes alleles that are not present anymore in the population from the "mutations" table, and measures nb of alleles and diversity at selected loci:
        
        if (gen % pas_c == 0)
        {
            le = 0; divsTot = 0;
            
            for (i = 0; i < nbSv; i++)  // for each locus
            {
                for (j = 0; j < mutations[i].size(); j++)
                    mutations[i][j][nv] = 0;
                
                for (j = 0; j < twoN; j++)
                    mutations[i][(*(pop[j])).sel[i]][nv] += 1.0;
                
                divs = 0;
                for (j = mutations[i].size() - 1; j >= 0; j--)
                {
                    if (mutations[i][j][nv] == 0)
                    {
                        delete [] mutations[i][j];
                        mutations[i].erase(mutations[i].begin() + j);
                        for (k = 0; k < twoN; k++)
                            if ((*(pop[k])).sel[i] > j)
                                (*(pop[k])).sel[i] -= 1;
                    }
                    else
                        divs += (double(mutations[i][j][nv]) / twoN) * (double(mutations[i][j][nv]) / twoN);
                }
                
                divsTot += 1 - divs;
                le += mutations[i].size();  // number of alleles present in the population at locus i
            }
            le /= nbSv; divsTot /= nbSv;
        }
        
        
        // fitnesses:
        
        lnwbar = 0;
        lnvarw = 0;
        wbar = 0;
        wmax = 0;
        varw = 0;
        cmpt = 0;
        for (i = 0; i < Nv; i++)
        {
            nb = twon * i;
            d = 0;
            for (j = 0; j < nv; j++)
            {
                x = pheno[nb + j] + pheno[nb + nv + j] + gasdev() * sqrn - opt[j];
                d += x * x; // "d" is the square of the distance to the optimum
            }

            // fitness of individual i:
            w = exp(-pow(d, hQ) / twoom2);
            Wtot[i] = w;

            wbar += w;
            varw += w * w;
            if (w > 0)
            {
                lnwbar += log(w);
                lnvarw += log(w) * log(w);
                cmpt++;
            }
            if (wmax < w)
                wmax = w;
        }
        wbar /= Nv;
        varw /= Nv;
        lnwbar /= cmpt;
        lnvarw /= cmpt;
        Wt+=wbar;
			
        // average trait values:
			
        for (i = 0; i < nv; i++)
        {
            m[i] = 0;
            ms[i] = 0;
        }
        for (i = 0; i < Nv; i++)
        {
            nb = twon * i;
            for (j = 0; j < nv; j++)
            {
                m[j] += (pheno[nb + j] + pheno[nb + nv + j]);
                ms[j] += (pheno[nb + j] + pheno[nb + nv + j]) * Wtot[i] / wbar;
            }
        }
        for (i = 0; i < nv; i++)
        {
            m[i] /= Nv;
            ms[i] /= Nv;
        }
			
        // variances:
			
        for (i = 0; i < nv; i++)
        {
            v[i] = 0;
            vs[i] = 0;
        }
        for (i = 0; i < Nv; i++)
        {
            nb = twon * i;
            for (j = 0; j < nv; j++)
            {
                v[j] += pow(pheno[nb + j] + pheno[nb + nv + j] - m[j], 2);
                vs[j] += pow(pheno[nb + j] + pheno[nb + nv + j] - ms[j], 2) * Wtot[i] / wbar;
            }
        }
        for (i = 0; i < nv; i++)
        {
            v[i] /= Nv;
            vs[i] /= Nv;
        }
			
        // third moments:
			
        for (i = 0; i < nv; i++)
            m3[i] = 0;
        for (i = 0; i < Nv; i++)
        {
            nb = twon * i;
            for (j = 0; j < nv; j++)
                m3[j] += pow(pheno[nb + j] + pheno[nb + nv + j] - m[j], 3);
        }
        for (i = 0; i < nv; i++)
            m3[i] /= Nv;
			
        // fourth moments:
			
        for (i = 0; i < nv; i++)
            m4[i] = 0;
        for (i = 0; i < Nv; i++)
        {
            nb = twon * i;
            for (j = 0; j < nv; j++)
                m4[j] += pow(pheno[nb + j] + pheno[nb + nv + j] - m[j], 4);
        }
        for (i = 0; i < nv; i++)
            m4[i] /= Nv;
			
        // fifth moments:
			
        for (i = 0; i < nv; i++)
            m5[i] = 0;
        for (i = 0; i < Nv; i++)
        {
            nb = twon * i;
            for (j = 0; j < nv; j++)
                m5[j] += pow(pheno[nb + j] + pheno[nb + nv + j] - m[j], 5);
        }
        for (i = 0; i < nv; i++)
            m5[i] /= Nv;
			
        // sixth moments:
			
        for (i = 0; i < nv; i++)
            m6[i] = 0;
        for (i = 0; i < Nv; i++)
        {
            nb = twon * i;
            for (j = 0; j < nv; j++)
                m6[j] += pow(pheno[nb + j] + pheno[nb + nv + j] - m[j], 6);
        }
        for (i = 0; i < nv; i++)
            m6[i] /= Nv;

		// sampling parents for the next generation:
        
		for (j = 0; j < Nv; j++)
		{
			nb = 2*j;
			
			do{
				chr1 = int(rnd.randExc() * twoN);
				prem1 = chr1 / 2;
                
			} while (rnd.rand() > Wtot[prem1] / wmax);
            
			par[nb] = 2 * prem1;
			
			rd = rnd.rand();
			
			// selfing:
			
			if (rd < sv)
				par[nb + 1] = 2 * prem1;
			
			// outcrossing:
			
			else
			{
				do{
					chr1 = int(rnd.randExc() * twoN);
					prem1 = chr1 / 2;
					
				} while (rnd.rand() > Wtot[prem1] / wmax);
				
				par[nb + 1] = 2 * prem1;
			}
		}
        
        free = 0;
        
        for (nb = 0; nb < twoN; nb++)
		{
			nb2 = nv * nb;
			
			nbCo = int(poisdev(Lv));
            mut = int(poisdev(Ut)); // number of mutations at selected loci
			muM = rnd.rand();

			// if no mutation or cross-over:

			if ((mut == 0) && (muM > mutMv) && (nbCo == 0))
			{
				rd = rnd.rand();
                if (rd < 0.5)
					parent = par[nb];
				else
					parent = par[nb] + 1;
				temp[nb] = pop[parent];
				nb3 = nv * parent;
				for (j = 0; j < nv; j++)
					phenoTemp[nb2 + j] = pheno[nb3 + j];
			}

			else
			{
                // looking for a free entry in "Chrm" table:

				while (Chrm[free].nbchr > 0)
					free++;

				// if no cross-over:

				if ((nbCo == 0) || (pop[par[nb]] == pop[par[nb] + 1]))
				{
					rd = rnd.rand();
                    if (rd < 0.5)
						parent = par[nb];
					else
						parent = par[nb] + 1;
                    for (i = 0; i < nbSv; i++)
                        Chrm[free].sel[i] = (*(pop[parent])).sel[i];
					Chrm[free].mod = (*(pop[parent])).mod;
					nb3 = nv * parent;
					for (j = 0; j < nv; j++)
						phenoTemp[nb2 + j] = pheno[nb3 + j];
				}

				// if cross-over:

				else
				{
                    rd = rnd.rand();
                    if (rd < 0.5)
						rec(Chrm[free], *(pop[par[nb]]), *(pop[par[nb] + 1]), nbCo, nbSv);
					else
						rec(Chrm[free], *(pop[par[nb] + 1]), *(pop[par[nb]]), nbCo, nbSv);
					
					for (j = 0; j < nv; j++)
					{
						d = 0;
						for (k = 0; k < nbSv; k++)
							d += mutations[k][Chrm[free].sel[k]][j];
						phenoTemp[nb2 + j] = d;
					}
				}

				// if mutation at neutral locus:

				if (muM < mutMv)
					Chrm[free].mod = rnd.rand();

				// mutations at selected loci:

				for (nm = 0; nm < mut; nm++)
				{
					site = int(nbSv * rnd.randExc());
                    
                    if (mutations[site].size() < uchmax1)
                    {
                        for (j = 0; j < mv; j++)
                            phenoTemp[nb2 + traits[site][j]] -= mutations[site][Chrm[free].sel[site]][traits[site][j]];
                    
                        float * mutant = new float [np1];
                        for (i = 0; i < np1; i++)
                            mutant[i] = 0;
                        for (j = 0; j < mv; j++)
                            mutant[traits[site][j]] = mutations[site][Chrm[free].sel[site]][traits[site][j]] + sigs * gasdev();
                        mutations[site].push_back(mutant);
                        Chrm[free].sel[site] = mutations[site].size() - 1;
                    
                        for (j = 0; j < mv; j++)
                            phenoTemp[nb2 + traits[site][j]] += mutations[site][Chrm[free].sel[site]][traits[site][j]];
                    }
				}

				temp[nb] = &Chrm[free];
				free++;
			}
		}
		
	if((step == 0)&&(gen % eqv == 0)){
        
        double blah = fabs(1 - (Wt/double(eqv))/(Wtb/double(eqv)));
        cout<<"checking "<<Wt<<" "<<Wtb<<" "<<blah<<endl;
        

        // Equation de fitness
        if (fabs(1 - (Wt/double(eqv))/(Wtb/double(eqv))) <= diffv){
             cout<<"EQ "<<endl;

            // Changement de l'optimum de certains traits

            // Mettre un marqueur pour dire que passage a rep à selection
            fout<<"opt "<<gen<<" ";

            if (nbopt == nv){
                for (i = 0; i < nv; i++)
                    opt[i] = Nnewoptv;
            }
            else{
                store.clear();

                for (i = 0; i < nbopt; i++){

                    do
                        k = int(nv * rnd.randExc());
                    while (find(store.begin(), store.end(), k) != store.end());
                    store.push_back(k);

                    opt[k] = Nnewoptv;
                    fout<<k<<" ";
                }
            }
        for (i = 0; i < nv; i++){
        cout<<opt[i]<<" ";}
        cout<<endl;
        
       
            // On prepare la deuxieme phase de simulation
            step = 1;
            gen=0;
            
            char nomFichier2[256];
            stringstream nomF2;
            nomF2 << "eq_N" << Nv << "_s" << sv << "_n" << nv << "_m" << mv << "_sig" << sigv << "_om2_" << om2v << "_Q" << Qv
            << "_U" << Uv << "_nbS" << nbSv << "_L" << Lv << "_u" << mutMv << "_nb"<<nbr<< "_gen"<< gen << "_rep" << rep<<".txt";
            nomF2 >> nomFichier2;
            ofstream fout2;
            fout2.open(nomFichier2);
            
            for (j = 0; j < twoN; j++)
            {
                for (i = 0; i < nbSv; i++){
                
                    fout2<< mutations[i][(*(pop[j])).sel[i]][0]<<" ";
                
                }
                fout2<<pheno[j]<<endl;
            }
            
            if((step == 1)&&(optionmut == 1)){
            Ut=0;
                }
        }
        else{
            Wtb = Wt;
            Wt = 0;
        }
    }	
		
	if((step == 1)&&(gen >= NbGenv)){
        step=2;
        
        char nomFichier2[256];
            stringstream nomF2;
            nomF2 << "FIN_N" << Nv << "_s" << sv << "_n" << nv << "_m" << mv << "_sig" << sigv << "_om2_" << om2v << "_Q" << Qv
            << "_U" << Uv << "_nbS" << nbSv << "_L" << Lv << "_u" << mutMv << "_nb"<<nbr<<"_"<<rep<< "_gen"<< gen << "_rep" << rep<<".txt";
            nomF2 >> nomFichier2;
            ofstream fout2;
            fout2.open(nomFichier2);
            
            
            for (j = 0; j < twoN; j++)
            {
                for (i = 0; i < nbSv; i++){
                
                    fout2<< mutations[i][(*(pop[j])).sel[i]][0]<<" ";
                
                }
                fout2<<pheno[j]<<endl;
            }
        
        break;
		}	
		
		// update population:

		for (i = 0; i < fourN; i++)
			Chrm[i].nbchr = 0;
		
		for (i = 0; i < twonN; i++)
			pheno[i] = phenoTemp[i];

		for (i = 0; i < twoN; i++)
		{
			pop[i] = temp[i];
			(*(pop[i])).nbchr++;
		}
			
		// measures diversity and writes in result file every "pasv" generations after environmental change:
		
		if ((gen % pasv == 0)&&(step == 1))
		{
			// inbreeding depression:
            
            // sampling parents for outcrossed offspring:
            
            for (i = 0; i < twoSmp; i++)
                parDep[i] = 2 * int(rnd.randExc() * Nv);
            
            // recombination:
            
            for (nb = 0; nb < twoSmp; nb++)
            {
                nb2 = nv * nb;
                nbCo = int(poisdev(Lv));
                
                // if no mutation or cross-over:
                
                if (nbCo == 0)
                {
                    rd = rnd.rand();
                    if (rd < 0.5)
                        parent = parDep[nb];
                    else
                        parent = parDep[nb] + 1;
                    nb3 = nv * parent;
                    for (j = 0; j < nv; j++)
                        phenoDep[nb2 + j] = pheno[nb3 + j];
                }
                
                else
                {
					rd = rnd.rand();
                    if (rd < 0.5)
                        rec(chrTemp, *(pop[parDep[nb]]), *(pop[parDep[nb] + 1]), nbCo, nbSv);
                    else
                        rec(chrTemp, *(pop[parDep[nb] + 1]), *(pop[parDep[nb]]), nbCo, nbSv);
                        
                    for (j = 0; j < nv; j++)
                    {
                        d = 0;
                        for (k = 0; k < nbSv; k++)
                            d += mutations[k][chrTemp.sel[k]][j];
                        phenoDep[nb2 + j] = d;
                    }
                }
            }
			
			wout = 0;
            lnwout = 0;
            cmpt = 0;
            for (i = 0; i < smpv; i++)
            {
                nb = twon * i;
                d = 0;
                for (j = 0; j < nv; j++)
                {
                    x = phenoDep[nb + j] + phenoDep[nb + nv + j] + gasdev() * sqrn -opt[j];
                    d += x * x; // "d" is the square of the distance to the optimum
                }
					
                // fitness of individual i:
                w = exp(-pow(d, hQ) / twoom2);
					
                wout += w;
                if (w > 0)
                {
                    lnwout += log(w);
                    cmpt++;
                }
            }
            wout /= smpv;
            lnwout /= cmpt;
				
            // average trait values:
				
            for (i = 0; i < nv; i++)
            {
                mout[i] = 0;
                vout[i] = 0;
            }
            for (i = 0; i < smpv; i++)
            {
                nb = twon * i;
                for (j = 0; j < nv; j++)
                    mout[j] += (phenoDep[nb + j] + phenoDep[nb + nv + j]);
            }
            for (i = 0; i < nv; i++)
                mout[i] /= smpv;
				
            // variances:
				
            for (i = 0; i < smpv; i++)
            {
                nb = twon * i;
                for (j = 0; j < nv; j++)
                    vout[j] += pow(phenoDep[nb + j] + phenoDep[nb + nv + j] - mout[j], 2);
            }
            for (i = 0; i < nv; i++)
                vout[i] /= smpv;
			
			// sampling parents for selfed offspring:
            
            for (i = 0; i < smpv; i++)
            {
                parDep[2*i] = 2 * int(rnd.randExc() * Nv);
                parDep[2*i+1] = parDep[2*i];
            }
            
            // recombination:
            
            for (nb = 0; nb < twoSmp; nb++)
            {
                nb2 = nv * nb;
                nbCo = int(poisdev(Lv));
                
                // if no mutation or cross-over:
                
                if (nbCo == 0)
                {
                    rd = rnd.rand();
                    if (rd < 0.5)
                        parent = parDep[nb];
                    else
                        parent = parDep[nb] + 1;
                    nb3 = nv * parent;
                    for (j = 0; j < nv; j++)
                        phenoDep[nb2 + j] = pheno[nb3 + j];
                }
                
                else
                {
					rd = rnd.rand();
                    if (rd < 0.5)
                        rec(chrTemp, *(pop[parDep[nb]]), *(pop[parDep[nb] + 1]), nbCo, nbSv);
                    else
                        rec(chrTemp, *(pop[parDep[nb] + 1]), *(pop[parDep[nb]]), nbCo, nbSv);
					
                    for (j = 0; j < nv; j++)
                    {
                        d = 0;
                        for (k = 0; k < nbSv; k++)
                            d += mutations[k][chrTemp.sel[k]][j];
                        phenoDep[nb2 + j] = d;
                    }
                }
            }
			
			wself = 0;
            lnwself = 0;
            cmpt = 0;
            for (i = 0; i < smpv; i++)
            {
                nb = twon * i;
                d = 0;
                for (j = 0; j < nv; j++)
                {
                    x = phenoDep[nb + j] + phenoDep[nb + nv + j] + gasdev() * sqrn -opt[j];
                    d += x * x; // "d" is the square of the distance to the optimum
                }
					
                // fitness of individual i:
                w = exp(-pow(d, hQ) / twoom2);
					
                wself += w;
                if (w > 0)
                {
                    lnwself += log(w);
                    cmpt++;
                }
            }
            wself /= smpv;
            lnwself /= cmpt;
				
            // average trait values:
				
            for (i = 0; i < nv; i++)
            {
                mself[i] = 0;
                vself[i] = 0;
            }
            for (i = 0; i < smpv; i++)
            {
                nb = twon * i;
                for (j = 0; j < nv; j++)
                    mself[j] += (phenoDep[nb + j] + phenoDep[nb + nv + j]);
            }
            for (i = 0; i < nv; i++)
                mself[i] /= smpv;
				
            // variances:
				
            for (i = 0; i < smpv; i++)
            {
                nb = twon * i;
                for (j = 0; j < nv; j++)
                    vself[j] += pow(phenoDep[nb + j] + phenoDep[nb + nv + j] - mself[j], 2);
            }
            for (i = 0; i < nv; i++)
                vself[i] /= smpv;
            
            // neutral diversity:
			
			Fr.clear();

			for (i = 0; i < fourN; i++)
				if (Chrm[i].nbchr > 0)
				{
					frs = Fr.size();
					for (j = 0; j < frs; j++)
                        if (Fr[j].all == Chrm[i].mod)
                        {
                            Fr[j].freq += Chrm[i].nbchr;
                            break;
                        }
                    if (j == frs)
                    {
                        allTemp.all = Chrm[i].mod;
                        allTemp.freq = Chrm[i].nbchr;
                        Fr.push_back(allTemp);
                    }
				}

			frs = Fr.size();
			div = 0;
			for (j = 0; j < frs; j++)
			{	
				pp = Fr[j].freq / twoN;
				div += pp * pp;
			}

			fout << 1 - div << " " << wbar << " " << varw - (wbar*wbar) << " " << lnwbar << " " << lnvarw - (lnwbar*lnwbar)
                 << " " << wself << " " << wout << " " << lnwself << " " << lnwout << " " << le << " " << divsTot;
				
            for (i = 0; i < nv; i++)
                fout << " " << m[i] << " " << v[i] << " " << m3[i] << " " << m4[i] << " "
                     << m5[i] << " " << m6[i] << " " << ms[i] << " " << vs[i]
                     << " " << mout[i] << " " << vout[i] << " " << mself[i] << " " << vself[i] ;
            fout << endl;
                   
            char nomFichier2[256];
            stringstream nomF2;
            nomF2 << "haplo_N" << Nv << "_s" << sv << "_n" << nv << "_m" << mv << "_sig" << sigv << "_om2_" << om2v << "_Q" << Qv
            << "_U" << Uv << "_nbS" << nbSv << "_L" << Lv << "_u" << mutMv << "_nb"<<nbr<<"_"<<rep<< "_gen"<<gen<< "_rep" << rep<<".txt";
            nomF2 >> nomFichier2;
            ofstream fout2;
            fout2.open(nomFichier2);
            
            for (j = 0; j < twoN; j++)
            {
                for (i = 0; i < nbSv; i++){
                
                    fout2<< mutations[i][(*(pop[j])).sel[i]][0]<<" ";
                
                }
                fout2<<pheno[j]<<endl;
            }
        
            
          /*  for (i = 0; i < nbSv; i++){
                for (j = 0; j < mv; j++)
                {
                cout<<int(traits[i][j])<<" ";
                }
            cout<<endl;
            }*/
        
        }
        
		
		/*if (gen % miniv == 0)
		{
			for (i = 0; i < nv; i++)
			{
				 Vg[i] = 0;
				 Vgintra[i] = 0;
			}
			for (k = 0; k < nbSv; k++)
			{
				nb = k * nv;
				pp = 0;
				hom = 0;
				for (i = 0; i < twoN; i++)
					if ((*(pop[i])).sel[k] == 1)
						pp = pp + 1.0;
				for (i = 0; i < Nv; i++)
					if (((*(pop[2*i])).sel[k] == 1) && ((*(pop[2*i+1])).sel[k] == 1))
						hom = hom + 1.0;
				pp /= twoN;
				hom /= Nv;
				for (i = 0; i < nv; i++)
				{
					Vg[i] += 2.0 * mutations[nb+i] * mutations[nb+i] * pp * (1 - pp);
					Vgintra[i] += 2.0 * mutations[nb+i] * mutations[nb+i] * (hom - pp*pp);
				}
			}
			
			for (i = 0; i < nv; i++)
                fout2 << Vg[i] << " " << Vgintra[i] << " ";
            fout2 << endl;
		}*/
		
	}while(step != 2); 
    
    /*for (k = 0; k < nbSv; k++)
    {
        pp = 0;
        for (i = 0; i < twoN; i++)
            if ((*(pop[i])).sel[k] == 1)
                pp = pp + 1.0;
        freqs[k] = pp / twoN;
    }
    for (k = 0; k < nbSv; k++)
        fout4 << freqs[k] << " ";*/
	
	fin = time(0);
	
    // writes in output file:
    fprintf(fichierS, "\n\nResultats dans fichier ");
    fprintf(fichierS, nomFichier);
    fprintf(fichierS, "\n");
	         
    // time length:
    int temps = int(difftime(fin, debut));
    fprintf(fichierS,
        "\n%d generations ont pris %d heure(s) %d minute(s) %d secondes\n",
        NbGenv, temps / 3600, (temps % 3600) / 60, temps % 60);
		
    // date and time:
    ptr=localtime(&fin);
    fprintf(fichierS, asctime(ptr));
    
    for (i = 0; i < fourN; i++)
        delete [] Chrm[i].sel;
    
    for (i = 0; i < nbSv; i++)
        for (j = 0; j < mutations[i].size(); j++)
            delete [] mutations[i][j];
    
    for (i = 0; i < nbSv; i++)
        delete [] traits[i];
    
	delete [] pop;
    delete [] temp;
	delete [] Chrm;
	delete [] pheno;
	delete [] phenoTemp;
	delete [] phenoDep;
	delete [] Wtot;
	delete [] par;
    delete [] parDep;
	delete [] m;
	delete [] v;
	delete [] ms;
	delete [] vs;
	delete [] mout;
	delete [] vout;
	delete [] mself;
	delete [] vself;
	delete [] m3;
	delete [] m4;
	delete [] m5;
	delete [] m6;
	//delete [] Vg;
	//delete [] Vgintra;
	}while (rep<nbrepv);
}
