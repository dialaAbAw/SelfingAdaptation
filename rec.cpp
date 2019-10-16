#include "fisher.h"
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

extern MTRand rnd;

// rec function: generates recombinant genome segment "res"
// from parental genome segments "c1" and "c2"
// nbCo is the number of cross-overs in the genome segment,
// nS the number of selected loci

void rec(chr &res, chr &c1, chr &c2, int nbCo, int nS)
{
    vector<int> pos;
    int i, j;
    double rd;
    bool test = 0;
    int hnS = nS / 2;
    if (nbCo == 0)
    {
        res.mod = c1.mod;
        for (i = 0; i < nS; i++)
            res.sel[i] = c1.sel[i];
    }
    else
    {
        int locus = 0;
        
        for (j = 0; j < nbCo; j++)
            pos.push_back(rnd.randInt(nS)); // positions of cross-overs
        sort(pos.begin(), pos.end());
        
        
        for (i = 0; i < pos.size(); i++)
        {
            if (i % 2 == 0)
                for (j = locus; j < pos[i]; j++)
                    res.sel[j] = c1.sel[j];
            else
                for (j = locus; j < pos[i]; j++)
                    res.sel[j] = c2.sel[j];
            locus = pos[i];
            
            if ((test == 0) && (pos[i] >= hnS))
            {
                test = 1;
                if (pos[i] == hnS)
                {
                    rd = rnd.rand();
                    if (rd < 0.5)
                        res.mod = c1.mod;
                    else
                        res.mod = c2.mod;
                }
                else
                {
                    if (i % 2 == 0)
                        res.mod = c2.mod;
                    else
                        res.mod = c1.mod;
                }
            }
        }
        
        if (pos.size() % 2 == 0)
            for (j = locus; j < nS; j++)
                res.sel[j] = c1.sel[j];
        else
            for (j = locus; j < nS; j++)
                res.sel[j] = c2.sel[j];
        
        if (test == 0)
        {
            if (nbCo % 2 == 0)
                res.mod = c1.mod;
            else
                res.mod = c2.mod;
        }
    }
}
