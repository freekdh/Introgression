/*=============================================================================================================
                                                IntrogressionSimulations.cpp
===============================================================================================================

 Simulations of genetic rescue

 C++-code accompanying:
		 
		(ms. in prep).

 Written by:
        F.J.H. de Haas
       	Theoretical Biology Group
        University of British Columbia
        the Netherlands

 Program version
		xx/xx/xxxx	:

=============================================================================================================*/

#include <assert.h>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <iostream>
#include "random.h"
#include "utils.h"
#include <Rcpp.h>

enum typenames {AB,Ab,aB,ab};
enum distribution {birth_AB,birth_Ab,birth_aB,birth_ab,death_AB,death_Ab,death_aB,death_ab};        

struct Parameters{
    Parameters(int in_AB0, int in_Ab0, int in_aB0, int in_ab0, double in_bA,double in_ba, double in_dA, double in_da, double in_r){
        AB0 = in_AB0;
        Ab0 = in_Ab0;
        aB0 = in_aB0;
        ab0 = in_ab0;
        bA = in_bA;
        ba = in_ba;
        dA = in_dA;
        da = in_da;
        r = in_r;
    }

    Parameters(Rcpp::List parslist){
        AB0 = parslist["AB0"];
        Ab0 = parslist["Ab0"];
        aB0 = parslist["aB0"];
        ab0 = parslist["ab0"];
        bA = parslist["bA"];
        ba = parslist["ba"];
        dA = parslist["dA"];
        da = parslist["da"];
        r = parslist["r"];
    }
    int AB0,Ab0,aB0,ab0;
    double bA,ba,dA,da;
    double r;
};

class BDPopulation{
    public:
    BDPopulation(const Parameters &pars) {
        type[AB] = pars.AB0;
        type[Ab] = pars.Ab0;
        type[aB] = pars.aB0;
        type[ab] = pars.ab0; 
    }
    bool Iterate(const Parameters &pars){
        static double f[4];
        assert(pars.da>=0.0);
        assert(pars.dA>=0.0);
        if(CalculateFrequencies(f)==false) {return false;};
        const double sumdeathrate = (f[AB]+f[Ab])*pars.dA+(f[aB]+f[ab])*pars.da;
        const double sumbirthrate = (f[AB]+f[Ab])*pars.bA+(f[aB]+f[ab])*pars.ba;
        const double Pdeath = sumdeathrate/(sumbirthrate+sumdeathrate);
        const double Pbirth = 1-Pdeath;
        const double rD = pars.r*(f[AB]*f[ab]-f[Ab]*f[aB]);
        const double PAB_death = f[AB]*pars.dA/sumdeathrate;
        const double PAb_death = f[Ab]*pars.dA/sumdeathrate;
        const double PaB_death = f[aB]*pars.da/sumdeathrate;
        const double Pab_death = f[ab]*pars.da/sumdeathrate;
        const double PAB_birth = f[AB]-rD;
        const double PAb_birth = f[Ab]+rD;
        const double PaB_birth = f[aB]+rD;
        const double Pab_birth = f[ab]-rD;

        rnd::discrete_distribution birthdeathevent(8);
        birthdeathevent[birth_AB] = Pbirth * PAB_birth;
        birthdeathevent[birth_Ab] = Pbirth * PAb_birth;
        birthdeathevent[birth_aB] = Pbirth * PaB_birth;
        birthdeathevent[birth_ab] = Pbirth * Pab_birth;
        
        birthdeathevent[death_AB] = Pdeath * PAB_death;
        birthdeathevent[death_Ab] = Pdeath * PAb_death;
        birthdeathevent[death_aB] = Pdeath * PaB_death;
        birthdeathevent[death_ab] = Pdeath * Pab_death;

        const int sample = birthdeathevent.sample();
        switch(sample) {
            case birth_AB:      ++type[AB]; break;
            case birth_Ab:      ++type[Ab]; break;
            case birth_aB:      ++type[aB]; break;
            case birth_ab:      ++type[ab]; break;
            case death_AB:      --type[AB]; break;
            case death_Ab:      --type[Ab]; break;
            case death_aB:      --type[aB]; break;
            case death_ab:      --type[ab]; break;
        }
        return true;
    }
    inline int returnTypes(int typenr){
        return type[typenr];
    }

    private:
    int type[4];
    bool CalculateFrequencies(double *f){
        const double sum=(double)(type[AB]+type[Ab]+type[aB]+type[ab]);
        if(sum<=0){return false;}
        f[AB] = (double)type[AB]/sum;
        f[Ab] = (double)type[Ab]/sum;
        f[aB] = (double)type[aB]/sum;
        f[ab] = (double)type[ab]/sum;
        return true;
    }
};

// [[Rcpp::export]]
Rcpp::List BDSim(int tend, Rcpp::List parslist){
    rnd::set_seed();
    Parameters pars(parslist);

    BDPopulation population(pars);

    int i = 0;
    while(population.Iterate(pars)==true && i!= tend){
        ++i;
    }
    
    return Rcpp::List::create(
        Rcpp::_["AB"] = population.returnTypes(0),
        Rcpp::_["Ab"] = population.returnTypes(1),
        Rcpp::_["aB"] = population.returnTypes(2),
        Rcpp::_["ab"] = population.returnTypes(3)
    );
}