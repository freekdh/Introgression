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
#include "random.h"
#include "utils.h"

enum typenames {AB,Ab,aB,ab};
enum distribution {birth_AB,birth_Ab,birth_aB,birth_ab,death_AB,death_Ab,death_aB,death_ab};        

struct BDPopulation{
    public:
    BDPopulation(const int &AB0, const int &Ab0, const int &aB0, const int ab0) {
        type[AB] = AB0;
        type[Ab] = Ab0;
        type[aB] = aB0;
        type[ab] = ab0; 
    }
    bool Iterate(const double &bA, const double &ba, const double &dA, const double &da, const double r){
        static double f[4];
        assert(da>=0.0);
        assert(dA>=0.0);
        if(CalculateFrequencies(f)==false) {return false;};
        const double sumdeathrate = (f[AB]+f[Ab])*dA+(f[aB]+f[ab])*da;
        const double sumbirthrate = (f[AB]+f[Ab])*bA+(f[aB]+f[ab])*ba;
        const double Pdeath = sumdeathrate/(sumbirthrate+sumdeathrate);
        const double Pbirth = 1-Pdeath;
        const double rD = r*(f[AB]*f[ab]-f[Ab]*f[aB]);
        const double PAB_death = f[AB]*dA/sumdeathrate;
        const double PAb_death = f[Ab]*dA/sumdeathrate;
        const double PaB_death = f[aB]*da/sumdeathrate;
        const double Pab_death = f[ab]*da/sumdeathrate;
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
 
    /*
        inline void returnTypes(int *typepointer){
            typepointer = type[];
        }
    */

    private:
    int type[4];
    bool CalculateFrequencies(double *f){
        const int sum=type[AB]+type[Ab]+type[aB]+type[ab];
        if(sum<=0.0){return false;}
        f[AB] = type[AB]/sum;
        f[Ab] = type[Ab]/sum;
        f[aB] = type[aB]/sum;
        f[ab] = type[ab]/sum;
        assert(f[AB]+f[Ab]+f[aB]+f[ab]==1.0);
        return true;
    }
};

int main(){

    BDPopulation population(0,1,9,0);
    int i = 0;
    while(population.Iterate(1.0,1.0,2.0-1.1,2.0-0.9,0.5)==true || i!= 100){
        ++i;
    }
    return 0;
}