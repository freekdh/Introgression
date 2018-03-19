
#include <boost/dynamic_bitset.hpp>
#include <assert.h>
#include "random.h"
#include <fstream>
#include <stdlib.h>

boost::dynamic_bitset<> r_global;
boost::dynamic_bitset<> m_global;

class Individual{
    public:
    Individual(const int NPLOIDY, const int NLOCI) {
			assert(NPLOIDY%2==0 && NPLOIDY>0 && NLOCI > 0);
			genome.resize(NPLOIDY);
			for(int i = 0; i < NPLOIDY; ++i){
				genome[i].resize(NLOCI);
			}
    }

    Individual(Individual parent[2]){
        const int NPLOIDY = parent[0].genome.size();
        const int NLOCI = parent[0].genome[0].size();

        assert(NPLOIDY == parent[1].genome.size());
        for(int i = 0; i < NPLOIDY; ++i){
            assert(NLOCI == parent[1].genome[i].size());
        }

        // Initialize individual
        assert(NPLOIDY%2==0 && NPLOIDY>0);
        genome.resize(NPLOIDY);
        for(int i = 0; i < NPLOIDY; ++i){
            genome[i].resize(NLOCI);
        }

        // Recombination & Mutation
        boost::dynamic_bitset<> temp1, temp2, temp3, r_localbit, m_localbit;

        r_localbit.resize(NLOCI);
        m_localbit.resize(NLOCI);

        assert(r_localbit.size() == NLOCI);
        assert(m_localbit.size() == NLOCI);

        for(int j = 0; j < NPLOIDY; j+=2){
            for(int p = 0; p < 2; ++p){
                Make_localbit(r_localbit,r_global);
                Make_localbit(m_localbit, m_global);

                temp1 = parent[p].genome[j] & r_localbit;
                temp2 = parent[p].genome[j+1] & r_localbit.flip();
                temp3 = temp1 | temp2;
            
                temp1 = temp3 & m_localbit;
                temp2 = temp3 | m_localbit;
                genome[j+p] = temp2 & temp1.flip();
            }
        }
    }

    double AdditiveViability(const std::vector<double> &scgenome){
        const int NPLOIDY = genome.size();
        const int NLOCI = genome[0].size();
        assert(scgenome.size() == genome[0].size());
        double output = 0.0;
        for(int i = 0; i < NPLOIDY; ++i){
            for(int j = 0; j < NLOCI; ++j){
                if(genome[i][j] == 1) {output += scgenome[j];}
            }
        }
        output /= (double)NLOCI;
        assert(0.0 <= output <= 1.0);
        return output;
    }

    void Flipbit(const int &chromosome, const int &locus){ genome[chromosome][locus].flip();}

    bool Genotype(const int &chromosome, const int &locus) {return genome[chromosome][locus];}

    private:
    void Make_localbit(boost::dynamic_bitset<> &local, const boost::dynamic_bitset<> &global){
        const int NLOCI = genome[0].size();
        const int global_size = global.size();

        assert(local.size() == NLOCI);
        assert(global_size > NLOCI);
        int start = rnd::integer(global_size-NLOCI);

        for(int i = 0; i < NLOCI; ++i){
            local[i] = global[start+i];
        }
    }

    std::vector<boost::dynamic_bitset<>> genome;    // genome of individual
};

std::vector<Individual*>::iterator it;

void ResetRGlobal(boost::dynamic_bitset<> &global, const int &GLOBALMAX, const double &RRATE){
    global.resize(GLOBALMAX);

    boost::dynamic_bitset<> a(1);
    a[0] = rnd::bernoulli(0.5);
    for(int i = 0; i < GLOBALMAX; ++i){
        if(rnd::bernoulli(RRATE) == true) a.flip();
        global[i] = a[0];
    }

}

void ResetMGlobal(boost::dynamic_bitset<> &global, const int &GLOBALMAX, const double &MRATE){
    global.resize(GLOBALMAX);
    
    for(int i = 0; i < GLOBALMAX; ++i){
        global[i] = rnd::bernoulli(MRATE);
    }
}

void ItteratePopulation(std::vector<Individual*> &population, const double &FERTILITY, const std::vector<double> &SC_GENOME){
    const int popsize = population.size();
    assert(popsize > 0);
    std::vector<Individual*> offspring;
    offspring.resize((double)popsize*FERTILITY);
    assert(offspring.size() > 0);

    for(it = offspring.begin(); it != offspring.end(); ++it){
        Individual c[2] = {*population[rnd::integer(popsize)], *population[rnd::integer(popsize)]};
        *it = new Individual(c);
    }

    for(it = population.begin(); it != population.end(); ++it){
        delete *it;
    }
    population.clear();

    for(it = offspring.begin(); it != offspring.end(); ++it){
        if(rnd::uniform() > (*it)->AdditiveViability(SC_GENOME)) {population.push_back(*it);} ;
    }
}

void StationaryPopulation(std::vector<Individual*> &population, const double &FERTILITY, const std::vector<double> &SC_GENOME){
    const int popsize = population.size();
    assert(popsize > 0);
    std::vector<Individual*> offspring(popsize);

    rnd::discrete_distribution viabilitydist(popsize);
    for(int i = 0; i < popsize; ++i){
        viabilitydist[i] = 1.0-population[i]->AdditiveViability(SC_GENOME);
    }
    
    for(it = offspring.begin(); it != offspring.end(); ++it){
        Individual c[2] = {*population[viabilitydist.sample()], *population[viabilitydist.sample()]};
        *it = new Individual(c);
    }

    for(it = population.begin(); it != population.end(); ++it){
        delete *it;
    }
    population.clear();

    population = offspring;
}

void CollectData(std::vector<Individual*> &population, std::ofstream &output, const int NPLOIDY, const int NLOCI){
    output.fill(',');
    std::vector<int> temp(NLOCI,0);
    for(it = population.begin(); it != population.end(); ++it){
        for(int i = 0; i < NPLOIDY; ++i){
            for(int j = 0; j < NLOCI; ++j){
                temp[j] += (int)(*it)->Genotype(i,j);
            }
        }
    }
    static int gen = 0;
    output << gen << output.fill() << population.size() << output.fill();
    for(int i = 0; i < NLOCI; ++i){
        output << (double)temp[i]/((double)NPLOIDY*(double)population.size()) << output.fill();
    }

    ++gen;
    
    output << std::endl;
}

int main(int argc, char *argv[]){

    rnd::set_seed();

    // Parameters
    const int GLOBALMAX = 100000;
    const double RRATE = 0.5;
    const double MRATE = 0.0001;
    const int FERTILITY = 1.1;

    const int NLOCI = atoi(argv[1]);        assert(NLOCI > 0);
    const int NPLOIDY = atoi(argv[2]);      assert(NPLOIDY%2 == 0 && NPLOIDY > 0);
    const int NPOPSIZE = atoi(argv[3]);     assert(NPOPSIZE > 0);
    const int NGEN = atoi(argv[4]);         assert(NGEN > 0);

    std::vector<double> SC_GENOME(NLOCI);
    for(int i = 0; i < NLOCI; ++i){
        SC_GENOME[i] = 0.0;
    }

    // Initialize global variables
    ResetRGlobal(r_global, GLOBALMAX, RRATE);
    ResetMGlobal(m_global, GLOBALMAX, MRATE);

    // Create population in mutation-selection-drift balance
    std::vector<Individual*> population(NPOPSIZE);
    for(it = population.begin(); it != population.end(); ++it){
        *it = new Individual(NPLOIDY,NLOCI);
    }

    std::ofstream output("output.csv");
    output.fill(',');
    output << "gen" << output.fill() << "popsize" << output.fill();
    for(int i = 0; i < NLOCI; ++i){
        output << "locus" << i << output.fill();
    }
    output << std::endl;

    for(int i = 0; i < NGEN; ++i){
        StationaryPopulation(population,FERTILITY,SC_GENOME);
        CollectData(population,output,NPLOIDY,NLOCI);
    }

    // Environmental change
    SC_GENOME[0] = 0.1;

    for(int i = 0; i < 100; ++i){
        ItteratePopulation(population,FERTILITY,SC_GENOME);
        CollectData(population,output,NPLOIDY,NLOCI);
    }

    return 0;
}