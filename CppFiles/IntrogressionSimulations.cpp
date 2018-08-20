#include <boost/dynamic_bitset.hpp>
#include <assert.h>
#include <stdlib.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <string>
#include "random.h"
#include "IntrogressionSimulations.h"

class Individual{
  public:
    Individual(const boost::dynamic_bitset<> &INITGENOME) : genome(INITGENOME) {;}

    Individual(Individual* parent1, Individual* parent2, const Parameters &pars){
        Individual* parent[2] = {parent1, parent2};
        const int NLOCI = parent[0]->genome[0].size();
        genome.resize(NLOCI);

        // Recombination & Mutation
        boost::dynamic_bitset<> temp1, temp2, temp3, r_localbit, m_localbit;
        r_localbit.resize(NLOCI);
        m_localbit.resize(NLOCI);

        Make_localbit(r_localbit, pars.r_global);
        Make_localbit(m_localbit, pars.m_global);

        temp1 = parent[0]->genome & r_localbit;
        temp2 = parent[1]->genome & r_localbit.flip();
        temp3 = temp1 | temp2;

        temp1 = temp3 & m_localbit;
        temp2 = temp3 | m_localbit;
        genome = temp2 & temp1.flip();
    }

    inline bool rescue(const Parameters &pars){ return genome[pars.index[0]] ? true : false;} 

    bool Genotype(const int &locus) {return genome[locus]; }

    int GenotypeCount() {return genome.count();}
    
    private:
    void Make_localbit(boost::dynamic_bitset<> &local, const boost::dynamic_bitset<> &global)
    {
        const int NLOCI = genome[0].size();
        const int global_size = global.size();

        assert(local.size() == NLOCI);
        assert(global_size > NLOCI);
        int start = rnd::integer(global_size - NLOCI);

        for (int i = 0; i < NLOCI; ++i)
        {
            local[i] = global[start + i];
        }
    };

    boost::dynamic_bitset<> genome; // haploids
};

void WriteToDataBlock(std::vector<Individual*> &population, const Parameters &pars, DataBlock* &SimData){
    // Popsize
    const int populationsize = population.size();
    SimData->popsize.push_back(populationsize);
    std::vector<int> count[2];
    count[0].resize(pars.NLOCI,0);
    count[1].resize(pars.NLOCI,0);
    int temp0 = 0, temp1 = 0, ind0 = 0, ind1 = 0;
    for(Individual* ind : population){
        if(ind->Genotype(pars.index[0])){
            ++ind1;
            temp1 += ind->GenotypeCount() - 1;
        }
        else{
            ++ind0;
            temp0 += ind->GenotypeCount(); 
        } 
        // Genome allele frequencies
        for(int j = 0; j < pars.NLOCI; ++j){
            ++count[ind->Genotype(j)][j];
        }
    }

    SimData->allele0.push_back(count[0]);
    SimData->allele1.push_back(count[1]);
    SimData->major0.push_back(ind0);
    SimData->major1.push_back(ind1);
    SimData->introgressed0.push_back(temp0/(double)((pars.NLOCI-1)*ind0));
    SimData->introgressed1.push_back(temp1/(double)((pars.NLOCI-1)*ind1));
}

bool ItteratePopulation(std::vector<Individual*> &population, const Parameters &pars){ 
    std::vector<Individual*>::iterator it;

    // Birth
    const int nparents = population.size();
    const double popgrowthrate = 1.0 + pars.BIRTHRATE * (1.0 - ((double)nparents / (double)pars.K));
    assert(popgrowthrate >= 0.0);
    const int noffspring = rnd::poisson((double)nparents * popgrowthrate);
    std::vector<Individual*> offspring(noffspring);

    // Random mating 
    for(it = offspring.begin(); it != offspring.end(); ++it){
        *it = new Individual(population[population[rnd::integer(nparents)]],population[rnd::integer(nparents)],pars);
    }

    // Death (selection)
    for(it = offspring.begin(); it != offspring.end(); ++it){
        it->Genotype(0,)
    }

    for(Individual* ind : offspring)
        ind->Genotype(0,pars.index[0]);

    // Cleanup
    for(Individual* ind : population) delete ind;

    population.clear();
    population = offspring;

    // Calculate if rescue?
    if(rescue == false){
            if((*it)->rescue(pars)==true) {rescue = true;}
        }

    return rescue;
}

bool RunSimulation(const Parameters &SimPars, SimData &SimulationData){

    // Set local parameters
    DataBlock* SimData = new DataBlock;

    // Initialize population
    std::vector<Individual*> population(SimPars.NINIT[0]+SimPars.NINIT[1]);
    for(int i = 0; i < SimPars.NINIT[0]; ++i){
        population[i] = new Individual(SimPars.INIT_GENOME[0]);
    }
    for(int i = SimPars.NINIT[0]; i < SimPars.NINIT[0]+SimPars.NINIT[1];++i){
        population[i] = new Individual(SimPars.INIT_GENOME[1]);
    }

    // Run simulation
    WriteToDataBlock(population, SimPars, SimData);

    for (int i = 0; i < SimPars.NGEN; ++i){
        if(ItteratePopulation(population, SimPars)==false){
            for (Individual* i: population) delete i;
            delete SimData;
            ++SimulationData.nofixcounter;
            return false;
            }
        else{
            WriteToDataBlock(population, SimPars, SimData);
            }
    }    
  
    // Cleanup and Write DataBlock
    for (Individual* i: population) delete i;
    SimulationData.push_back_protect(SimData);
    return true;
}

#ifdef SHINYFUNCTION_H
Parameters::Parameters(const Rcpp::List &parslist){
    MUTATIONRATE = 0.0;
    BIRTHRATE = parslist["b"];
    DEATHRATEA = parslist["dA"];
    DEATHRATEa = parslist["da"];
    NLOCI = parslist["nloci"];
    NINIT[0] = parslist["ninit0"];
    NINIT[1] = parslist["ninit1"];
    NGEN = parslist["ngen"];
    NREP = parslist["nrep"];
    RECOMBINATIONRATE = parslist["rec"];
    K = parslist["k"];

    Initialize();
}
#endif

Parameters::Parameters(int argc, char *argv[]){
        MUTATIONRATE = 0.0;
        BIRTHRATE = std::atof(argv[1]);
        DEATHRATEA = std::atof(argv[2]);
        DEATHRATEa = std::atof(argv[3]);
        NLOCI = std::atoi(argv[4]);
        NPLOIDY = std::atoi(argv[5]);
        NINIT[0] = std::atoi(argv[6]);
        NINIT[1] = std::atoi(argv[7]);
        NGEN = std::atoi(argv[8]);
        NREP = std::atoi(argv[9]);
        RECOMBINATIONRATE = std::atof(argv[10]);
        K = std::atoi(argv[11]);

        Initialize();
}

void Parameters::Initialize(){
    boost::dynamic_bitset<>> INIT_GENOME0(NLOCI));
    boost::dynamic_bitset<> INIT_GENOME1(NLOCI).set();
    INIT_GENOME[0] = INIT_GENOME0;
    INIT_GENOME[1] = INIT_GENOME1;

    //Genetic architecture of selection
    index.clear();
    index.resize(NLOCI);
    index[0] = std::floor((double)NLOCI/2.0);

    //Recombination and mutation
    const int GLOBALMAX = 100000;
    r_global.resize(GLOBALMAX);
    boost::dynamic_bitset<> a(1);
    a[0] = rnd::bernoulli(0.5);
    for (int i = 0; i < GLOBALMAX; ++i)
    {
        if (rnd::bernoulli(RECOMBINATIONRATE) == true)
        a.flip();
        r_global[i] = a[0];
    }
    m_global.resize(GLOBALMAX);
    for (int i = 0; i < GLOBALMAX; ++i)
    {
        m_global[i] = rnd::bernoulli(MUTATIONRATE);
    }
}

/*
double AdditiveFitness(const Parameters &pars){
    double viability = 1.0;
    for(int i = 0; i < pars.NPLOIDY; ++i){
        for(int j = 0; j < pars.NLOCAL_ADAPTED_LOCI+1; ++j){
            if(genome[i][pars.index[j]] == true){viability += pars.SC_GENOME[pars.index[j]];}
        }
    }
    return viability;
}
*/