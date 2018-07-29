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
    Individual(const std::vector<boost::dynamic_bitset<>> &INITGENOME) : genome(INITGENOME) {;}

    Individual(Individual* parent1, Individual* parent2, const Parameters &pars){
        Individual* parent[2] = {parent1, parent2};
        const int NPLOIDY = parent[0]->genome.size();
        const int NLOCI = parent[0]->genome[0].size();

        assert(NPLOIDY == parent[1]->genome.size());
        for (int i = 0; i < NPLOIDY; ++i){
            assert(NLOCI == parent[1]->genome[i].size());
        }

        // Initialize individual
        assert(NPLOIDY % 2 == 0 && NPLOIDY > 0);
        genome.resize(NPLOIDY);
        for (int i = 0; i < NPLOIDY; ++i)
        {
            genome[i].resize(NLOCI);
        }

        // Recombination & Mutation
        boost::dynamic_bitset<> temp1, temp2, temp3, r_localbit, m_localbit;

        r_localbit.resize(NLOCI);
        m_localbit.resize(NLOCI);

        assert(r_localbit.size() == NLOCI);
        assert(m_localbit.size() == NLOCI);

        for (int j = 0; j < NPLOIDY; j += 2)
        {
            for (int p = 0; p < 2; ++p)
            {
                Make_localbit(r_localbit, pars.r_global);
                Make_localbit(m_localbit, pars.m_global);

                temp1 = parent[p]->genome[j] & r_localbit;
                temp2 = parent[p]->genome[j + 1] & r_localbit.flip();
                temp3 = temp1 | temp2;

                temp1 = temp3 & m_localbit;
                temp2 = temp3 | m_localbit;
                genome[j + p] = temp2 & temp1.flip();
            }
        }
    }

    double AdditiveFitness(const Parameters &pars){
        double viability = 1.0;
        for(int i = 0; i < pars.NPLOIDY; ++i){
            for(int j = 0; j < pars.NLOCAL_ADAPTED_LOCI+1; ++j){
                if(genome[i][pars.index[j]] == true){viability += pars.SC_GENOME[pars.index[j]];}
            }
        }
        return viability;
    }

    bool rescue(const Parameters &pars){
        for(int i = 0; i < pars.NPLOIDY; ++i)
            if(genome[i][pars.index[0]]) return true;
        
        return false;
    }

    bool Genotype(const int &chromosome, const int &locus) { return genome[chromosome][locus]; }

    int GenotypeCount(const int &chromosome) {return genome[chromosome].count();}
    
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

    std::vector<boost::dynamic_bitset<>> genome; // circular genome of individual
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
        for(int i = 0; i < pars.NPLOIDY; ++i){
            if(ind->Genotype(i,pars.index[0])){
                ++ind1;
                temp1 += ind->GenotypeCount(i) - 1;
            }
            else{
                ++ind0;
                temp0 += ind->GenotypeCount(i); 
            } 
            // Genome allele frequencies
            for(int j = 0; j < pars.NLOCI; ++j){
                ++count[ind->Genotype(i,j)][j];
            }
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

    // Calculate offspring population size
    const int nparents = population.size();
    const double popgrowthrate = 1.0 + pars.INTRINSIC_GROWTHRATE * (1.0 - ((double)nparents / (double)pars.K));
    assert(popgrowthrate >= 0.0);
    const int noffspring = rnd::poisson((double)nparents * popgrowthrate);
    std::vector<Individual*> offspring(noffspring);

    // Selection 
    bool rescue = false;
    rnd::discrete_distribution fitnessdist((int)nparents);
    for(int i = 0; i < nparents; ++i)
        fitnessdist[i] = population[i]->AdditiveFitness(pars);
    for(it = offspring.begin(); it != offspring.end(); ++it){
        *it = new Individual(population[fitnessdist.sample()],population[fitnessdist.sample()],pars);
        if(rescue == false){
            if((*it)->rescue(pars)==true) {rescue = true;}
        }
    }

    for(Individual* ind : offspring)
        ind->Genotype(0,pars.index[0]);

    // Cleanup
    for(Individual* ind : population) delete ind;

    population.clear();
    population = offspring;

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

Parameters::Parameters(const Rcpp::List &parslist){
    MUTATIONRATE = 0.0;
    INTRINSIC_GROWTHRATE = parslist["r"];
    NLOCI = parslist["nloci"];
    NPLOIDY = parslist["nploidy"];
    NINIT[0] = parslist["ninit0"];
    NINIT[1] = parslist["ninit1"];
    NLOCAL_ADAPTED_LOCI = 1;
    DISTLOCAL = parslist["distlocal"];
    SC_MAJOR = parslist["scmajor"];
    SC_LOCAL = parslist["sclocal"];
    NGEN = parslist["ngen"];
    NREP = parslist["nrep"];
    RECOMBINATIONRATE = parslist["rec"];
    K = parslist["k"];

    Initialize();
}

void Parameters::Initialize(){
    std::vector<boost::dynamic_bitset<>> INIT_GENOME0(NPLOIDY, boost::dynamic_bitset<>(NLOCI));
    std::vector<boost::dynamic_bitset<>> INIT_GENOME1(NPLOIDY, boost::dynamic_bitset<>(NLOCI).set());
    INIT_GENOME[0] = INIT_GENOME0;
    INIT_GENOME[1] = INIT_GENOME1;

    SC_GENOME.clear();
    SC_GENOME.resize(NLOCI, 0.0);

    //Genetic architecture of selection
    index.clear();
    index.resize(NLOCI);
    index[0] = std::floor((double)NLOCI/2.0);
    assert((DISTLOCAL+index[0]) < NLOCI);
    index[1] = index[0]+DISTLOCAL;

    SC_GENOME.clear();
    SC_GENOME.resize(NLOCI, 0.0);
    SC_GENOME[index[0]] = SC_MAJOR;
    SC_GENOME[index[1]] = SC_LOCAL;

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