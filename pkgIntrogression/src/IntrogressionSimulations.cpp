#include "IntrogressionSimulations.h"

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

Parameters::Initialize(){
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