#include <boost/dynamic_bitset.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <assert.h>
#include "random.h"
#include <stdlib.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <string>
#include <atomic>
#include <mutex>
#include <Rcpp.h>
#include <progress.hpp>
#ifdef _OPENMP
    #include <omp.h>
#endif

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(BH)]]

using namespace boost::accumulators;

std::mutex mu_datablock;

struct Parameters{
    public:
    Parameters(double r, int nloci, int nploidy, int ninit0, int ninit1, int distlocal, double scmajor, double sclocal, int ngen, int nrep, double rec, int k, int threads = 0){
        MUTATIONRATE = 0.0;
        INTRINSIC_GROWTHRATE = r;
        NLOCI = nloci;
        NPLOIDY = nploidy;
        NINIT[0] = ninit0;
        NINIT[1] = ninit1;
        NLOCAL_ADAPTED_LOCI = 1;
        DISTLOCAL = distlocal;
        SC_MAJOR = scmajor;
        SC_LOCAL = sclocal;
        NGEN = ngen;
        NREP = nrep;
        RECOMBINATIONRATE = rec;
        K = k;

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
    
    // Can be public bc no problem with multithread access.
    double RECOMBINATIONRATE;
    double MUTATIONRATE = 0.0;
    double INTRINSIC_GROWTHRATE;
    int NGEN;
    int NLOCI;
    int DISTLOCAL;
    int NPLOIDY;
    int NREP;
    int NINIT[2];
    int K;
    double SC_MAJOR;
    double SC_LOCAL;
    int NLOCAL_ADAPTED_LOCI;

    std::vector<double> SC_GENOME;
    std::vector<int> index;
    std::vector<boost::dynamic_bitset<>> INIT_GENOME[2];
    boost::dynamic_bitset<> r_global;
    boost::dynamic_bitset<> m_global;
};

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

struct DataBlock{
    public:
    std::vector<int> popsize;
    std::vector<std::vector<int>> allele0;
    std::vector<std::vector<int>> allele1;
    std::vector<int> major0;
    std::vector<int> major1;
    std::vector<double> introgressed0;
    std::vector<double> introgressed1;
};

// Public struct for multithreading
struct SimData{  
    SimData(){
        nofixcounter = 0;
    }

    void push_back_protect(DataBlock* &datablock){
        mu_datablock.lock();
        DataSet.push_back(datablock);
        mu_datablock.unlock();
    };
 
    std::vector<DataBlock*> DataSet;
    std::atomic<int> nofixcounter; // for thread safety incrementation
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

Rcpp::List WriteOutput(const Parameters &GlobalPars, SimData &SimulationData){

   // Parameters
    Rcpp::DataFrame parsdata =  Rcpp::DataFrame::create(
        Rcpp::_["NINIT0"] = GlobalPars.NINIT[0], 
        Rcpp::_["NINIT1"] = GlobalPars.NINIT[1],
        Rcpp::_["NGEN"] = GlobalPars.NGEN,
        Rcpp::_["NREP"] = GlobalPars.NREP,
        Rcpp::_["NLOCI"] = GlobalPars.NLOCI,
        Rcpp::_["NPLOIDY"] = GlobalPars.NPLOIDY,
        Rcpp::_["MUTATIONRATE"] = GlobalPars.MUTATIONRATE,
        Rcpp::_["RECOMBINATIONRATE"] = GlobalPars.RECOMBINATIONRATE,
        Rcpp::_["GROWTHRATE"] = GlobalPars.INTRINSIC_GROWTHRATE,
        Rcpp::_["CARRYINGCAPACITY"] = GlobalPars.K,
        Rcpp::_["SC_MAJOR"] = GlobalPars.SC_MAJOR,
        Rcpp::_["SC_LOCAL"] = GlobalPars.SC_LOCAL,
        Rcpp::_["INDEX_MAJOR"] = GlobalPars.index[0],
        Rcpp::_["INDEX_LOCAL"] = GlobalPars.index[1],
        Rcpp::_["DISTLOCAL"] = GlobalPars.DISTLOCAL
    );
    
    // Write outputfiles
    Rcpp::NumericVector generation(GlobalPars.NGEN);    
    Rcpp::NumericVector popsizevecm(GlobalPars.NGEN);
    Rcpp::NumericVector major0vecm(GlobalPars.NGEN);
    Rcpp::NumericVector major1vecm(GlobalPars.NGEN);
    Rcpp::NumericVector introgressed0vecm(GlobalPars.NGEN);
    Rcpp::NumericVector introgressed1vecm(GlobalPars.NGEN);
    Rcpp::NumericVector popsizevecv(GlobalPars.NGEN);
    Rcpp::NumericVector major0vecv(GlobalPars.NGEN);
    Rcpp::NumericVector major1vecv(GlobalPars.NGEN);
    Rcpp::NumericVector introgressed0vecv(GlobalPars.NGEN);
    Rcpp::NumericVector introgressed1vecv(GlobalPars.NGEN);
    
    for(int i = 0; i < GlobalPars.NGEN; ++i){
        accumulator_set<int, stats<tag::mean, tag::variance > > popsize;
        accumulator_set<double, stats<tag::mean, tag::variance > > major0;
        accumulator_set<double, stats<tag::mean, tag::variance > > major1;
        accumulator_set<double, stats<tag::mean, tag::variance > > introgressed0;
        accumulator_set<double, stats<tag::mean, tag::variance > > introgressed1;

        for(int j = 0; j < GlobalPars.NREP; ++j){
            popsize(SimulationData.DataSet[j]->popsize[i]);
            major0((double)SimulationData.DataSet[j]->major0[i]/(double)GlobalPars.NPLOIDY);
            major1((double)SimulationData.DataSet[j]->major1[i]/(double)GlobalPars.NPLOIDY);
            introgressed0(SimulationData.DataSet[j]->introgressed0[i]);
            introgressed1(SimulationData.DataSet[j]->introgressed1[i]);
        }
        generation[i] = i;
        popsizevecm[i] = mean(popsize);
        major0vecm[i] = mean(major0);
        major1vecm[i] = mean(major1);
        introgressed0vecm[i] = mean(introgressed0);
        introgressed1vecm[i] = mean(introgressed1);
        popsizevecv[i] = variance(popsize);
        major0vecv[i] = variance(major0);
        major1vecv[i] = variance(major1);
        introgressed0vecv[i] = variance(introgressed0);
        introgressed1vecv[i] = variance(introgressed1);
    };
    
    Rcpp::DataFrame data =  Rcpp::DataFrame::create(
        Rcpp::_["Generation"] = generation,
        Rcpp::_["Popsize_avg"] = (popsizevecm), 
        Rcpp::_["Major0_avg"] = (major0vecm), 
        Rcpp::_["Major1_avg"] = (major1vecm), 
        Rcpp::_["Introgressed0_avg"] = (introgressed0vecm), 
        Rcpp::_["Introgressed1_avg"] = (introgressed1vecm),
        Rcpp::_["Popsize_var"] = (popsizevecv), 
        Rcpp::_["Major0_var"] = (major0vecv), 
        Rcpp::_["Major1_var"] = (major1vecv), 
        Rcpp::_["Introgressed0_var"] = (introgressed0vecv), 
        Rcpp::_["Introgressed1_var"] = (introgressed1vecv)
    );

    std::vector<Rcpp::NumericVector> allelefrequencymean(GlobalPars.NLOCI);
    std::vector<Rcpp::NumericVector> allelefrequencyvar(GlobalPars.NLOCI);

    for(int i = 0; i < GlobalPars.NLOCI; ++i){
        allelefrequencymean[i] = Rcpp::NumericVector(GlobalPars.NGEN); 
        allelefrequencyvar[i] = Rcpp::NumericVector(GlobalPars.NGEN); 
    }

    for(int i = 0; i < GlobalPars.NGEN; ++i){
        for(int l = 0; l < GlobalPars.NLOCI; ++l)
        {
            accumulator_set<double, stats<tag::mean, tag::variance > > locus;
            for(int r = 0; r < GlobalPars.NREP; ++r)
            {
                locus((double)SimulationData.DataSet[r]->allele0[i][l] / ((double)SimulationData.DataSet[r]->popsize[i] * (double)GlobalPars.NPLOIDY));
            }
            allelefrequencymean[l][i] = mean(locus);
            allelefrequencyvar[l][i] = variance(locus);
        }
    }   

    Rcpp::DataFrame alleledatamean, alleledatavar;
    std::string avglocus = "avglocus", varlocus = "varlocus";
    alleledatamean.push_back(generation,"Generation");
    alleledatavar.push_back(generation,"Generation");
    for(int i = 0; i < GlobalPars.NLOCI; ++i){
        alleledatamean.push_back((allelefrequencymean[i]), avglocus+std::to_string(i));
        alleledatavar.push_back((allelefrequencyvar[i]), varlocus+std::to_string(i));
    }

    double fixation = (double)GlobalPars.NREP/(double(GlobalPars.NREP+SimulationData.nofixcounter.load()));
        
    // Cleanup
    for(int i = 0; i < GlobalPars.NREP; ++i){
        delete SimulationData.DataSet[i];
    }

    return Rcpp::List::create(
        Rcpp::_["pars"] = (parsdata),
        Rcpp::_["data"] = (data),
        Rcpp::_["allelefavg"] = (alleledatamean),
        Rcpp::_["allelefvar"] = (alleledatavar),
        Rcpp::_["fixation"] = fixation
    );

}

// [[Rcpp::export]]
Rcpp::List IntrogressionSimulation(double r, int nloci, int nploidy, int ninit0, int ninit1, int distlocal, double scmajor, double sclocal, int ngen, int nrep, double rec, int k, int threads = 0){
    // Prepare for simulation
    rnd::set_seed();
    const Parameters GlobalPars(r, nloci, nploidy, ninit0, ninit1, distlocal, scmajor, sclocal, ngen, nrep, rec, k);
    SimData SimulationData;

    // Run nrep successful simulations
    #ifdef _OPENMP
        if(threads>0) omp_set_num_threads(threads);
        REprintf("Parallel activated : Number of threads=%i\n",omp_get_max_threads());   
    #endif
    Progress p(nrep, true);

    #pragma omp parallel for schedule(static)
    for (int task = 0; task < nrep; ++task){
        while(RunSimulation(GlobalPars, SimulationData)==false);
        p.increment();
    }

    return WriteOutput(GlobalPars, SimulationData);
}
