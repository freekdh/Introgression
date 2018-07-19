#ifndef GLOBALSTRUCT_H
#define GLOBALSTRUCT_H

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
#include <mutex>
#include <atomic>
#include <vector>

std::mutex mu_datablock;

struct Parameters{
    Parameters(){};
    #ifdef RCPPFUNCTION_H
        Parameters(const Rcpp::List &);
    #endif
    #ifdef CSVFUNCTION_H
        Parameters(double r, int nloci, int nploidy, int ninit0, int ninit1, int distlocal, double scmajor, double sclocal, int ngen, int nrep, double rec, int k, int threads);    
    #endif
    
    double RECOMBINATIONRATE;
    double MUTATIONRATE = 0.0   ;
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
    boost::dynamic_bitset<> r_global, m_global;
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

struct SimData{  
    SimData(){nofixcounter = 0;}
    void push_back_protect(DataBlock* &datablock){
        mu_datablock.lock();
        DataSet.push_back(datablock);
        mu_datablock.unlock();
    };
 
    std::vector<DataBlock*> DataSet;
    std::atomic<int> nofixcounter;
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

void WriteToDataBlock(std::vector<Individual*> &population, const Parameters &pars, DataBlock* &SimData){};

bool ItteratePopulation(std::vector<Individual*> &population, const Parameters &pars){};

bool RunSimulation(const Parameters &SimPars, SimData &SimulationData){};

#endif
