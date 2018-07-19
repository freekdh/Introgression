#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <assert.h>
#include "random.h"
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#ifdef _OPENMP
    #include <omp.h>
#endif

using namespace boost::accumulators;
enum nameofstream {Parametersoff, Dataoff, AlleleFrequencyoff_mean, AlleleFrequencyoff_var};

boost::dynamic_bitset<> r_global;
boost::dynamic_bitset<> m_global;

struct Parameters{
    double RECOMBINATIONRATE = 0.5;
    double MUTATIONRATE = 0.0;
    double INTRINSIC_GROWTHRATE = 0.01;
    int NGEN = 20;
    int NLOCI = 20;
    int DISTLOCAL = 1;
    int NPLOIDY = 2;
    int NREP = 1;
    int NINIT[2] = {1,10};
    int K = 100;
    std::vector<double> SC_GENOME;
    std::vector<int> index;    
    std::vector<int> v;
    int threads = 0;

    std::vector<boost::dynamic_bitset<>> INIT_GENOME[2];

    double SC_MAJOR;
    double SC_LOCAL;
    int NLOCAL_ADAPTED_LOCI;
};

class Individual{
  public:
    Individual(const std::vector<boost::dynamic_bitset<>> &INITGENOME) : genome(INITGENOME) {;}

    Individual(Individual* parent1, Individual* parent2){
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
                Make_localbit(r_localbit, r_global);
                Make_localbit(m_localbit, m_global);

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
    std::vector<int> introgressed0;
    std::vector<int> introgressed1;
};

std::vector<DataBlock*> DataSet; 

void WriteToDataBlock(std::vector<Individual*> &population, const Parameters &pars, DataBlock* &SimData){

    // Popsize
    SimData->popsize.push_back(population.size());
    std::vector<int> count[2];
    count[0].resize(pars.NLOCI,0);
    count[1].resize(pars.NLOCI,0);
    int temp0 = 0, temp1 = 0, ind0 = 0, ind1 = 0;
    for(Individual* ind : population){
        for(int i = 0; i < pars.NPLOIDY; ++i){
            if(ind->Genotype(i,pars.index[0])){
                ++ind0;
                temp0 += ind->GenotypeCount(i) - 1;
            }
            else{
                ++ind1;
                temp1 += ind->GenotypeCount(i);
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
    SimData->introgressed0.push_back(temp0);
    SimData->introgressed1.push_back(temp1);
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
        *it = new Individual(population[fitnessdist.sample()],population[fitnessdist.sample()]);
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

bool RunSimulation(const Parameters &SimPars, const int &tasknr){
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
            return false;
            }
        else{
            WriteToDataBlock(population, SimPars, SimData);
            }
    }    
  
    // Cleanup and Write DataBlock
    for (Individual* i: population) delete i;
    DataSet[tasknr] = SimData;
    return true;
}

void InitializeGlobalEnvironment(const Parameters &GlobalPars){
    rnd::set_seed();
    //srand(static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    const int GLOBALMAX = 100000;

    r_global.resize(GLOBALMAX);
    boost::dynamic_bitset<> a(1);
    a[0] = rnd::bernoulli(0.5);
    for (int i = 0; i < GLOBALMAX; ++i)
    {
        if (rnd::bernoulli(GlobalPars.RECOMBINATIONRATE) == true)
        a.flip();
        r_global[i] = a[0];
    }

    m_global.resize(GLOBALMAX);
    for (int i = 0; i < GLOBALMAX; ++i)
    {
        m_global[i] = rnd::bernoulli(GlobalPars.MUTATIONRATE);
    }
}

void CollectParameters(int argc, char* argv[], Parameters &GlobalPars){
    if (argc != 13) {std::cout << "Argc != 13" << std::endl;}
    
    //General 
    GlobalPars.MUTATIONRATE = 0.0;
    GlobalPars.INTRINSIC_GROWTHRATE = 0.1;
    GlobalPars.NLOCI = atoi(argv[1]);
    GlobalPars.NPLOIDY = atoi(argv[2]);
    GlobalPars.NINIT[0] = atoi(argv[3]);
    GlobalPars.NINIT[1] = atoi(argv[4]);
    GlobalPars.NLOCAL_ADAPTED_LOCI = 1;
    GlobalPars.DISTLOCAL = atoi(argv[5]);
    GlobalPars.SC_MAJOR = atof(argv[6]);
    GlobalPars.SC_LOCAL = atof(argv[7]);
    GlobalPars.NGEN = atoi(argv[8]);
    GlobalPars.NREP = atoi(argv[9]);
    GlobalPars.RECOMBINATIONRATE = atof(argv[10]);
    GlobalPars.K = atoi(argv[11]);
    GlobalPars.threads = atoi(argv[12]);
    GlobalPars.SC_GENOME.clear();
    GlobalPars.SC_GENOME.resize(GlobalPars.NLOCI, 0.0);
    std::vector<boost::dynamic_bitset<>> INIT_GENOME0(GlobalPars.NPLOIDY, boost::dynamic_bitset<>(GlobalPars.NLOCI));
    std::vector<boost::dynamic_bitset<>> INIT_GENOME1(GlobalPars.NPLOIDY, boost::dynamic_bitset<>(GlobalPars.NLOCI).set());
    GlobalPars.INIT_GENOME[0] = INIT_GENOME0;
    GlobalPars.INIT_GENOME[1] = INIT_GENOME1;

    //Genetic architecture of selection
    GlobalPars.index.clear();
    GlobalPars.index.resize(GlobalPars.NLOCI);
    GlobalPars.index[0] = std::floor((double)GlobalPars.NLOCI/2.0);
    assert((GlobalPars.DISTLOCAL+GlobalPars.index[0]) < GlobalPars.NLOCI);
    GlobalPars.index[1] = GlobalPars.index[0]+GlobalPars.DISTLOCAL;

    GlobalPars.SC_GENOME.clear();
    GlobalPars.SC_GENOME.resize(GlobalPars.NLOCI, 0.0);
    GlobalPars.SC_GENOME[GlobalPars.index[0]] = GlobalPars.SC_MAJOR;
    GlobalPars.SC_GENOME[GlobalPars.index[1]] = GlobalPars.SC_LOCAL;

    // Make sure all parameters are there
    assert(GlobalPars.NREP > 0);
    assert(GlobalPars.NLOCI > 1);
    assert(GlobalPars.NLOCI >= GlobalPars.NLOCAL_ADAPTED_LOCI + 1);// nlocaladaptedloci + major locus
    assert(GlobalPars.NPLOIDY % 2 == 0 && GlobalPars.NPLOIDY > 0);
    assert(GlobalPars.NINIT[0] > 0);
    assert(GlobalPars.NINIT[1] > 0);
    assert(GlobalPars.NGEN > 0);
    assert(1.0 >= GlobalPars.SC_LOCAL >= 0.0);
    assert(1.0 >= GlobalPars.SC_MAJOR >= 0.0);
    assert(GlobalPars.NLOCAL_ADAPTED_LOCI >= 0);
    assert(GlobalPars.NLOCAL_ADAPTED_LOCI + 1 < GlobalPars.NLOCI);
}

void WriteOutput(std::ofstream arrayofstream[4], Parameters &pars){
    // Parameters output
    arrayofstream[Parametersoff] << "SC_MAJOR" << arrayofstream[Parametersoff].fill() << pars.SC_MAJOR << std::endl;
    arrayofstream[Parametersoff] << "SC_LOCAL" << arrayofstream[Parametersoff].fill() << pars.SC_LOCAL << std::endl;
    arrayofstream[Parametersoff] << "NLOCALLOCI" << arrayofstream[Parametersoff].fill() << pars.NLOCAL_ADAPTED_LOCI << std::endl;
    for (int i = 0; i < 2; ++i)
    {
        arrayofstream[Parametersoff] << "NINIT" << i << arrayofstream[Parametersoff].fill() << pars.NINIT[i] << std::endl;
    }
    arrayofstream[Parametersoff] << "NGEN" << arrayofstream[Parametersoff].fill() << pars.NGEN << std::endl;
    arrayofstream[Parametersoff] << "NLOCI" << arrayofstream[Parametersoff].fill() << pars.NLOCI << std::endl;
    arrayofstream[Parametersoff] << "NPLOIDY" << arrayofstream[Parametersoff].fill() << pars.NPLOIDY << std::endl;
    arrayofstream[Parametersoff] << "NREP" << arrayofstream[Parametersoff].fill() << pars.NREP << std::endl;
    arrayofstream[Parametersoff] << "MUTATIONRATE" << arrayofstream[Parametersoff].fill() << pars.MUTATIONRATE << std::endl;
    arrayofstream[Parametersoff] << "RECOMBINATIONRATE" << arrayofstream[Parametersoff].fill() << pars.RECOMBINATIONRATE << std::endl;
    arrayofstream[Parametersoff] << "INTRINSICGROWTHRATE" << arrayofstream[Parametersoff].fill() << pars.INTRINSIC_GROWTHRATE << std::endl;
    arrayofstream[Parametersoff] << "CARRYINGCAPACITY" << arrayofstream[Parametersoff].fill() << pars.K << std::endl;
    arrayofstream[Parametersoff] << "DISTANCEFROMMAJOR" << arrayofstream[Parametersoff].fill() << pars.DISTLOCAL << std::endl;
    arrayofstream[Parametersoff] << "SC_GENOME" << std::endl;

    for (int j = 0; j < pars.NLOCI; ++j) {
            arrayofstream[Parametersoff] << pars.SC_GENOME[j] << " ";
        }
    arrayofstream[Parametersoff] << std::endl;

    //Simulation output
    arrayofstream[Dataoff] 
    << "Generation" << arrayofstream[Dataoff].fill() 
    << "AVG_popsize" << arrayofstream[Dataoff].fill()
    << "AVG_major0" << arrayofstream[Dataoff].fill() 
    << "AVG_major1" << arrayofstream[Dataoff].fill() 
    << "AVG_introgressed0" << arrayofstream[Dataoff].fill()
    << "AVG_introgressed1" << arrayofstream[Dataoff].fill() 
    << "VAR_popsize" << arrayofstream[Dataoff].fill()
    << "VAR_major0" << arrayofstream[Dataoff].fill() 
    << "VAR_major1" << arrayofstream[Dataoff].fill()
    << "VAR_introgressed0" << arrayofstream[Dataoff].fill() 
    << "VAR_introgressed1" << std::endl;

    for(int i = 0; i < pars.NGEN; ++i){
        accumulator_set<int, stats<tag::mean, tag::variance > > popsize;
        accumulator_set<int, stats<tag::mean, tag::variance > > major0;
        accumulator_set<int, stats<tag::mean, tag::variance > > major1;
        accumulator_set<int, stats<tag::mean, tag::variance > > introgressed0;
        accumulator_set<int, stats<tag::mean, tag::variance > > introgressed1;

        ///FINISH ANALYSIS OF THE GENOTYPE RESCUE PER LOCUS PART. 
        for(int j = 0; j < pars.NREP; ++j){
            popsize(DataSet[j]->popsize[i]);
            major0(DataSet[j]->major0[i]);
            major1(DataSet[j]->major1[i]);
            introgressed0(DataSet[j]->introgressed0[i]);
            introgressed1(DataSet[j]->introgressed1[i]);
        }

        arrayofstream[Dataoff]
        << i << arrayofstream[Dataoff].fill() 
        << mean(popsize) << arrayofstream[Dataoff].fill()
        << (double)mean(major0)/(double)pars.NPLOIDY << arrayofstream[Dataoff].fill()
        << (double)mean(major1)/(double)pars.NPLOIDY << arrayofstream[Dataoff].fill() 
        << mean(introgressed0)/(double)(pars.NPLOIDY*(pars.NLOCI-1)) << arrayofstream[Dataoff].fill()
        << mean(introgressed1)/(double)(pars.NPLOIDY*(pars.NLOCI-1)) << arrayofstream[Dataoff].fill()

        << variance(popsize) << arrayofstream[Dataoff].fill()
        << (double)variance(major0)/((double)pars.NPLOIDY*(double)pars.NPLOIDY) << arrayofstream[Dataoff].fill()
        << (double)variance(major1)/((double)pars.NPLOIDY*(double)pars.NPLOIDY) << arrayofstream[Dataoff].fill()
        << variance(introgressed0)/(double)(pars.NPLOIDY*(pars.NLOCI-1)*pars.NPLOIDY*(pars.NLOCI-1)) << arrayofstream[Dataoff].fill()
        << variance(introgressed1)/(double)(pars.NPLOIDY*(pars.NLOCI-1)*pars.NPLOIDY*(pars.NLOCI-1)) << arrayofstream[Dataoff].fill() << std::endl;
    };
    
    for(int i = 0; i < pars.NGEN; ++i){
        arrayofstream[AlleleFrequencyoff_mean] << i << arrayofstream[AlleleFrequencyoff_mean].fill(); 
        arrayofstream[AlleleFrequencyoff_var] << i << arrayofstream[AlleleFrequencyoff_var].fill();
        for(int l = 0; l < pars.NLOCI; ++l)
        {
            accumulator_set<double, stats<tag::mean, tag::variance > > locus;
            for(int r = 0; r < pars.NREP; ++r)
            {
                locus((double)DataSet[r]->allele0[i][l] / ((double)DataSet[r]->popsize[i] * (double)pars.NPLOIDY));
            }
            arrayofstream[AlleleFrequencyoff_mean] << mean(locus)  << arrayofstream[AlleleFrequencyoff_mean].fill();
            arrayofstream[AlleleFrequencyoff_var] << variance(locus) << arrayofstream[AlleleFrequencyoff_var].fill(); 
        }
        arrayofstream[AlleleFrequencyoff_mean] << std::endl;
        arrayofstream[AlleleFrequencyoff_var] << std::endl;
    }
}

std::string CreateOutputStreams(std::ofstream arrayofstream[4]){
    
    std::string CurrentWorkingDirectory = boost::filesystem::current_path().c_str();
    CurrentWorkingDirectory.append("/Data");

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
    std::string str = oss.str();

    std::string MainOutputFolder = CurrentWorkingDirectory;
    MainOutputFolder.append("/");
    MainOutputFolder.append(str);

    boost::filesystem::create_directories(MainOutputFolder.c_str());

    std::string ParameterOutput = MainOutputFolder;
    ParameterOutput.append("/Parameters.csv");

    arrayofstream[Parametersoff].open(ParameterOutput);
    arrayofstream[Parametersoff].fill(',');

    std::string DataOutput = MainOutputFolder;
    DataOutput.append("/Data.csv");

    std::string AlleleFOutput_mean = MainOutputFolder;
    std::string AlleleFOutput_var = MainOutputFolder;

    AlleleFOutput_mean.append("/AlleleF_mean.csv");
    AlleleFOutput_var.append("/AlleleF_var.csv");

    arrayofstream[Dataoff].open(DataOutput);
    arrayofstream[AlleleFrequencyoff_mean].open(AlleleFOutput_mean);
    arrayofstream[AlleleFrequencyoff_var].open(AlleleFOutput_var);
    arrayofstream[Dataoff].fill(',');
    arrayofstream[AlleleFrequencyoff_mean].fill(',');
    arrayofstream[AlleleFrequencyoff_var].fill(',');

    return MainOutputFolder;
}

int main(int argc, char *argv[]){

    // Initialize simulation
    Parameters GlobalPars;
    CollectParameters(argc,argv,GlobalPars);
    InitializeGlobalEnvironment(GlobalPars);

    // Run nrep successful simulations
    auto start = std::chrono::high_resolution_clock::now();
    #ifdef _OPENMP
        if(GlobalPars.threads > 0)
            omp_set_num_threads(GlobalPars.threads);
        std::cout << "Parallel activated : Number of threads= " << omp_get_max_threads() << std::endl;   
    #endif
    #pragma omp parallel for
    for (int task = 0; task < GlobalPars.NREP; ++task){
        while(RunSimulation(GlobalPars, task)==false); 
    }
    auto finish = std::chrono::high_resolution_clock::now();
   
    // Write outputfiles
    std::ofstream arrayofstream[4]; 
    CreateOutputStreams(arrayofstream); 
	WriteOutput(arrayofstream, GlobalPars);

    // Cleanup
    for(DataBlock* i : DataSet) delete i;
    DataSet.clear();

    // Show Time elapsed
    std::chrono::duration<double> elapsed = finish-start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
    return 0;
}
