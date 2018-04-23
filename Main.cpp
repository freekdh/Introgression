
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

using namespace boost::accumulators;

boost::dynamic_bitset<> r_global;
boost::dynamic_bitset<> m_global;

struct Parameters
{
    //default parameters
    double RECOMBINATIONRATE;
    double MUTATIONRATE;
    double INTRINSIC_GROWTHRATE;
    int NGEN;
    int NLOCI;
    int NPLOIDY;
    int NREP;
    int NINIT[2];
    int K;
    std::vector<double> SC_GENOME;
    std::vector<int> index;    
    std::vector<int> v;

    std::vector<boost::dynamic_bitset<>> INIT_GENOME[2];

    double SC_MAJOR;
    double SC_LOCAL;
    int NLOCAL_ADAPTED_LOCI;
};

class Individual
{
  public:
    Individual(const std::vector<boost::dynamic_bitset<>> &INITGENOME) : genome(INITGENOME) { ; }

    Individual(Individual parent[2])
    {
        const int NPLOIDY = parent[0].genome.size();
        const int NLOCI = parent[0].genome[0].size();

        assert(NPLOIDY == parent[1].genome.size());
        for (int i = 0; i < NPLOIDY; ++i)
        {
            assert(NLOCI == parent[1].genome[i].size());
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

                temp1 = parent[p].genome[j] & r_localbit;
                temp2 = parent[p].genome[j + 1] & r_localbit.flip();
                temp3 = temp1 | temp2;

                temp1 = temp3 & m_localbit;
                temp2 = temp3 | m_localbit;
                genome[j + p] = temp2 & temp1.flip();
            }
        }
    }

    double MultiplicativeViability(const Parameters &pars)
    {
        // Sally? what is going on here. Don't understand this multiplicative viability thing.
        double viability = 0.0;
        for(int i = 0; i < pars.NPLOIDY; ++i){
            for(int j = 0; j < pars.NLOCAL_ADAPTED_LOCI+1; ++j){
                if(genome[i][pars.index[j]] == true){viability *= pars.SC_GENOME[pars.index[j]];}
            }
        }

        viability = 1.0-viability;

        assert(0.0 <= viability <= 1.0);
        return viability;
    }

    void Flipbit(const int &chromosome, const int &locus) { genome[chromosome][locus].flip(); }

    bool Genotype(const int &chromosome, const int &locus) { return genome[chromosome][locus]; }

    int Neutralloci(const int &chromosome, const Parameters &pars){
        // v[0] = location of major locus
        // v[1:Nlocal] = location of local loci
        // v[Nlocal+1 : v.end()] = location of neutral loci
        int counter = 0;
        counter += genome[chromosome].count();
        for(int i = 0; i < pars.NLOCAL_ADAPTED_LOCI + 1; ++i){
            counter -= genome[chromosome][pars.index[i]];
        }
        assert(counter>=0);
        return counter;
    }
    
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

std::vector<Individual*>::iterator it;

struct DataBlock{
    public:

    std::vector<double> SC_GENOME;
    std::vector<double> p;
    std::vector<double> q;
    std::vector<int> popsize;
    std::vector<int> rescue;
};

std::vector<DataBlock*> DataSet;    // Store replicates

void CalculateMCMC(std::vector<Individual*> &population, const Parameters &pars, DataBlock* &SimData){
    int n[2][2];
    n[0][0] = 0;
    n[1][0] = 0;
    n[0][1] = 0;
    n[1][1] = 0;
    int nmajor = 0;
    for(it = population.begin(); it != population.end(); ++it)
        for(int i = 0; i < pars.NPLOIDY; ++i){
            nmajor += (*it)->Genotype(i,pars.index[0]);
            for(int j = 0; j < pars.NLOCI-1; ++j)
                ++n[(*it)->Genotype(i,j)][(*it)->Genotype(i,j+1)];
        }
            
    double pmax = double(n[0][1]) / double(n[0][1]+n[0][0]);
    double qmax = double(n[1][0]) / double(n[1][0]+n[1][1]);

    SimData->p.push_back(pmax);
    SimData->q.push_back(qmax);
    SimData->popsize.push_back(population.size());
}

void ResetRGlobal(boost::dynamic_bitset<> &global, const int &GLOBALMAX, const double &RRATE)
{
    global.resize(GLOBALMAX);

    boost::dynamic_bitset<> a(1);
    a[0] = rnd::bernoulli(0.5);
    for (int i = 0; i < GLOBALMAX; ++i)
    {
        if (rnd::bernoulli(RRATE) == true)
            a.flip();
        global[i] = a[0];
    }
}

void ResetMGlobal(boost::dynamic_bitset<> &global, const int &GLOBALMAX, const double &MRATE)
{
    global.resize(GLOBALMAX);

    for (int i = 0; i < GLOBALMAX; ++i)
    {
        global[i] = rnd::bernoulli(MRATE);
    }
}

void ItteratePopulation(std::vector<Individual*> &population, const Parameters &pars)
{
    // Density dependent birth rate (fertility)
    const int popsize = population.size();
    std::vector<Individual*> offspring(rnd::poisson(
        (double)popsize * (1.0 + pars.INTRINSIC_GROWTHRATE * (1.0 - ((double)popsize / (double)pars.K)))
        ));

    // Mating + Create offspring
    for (it = offspring.begin(); it != offspring.end(); ++it)
    {
        Individual c[2] = {*population[rnd::integer(popsize)], *population[rnd::integer(popsize)]};
        *it = new Individual(c);
    }

    // Parents die
    for (it = population.begin(); it != population.end(); ++it)
    {
        delete *it;
    }
    population.clear();

    // Viability selection on offspring
    for (it = offspring.begin(); it != offspring.end(); ++it)
    {
        if (rnd::uniform() < (*it)->MultiplicativeViability(pars))
        {
            population.push_back(*it);
        }
        else {delete *it;} 
    }

}

std::string ReturnTimeStamp(const std::string &CurrentDirectory)
{
    std::string out;

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
    std::string str = oss.str();

    out = CurrentDirectory;
    out.append("/");
    out.append(str);

    return out;
}

void OutputParameters(std::ofstream &ofstream, const Parameters &pars)
{

    ofstream.fill(',');

    ofstream << "SC_MAJOR" << ofstream.fill() << pars.SC_MAJOR << std::endl;
    ofstream << "SC_LOCAL" << ofstream.fill() << pars.SC_LOCAL << std::endl;
    ofstream << "NLOCALLOCI" << ofstream.fill() << pars.NLOCAL_ADAPTED_LOCI << std::endl;

    for (int i = 0; i < 2; ++i)
    {
        ofstream << "NINIT" << i << ofstream.fill() << pars.NINIT[i] << std::endl;
    }
    ofstream << "NGEN" << ofstream.fill() << pars.NGEN << std::endl;
    ofstream << "NLOCI" << ofstream.fill() << pars.NLOCI << std::endl;
    ofstream << "NPLOIDY" << ofstream.fill() << pars.NPLOIDY << std::endl;
    ofstream << "NREP" << ofstream.fill() << pars.NREP << std::endl;
    ofstream << "MUTATIONRATE" << ofstream.fill() << pars.MUTATIONRATE << std::endl;
    ofstream << "RECOMBINATIONRATE" << ofstream.fill() << pars.RECOMBINATIONRATE << std::endl;
    ofstream << "INTRINSICGROWTHRATE" << ofstream.fill() << pars.INTRINSIC_GROWTHRATE << std::endl;
    ofstream << "CARRYINGCAPACITY" << ofstream.fill() << pars.K << std::endl;

    ofstream << "SC_GENOME" << std::endl;
    for (int i = 0; i < DataSet.size(); ++i){
        for (int j = 0; j < pars.NLOCI; ++j) {
            ofstream << DataSet[i]->SC_GENOME[j] << ofstream.fill();
        }
        ofstream << std::endl;
    }
    ofstream << std::endl;
}

std::string CreateOutputStreams(std::ofstream &ParametersOfstream, std::ofstream &DataOfstream)
{
    std::string CurrentWorkingDirectory = "/home/freek/Introgression";
    //std::string Current = boost::filesystem::current_path();
    std::string MainOutputFolder = ReturnTimeStamp(CurrentWorkingDirectory);
    boost::filesystem::create_directories(MainOutputFolder.c_str());

    std::string ParameterOutput = MainOutputFolder;
    ParameterOutput.append("/Parameters.csv");

    ParametersOfstream.open(ParameterOutput);
    ParametersOfstream.fill(',');

    std::string DataOutput = MainOutputFolder;
    DataOutput.append("/Data.csv");

    DataOfstream.open(DataOutput);
    DataOfstream.fill(',');
    
    return MainOutputFolder;
}

bool RescuePopulation(std::vector<Individual*> population, const Parameters &pars){
    for(it = population.begin(); it != population.end(); ++it){
        for(int i = 0; i < pars.NPLOIDY; ++i){
            if((*it)->Genotype(i,pars.index[0])) return true;
        }
    }
    
    return false;
}

void RunSimulation(const Parameters &GlobalPars)
{
    // Set local parameters
    Parameters SimPars = GlobalPars;
    DataBlock* SimData = new DataBlock;

    SimPars.index.clear();
    SimPars.index.resize(SimPars.NLOCI);
    std::iota(std::begin(SimPars.index), std::end(SimPars.index), 0);
    std::random_shuffle(SimPars.index.begin(), SimPars.index.end());

    SimPars.SC_GENOME.clear();
    SimPars.SC_GENOME.resize(SimPars.NLOCI, 0.0);
    SimPars.SC_GENOME[SimPars.index[0]] = SimPars.SC_MAJOR;          // Major allele is index = 0;
    for (int i = 0; i < SimPars.NLOCAL_ADAPTED_LOCI; ++i)
    {
        SimPars.SC_GENOME[SimPars.index[i + 1]] = SimPars.SC_LOCAL;  // Local adapted alleles index 1 - nlocaladaptedloci;
    }

    SimData->SC_GENOME = SimPars.SC_GENOME;
    
    // Initialize population
    std::vector<Individual*> population(SimPars.NINIT[0] + SimPars.NINIT[1]);
    for (it = population.begin(); it != population.begin() + SimPars.NINIT[0]; ++it)
    {
        *it = new Individual(SimPars.INIT_GENOME[0]);
    }
    for (it = population.begin() + SimPars.NINIT[0]; it != population.end(); ++it)
    {
        *it = new Individual(SimPars.INIT_GENOME[1]);
    }

    // Itterate population
    for (int i = 0; i < SimPars.NGEN; ++i)
    {
        CalculateMCMC(population, SimPars, SimData);
        ItteratePopulation(population, SimPars);
    }    

    for (it = population.begin(); it != population.end(); ++it)
    {
        delete *it;
    }

    DataSet.push_back(SimData);
}

void AssertParameters(Parameters &pars){
    // Make sure all parameters are there
    assert(pars.NREP > 0);
    assert(pars.NLOCI > 1);
    assert(pars.NLOCI >= pars.NLOCAL_ADAPTED_LOCI + 1);// nlocaladaptedloci + major locus
    assert(pars.NPLOIDY % 2 == 0 && pars.NPLOIDY > 0);
    assert(pars.NINIT[0] > 0);
    assert(pars.NINIT[1] > 0);
    assert(pars.NGEN > 0);
    assert(1.0 >= pars.SC_LOCAL >= 0.0);
    assert(1.0 >= pars.SC_MAJOR >= 0.0);
    assert(pars.NLOCAL_ADAPTED_LOCI >= 0);
    assert(pars.NLOCAL_ADAPTED_LOCI + 1 < pars.NLOCI);
}

void WriteOutput(std::ofstream &output, Parameters &pars){

    output 
    << "Generation" << output.fill() 
    << "AVG_p" << output.fill() 
    << "VAR_p" << output.fill() 
    << "AVG_q" << output.fill() 
    << "VAR_q" << output.fill() 
    << "AVG_size" << output.fill()
    << "VAR_size" << std::endl; 
    for(int i = 0; i < pars.NGEN; ++i){

        // Analysis
        accumulator_set<double, stats<tag::mean, tag::variance > > pGlobal;
        accumulator_set<double, stats<tag::mean, tag::variance > > qGlobal;
        accumulator_set<int, stats<tag::mean, tag::variance > > popsize;

        for(int j = 0; j < pars.NREP; ++j){
            pGlobal(DataSet[j]->p[i]);
            qGlobal(DataSet[j]->q[i]);
            popsize(DataSet[j]->popsize[i]);
        }

        output
    
        << i << output.fill() << mean(pGlobal) << output.fill()
        << variance(pGlobal) << output.fill()
        << mean(qGlobal) << output.fill() 
        << variance(qGlobal) << output.fill()
        << mean(popsize) << output.fill()
        << variance(popsize) << std::endl;
        
    }
}

int main(int argc, char *argv[])
{
    rnd::set_seed();
    srand(static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));

    // Parameters
    if (argc != 12) {std::cout << "Argc != 12" << std::endl; return -1;}
    const int GLOBALMAX = 100000;
    Parameters GlobalPars;
    {
        GlobalPars.MUTATIONRATE = 0.0;
        GlobalPars.INTRINSIC_GROWTHRATE = 0.1;

        GlobalPars.NLOCI = atoi(argv[1]);
        GlobalPars.NPLOIDY = atoi(argv[2]);
        GlobalPars.NINIT[0] = atoi(argv[3]);
        GlobalPars.NINIT[1] = atoi(argv[4]);
        GlobalPars.NLOCAL_ADAPTED_LOCI = atoi(argv[5]);
        GlobalPars.SC_MAJOR = atof(argv[6]);
        GlobalPars.SC_LOCAL = atof(argv[7]);
        GlobalPars.NGEN = atoi(argv[8]);
        GlobalPars.NREP = atoi(argv[9]);
        GlobalPars.RECOMBINATIONRATE = atof(argv[10]);
        GlobalPars.K = atoi(argv[11]);

        GlobalPars.SC_GENOME.clear();
        GlobalPars.SC_GENOME.resize(GlobalPars.NLOCI, 0.0);
        std::vector<boost::dynamic_bitset<>> INIT_GENOME0(GlobalPars.NPLOIDY, boost::dynamic_bitset<>(GlobalPars.NLOCI));
        std::vector<boost::dynamic_bitset<>> INIT_GENOME1(GlobalPars.NPLOIDY, boost::dynamic_bitset<>(GlobalPars.NLOCI).set());

        GlobalPars.INIT_GENOME[0] = INIT_GENOME0;
        GlobalPars.INIT_GENOME[1] = INIT_GENOME1;
    }
    AssertParameters(GlobalPars);
    ResetRGlobal(r_global, GLOBALMAX, GlobalPars.RECOMBINATIONRATE);
    ResetMGlobal(m_global, GLOBALMAX, GlobalPars.MUTATIONRATE);


    // Simulations
    for (int i = 0; i < GlobalPars.NREP; ++i)
    {
        std::cout << "Replicate: " << i << std::endl;
        RunSimulation(GlobalPars);
    }
   
    // Output
    std::ofstream Parametersoff, Dataoff;
    std::string MainFolder;
    MainFolder = CreateOutputStreams(Parametersoff,Dataoff); 

    OutputParameters(Parametersoff, GlobalPars);
	WriteOutput(Dataoff, GlobalPars);

    return 0;
}

/*
struct DataSetje
{
    DataSet(Parameters *parspointer) : pars(parspointer) {
        data.resize(pars->NPLOIDY, std::vector<int>(pars->NLOCI, 0));
        NeutralCount[0].resize(pars->NGEN,0);
        NeutralCount[1].resize(pars->NGEN,0);
        TotalPopulationSizeData[0].resize(pars->NGEN,0);
        TotalPopulationSizeData[1].resize(pars->NGEN,0);
    }

    void AddCount(const std::vector<int> data[2], const std::vector<int> popsize[2]) {
        assert(data[0].size() == data[1].size());
        assert(data[0].size() == NeutralCount[0].size());
        for(int i = 0; i < pars->NGEN; ++i){
            NeutralCount[0][i] += data[0][i];
            NeutralCount[1][i] += data[1][i];
            TotalPopulationSizeData[0][i] += popsize[0][i];
            TotalPopulationSizeData[1][i] += popsize[1][i];
       }
    }

    void AddSC_GENOME(const std::vector<double> SC_GENOME){
        SC_GENOME_data.push_back(SC_GENOME);
    }

    void WriteOutput(std::ofstream &output){
        output 
        << "Generation" << output.fill() 
        << "NeutralCount0" << output.fill() 
        << "NeutralCount1" << output.fill() 
        << "PopSize0" << output.fill() 
        << "PopSize1" << output.fill() << std::endl; 
        for(int i = 0; i < pars->NGEN; ++i){
            output
            << i << output.fill() 
            << NeutralCount[0][i] / ((double)pars->NREP * (double)pars->NPLOIDY * (double)(pars->NLOCI - pars->NLOCAL_ADAPTED_LOCI-1)) << output.fill() 
            << NeutralCount[1][i] / ((double)pars->NREP * (double)pars->NPLOIDY * (double)(pars->NLOCI - pars->NLOCAL_ADAPTED_LOCI-1)) << output.fill()
            << (double)TotalPopulationSizeData[0][i] / (double)(pars->NREP) << output.fill()
            << (double)TotalPopulationSizeData[1][i] / (double)(pars->NREP) << output.fill() << std::endl;
        }

        output << "\n \n \n \n" << std::endl;
    }

    std::vector<std::vector<double>> SC_GENOME_data;

    private:
    std::vector<std::vector<int>> data; // data[NGEN][NLOCI]
    std::vector<int> NeutralCount[2];
    std::vector<int> TotalPopulationSizeData[2];
    const Parameters *pars;
    int counter = 0;
};
*/
