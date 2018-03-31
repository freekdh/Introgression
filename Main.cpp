
#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <assert.h>
#include "random.h"
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <iomanip>

/// make genome circular

boost::dynamic_bitset<> r_global;
boost::dynamic_bitset<> m_global;

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

    double MultiplicativeViability(const std::vector<double> &scgenome)
    {
        const int NPLOIDY = genome.size();
        const int NLOCI = genome[0].size();
        assert(scgenome.size() == genome[0].size());
        double output = 1.0;
        for (int i = 0; i < NPLOIDY; ++i)
        {
            for (int j = 0; j < NLOCI; ++j)
            {
                if (genome[i][j] == 1)
                {
                    output *= (1.0 - scgenome[j]);
                }
            }
        }
        output /= (double)NLOCI;
        assert(0.0 <= output <= 1.0);
        return output;
    }

    void Flipbit(const int &chromosome, const int &locus) { genome[chromosome][locus].flip(); }

    bool Genotype(const int &chromosome, const int &locus) { return genome[chromosome][locus]; }
    

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

struct Parameters
{
    //default parameters
    double RECOMBINATIONRATE;
    double MUTATIONRATE;
    double FERTILITY;
    int NGEN;
    int NLOCI;
    int NPLOIDY;
    int NREP;
    int NINIT[2];
    std::vector<double> SC_GENOME;

    std::vector<boost::dynamic_bitset<>> INIT_GENOME[2];

    double SC_MAJOR;
    double SC_LOCAL;
    int NLOCAL_ADAPTED_LOCI;
};

struct DataSet
{
    DataSet(Parameters *parspointer) : pars(parspointer) {
        data.resize(pars->NPLOIDY, std::vector<int>(pars->NLOCI, 0));
    }

    void AddCount(std::vector<Individual*> &population, const int &gen) {
        for(it = population.begin(); it != population.end(); ++it){
            for(int k = 0; k < pars->NPLOIDY; ++k){
                for (int i = 0; i < pars->NLOCI; ++i){
                    data[gen][i] += (int)(*it)->Genotype(k,i);
                }
            }
        }
    }

    void Analysis(std::ofstream &output){

    }

    private:
    std::vector<std::vector<int>> data; // data[NGEN][NLOCI]
    const Parameters *pars;
    int counter = 0;
};

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
    // Constant fertility 
    const int popsize = population.size();
    std::vector<Individual*> offspring(rnd::poisson((double)popsize * pars.FERTILITY));

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
        if (rnd::uniform() < (*it)->MultiplicativeViability(pars.SC_GENOME))
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

    ofstream << "SC_GENOME" << ofstream.fill();
    for (int i = 0; i < pars.NLOCI; ++i)
    {
        ofstream << pars.SC_GENOME[i] << ofstream.fill();
    }
    ofstream << std::endl;

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
    ofstream << "FERTILITY" << ofstream.fill() << pars.FERTILITY << std::endl;
}

std::string CreateOutputStreams(std::ofstream &ParametersOfstream, std::ofstream &DataOfstream)
{
    std::string CurrentWorkingDirectory = "/home/freek/Introgression";
    //std::string Current = boost::filesystem::current_path();
    std::string MainOutputFolder = ReturnTimeStamp(CurrentWorkingDirectory);
    boost::filesystem::create_directories(MainOutputFolder.c_str());

    std::string ParameterOutput = MainOutputFolder.append("/Parameters.csv");

    ParametersOfstream.open(ParameterOutput);
    ParametersOfstream.fill(',');

    std::string DataOutput = MainOutputFolder.append("/Data.csv");

    DataOfstream.open(DataOutput);
    DataOfstream.fill(',');
    
    return MainOutputFolder;
}

void RunSimulation(Parameters &pars, DataSet &data)
{
    // Decide on major rescue locus (element 0) and element 1:n are locally favored loci
    std::vector<int> v(pars.NLOCI);
    std::iota(std::begin(v), std::end(v), 0);
    std::random_shuffle(v.begin(), v.end());

    // Set SC_GENOME for this simulation
    pars.SC_GENOME.clear();
    pars.SC_GENOME.resize(pars.NLOCI, 0.0);
    pars.SC_GENOME[v[0]] = pars.SC_MAJOR;
    for (int i = 0; i < pars.NLOCAL_ADAPTED_LOCI; ++i)
    {
        pars.SC_GENOME[v[i + 1]] = pars.SC_LOCAL;
    }

    // Initialize population
    std::vector<Individual*> population(pars.NINIT[0] + pars.NINIT[1]);
    for (it = population.begin(); it != population.begin() + pars.NINIT[0]; ++it)
    {
        *it = new Individual(pars.INIT_GENOME[0]);
    }
    for (it = population.begin() + pars.NINIT[0]; it != population.end(); ++it)
    {
        *it = new Individual(pars.INIT_GENOME[1]);
    }

    for (int i = 0; i < pars.NGEN; ++i)
    {
        data.AddCount(population,i);
        ItteratePopulation(population, pars);
    }

    for (it = population.begin(); it != population.end(); ++it)
    {
        delete *it;
    }
}

void AssertParameters(Parameters &pars){
    // Make sure all parameters are there
    assert(pars.NREP > 0);
    assert(pars.NLOCI > 0);
    assert(pars.NPLOIDY % 2 == 0 && pars.NPLOIDY > 0);
    assert(pars.NINIT[0] > 0);
    assert(pars.NINIT[1] > 0);
    assert(pars.NGEN > 0);
    assert(1.0 >= pars.SC_LOCAL >= 0.0);
    assert(1.0 >= pars.SC_MAJOR >= 0.0);
    assert(pars.NLOCAL_ADAPTED_LOCI >= 0);
    assert(pars.NLOCAL_ADAPTED_LOCI + 1 < pars.NLOCI);
}

int main(int argc, char *argv[])
{
    if (argc != 10) {std::cout << "Argc != 10" << std::endl; return -1;}

    const int GLOBALMAX = 100000;

    // Parameters
    Parameters pars;
    pars.RECOMBINATIONRATE = 0.5;
    pars.MUTATIONRATE = 0.00001;
    pars.FERTILITY = 1.1;

    pars.NLOCI = atoi(argv[1]);
    pars.NPLOIDY = atoi(argv[2]);
    pars.NINIT[0] = atoi(argv[3]);
    pars.NINIT[1] = atoi(argv[4]);
    pars.NLOCAL_ADAPTED_LOCI = atoi(argv[5]);
    pars.SC_MAJOR = atof(argv[6]);
    pars.SC_LOCAL = atof(argv[7]);
    pars.NGEN = atoi(argv[8]);
    pars.NREP = atoi(argv[9]);

    pars.SC_GENOME.clear();
    pars.SC_GENOME.resize(pars.NLOCI, 0.0);
    std::vector<boost::dynamic_bitset<>> INIT_GENOME0(pars.NPLOIDY, boost::dynamic_bitset<>(pars.NLOCI, false));
    std::vector<boost::dynamic_bitset<>> INIT_GENOME1(pars.NPLOIDY, boost::dynamic_bitset<>(pars.NLOCI, true));
    pars.INIT_GENOME[0] = INIT_GENOME0;
    pars.INIT_GENOME[1] = INIT_GENOME1;

    AssertParameters(pars);

    // Initialize global variables
    rnd::set_seed();
    ResetRGlobal(r_global, GLOBALMAX, pars.RECOMBINATIONRATE);
    ResetMGlobal(m_global, GLOBALMAX, pars.MUTATIONRATE);

    DataSet data(&pars);
    for (int i = 0; i < pars.NREP; ++i)
    {
        std::cout << "Replicate: " << i << std::endl;
        RunSimulation(pars, data); // input, output
    }

    // DATA OFFSTREAM
    std::ofstream Parametersoff, Dataoff;
    std::string MainFolder;
    MainFolder = CreateOutputStreams(Parametersoff,Dataoff); 

    OutputParameters(Parametersoff,pars);
	data.Analysis(Dataoff);

    return 0;
}
