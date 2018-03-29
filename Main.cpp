
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
    Individual(const std::vector<boost::dynamic_bitset<>> &INITGENOME) : genome(INITGENOME) {}

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
    }

    std::vector<boost::dynamic_bitset<>> genome; // circular genome of individual
};

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
    DataSet(Parameters *parspointer) : pars(parspointer)
    {
    }

    void AddCount()
    {
    }

    void Analysis()
    {
    }

  private:
    const Parameters *pars;
    int counter = 0;
};

std::vector<Individual *>::iterator it;

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

void ItteratePopulation(std::vector<Individual *> &population, const double &FERTILITY, const std::vector<double> &SC_GENOME)
{
    const int popsize = population.size();
    assert(popsize > 0);
    std::vector<Individual *> offspring;
    offspring.resize((double)popsize * FERTILITY);
    assert(offspring.size() > 0);

    for (it = offspring.begin(); it != offspring.end(); ++it)
    {
        Individual c[2] = {*population[rnd::integer(popsize)], *population[rnd::integer(popsize)]};
        *it = new Individual(c);
    }

    for (it = population.begin(); it != population.end(); ++it)
    {
        delete *it;
    }
    population.clear();

    for (it = offspring.begin(); it != offspring.end(); ++it)
    {
        if (rnd::uniform() < (*it)->MultiplicativeViability(SC_GENOME))
        {
            population.push_back(*it);
        };
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

std::string CreateOutputStreams(std::ofstream &ParametersOfstream, std::vector<std::ofstream> &MetaOfstream, std::vector<std::vector<std::ofstream>> &SubOfStream)
{
    std::string CurrentWorkingDirectory = "/home/freek/PopGenDensityDependence";
    //std::string Current = boost::filesystem::current_path();
    std::string MainOutputFolder = ReturnTimeStamp(CurrentWorkingDirectory);
    std::string returnmainfolder = MainOutputFolder;
    std::string MetaOutputFolder = MainOutputFolder;
    std::string SubOutputFolder = MainOutputFolder;
    MetaOutputFolder.append("/MetaPopulation");
    SubOutputFolder.append("/SubPopulations");
    boost::filesystem::create_directories(MainOutputFolder.c_str());
    boost::filesystem::create_directories(MetaOutputFolder.c_str());
    boost::filesystem::create_directories(SubOutputFolder.c_str());

    ParametersOfstream.open(MainOutputFolder.append("/Parameters.csv"));

    // MetaFiles:

    for (unsigned int i = 0; i < MetaOfstream.size(); ++i)
    {
        std::string focallocus = MetaOutputFolder;
        focallocus.append("/Locus");
        focallocus = focallocus + std::to_string(i);
        focallocus = focallocus.append(".csv");
        MetaOfstream[i].open(focallocus);
        assert(MetaOfstream[i].is_open());
        MetaOfstream[i].fill(',');
    }

    // SubFiles:
    for (unsigned int sub = 0; sub < SubOfStream.size(); ++sub)
    {
        std::string focal = SubOutputFolder;
        focal.append("/Sub");
        focal = focal + std::to_string(sub);
        boost::filesystem::create_directories(focal.c_str());
        for (unsigned loc = 0; loc < SubOfStream[sub].size(); ++loc)
        {
            std::string temp = focal;
            temp.append("/Locus");
            temp = temp + std::to_string(loc);
            temp = temp.append(".csv");
            SubOfStream[sub][loc].open(temp);
            assert(SubOfStream[sub][loc].is_open());
            SubOfStream[sub][loc].fill(',');
        }
    }

    return returnmainfolder;
}

void RunSimulation(Parameters &pars, DataSet &data)
{
    // Decide on major rescue gene:
    std::vector<int> v(pars.NLOCI);
    std::iota(std::begin(v), std::end(v), 0);
    std::random_shuffle(v.begin(), v.end());

    pars.SC_GENOME[v[0]] = pars.SC_MAJOR;
    for (int i = 0; i < pars.NLOCAL_ADAPTED_LOCI; ++i)
    {
        pars.SC_GENOME[v[i + 1]] = pars.SC_LOCAL;
    }

    // Initialize population
    std::vector<Individual *> population(pars.NINIT[0] + pars.NINIT[1]);
    for (it = population.begin(); it != population.begin() + pars.NINIT[0]; ++it)
    {
        *it = new Individual(pars.INIT_GENOME[0]);
    }
    for (it = population.begin() + pars.NINIT[0]; it != population.end(); ++it)
    {
        *it = new Individual(pars.INIT_GENOME[1]);
    }


}

int main(int argc, char *argv[])
{

    const int GLOBALMAX = 100000;

    // Parameters
    Parameters pars;
    pars.RECOMBINATIONRATE = 0.5;
    pars.MUTATIONRATE = 0.00001;
    pars.FERTILITY = 1.1;

    pars.NLOCI = atoi(argv[1]);
    assert(pars.NLOCI > 0);
    pars.NPLOIDY = atoi(argv[2]);
    assert(pars.NPLOIDY % 2 == 0 && pars.NPLOIDY > 0);
    pars.NINIT[0] = atoi(argv[3]);
    assert(pars.NINIT[0] > 0);
    pars.NINIT[1] = atoi(argv[4]);
    assert(pars.NINIT[1] > 0);
    pars.NLOCAL_ADAPTED_LOCI = atoi(argv[5]);
    assert(pars.NLOCAL_ADAPTED_LOCI >= 0);
    assert(pars.NLOCAL_ADAPTED_LOCI + 1 < pars.NLOCI);
    pars.SC_MAJOR = atof(argv[6]);
    assert(1.0 >= pars.SC_MAJOR >= 0.0);
    pars.SC_LOCAL = atof(argv[7]);
    assert(1.0 >= pars.SC_LOCAL >= 0.0);
    pars.NGEN = atoi(argv[8]);
    assert(pars.NGEN > 0);
    pars.NREP = atoi(argv[9]);
    assert(pars.NREP > 0);

    pars.SC_GENOME.clear();
    pars.SC_GENOME.resize(pars.NLOCI, 0.0);

    std::vector<boost::dynamic_bitset<>> INIT_GENOME0(pars.NPLOIDY, boost::dynamic_bitset<>(pars.NLOCI, false));
    std::vector<boost::dynamic_bitset<>> INIT_GENOME1(pars.NPLOIDY, boost::dynamic_bitset<>(pars.NLOCI, true));
    pars.INIT_GENOME[0] = INIT_GENOME0;
    pars.INIT_GENOME[1] = INIT_GENOME1;

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
    std::ofstream output("output.csv");
    output.fill(',');
    output << "gen" << output.fill() << "popsize" << output.fill();
    for (int i = 0; i < pars.NLOCI; ++i)
    {
        output << "locus" << i << output.fill();
    }
    output << std::endl;

    return 0;
}

/* 
GRAVEYARD:

void StationaryPopulation(std::vector<Individual *> &population, const double &FERTILITY, const std::vector<double> &SC_GENOME)
{
    const int popsize = population.size();
    assert(popsize > 0);
    std::vector<Individual *> offspring(popsize);

    rnd::discrete_distribution viabilitydist(popsize);
    for (int i = 0; i < popsize; ++i)
    {
        viabilitydist[i] = population[i]->MultiplicativeViability(SC_GENOME);
    }

    for (it = offspring.begin(); it != offspring.end(); ++it)
    {
        Individual c[2] = {*population[viabilitydist.sample()], *population[viabilitydist.sample()]};
        *it = new Individual(c);
    }

    for (it = population.begin(); it != population.end(); ++it)
    {
        delete *it;
    }
    population.clear();

    population = offspring;
}
*/
