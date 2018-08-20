#ifndef CSVFUNCTION_H
#define CSVFUNCTION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include "IntrogressionSimulations.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/filesystem.hpp>
#include <stdio.h>
#include <iomanip>
#include <chrono>
#include "random.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

enum nameofstream {Parametersoff, Dataoff, AlleleFrequencyoff_mean, AlleleFrequencyoff_var};

void CSV_WriteOutput(std::ofstream arrayofstream[4], const Parameters &pars, SimData &SimulationData){

    using namespace boost::accumulators;
    // Parameters output
    for (int i = 0; i < 2; ++i)
    {
        arrayofstream[Parametersoff] << "NINIT" << i << arrayofstream[Parametersoff].fill() << pars.NINIT[i] << std::endl;
    }
    arrayofstream[Parametersoff] << "NGEN" << arrayofstream[Parametersoff].fill() << pars.NGEN << std::endl;
    arrayofstream[Parametersoff] << "NLOCI" << arrayofstream[Parametersoff].fill() << pars.NLOCI << std::endl;
    arrayofstream[Parametersoff] << "NREP" << arrayofstream[Parametersoff].fill() << pars.NREP << std::endl;
    arrayofstream[Parametersoff] << "MUTATIONRATE" << arrayofstream[Parametersoff].fill() << pars.MUTATIONRATE << std::endl;
    arrayofstream[Parametersoff] << "RECOMBINATIONRATE" << arrayofstream[Parametersoff].fill() << pars.RECOMBINATIONRATE << std::endl;
    arrayofstream[Parametersoff] << "CARRYINGCAPACITY" << arrayofstream[Parametersoff].fill() << pars.K << std::endl;
    arrayofstream[Parametersoff] << "BIRTHRATE" << arrayofstream[Parametersoff].fill() << pars.BIRTHRATE << std::endl;
    arrayofstream[Parametersoff] << "DEATHRATEA" << arrayofstream[Parametersoff].fill() << pars.DEATHRATEA << std::endl;
    arrayofstream[Parametersoff] << "DEATHRATEa" << arrayofstream[Parametersoff].fill() << pars.DEATHRATEa << std::endl;
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
            popsize(SimulationData.DataSet[j]->popsize[i]);
            major0(SimulationData.DataSet[j]->major0[i]);
            major1(SimulationData.DataSet[j]->major1[i]);
            introgressed0(SimulationData.DataSet[j]->introgressed0[i]);
            introgressed1(SimulationData.DataSet[j]->introgressed1[i]);
        }

        arrayofstream[Dataoff]
        << i << arrayofstream[Dataoff].fill() 
        << mean(popsize) << arrayofstream[Dataoff].fill()
        << (double)mean(major0) << arrayofstream[Dataoff].fill()
        << (double)mean(major1) << arrayofstream[Dataoff].fill() 
        << mean(introgressed0)/(double)(pars.NLOCI-1) << arrayofstream[Dataoff].fill()
        << mean(introgressed1)/(double)(pars.NLOCI-1) << arrayofstream[Dataoff].fill()

        << variance(popsize) << arrayofstream[Dataoff].fill()
        << (double)variance(major0) << arrayofstream[Dataoff].fill()
        << (double)variance(major1) << arrayofstream[Dataoff].fill()
        << variance(introgressed0)/(double)((pars.NLOCI-1)*(pars.NLOCI-1)) << arrayofstream[Dataoff].fill()
        << variance(introgressed1)/(double)((pars.NLOCI-1)*(pars.NLOCI-1)) << arrayofstream[Dataoff].fill() << std::endl;
    };
    
    for(int i = 0; i < pars.NGEN; ++i){
        arrayofstream[AlleleFrequencyoff_mean] << i << arrayofstream[AlleleFrequencyoff_mean].fill(); 
        arrayofstream[AlleleFrequencyoff_var] << i << arrayofstream[AlleleFrequencyoff_var].fill();
        for(int l = 0; l < pars.NLOCI; ++l)
        {
            accumulator_set<double, stats<tag::mean, tag::variance > > locus;
            for(int r = 0; r < pars.NREP; ++r)
            {
                locus((double)SimulationData.DataSet[r]->allele0[i][l] / ((double)SimulationData.DataSet[r]->popsize[i]));
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
    CurrentWorkingDirectory.append("/CppFiles/Data");

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
    rnd::set_seed();
    const Parameters GlobalPars(argc, argv);

    SimData SimulationData;

    // Run nrep successful simulations
    // auto start = std::chrono::high_resolution_clock::now();
    #ifdef _OPENMP
        printf("Parallel activated : Number of threads=%i\n",omp_get_max_threads());   
    #endif

    #pragma omp parallel for schedule(static)
    for(int task = 0; task < GlobalPars.NREP; ++task){
        while(RunSimulation(GlobalPars, SimulationData)==false);
    }
    // auto finish = std::chrono::high_resolution_clock::now();
   
    // Write outputfiles
    std::ofstream arrayofstream[4]; 
    CreateOutputStreams(arrayofstream); 
	CSV_WriteOutput(arrayofstream, GlobalPars, SimulationData);

    // Cleanup
    for(int i = 0; i < GlobalPars.NREP; ++i){
        delete SimulationData.DataSet[i];
    }

    // Show Time elapsed
    /*
    std::chrono::duration<double> elapsed = finish-start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
    */
    return 0;
}

#endif
