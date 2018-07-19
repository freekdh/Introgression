
#include "CSVfunctions.h"

enum nameofstream {Parametersoff, Dataoff, AlleleFrequencyoff_mean, AlleleFrequencyoff_var};

Parameters::Parameters(double r, int nloci, int nploidy, int ninit0, int ninit1, int distlocal, double scmajor, double sclocal, int ngen, int nrep, double rec, int k, int threads){
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

        Initialize();
}

void CSV_WriteOutput(const Parameters &pars, SimData &SimulationData, std::ofstream arrayofstream[4]){
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

int CSV_IntrogressionSimulation(int argc, char *argv[]){

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
	WriteOutputCSV(arrayofstream, GlobalPars);

    // Cleanup
    for(DataBlock* i : DataSet) delete i;
    DataSet.clear();

    // Show Time elapsed
    std::chrono::duration<double> elapsed = finish-start;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
    return 0;
}

