// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <progress.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include "IntrogressionSimulations.h"
#include "random.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

Rcpp::List Rcpp_WriteOutput(const Parameters &GlobalPars, SimData &SimulationData){
    using namespace boost::accumulators;

    Rcpp::DataFrame parsdata =  Rcpp::DataFrame::create(
        Rcpp::_["NINIT0"] = GlobalPars.NINIT[0], 
        Rcpp::_["NINIT1"] = GlobalPars.NINIT[1],
        Rcpp::_["NGEN"] = GlobalPars.NGEN,
        Rcpp::_["NREP"] = GlobalPars.NREP,
        Rcpp::_["NLOCI"] = GlobalPars.NLOCI,
        Rcpp::_["RECOMBINATIONRATE"] = GlobalPars.RECOMBINATIONRATE,
        Rcpp::_["GROWTHRATE"] = GlobalPars.BIRTHRATE,
        Rcpp::_["CARRYINGCAPACITY"] = GlobalPars.K,
        Rcpp::_["INDEX_MAJOR"] = GlobalPars.index[0]
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
            major0((double)SimulationData.DataSet[j]->major0[i]);
            major1((double)SimulationData.DataSet[j]->major1[i]);
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
                locus((double)SimulationData.DataSet[r]->allele0[i][l] / ((double)SimulationData.DataSet[r]->popsize[i]));
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
Rcpp::List RcppIntrogressionSimulation(Rcpp::List parslist, int setthreads = 0, bool progressbar = false){
    
    // Prepare for simulation
    rnd::set_seed();
    const Parameters GlobalPars(parslist);
    SimData SimulationData;

    // Run nrep successful simulations
    #ifdef _OPENMP
        const static int maxthreads = omp_get_max_threads();
        if(setthreads>0) omp_set_num_threads(setthreads);
        else omp_set_num_threads(maxthreads);
        REprintf("Parallel activated : Number of threads=%i\n",omp_get_max_threads());   
    #endif
    //if(progressbar == true) Progress p(GlobalPars.NREP, true);

    #pragma omp parallel for schedule(static)
    for (int task = 0; task < GlobalPars.NREP; ++task){
        while(RunSimulation(GlobalPars, SimulationData)==false);
        //if(progressbar == true) p.increment();
    }

    return Rcpp_WriteOutput(GlobalPars, SimulationData);
}
