// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <progress.hpp>
#include "IntrogressionSimulations.h"

// [[Rcpp::export]]
void Shiny_InitializeSimulation(double r, int nloci, int nploidy, int ninit0, int ninit1, int distlocal, double scmajor, double sclocal, int ngen, double rec, int k){
    if(RUNSIMULATION_FUN == true) {std::cerr << "First write output or clear dataset : WriteOutputandCleanupt()" << std::endl;}
    else{
        // Prepare for simulations
        CollectParameters( r,  nloci,  nploidy,  ninit0,  ninit1,  distlocal,  scmajor,  sclocal,  ngen,  rec,  k);
        InitializeGlobalEnvironment(GlobalPars);
        INITIALIZESIMULATION_FUN = true;
    } 
}

// [[Rcpp::export]]
void Shiny_RunSimulation(){
    if(INITIALIZESIMULATION_FUN == true){
        // Run nrep successful simulations
        while(RunSimulation(GlobalPars)==false);        
        RUNSIMULATION_FUN = true;
    }
    else{std::cerr << "First initialize simulation: InitializeSimulation()" << std::endl;}
}

// [[Rcpp::export]]
Rcpp::List Shiny_WriteOutputandCleanup(){

    if(INITIALIZESIMULATION_FUN == true && RUNSIMULATION_FUN == true){
        // Parameters
        const int NREP = DataSet.size();
        Rcpp::DataFrame parsdata =  Rcpp::DataFrame::create(
            Rcpp::_["NINIT0"] = GlobalPars.NINIT[0], 
            Rcpp::_["NINIT1"] = GlobalPars.NINIT[1],
            Rcpp::_["NGEN"] = GlobalPars.NGEN,
            Rcpp::_["NREP"] = NREP,
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

            for(int j = 0; j < NREP; ++j){
                popsize(DataSet[j]->popsize[i]);
                major0((double)DataSet[j]->major0[i]/(double)GlobalPars.NPLOIDY);
                major1((double)DataSet[j]->major1[i]/(double)GlobalPars.NPLOIDY);
                introgressed0(DataSet[j]->introgressed0[i]);
                introgressed1(DataSet[j]->introgressed1[i]);
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
                for(int r = 0; r < NREP; ++r)
                {
                    locus((double)DataSet[r]->allele0[i][l] / ((double)DataSet[r]->popsize[i] * (double)GlobalPars.NPLOIDY));
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

        // Fixation probability:
        double fixation = (double)fix[true]/(double)(fix[true]+fix[false]);

        // Cleanup
        for(int i = 0; i < NREP; ++i){
            delete DataSet[i];
        }
        DataSet.clear();
        fix[true] = 0;
        fix[true] = 0;

        INITIALIZESIMULATION_FUN = false;
        RUNSIMULATION_FUN = false;
        
        return Rcpp::List::create(
            Rcpp::_["pars"] = (parsdata),
            Rcpp::_["data"] = (data),
            Rcpp::_["allelefavg"] = (alleledatamean),
            Rcpp::_["allelefvar"] = (alleledatavar),
            Rcpp::_["fixation"] = fixation
        );
    }

    else{
        std::cerr << "First initialize and run a simulation: InitializeSimulation() and RunSimulation()" << std::endl;
        return NULL;
    }
}
