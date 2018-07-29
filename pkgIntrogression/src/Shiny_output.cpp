// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(BH)]]

#ifndef SHINYFUNCTION_H
#define SHINYFUNCTION_H

#include "IntrogressionSimulations.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <progress.hpp>
#include <Rcpp.h>

static bool RUNSIMULATION_FUN = false, INITIALIZESIMULATION_FUN = false;
static Parameters Shiny_Pars;
static SimData Shiny_Data; 

// [[Rcpp::export]]
void ShinyInitializeSimulation(const Rcpp::List &parslist){
    if(RUNSIMULATION_FUN == true) {std::cerr << "First write output or clear dataset : WriteOutputandCleanupt()" << std::endl;}
    else{
        // Prepare for simulations
        Parameters newpars(parslist);
        Shiny_Pars = newpars;
        INITIALIZESIMULATION_FUN = true;
    } 
}

// [[Rcpp::export]]
void ShinyRunSimulation(){
    if(INITIALIZESIMULATION_FUN == true){
        // Run nrep successful simulations
        while(RunSimulation(Shiny_Pars, Shiny_Data)==false);        
        RUNSIMULATION_FUN = true;
    }
    else{std::cerr << "First initialize simulation: InitializeSimulation()" << std::endl;}
}

// [[Rcpp::export]]
Rcpp::List ShinyWriteOutputandCleanup(){
    if(INITIALIZESIMULATION_FUN == true && RUNSIMULATION_FUN == true){
        // Parameters
        using namespace boost::accumulators;
        const int NREP = Shiny_Data.DataSet.size();
        Rcpp::DataFrame parsdata =  Rcpp::DataFrame::create(
            Rcpp::_["NINIT0"] = Shiny_Pars.NINIT[0], 
            Rcpp::_["NINIT1"] = Shiny_Pars.NINIT[1],
            Rcpp::_["NGEN"] = Shiny_Pars.NGEN,
            Rcpp::_["NREP"] = NREP,
            Rcpp::_["NLOCI"] = Shiny_Pars.NLOCI,
            Rcpp::_["NPLOIDY"] = Shiny_Pars.NPLOIDY,
            Rcpp::_["MUTATIONRATE"] = Shiny_Pars.MUTATIONRATE,
            Rcpp::_["RECOMBINATIONRATE"] = Shiny_Pars.RECOMBINATIONRATE,
            Rcpp::_["GROWTHRATE"] = Shiny_Pars.INTRINSIC_GROWTHRATE,
            Rcpp::_["CARRYINGCAPACITY"] = Shiny_Pars.K,
            Rcpp::_["SC_MAJOR"] = Shiny_Pars.SC_MAJOR,
            Rcpp::_["SC_LOCAL"] = Shiny_Pars.SC_LOCAL,
            Rcpp::_["INDEX_MAJOR"] = Shiny_Pars.index[0],
            Rcpp::_["INDEX_LOCAL"] = Shiny_Pars.index[1],
            Rcpp::_["DISTLOCAL"] = Shiny_Pars.DISTLOCAL
        );
    
        // Write outputfiles
        Rcpp::NumericVector generation(Shiny_Pars.NGEN);    
        Rcpp::NumericVector popsizevecm(Shiny_Pars.NGEN);
        Rcpp::NumericVector major0vecm(Shiny_Pars.NGEN);
        Rcpp::NumericVector major1vecm(Shiny_Pars.NGEN);
        Rcpp::NumericVector introgressed0vecm(Shiny_Pars.NGEN);
        Rcpp::NumericVector introgressed1vecm(Shiny_Pars.NGEN);
        Rcpp::NumericVector popsizevecv(Shiny_Pars.NGEN);
        Rcpp::NumericVector major0vecv(Shiny_Pars.NGEN);
        Rcpp::NumericVector major1vecv(Shiny_Pars.NGEN);
        Rcpp::NumericVector introgressed0vecv(Shiny_Pars.NGEN);
        Rcpp::NumericVector introgressed1vecv(Shiny_Pars.NGEN);
        

        for(int i = 0; i < Shiny_Pars.NGEN; ++i){
            accumulator_set<int, stats<tag::mean, tag::variance > > popsize;
            accumulator_set<double, stats<tag::mean, tag::variance > > major0;
            accumulator_set<double, stats<tag::mean, tag::variance > > major1;
            accumulator_set<double, stats<tag::mean, tag::variance > > introgressed0;
            accumulator_set<double, stats<tag::mean, tag::variance > > introgressed1;

            for(int j = 0; j < NREP; ++j){
                popsize(Shiny_Data.DataSet[j]->popsize[i]);
                major0((double)Shiny_Data.DataSet[j]->major0[i]/(double)Shiny_Pars.NPLOIDY);
                major1((double)Shiny_Data.DataSet[j]->major1[i]/(double)Shiny_Pars.NPLOIDY);
                introgressed0(Shiny_Data.DataSet[j]->introgressed0[i]);
                introgressed1(Shiny_Data.DataSet[j]->introgressed1[i]);
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

        std::vector<Rcpp::NumericVector> allelefrequencymean(Shiny_Pars.NLOCI);
        std::vector<Rcpp::NumericVector> allelefrequencyvar(Shiny_Pars.NLOCI);

        for(int i = 0; i < Shiny_Pars.NLOCI; ++i){
            allelefrequencymean[i] = Rcpp::NumericVector(Shiny_Pars.NGEN); 
            allelefrequencyvar[i] = Rcpp::NumericVector(Shiny_Pars.NGEN); 
        }

        for(int i = 0; i < Shiny_Pars.NGEN; ++i){
            for(int l = 0; l < Shiny_Pars.NLOCI; ++l)
            {
                accumulator_set<double, stats<tag::mean, tag::variance > > locus;
                for(int r = 0; r < NREP; ++r)
                {
                    locus((double)Shiny_Data.DataSet[r]->allele0[i][l] / ((double)Shiny_Data.DataSet[r]->popsize[i] * (double)Shiny_Pars.NPLOIDY));
                }
                allelefrequencymean[l][i] = mean(locus);
                allelefrequencyvar[l][i] = variance(locus);
            }
        }   

        Rcpp::DataFrame alleledatamean, alleledatavar;
        std::string avglocus = "avglocus", varlocus = "varlocus";
        alleledatamean.push_back(generation,"Generation");
        alleledatavar.push_back(generation,"Generation");
        for(int i = 0; i < Shiny_Pars.NLOCI; ++i){
            alleledatamean.push_back((allelefrequencymean[i]), avglocus+std::to_string(i));
            alleledatavar.push_back((allelefrequencyvar[i]), varlocus+std::to_string(i));
        }

        // Fixation probability:
        double fixation = (double)Shiny_Pars.NREP/(double(Shiny_Pars.NREP+Shiny_Data.nofixcounter.load()));

        // Cleanup
        for(int i = 0; i < NREP; ++i){
            delete Shiny_Data.DataSet[i];
        }
        Shiny_Data.DataSet.clear();
        Shiny_Data.nofixcounter = 0;

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

#endif