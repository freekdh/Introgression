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

extern Rcpp::List Rcpp_WriteOutput(const Parameters &GlobalPars, SimData &SimulationData);

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
        Rcpp_WriteOutput(Shiny_Pars,Shiny_Data);
    }
    else{
        std::cerr << "First initialize and run a simulation: InitializeSimulation() and RunSimulation()" << std::endl;
        return NULL;
    }
}

#endif