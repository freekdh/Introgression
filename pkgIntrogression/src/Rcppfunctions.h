#ifndef RCPPFUNCTION_H
#define RCPPFUNCTION_H

#include <Rcpp.h>
#include <progress.hpp>
#include "IntrogressionSimulations.h"

Parameters::Parameters(const Rcpp::List &parslist){};

Rcpp::List Rcpp_WriteOutput(const Parameters &GlobalPars, SimData &SimulationData){};

Rcpp::List Rcpp_IntrogressionSimulation(Rcpp::List parslist, int setthreads = 0, bool progressbar = true){};

#endif
