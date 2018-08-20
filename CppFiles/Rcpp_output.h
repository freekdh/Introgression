#ifndef RCPPOUTPUT_H
#define RCPPOUTPUT_H

#include <Rcpp.h>
#include "IntrogressionSimulations.h"

Rcpp::List Rcpp_WriteOutput(const Parameters &GlobalPars, SimData &SimulationData);

#endif

