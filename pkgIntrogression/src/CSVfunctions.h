#ifndef CSVFUNCTION_H
#define CSVFUNCTION_H

#include <iostream>
#include <fstream>
#include "IntrogressionSimulations.h"

void CSV_WriteOutput(const Parameters &pars, SimData &SimulationData, std::ofstream arrayofstream[4]){};

std::string CreateOutputStreams(std::ofstream arrayofstream[4]){};

int CSV_IntrogressionSimulation(int argc, char *argv[]) {};

Parameters::Parameters(double r, int nloci, int nploidy, int ninit0, int ninit1, int distlocal, double scmajor, double sclocal, int ngen, int nrep, double rec, int k, int threads = 0){};

#endif
