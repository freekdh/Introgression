#ifndef SHINYFUNCTION_H
#define SHINYFUNCTION_H

void Shiny_InitializeSimulation(double r, int nloci, int nploidy, int ninit0, int ninit1, int distlocal, double scmajor, double sclocal, int ngen, double rec, int k){};

void Shiny_RunSimulation(){};

Rcpp::List Shiny_WriteOutputandCleanup(){};

#endif
