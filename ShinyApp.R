library(shiny)
library(Rcpp)
library(BH)

Rcpp.package.skeleton(
    name = "IntrogressionRcpp1", 
    cpp_files = c("MainRcpp.cpp", "random.h", "random.cpp", "utils.h", "utils.cpp"), 
    example_code = FALSE,
    author = "F.J.H. de Haas",
    email = "dehaas@zoology.ubc.ca")


library(IntrogressionRcpp1)
#                   growthrate, nloci,  nploidy,    ninit0,     ninit1,     distlocal,  scmajor,    sclocal,    ngen,   nrep,   rec,    k
data=RunSimulation( 0.2,        5,      2,          5,          5,          1,          0.0,        0.0,        50,     100,   0.5,    50)

View(data[[2]])
data[[2]][[7]] = sqrt(data[[2]][[7]])
plot(data[[2]][[1]],data[[2]][[2]])