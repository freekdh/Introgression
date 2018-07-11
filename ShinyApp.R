library(shiny)
library(Rcpp)
library(BH)

Rcpp.package.skeleton(
    name = "mypackage4", 
    cpp_files = c("MainRcpp.cpp", "random.h", "random.cpp", "utils.h", "utils.cpp"), 
    example_code = FALSE,
    author = "F.J.H. de Haas",
    email = "dehaas@zoology.ubc.ca")


library(mypackage4)
a=RunSimulation(0.1, 10, 2, 20, 20, 1, 0.0, 0.0, 10, 10, 0.5, 100)
a[[3]]