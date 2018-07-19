library(Rcpp)

# Make R package
Rcpp.package.skeleton(
    name = "IntrogressionRcpp", 
    cpp_files = c("MainRcpp_Multithread.cpp", "random.h", "random.cpp", "utils.h", "utils.cpp"), 
    example_code = FALSE,
    author = "F.J.H. de Haas",
    email = "dehaas@zoology.ubc.ca")
#Add to the DESCRIPTION file:
#LinkingTo: Rcpp, BH, RcppProgress
#Depends: BH, RcppProgress
# Something