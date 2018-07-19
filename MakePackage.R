library(Rcpp)
paste(getwd(),"/Cppfiles", sep="")
# Make R package
Rcpp.package.skeleton(
    name = "pkgIntrogression2", 
    cpp_files = c("Cppfiles/MainRcpp_Multithread.cpp", "Cppfiles/random.h", "Cppfiles/random.cpp", "Cppfiles/utils.h", "Cppfiles/utils.cpp"), 
    example_code = FALSE,
    author = "F.J.H. de Haas",
    email = "dehaas@zoology.ubc.ca")
#Add to the DESCRIPTION file:
#LinkingTo: Rcpp, BH, RcppProgress
#Depends: BH, RcppProgress
# Something