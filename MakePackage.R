library(Rcpp)
paste(getwd(),"/Cppfiles", sep="")
# Make R package
Rcpp.package.skeleton(
    name = "pkgIntrogression", 
    cpp_files = c("sourcefiles/IntrogressionSimulations.cpp",
    "sourcefiles/IntrogressionSimulations.h",
    "sourcefiles/Rcpp_output.cpp",
    "sourcefiles/Shiny_output.cpp", 
    "sourcefiles/random.h", 
    "sourcefiles/random.cpp", 
    "sourcefiles/utils.h", 
    "sourcefiles/utils.cpp"), 
    example_code = FALSE,
    author = "F.J.H. de Haas",
    email = "dehaas@zoology.ubc.ca")
    
#Add to the DESCRIPTION file:
#Imports: Rcpp (>= 0.12.17), BH, RcppProgress
#LinkingTo: Rcpp, BH, RcppProgress

#And add the MakeVars file with the openmp flag