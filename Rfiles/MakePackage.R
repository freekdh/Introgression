library(Rcpp)
paste(getwd(),"/Cppfiles", sep="")
# Make R package
Rcpp.package.skeleton(
    name = "pkgIntrogressionv2`", 
    code_files = c("Rfiles/ShinyExport.R", "Rfiles/ShinyApp/server.r", "Rfiles/ShinyApp/ui.r"),
    cpp_files = c("CppFiles/IntrogressionSimulations.cpp",
    "CppFiles/IntrogressionSimulations.h",
    "CppFiles/Rcpp_output.cpp",
    "CppFiles/Rcpp_output.h",
    "CppFiles/Shiny_output.cpp", 
    "CppFiles/random.h", 
    "CppFiles/random.cpp", 
    "CppFiles/utils.h", 
    "CppFiles/utils.cpp",
    "CppFiles/Makevars"), 
    example_code = FALSE,
    author = "F.J.H. de Haas",
    email = "dehaas@zoology.ubc.ca")
    
#Add to the DESCRIPTION file:
#Imports: Rcpp (>= 0.12.17), BH, RcppProgress
#LinkingTo: Rcpp, BH, RcppProgress

