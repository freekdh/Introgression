library(Rcpp)
paste(getwd(),"/Cppfiles", sep="")
# Make R package
Rcpp.package.skeleton(
    name = "pkgIntrogression", 
    code_files = c("ShinyApp/ShinyExport.R", "ShinyApp/server.r", "ShinyApp/ui.r"),
    cpp_files = c("CppFiles/IntrogressionSimulations.cpp",
    "CppFiles/IntrogressionSimulations.h",
    "CppFiles/Rcpp_output.cpp",
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

