library(Rcpp)
setwd("/home/freek/Documents/temp2/Introgression/CppFiles/")
sourceCpp("BDSimulations.cpp")

testpars <- list(
tend = 10000,
AB0 = 0,
Ab0 = 1,
aB0 = 10, 
ab0 = 0,
bA = 1,
ba = 1,
dA = 2-1.1,
da = 2-0.9,
r = 0.5
)

data = c()
for (i in 1:10000) {
    res <- BDSim(tend, testpars)
    if(res$AB+res$Ab+res$aB+res$ab!=0){
        data <- c(data, (res$AB/(res$AB+res$Ab)))
    }
}
mean(data)