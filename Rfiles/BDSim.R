library(Rcpp)
library(parallel)
library(plot3D)

setwd("/home/freek/Documents/temp2/Introgression/CppFiles/")
sourceCpp("BDSimulations.cpp")
 
# Calculate the number of cores
no_cores <- detectCores() - 1
 
# Initiate cluster
cl <- makeCluster(no_cores)

# Test parameters
testpars <- list(
tend = 1000,
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

mymatRescue = matrix(0,10,10)
mymatFix = matrix(0,10,10)
for(Ab in 1:9){
    for(aB in 0:9){
        testpars$Ab0 <- Ab
        testpars$aB0 <- aB
        data = c()
        fixcounter=0
        for (i in 1:1000) {
            res <- BDSim(tend, testpars)
            if(res$AB+res$Ab+res$aB+res$ab!=0){
                data <- c(data, (res$AB/(res$AB+res$Ab)))
            }
            else{
                fixcounter = fixcounter + 1
            }
        }
        mymatRescue[Ab+1,aB+1] <- mean(data)
        mymatFix[Ab+1,aB+1] <- fixcounter/1000
        print("done")
        print(Ab)
        print(aB)
    }
}

write.table(mymatRescue, "/home/freek/Documents/temp2/Introgression/Rfiles/data/BDSim.csv", row.names = FALSE, col.names = FALSE, sep = ',')
write.table(mymatFix, "/home/freek/Documents/temp2/Introgression/Rfiles/data/BDSimFix.csv", row.names = FALSE, col.names = FALSE, sep = ',')

