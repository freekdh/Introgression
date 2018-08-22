library(pkgIntrogression)

testpars <- list(
    b=1.1,
    dA=1,
    da=0.8,a
    r=0.05,
    nloci=100,
    ninit0=20,
    ninit1=2,
    ngen=20,
    nrep=100,
    rec=0.05,
    k=100
)

out <- matrix(0,nrow = 10, ncol = 10)
for(i in 1:10){
    for(j in 1:10){
        testpars$ninit0 <- i
        testpars$ninit1 <- j
        testdata <- RcppIntrogressionSimulation(testpars,0)
        out[i,j]<-1-tail(testdata$data$Introgressed1_avg,n=1)
    }
}

write.table(out,file = sprintf("./Rfiles/data/SimData_r%s.csv", 100*testpars$r),row.names=FALSE,col.names=FALSE,sep=",")


### Increased growth plot simulations

testpars <- list(
    b=1.1,
    dA=1,
    da=0.8,
    r=0.05,
    nloci=100,
    ninit0=20,
    ninit1=1,
    ngen=50,
    nrep=1000,
    rec=0.5,
    k=5000
)

testdata <- RcppIntrogressionSimulation(testpars,0)
testdata$data$Major1_avg

plot(log(testdata$data$Major1_avg))

write.table(log(testdata$data$Major1_avg),"./Rfiles/data/loggrowthA.csv",row.names=FALSE,col.names=FALSE,sep=",")


### Fixation probabilities data

testpars <- list(
    b=1.1,
    dA=1,
    da=0.8,
    nloci=100,
    ninit0=50,
    ninit1=1,
    ngen=50,
    nrep=1000,
    rec=0.5,
    k=500
)

out <- matrix(0,nrow = 20, ncol = 2)
index = 1
for(i in 1:20){
    testpars$ninit0 <- 20
    testpars$ninit1 <- i
    testpars$b <- 1.01
    testdata <- RcppIntrogressionSimulation(testpars,0)
    out[index,1]<-i
    out[index,2]<-testdata$fixation
    index = index+1
}
out
write.table(out,"./Rfiles/data/fixationsimulations.csv",row.names=FALSE,col.names=FALSE,sep=",")
