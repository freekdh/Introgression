library(pkgIntrogression)
testpars <- list(
    b=1.1,
    dA=1,
    da=0.8,
    r=0.0005,
    nloci=100,
    ninit0=20,
    ninit1=2,
    ngen=20,
    nrep=1000,
    rec=0.5,
    k=20)

ShinyInitializeSimulation(testpars)
ShinyRunSimulation()
ShinyWriteOutputandCleanup()

testdata <- RcppIntrogressionSimulation(testpars,0)
tail(testdata$data$Introgressed1_avg,n=1)

plot(testdata$data$Major0_avg)
plot(testdata$data$Major1_avg)
plot(testdata$data$Popsize_avg)
## shiny stuff

shiny::runApp("/home/freek/Documents/VisualCode/C++/Introgression/pkgIntrogression/inst/ShinyApp", display.mode = "normal")

