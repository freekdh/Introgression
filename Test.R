library(pkgIntrogression)
testpars <- list(
    r=0.1,
    nloci=10,
    nploidy=2,
    ninit0=10,
    ninit1=10,
    distlocal=1,
    scmajor=0.0,
    sclocal=0.0,
    ngen=100,
    nrep=1000,
    rec=0.5,
    k=100)

ShinyInitializeSimulation(testpars)
ShinyRunSimulation()
ShinyWriteOutputandCleanup()

hello <- RcppIntrogressionSimulation(testpars,4)

names(hello)
plot(hello$data$Popsize_avg)
## shiny stuff

shiny::runApp("/home/freek/Documents/VisualCode/C++/Introgression/pkgIntrogression/ShinyApp", display.mode = "normal")

