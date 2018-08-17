library(pkgIntrogression)
testpars <- list(
    r=0.1,
    nloci=10,
    nploidy=2,
    ninit0=50,
    ninit1=2,
    distlocal=1,
    scmajor=0.0,
    sclocal=0.0,
    ngen=50,
    nrep=100,
    rec=0.5,
    k=50)

ShinyInitializeSimulation(testpars)
ShinyRunSimulation()
ShinyWriteOutputandCleanup()


hello <- RcppIntrogressionSimulation(testpars,0)

names(hello)
plot(hello$data$Popsize_avg)
## shiny stuff

shiny::runApp("/home/freek/Documents/VisualCode/C++/Introgression/pkgIntrogression/inst/ShinyApp", display.mode = "normal")

