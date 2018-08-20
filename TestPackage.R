library(pkgIntrogression)
testpars <- list(
    r=0.5,
    nloci=100,
    nploidy=2,
    ninit0=10,
    ninit1=2,
    distlocal=1,
    scmajor=0.22,
    sclocal=0.0,
    ngen=20,
    nrep=1000,
    rec=0.5,
    k=100)

ShinyInitializeSimulation(testpars)
ShinyRunSimulation()
ShinyWriteOutputandCleanup()


hello <- RcppIntrogressionSimulation(testpars,0)
names(hello)
hello$data

plot(hello$data$Popsize_avg)
## shiny stuff

shiny::runApp("/home/freek/Documents/VisualCode/C++/Introgression/pkgIntrogression/inst/ShinyApp", display.mode = "normal")

