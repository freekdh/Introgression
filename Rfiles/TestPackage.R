library(pkgIntrogressionv2)
testpars <- list(
    b=1.1,
    dA=1,
    da=0.8,
    r=0.5,
    nloci=100,
    ninit0=10,
    ninit1=2,
    ngen=20,
    nrep=1000,
    rec=0.5,
    k=100)

ShinyInitializeSimulation(testpars)
ShinyRunSimulation()
ShinyWriteOutputandCleanup()

hello <- RcppIntrogressionSimulation(testpars,0)
names(hello)
hello$pars

plot(hello$data$Popsize_avg)
## shiny stuff

shiny::runApp("/home/freek/Documents/VisualCode/C++/Introgression/pkgIntrogression/inst/ShinyApp", display.mode = "normal")

