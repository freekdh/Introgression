library(shiny)
library(Rcpp)
library(BH)
library(RcppProgress)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(plotly)

sourceCpp("/home/freek/Documents/VisualCode/C++/Introgression/MainRcpp.cpp")

testdata=RunSimulation( 
r = 0.1,    #growthrate
nloci = 100,    #nloci
nploidy = 2,      #nploidy
ninit0 = 190,      #ninit0
ninit1 = 10,      #ninit1
distlocal = 1,      #distlocal
scmajor = 0.1,    #scmajor
sclocal = 0.0,    #sclocal
ngen = 50,     #ngen
nrep = 100,    #nrep  
rec= 0.01,    #rec
k = 200,        #k
threads = 1)      #threads

#Population size and type0 & type1
plotPopulation <- function(data){
ggplot(data[[2]]) + 
geom_line(aes(x=Generation, y=Popsize_avg)) + geom_ribbon(aes(x=Generation, ymin=Popsize_avg-sqrt(Popsize_var), ymax=Popsize_avg+sqrt(Popsize_var)), fill = "blue", alpha = 0.3) +
geom_line(aes(x=Generation, y=Major0_avg)) + geom_ribbon(aes(x=Generation, ymin=Major0_avg-sqrt(Major0_avg), ymax=Major0_avg+sqrt(Major0_avg)), fill = "red", alpha = 0.3) +
geom_line(aes(x=Generation, y=Major1_avg)) + geom_ribbon(aes(x=Generation, ymin=Major1_avg-sqrt(Major1_avg), ymax=Major1_avg+sqrt(Major1_avg)), fill = "green", alpha = 0.3) +
theme_minimal() + xlab("Generation") + ylab("Populationsize")
}


#Introgression 
plotIntrogression <- function(data){
ggplot(data[[2]]) + 
geom_line(aes(x=Generation, y=Introgressed0_avg)) +
geom_ribbon(aes(x=Generation, ymin=Introgressed0_avg-sqrt(Introgressed0_var), ymax=Introgressed0_avg+sqrt(Introgressed0_var)), fill = "red", alpha = 0.3) +
geom_line(aes(x=Generation, y=Introgressed1_avg)) +
geom_ribbon(aes(x=Generation, ymin=Introgressed1_avg-sqrt(Introgressed1_var), ymax=Introgressed1_avg+sqrt(Introgressed1_var)), fill = "green", alpha = 0.3) +
theme_minimal() + 
xlab("Generation") + 
ylab("Introgression")
}


grid.arrange(plotPopulation(testdata), plotIntrogression(testdata),ncol = 1)

#allelefrequencies
dataallelemean <-tidyr::gather(as.data.frame(testdata[[3]]), locus, value, 2:ncol(as.data.frame(testdata[[3]])))
dataallelevar <-tidyr::gather(as.data.frame(testdata[[4]]), locus, value, 2:ncol(as.data.frame(testdata[[4]])))
dataallele <- dplyr::bind_cols(dataallelemean,dataallelevar[3])
dataallele <- dplyr::mutate(dataallele, col=factor(locus,locus))
allelefrequencyplot <- ggplot(dataallele, aes(col, value, ymin=value-sqrt(value1),ymax=value+sqrt(value1))) + 
geom_pointrange(aes(frame = Generation), shape=22) + 
theme_minimal() + 
xlab("Locus") +
ylab("Allelefrequency") + 
theme(axis.text.x = element_blank())

ggplotly(allelefrequencyplot)


shinyApp(ui,server)

Rcpp.package.skeleton(
    name = "IntrogressionRcpp7", 
    cpp_files = c("MainRcpp.cpp", "random.h", "random.cpp", "utils.h", "utils.cpp"), 
    example_code = FALSE,
    author = "F.J.H. de Haas",
    email = "dehaas@zoology.ubc.ca")
#Add to the DESCRIPTION file:
#LinkingTo: Rcpp, BH, RcppProgress
#Depends: BH, RcppProgress
# Something


