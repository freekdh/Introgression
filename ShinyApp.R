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
0.1,    #growthrate
100,     #nloci
2,      #nploidy
190,      #ninit0
10,      #ninit1
1,      #distlocal
0.1,    #scmajor
0.0,    #sclocal
50,     #ngen
100,    #nrep  
0.01,    #rec
200,    #k
4)      #threads

#Population size and type0 & type1
populationsizeplot <- ggplot(testdata[[2]]) + 
geom_line(aes(x=Generation, y=Popsize_avg)) + 
geom_ribbon(aes(x=Generation, ymin=Popsize_avg-sqrt(Popsize_var), ymax=Popsize_avg+sqrt(Popsize_var)), fill = "blue", alpha = 0.3) +
geom_line(aes(x=Generation, y=Major0_avg)) + 
geom_ribbon(aes(x=Generation, ymin=Major0_avg-sqrt(Major0_avg), ymax=Major0_avg+sqrt(Major0_avg)), fill = "red", alpha = 0.3) +
geom_line(aes(x=Generation, y=Major1_avg)) +
geom_ribbon(aes(x=Generation, ymin=Major1_avg-sqrt(Major1_avg), ymax=Major1_avg+sqrt(Major1_avg)), fill = "green", alpha = 0.3) +
theme_minimal() + 
xlab("Generation") + 
ylab("Populationsize")

#Introgression 
introgressionplot <- ggplot(testdata[[2]]) + 
geom_line(aes(x=Generation, y=Introgressed0_avg)) +
geom_ribbon(aes(x=Generation, ymin=Introgressed0_avg-sqrt(Introgressed0_var), ymax=Introgressed0_avg+sqrt(Introgressed0_var)), fill = "red", alpha = 0.3) +
geom_line(aes(x=Generation, y=Introgressed1_avg)) +
geom_ribbon(aes(x=Generation, ymin=Introgressed1_avg-sqrt(Introgressed1_var), ymax=Introgressed1_avg+sqrt(Introgressed1_var)), fill = "green", alpha = 0.3) +
theme_minimal() + 
xlab("Generation") + 
ylab("Introgression")

grid.arrange(populationsizeplot, introgressionplot,ncol = 1)

#allelefrequencies
dataallelemean <-tidyr::gather(as.data.frame(testdata[[3]]), locus, value, 2:ncol(as.data.frame(testdata[[3]])))
dataallelevar <-tidyr::gather(as.data.frame(testdata[[4]]), locus, value, 2:ncol(as.data.frame(testdata[[4]])))
dataallele <- dplyr::bind_cols(dataallelemean,dataallelevar[3])
dataallele <- dplyr::mutate(dataallele, col=factor(locus,locus))
allelefrequencyplot <- ggplot(dataallele, aes(col, value, ymin=value-sqrt(value1),ymax=value+sqrt(value1))) + 
geom_pointrange(aes(frame = Generation)) + 
theme_minimal() + 
xlab("Locus") +
ylab("Allelefrequency") + 
theme(axis.text.x = element_blank())

ggplotly(allelefrequencyplot)



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
