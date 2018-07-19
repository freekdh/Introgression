library(shiny)
library(Rcpp)
library(BH)
library(RcppProgress)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(rbenchmark)
library(gridExtra)

sourceCpp("/home/freek/Documents/VisualCode/C++/Introgression/MainRcpp_Multithread.cpp")

### BASELINE PARAMETERS ###
pars <- list(
"r" = 0.1,          #growthrate
"nloci" = 10,       #nloci
"nploidy" = 2,      #nploidy
"ninit0" = 100,      #ninit0
"ninit1" = 1,       #ninit1
"distlocal" = 1,    #distlocal
"scmajor" = 0.1,    #scmajor
"sclocal" = 0.0,    #sclocal
"ngen" = 50,        #ngen
"rec"= 0.5,         #rec
"k" = 100,           #k
"nrep" = 500,       #nrep
"threads" = 0       #multithreading
)

start_time <- Sys.time()
testdata <- do.call(IntrogressionSimulation, pars)
end_time <- Sys.time()
end_time - start_time

#################### FIXATION PROBABILITY ~ FIGURE 1 ######################

#Collect data
fig1data <- data.frame(fixation = numeric(), ninit1 = numeric(), scmajor = numeric())
localpars <- pars
for(sc in seq(0,0.1,0.02))
for(i in c(1:10)){
    #set pars
    localpars$scmajor <- sc
    localpars$ninit1 <- i
    localpars$ninit0 <- pars$k - i
    #do simulation
    vx <- do.call(IntrogressionSimulation, localpars)$fixation
    #store data:
    fig1data <- rbind(fig1data, data.frame(fixation = vx, ninit1 = localpars$ninit1, scmajor = localpars$scmajor))
}

# on 1 thread: 1.351944 min !!
# on 4 threads: 37.104 secs !!
# on 8 threads: 34.766 secs !!

#Plot data
fig1data$scmajor <- as.factor(fig1data$scmajor)
theme_set(theme_gray(base_size = 20))
FixationPlot <- ggplot(fig1data) + 
geom_line(aes(x=ninit1, y = fixation, colour = scmajor),size = 2) +
xlab("# rescuetype individuals") + ylab("fixation probability") + ylim(0,1)
ggtitle("Figure 1") +
theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks=seq(1,10,3)) + scale_y_continuous(breaks=seq(0.0,1,0.25))

#Save plot
FixationPlot
ggsave("./Figures/Fig1.pdf", width = 20, height = 20, units = "cm")

#################### INTROGRESSION ~ FIGURE 2 ######################

#Collect data
fig2data <- data.frame(introgressed_avg = numeric(), introgressed_var = numeric(), ninit1 = numeric(), scmajor = numeric())
localpars <- pars
for(sc in seq(0,0.1,0.02))
for(i in c(1:10)){
    #set pars
    localpars$ninit1 <- i
    localpars$ninit0 <- pars$k - i
    localpars$scmajor <- sc
    #do simulation
    vx <- do.call(IntrogressionSimulation, localpars)
    introgressed_value_avg <- tail(vx$data$Introgressed1_avg , n=1)
    introgressed_value_var <- tail(vx$data$Introgressed1_var , n=1)
    #store data:
    fig2data <- rbind(fig2data, data.frame(introgressed_avg = introgressed_value_avg, introgressed_var = introgressed_value_var, ninit1 = localpars$ninit1, scmajor = localpars$scmajor))
}

fig2data$introgressed_avg <- 1-fig2data$introgressed_avg
fig2data$scmajor <- as.factor(fig2data$scmajor)

#Plot data
theme_set(theme_gray(base_size = 20))
IntrogressedPlot <- ggplot(fig2data) + 
geom_line(aes(x=ninit1, y = introgressed_avg, colour = scmajor),size = 2) +
#geom_ribbon(alpha = 0.1, aes(x=ninit1, ymin=introgressed_avg-sqrt(introgressed_var), ymax=introgressed_avg+sqrt(introgressed_var),colour = scmajor)) +
xlab("# rescuetype individuals") + ylab("Introgression") + ylim(0,1)
ggtitle("Figure 2") +
theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks=seq(1,10,3)) + scale_y_continuous(breaks=seq(0.0,1,0.1))

#Save plot
IntrogressedPlot
ggsave("./Figures/Fig2.pdf", width = 20, height = 20, units = "cm")

#################### RECOMBINATION ~ FIGURE 3 ######################

#Collect data
fig3data <- data.frame(introgressed_avg = numeric(), fixation = numeric(), scmajor = numeric(), recombination = numeric())
localpars <- pars
for(r in seq(0,0.5,0.01))
for(sc in seq(0,0.2,0.2)){
    #set pars
    localpars$rec <- r
    localpars$scmajor <- sc
    #do simulation
    vx <- do.call(IntrogressionSimulation, localpars)
    introgressed_value_avg <- tail(vx$data$Introgressed1_avg , n=1)
    #store data:
    fig3data <- rbind(fig3data, data.frame(introgressed_avg = introgressed_value_avg, scmajor = localpars$scmajor,recombination = localpars$rec))
}

fig3data
fig3data$introgressed_avg <- 1-fig3data$introgressed_avg
fig3data$scmajor <- as.factor(fig3data$scmajor)

#Plot data
theme_set(theme_gray(base_size = 20))
Fig3Plot <- ggplot(fig3data) + 
geom_line(aes(x=recombination, y = introgressed_avg, colour = scmajor),size = 2) +
#geom_ribbon(alpha = 0.1, aes(x=ninit1, ymin=introgressed_avg-sqrt(introgressed_var), ymax=introgressed_avg+sqrt(introgressed_var),colour = scmajor)) +
xlab("recombination rate") + ylab("Introgression") + 
ggtitle("Figure 3") +
theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks=seq(0,0.5,0.1)) + scale_y_continuous(breaks=seq(0.0,1,0.1))

#Save plot
Fig3Plot
ggsave("./Figures/Fig3.pdf", width = 20, height = 20, units = "cm")

grid.arrange(FixationPlot, IntrogressedPlot, Fig3Plot)

#################### REST ######################


#Population size and type0 & type1
ggplot(data[[2]]) + 
geom_line(aes(x=Generation, y=Popsize_avg)) + 
geom_ribbon(aes(x=Generation, ymin=Popsize_avg-sqrt(Popsize_var), ymax=Popsize_avg+sqrt(Popsize_var)), fill = "blue", alpha = 0.3) +
geom_line(aes(x=Generation, y=Major0_avg)) + 
geom_ribbon(aes(x=Generation, ymin=Major0_avg-sqrt(Major0_avg), ymax=Major0_avg+sqrt(Major0_avg)), fill = "red", alpha = 0.3) +
geom_line(aes(x=Generation, y=Major1_avg)) + 
geom_ribbon(aes(x=Generation, ymin=Major1_avg-sqrt(Major1_avg), ymax=Major1_avg+sqrt(Major1_avg)), fill = "green", alpha = 0.3) +
theme_minimal() + xlab("Generation") + ylab("Populationsize")

#Introgression 
ggplot(data[[2]]) + 
geom_line(aes(x=Generation, y=Introgressed0_avg)) +
geom_ribbon(aes(x=Generation, ymin=Introgressed0_avg-sqrt(Introgressed0_var), ymax=Introgressed0_avg+sqrt(Introgressed0_var)), fill = "red", alpha = 0.3) +
geom_line(aes(x=Generation, y=Introgressed1_avg)) +
geom_ribbon(aes(x=Generation, ymin=Introgressed1_avg-sqrt(Introgressed1_var), ymax=Introgressed1_avg+sqrt(Introgressed1_var)), fill = "green", alpha = 0.3) +
theme_minimal() + 
xlab("Generation") + 
ylab("Introgression")

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
