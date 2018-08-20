library(pkgIntrogression)
library(ggplot2)
library(assertthat)

paramNames <- c("r", "nloci", "nploidy", "ninit0" ,"ninit1", 
	"distlocal", "scmajor", "sclocal", "ngen", "nrep", "rec", "k"
)

plotfixation <- function(pars, xvar, xvec, groupvar, groupvec) {
    assert_that(length(xvec)>2)
    assert_that(length(groupvec)>2)
    assert_that(is.element(xvar, paramNames))
    assert_that(is.element(groupvar, paramNames))

    #Collect data
    fixdata <- data.frame(xvar = numeric(), fixation = numeric(), groupvar = numeric())
    for(x in xvec)
    for(group in groupvec){
        #set pars
        pars$xvar <- x
        pars$groupvar <- group
        #do simulation
        data <- RcppIntrogressionSimulation(pars)$fixation
        #store data:
        fixdata <- rbind(fixdata, data.frame(xvar = x, fixation =  data, groupvar = group))
    }

    #Plot data
    fixdata$scmajor <- as.factor(fixdata$groupvar)
    theme_set(theme_gray(base_size = 20))
    FixationPlot <- ggplot(fixdata) + 
    geom_line(aes(x=xvar, y = fixation, colour = scmajor),size = 2) +
    xlab(xvar) + ylab("fixation probability") + ylim(0,1)
    ggtitle("Figure 1")

    #Save plot
    FixationPlot
#    ggsave("./Figures/Fig1.pdf", width = 20, height = 20, units = "cm")
}
#

#plotfixation(testpars, "nloci", c(10,20,30), "ninit0" , c(0.01,0.02,0.03,0.04,0.05,2))
