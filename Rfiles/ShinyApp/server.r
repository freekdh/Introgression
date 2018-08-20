
paramNames <- c("r", "nloci", "nploidy", "ninit0" ,"ninit1", 
	"distlocal", "scmajor", "sclocal", "ngen", "nrep", "rec", "k"
)

plot_pop <- function(nav) {
  ggplot(nav[[2]]) + 
  geom_line(aes(x=Generation, y=Popsize_avg)) + geom_ribbon(aes(x=Generation, ymin=Popsize_avg-sqrt(Popsize_var), ymax=Popsize_avg+sqrt(Popsize_var)), fill = "blue", alpha = 0.3) +
  geom_line(aes(x=Generation, y=Major0_avg)) + geom_ribbon(aes(x=Generation, ymin=Major0_avg-sqrt(Major0_avg), ymax=Major0_avg+sqrt(Major0_avg)), fill = "red", alpha = 0.3) +
  geom_line(aes(x=Generation, y=Major1_avg)) + geom_ribbon(aes(x=Generation, ymin=Major1_avg-sqrt(Major1_avg), ymax=Major1_avg+sqrt(Major1_avg)), fill = "green", alpha = 0.3) +
  theme_minimal() + xlab("Generation") + ylab("Populationsize")
}

plot_intro <- function(nav) {
  ggplot(nav[[2]]) + 
  geom_line(aes(x=Generation, y=Introgressed0_avg)) +
  geom_ribbon(aes(x=Generation, ymin=Introgressed0_avg-sqrt(Introgressed0_var), ymax=Introgressed0_avg+sqrt(Introgressed0_var)), fill = "red", alpha = 0.3) +
  geom_line(aes(x=Generation, y=Introgressed1_avg)) +
  geom_ribbon(aes(x=Generation, ymin=Introgressed1_avg-sqrt(Introgressed1_var), ymax=Introgressed1_avg+sqrt(Introgressed1_var)), fill = "green", alpha = 0.3) +
  theme_minimal() + 
  xlab("Generation") + 
  ylab("Introgression")
}

RunAllSimulation <- function(pars){
  pkgIntrogression::ShinyInitializeSimulation(pars)
  
  withProgress(message="Running simulations",value = 0,{
    for(x in c(1:pars$nrep)){
      incProgress(1/pars$nrep, detail = paste("Doing replicate", x))
      pkgIntrogression::ShinyRunSimulation()
    }
  })

  pkgIntrogression::ShinyWriteOutputandCleanup()
}

function(input, output, session) {

  getParams <- function(prefix) {
    input[[paste0(prefix, "_recalc")]]

    params <- list()
    length(params) <- 12
    params <- lapply(paramNames, function(p) {input[[paste0(prefix, "_", p)]]})
    names(params) <- paramNames
    as.list(params)
  }

  navA <- eventReactive(input$a_recalc, RunAllSimulation(getParams("a")) ) # 

  output$a_PopulationPlot <- renderPlot({
    suppressWarnings(plot_pop(navA()))
  })

  output$a_IntrogressionPlot <- renderPlot({
    suppressWarnings(plot_intro(navA()))
  })

  output$fixation_output <- renderText({ 
    paste("The fixation probability is: ", navA()$fixation)
  })

  output$a_lociPlot <- renderPlotly({
  #allelefrequencies
  dataallelemean <-tidyr::gather(as.data.frame(navA()[[3]]), locus, value, 2:ncol(as.data.frame(navA()[[3]])))
  dataallelevar <-tidyr::gather(as.data.frame(navA()[[4]]), locus, value, 2:ncol(as.data.frame(navA()[[4]])))
  dataallele <- dplyr::bind_cols(dataallelemean,dataallelevar[3])
  dataallele <- dplyr::mutate(dataallele, col=factor(locus,locus))
  allelefrequencyplot <- ggplot(dataallele, aes(col, value, ymin=value-sqrt(value1),ymax=value+sqrt(value1))) + 
  geom_pointrange(aes(frame = Generation), shape=22) + 
  theme_minimal() + 
  xlab("Locus") +
  ylab("Allelefrequency") + 
  theme(axis.text.x = element_blank())

  suppressWarnings(ggplotly(allelefrequencyplot))
})

}