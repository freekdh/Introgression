
paramNames <- c("r", "nloci", "nploidy", "ninit0" ,"ninit1", 
	"distlocal", "scmajor", "sclocal", "ngen", "nrep", "rec", "k")

sourceCpp("/home/freek/Documents/VisualCode/C++/Introgression/MainRcpp.cpp")

plot_pop <- function(nav) {
  plotPopulation(nav)
}

plot_intro <- function(nav) {
  plotIntrogression(nav)
}

# Define server logic required to generate and plot a random distribution
#
# Idea and original code by Pierre Chretien
# Small updates by Michael Kapler
#
function(input, output, session) {

	getParams <- function(prefix) {
		input[[paste0(prefix, "_recalc")]]

		params <- lapply(paramNames, function(p) {
			input[[paste0(prefix, "_", p)]]
		})
		names(params) <- paramNames
		params
	}

  # Function that generates scenarios and computes NAV. The expression
  # is wrapped in a call to reactive to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically
  #     re-executed when inputs change
  #
  navA <- reactive(do.call(RunSimulation, getParams("a"))) # 

  # Expression that plot NAV paths. The expression
  # is wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically
  #     re-executed when inputs change
  #  2) Its output type is a plot
  #

output$a_PopulationPlot <- renderPlot({
    plot_pop(navA())
  })

output$a_IntrogressionPlot <- renderPlot({
    plot_intro(navA())
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

  ggplotly(allelefrequencyplot)
})

}