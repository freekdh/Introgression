
paramNames <- c("r", "nloci", "nploidy", "ninit0" ,"ninit1", 
	"distlocal", "scmajor", "sclocal", "ngen", "nrep", "rec", "k", "threads")

plot_nav <- function(nav) {

  layout(matrix(c(1,2,1,3),2,2))

  palette(c("black", "grey50", "grey30", "grey70", "#d9230f"))

  # plot all scenarios
  matplot(nav,
    type = 'l', lwd = 0.5, lty = 1, col = 1:5,
    xlab = 'Months', ylab = 'Millions',
    main = 'Projected Value of Initial Capital')

  # plot % of scenarios that are still paying
  p.alive = 1 - rowSums(is.na(nav)) / ncol(nav)

  plot(100 * p.alive, las = 1, xlab = 'Months', ylab = 'Percentage Paying',
    main = 'Percentage of Paying Scenarios', ylim=c(0,100))
  grid()


  last.period = nrow(nav)

  # plot distribution of final wealth
  final.nav = nav[last.period, ]
  final.nav = final.nav[!is.na(final.nav)]

  if(length(final.nav) ==  0) return()

  plot(density(final.nav, from=0, to=max(final.nav)), las = 1, xlab = 'Final Capital',
    main = paste0('Distribution of Final Capital\n', 100 * p.alive[last.period], '% are still paying'))
  grid()
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
  navA <- reactive(do.call(RunSimulation, getParams("a")))
  navB <- reactive(do.call(RunSimulation, getParams("b")))

  # Expression that plot NAV paths. The expression
  # is wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically
  #     re-executed when inputs change
  #  2) Its output type is a plot
  #
  output$a_distPlot <- renderPlot({
  	plot_nav(navA())
  })
  output$b_distPlot <- renderPlot({
  	plot_nav(navB())
  })

}