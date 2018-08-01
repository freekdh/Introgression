library(shiny)
library(plotly)

renderInputs <- function(prefix) {
  wellPanel(
    fluidRow(
        column(3,
        sliderInput(paste0(prefix, "_", "ngen"), "Number of generations:", min = 10, max = 100, value = 100,step = 20),
        sliderInput(paste0(prefix, "_", "nloci"), "Number of loci:", min = 50, max = 200, value = 50, step = 25),
        sliderInput(paste0(prefix, "_", "nrep"), "Number of repetitions:", min = 1, max = 100, value = 50, step = 10)
        ),
        column(3,
        sliderInput(paste0(prefix, "_", "nploidy"), "Number of ploidy:", min = 2, max = 6, value = 2, step = 2),
        sliderInput(paste0(prefix, "_", "ninit0"), "Number of wildtype individuals:", min = 10, max = 100, value = 50, step = 10),
        sliderInput(paste0(prefix, "_", "ninit1"), "Number of rescuetype individuals:", min = 1, max = 100, value = 5)
        ),
        column(3,
        sliderInput(paste0(prefix, "_", "k"), "Carrying capacity of the population:", min = 10, max = 150, value = 100, step = 20),
        sliderInput(paste0(prefix, "_", "r"), "Intrinsic growthrate of the population:", min = 0.0, max = 0.5, value = 0.1, step = 0.05),
        sliderInput(paste0(prefix, "_", "distlocal"), "Distance of locally adapted locus from major locus:", min = 1, max = 10, value = 1)
        ),
        column(3,
        sliderInput(paste0(prefix, "_", "scmajor"), "Major locus additive fitness:", min = 0.0, max = 0.5, value = 0.05, step = 0.05),
        sliderInput(paste0(prefix, "_", "sclocal"), "Locally adapted locus additive fitness:", min = -0.1, max = 0.0, value = -0.1, step = 0.01),
        sliderInput(paste0(prefix, "_", "rec"), "Recombination rate:", min = 0.0, max = 0.5, value = 0.5, step = 0.05)
        )
    ),
    p(actionButton(paste0(prefix, "_", "recalc"),
      "Run simulation", icon("random")
    ))
  )
}

# Define UI for application that plots random distributions
fluidPage(theme="simplex.min.css",
  tags$style(type="text/css",
    "label {font-size: 12px;}",
    ".recalculating {opacity: 1.0;}"
  ),

  # Application title
  tags$h2("Evolutinar rescue due to introgressive hybridization"),
  p("Simulations implemented in c++"),
  hr(),

  fluidRow(
    column(6,
    renderInputs("a")
    ),
    column(3,
    plotOutput("a_PopulationPlot")
    ),
    column(3,
    plotOutput("a_IntrogressionPlot")
    )
  ),
  fluidRow(
    textOutput("fixation_output")
  ),
  fluidRow(
      plotlyOutput("a_lociPlot")
  )
)