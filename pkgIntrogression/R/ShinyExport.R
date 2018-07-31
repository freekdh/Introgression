
runExample <- function() {
    require(shiny)
    shinyApp(ui = "ui.r", server = "server.r")
}

runExample <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "mypackage")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}