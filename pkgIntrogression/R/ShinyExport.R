
IntrogressionShiny <- function() {
  appDir <- system.file("ShinyApp", package = "pkgIntrogression")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
