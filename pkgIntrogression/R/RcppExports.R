# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

RcppIntrogressionSimulation <- function(parslist, setthreads = 0L, progressbar = FALSE) {
    .Call('_pkgIntrogression_RcppIntrogressionSimulation', PACKAGE = 'pkgIntrogression', parslist, setthreads, progressbar)
}

ShinyInitializeSimulation <- function(parslist) {
    invisible(.Call('_pkgIntrogression_ShinyInitializeSimulation', PACKAGE = 'pkgIntrogression', parslist))
}

ShinyRunSimulation <- function() {
    invisible(.Call('_pkgIntrogression_ShinyRunSimulation', PACKAGE = 'pkgIntrogression'))
}

ShinyWriteOutputandCleanup <- function() {
    .Call('_pkgIntrogression_ShinyWriteOutputandCleanup', PACKAGE = 'pkgIntrogression')
}

