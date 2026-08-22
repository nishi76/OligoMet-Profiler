# =============================================================================
# run_app.R -- launcher for the bundled Shiny dashboard
# =============================================================================

# Launch the OligoMet Profiler Shiny dashboard from the installed package.
# Extra arguments are passed to shiny::runApp() (launch.browser, port, host).
run_app <- function(...) {
  missing <- c("shiny", "DT")[!vapply(c("shiny", "DT"), requireNamespace,
                                      logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("The dashboard needs the following package(s): ",
         paste(missing, collapse = ", "), ". Install with install.packages(c(",
         paste(sprintf('"%s"', missing), collapse = ", "), ")).",
         call. = FALSE)
  }

  app_dir <- system.file("app", package = "OligoMetProfiler")
  if (!nzchar(app_dir)) {
    stop("Could not find the app directory in the installed package. ",
         "Reinstall with remotes::install_github(\"nishi76/OligoMet-Profiler\").",
         call. = FALSE)
  }

  shiny::runApp(app_dir, ...)
}
