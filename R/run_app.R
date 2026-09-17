# =============================================================================
# run_app.R -- launcher for the bundled Shiny dashboard
# =============================================================================

#' Launch the OligoMet Profiler Shiny dashboard
#'
#' Starts the bundled Shiny app for interactive sequence entry, library
#' generation, MS matching, and batch processing. For headless/scripted
#' use with no Shiny dependency at all, see `run_custom_oligo.R`
#' (single sequence/file) or `run_batch_ms.R` (multi-file batch) in the
#' package/repository root instead -- copy one, edit its CONFIG block,
#' and `Rscript` it.
#'
#' @param ... Passed straight through to `shiny::runApp()` (e.g.
#'   `launch.browser`, `port`, `host`).
#' @return Does not return under normal use -- `shiny::runApp()` blocks
#'   until the app is stopped.
#' @export
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
