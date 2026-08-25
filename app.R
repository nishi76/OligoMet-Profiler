# =============================================================================
# app.R -- OligoMet Profiler (repository-root launcher)
#
# Author and developer: Nishikant Wase, PhD <nishikant.wase@gmail.com>
# Research Scientist, Thermo Fisher Scientific. An independent personal
# project, not a Thermo Fisher Scientific product.
#
# FOR RESEARCH USE ONLY. Not for diagnostic, clinical, or regulatory
# submission use. Provided without warranty; the author accepts no
# liability for its use. See DISCLAIMER.md, or run oligomet_about().
#
# The dashboard itself lives in inst/app/app.R so that it ships with the
# installed package; installed users launch it with
# OligoMetProfiler::run_app(). This file keeps the repository root working
# as a Shiny app directory, for shiny::runApp("."), shiny::runGitHub(), and
# RStudio's "Run App" button.
# =============================================================================

source(file.path("inst", "app", "app.R"))$value
