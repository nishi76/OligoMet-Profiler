# =============================================================================
# app.R -- OligoMet Profiler (repository-root launcher)
# 
# Author: Nishikant Wase, PhD 
# email: nishikant.wase@gmail.com
# Research Scientist
# Thermofisher Scientific
# 
# 
# 
# The dashboard itself lives in inst/app/app.R so that it ships with the
# installed package; installed users launch it with
# OligoMetProfiler::run_app(). This file keeps the repository root working
# as a Shiny app directory, for shiny::runApp("."), shiny::runGitHub(), and
# RStudio's "Run App" button.
# =============================================================================

source(file.path("inst", "app", "app.R"))$value
