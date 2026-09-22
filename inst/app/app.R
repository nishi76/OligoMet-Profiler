# =============================================================================
# app.R -- OligoMet Profiler
# Shiny app wrapper for the oligonucleotide metabolite identification
# pipeline. Runs on any oligonucleotide -- standard chemistry works out of
# the box, and the custom chemistry table below covers anything else.
#
# Run with:
#   OligoMetProfiler::run_app()          # from the installed package
#   shiny::runApp(".", launch.browser=TRUE)   # from a repository checkout
#   or in RStudio: open app.R in the repository root and click "Run App"
#
# The app takes its pipeline functions from the installed package, or from
# R/ when run inside a checkout, and provides an interactive interface for:
#   - Sequence input (triplet or OligoDistiller notation, auto-detected)
#   - Custom chemistry overrides (editable table)
#   - Full parameter control (truncation depth, charge range, oxidation, etc.)
#   - Optional MS data upload for matching
#   - Summary dashboard with 4 plots + download buttons
# =============================================================================

## ---- Module sourcing -------------------------------------------------------
# The app runs in two situations:
#   1. From an installed OligoMetProfiler package (OligoMetProfiler::run_app()),
#      where the pipeline functions already live in the package namespace.
#   2. From a repository checkout, where they must be source()d from R/.
# A checkout wins when one is found, so edits take effect without reinstalling.
.find_module_dir <- function() {
  candidates <- c(getwd(), file.path(getwd(), "..", ".."))
  file_arg <- sub("^--file=", "",
                  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) > 0) {
    here <- dirname(normalizePath(file_arg[1], mustWork = FALSE))
    candidates <- c(candidates, here, file.path(here, "..", ".."))
  }
  for (d in candidates) {
    if (file.exists(file.path(d, "R", "chemistry_dict.R")))
      return(normalizePath(d))
  }
  NULL
}

.module_dir <- .find_module_dir()

if (!is.null(.module_dir)) {
  for (.f in c("about.R", "default_params.R", "progress_utils.R",
               "chemistry_dict.R", "oligo_io.R",
               "metabolites.R", "mass_isotope.R", "fragments.R",
               "ms_matching.R", "spectra_io.R", "batch_ms_processing.R", "statistics.R",
               "degradation.R", "multivariate.R", "build_workbook.R", "build_report.R",
               "export_acquisition.R", "export_spectral.R", "mirror_plot.R",
               "agent_tools.R", "agent_core.R")) {
    source(file.path(.module_dir, "R", .f))
  }
} else if (requireNamespace("OligoMetProfiler", quietly = TRUE)) {
  library(OligoMetProfiler)
} else {
  stop("Cannot locate the pipeline modules. Either install the package -- ",
       "remotes::install_github(\"nishi76/OligoMet-Profiler\"), then launch ",
       "with OligoMetProfiler::run_app() -- or run shiny::runApp() on a ",
       "repository checkout.")
}

if (!requireNamespace("DT", quietly = TRUE)) {
  stop("The 'DT' package is required by the dashboard. ",
       "Install it with install.packages(\"DT\").")
}

# Optional: without shinyFiles the "Save to folder" box still works, it just
# loses its folder-browser button.
.have_shinyfiles <- requireNamespace("shinyFiles", quietly = TRUE)

library(shiny)
library(DT)

# Shiny's default per-request upload cap is 5MB -- fine for the single-file
# "MS Matching" section, but batch MS files (esp. HRMS) routinely run from
# tens of MB to multiple GB, so the unmodified default silently rejects
# them before batch_files' server logic ever runs. Raised generously here;
# a hosted deployment (shinyapps.io/Posit Connect) may still enforce its
# own front-end request-size ceiling independent of this option.
options(shiny.maxRequestSize = 20 * 1024^3)  # 20 GB

# Best-effort hosted-deployment detection, for advisory messaging only --
# Posit Connect sets CONNECT_SERVER for content it runs; shinyapps.io and
# other managed Shiny hosts set SHINY_PORT (not exclusive to them, but a
# reasonable signal in combination with the absence of an interactive
# session). Neither is a guarantee, and there is no way to read the actual
# admin-configured upload cap from inside the app -- it's enforced by the
# reverse proxy in front of R, independent of the option above. This only
# gates a warning telling the user where to look if uploads fail; it never
# blocks anything itself.
.is_hosted_deployment <- function() {
  nzchar(Sys.getenv("CONNECT_SERVER")) || nzchar(Sys.getenv("SHINY_PORT"))
}

## ---- Example sequences (illustrative, not the assumed input) --------------
# A generic gapmer built entirely from standard dictionary codes, plus four
# approved oligonucleotide therapeutics -- one per modality class in
# Takakusa et al. (2023) Table 1. The two antisense examples reproduce
# their published molecular formulas exactly (see validate_reference() in
# chemistry_dict.R). None is required; paste any sequence into the box.
# Duplex drugs are given as their sense strand -- run each strand
# separately.
.EXAMPLE_SEQS <- c(
  list("Generic 2'MOE/DNA gapmer" = list(
    seq = "Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm",
    conj5 = "none", conj3 = "none")),
  setNames(
    lapply(REFERENCE_OLIGOS, function(r)
      list(seq = r$triplet, conj5 = r$conj5, conj3 = r$conj3)),
    vapply(REFERENCE_OLIGOS,
           function(r) sprintf("%s - %s", r$modality, r$name), ""))
)

## ---- Conjugate options -----------------------------------------------------
.fatty_acid_choices <- c("myristoyl", "palmitoyl", "stearoyl", "docosanoyl")
.conj5_choices <- c("none", "5'-phosphate", "5'-thiophosphate", "biotin",
                    "cAG_cap", "cAU_cap", "ARCA_cap", "mCAP",
                    "GalNAc", "GalNAc3",
                    "cholesterol", "C6", "C12", "TEG", "FAM", "Cy3",
                    .fatty_acid_choices)
.conj3_choices <- c("none", "3'-phosphate", "3'-cyclophos", "3'-thiophosphate",
                    "GalNAc", "GalNAc3", "GalNAc3_triantennary",
                    "cholesterol", "C6", "C12", "TEG", "FAM", "Cy3",
                    .fatty_acid_choices)

## ---- Custom chemistry table initial data -----------------------------------
.custom_chem_init <- data.frame(
  Code    = c("", "", ""),
  Formula = c("", "", ""),
  Name    = c("", "", ""),
  Type    = c("base", "sugar", "linkage"),
  Attach  = c("add", "add", "add"),
  stringsAsFactors = FALSE
)

## ---- Custom chemistry table -> build_dictionary() overrides ----------------
# Shared by the Run handler and by manual sequence entry, so a code typed
# into the Custom Chemistry table validates the same way in both.
# Rows with an invalid formula are skipped here (never reach
# build_dictionary()) -- see .chem_row_status()/.chem_table_errors() below,
# which are what block Run and surface the problem to the user instead of
# it silently producing a zero-mass residue.
.overrides_from_table <- function(cd) {
  overrides <- list()
  for (i in seq_len(nrow(cd))) {
    code <- trimws(cd$Code[i])
    formula_str <- trimws(cd$Formula[i])
    if (nchar(code) == 0 || nchar(formula_str) == 0) next
    if (!is_valid_formula_string(formula_str)) next
    entry <- list(
      formula = formula_str,
      name = if (nchar(trimws(cd$Name[i])) > 0) trimws(cd$Name[i]) else code
    )
    if (cd$Type[i] == "conjugate") entry$attach <- cd$Attach[i]
    entry$kind <- cd$Type[i]
    overrides[[code]] <- entry
  }
  overrides
}

## ---- Custom chemistry table validation -------------------------------------
# A blank row is intentionally ignored (see .overrides_from_table() above);
# any other row must have a code AND a formula that parses to a valid
# elemental composition, or it is flagged here rather than silently
# reaching build_dictionary() as an all-zero-mass entry.
.chem_row_status <- function(cd) {
  vapply(seq_len(nrow(cd)), function(i) {
    code <- trimws(cd$Code[i]); formula_str <- trimws(cd$Formula[i])
    if (!nzchar(code) && !nzchar(formula_str)) return("")
    if (!nzchar(code)) return("✗ missing code")
    if (!nzchar(formula_str)) return("✗ missing formula")
    if (!is_valid_formula_string(formula_str)) return("✗ invalid formula")
    "✓ valid"
  }, character(1))
}
.chem_table_errors <- function(cd) {
  st <- .chem_row_status(cd)
  bad <- which(startsWith(st, "✗"))
  if (length(bad) == 0) return(character(0))
  sprintf("row %d (code '%s'): %s", bad, trimws(cd$Code[bad]), st[bad])
}

## ---- Small inline help icon -------------------------------------------------
# A native title-attribute tooltip -- no JS dependency, works everywhere --
# for surfacing the tradeoffs behind a parameter (e.g. "20% NCE is a
# starting point, not validated") right next to its input instead of only
# in the README.
.info_icon <- function(text) {
  tags$span("ⓘ", title = text,
            style = "cursor: help; color: #6c757d; margin-left: 4px; font-size: 12px;")
}
.with_info <- function(label, ...) tagList(label, .info_icon(paste0(...)))

## ---- Session save/load ------------------------------------------------------
# "What exactly produced this inclusion list" is a reproducibility question
# for a tool whose outputs feed instrument methods -- a saved session is the
# sequence, every parameter below, and the Custom Chemistry table, as one
# restorable .json file. Uploaded MS/batch files are deliberately excluded
# (too large to round-trip through a small JSON, and re-uploading is the
# normal Shiny flow anyway); the file input's own name is not restorable.
.session_input_ids <- c(
  "seq", "oligo_name", "output_prefix", "output_dir", "conj5", "conj3",
  "max_3p", "max_5p", "endo", "endo_sites", "min_frag_len",
  "z_min", "z_max", "n_iso", "max_oxid", "h_offset", "use_envipat",
  "method_length", "ms2_z_min", "ms2_z_max", "hcd_nce",
  "ms1_target_cap", "ms2_target_cap",
  "enable_ms", "ppm_tol", "noise_mode", "min_intensity", "sn_threshold",
  "adducts", "frag_tol_ppm", "frag_z_max",
  "enable_batch", "batch_run_ms2", "batch_n_workers", "batch_deconv_ppm",
  "batch_noise_mode", "batch_min_intensity", "batch_sn_threshold", "batch_dir",
  "man_bases", "man_sugars", "man_linkages"
)
# One update*Input() call per id above -- dispatch table so loading a
# session doesn't need a long if/else chain keyed on input type.
.session_update <- list(
  seq            = function(s, v) updateTextAreaInput(s, "seq", value = v),
  oligo_name     = function(s, v) updateTextInput(s, "oligo_name", value = v),
  output_prefix  = function(s, v) updateTextInput(s, "output_prefix", value = v),
  output_dir     = function(s, v) updateTextInput(s, "output_dir", value = v),
  conj5          = function(s, v) updateSelectInput(s, "conj5", selected = v),
  conj3          = function(s, v) updateSelectInput(s, "conj3", selected = v),
  max_3p         = function(s, v) updateNumericInput(s, "max_3p", value = v),
  max_5p         = function(s, v) updateNumericInput(s, "max_5p", value = v),
  endo           = function(s, v) updateCheckboxInput(s, "endo", value = v),
  endo_sites     = function(s, v) updateRadioButtons(s, "endo_sites", selected = v),
  min_frag_len   = function(s, v) updateNumericInput(s, "min_frag_len", value = v),
  z_min          = function(s, v) updateNumericInput(s, "z_min", value = v),
  z_max          = function(s, v) updateNumericInput(s, "z_max", value = v),
  n_iso          = function(s, v) updateNumericInput(s, "n_iso", value = v),
  max_oxid       = function(s, v) updateNumericInput(s, "max_oxid", value = v),
  h_offset       = function(s, v) updateNumericInput(s, "h_offset", value = v),
  use_envipat    = function(s, v) updateCheckboxInput(s, "use_envipat", value = v),
  method_length  = function(s, v) updateNumericInput(s, "method_length", value = v),
  ms2_z_min      = function(s, v) updateNumericInput(s, "ms2_z_min", value = v),
  ms2_z_max      = function(s, v) updateNumericInput(s, "ms2_z_max", value = v),
  hcd_nce        = function(s, v) updateNumericInput(s, "hcd_nce", value = v),
  ms1_target_cap = function(s, v) updateNumericInput(s, "ms1_target_cap", value = v),
  ms2_target_cap = function(s, v) updateNumericInput(s, "ms2_target_cap", value = v),
  enable_ms      = function(s, v) updateCheckboxInput(s, "enable_ms", value = v),
  ppm_tol        = function(s, v) updateNumericInput(s, "ppm_tol", value = v),
  noise_mode     = function(s, v) updateRadioButtons(s, "noise_mode", selected = v),
  min_intensity  = function(s, v) updateNumericInput(s, "min_intensity", value = v),
  sn_threshold   = function(s, v) updateNumericInput(s, "sn_threshold", value = v),
  adducts        = function(s, v) updateCheckboxGroupInput(s, "adducts", selected = v),
  frag_tol_ppm   = function(s, v) updateNumericInput(s, "frag_tol_ppm", value = v),
  frag_z_max     = function(s, v) updateNumericInput(s, "frag_z_max", value = v),
  enable_batch   = function(s, v) updateCheckboxInput(s, "enable_batch", value = v),
  batch_run_ms2  = function(s, v) updateCheckboxInput(s, "batch_run_ms2", value = v),
  batch_n_workers  = function(s, v) updateNumericInput(s, "batch_n_workers", value = v),
  batch_deconv_ppm = function(s, v) updateNumericInput(s, "batch_deconv_ppm", value = v),
  batch_noise_mode = function(s, v) updateRadioButtons(s, "batch_noise_mode", selected = v),
  batch_min_intensity = function(s, v) updateNumericInput(s, "batch_min_intensity", value = v),
  batch_sn_threshold = function(s, v) updateNumericInput(s, "batch_sn_threshold", value = v),
  batch_dir      = function(s, v) updateTextInput(s, "batch_dir", value = v),
  man_bases      = function(s, v) updateTextInput(s, "man_bases", value = v),
  man_sugars     = function(s, v) updateTextInput(s, "man_sugars", value = v),
  man_linkages   = function(s, v) updateTextInput(s, "man_linkages", value = v)
)

## ---- Help documents ---------------------------------------------------------
# The guides live in inst/help/, which means they resolve two ways: through
# system.file() when the package is installed, and relative to the checkout
# root when the app is run from a clone (where .module_dir is set and the
# package may not be installed at all).
.help_file <- function(name) {
  candidates <- character(0)
  if (!is.null(.module_dir))
    candidates <- c(candidates, file.path(.module_dir, "inst", "help", name))
  p <- tryCatch(system.file("help", name, package = "OligoMetProfiler"),
                error = function(e) "")
  if (nzchar(p)) candidates <- c(candidates, p)
  # Running the app directly out of inst/app/ in a checkout.
  candidates <- c(candidates, file.path("..", "help", name))
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) NULL else found[1]
}

# Render one guide, degrading to a link rather than an error: markdown is a
# Suggests-level dependency and the file is missing in an odd layout.
.help_ui <- function(name) {
  path <- .help_file(name)
  url <- paste0("https://github.com/nishi76/OligoMet-Profiler/blob/main/inst/help/",
                name)
  fallback <- function(why) {
    tags$div(class = "help-doc",
      tags$p(why),
      tags$p(tags$a(href = url, target = "_blank", rel = "noopener",
                    paste("Read", name, "on GitHub"))))
  }
  if (is.null(path))
    return(fallback(paste0("Could not find ", name, " in this installation.")))
  if (!requireNamespace("markdown", quietly = TRUE))
    return(fallback(paste0("Rendering the guides in-app needs the 'markdown' ",
                           "package -- install.packages(\"markdown\").")))
  tags$div(class = "help-doc", withMathJax(includeMarkdown(path)))
}

## =============================================================================
## UI
## =============================================================================
ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly",
                          primary = "#C8912E", secondary = "#2c3e50"),

  titlePanel("OligoMet Profiler"),

  ## Research-use-only banner. Deliberately above the fold and not
  ## dismissible: every number this app produces is a prediction, and the
  ## outputs get shared as files that leave the app. Kept as persistent
  ## chrome above the tab bar (not tab-scoped) -- see the About sub-tab
  ## under Help & Quick Start for the full disclaimer text.
  tags$div(class = "ruo-banner",
    tags$span(class = "ruo-tag", "RESEARCH USE ONLY"),
    tags$span("Not for diagnostic, clinical, or regulatory submission use. ",
              "All values are computed predictions, not measurements -- ",
              "confirm every assignment experimentally."),
    tags$a(href = "#", "See Help \u2192 About for the full disclaimer")
  ),

  tags$head(tags$style(HTML("
    .sidebar-section { margin-bottom: 18px; }
    .sidebar-section h5 { font-weight: 600; color: #2c3e50;
                          border-bottom: 1px solid #ecf0f1; padding-bottom: 4px; }
    .metric-card { background: #f8f9fa; border-radius: 6px; padding: 10px 14px;
                   text-align: center; border: 1px solid #dee2e6; }
    .metric-card .label { font-size: 11px; color: #6c757d; text-transform: uppercase;
                          letter-spacing: 0.5px; }
    .metric-card .value { font-size: 18px; font-weight: 600; color: #C8912E; }
    .manual-entry { background: #f8f9fa; border: 1px solid #dee2e6;
                    border-radius: 6px; padding: 14px 16px 6px;
                    margin-bottom: 14px; }
    .manual-entry h5 { font-weight: 600; color: #2c3e50; }
    .manual-entry .hint { font-size: 12px; color: #6c757d; }
    .manual-entry .form-group { margin-bottom: 8px; }
    .man-ok { color: #18632f; font-size: 13px; }
    .man-err { color: #a3231b; font-size: 13px; }
    .man-seq { font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
               font-size: 12px; word-break: break-all; }
    .help-panel { border: 1px solid #dee2e6; border-radius: 6px;
                  margin-bottom: 14px; background: #fff; }
    .help-panel > summary { cursor: pointer; padding: 10px 16px;
                            font-weight: 600; color: #2c3e50;
                            list-style: revert; }
    .help-panel[open] > summary { border-bottom: 1px solid #dee2e6; }
    .help-body { padding: 10px 16px 4px; }
    .help-doc { max-height: 60vh; overflow-y: auto; padding: 12px 4px 0;
                font-size: 13px; }
    .help-doc h1 { font-size: 20px; }
    .help-doc h2 { font-size: 17px; margin-top: 18px; }
    .help-doc h3 { font-size: 15px; margin-top: 14px; }
    .help-doc table { font-size: 12px; margin-bottom: 12px;
                      border-collapse: collapse; }
    .help-doc th, .help-doc td { border: 1px solid #dee2e6;
                                 padding: 3px 8px; }
    .help-doc pre { background: #f8f9fa; border: 1px solid #e9ecef;
                    border-radius: 4px; padding: 8px; font-size: 12px;
                    overflow-x: auto; }
    .help-doc code { font-size: 12px; }
    .help-doc img { max-width: 100%; }
    .ruo-banner { background: #fff8e6; border: 1px solid #f0d68a;
                  border-radius: 6px; padding: 8px 14px; margin-bottom: 10px;
                  font-size: 12px; color: #5c4813; line-height: 1.5; }
    .ruo-tag { display: inline-block; background: #a3231b; color: #fff;
               font-weight: 700; font-size: 10.5px; letter-spacing: 0.6px;
               border-radius: 3px; padding: 1px 7px; margin-right: 8px; }
    .about-block h6 { font-weight: 600; color: #2c3e50; font-size: 13px;
                      margin-bottom: 4px; }
    .about-block p { margin-bottom: 8px; font-size: 12px; color: #6c757d; }
    .adv-panel { border: 1px solid #dee2e6; border-radius: 6px;
                 margin-bottom: 14px; background: #fff; }
    .adv-panel > summary { cursor: pointer; padding: 8px 12px;
                           font-weight: 600; color: #2c3e50;
                           list-style: revert; font-size: 14px; }
    .adv-panel[open] > summary { border-bottom: 1px solid #dee2e6; }
    .adv-body { padding: 10px 12px 4px; }
    .session-row .btn { font-size: 12px; }
    .chem-hint { font-size: 12px; color: #6c757d; }
    /* Gold active-tab fill, matching the reconciled mockup */
    .nav-tabs .nav-link.active { background-color: #C8912E !important;
                                 color: #fff !important; font-weight: 700;
                                 border-color: #C8912E; }
    .nav-tabs .nav-link { color: #4B5563; }
    .workflow-status-row { display: flex; align-items: center; gap: 9px;
                           font-size: 12px; color: #374151; margin-bottom: 10px; }
    .status-dot { width: 20px; height: 20px; border-radius: 50%; display: flex;
                 align-items: center; justify-content: center; font-size: 11px;
                 color: #fff; flex-shrink: 0; font-weight: 700; }
    .status-dot.done { background: #22C55E; }
    .status-dot.pending { background: #D1D5DB; }
    .oligo-box { background: #FAFAFA; border: 1px solid #dee2e6; border-radius: 6px;
                padding: 10px 12px; font-size: 11.5px; color: #374151; }
    .stepper { display: flex; align-items: center; gap: 10px; margin-bottom: 20px; }
    .wf-step { flex: 1; background: #F1F1F1; border-radius: 8px; padding: 14px 6px;
              text-align: center; }
    .wf-step.wf-active { background: #C8912E; }
    .wf-step .wf-circ { width: 32px; height: 32px; border-radius: 50%; background: #fff;
                        color: #9CA3AF; display: flex; align-items: center;
                        justify-content: center; margin: 0 auto 6px; font-weight: 800; }
    .wf-step.wf-active .wf-circ { color: #C8912E; }
    .wf-step .wf-lbl { font-size: 11.5px; font-weight: 700; color: #6B7280; }
    .wf-step.wf-active .wf-lbl { color: #fff; }
    .wf-arrow { color: #C0C0C0; font-size: 16px; }
    .dashboard-callout { background: #FDFBF3; border-left: 4px solid #C8912E;
                         border-radius: 6px; padding: 14px 18px; }
    .dashboard-callout p { font-size: 13px; margin: 0 0 10px; color: #374151; }
    .placeholder-note { font-size: 11px; color: #B45309; background: #FFF7ED;
                        border: 1px solid #FED7AA; border-radius: 5px;
                        padding: 6px 10px; margin-top: 6px; }
  "))),

  tabsetPanel(id = "main_nav", type = "tabs",

    ## =========================================================================
    ## TAB: Dashboard
    ## =========================================================================
    tabPanel("Dashboard",
      tags$div(style = "height: 14px;"),
      fluidRow(
        column(4,
          tags$div(class = "sidebar-section",
            tags$h5("Session"),
            downloadButton("dl_session", "Save Session (JSON)", class = "btn-primary w-100"),
            tags$div(style = "height: 8px;"),
            fileInput("load_session_file", NULL, accept = ".json",
                      buttonLabel = "Load Session (JSON)...", placeholder = ""),
            tags$p(style = "font-size: 10.5px; color: #6c757d; margin-top: -8px;",
                   "Saves/restores the sequence, custom chemistry, and all ",
                   "parameters as a .json file -- not uploaded MS/batch files, ",
                   "and not results (see Analysis State below)."),
            tags$hr(style = "margin: 10px 0;"),
            downloadButton("dl_analysis_state", "Save Analysis State (.rds)", class = "btn-outline-primary w-100"),
            tags$div(style = "height: 8px;"),
            fileInput("load_analysis_state_file", NULL, accept = ".rds",
                      buttonLabel = "Load Analysis State (.rds)...", placeholder = ""),
            tags$p(style = "font-size: 10.5px; color: #6c757d; margin-top: -8px;",
                   "Saves/restores the library, batch results, statistics, and ",
                   "quantification tables, so reloading doesn't require re-running ",
                   "the batch. Raw mzML/raw files and per-hit MS2 spectra (needed only ",
                   "for redrawing mirror plots) are never included -- re-run with MS2 ",
                   "confirmation on if you need those back.")
          ),
          tags$div(class = "sidebar-section",
            tags$h5("Workflow Status"),
            uiOutput("dashboard_workflow_status")
          ),
          tags$div(class = "sidebar-section",
            tags$h5("Active Oligo"),
            uiOutput("dashboard_active_oligo")
          )
        ),
        column(8,
          uiOutput("dashboard_stepper"),
          tags$div(class = "dashboard-callout",
            tags$p(tags$b("Step 1: "), "Enter your oligonucleotide sequence and generate the MS1 library (Library Generation)."),
            tags$p(tags$b("Step 2: "), "Upload raw MS2 data to build an empirical fragment library (Empirical MS2 Library)."),
            tags$p(tags$b("Step 3: "), "Process batch samples -- deconvolution, matching, and calibration/quantification (Batch Processing)."),
            tags$p(style = "margin-bottom: 0;", tags$b("Step 4: "), "Run statistical analysis -- data matrix, univariate statistics, class comparison (Statistical Analysis).")
          ),
          tags$div(style = "height: 14px;"),
          verbatimTextOutput("status", placeholder = TRUE)
        )
      )
    ),

    ## =========================================================================
    ## TAB: Library Generation
    ## =========================================================================
    tabPanel("Library Generation",
      tags$div(style = "height: 14px;"),
      fluidRow(
        column(4,
          tags$div(class = "sidebar-section",
            tags$h5("Input"),
            textAreaInput("seq",
                          label = .with_info("Sequence (triplet, OligoDistiller, or FASTA)",
                            "Notation is auto-detected: starts with '>' = FASTA ",
                            "(paste a BioPharma Finder-exported record directly), ",
                            "starts with 'OH-' = OligoDistiller, otherwise triplet."),
                          value = .EXAMPLE_SEQS[[1]]$seq, rows = 3,
                          placeholder = "e.g. Te-sSe-sAe-sSe-... or OH-Am*-Gm*-...-OH or a pasted FASTA record"),
            fluidRow(
              column(8, selectInput("example_seq", NULL,
                        choices = c("Choose an example..." = "", names(.EXAMPLE_SEQS)),
                        selected = "")),
              column(4, actionButton("load_example", "Load", class = "btn-sm btn-outline-secondary w-100"))
            ),
            textInput("oligo_name", "Oligo name", value = "my_oligo"),
            textInput("output_prefix", "Output prefix", value = "my_oligo_metabolite"),
            tags$label("Save to folder (optional)", style = "font-size: 14px; font-weight: 500;"),
            fluidRow(
              column(8, textInput("output_dir", NULL, value = "",
                                  placeholder = "e.g. C:/Users/you/Documents/results")),
              column(4, if (.have_shinyfiles)
                shinyFiles::shinyDirButton("browse_output_dir", "Browse...", "Choose a folder",
                                           class = "btn-sm btn-outline-secondary w-100"))
            ),
            tags$p(style = "font-size: 11px; color: #6c757d; margin-top: -6px;",
                   "Leave blank to only use the download buttons below. If set, ",
                   "the workbook, report, PRM list, and acquisition method lists ",
                   "are also written directly to this folder when you click Generate Library."),
            selectInput("conj5", "5' conjugate", choices = .conj5_choices, selected = "none"),
            selectInput("conj3", "3' conjugate", choices = .conj3_choices, selected = "none")
          ),
          tags$div(class = "sidebar-section",
            tags$h5("Metabolite Generation"),
            fluidRow(
              column(6, numericInput("max_3p", "Max 3' trunc.",
                                     value = DEFAULT_PIPELINE_PARAMS$max_3p, min = 0, max = 50)),
              column(6, numericInput("max_5p", "Max 5' trunc.",
                                     value = DEFAULT_PIPELINE_PARAMS$max_5p, min = 0, max = 50))
            ),
            checkboxInput("endo", "Include endonuclease fragments",
                          value = DEFAULT_PIPELINE_PARAMS$endo),
            radioButtons("endo_sites", "Endo cleavage sites",
                         choices = c("All positions" = "all", "DNA gap only" = "gap"),
                         selected = DEFAULT_PIPELINE_PARAMS$endo_sites, inline = TRUE),
            numericInput("min_frag_len", "Min fragment length (nt)",
                        value = DEFAULT_PIPELINE_PARAMS$min_frag_len, min = 1, max = 20)
          ),
          tags$div(class = "sidebar-section",
            tags$h5("Mass & Isotope"),
            fluidRow(
              column(6, numericInput("z_min", "Min charge z",
                                     value = DEFAULT_PIPELINE_PARAMS$z_min, min = 1, max = 20)),
              column(6, numericInput("z_max", "Max charge z",
                                     value = DEFAULT_PIPELINE_PARAMS$z_max, min = 1, max = 30))
            ),
            fluidRow(
              column(6, numericInput("n_iso", "Isotope peaks",
                                     value = DEFAULT_PIPELINE_PARAMS$n_iso, min = 1, max = 20)),
              column(6, numericInput("max_oxid",
                        label = .with_info("Max PS oxid.",
                          "Trades completeness for target count: a higher cap covers deeper ",
                          "oxidation states but multiplies Charge Envelope rows and PRM/MS1 ",
                          "targets -- each metabolite is modeled at every oxidation level ",
                          "from 0 up to this cap."),
                        value = DEFAULT_PIPELINE_PARAMS$max_oxid, min = 0, max = 30))
            ),
            numericInput("h_offset",
                        label = .with_info("Envelope offset (Da)",
                          "0 = standard [M-zH]^z- charge envelope. Nonzero only to ",
                          "reproduce a legacy or lab-specific mass convention."),
                        value = DEFAULT_PIPELINE_PARAMS$h_offset, step = 0.001),
            checkboxInput("use_envipat", "Use enviPat for isotopes",
                          value = DEFAULT_PIPELINE_PARAMS$use_envipat)
          ),
          tags$details(class = "adv-panel",
            tags$summary("Orbitrap Acquisition Method"),
            tags$div(class = "adv-body",
              numericInput("method_length",
                          label = .with_info("Method length (min)",
                            "Total LC-MS run time; sets the RT window end for every row ",
                            "in the MS1 inclusion and MS2 PRM target lists."),
                          value = DEFAULT_PIPELINE_PARAMS$method_length, min = 1, max = 999),
              fluidRow(
                column(6, numericInput("ms2_z_min", "MS2 charge z min",
                                       value = DEFAULT_PIPELINE_PARAMS$ms2_z_min, min = 1, max = 30)),
                column(6, numericInput("ms2_z_max", "MS2 charge z max",
                                       value = DEFAULT_PIPELINE_PARAMS$ms2_z_max, min = 1, max = 30))
              ),
              numericInput("hcd_nce",
                          label = .with_info("HCD NCE (%)",
                            "20% is a PS-backbone starting point, not a validated instrument ",
                            "parameter -- optimize per method and per instrument."),
                          value = DEFAULT_PIPELINE_PARAMS$hcd_nce, min = 1, max = 200),
              fluidRow(
                column(6, numericInput("ms1_target_cap", "Max MS1 targets",
                          value = DEFAULT_PIPELINE_PARAMS$ms1_target_cap, min = 1, max = 150000)),
                column(6, numericInput("ms2_target_cap", "Max MS2 targets",
                          value = DEFAULT_PIPELINE_PARAMS$ms2_target_cap, min = 1, max = 150000))
              )
            )
          ),
          tags$hr(),
          actionButton("run_phase1", "Generate Library",
                       class = "btn-primary btn-lg w-100"),
          tags$p(style = "font-size: 11px; color: #6c757d; margin: 4px 0 10px;",
                 "Builds the MS1-only library -- .mgf/.msp, inclusion lists, and the ",
                 "Excel workbook. No MS2 library yet: that comes from real data ",
                 "in the Empirical MS2 Library tab.")
        ),

        column(8,
          tags$div(class = "manual-entry",
            tags$h5("Manual sequence entry"),
            tags$p(class = "hint",
                   "Type the three lines from your chemical analysis file. New to ",
                   "this? Open ", tags$strong("Help & Quick Start \u2192 Sequence Guide"), "."),
            fluidRow(
              column(4, textInput("man_bases", "Bases (5'->3')", value = "",
                                  placeholder = "TSASTTTSATAATGSTGG")),
              column(4, textInput("man_sugars", "Sugars", value = "",
                                  placeholder = "eeeeeeeeeeeeeeeeee")),
              column(4, textInput("man_linkages", "Linkages", value = "",
                                  placeholder = "sssssssssssssssss"))
            ),
            fluidRow(
              column(3, actionButton("man_submit", "Submit", class = "btn-primary w-100")),
              column(3, actionButton("man_example", "Fill example", class = "btn-outline-secondary w-100")),
              column(3, actionButton("man_clear", "Clear", class = "btn-outline-secondary w-100")),
              column(3, downloadButton("dl_fasta", "BPF FASTA", class = "btn-outline-secondary w-100"))
            ),
            htmlOutput("man_feedback")
          ),
          tags$details(class = "help-panel",
            tags$summary("Custom Chemistry (advanced)"),
            tags$div(class = "help-body",
              tags$p(class = "chem-hint",
                     "Add custom base/sugar/linkage/conjugate entries not in the ",
                     "standard dictionary."),
              DT::dataTableOutput("custom_chem"),
              tags$div(style = "height: 6px;"),
              fluidRow(
                column(3, actionButton("add_row", "Add Row", class = "btn-sm btn-outline-primary w-100")),
                column(3, actionButton("remove_row", "Remove Row", class = "btn-sm btn-outline-secondary w-100"))
              )
            )
          ),
          conditionalPanel(
            condition = "output.status_ready == 'true'",
            tags$hr(),
            fluidRow(
              column(3, tags$div(class = "metric-card",
                tags$div(class = "label", "Formula"), tags$div(class = "value", textOutput("m_formula")))),
              column(2, tags$div(class = "metric-card",
                tags$div(class = "label", "Mono Mass (Da)"), tags$div(class = "value", textOutput("m_mono_mass")))),
              column(2, tags$div(class = "metric-card",
                tags$div(class = "label", "Avg Mass (Da)"), tags$div(class = "value", textOutput("m_avg_mass")))),
              column(1, tags$div(class = "metric-card",
                tags$div(class = "label", "Length"), tags$div(class = "value", textOutput("m_length")))),
              column(2, tags$div(class = "metric-card",
                tags$div(class = "label", "Metabolites"), tags$div(class = "value", textOutput("m_n_mets")))),
              column(2, tags$div(class = "metric-card",
                tags$div(class = "label", "PRM Entries"), tags$div(class = "value", textOutput("m_n_prm"))))
            ),
            tags$div(style = "height: 10px;"),
            tabsetPanel(
              tabPanel("Charge Envelope", plotOutput("plot_envelope", height = "380px")),
              tabPanel("Truncation Series", plotOutput("plot_truncation", height = "380px")),
              tabPanel("Isotope Pattern", plotOutput("plot_isotope", height = "380px")),
              tabPanel("Oxidation Series", plotOutput("plot_oxidation", height = "380px")),
              tabPanel("MS2 Explorer (predicted)",
                tags$div(style = "padding-top: 12px;",
                  fluidRow(
                    column(8, tags$p(style = "font-size: 12px; color: #6c757d;",
                      "Predicted (theoretical) MS2 library for interactive browsing --",
                      " see Empirical MS2 Library for spectra built from real data.")),
                    column(4, actionButton("build_ms2_explorer", "Build / Refresh Library",
                      class = "btn-sm btn-outline-primary w-100"))
                  ),
                  conditionalPanel(
                    condition = "output.ms2_explorer_ready == 'true'",
                    DT::dataTableOutput("ms2_explorer_table"),
                    tags$hr(),
                    uiOutput("ms2_explorer_title"),
                    fluidRow(
                      column(7, plotOutput("ms2_explorer_plot", height = "300px")),
                      column(5, DT::dataTableOutput("ms2_explorer_peaks"))
                    )
                  )
                )
              )
            ),
            tags$div(style = "height: 10px;"),
            tags$div(class = "sidebar-section",
              tags$h5("Downloads -- Library (MS1-only)"),
              fluidRow(
                column(4, downloadButton("dl_workbook", "Excel Workbook (.xlsx)", class = "btn-success w-100")),
                column(4, downloadButton("dl_report", "HTML Report (.html)", class = "btn-info w-100")),
                column(4, downloadButton("dl_prm", "PRM Inclusion List (.csv)", class = "btn-warning w-100"))
              ),
              tags$div(style = "height: 8px;"),
              fluidRow(
                column(3, downloadButton("dl_ms1_inclusion", "MS1 Inclusion List (.csv)", class = "btn-outline-secondary w-100")),
                column(3, downloadButton("dl_ms2_prm", "MS2 PRM Target List (.csv)", class = "btn-outline-secondary w-100")),
                column(3, downloadButton("dl_frag_ref", "MS2 Fragment Reference (.csv)", class = "btn-outline-secondary w-100")),
                column(3, downloadButton("dl_all", "Download All (.zip)", class = "btn-dark w-100"))
              ),
              tags$div(style = "height: 8px;"),
              fluidRow(
                column(6, downloadButton("dl_ms1_mgf", "MS1 library (.mgf)", class = "btn-outline-primary w-100")),
                column(6, downloadButton("dl_ms1_msp", "MS1 library (.msp)", class = "btn-outline-primary w-100"))
              )
            )
          )
        )
      )
    ),

    ## =========================================================================
    ## TAB: Empirical MS2 Library
    ## =========================================================================
    tabPanel("Empirical MS2 Library",
      tags$div(style = "height: 14px;"),
      fluidRow(
        column(4,
          tags$div(class = "sidebar-section",
            tags$h5("Acquisition Mode"),
            radioButtons("ms2_acquisition_mode", NULL,
                         choices = c("DDA" = "dda", "PRM" = "prm", "DIA" = "dia", "AcquireX" = "acquirex"),
                         selected = "dda", inline = TRUE),
            tags$p(style = "font-size: 10px; color: #B45309; background: #FFF7ED; border: 1px solid #FED7AA; border-radius: 5px; padding: 5px 8px; margin-top: 4px;",
                   "DIA: wide-isolation precursor deconvolution is not implemented yet -- ",
                   "matching below assumes a 1:1 precursor\u2194MS2 scan (DDA/PRM/AcquireX). ",
                   "This selector is currently informational only.")
          ),
          tags$div(class = "sidebar-section",
            tags$h5("MS Matching (single file)"),
            checkboxInput("enable_ms", "Enable MS matching",
                          value = DEFAULT_PIPELINE_PARAMS$enable_ms),
            conditionalPanel(
              condition = "input.enable_ms == true",
              fileInput("ms_file", "Upload MS file (.mzML, .mzXML, .raw, .csv)",
                        accept = c(".mzML", ".mzXML", ".mzml", ".mzxml", ".raw", ".csv", ".txt")),
              numericInput("ppm_tol", "MS1 tolerance (ppm)",
                          value = DEFAULT_PIPELINE_PARAMS$ppm_tol, min = 1, max = 50),
              radioButtons("noise_mode", "Background/noise threshold",
                           choices = c("Signal-to-noise multiple" = "sn",
                                       "Fixed intensity value" = "fixed"),
                           selected = "sn", inline = TRUE),
              conditionalPanel(
                condition = "input.noise_mode == 'sn'",
                numericInput("sn_threshold", "S/N threshold (x noise level)",
                            value = DEFAULT_PIPELINE_PARAMS$sn_threshold, min = 0, step = 0.5)
              ),
              conditionalPanel(
                condition = "input.noise_mode == 'fixed'",
                numericInput("min_intensity", "Fixed intensity threshold",
                            value = DEFAULT_PIPELINE_PARAMS$min_intensity, min = 0)
              ),
              checkboxGroupInput("adducts", "Adducts",
                                 choices = c("H", "Na", "K", "NH4"),
                                 selected = DEFAULT_PIPELINE_PARAMS$adducts, inline = TRUE),
              fluidRow(
                column(6, numericInput("frag_tol_ppm", "Fragment tol (ppm)",
                          value = DEFAULT_PIPELINE_PARAMS$frag_tol_ppm, min = 5, max = 100)),
                column(6, numericInput("frag_z_max", "Fragment max z",
                          value = DEFAULT_PIPELINE_PARAMS$frag_z_max, min = 1, max = 5))
              )
            )
          ),
          tags$div(class = "sidebar-section",
            tags$h5("Batch Raw Files (multi-file)"),
            checkboxInput("enable_batch", "Enable batch processing",
                          value = DEFAULT_PIPELINE_PARAMS$enable_batch),
            conditionalPanel(
              condition = "input.enable_batch == true",
              fileInput("batch_files", "Upload raw files (.mzML, .mzXML, .raw)", multiple = TRUE,
                        accept = c(".mzML", ".mzml", ".mzXML", ".mzxml", ".raw")),
              uiOutput("batch_upload_notice"),
              tags$label("...or point at a local folder (optional)",
                         style = "font-size: 13px; font-weight: 500;"),
              fluidRow(
                column(8, textInput("batch_dir", NULL, value = "",
                                    placeholder = "e.g. /path/to/mzml_folder")),
                column(4, if (.have_shinyfiles)
                  shinyFiles::shinyDirButton("browse_batch_dir", "Browse...", "Choose a folder",
                                             class = "btn-sm btn-outline-secondary w-100"))
              ),
              checkboxInput("batch_run_ms2", "Confirm hits with MS2 (build empirical library)",
                            value = TRUE),
              fluidRow(
                column(6, numericInput("batch_n_workers", "Parallel workers",
                          value = DEFAULT_PIPELINE_PARAMS$batch_n_workers, min = 1, max = 64)),
                column(6, numericInput("batch_deconv_ppm", "Deconv mass tol (ppm)",
                          value = DEFAULT_PIPELINE_PARAMS$batch_deconv_ppm, min = 1, max = 100))
              ),
              radioButtons("batch_noise_mode", "Background/noise threshold",
                           choices = c("Signal-to-noise multiple" = "sn",
                                       "Fixed intensity value" = "fixed"),
                           selected = "sn", inline = TRUE),
              conditionalPanel(
                condition = "input.batch_noise_mode == 'sn'",
                numericInput("batch_sn_threshold", "S/N threshold (x noise level)",
                            value = DEFAULT_PIPELINE_PARAMS$sn_threshold, min = 0, step = 0.5)
              ),
              conditionalPanel(
                condition = "input.batch_noise_mode == 'fixed'",
                numericInput("batch_min_intensity", "Fixed intensity threshold",
                            value = DEFAULT_PIPELINE_PARAMS$min_intensity, min = 0)
              ),
              tags$p(style = "font-size: 10.5px; color: #6c757d;",
                     "Sample metadata (group/timepoint/type) is set on the Batch ",
                     "Processing tab -- not needed just to build the empirical library.")
            )
          ),
          conditionalPanel(
            condition = "output.library_ready == 'true'",
            actionButton("run_ms2_library", "Generate Empirical MS2 Library",
                         class = "btn-primary btn-lg w-100")
          ),
          conditionalPanel(
            condition = "output.library_ready != 'true'",
            tags$p(style = "font-size: 11px; color: #6c757d;",
                   "Run \"Generate Library\" on the Library Generation tab first.")
          )
        ),

        column(8,
          conditionalPanel(
            condition = "input.enable_ms == true && output.status_ready == 'true'",
            fluidRow(
              column(4, tags$div(class = "metric-card",
                tags$div(class = "label", "MS1 Matches (single file)"), tags$div(class = "value", textOutput("m_ms1_matches")))),
              column(4, tags$div(class = "metric-card",
                tags$div(class = "label", "Annotated Mets"), tags$div(class = "value", textOutput("m_annotated")))),
              column(4, tags$div(class = "metric-card",
                tags$div(class = "label", "Putative IDs"), tags$div(class = "value", textOutput("m_confident"))))
            ),
            tags$div(style = "height: 10px;")
          ),
          conditionalPanel(
            condition = "output.batch_ready == 'true'",
            fluidRow(
              column(4, tags$div(class = "metric-card",
                tags$div(class = "label", "MS2 Confirmed Hits"), tags$div(class = "value", textOutput("ms2lib_n_matched")))),
              column(4, tags$div(class = "metric-card",
                tags$div(class = "label", "MS2 Spectra Extracted"), tags$div(class = "value", textOutput("ms2lib_n_spectra")))),
              column(4, tags$div(class = "metric-card",
                tags$div(class = "label", "Consensus Spectra Built"), tags$div(class = "value", textOutput("ms2lib_n_consensus"))))
            ),
            tags$div(style = "height: 10px;")
          ),
          tabsetPanel(
            tabPanel("MS2 Mirror Plot (single file)",
              tags$div(style = "padding-top: 12px;",
                conditionalPanel(
                  condition = "output.ms2_mirror_ready == 'true'",
                  DT::dataTableOutput("ms2_mirror_table"),
                  tags$hr(),
                  plotOutput("ms2_mirror_plot", height = "400px")
                ),
                conditionalPanel(
                  condition = "output.ms2_mirror_ready != 'true'",
                  tags$p(style = "color: #6c757d;",
                    "Enable MS matching, upload an MS2-containing file, and run to see the mirror plot here.")
                )
              )
            ),
            tabPanel("Empirical Library",
              tags$div(style = "padding-top: 12px;",
                conditionalPanel(
                  condition = "output.batch_ready == 'true'",
                  DT::dataTableOutput("empirical_ms2_summary_table"),
                  tags$div(style = "height: 8px;"),
                  fluidRow(
                    column(3, downloadButton("dl_empirical_ms2_msp", "Empirical MS2 (.msp)", class = "btn-outline-primary w-100")),
                    column(3, downloadButton("dl_empirical_ms2_summary", "Summary (.csv)", class = "btn-outline-primary w-100")),
                    column(3, downloadButton("dl_batch_annotated_msp", "Annotated MS2 (.msp)", class = "btn-outline-primary w-100")),
                    column(3, downloadButton("dl_batch_mirror_pdf", "Mirror Plots (.pdf)", class = "btn-outline-primary w-100"))
                  )
                ),
                conditionalPanel(
                  condition = "output.batch_ready != 'true'",
                  tags$p(style = "color: #6c757d;",
                    "Upload raw files above and click \"Generate Empirical MS2 Library\" to see results here.")
                ),
                tags$hr(),
                tags$h6("Predicted MS2 Library (theoretical)"),
                fluidRow(
                  column(6, tagAppendAttributes(
                    downloadButton("dl_ms2_mgf", "Predicted MS2 (.mgf) \u26a0", class = "btn-outline-secondary w-100"),
                    title = "MS2 intensities are placeholders (rule-based heuristic) -- match on m/z only.")),
                  column(6, tagAppendAttributes(
                    downloadButton("dl_ms2_msp", "Predicted MS2 (.msp) \u26a0", class = "btn-outline-secondary w-100"),
                    title = "MS2 intensities are placeholders (rule-based heuristic) -- match on m/z only."))
                )
              )
            )
          )
        )
      )
    ),

    ## =========================================================================
    ## TAB: Batch Processing
    ## =========================================================================
    tabPanel("Batch Processing",
      tags$div(style = "height: 14px;"),
      fluidRow(
        column(4,
          tags$div(class = "sidebar-section",
            tags$h5("Sample Metadata"),
            fluidRow(
              column(6, downloadButton("dl_batch_meta_template", "Download CSV template",
                                        class = "btn-outline-secondary btn-sm w-100")),
              column(6, fileInput("batch_meta_csv", NULL, accept = ".csv",
                                   placeholder = "Upload filled-in sample info CSV..."))
            ),
            tags$p(style = "font-size: 10.5px; color: #6c757d; margin-top: -10px;",
                   "Columns: sample, group, timepoint, sample_type, concentration. ",
                   "Timepoint must be numeric (0, 4, 24) -- flagged live below if not."),
            uiOutput("batch_meta_upload_status"),
            DT::DTOutput("sample_meta_table"),
            uiOutput("sample_meta_timepoint_warn"),
            uiOutput("sample_meta_type_warn"),
            tags$p(style = "font-size: 11px; color: #6c757d; margin-top: 4px;",
                   "One row per uploaded file (see Empirical MS2 Library tab to upload ",
                   "raw files). Fill in Group (2+ groups) or Timepoint (time series) ",
                   "before running -- leave both blank to only extract and match ",
                   "features, with no statistics.")
          ),
          tags$div(class = "sidebar-section",
            tags$h5("Processing Options"),
            tags$p(style = "font-size: 10.5px; color: #6c757d; margin-top: -6px;",
                   "Re-uses the same noise-threshold and confirmation settings as ",
                   "Empirical MS2 Library -- set there before running here.")
          ),
          tags$details(class = "adv-panel", open = NA,
            tags$summary("Calibration & Quantification"),
            tags$div(class = "adv-body",
              selectizeInput("absolute_quant_mets",
                             "Absolute-quantify these metabolites (calibration curve)",
                             choices = character(0), multiple = TRUE,
                             options = list(placeholder = "None selected -- every metabolite uses relative quantification")),
              fluidRow(
                column(6, selectInput("calibration_weighting", "Calibration curve weighting",
                          choices = c("1/x² weighted (recommended)" = "1/x2",
                                      "1/x weighted" = "1/x",
                                      "Unweighted (OLS)" = "none"))),
                column(6, selectInput("control_group", "Control group (relative quant)",
                          choices = character(0)))
              ),
              tags$p(style = "font-size: 10.5px; color: #6c757d; margin-top: -10px;",
                     "Selected metabolites are absolute-quantified from their own ",
                     "calibration curve (Sample Type = standard rows with a ",
                     "Concentration value; linear regression). Every other ",
                     "metabolite is relative-quantified.")
            )
          ),
          conditionalPanel(
            condition = "output.library_ready == 'true'",
            actionButton("run_phase2", "Run Batch Processing",
                         class = "btn-primary btn-lg w-100")
          ),
          conditionalPanel(
            condition = "output.library_ready != 'true'",
            tags$p(style = "font-size: 11px; color: #6c757d;",
                   "Run \"Generate Library\" on the Library Generation tab first.")
          )
        ),

        column(8,
          tabsetPanel(
            tabPanel("Batch Results",
              conditionalPanel(
                condition = "output.batch_ready == 'true'",
                tags$div(style = "padding-top: 12px;",
                  DT::DTOutput("batch_matches_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_batch_tsv", "Download combined features (.tsv)", class = "btn-outline-primary"),
                  tags$hr(),
                  tags$h6("MS2 Mirror Plot"),
                  tags$p(style = "font-size: 12px; color: #6c757d;",
                    "Select a row above with MS2 confirmation data (n_ms2_peaks > 0) to plot it."),
                  plotOutput("batch_mirror_plot", height = "380px")
                )
              ),
              conditionalPanel(
                condition = "output.batch_ready != 'true'",
                tags$p(style = "padding-top: 12px; color: #6c757d;",
                       "Upload files (Empirical MS2 Library tab) and click Run Batch Processing to see results here.")
              )
            ),
            tabPanel("Unidentified Peaks",
              conditionalPanel(
                condition = "output.batch_ready == 'true'",
                tags$div(style = "padding-top: 12px;",
                  DT::DTOutput("unmatched_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_unmatched_csv", "Download unidentified peaks (.csv)", class = "btn-outline-primary")
                )
              )
            ),
            tabPanel("Calibration Curves",
              conditionalPanel(
                condition = "output.quant_ready == 'true'",
                tags$div(style = "padding-top: 12px;",
                  tags$h6("Calibration curves"),
                  tags$p(style = "font-size: 12px; color: #6c757d;",
                    "One row per metabolite selected for absolute quantification. ",
                    "\"note\" explains why a curve is missing instead of just disappearing."),
                  DT::DTOutput("calibration_curves_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_calibration_curves_csv", "Download calibration curves (.csv)", class = "btn-outline-primary"),
                  tags$hr(),
                  tags$h6("Absolute quantification"),
                  DT::DTOutput("quant_absolute_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_quant_absolute_csv", "Download absolute quantification (.csv)", class = "btn-outline-primary"),
                  tags$hr(),
                  tags$h6("Relative quantification (fold-change)"),
                  DT::DTOutput("quant_relative_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_quant_relative_csv", "Download relative quantification (.csv)", class = "btn-outline-primary")
                )
              ),
              conditionalPanel(
                condition = "output.quant_ready != 'true'",
                uiOutput("quant_not_ready_note")
              )
            ),
            tabPanel("Degradation Summary",
              conditionalPanel(
                condition = "output.batch_ready == 'true'",
                tags$div(style = "padding-top: 12px;",
                  tags$h6("% Degradation per sample"),
                  DT::DTOutput("degradation_per_sample_table"),
                  tags$div(style = "height: 8px;"),
                  tags$h6("Composition by class"),
                  plotOutput("plot_degradation_composition", height = "300px"),
                  DT::DTOutput("degradation_composition_table"),
                  tags$div(style = "height: 8px;"),
                  tags$h6("Top degradant species"),
                  DT::DTOutput("degradation_top_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_degradation_csv", "Download degradation summary (.csv)", class = "btn-outline-primary")
                )
              )
            )
          )
        )
      )
    ),

    ## =========================================================================
    ## TAB: Statistical Analysis
    ## =========================================================================
    tabPanel("Statistical Analysis",
      tags$div(style = "height: 14px;"),
      fluidRow(
        column(4,
          tags$div(class = "sidebar-section",
            tags$h5("Experimental Design"),
            uiOutput("stats_design_summary"),
            fluidRow(
              column(6, selectInput("stats_group_a", "Group A (optional override)", choices = character(0))),
              column(6, selectInput("stats_group_b", "Group B (optional override)", choices = character(0)))
            ),
            tags$p(style = "font-size: 10px; color: #6c757d; margin-top: -8px;",
                   "Leave both blank to auto-pick the first two Group values found. ",
                   "Only used for a 2-group design; ignored for 3+ groups or time-course."),
            fluidRow(
              column(6, selectInput("stats_padjust", "P-value adjustment",
                        choices = c("Benjamini-Hochberg" = "BH", "Bonferroni" = "bonferroni",
                                    "Holm" = "holm", "None" = "none"), selected = "BH")),
              column(6, numericInput("stats_log2fc_threshold", "Log2FC highlight threshold",
                        value = 1, min = 0, step = 0.1))
            )
          ),
          tags$details(class = "adv-panel",
            tags$summary("Charge Grouping (Advanced)"),
            tags$div(class = "adv-body",
              tags$p(style = "font-size: 9.5px; color: #6c757d; font-style: italic;",
                     "Groups co-eluting peaks by neutral mass agreement (charge_group.py logic)."),
              fluidRow(
                column(6, numericInput("cg_rt_tol", "RT tolerance (min)", value = 0.15, min = 0, step = 0.01)),
                column(6, numericInput("cg_mass_tol_ppm", "Mass tolerance (ppm)", value = 20, min = 1))
              ),
              numericInput("cg_min_charge_states", "Min charge states", value = 2, min = 1, max = 10),
              actionButton("cg_rerun", "Re-run Aggregation", class = "btn-sm w-100",
                          style = "background:#2F6FED;color:#fff;"),
              uiOutput("cg_rerun_status")
            )
          ),
          actionButton("run_stats", "Run Statistical Analysis", class = "btn-primary btn-lg w-100"),
          tags$p(style = "font-size: 10.5px; color: #6c757d;",
                 "Statistics already run automatically as part of Batch Processing when ",
                 "Group/Timepoint is filled in -- this button just re-applies the ",
                 "P-adjustment/threshold settings above to the existing results.")
        ),

        column(8,
          tabsetPanel(
            tabPanel("Data Matrix",
              conditionalPanel(
                condition = "output.batch_ready == 'true'",
                tags$div(style = "padding-top: 12px;",
                  tags$p(style = "font-size: 12px; color: #6c757d;",
                    "One row per metabolite (kind/z), one column per sample -- ",
                    "max-intensity match per metabolite per sample, matching how the ",
                    "batch pipeline already collapses charge states before matching. ",
                    "A genuinely separate \"raw, one row per charge state\" matrix and a ",
                    "true charge_group.py-style post-hoc re-aggregation are not implemented ",
                    "yet -- see the Charge Grouping panel note."),
                  fluidRow(
                    column(6, downloadButton("dl_data_matrix_wide", "Data Matrix -- wide (.csv)", class = "btn-outline-primary w-100")),
                    column(6, downloadButton("dl_data_matrix_long", "Data Matrix -- long (.csv)", class = "btn-outline-primary w-100"))
                  )
                )
              ),
              conditionalPanel(
                condition = "output.batch_ready != 'true'",
                tags$p(style = "padding-top: 12px; color: #6c757d;",
                       "Run Batch Processing to see the data matrix here.")
              )
            ),
            tabPanel("Univariate Statistics",
              conditionalPanel(
                condition = "output.stats_ready == 'true'",
                tags$div(style = "padding-top: 12px;",
                  uiOutput("stats_met_selector"),
                  plotOutput("plot_stats_main", height = "340px"),
                  DT::DTOutput("stats_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_stats_csv", "Download statistics table (.csv)", class = "btn-outline-primary")
                )
              ),
              conditionalPanel(
                condition = "output.stats_ready != 'true'",
                tags$p(style = "padding-top: 12px; color: #6c757d;",
                       "Fill in Group or Timepoint in the batch sample table to see comparisons here. ",
                       "Time-course vs. group comparison is auto-detected from that column.")
              )
            ),
            tabPanel("Class Comparison",
              conditionalPanel(
                condition = "output.kind_stats_ready == 'true'",
                tags$div(style = "padding-top: 12px;",
                  tags$p(style = "font-size: 12px; color: #6c757d;",
                    "Same comparison, rolled up by metabolite class (parent / 5' exonuclease / ",
                    "3' exonuclease / endonuclease) instead of per-metabolite."),
                  DT::DTOutput("kind_stats_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_kind_stats_csv", "Download class comparison table (.csv)", class = "btn-outline-primary")
                )
              )
            ),
            tabPanel("Time Course",
              conditionalPanel(
                condition = "output.batch_ready == 'true'",
                tags$div(style = "padding-top: 12px;",
                  tags$p(style = "font-size: 12px; color: #6c757d;",
                    "Per-metabolite time-course comparison already runs automatically under ",
                    "Univariate Statistics when Timepoint is filled in -- this is the same ",
                    "data, several metabolites overlaid on one shared scale so their kinetics ",
                    "can be compared directly. Calibration standards/QC/blanks are excluded."),
                  selectizeInput("tc_met_ids", "Metabolites to plot", choices = character(0),
                                 multiple = TRUE,
                                 options = list(placeholder = "Defaults to the top few by signal")),
                  checkboxInput("tc_normalize",
                                "Normalize to earliest timepoint (recommended -- puts metabolites of very different abundance on one scale)",
                                value = TRUE),
                  uiOutput("time_course_note"),
                  plotOutput("plot_time_course_trend", height = "340px"),
                  DT::DTOutput("time_course_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_time_course_csv", "Download time-course summary (.csv)", class = "btn-outline-primary")
                )
              ),
              conditionalPanel(
                condition = "output.batch_ready != 'true'",
                tags$p(style = "padding-top: 12px; color: #6c757d;",
                       "Run Batch Processing with Timepoint filled in to see a multi-metabolite trend here.")
              )
            ),
            tabPanel("Multivariate",
              conditionalPanel(
                condition = "output.batch_ready == 'true'",
                tags$div(style = "padding-top: 12px;",
                  tags$p(style = "font-size: 12px; color: #6c757d;",
                    "PCA and hierarchical clustering of samples by metabolite profile. ",
                    "Calibration standards, QC, and blanks are excluded (Sample Type on the ",
                    "Batch Processing tab), same as everywhere else in this tab."),
                  selectizeInput("mv_met_ids", "Metabolites to include", choices = character(0),
                                 multiple = TRUE, options = list(placeholder = "All metabolites (default)")),
                  fluidRow(
                    column(6, checkboxInput("mv_log", "Log2-transform", value = TRUE)),
                    column(6, checkboxInput("mv_scale", "Unit-variance scale", value = TRUE))
                  ),
                  uiOutput("multivariate_note"),
                  tags$h6("PCA -- sample scores"),
                  plotOutput("plot_pca_scores", height = "320px"),
                  tags$p(style = "font-size: 11px; color: #6c757d; margin-top: -4px;",
                         "Loadings: which metabolites drive each component, ranked by |PC1|."),
                  DT::DTOutput("pca_loadings_table"),
                  fluidRow(
                    column(6, downloadButton("dl_pca_scores_csv", "Download PCA scores (.csv)", class = "btn-outline-primary w-100")),
                    column(6, downloadButton("dl_pca_loadings_csv", "Download PCA loadings (.csv)", class = "btn-outline-primary w-100"))
                  ),
                  tags$hr(),
                  tags$h6("Hierarchical clustering"),
                  numericInput("mv_k", "Number of clusters", value = 2, min = 2, step = 1),
                  plotOutput("plot_dendrogram", height = "320px"),
                  DT::DTOutput("hclust_clusters_table"),
                  tags$div(style = "height: 8px;"),
                  downloadButton("dl_hclust_clusters_csv", "Download cluster assignments (.csv)", class = "btn-outline-primary")
                )
              ),
              conditionalPanel(
                condition = "output.batch_ready != 'true'",
                tags$p(style = "padding-top: 12px; color: #6c757d;",
                       "Run Batch Processing to see PCA / clustering here.")
              )
            )
          )
        )
      )
    ),

    ## =========================================================================
    ## TAB: Ask OligoMet
    ## =========================================================================
    tabPanel("Ask OligoMet",
      tags$div(style = "height: 14px; padding: 0 15px;"),
      fluidRow(
        column(4,
          tags$div(class = "sidebar-section",
            tags$h5("LLM Backend"),
            radioButtons("agent_backend", NULL,
                         choices = c("Anthropic Claude" = "anthropic", "OpenAI" = "openai"),
                         inline = TRUE),
            textOutput("agent_key_status"),
            tags$p(style = "font-size: 11px; color: #6c757d; margin-top: 6px;",
                   "API keys are read from the ANTHROPIC_API_KEY / OPENAI_API_KEY ",
                   "environment variables -- never typed into this page.")
          ),
          tags$div(class = "sidebar-section",
            tags$h5("What it can see"),
            tags$p(style = "font-size: 12px; color: #6c757d;",
                   "Whatever sequence / library / batch results are already loaded in ",
                   "this session, across all tabs. It does not write back into the main ",
                   "tabs -- re-run the relevant tab yourself if you want a suggestion ",
                   "it made reflected in the main view.")
          )
        ),
        column(8,
          tags$div(style = paste("border:1px solid #ddd; border-radius:8px;",
                                 "padding:10px; min-height:340px; max-height:480px;",
                                 "overflow-y:auto; background:#fafafa;"),
            uiOutput("agent_chat_html")
          ),
          tags$div(style = "height: 8px;"),
          fluidRow(
            column(9, textAreaInput("agent_input", NULL, width = "100%", rows = 2,
                     placeholder = "Ask about your sequence, a match, or what's degrading...")),
            column(3,
              actionButton("agent_send", "Send", class = "btn-primary w-100"),
              tags$div(style = "height: 4px;"),
              actionButton("agent_clear", "Clear", class = "btn-outline-secondary btn-sm w-100"))
          ),
          conditionalPanel(condition = "output.agent_busy_flag == 'true'",
            tags$p(style = "color: #6c757d; font-size: 12px;", "Thinking..."))
        )
      )
    ),

    ## =========================================================================
    ## TAB: Help & Quick Start
    ## =========================================================================
    tabPanel("Help & Quick Start",
      tags$div(style = "height: 14px;"),
      tabsetPanel(
        id = "help_tabs",
        tabPanel("Quick Start", tags$div(style = "padding-top: 12px;", uiOutput("help_quickstart"))),
        tabPanel("Sequence Guide", tags$div(style = "padding-top: 12px;", uiOutput("help_sequence"))),
        tabPanel("Modifications", tags$div(style = "padding-top: 12px;", uiOutput("help_modifications"))),
        tabPanel("No-Shiny (CLI)", tags$div(style = "padding-top: 12px;", uiOutput("help_quickstart_cli"))),
        tabPanel("About",
          tags$div(class = "about-block", style = "padding-top: 16px;",
            tags$h6("About"),
            tags$p(
              tags$strong("OligoMet Profiler"), HTML("&mdash;"),
              textOutput("about_version", inline = TRUE), tags$br(),
              paste0(OLIGOMET_AUTHOR_ROLE, ": ", OLIGOMET_AUTHOR),
              HTML(paste0("(<a href=\"mailto:", OLIGOMET_AUTHOR_EMAIL, "\">",
                          OLIGOMET_AUTHOR_EMAIL, "</a>)")), tags$br(),
              paste0(OLIGOMET_AUTHOR_TITLE, ", ", OLIGOMET_AUTHOR_AFFILIATION,
                     " -- an independent personal project, not a ",
                     OLIGOMET_AUTHOR_AFFILIATION, " product."), tags$br(),
              tags$a(href = OLIGOMET_URL, target = "_blank", rel = "noopener", OLIGOMET_URL),
              HTML("&mdash; released under the MIT licence.")),
            tags$h6("Disclaimer"),
            lapply(OLIGOMET_DISCLAIMER, tags$p)
          )
        )
      )
    )
  )
)

## =============================================================================
## SERVER
## =============================================================================
server <- function(input, output, session) {

  ## ---- Reactive values -----------------------------------------------------
  rv <- reactiveValues(
    spec = NULL, mets = NULL, dict = NULL, prm = NULL,
    ms_results = NULL, wb_path = NULL, report_path = NULL,
    plots = list(), ready = FALSE,
    status_text = "Enter a sequence and click \"Generate Library\" on the Library Generation tab.\n",
    batch_features = NULL, batch_ms_results = NULL,
    sample_meta = NULL, stats_results = NULL, kind_stats_results = NULL,
    quant_results = NULL,
    library_ready = FALSE,
    agent_messages = list(), agent_ctx = list(), agent_busy = FALSE
  )

  ## ---- Dashboard tab ---------------------------------------------------------
  # Pure display, computed off rv$* state that's already set elsewhere --
  # this tab adds no new pipeline logic, just an at-a-glance view of it.
  output$dashboard_workflow_status <- renderUI({
    step2_done <- !is.null(rv$batch_ms_results) &&
      (length(rv$batch_ms_results$ms2_confirmations) > 0 && nrow(rv$batch_ms_results$ms2_confirmations) > 0)
    step3_done <- !is.null(rv$batch_ms_results) && !is.null(rv$batch_ms_results$ms1_matches) &&
      nrow(rv$batch_ms_results$ms1_matches) > 0
    step4_done <- !is.null(rv$stats_results) || !is.null(rv$quant_results)
    steps <- list(
      list("1. Library Generation", isTRUE(rv$library_ready)),
      list("2. Empirical MS2 Library", step2_done),
      list("3. Batch Processing", step3_done),
      list("4. Statistical Analysis", step4_done),
      list("5. Complete", isTRUE(rv$library_ready) && step3_done && step4_done)
    )
    tagList(lapply(steps, function(s) {
      tags$div(class = "workflow-status-row",
        tags$div(class = paste("status-dot", if (s[[2]]) "done" else "pending"),
                 if (s[[2]]) "✓" else ""),
        s[[1]])
    }))
  })

  output$dashboard_active_oligo <- renderUI({
    if (is.null(rv$spec)) {
      return(tags$div(class = "oligo-box", "No sequence loaded yet."))
    }
    tags$div(class = "oligo-box",
      sprintf("%s (%d-mer)", input$oligo_name %||% "my_oligo", rv$spec$n %||% 0))
  })

  output$dashboard_stepper <- renderUI({
    step2_done <- !is.null(rv$batch_ms_results) &&
      (length(rv$batch_ms_results$ms2_confirmations) > 0 && nrow(rv$batch_ms_results$ms2_confirmations) > 0)
    step3_done <- !is.null(rv$batch_ms_results) && !is.null(rv$batch_ms_results$ms1_matches) &&
      nrow(rv$batch_ms_results$ms1_matches) > 0
    step4_done <- !is.null(rv$stats_results) || !is.null(rv$quant_results)
    active_idx <- if (!isTRUE(rv$library_ready)) 1L
      else if (!step2_done) 2L else if (!step3_done) 3L else if (!step4_done) 4L else 5L
    labels <- c("1. Library", "2. MS2 Library", "3. Batch Processing", "4. Statistics", "5. Complete")
    tags$div(class = "stepper",
      tagList(lapply(seq_along(labels), function(i) {
        tagList(
          tags$div(class = paste("wf-step", if (i == active_idx) "wf-active" else ""),
            tags$div(class = "wf-circ", if (i < active_idx) "✓" else as.character(i)),
            tags$div(class = "wf-lbl", labels[i])),
          if (i < length(labels)) tags$span(class = "wf-arrow", "→")
        )
      }))
    )
  })

  ## ---- Batch sample metadata table (group/timepoint/sample_type/concentration)
  # Seeded from the uploaded batch_files' original filenames; edited in place
  # via DT's edit feature, or in bulk via the sample info CSV upload below.
  # Blank group/timepoint columns mean "extract and match only, no
  # statistics" -- see the Run handler below. sample_type defaults to
  # "unknown" (one of its own controlled-vocabulary values, not blank) since
  # every row IS some kind of sample even before the user has classified it.
  # concentration is the nominal/spiked concentration for a calibration
  # standard or QC sample (blank for study "unknown" samples) -- stored as
  # character like every other column here so DT cell-edits and the CSV
  # merge don't need special-casing; parsed to numeric only where consumed
  # (fit_calibration_curve() in R/statistics.R).
  .SAMPLE_TYPE_LEVELS <- c("unknown", "standard", "quality_control",
                           "reagent_blank", "matrix_blank")

  # Case/spacing/punctuation-insensitive mapping of a typed sample_type
  # value against the controlled vocabulary above ("Quality Control"/"QC"/
  # "qc " all map to "quality_control") -- NA for anything that doesn't
  # map, so a caller can flag it rather than silently keep or drop it.
  # Shared by every place a sample_type value can reach batch_meta_data():
  # the CSV merge below AND direct cell-editing in the table (real-user
  # report: typing "Standard" into the table, instead of the exact lower-
  # case "standard" fit_calibration_curve()/.is_study_sample() match
  # against, silently broke calibration-curve fitting with no error --
  # the CSV path already normalized this, direct cell edits didn't).
  .map_sample_type <- function(x) {
    norm <- gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(x))))
    canon <- c(unknown = "unknown", standard = "standard",
               quality_control = "quality_control", qc = "quality_control",
               reagent_blank = "reagent_blank", blank = "reagent_blank",
               matrix_blank = "matrix_blank")
    unname(canon[norm])
  }

  batch_meta_data <- reactiveVal(data.frame(
    sample = character(0), group = character(0), timepoint = character(0),
    sample_type = character(0), concentration = character(0), stringsAsFactors = FALSE))

  observeEvent(input$batch_files, {
    samples <- tools::file_path_sans_ext(input$batch_files$name)
    batch_meta_data(data.frame(sample = samples, group = "", timepoint = "",
                                sample_type = "unknown", concentration = "",
                                stringsAsFactors = FALSE))
  })

  # Local-folder alternative to the upload above (see the UI section and
  # the shinyDirChoose wiring further down) -- skips the browser entirely,
  # so it has no upload-size limit. Triggers on any change to the batch_dir
  # text field, whether set via "Browse..." or typed directly.
  observeEvent(input$batch_dir, {
    dir <- input$batch_dir
    if (!is.null(dir) && nzchar(dir) && dir.exists(dir)) {
      paths <- list.files(dir, pattern = "\\.(mzML|mzXML|raw|wiff|baf|yep)$",
                          full.names = TRUE, ignore.case = TRUE)
      if (length(paths) > 0) {
        samples <- tools::file_path_sans_ext(basename(paths))
        batch_meta_data(data.frame(sample = samples, group = "", timepoint = "",
                                    sample_type = "unknown", concentration = "",
                                    stringsAsFactors = FALSE))
      }
    }
  }, ignoreInit = TRUE)

  # Unified batch-file source: a local folder (input$batch_dir) takes
  # precedence over the fileInput upload when both are set, since it's
  # strictly better when available (no size limit, no upload time). Always
  # returns the same shape fileInput does (a data.frame with $datapath and
  # $name), so downstream code (the Run handler) doesn't need to know which
  # source it came from.
  .batch_input_files <- function() {
    dir <- input$batch_dir
    if (!is.null(dir) && nzchar(dir) && dir.exists(dir)) {
      paths <- list.files(dir, pattern = "\\.(mzML|mzXML|raw|wiff|baf|yep)$",
                          full.names = TRUE, ignore.case = TRUE)
      if (length(paths) > 0) {
        return(data.frame(datapath = paths, name = basename(paths), stringsAsFactors = FALSE))
      }
    }
    input$batch_files
  }

  # Advisory only -- see .is_hosted_deployment()'s comment on why the actual
  # server-enforced cap can't be read from inside the app. Shows a general
  # note as soon as a hosted deployment is detected (before any files are
  # even chosen), and escalates to a size-specific warning once the
  # uploaded files' reported sizes (input$batch_files$size, in bytes) cross
  # a conservative heuristic threshold -- not a known limit, just a point
  # past which a rejection becomes plausible on a typical hosted setup.
  output$batch_upload_notice <- renderUI({
    hosted <- .is_hosted_deployment()
    files <- input$batch_files
    total_mb <- if (!is.null(files)) sum(files$size) / 1024^2 else NA_real_
    large <- !is.na(total_mb) && total_mb > 200
    no_python <- is.na(find_python())

    if (!hosted && !large) return(NULL)

    msg <- if (large) {
      sprintf(paste0(
        "⚠ Selected files total %.0f MB. %s can enforce its own upload ",
        "limit independent of this app's 20GB setting, and large HRMS files ",
        "may still be rejected. Running this app locally in RStudio (see ",
        "Quick Start in the README, or `git clone` + `shiny::runApp()`) ",
        "removes that limit entirely -- or use the local-folder input below ",
        "if the files already live on this machine."),
        total_mb, if (hosted) "This hosted deployment" else "A hosted Shiny deployment")
    } else {
      paste0(
        "⚠ This looks like a hosted deployment, which can enforce its ",
        "own upload limit independent of this app's 20GB setting. If a ",
        "large upload gets rejected, running this app locally in RStudio ",
        "removes that limit entirely -- or use the local-folder input below.")
    }
    if (hosted && no_python) {
      msg <- paste(msg,
        "Batch (multi-file, parallel) MS processing also needs a Python 3 ",
        "interpreter on PATH, which this hosted deployment does not have -- ",
        "batch processing will fail here. Run the app locally in RStudio ",
        "for batch MS processing.")
    }
    tags$p(style = "font-size: 11px; color: #a3231b; font-weight: 600; margin: 2px 0 6px;", msg)
  })

  # A timepoint value must be numeric-coercible or it's silently dropped
  # several steps downstream (compare_time_series()/quantify_relative() in
  # R/statistics.R both do as.numeric(timepoint) and filter out the NAs with
  # no error) -- flagging it right here, at entry, is the only place a typo
  # like "Day 1" is cheap to catch. .tp_invalid is appended as a hidden
  # extra column purely for DT's formatStyle to key off of; it sits AFTER
  # the 5 real columns, so it never shifts sample_meta_table_cell_edit's
  # $col indices for the columns a user can actually click into.
  .invalid_timepoints <- function(meta) {
    nzchar(meta$timepoint) & is.na(suppressWarnings(as.numeric(meta$timepoint)))
  }

  # A sample_type value that doesn't map onto the controlled vocabulary
  # (even after .map_sample_type()'s case/spacing/punctuation-insensitive
  # matching) breaks every exact sample_type == "standard"/"unknown" check
  # downstream (fit_calibration_curve(), .is_study_sample(), ...) with no
  # error -- same "flag it right here, at entry" reasoning as
  # .invalid_timepoints() above.
  .invalid_sample_types <- function(meta) {
    nzchar(meta$sample_type) & is.na(.map_sample_type(meta$sample_type))
  }

  output$sample_meta_table <- DT::renderDT({
    df <- batch_meta_data()
    df$.tp_invalid <- .invalid_timepoints(df)
    df$.st_invalid <- .invalid_sample_types(df)
    dt <- DT::datatable(df, editable = TRUE, rownames = FALSE,
                   options = list(dom = "t", paging = FALSE, scrollY = "180px",
                     columnDefs = list(list(visible = FALSE,
                       targets = which(names(df) %in% c(".tp_invalid", ".st_invalid")) - 1))))
    dt <- DT::formatStyle(dt, "timepoint", valueColumns = ".tp_invalid",
                     backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#fdecea", "white")),
                     color = DT::styleEqual(c(TRUE, FALSE), c("#a3231b", "inherit")),
                     fontWeight = DT::styleEqual(c(TRUE, FALSE), c("600", "normal")))
    DT::formatStyle(dt, "sample_type", valueColumns = ".st_invalid",
                     backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#fdecea", "white")),
                     color = DT::styleEqual(c(TRUE, FALSE), c("#a3231b", "inherit")),
                     fontWeight = DT::styleEqual(c(TRUE, FALSE), c("600", "normal")))
  })
  observeEvent(input$sample_meta_table_cell_edit, {
    df <- batch_meta_data()
    edit <- input$sample_meta_table_cell_edit
    val <- edit$value
    # sample_type is free-text-editable in this table, same as every other
    # cell -- canonicalize it the same way the CSV-merge path already does
    # (.map_sample_type()), so typing "Standard"/"QC"/"Quality Control"
    # lands on the exact lowercase controlled-vocabulary string
    # fit_calibration_curve()/.is_study_sample() match against, instead of
    # silently breaking calibration-curve fitting with no error. Anything
    # that still doesn't map is kept as typed and flagged red (see
    # .invalid_sample_types() above), the same "keep as typed, flag it"
    # policy the CSV path uses -- never silently coerced to "unknown".
    if (identical(names(df)[edit$col + 1], "sample_type") && nzchar(val)) {
      mapped <- .map_sample_type(val)
      if (!is.na(mapped)) val <- mapped
    }
    df[edit$row, edit$col + 1] <- val
    batch_meta_data(df)
  })

  # Live, human-readable companion to the red cell highlighting above --
  # updates on every cell edit and every CSV merge, since both write
  # through batch_meta_data().
  output$sample_meta_type_warn <- renderUI({
    meta <- batch_meta_data()
    if (nrow(meta) == 0) return(NULL)
    bad <- .invalid_sample_types(meta)
    if (!any(bad)) return(NULL)
    tags$p(style = "font-size: 11px; color: #a3231b; font-weight: 600; margin: 4px 0;",
      sprintf("Sample Type must be one of {%s} -- not free text. Fix: %s.",
              paste(.SAMPLE_TYPE_LEVELS, collapse = ", "),
              paste0(meta$sample[bad], " (\"", meta$sample_type[bad], "\")", collapse = ", ")))
  })

  # Live, human-readable companion to the red cell highlighting above --
  # updates on every cell edit and every CSV merge, same as the table does,
  # since both write through batch_meta_data().
  output$sample_meta_timepoint_warn <- renderUI({
    meta <- batch_meta_data()
    if (nrow(meta) == 0) return(NULL)
    bad <- .invalid_timepoints(meta)
    if (!any(bad)) return(NULL)
    tags$p(style = "font-size: 11px; color: #a3231b; font-weight: 600; margin: 4px 0;",
      sprintf("Timepoint must be numeric (e.g. 0, 4, 24) in one consistent unit -- not text. Fix: %s.",
              paste0(meta$sample[bad], " (\"", meta$timepoint[bad], "\")", collapse = ", ")))
  })

  ## ---- Sample info CSV: template download + bulk upload ---------------------
  # Template is pre-filled with whatever sample list is already on the
  # table (from uploaded files or a local folder), so round-tripping through
  # Excel doesn't require retyping/copy-pasting sample names by hand.
  output$dl_batch_meta_template <- downloadHandler(
    filename = function() "sample_info_template.csv",
    content = function(file) {
      df <- batch_meta_data()
      if (nrow(df) == 0) {
        df <- data.frame(sample = character(0), group = character(0),
                          timepoint = character(0), sample_type = character(0),
                          concentration = character(0))
      }
      utils::write.csv(df, file, row.names = FALSE)
    }
  )

  batch_meta_upload_status <- reactiveVal(NULL)
  output$batch_meta_upload_status <- renderUI({
    msg <- batch_meta_upload_status()
    if (is.null(msg)) return(NULL)
    color <- if (startsWith(msg, "WARNING")) "#a3231b" else "#2e7d32"
    tags$p(style = paste0("font-size: 11px; color: ", color, "; margin: -6px 0 6px;"), msg)
  })

  # Merges by `sample` name, not row position -- the sample list already on
  # the table (derived from the actual uploaded/local files) stays
  # authoritative for WHICH samples exist; the CSV only supplies values for
  # group/timepoint/sample_type, so it can be partial or reordered safely.
  observeEvent(input$batch_meta_csv, {
    current <- batch_meta_data()
    if (nrow(current) == 0) {
      batch_meta_upload_status(
        "WARNING: upload batch files (or set a local folder) first, so there's a sample list to merge the CSV onto.")
      return()
    }
    csv <- tryCatch(
      utils::read.csv(input$batch_meta_csv$datapath, stringsAsFactors = FALSE, colClasses = "character"),
      error = function(e) NULL)
    if (is.null(csv)) {
      batch_meta_upload_status("WARNING: could not read that file as CSV.")
      return()
    }
    names(csv) <- tolower(trimws(names(csv)))
    if (!"sample" %in% names(csv)) {
      batch_meta_upload_status("WARNING: CSV needs a 'sample' column matching the uploaded file names.")
      return()
    }
    csv$sample <- trimws(csv$sample)

    # Canonicalize sample_type against the controlled vocabulary (see
    # .map_sample_type() above) so a human-typed CSV isn't rejected over
    # formatting differences. Anything that still doesn't map is kept as
    # typed and flagged, rather than silently coerced to "unknown".
    invalid_types <- character(0)
    if ("sample_type" %in% names(csv)) {
      mapped <- .map_sample_type(csv$sample_type)
      unmapped <- is.na(mapped) & nzchar(csv$sample_type)
      invalid_types <- unique(csv$sample_type[unmapped])
      csv$sample_type[!is.na(mapped)] <- mapped[!is.na(mapped)]
    }

    # concentration must parse as a plain number (blank stays blank, meaning
    # "not a standard/QC" or "no target set") -- a non-numeric entry is kept
    # as typed and flagged, same treatment as an unmapped sample_type, since
    # fit_calibration_curve() would otherwise silently drop it with no
    # explanation of why that sample never made it onto the curve.
    invalid_conc <- character(0)
    if ("concentration" %in% names(csv)) {
      non_numeric <- nzchar(csv$concentration) & is.na(suppressWarnings(as.numeric(csv$concentration)))
      invalid_conc <- unique(csv$concentration[non_numeric])
    }

    # Same treatment as concentration above -- kept as typed (so the CSV
    # merge doesn't silently drop the value), flagged here in the summary,
    # and flagged again live via .invalid_timepoints()/the red cell
    # highlighting once it lands on the table.
    invalid_time <- character(0)
    if ("timepoint" %in% names(csv)) {
      non_numeric <- nzchar(csv$timepoint) & is.na(suppressWarnings(as.numeric(csv$timepoint)))
      invalid_time <- unique(csv$timepoint[non_numeric])
    }

    matched <- intersect(current$sample, csv$sample)
    unmatched_csv <- setdiff(csv$sample, current$sample)
    for (col in intersect(c("group", "timepoint", "sample_type", "concentration"), names(csv))) {
      for (s in matched) {
        val <- csv[[col]][csv$sample == s][1]
        if (!is.na(val) && nzchar(val)) current[current$sample == s, col] <- val
      }
    }
    batch_meta_data(current)

    msg <- sprintf("Sample info CSV applied: %d/%d uploaded samples matched.",
                    length(matched), nrow(current))
    if (length(unmatched_csv) > 0) {
      msg <- paste0(msg, " Ignored ", length(unmatched_csv),
                    " CSV row(s) with no matching uploaded file: ",
                    paste(unmatched_csv, collapse = ", "), ".")
    }
    if (length(invalid_types) > 0) {
      msg <- paste0(msg, " Sample Type value(s) not in {",
                    paste(.SAMPLE_TYPE_LEVELS, collapse = ", "),
                    "} kept as typed: ", paste(invalid_types, collapse = ", "), ".")
    }
    if (length(invalid_conc) > 0) {
      msg <- paste0(msg, " Concentration value(s) not numeric, kept as typed: ",
                    paste(invalid_conc, collapse = ", "), ".")
    }
    if (length(invalid_time) > 0) {
      msg <- paste0(msg, " Timepoint value(s) not numeric, kept as typed: ",
                    paste(invalid_time, collapse = ", "), ".")
    }
    batch_meta_upload_status(msg)
  })

  ## ---- Quantification controls: keep dropdown choices in sync --------------
  # Absolute-quant metabolite choices come from the library just built
  # (rv$mets), not from any matches yet -- the selection has to exist
  # BEFORE the batch run that will act on it. Re-populates whenever a new
  # library is generated; existing selections that are still valid met_ids
  # are preserved across a re-generation rather than silently cleared.
  observeEvent(rv$mets, {
    mets <- rv$mets
    if (is.null(mets) || length(mets) == 0) {
      updateSelectizeInput(session, "absolute_quant_mets", choices = character(0))
      return()
    }
    ids <- vapply(mets, function(m) m$id, character(1))
    labels <- vapply(mets, function(m) paste0(m$name, " (", m$id, ")"), character(1))
    choices <- stats::setNames(ids, labels)
    keep_selected <- intersect(isolate(input$absolute_quant_mets), ids)
    updateSelectizeInput(session, "absolute_quant_mets", choices = choices, selected = keep_selected)
  }, ignoreNULL = FALSE)

  # Control-group choices come from whatever Group values are actually on
  # the sample metadata table right now, so the dropdown never offers a
  # group that doesn't exist -- updates live as the DT table is edited or a
  # CSV is merged in.
  observeEvent(batch_meta_data(), {
    meta <- batch_meta_data()
    groups <- if ("group" %in% names(meta)) unique(meta$group[nzchar(meta$group)]) else character(0)
    keep_selected <- if (isolate(input$control_group) %in% groups) isolate(input$control_group) else NULL
    updateSelectInput(session, "control_group", choices = groups, selected = keep_selected)
  }, ignoreNULL = FALSE)

  ## ---- MS2 library explorer --------------------------------------------------
  # Cached separately from rv$ms_results etc.: building the full MS2 library
  # is expensive (see .ms2_library() below), so it's built once on demand via
  # the "Build / Refresh Library" button rather than on every Run, and stays
  # cached until the next Run invalidates it or the user rebuilds it.
  ms2_lib_rv <- reactiveVal(NULL)

  ## ---- Custom chemistry table ----------------------------------------------
  custom_chem_data <- reactiveVal(.custom_chem_init)

  ## ---- Manual three-line sequence entry ------------------------------------
  # Bases / sugars / linkages as three plain strings -- the layout a chemical
  # analysis file or a BioPharma Finder sequence entry uses. Submit converts
  # them to triplet notation and loads that into the sidebar's sequence box,
  # so the rest of the app sees an ordinary sequence and nothing downstream
  # needs to know which way it was entered.
  man_spec <- reactiveVal(NULL)

  observeEvent(input$man_example, {
    # The 18-mer worked example from inst/help/SEQUENCE_GUIDE.md (nusinersen:
    # uniform 2'-MOE, fully phosphorothioate, 5-methyl-C written as S).
    updateTextInput(session, "man_bases", value = "TSASTTTSATAATGSTGG")
    updateTextInput(session, "man_sugars", value = "eeeeeeeeeeeeeeeeee")
    updateTextInput(session, "man_linkages", value = "sssssssssssssssss")
  })

  observeEvent(input$man_clear, {
    for (id in c("man_bases", "man_sugars", "man_linkages")) {
      updateTextInput(session, id, value = "")
    }
    man_spec(NULL)
    output$man_feedback <- renderUI(NULL)
  })

  observeEvent(input$man_submit, {
    fields <- lapply(c("man_bases", "man_sugars", "man_linkages"),
                     function(id) trimws(input[[id]] %||% ""))
    names(fields) <- c("bases", "sugars", "linkages")
    empty <- names(fields)[!nzchar(unlist(fields))]
    if (length(empty) > 0) {
      man_spec(NULL)
      output$man_feedback <- renderUI(tags$p(class = "man-err",
        paste0("Fill in all three fields -- still empty: ",
               paste(empty, collapse = ", "), ".")))
      return()
    }

    # Build against the same dictionary the run will use, so a custom code
    # typed into the Custom Chemistry table validates here too.
    dict <- tryCatch(build_dictionary(overrides = .overrides_from_table(
      custom_chem_data())), error = function(e) STANDARD_DICT)

    spec <- tryCatch(
      parse_three_line(fields$bases, fields$sugars, fields$linkages,
                       conj5 = input$conj5, conj3 = input$conj3, dict = dict),
      error = function(e) structure(conditionMessage(e), class = "man_error"))

    if (inherits(spec, "man_error")) {
      man_spec(NULL)
      output$man_feedback <- renderUI(tags$p(class = "man-err",
        paste0("Could not build the sequence: ", as.character(spec))))
      return()
    }

    triplet <- format_triplet(spec)
    man_spec(spec)
    updateTextAreaInput(session, "seq", value = triplet)

    info <- tryCatch(metabolite_mass_info(spec, dict), error = function(e) NULL)
    output$man_feedback <- renderUI(tagList(
      tags$p(class = "man-ok",
             sprintf("Built a %d-mer (%d linkages) and loaded it into the sequence box. Click \"1. Generate Library\".",
                     spec$n, spec$n - 1L)),
      tags$p(class = "man-seq", triplet),
      if (!is.null(info)) tags$p(class = "hint",
        sprintf("Formula %s -- monoisotopic %.4f Da, average %.2f Da.",
                info$formula_str, info$mono_mass, info$avg_mass))
    ))
  })

  output$dl_fasta <- downloadHandler(
    filename = function() paste0(input$oligo_name %||% "oligo", ".fasta"),
    content = function(file) {
      # Prefer whatever is in the sequence box, so this works whether the
      # sequence came from manual entry, an example, or was typed directly.
      spec <- man_spec()
      if (is.null(spec)) {
        spec <- tryCatch(parse_input(trimws(input$seq %||% "")),
                         error = function(e) NULL)
        if (!is.null(spec)) {
          spec$conj5 <- input$conj5
          spec$conj3 <- input$conj3
        }
      }
      validate(need(!is.null(spec),
                    "Enter a valid sequence before downloading a FASTA."))
      writeLines(format_biopharma_fasta(spec, input$oligo_name %||% "oligo"),
                 file)
    }
  )

  # Load a selected example sequence into the input box. Terminal
  # conjugates travel with the example (triplet notation can't express
  # them), so the conjugate dropdowns are updated alongside the sequence --
  # otherwise loading the GalNAc-siRNA example would silently drop its
  # GalNAc and compute the wrong parent mass.
  observeEvent(input$load_example, {
    sel <- input$example_seq
    if (nzchar(sel) && sel %in% names(.EXAMPLE_SEQS)) {
      ex <- .EXAMPLE_SEQS[[sel]]
      updateTextAreaInput(session, "seq", value = ex$seq)
      updateSelectInput(session, "conj5", selected = ex$conj5 %||% "none")
      updateSelectInput(session, "conj3", selected = ex$conj3 %||% "none")
    }
  })

  # Folder picker for "Save to folder", via shinyFiles. Deliberately not
  # using tcltk::tk_choose.dir()/utils::choose.dir() here: those open a
  # blocking native OS dialog, and Shiny runs on a single R thread -- if
  # that dialog fails to initialize properly (seen intermittently on
  # Windows depending on how the R process was spawned) it can freeze the
  # entire app, Run button included, not just the picker. shinyFiles
  # renders an HTML folder browser inside the Shiny UI instead, so it can't
  # block the R process the way a native dialog can.
  if (.have_shinyfiles) {
    volumes <- c(Home = path.expand("~"), shinyFiles::getVolumes()())
    shinyFiles::shinyDirChoose(input, "browse_output_dir", roots = volumes,
                               session = session)
    observeEvent(input$browse_output_dir, {
      sel <- input$browse_output_dir
      if (is.list(sel) && !is.null(sel$path)) {
        path <- shinyFiles::parseDirPath(volumes, sel)
        if (length(path) == 1 && nzchar(path)) {
          updateTextInput(session, "output_dir", value = path)
        }
      }
    })

    shinyFiles::shinyDirChoose(input, "browse_batch_dir", roots = volumes, session = session)
    observeEvent(input$browse_batch_dir, {
      sel <- input$browse_batch_dir
      if (is.list(sel) && !is.null(sel$path)) {
        path <- shinyFiles::parseDirPath(volumes, sel)
        if (length(path) == 1 && nzchar(path)) {
          updateTextInput(session, "batch_dir", value = path)
        }
      }
    })
  }

  output$custom_chem <- DT::renderDataTable({
    d <- custom_chem_data()
    disp <- d
    disp$Status <- .chem_row_status(d)
    # Status (column index 5, 0-based) is display-only -- disabling it here
    # keeps DT's cell-edit column indices for Code..Attach (0-4) unchanged,
    # so the cell_edit handler below needs no adjustment for it.
    DT::datatable(disp, rownames = FALSE,
                  editable = list(target = "cell", disable = list(columns = 5)),
                  options = list(dom = "t", paging = FALSE, ordering = FALSE,
                                 autoWidth = TRUE,
                                 columnDefs = list(list(width = "60px", targets = 0),
                                                   list(width = "100px", targets = 1),
                                                   list(width = "80px", targets = 3),
                                                   list(width = "110px", targets = 5))),
                  selection = "none") |>
      DT::formatStyle("Status", color = DT::JS(
        "value.indexOf(String.fromCharCode(10003)) === 0 ? '#18632f' : (value === '' ? '#adb5bd' : '#a3231b')"))
  })

  # Capture cell edits
  observeEvent(input$custom_chem_cell_edit, {
    info <- input$custom_chem_cell_edit
    d <- custom_chem_data()
    d[info$row, info$col + 1] <- info$value  # +1 because rownames=FALSE shifts col index
    custom_chem_data(d)
  }, priority = 1000)

  # Add row
  observeEvent(input$add_row, {
    d <- custom_chem_data()
    d <- rbind(d, data.frame(Code = "", Formula = "", Name = "",
                             Type = "base", Attach = "add",
                             stringsAsFactors = FALSE))
    custom_chem_data(d)
  })

  # Remove row (keep at least 1)
  observeEvent(input$remove_row, {
    d <- custom_chem_data()
    if (nrow(d) > 1) {
      d <- d[-nrow(d), , drop = FALSE]
      custom_chem_data(d)
    }
  })

  ## ---- Session save/load -----------------------------------------------------
  output$dl_session <- downloadHandler(
    filename = function() paste0(input$oligo_name %||% "oligomet", "_session.json"),
    content = function(file) {
      vals <- stats::setNames(lapply(.session_input_ids, function(id) input[[id]]),
                               .session_input_ids)
      session_obj <- list(
        app = "OligoMetProfiler", format_version = 1L,
        saved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        inputs = vals,
        custom_chemistry = custom_chem_data()
      )
      writeLines(jsonlite::toJSON(session_obj, auto_unbox = TRUE, pretty = TRUE,
                                  na = "null"), file)
    }
  )

  observeEvent(input$load_session_file, {
    req(input$load_session_file)
    parsed <- tryCatch(
      jsonlite::fromJSON(input$load_session_file$datapath, simplifyVector = TRUE),
      error = function(e) NULL)
    if (is.null(parsed) || is.null(parsed$inputs)) {
      rv$status_text <- paste0(
        "ERROR: could not read '", input$load_session_file$name, "' -- not a ",
        "valid OligoMet Profiler session .json file.\n")
      return()
    }
    for (id in names(parsed$inputs)) {
      upd <- .session_update[[id]]
      v <- parsed$inputs[[id]]
      if (!is.null(upd) && !is.null(v) && length(v) > 0) {
        tryCatch(upd(session, v), error = function(e) NULL)
      }
    }
    if (!is.null(parsed$custom_chemistry)) {
      cc <- tryCatch(as.data.frame(parsed$custom_chemistry, stringsAsFactors = FALSE),
                     error = function(e) NULL)
      need_cols <- c("Code", "Formula", "Name", "Type", "Attach")
      if (!is.null(cc) && all(need_cols %in% names(cc)))
        custom_chem_data(cc[, need_cols])
    }
    # A restored session only has form parameters (see the module header
    # comment above) -- any library from a PRIOR run in this same session is
    # now stale relative to the just-restored parameters, so the Phase 2
    # button must go back behind Phase 1 until it's re-run.
    rv$library_ready <- FALSE
    rv$status_text <- paste0(
      "Loaded session from '", input$load_session_file$name, "'. Review the ",
      "restored parameters, then click \"Generate Library\" (Library Generation ",
      "tab) to regenerate the library (fast, deterministic -- reproduces exactly ",
      "what the saved session would have produced). Re-upload any MS/batch files ",
      "(not saved in a session), then run batch processing again.\n")
  })

  ## ---- Analysis state save/load ----------------------------------------------
  # A session .json (above) only ever held UI parameters -- fast/deterministic
  # to reproduce, so there was never a reason to persist their OUTPUT. Batch
  # results/stats/quantification are a different story: a real study batch
  # can take many minutes (deconvolution + matching across every file), and
  # losing that on a browser refresh means re-running the whole thing. RDS
  # (not JSON) because this is arbitrary R objects -- data.frames and nested
  # lists round-trip natively, with binary+xz compression keeping the file
  # far smaller than an equivalent JSON would be, and because the JSON writer
  # has no clean way to serialize an lm() fit object at all (see below).
  #
  # Deliberately NOT included, both to keep the file small and because
  # neither round-trips usefully:
  #   - ms2_spectra / ms_results$ms2_results[[.]]$best_spec: raw per-hit peak
  #     lists, needed only to redraw a specific mirror plot on demand -- by
  #     far the largest contributor for a batch with MS2 confirmation on.
  #   - quant_results$calibration_curves[[.]]$model: the raw lm() fit object.
  #     Its intercept/slope/r_squared/points/note are already extracted into
  #     the same list entry (everything quant_absolute_table/
  #     calibration_curves_table actually display) -- the fit object itself
  #     carries its own environment and isn't needed to show any of that.
  .strip_for_analysis_state <- function(ms_results, batch_ms_results, quant_results) {
    if (!is.null(ms_results) && length(ms_results$ms2_results) > 0) {
      ms_results$ms2_results <- lapply(ms_results$ms2_results, function(r) {
        r$best_spec <- NULL
        r
      })
    }
    if (!is.null(batch_ms_results)) batch_ms_results$ms2_spectra <- NULL
    if (!is.null(quant_results) && length(quant_results$calibration_curves) > 0) {
      quant_results$calibration_curves <- lapply(quant_results$calibration_curves, function(c) {
        c$model <- NULL
        c
      })
    }
    list(ms_results = ms_results, batch_ms_results = batch_ms_results, quant_results = quant_results)
  }

  output$dl_analysis_state <- downloadHandler(
    filename = function() paste0(input$oligo_name %||% "oligomet", "_analysis_state.rds"),
    content = function(file) {
      stripped <- .strip_for_analysis_state(rv$ms_results, rv$batch_ms_results, rv$quant_results)
      state <- list(
        app = "OligoMetProfiler", format_version = 1L,
        saved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        spec = rv$spec, mets = rv$mets, dict = rv$dict, prm = rv$prm,
        batch_features = rv$batch_features, sample_meta = rv$sample_meta,
        ms_results = stripped$ms_results,
        batch_ms_results = stripped$batch_ms_results,
        stats_results = rv$stats_results, kind_stats_results = rv$kind_stats_results,
        quant_results = stripped$quant_results
      )
      saveRDS(state, file, compress = "xz")
    }
  )

  observeEvent(input$load_analysis_state_file, {
    req(input$load_analysis_state_file)
    state <- tryCatch(readRDS(input$load_analysis_state_file$datapath), error = function(e) NULL)
    if (is.null(state) || is.null(state$format_version)) {
      rv$status_text <- paste0(
        "ERROR: could not read '", input$load_analysis_state_file$name, "' -- not a ",
        "valid OligoMet Profiler analysis state .rds file.\n")
      return()
    }
    rv$spec <- state$spec; rv$mets <- state$mets; rv$dict <- state$dict; rv$prm <- state$prm
    rv$batch_features <- state$batch_features
    rv$sample_meta <- state$sample_meta
    rv$ms_results <- state$ms_results
    rv$batch_ms_results <- state$batch_ms_results
    rv$stats_results <- state$stats_results
    rv$kind_stats_results <- state$kind_stats_results
    rv$quant_results <- state$quant_results
    if (!is.null(state$sample_meta)) batch_meta_data(state$sample_meta)
    rv$library_ready <- !is.null(rv$mets)
    rv$ready <- !is.null(rv$mets)
    rv$status_text <- paste0(
      "Loaded analysis state from '", input$load_analysis_state_file$name,
      "' (saved ", state$saved_at %||% "unknown time", "). Library, batch results, ",
      "statistics, and quantification are restored -- no need to re-run unless ",
      "you're adding new data. Mirror plots and the empirical MS2 library need ",
      "MS2 spectra, which this file doesn't carry -- re-run with \"Confirm hits ",
      "with MS2\" on if you need those back. Raw mzML/raw files are never stored ",
      "here either way.\n")
  })

  ## ---- About ---------------------------------------------------------------
  output$about_version <- renderText({
    v <- tryCatch(as.character(utils::packageVersion("OligoMetProfiler")),
                  error = function(e) NA_character_)
    if (is.na(v)) "running from a repository checkout" else paste("version", v)
  })

  ## ---- Help documents ------------------------------------------------------
  # Rendered lazily on first view and then cached by Shiny, so the markdown
  # is not converted on every page load of a collapsed panel.
  output$help_quickstart    <- renderUI(.help_ui("QUICKSTART.md"))
  output$help_quickstart_cli <- renderUI(.help_ui("QUICKSTART_CLI.md"))
  output$help_sequence      <- renderUI(.help_ui("SEQUENCE_GUIDE.md"))
  output$help_modifications <- renderUI(.help_ui("MODIFICATIONS.md"))

  ## ---- Status output -------------------------------------------------------
  output$status <- renderText({ rv$status_text })

  # Hidden flag for conditionalPanel
  output$status_ready <- reactive({ if (rv$ready) "true" else "false" })
  outputOptions(output, "status_ready", suspendWhenHidden = FALSE)

  # Hidden flag gating the Phase 2 button -- TRUE once Phase 1 has produced
  # a library (rv$spec/mets/dict/prm) this session.
  output$library_ready <- reactive({ if (rv$library_ready) "true" else "false" })
  outputOptions(output, "library_ready", suspendWhenHidden = FALSE)

  ## ---- Run pipeline --------------------------------------------------------
  ## ---- Shared: save copies directly to a user-chosen local folder ----------
  # Called at the end of BOTH Phase 1 (library-only) and Phase 2 (rebuilds
  # the workbook/report with MS results, so the local-folder copy should be
  # refreshed too). PRM/acquisition-method CSVs and spectral libraries only
  # ever depend on mets/dict/prm, not on whether MS data has been processed
  # yet, so this one helper covers both phases unmodified.
  .save_to_local_folder <- function(mets, dict, prm, z_range) {
    out_dir <- trimws(input$output_dir %||% "")
    if (!nzchar(out_dir)) return(NULL)
    tryCatch({
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      if (!dir.exists(out_dir)) stop("could not create folder")
      wb_dest <- file.path(out_dir, paste0(input$output_prefix, "_library.xlsx"))
      rp_dest <- file.path(out_dir, paste0(input$output_prefix, "_report.html"))
      prm_dest <- file.path(out_dir, paste0(input$output_prefix, "_prm_list.csv"))
      file.copy(rv$wb_path, wb_dest, overwrite = TRUE)
      if (!is.null(rv$report_path) && file.exists(rv$report_path)) {
        file.copy(rv$report_path, rp_dest, overwrite = TRUE)
      }
      utils::write.csv(prm, prm_dest, row.names = FALSE)
      ms1_dest <- file.path(out_dir, paste0(input$output_prefix, "_MS1_inclusion_list.csv"))
      ms2_dest <- file.path(out_dir, paste0(input$output_prefix, "_MS2_PRM_target_list.csv"))
      frag_dest <- file.path(out_dir, paste0(input$output_prefix, "_MS2_fragment_reference.csv"))
      z2 <- min(input$ms2_z_min, input$ms2_z_max):max(input$ms2_z_min, input$ms2_z_max)
      utils::write.csv(
        thermo_ms1_inclusion_list(mets, dict, z_range = z_range,
          h_offset = input$h_offset, max_oxid = input$max_oxid,
          rt_end = input$method_length, max_targets = input$ms1_target_cap),
        ms1_dest, row.names = FALSE)
      utils::write.csv(
        thermo_ms2_prm_target_list(mets, dict, z_range = z2,
          h_offset = input$h_offset, max_oxid = input$max_oxid,
          rt_end = input$method_length, nce = input$hcd_nce,
          max_targets = input$ms2_target_cap),
        ms2_dest, row.names = FALSE)
      utils::write.csv(ms2_fragment_reference(mets, dict, z_range = z2),
                       frag_dest, row.names = FALSE)
      export_spectral_libraries(
        mets, dict, out_dir = out_dir, prefix = input$output_prefix,
        z_range = z_range, n_iso = input$n_iso, max_oxid = input$max_oxid,
        precursor_z_range = z2, frag_z_range = 1:input$frag_z_max,
        h_offset = input$h_offset, use_envipat = input$use_envipat,
        oligo_name = input$oligo_name)
      normalizePath(out_dir)
    }, error = function(e) {
      paste0("WARNING: could not save to '", out_dir, "': ", conditionMessage(e))
    })
  }

  ## ---- Phase 1: generate the predicted library ------------------------------
  # Steps 1-6 (dictionary, sequence parsing, metabolite generation, masses,
  # fragments, PRM list) plus a library-only workbook/report build -- no MS
  # data needed. Sets rv$library_ready, which gates the Phase 2 button.
  observeEvent(input$run_phase1, {
    # Reset
    rv$ready <- FALSE
    rv$library_ready <- FALSE
    rv$status_text <- "Running...\n"
    ms2_lib_rv(NULL)

    # Validate inputs
    seq_str <- trimws(input$seq)
    if (nchar(seq_str) == 0) {
      rv$status_text <- "ERROR: Sequence is empty. Enter a triplet or OligoDistiller sequence.\n"
      return()
    }
    if (input$z_min >= input$z_max) {
      rv$status_text <- "ERROR: Min charge must be less than max charge.\n"
      return()
    }
    if (input$max_3p < 0 || input$max_5p < 0) {
      rv$status_text <- "ERROR: Truncation counts cannot be negative.\n"
      return()
    }
    chem_errors <- .chem_table_errors(custom_chem_data())
    if (length(chem_errors) > 0) {
      rv$status_text <- paste0(
        "ERROR: Custom Chemistry table has invalid row(s):\n  ",
        paste(chem_errors, collapse = "\n  "),
        "\nFix or clear these rows before running -- an invalid formula ",
        "would otherwise silently produce a zero-mass entry.\n")
      return()
    }

    # Wrap the whole phase in a top-level tryCatch: any error anywhere below
    # (not just the points with their own tryCatch) gets reported here
    # instead of silently stopping the reactive with nothing visible in the
    # status panel.
    tryCatch({

    withProgress(message = "Generating library...", value = 0, {

      prog <- progress_tracker(c(
        "Building dictionary" = 1,
        "Parsing sequence" = 1,
        "Generating metabolites" = 2,
        "Computing masses" = 1,
        "Generating fragment ions" = 2,
        "Generating PRM list" = 1,
        "Building Excel workbook" = 15,
        "Building report" = 4,
        "Generating plots" = 2
      ))

      # Step 1: Build dictionary with custom overrides
      incProgress(0.05, detail = "Building dictionary")
      progress_next(prog)
      dict <- build_dictionary(overrides =
                                 .overrides_from_table(custom_chem_data()))
      rv$dict <- dict

      # Step 2: Parse input
      incProgress(0.10, detail = "Parsing sequence")
      progress_next(prog)
      spec <- tryCatch(
        parse_input(seq_str, dict = dict),
        error = function(e) {
          rv$status_text <- paste0("ERROR parsing sequence: ", conditionMessage(e), "\n")
          NULL
        }
      )
      if (is.null(spec)) return()
      # Override conjugates from dropdowns
      spec$conj5 <- input$conj5
      spec$conj3 <- input$conj3
      rv$spec <- spec

      # Step 3: Generate metabolites
      incProgress(0.20, detail = "Generating metabolites")
      progress_next(prog)
      met_opts <- list(
        oligo_name = input$oligo_name,
        max_3p = input$max_3p,
        max_5p = input$max_5p,
        endo = input$endo,
        endo_sites = input$endo_sites,
        min_frag_len = input$min_frag_len
      )
      mets <- tryCatch(
        generate_metabolites(spec, opts = met_opts, dict = dict),
        error = function(e) {
          rv$status_text <- paste0("ERROR generating metabolites: ", conditionMessage(e), "\n")
          NULL
        }
      )
      if (is.null(mets)) return()
      rv$mets <- mets

      # Step 4: Compute parent mass/formula
      incProgress(0.30, detail = "Computing masses")
      progress_next(prog)
      parent_info <- metabolite_mass_info(mets[[1]], dict)

      # Step 5: Generate fragments
      incProgress(0.40, detail = "Generating fragment ions")
      progress_next(prog)
      frags <- generate_fragments(mets[[1]], dict)
      ifrags <- generate_internal_fragments(mets[[1]], dict)
      n_frags <- length(frags) + length(ifrags)

      # Step 6: Generate PRM inclusion list
      incProgress(0.50, detail = "Generating PRM list")
      progress_next(prog)
      z_range <- input$z_min:input$z_max
      prm <- prm_inclusion_list(mets, dict, z_range = z_range,
                                 max_oxid = input$max_oxid,
                                 h_offset = input$h_offset)
      rv$prm <- prm

      # A fresh (or regenerated) library invalidates any MS/batch results
      # from a previous Phase 2 run in this session -- they were matched
      # against the OLD library and are no longer valid to show alongside it.
      rv$ms_results <- NULL
      rv$batch_features <- NULL
      rv$batch_ms_results <- NULL
      rv$sample_meta <- NULL
      rv$stats_results <- NULL
      rv$kind_stats_results <- NULL
      rv$quant_results <- NULL

      # Library-only workbook/report build (ms_results = NULL, no batch
      # results) -- this is the "save the library" deliverable; Phase 2
      # rebuilds both with MS-matching results once data is imported.
      incProgress(0.60, detail = "Building Excel workbook")
      progress_next(prog)
      build_opts <- list(
        z_range = z_range, n_iso = input$n_iso,
        max_oxid = input$max_oxid, h_offset = input$h_offset,
        use_envipat = input$use_envipat, include_internal = FALSE,
        max_3p = input$max_3p, max_5p = input$max_5p,
        endo = input$endo,
        adducts = if (input$enable_ms && !is.null(input$adducts)) input$adducts else c("H","Na","K","NH4"),
        ppm_tol = if (input$enable_ms) input$ppm_tol else 10
      )
      wb_file <- paste0(input$output_prefix, "_library.xlsx")
      wb_path <- tryCatch({
        build_workbook(spec, mets, dict, NULL, NULL, NULL,
                        build_opts, wb_file,  # build_workbook writes to a scratch dir; download handler reads rv$wb_path
                        progress = function(msg) incProgress(0, detail = msg),
                        console_tracker = prog)
      }, error = function(e) {
        rv$status_text <- paste0("ERROR building workbook: ", conditionMessage(e), "\n")
        NULL
      })
      if (is.null(wb_path)) return()
      rv$wb_path <- wb_path

      incProgress(0.75, detail = "Building report")
      progress_next(prog)
      report_file <- paste0(input$output_prefix, "_report")
      report_path <- tryCatch({
        build_report(spec, mets, dict, NULL, NULL, NULL,
                      build_opts, output_file = report_file,
                      output_format = "html", plot_dir = tempdir())
      }, error = function(e) {
        file.path(tempdir(), paste0(report_file, ".html"))  # fallback
      })
      rv$report_path <- report_path

      # Plots (theoretical only -- don't depend on MS data)
      incProgress(0.85, detail = "Generating plots")
      progress_next(prog)
      z_rep <- min(max(round(median(z_range)), z_range[1]), z_range[length(z_range)])
      plots <- tryCatch({
        list(
          envelope = plot_charge_envelope(mets[[1]], dict, z_range = z_range,
                                           h_offset = input$h_offset,
                                           max_oxid = min(input$max_oxid, 3)),
          truncation = plot_truncation_series(mets, dict),
          isotope = plot_isotope_pattern(mets[[1]], dict, z = z_rep,
                                          n_iso = max(input$n_iso, 8),
                                          h_offset = input$h_offset,
                                          use_envipat = input$use_envipat),
          oxidation = plot_oxidation_series(mets[[1]], dict,
                                             max_oxid = input$max_oxid)
        )
      }, error = function(e) {
        rv$status_text <- paste0(rv$status_text,
          "WARNING: Some plots failed: ", conditionMessage(e), "\n")
        list()
      })
      rv$plots <- plots

      # Optional: save copies directly to a user-chosen local folder.
      saved_to <- .save_to_local_folder(mets, dict, prm, z_range)

      # Done
      incProgress(1.0, detail = "Complete")
      progress_finish(prog)

      # Build status summary
      n_mets <- length(mets)
      n_prm <- nrow(prm)
      status <- paste0(
        "Library generated.\n",
        "  Formula: ", parent_info$formula_str, "\n",
        "  Mono mass: ", sprintf("%.6f Da", parent_info$mono_mass), "\n",
        "  Avg mass: ", sprintf("%.4f Da", parent_info$avg_mass), "\n",
        "  Length: ", spec$n, " nucleotides\n",
        "  Metabolites: ", n_mets, "\n",
        "  Fragment ions: ", n_frags, "\n",
        "  PRM entries: ", n_prm, "\n"
      )
      if (!is.null(saved_to)) {
        status <- paste0(status,
          if (startsWith(saved_to, "WARNING")) paste0("  ", saved_to, "\n")
          else paste0("  Saved to: ", saved_to, "\n"))
      }
      status <- paste0(status,
        "\nDownload the library now, or upload MS data above and click ",
        "\"2. Import & Process MS Data\" to identify metabolites.\n")
      # Prepend rather than overwrite: rv$status_text may already carry
      # WARNING lines appended earlier in this same run (e.g. some plots
      # failed) -- losing those here would hide a real problem behind a
      # clean-looking "Library generated" message.
      warnings_so_far <- sub("^Running\\.\\.\\.\n", "", rv$status_text)
      rv$status_text <- paste0(warnings_so_far, status)
      rv$ready <- TRUE
      rv$library_ready <- TRUE

    })  # end withProgress

    }, error = function(e) {
      rv$status_text <- paste0(
        "ERROR: unexpected failure -- ", conditionMessage(e), "\n",
        "Check the R console for a full traceback.\n")
      rv$ready <- FALSE
      rv$library_ready <- FALSE
    })
  })  # end observeEvent(input$run_phase1)

  ## ---- Comparison stats + quantification (Statistical Analysis tab) --------
  # Re-reads the LIVE batch sample metadata table (batch_meta_data(), the
  # thing the user actually edits -- Group/Timepoint/Sample Type/
  # concentration) rather than rv$sample_meta, which is otherwise only ever
  # set once, at the end of the (slow, 30+ minute) full batch MS run --
  # editing Sample Type/concentration afterward to mark a calibration curve
  # and clicking "Run Statistical Analysis" used to silently do nothing
  # because this read the stale snapshot instead. Computes
  # rv$stats_results/rv$kind_stats_results/rv$quant_results/
  # rv$batch_ms_results$degradation from the CURRENT settings -- called
  # automatically right after a batch run finishes, and again from "Run
  # Statistical Analysis" when the user only wants to change metadata / the
  # P-adjustment method / Group A-B override / calibration settings without
  # re-running the whole batch pipeline.
  .recompute_stats_and_quant <- function() {
    meta <- batch_meta_data()
    bres <- rv$batch_ms_results
    if (is.null(meta) || nrow(meta) == 0 || is.null(bres) || nrow(bres$ms1_matches) == 0) {
      return(invisible(NULL))
    }
    rv$sample_meta <- meta

    has_group <- !is.null(meta$group) && any(nzchar(meta$group))
    has_time <- !is.null(meta$timepoint) && any(nzchar(meta$timepoint))

    p_adjust_method <- input$stats_padjust %||% "BH"
    # An explicit Group A/B override (Statistical Analysis tab) takes
    # priority over the first two Group values found -- lets a 3+-group
    # design still get a specific pairwise 2-group test without retyping
    # the sample metadata.
    .pick_pair <- function(groups) {
      ga <- input$stats_group_a; gb <- input$stats_group_b
      if (!is.null(ga) && !is.null(gb) && nzchar(ga) && nzchar(gb) &&
          ga %in% groups && gb %in% groups && ga != gb) c(ga, gb) else groups[1:2]
    }

    abund <- build_abundance_matrix(bres$ms1_matches)
    stats_res <- tryCatch({
      if (has_time) {
        sm <- meta[, c("sample", "timepoint")]
        sm$timepoint <- suppressWarnings(as.numeric(sm$timepoint))
        long <- abundance_long(abund, sm[!is.na(sm$timepoint), ])
        list(mode = "time_series", result = compare_time_series(long, p_adjust_method = p_adjust_method))
      } else {
        sm <- meta[nzchar(meta$group), c("sample", "group")]
        groups <- unique(sm$group)
        long <- abundance_long(abund, sm)
        if (length(groups) == 2) {
          pair <- .pick_pair(groups)
          list(mode = "two_group", result = compare_two_groups(long, pair[1], pair[2], p_adjust_method = p_adjust_method))
        } else if (length(groups) > 2) {
          list(mode = "multi_group", result = compare_multi_groups(long, groups, p_adjust_method = p_adjust_method))
        } else NULL
      }
    }, error = function(e) {
      rv$status_text <- paste0(rv$status_text,
        "WARNING: statistical comparison failed: ", conditionMessage(e), "\n")
      NULL
    })
    rv$stats_results <- stats_res

    # Same comparison, but on kind-level (composition-class) totals instead
    # of per-metabolite abundance -- reuses the same compare_*() functions
    # unmodified (see build_kind_abundance_matrix() in R/statistics.R).
    signal_col <- if ("area" %in% names(bres$ms1_matches) &&
                       any(!is.na(bres$ms1_matches$area))) "area" else "intensity"
    kind_abund <- build_kind_abundance_matrix(bres$ms1_matches, signal_col = signal_col)
    kind_stats_res <- tryCatch({
      if (has_time) {
        sm <- meta[, c("sample", "timepoint")]
        sm$timepoint <- suppressWarnings(as.numeric(sm$timepoint))
        long <- kind_abundance_long(kind_abund, sm[!is.na(sm$timepoint), ])
        list(mode = "time_series", result = compare_time_series(long, p_adjust_method = p_adjust_method))
      } else {
        sm <- meta[nzchar(meta$group), c("sample", "group")]
        groups <- unique(sm$group)
        long <- kind_abundance_long(kind_abund, sm)
        if (length(groups) == 2) {
          pair <- .pick_pair(groups)
          list(mode = "two_group", result = compare_two_groups(long, pair[1], pair[2], p_adjust_method = p_adjust_method))
        } else if (length(groups) > 2) {
          list(mode = "multi_group", result = compare_multi_groups(long, groups, p_adjust_method = p_adjust_method))
        } else NULL
      }
    }, error = function(e) {
      rv$status_text <- paste0(rv$status_text,
        "WARNING: kind-level statistical comparison failed: ", conditionMessage(e), "\n")
      NULL
    })
    rv$kind_stats_results <- kind_stats_res

    # Absolute quantification (calibration curve) for the user-selected
    # metabolites, relative quantification (fold-change vs pre-dose/time-0,
    # or vs the Control group) for every other metabolite -- see
    # quantify_metabolites() in R/statistics.R. Uses the FULL sample_meta
    # (not the sample+group/timepoint slice above), since it also needs
    # sample_type/concentration for the calibration curve half.
    quant_res <- tryCatch({
      quantify_metabolites(
        bres$ms1_matches, meta,
        absolute_met_ids = input$absolute_quant_mets,
        mode = if (has_time) "time_series" else "group",
        control_group = if (has_time) NULL else input$control_group,
        weighting = input$calibration_weighting)
    }, error = function(e) {
      rv$status_text <- paste0(rv$status_text,
        "WARNING: quantification failed: ", conditionMessage(e), "\n")
      NULL
    })
    rv$quant_results <- quant_res

    # Degradation summary (per-sample %, composition-by-class) also needs
    # to reflect the CURRENT Sample Type/Group/Timepoint edits -- it's
    # computed once inside annotate_metabolites_batch() at batch-run time,
    # before the user has necessarily finished marking calibration
    # standards/QC, so it's recomputed here from the live metadata too
    # (excludes standards/QC/blanks, joins group/timepoint -- see
    # degradation_summary()'s sample_meta argument).
    bres$degradation <- tryCatch(
      degradation_summary(bres$ms1_matches, sample_meta = meta),
      error = function(e) {
        rv$status_text <- paste0(rv$status_text,
          "WARNING: degradation summary failed: ", conditionMessage(e), "\n")
        NULL
      })
    rv$batch_ms_results <- bres
    invisible(NULL)
  }

  observeEvent(input$run_stats, { .recompute_stats_and_quant() })

  # Group A/B override choices (Statistical Analysis tab) -- same source and
  # update pattern as control_group above, kept in sync live as the sample
  # metadata table is edited or a CSV is merged in.
  observeEvent(batch_meta_data(), {
    meta <- batch_meta_data()
    groups <- if ("group" %in% names(meta)) unique(meta$group[nzchar(meta$group)]) else character(0)
    ch <- c("Auto (first two found)" = "", groups)
    keep_a <- if (isolate(input$stats_group_a) %in% groups) isolate(input$stats_group_a) else NULL
    keep_b <- if (isolate(input$stats_group_b) %in% groups) isolate(input$stats_group_b) else NULL
    updateSelectInput(session, "stats_group_a", choices = ch, selected = keep_a)
    updateSelectInput(session, "stats_group_b", choices = ch, selected = keep_b)
  }, ignoreNULL = FALSE)

  output$stats_design_summary <- renderUI({
    meta <- rv$sample_meta
    if (is.null(meta)) {
      return(tags$p(style = "font-size: 10.5px; color: #6c757d;",
                    "Run Batch Processing with Group or Timepoint filled in to see the detected design here."))
    }
    has_time <- !is.null(meta$timepoint) && any(nzchar(meta$timepoint))
    mode_txt <- if (has_time) "Time-course (Timepoint filled)" else "Group comparison (Group filled)"
    tags$p(style = "font-size: 10.5px; color: #1f2430;",
           "Auto-detected: ", tags$b(mode_txt))
  })

  # Post-hoc, R-side charge-state re-aggregation (mirroring charge_group.py's
  # neutral-mass-consensus grouping) is NOT implemented yet -- this handler
  # says so plainly rather than leaving the button silently do nothing.
  cg_rerun_msg <- reactiveVal(NULL)
  observeEvent(input$cg_rerun, {
    cg_rerun_msg(paste0(
      "Post-hoc R-side charge-grouping (a charge_group.py port applied to ",
      "the already-matched batch features) is not implemented yet -- this ",
      "panel is a placeholder for that future work. RT tol=", input$cg_rt_tol,
      " min, mass tol=", input$cg_mass_tol_ppm, " ppm, min charge states=",
      input$cg_min_charge_states, " were noted but not applied."))
  })
  output$cg_rerun_status <- renderUI({
    msg <- cg_rerun_msg()
    if (is.null(msg)) return(NULL)
    tags$p(class = "placeholder-note", msg)
  })

  ## ---- Phase 2: import & process MS data ------------------------------------
  # Steps 7/7b (single-file and/or batch MS matching) against the library
  # Phase 1 built, then rebuilds the workbook/report with those results
  # included. Requires rv$library_ready (Phase 1 having run this session).
  # Extracted to a plain function so it can be triggered from either of two
  # buttons now that the UI is split across tabs: "Run Batch Processing"
  # (Batch Processing tab) and "Generate Empirical MS2 Library" (Empirical
  # MS2 Library tab) -- both run the exact same pipeline; which downstream
  # results end up populated just depends on which checkboxes/files are set
  # (enable_ms/enable_batch), same as before this tab split.
  .run_phase2_now <- function() {
    if (!isTRUE(rv$library_ready)) {
      rv$status_text <- "ERROR: generate the library first (\"1. Generate Library\").\n"
      return()
    }

    rv$ready <- FALSE
    rv$status_text <- "Running...\n"

    tryCatch({

    withProgress(message = "Importing and processing MS data...", value = 0, {

      prog <- progress_tracker(c(
        "MS data import and matching" = 3,
        "Batch MS processing (parallel)" = 10,
        "Building Excel workbook" = 15,
        "Building report" = 4
      ))

      # Re-derive the (cheap, deterministic) locals Phase 1 already computed
      # and stored on rv -- steps 7/7b need them, but each observeEvent is
      # its own closure so they aren't visible here directly.
      spec <- rv$spec
      mets <- rv$mets
      dict <- rv$dict
      prm <- rv$prm
      z_range <- input$z_min:input$z_max
      parent_info <- metabolite_mass_info(mets[[1]], dict)
      frags <- generate_fragments(mets[[1]], dict)
      ifrags <- generate_internal_fragments(mets[[1]], dict)
      n_frags <- length(frags) + length(ifrags)

      # Step 7: Optional MS matching
      progress_next(prog)
      ms_results <- NULL
      if (input$enable_ms && !is.null(input$ms_file)) {
        incProgress(0.10, detail = "Importing MS data")
        ms_results <- tryCatch({
          # Vendor .raw (or other vendor formats) is converted to .mzML via
          # msconvert first, if needed -- resolve_ms_input_file() passes
          # .mzML/.mzXML/.csv/.txt straight through unchanged.
          ms_path <- resolve_ms_input_file(input$ms_file$datapath)
          if (grepl("\\.mzML$|\\.mzml$|\\.mzXML$|\\.mzxml$", ms_path, ignore.case = TRUE)) {
            ms_data <- read_ms_file(ms_path)
            if (isTRUE(ms_data$info$profile_mode)) {
              rv$status_text <- paste0(rv$status_text,
                "WARNING: this mzML file appears to be PROFILE mode (not centroided). ",
                "MS1 feature extraction is designed for centroided peaks and may run slowly ",
                "and/or produce poor results -- convert with `msconvert --centroid` first for ",
                "reliable results.\n")
            }
            ms1_features <- if (input$noise_mode == "sn") {
              extract_ms1_features(ms_data$ms1, ppm = input$ppm_tol, sn_threshold = input$sn_threshold)
            } else {
              extract_ms1_features(ms_data$ms1, ppm = input$ppm_tol, min_intensity = input$min_intensity)
            }
            adducts <- input$adducts
            if (is.null(adducts)) adducts <- "H"
            annotate_metabolites(mets, ms1_features, ms_data$ms2,
              dict = dict, ppm_tol = input$ppm_tol, z_range = z_range,
              adducts = adducts, max_oxid = input$max_oxid,
              h_offset = input$h_offset, frag_tol_ppm = input$frag_tol_ppm,
              frag_z_range = 1:input$frag_z_max,
              n_iso = input$n_iso, use_envipat = input$use_envipat)
          } else {
            ms1_features <- import_peak_list(ms_path)
            adducts <- input$adducts
            if (is.null(adducts)) adducts <- "H"
            annotate_metabolites(mets, ms1_features, NULL,
              dict = dict, ppm_tol = input$ppm_tol, z_range = z_range,
              adducts = adducts, max_oxid = input$max_oxid,
              h_offset = input$h_offset,
              n_iso = input$n_iso, use_envipat = input$use_envipat)
          }
        }, error = function(e) {
          rv$status_text <- paste0(rv$status_text,
            "WARNING: MS matching failed: ", conditionMessage(e), "\n")
          NULL
        })
      }
      rv$ms_results <- ms_results

      # Step 7b: Optional batch MS processing (parallel, multi-file)
      progress_next(prog)
      rv$batch_features <- NULL
      rv$batch_ms_results <- NULL
      rv$sample_meta <- NULL
      rv$stats_results <- NULL
      rv$kind_stats_results <- NULL
      rv$quant_results <- NULL
      batch_files_df <- .batch_input_files()
      if (input$enable_batch && !is.null(batch_files_df) && nrow(batch_files_df) > 0) {
        incProgress(0.25, detail = "Batch deconvolution (parallel)")
        batch_out <- tryCatch({
          adducts <- input$adducts
          if (is.null(adducts)) adducts <- "H"

          watchlist_path <- NULL
          if (isTRUE(input$batch_run_ms2)) {
            watchlist_path <- tempfile(fileext = ".txt")
            write_precursor_watchlist(mets, dict, z_range = z_range,
                                       max_oxid = input$max_oxid, h_offset = input$h_offset,
                                       out_path = watchlist_path)
          }

          # Resolve each input file to an .mzML path first (vendor formats
          # go through the msconvert bridge). Each conversion is wrapped
          # individually so one bad/unconvertible file doesn't abort the
          # whole batch -- it's dropped with a warning instead.
          resolved <- vapply(seq_len(nrow(batch_files_df)), function(i) {
            tryCatch(resolve_ms_input_file(batch_files_df$datapath[i]),
                     error = function(e) {
                       rv$status_text <- paste0(rv$status_text,
                         "WARNING: skipping '", batch_files_df$name[i], "': ",
                         conditionMessage(e), "\n")
                       NA_character_
                     })
          }, character(1))
          ok <- !is.na(resolved)
          if (!any(ok)) stop("no batch files could be read (all failed vendor conversion)")
          resolved_paths <- resolved[ok]
          resolved_names <- batch_files_df$name[ok]

          deconv <- if (input$batch_noise_mode == "sn") {
            run_batch_deconvolution(
              resolved_paths, precursor_watchlist = watchlist_path,
              mass_tol_ppm = input$batch_deconv_ppm, n_workers = input$batch_n_workers,
              sn_threshold = input$batch_sn_threshold,
              progress = function(msg) incProgress(0, detail = msg))
          } else {
            run_batch_deconvolution(
              resolved_paths, precursor_watchlist = watchlist_path,
              mass_tol_ppm = input$batch_deconv_ppm, n_workers = input$batch_n_workers,
              min_intensity = input$batch_min_intensity,
              progress = function(msg) incProgress(0, detail = msg))
          }

          if (!is.null(deconv$profile_mode_files) && nrow(deconv$profile_mode_files) > 0) {
            rv$status_text <- paste0(rv$status_text,
              "WARNING: ", nrow(deconv$profile_mode_files), " batch file(s) appear to be PROFILE ",
              "mode (not centroided): ", paste(deconv$profile_mode_files$sample, collapse = ", "),
              ". ROI/charge-envelope detection is designed for centroided peaks and may run ",
              "slowly and/or produce poor results -- convert with `msconvert --centroid` first ",
              "for reliable results.\n")
          }

          # Only present in S/N noise mode -- the same sn_threshold multiple
          # produces a different absolute cutoff per file (each file has its
          # own noise level), so worth surfacing what actually got applied
          # rather than leaving it invisible inside the run.
          if (!is.null(deconv$noise_thresholds) && nrow(deconv$noise_thresholds) > 0) {
            nt <- deconv$noise_thresholds
            summary_lines <- sprintf("%s: noise=%.0f -> threshold=%.0f",
                                      nt$sample, nt$noise_level, nt$effective_min_intensity)
            rv$status_text <- paste0(rv$status_text,
              "S/N noise threshold applied per file (", nrow(nt), " x ", nt$sn_threshold[1], "):\n  ",
              paste(summary_lines, collapse = "\n  "), "\n")
          }

          # Shiny renames uploads to random tmp paths (local-folder paths
          # are already their own basename, so this is a no-op there);
          # recover the sample name from the ORIGINAL filename, not the
          # (possibly msconvert-generated) resolved datapath basename.
          name_map <- stats::setNames(tools::file_path_sans_ext(resolved_names),
                                       tools::file_path_sans_ext(basename(resolved_paths)))
          feats <- read_batch_features(deconv$features_path)
          feats$sample <- unname(name_map[feats$sample])
          ms2 <- if (!is.null(deconv$ms2_path)) read_batch_ms2(deconv$ms2_path) else NULL
          if (!is.null(ms2) && nrow(ms2) > 0) ms2$sample <- unname(name_map[ms2$sample])

          incProgress(0, detail = "Matching + MS2 confirmation")
          batch_results <- annotate_metabolites_batch(
            mets, feats, ms2, dict = dict, ppm_tol = input$ppm_tol, z_range = z_range,
            adducts = adducts, max_oxid = input$max_oxid, h_offset = input$h_offset,
            n_iso = input$n_iso, use_envipat = input$use_envipat,
            frag_tol_ppm = input$frag_tol_ppm, frag_z_range = 1:input$frag_z_max,
            sample_meta = batch_meta_data())

          list(features = feats, results = batch_results, meta = batch_meta_data())
        }, error = function(e) {
          # conditionMessage(e) carries the full captured Python stdout/stderr
          # (see run_batch_deconvolution()'s stop() call in
          # R/batch_ms_processing.R) -- printing it to the console, not just
          # into rv$status_text, means a batch failure is visible immediately
          # in the R/RStudio console instead of only in the app's Status
          # panel, which the rest of the pipeline (workbook/report/plots)
          # completes past silently since they all tolerate a NULL batch_out.
          message("Batch MS processing failed:\n", conditionMessage(e))
          rv$status_text <- paste0(rv$status_text,
            "WARNING: Batch MS processing failed: ", conditionMessage(e), "\n")
          NULL
        })

        if (!is.null(batch_out)) {
          rv$batch_features <- batch_out$features
          rv$batch_ms_results <- batch_out$results
          rv$sample_meta <- batch_out$meta

          .recompute_stats_and_quant()
        }
      }

      # Step 8: Rebuild the workbook, now WITH MS-matching results
      incProgress(0.60, detail = "Building Excel workbook")
      progress_next(prog)
      build_opts <- list(
        z_range = z_range, n_iso = input$n_iso,
        max_oxid = input$max_oxid, h_offset = input$h_offset,
        use_envipat = input$use_envipat, include_internal = FALSE,
        max_3p = input$max_3p, max_5p = input$max_5p,
        endo = input$endo,
        adducts = if (input$enable_ms && !is.null(input$adducts)) input$adducts else c("H","Na","K","NH4"),
        ppm_tol = if (input$enable_ms) input$ppm_tol else 10
      )
      wb_file <- paste0(input$output_prefix, "_library.xlsx")
      wb_path <- tryCatch({
        build_workbook(spec, mets, dict, NULL, ms_results,
                        if (input$enable_ms) list(n_ms1 = 0, n_ms2 = 0) else NULL,
                        build_opts, wb_file,  # build_workbook writes to a scratch dir; download handler reads rv$wb_path
                        progress = function(msg) incProgress(0, detail = msg),
                        console_tracker = prog,
                        batch_ms_results = rv$batch_ms_results,
                        stats_results = rv$stats_results,
                        kind_stats_results = rv$kind_stats_results)
      }, error = function(e) {
        rv$status_text <- paste0("ERROR building workbook: ", conditionMessage(e), "\n")
        NULL
      })
      if (is.null(wb_path)) return()
      rv$wb_path <- wb_path

      # Step 9: Rebuild the report
      incProgress(0.80, detail = "Building report")
      progress_next(prog)
      report_file <- paste0(input$output_prefix, "_report")
      report_path <- tryCatch({
        build_report(spec, mets, dict, NULL, ms_results,
                      if (input$enable_ms) list(n_ms1 = 0, n_ms2 = 0) else NULL,
                      build_opts, output_file = report_file,
                      output_format = "html", plot_dir = tempdir())
      }, error = function(e) {
        file.path(tempdir(), paste0(report_file, ".html"))  # fallback
      })
      rv$report_path <- report_path

      # Optional: refresh copies in the user-chosen local folder.
      saved_to <- .save_to_local_folder(mets, dict, prm, z_range)

      # Done
      incProgress(1.0, detail = "Complete")
      progress_finish(prog)

      # Build status summary
      n_mets <- length(mets)
      n_prm <- nrow(prm)
      status <- paste0(
        "Pipeline complete.\n",
        "  Formula: ", parent_info$formula_str, "\n",
        "  Mono mass: ", sprintf("%.6f Da", parent_info$mono_mass), "\n",
        "  Avg mass: ", sprintf("%.4f Da", parent_info$avg_mass), "\n",
        "  Length: ", spec$n, " nucleotides\n",
        "  Metabolites: ", n_mets, "\n",
        "  Fragment ions: ", n_frags, "\n",
        "  PRM entries: ", n_prm, "\n"
      )
      if (!is.null(ms_results) && !is.null(ms_results$summary) && nrow(ms_results$summary) > 0) {
        n_ms1 <- nrow(ms_results$ms1_matches)
        n_annot <- nrow(ms_results$summary)
        n_conf <- sum(ms_results$summary$confident, na.rm = TRUE)
        status <- paste0(status,
          "  MS1 matches: ", n_ms1, "\n",
          "  Annotated metabolites: ", n_annot, "\n",
          "  Putative IDs (ppm match): ", n_conf, "\n")
      }
      if (!is.null(saved_to)) {
        status <- paste0(status,
          if (startsWith(saved_to, "WARNING")) paste0("  ", saved_to, "\n")
          else paste0("  Saved to: ", saved_to, "\n"))
      }
      # Prepend rather than overwrite: rv$status_text may already carry
      # WARNING lines appended earlier in this same run (MS matching
      # failures, skipped batch files, profile-mode notices, ...) -- losing
      # those here would silently hide a real problem behind a clean-looking
      # "Pipeline complete" message.
      warnings_so_far <- sub("^Running\\.\\.\\.\n", "", rv$status_text)
      rv$status_text <- paste0(warnings_so_far, status)
      rv$ready <- TRUE

    })  # end withProgress

    }, error = function(e) {
      rv$status_text <- paste0(
        "ERROR: unexpected failure -- ", conditionMessage(e), "\n",
        "Check the R console for a full traceback.\n")
      rv$ready <- FALSE
    })
  }  # end .run_phase2_now()
  observeEvent(input$run_phase2, { .run_phase2_now() })
  observeEvent(input$run_ms2_library, { .run_phase2_now() })

  ## ---- Summary metrics outputs ---------------------------------------------
  output$m_formula <- renderText({
    if (rv$ready && !is.null(rv$mets)) {
      metabolite_mass_info(rv$mets[[1]], rv$dict)$formula_str
    }
  })
  output$m_mono_mass <- renderText({
    if (rv$ready && !is.null(rv$mets))
      sprintf("%.4f", metabolite_mass_info(rv$mets[[1]], rv$dict)$mono_mass)
  })
  output$m_avg_mass <- renderText({
    if (rv$ready && !is.null(rv$mets))
      sprintf("%.4f", metabolite_mass_info(rv$mets[[1]], rv$dict)$avg_mass)
  })
  output$m_length <- renderText({
    if (rv$ready && !is.null(rv$spec)) as.character(rv$spec$n)
  })
  output$m_n_mets <- renderText({
    if (rv$ready && !is.null(rv$mets)) as.character(length(rv$mets))
  })
  output$m_n_prm <- renderText({
    if (rv$ready && !is.null(rv$prm)) as.character(nrow(rv$prm))
  })
  output$m_ms1_matches <- renderText({
    if (rv$ready && !is.null(rv$ms_results) && !is.null(rv$ms_results$ms1_matches))
      as.character(nrow(rv$ms_results$ms1_matches)) else "0"
  })
  output$m_annotated <- renderText({
    if (rv$ready && !is.null(rv$ms_results) && !is.null(rv$ms_results$summary))
      as.character(nrow(rv$ms_results$summary)) else "0"
  })
  output$m_confident <- renderText({
    if (rv$ready && !is.null(rv$ms_results) && !is.null(rv$ms_results$summary))
      as.character(sum(rv$ms_results$summary$confident, na.rm = TRUE)) else "0"
  })

  ## ---- Empirical MS2 Library tab metrics ------------------------------------
  output$ms2lib_n_matched <- renderText({
    ms2c <- rv$batch_ms_results$ms2_confirmations
    if (!is.null(ms2c)) as.character(nrow(ms2c)) else "0"
  })
  output$ms2lib_n_spectra <- renderText({
    sp <- rv$batch_ms_results$ms2_spectra
    if (!is.null(sp)) as.character(length(sp)) else "0"
  })
  output$ms2lib_n_consensus <- renderText({
    tryCatch(as.character(nrow(.empirical_ms2_library()$summary)),
             error = function(e) "0")
  })

  ## ---- Plot outputs --------------------------------------------------------
  output$plot_envelope <- renderPlot({
    if (rv$ready && !is.null(rv$plots$envelope)) rv$plots$envelope
  })
  output$plot_truncation <- renderPlot({
    if (rv$ready && !is.null(rv$plots$truncation)) rv$plots$truncation
  })
  output$plot_isotope <- renderPlot({
    if (rv$ready && !is.null(rv$plots$isotope)) rv$plots$isotope
  })
  output$plot_oxidation <- renderPlot({
    if (rv$ready && !is.null(rv$plots$oxidation)) rv$plots$oxidation
  })

  ## ---- Download handlers ---------------------------------------------------
  output$dl_workbook <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_library.xlsx"),
    content = function(file) {
      if (!is.null(rv$wb_path) && file.exists(rv$wb_path))
        file.copy(rv$wb_path, file)
    }
  )
  output$dl_report <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_report.html"),
    content = function(file) {
      if (!is.null(rv$report_path) && file.exists(rv$report_path))
        file.copy(rv$report_path, file)
    }
  )
  output$dl_prm <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_prm_list.csv"),
    content = function(file) {
      if (!is.null(rv$prm)) write.csv(rv$prm, file, row.names = FALSE)
    }
  )
  output$dl_ms1_inclusion <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_MS1_inclusion_list.csv"),
    content = function(file) {
      req(rv$ready, rv$mets, rv$dict)
      lst <- thermo_ms1_inclusion_list(
        rv$mets, rv$dict, z_range = input$z_min:input$z_max,
        h_offset = input$h_offset, max_oxid = input$max_oxid,
        rt_start = 0, rt_end = input$method_length,
        max_targets = input$ms1_target_cap)
      write.csv(lst, file, row.names = FALSE)
    }
  )
  output$dl_ms2_prm <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_MS2_PRM_target_list.csv"),
    content = function(file) {
      req(rv$ready, rv$mets, rv$dict)
      z2 <- min(input$ms2_z_min, input$ms2_z_max):max(input$ms2_z_min, input$ms2_z_max)
      lst <- thermo_ms2_prm_target_list(
        rv$mets, rv$dict, z_range = z2,
        h_offset = input$h_offset, max_oxid = input$max_oxid,
        rt_start = 0, rt_end = input$method_length,
        nce = input$hcd_nce, max_targets = input$ms2_target_cap)
      write.csv(lst, file, row.names = FALSE)
    }
  )
  output$dl_frag_ref <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_MS2_fragment_reference.csv"),
    content = function(file) {
      req(rv$ready, rv$mets, rv$dict)
      z2 <- min(input$ms2_z_min, input$ms2_z_max):max(input$ms2_z_min, input$ms2_z_max)
      ref <- ms2_fragment_reference(rv$mets, rv$dict, z_range = z2)
      write.csv(ref, file, row.names = FALSE)
    }
  )

  ## ---- Spectral library downloads ------------------------------------------
  # Built on demand rather than during the run: the MS2 library in
  # particular is large (one spectrum per metabolite per precursor charge,
  # each holding every fragment ion), and most runs never download it.
  .ms1_library <- function() {
    build_ms1_library(rv$mets, rv$dict, z_range = input$z_min:input$z_max,
                      n_iso = input$n_iso, max_oxid = input$max_oxid,
                      h_offset = input$h_offset,
                      use_envipat = input$use_envipat,
                      oligo_name = input$oligo_name)
  }
  .ms2_library <- function() {
    z2 <- min(input$ms2_z_min, input$ms2_z_max):max(input$ms2_z_min, input$ms2_z_max)
    build_ms2_library(rv$mets, rv$dict, precursor_z_range = z2,
                      frag_z_range = 1:input$frag_z_max,
                      h_offset = input$h_offset,
                      oligo_name = input$oligo_name)
  }
  .spectral_download <- function(build, writer, label) {
    downloadHandler(
      filename = function() paste0(input$output_prefix, label),
      content = function(file) {
        req(rv$ready, rv$mets, rv$dict)
        withProgress(message = paste0("Building ", label, " ..."), value = 0.5, {
          writer(build(), file)
        })
      }
    )
  }
  output$dl_ms1_mgf <- .spectral_download(.ms1_library, write_mgf,
                                          "_MS1_library.mgf")
  output$dl_ms1_msp <- .spectral_download(.ms1_library, write_msp,
                                          "_MS1_library.msp")
  output$dl_ms2_mgf <- .spectral_download(.ms2_library, write_mgf,
                                          "_MS2_library.mgf")
  output$dl_ms2_msp <- .spectral_download(.ms2_library, write_msp,
                                          "_MS2_library.msp")

  ## ---- MS2 library explorer -------------------------------------------------
  # Builds the MS2 library once (via the button below) and caches it in
  # ms2_lib_rv, so filtering the table afterwards is instant instead of
  # rebuilding the library on every click.
  observeEvent(input$build_ms2_explorer, {
    req(rv$ready, rv$mets, rv$dict)
    withProgress(message = "Building MS2 library for exploration...", value = 0.3, {
      ms2_lib_rv(.ms2_library())
    })
  })

  output$ms2_explorer_ready <- reactive({
    if (!is.null(ms2_lib_rv()) && length(ms2_lib_rv()) > 0) "true" else "false"
  })
  outputOptions(output, "ms2_explorer_ready", suspendWhenHidden = FALSE)

  # One row per precursor spectrum (metabolite x charge); the row order here
  # matches ms2_lib_rv() exactly, so a DT row-selection index can be used to
  # index straight into the library list below.
  ms2_explorer_summary <- reactive({
    lib <- ms2_lib_rv()
    req(lib, length(lib) > 0)
    data.frame(
      Metabolite = vapply(lib, function(r)
        sub(" \\[z=.*\\] MS2$", "", r$name), ""),
      Charge = factor(vapply(lib, function(r) r$z, numeric(1))),
      `Precursor m/z` = round(vapply(lib, function(r) r$precursor_mz, numeric(1)), 4),
      Formula = vapply(lib, function(r) r$formula, ""),
      `Mono mass (Da)` = round(vapply(lib, function(r) r$mono_mass, numeric(1)), 4),
      `N fragments` = vapply(lib, function(r) nrow(r$peaks), numeric(1)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })

  output$ms2_explorer_table <- DT::renderDataTable({
    DT::datatable(ms2_explorer_summary(), filter = "top", rownames = FALSE,
                  selection = "single",
                  options = list(pageLength = 10,
                                 order = list(list(2, "asc"))))
  })

  output$ms2_explorer_title <- renderUI({
    sel <- input$ms2_explorer_table_rows_selected
    if (is.null(sel))
      return(tags$p(class = "hint",
        "Select a row above to plot its predicted MS2 spectrum."))
    tags$h6(ms2_lib_rv()[[sel]]$name)
  })

  # Basic stick plot: every fragment ion at its m/z, height set by the
  # rule-based intensity heuristic (fragment_intensity_weight() in
  # R/fragments.R) -- a coarse relative ranking, not a predicted abundance
  # (see that function's header for exactly what it does and does not
  # encode). Exact per-peak annotations are in the paired table, not on the
  # plot itself, since a precursor can carry hundreds to thousands of
  # internal-fragment peaks that would make on-plot labels unreadable.
  output$ms2_explorer_plot <- renderPlot({
    sel <- input$ms2_explorer_table_rows_selected
    req(sel)
    r <- ms2_lib_rv()[[sel]]
    ggplot2::ggplot(r$peaks, ggplot2::aes(x = mz, y = intensity)) +
      ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0),
                            color = "#2c3e50", linewidth = 0.3, alpha = 0.7) +
      ggplot2::labs(x = "m/z", y = "Relative intensity (rule-based heuristic)",
                    title = r$name,
                    subtitle = paste0(nrow(r$peaks), " predicted fragment ions")) +
      ggplot2::ylim(0, 110) +
      ggplot2::theme_bw()
  })

  output$ms2_explorer_peaks <- DT::renderDataTable({
    sel <- input$ms2_explorer_table_rows_selected
    req(sel)
    pk <- ms2_lib_rv()[[sel]]$peaks
    pk <- pk[order(pk$mz), ]
    DT::datatable(
      data.frame(`m/z` = round(pk$mz, 4), Annotation = pk$annotation,
                check.names = FALSE, stringsAsFactors = FALSE),
      rownames = FALSE, options = list(pageLength = 15, dom = "tip"))
  })

  ## ---- MS2 mirror plot (single-file mode) ------------------------------------
  # Acquired-vs-theoretical mirror plot for metabolites MS2-confirmed during
  # single-file MS matching (rv$ms_results, built by annotate_metabolites()
  # in R/ms_matching.R -- distinct from the pure-theoretical MS2 Explorer
  # above, which never sees acquired data).
  output$ms2_mirror_ready <- reactive({
    if (!is.null(rv$ms_results) && length(rv$ms_results$ms2_results) > 0) "true" else "false"
  })
  outputOptions(output, "ms2_mirror_ready", suspendWhenHidden = FALSE)

  .ms2_mirror_hits <- reactive({
    req(rv$ms_results)
    s <- rv$ms_results$summary
    req(nrow(s) > 0)
    s[s$has_ms2, , drop = FALSE]
  })

  # .find_met() also exists in R/mirror_plot.R, but as an internal
  # (non-exported) package function -- app.R runs OUTSIDE the package
  # namespace once installed normally (OligoMetProfiler::run_app()), so it
  # can only ever see EXPORTED names there. A local copy is the same fix
  # pattern as every other in-app .xxx helper already defined directly in
  # this file, rather than the first `:::` reference in app.R.
  .find_met <- function(mets, met_id) {
    idx <- which(vapply(mets, function(m) identical(m$id, met_id), logical(1)))
    if (length(idx) == 0) return(NULL)
    mets[[idx[1]]]
  }

  output$ms2_mirror_table <- DT::renderDataTable({
    hits <- .ms2_mirror_hits()
    req(nrow(hits) > 0)
    df <- data.frame(
      Metabolite = hits$met_name, Kind = hits$kind, n = hits$n,
      `MS2 score` = hits$ms2_score, `Coverage %` = round(100 * hits$ms2_coverage, 1),
      `N frag. matches` = hits$n_ms2_frags, Confident = hits$confident,
      check.names = FALSE, stringsAsFactors = FALSE)
    DT::datatable(df, rownames = FALSE, selection = "single",
                  options = list(pageLength = 10, dom = "tip"))
  })

  output$ms2_mirror_plot <- renderPlot({
    hits <- .ms2_mirror_hits()
    sel <- input$ms2_mirror_table_rows_selected
    req(sel)
    met_id <- hits$met_id[sel]
    r <- rv$ms_results$ms2_results[[met_id]]
    req(r, r$best_spec)
    met <- .find_met(rv$mets, met_id)
    req(met)
    spec <- mirror_spectrum_data(met, r$best_spec, rv$dict,
                                  tol_ppm = input$frag_tol_ppm,
                                  z_range = 1:input$frag_z_max)
    plot_mirror_spectrum(spec, title = paste0(met$name, " -- MS2 mirror plot"))
  })

  ## ---- Batch MS Processing / Unidentified Peaks / Statistics tabs -----------
  output$batch_ready <- reactive({
    if (!is.null(rv$batch_ms_results) && nrow(rv$batch_ms_results$ms1_matches) >= 0) "true" else "false"
  })
  outputOptions(output, "batch_ready", suspendWhenHidden = FALSE)

  output$stats_ready <- reactive({
    if (!is.null(rv$stats_results) && !is.null(rv$stats_results$result)) "true" else "false"
  })
  outputOptions(output, "stats_ready", suspendWhenHidden = FALSE)

  output$kind_stats_ready <- reactive({
    if (!is.null(rv$kind_stats_results) && !is.null(rv$kind_stats_results$result)) "true" else "false"
  })
  outputOptions(output, "kind_stats_ready", suspendWhenHidden = FALSE)

  ## ---- Quantification tab ----------------------------------------------------
  output$quant_ready <- reactive({
    q <- rv$quant_results
    # Also true when a calibration curve was ATTEMPTED but failed to fit --
    # calibration_curves_table exists specifically to surface WHY (too few
    # standard points, no standard-type rows, ...) via its own "note"
    # column, so gating the whole panel on absolute/relative having actual
    # rows used to hide that explanation right when it was most needed.
    !is.null(q) && (nrow(q$absolute) > 0 || nrow(q$relative) > 0 || length(q$calibration_curves) > 0)
  })
  outputOptions(output, "quant_ready", suspendWhenHidden = FALSE)

  # A calibration curve needs standard-type rows/concentrations in the
  # metadata table AND the metabolite selected under "Calibration &
  # Quantification (Advanced)" in the sidebar (absolute_quant_mets) -- two
  # separate steps, easy to do one and forget the other. The generic
  # "run batch processing..." message doesn't distinguish those cases, so
  # a user who filled in Sample Type = standard but never picked a
  # metabolite there sees the exact same text as someone who hasn't
  # touched the metadata table at all -- this fills in which one it is.
  output$quant_not_ready_note <- renderUI({
    meta <- rv$sample_meta
    n_std <- if (!is.null(meta) && "sample_type" %in% names(meta)) {
      sum(.map_sample_type(meta$sample_type) == "standard", na.rm = TRUE)
    } else 0L
    msg <- if (n_std > 0 && length(input$absolute_quant_mets %||% character(0)) == 0) {
      sprintf(paste("%d sample(s) are marked Sample Type = standard, but no metabolite is",
                     "selected to build a calibration curve from -- pick one under",
                     "\"Calibration & Quantification (Advanced)\" in the sidebar,",
                     "then click Run Statistical Analysis."), n_std)
    } else {
      "Run batch processing with Group/Timepoint (and Sample Type = standard for calibration) filled in to see this here."
    }
    tags$p(style = "padding-top: 12px; color: #6c757d;", msg)
  })

  # Flattens rv$quant_results$calibration_curves (a named list of
  # fit_calibration_curve() results, one per absolute-quant metabolite) into
  # one row per metabolite -- surfaces the "note" field that explains WHY a
  # curve is missing/unreliable (too few standard points, no standard rows
  # at all, ...) instead of that metabolite just silently having no rows in
  # quant_absolute_table below.
  .calibration_curves_display <- reactive({
    req(rv$quant_results)
    cc <- rv$quant_results$calibration_curves
    req(length(cc) > 0)
    do.call(rbind, lapply(cc, function(c) {
      mean_abs_re <- if (!is.null(c$points) && nrow(c$points) > 0 && "percent_re" %in% names(c$points))
        mean(abs(c$points$percent_re), na.rm = TRUE) else NA_real_
      data.frame(met_id = c$met_id, weighting = c$weighting,
                 r_squared = c$r_squared, n_points = c$n_points,
                 mean_abs_percent_re = mean_abs_re,
                 note = if (nzchar(c$note)) c$note else "ok",
                 stringsAsFactors = FALSE)
    }))
  })
  output$calibration_curves_table <- DT::renderDT({
    DT::datatable(.calibration_curves_display(), rownames = FALSE,
                  options = list(pageLength = 15, scrollX = TRUE)) |>
      DT::formatRound(c("r_squared", "mean_abs_percent_re"), digits = 3)
  })
  output$dl_calibration_curves_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_calibration_curves.csv"),
    content = function(file) utils::write.csv(.calibration_curves_display(), file, row.names = FALSE)
  )

  output$quant_absolute_table <- DT::renderDT({
    req(rv$quant_results)
    DT::datatable(rv$quant_results$absolute, rownames = FALSE,
                   options = list(pageLength = 15, scrollX = TRUE)) |>
      DT::formatRound(intersect(c("signal", "concentration_calc", "nominal_concentration",
                                   "percent_re", "curve_r_squared"),
                                 names(rv$quant_results$absolute)), digits = 3)
  })
  output$dl_quant_absolute_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_absolute_quantification.csv"),
    content = function(file) utils::write.csv(rv$quant_results$absolute, file, row.names = FALSE)
  )

  output$quant_relative_table <- DT::renderDT({
    req(rv$quant_results)
    rel <- rv$quant_results$relative
    disp_cols <- setdiff(names(rel), c(".time", ".arm"))
    DT::datatable(rel[, disp_cols, drop = FALSE], rownames = FALSE,
                   options = list(pageLength = 15, scrollX = TRUE)) |>
      DT::formatRound(intersect(c("signal", "baseline_signal", "relative_signal"), disp_cols), digits = 3)
  })
  output$dl_quant_relative_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_relative_quantification.csv"),
    content = function(file) utils::write.csv(rv$quant_results$relative, file, row.names = FALSE)
  )

  ## ---- Ask OligoMet: in-app LLM chat assistant ------------------------------
  # See R/agent_core.R (run_agent_turn(), the provider-agnostic ReAct loop)
  # and R/agent_tools.R (the tool registry it calls). This tab is the
  # "naive researcher" front end from /root/.claude/plans/staged-mapping-
  # blanket.md -- Track A; the MCP server (Track B) shares the same tool
  # registry from a separate process, see inst/mcp_server/.
  output$agent_key_status <- renderText({
    a <- nzchar(Sys.getenv("ANTHROPIC_API_KEY", ""))
    o <- nzchar(Sys.getenv("OPENAI_API_KEY", ""))
    paste0("Anthropic key: ", if (a) "configured" else "not set",
          "  |  OpenAI key: ", if (o) "configured" else "not set")
  })

  # Turns the canonical message list (content BLOCKS -- text/tool_use/
  # tool_result) into flat (role, text) entries for display. Tool calls get
  # their own compact "reasoning trace" line (matching MSAgent's own
  # transparency design) rather than being hidden; tool_result blocks
  # (the raw JSON a tool returned) are not shown directly -- the
  # assistant's own next text turn is what narrates them.
  .agent_display_entries <- function(messages) {
    entries <- list()
    for (m in messages) {
      for (b in m$content) {
        if (identical(b$type, "text") && nzchar(b$text %||% "")) {
          entries[[length(entries) + 1]] <- list(role = m$role, text = b$text)
        } else if (identical(b$type, "tool_use")) {
          entries[[length(entries) + 1]] <- list(
            role = "tool",
            text = paste0("Called ", b$name, "(",
                         as.character(jsonlite::toJSON(b$input, auto_unbox = TRUE, na = "null")), ")"))
        }
      }
    }
    entries
  }

  .agent_chat_bubble_html <- function(role, text) {
    esc <- gsub("\n", "<br/>", htmltools::htmlEscape(text), fixed = TRUE)
    style <- switch(role,
      user = "background:#0279EE; color:#fff; margin-left:auto; margin-right:0;",
      tool = paste("background:#fff; color:#6c757d; border:1px dashed #ccc;",
                   "font-family:monospace; font-size:11px; margin-right:auto;"),
      "background:#fff; border:1px solid #ddd; margin-right:auto;")
    label <- switch(role, user = "You", tool = "Tool call", "Assistant")
    sprintf(paste0('<div style="max-width:80%%; padding:8px 12px; margin:6px 0;',
                  'border-radius:8px; %s">',
                  '<div style="font-size:10px; opacity:0.7; margin-bottom:2px;">%s</div>',
                  '%s</div>'),
           style, label, esc)
  }

  output$agent_chat_html <- renderUI({
    entries <- .agent_display_entries(rv$agent_messages)
    if (length(entries) == 0) {
      return(tags$p(style = "color: #6c757d;",
                    "Ask a question about your sequence, a match, or what's degrading."))
    }
    HTML(paste(vapply(entries, function(e) .agent_chat_bubble_html(e$role, e$text), character(1)),
              collapse = ""))
  })

  output$agent_busy_flag <- reactive(if (isTRUE(rv$agent_busy)) "true" else "false")
  outputOptions(output, "agent_busy_flag", suspendWhenHidden = FALSE)

  observeEvent(input$agent_clear, {
    rv$agent_messages <- list()
    rv$agent_ctx <- list()
  })

  observeEvent(input$agent_send, {
    txt <- trimws(input$agent_input %||% "")
    req(nzchar(txt))
    updateTextAreaInput(session, "agent_input", value = "")
    rv$agent_busy <- TRUE

    # The agent's own accumulated context (whatever it has already parsed/
    # built within this chat) takes priority; the main app's currently
    # loaded sequence/library/batch results are only the DEFAULT for
    # anything the chat hasn't touched yet.
    base_ctx <- list(dict = rv$dict %||% STANDARD_DICT)
    if (!is.null(rv$spec)) base_ctx$spec <- rv$spec
    if (!is.null(rv$mets)) base_ctx$mets <- rv$mets
    if (!is.null(rv$batch_ms_results) && !is.null(rv$batch_ms_results$ms1_matches) &&
        nrow(rv$batch_ms_results$ms1_matches) > 0) {
      base_ctx$ms1_matches <- rv$batch_ms_results$ms1_matches
    }
    turn_ctx <- utils::modifyList(base_ctx, rv$agent_ctx %||% list())
    messages <- c(rv$agent_messages, list(user_turn(txt)))

    result <- tryCatch(
      run_agent_turn(messages, ctx = turn_ctx, backend = input$agent_backend %||% "anthropic"),
      error = function(e) list(
        messages = c(messages, list(list(role = "assistant",
                                         content = list(list(type = "text",
                                                             text = paste("Error:", conditionMessage(e))))))),
        ctx = turn_ctx)
    )
    rv$agent_messages <- result$messages
    rv$agent_ctx <- result$ctx
    rv$agent_busy <- FALSE
  })

  # Merges in the MS2 confirmation columns (n_ms2_peaks, coverage,
  # confirmation_score, confident, ...) when batch MS2 confirmation ran --
  # shared between the table render and the mirror-plot row lookup below so
  # a DT row-selection index always resolves against the same row order.
  .batch_matches_display <- reactive({
    req(rv$batch_ms_results)
    m <- rv$batch_ms_results$ms1_matches
    ms2c <- rv$batch_ms_results$ms2_confirmations
    if (!is.null(ms2c) && nrow(ms2c) > 0) {
      ms2c$met_name <- NULL
      m <- merge(m, ms2c, by = c("sample", "met_id", "k_oxid", "z", "adduct"),
                 all.x = TRUE, sort = FALSE)
    }
    m
  })

  output$batch_matches_table <- DT::renderDT({
    m <- .batch_matches_display()
    DT::datatable(m, filter = "top", rownames = FALSE, selection = "single",
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$batch_mirror_plot <- renderPlot({
    m <- .batch_matches_display()
    sel <- input$batch_matches_table_rows_selected
    req(sel)
    row <- m[sel, ]
    req(!is.null(row$n_ms2_peaks), !is.na(row$n_ms2_peaks), row$n_ms2_peaks > 0)
    key <- paste(row$sample, row$met_id, row$k_oxid, row$z, row$adduct, sep = "|")
    best_spec <- rv$batch_ms_results$ms2_spectra[[key]]
    req(best_spec)
    met <- .find_met(rv$mets, row$met_id)
    req(met)
    spec <- mirror_spectrum_data(met, best_spec, rv$dict,
                                  tol_ppm = input$frag_tol_ppm,
                                  z_range = 1:input$frag_z_max)
    plot_mirror_spectrum(spec, title = paste0(row$met_name, " (", row$sample,
                                               ") -- MS2 mirror plot"))
  })

  output$dl_batch_mirror_pdf <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_batch_MS2_mirror_plots.pdf"),
    content = function(file) {
      req(rv$batch_ms_results, rv$mets, rv$dict)
      batch_mirror_plots_pdf(rv$mets, rv$dict,
                              rv$batch_ms_results$ms2_confirmations,
                              rv$batch_ms_results$ms2_spectra,
                              file = file, tol_ppm = input$frag_tol_ppm,
                              z_range = 1:input$frag_z_max, h_offset = input$h_offset)
    }
  )

  output$dl_batch_annotated_msp <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_batch_MS2_annotated.msp"),
    content = function(file) {
      req(rv$batch_ms_results, rv$mets, rv$dict)
      recs <- batch_annotated_msp_records(rv$mets, rv$dict,
                                           rv$batch_ms_results$ms2_confirmations,
                                           rv$batch_ms_results$ms2_spectra,
                                           tol_ppm = input$frag_tol_ppm,
                                           z_range = 1:input$frag_z_max,
                                           h_offset = input$h_offset)
      write_msp(recs, file, measured = TRUE)
    }
  )

  # One consensus spectrum per (metabolite, oxidation, charge, adduct),
  # pooled across every sample/replicate that confirmed it -- see
  # build_empirical_ms2_library() in R/export_spectral.R for why this is
  # not gated on confirmation_score/coverage. Rebuilt on demand, same
  # reasoning as .ms1_library()/.ms2_library() above.
  .empirical_ms2_library <- function() {
    req(rv$batch_ms_results, rv$mets, rv$dict)
    build_empirical_ms2_library(rv$mets, rv$dict,
                                rv$batch_ms_results$ms2_confirmations,
                                rv$batch_ms_results$ms2_spectra,
                                frag_z_range = 1:input$frag_z_max,
                                h_offset = input$h_offset,
                                fragment_ppm = input$frag_tol_ppm)
  }
  output$dl_empirical_ms2_msp <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_MS2_library_EMPIRICAL.msp"),
    content = function(file) write_msp(.empirical_ms2_library()$records, file, measured = TRUE)
  )
  output$dl_empirical_ms2_summary <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_empirical_MS2_library_summary.csv"),
    content = function(file) write.csv(.empirical_ms2_library()$summary, file, row.names = FALSE)
  )
  output$empirical_ms2_summary_table <- DT::renderDT({
    tryCatch({
      s <- .empirical_ms2_library()$summary
      req(nrow(s) > 0)
      DT::datatable(s, filter = "top", rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE))
    }, error = function(e) DT::datatable(data.frame(), rownames = FALSE))
  })

  output$unmatched_table <- DT::renderDT({
    req(rv$batch_ms_results)
    DT::datatable(rv$batch_ms_results$unmatched, filter = "top", rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  ## ---- Degradation tab (M4: R/degradation.R) --------------------------------
  output$degradation_per_sample_table <- DT::renderDT({
    req(rv$batch_ms_results$degradation)
    DT::datatable(rv$batch_ms_results$degradation$per_sample, rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$plot_degradation_composition <- renderPlot({
    req(rv$batch_ms_results$degradation)
    plot_degradation_composition(rv$batch_ms_results$degradation)
  })

  output$degradation_composition_table <- DT::renderDT({
    req(rv$batch_ms_results$degradation)
    DT::datatable(rv$batch_ms_results$degradation$composition, rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$degradation_top_table <- DT::renderDT({
    req(rv$batch_ms_results$degradation)
    DT::datatable(rv$batch_ms_results$degradation$top_degradants, rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$dl_degradation_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_degradation_summary.csv"),
    content = function(file) {
      req(rv$batch_ms_results$degradation)
      utils::write.csv(rv$batch_ms_results$degradation$per_sample, file, row.names = FALSE)
    }
  )

  # Flattens compare_multi_groups()'s list(omnibus, posthoc) to a single
  # table for display; two_group/time_series results are already flat.
  .stats_display_table <- reactive({
    sr <- rv$stats_results
    req(sr)
    if (sr$mode == "multi_group") sr$result$omnibus else sr$result
  })

  output$stats_met_selector <- renderUI({
    df <- .stats_display_table()
    req(nrow(df) > 0)
    selectInput("stats_met_select", "Metabolite", choices = stats::setNames(df$met_id, df$met_name))
  })

  # .is_study_sample() also exists in R/chemistry_dict.R, but as an
  # internal (non-exported) package function -- app.R runs OUTSIDE the
  # package namespace once installed normally (OligoMetProfiler::run_app()),
  # so it can only ever see EXPORTED names there (confirmed by a real user
  # report: "could not find function \".is_study_sample\""). A local copy
  # is the same fix pattern as .find_met() above and every other in-app
  # .xxx helper already defined directly in this file.
  .is_study_sample <- function(sample_type) {
    is.na(sample_type) | !nzchar(sample_type) | sample_type == "unknown"
  }

  # Calibration standards/QC/blanks are for assay performance, not the
  # study design -- excluded here (same rule as quantify_relative()/
  # degradation_summary()) so they don't show up as extra points with no
  # real x-axis value (blank Group/Timepoint) on the trend/boxplot below.
  .study_sample_meta <- reactive({
    meta <- rv$sample_meta
    req(meta)
    if ("sample_type" %in% names(meta)) meta[.is_study_sample(meta$sample_type), ] else meta
  })

  output$plot_stats_main <- renderPlot({
    sr <- rv$stats_results
    req(sr, input$stats_met_select)
    if (sr$mode == "two_group") {
      plot_volcano(sr$result)
    } else if (sr$mode == "time_series") {
      abund <- build_abundance_matrix(rv$batch_ms_results$ms1_matches)
      long <- abundance_long(abund, .study_sample_meta()[, c("sample", "timepoint")])
      plot_trend(long, input$stats_met_select)
    } else {
      abund <- build_abundance_matrix(rv$batch_ms_results$ms1_matches)
      long <- abundance_long(abund, .study_sample_meta()[, c("sample", "group")])
      plot_group_boxplot(long, input$stats_met_select)
    }
  })

  output$stats_table <- DT::renderDT({
    df <- .stats_display_table()
    dt <- DT::datatable(df, rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
    # Bold+color any row whose |log2fc| clears the sidebar threshold --
    # two_group mode only, since multi_group/time_series results don't have
    # a log2fc column at all.
    thr <- input$stats_log2fc_threshold %||% 1
    if ("log2fc" %in% names(df) && is.finite(thr) && thr > 0) {
      dt <- DT::formatStyle(dt, "log2fc",
        fontWeight = DT::styleInterval(c(-thr, thr), c("bold", "normal", "bold")),
        color = DT::styleInterval(c(-thr, thr), c("#a3231b", "inherit", "#18632f")))
    }
    dt
  })

  # Same flatten as .stats_display_table(), for the kind-level (composition-
  # class) comparison -- its met_id/met_name columns hold "exo_3p"/
  # "endo_5frag"/etc rather than real metabolite IDs (see
  # build_kind_abundance_matrix() in R/statistics.R for why that's fine).
  .kind_stats_display_table <- reactive({
    sr <- rv$kind_stats_results
    req(sr)
    if (sr$mode == "multi_group") sr$result$omnibus else sr$result
  })

  output$kind_stats_table <- DT::renderDT({
    df <- .kind_stats_display_table()
    # met_id and met_name are identical here (both hold the kind string --
    # see build_kind_abundance_matrix()) -- show one column, labeled "kind".
    df$met_name <- NULL
    names(df)[names(df) == "met_id"] <- "kind"
    DT::datatable(df, rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$dl_kind_stats_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_class_comparison.csv"),
    content = function(file) {
      utils::write.csv(.kind_stats_display_table(), file, row.names = FALSE)
    }
  )

  ## ---- Time Course tab: multi-metabolite trend --------------------------------
  # Metabolite choices for both this and the Multivariate tab below come
  # straight from the batch matches, not from any stats result -- unlike
  # stats_met_selector (which only lists metabolites a compare_*() call
  # already succeeded on), these should offer every matched metabolite even
  # before/without a two-group or 3+-group comparison being meaningful.
  observeEvent(rv$batch_ms_results, {
    bres <- rv$batch_ms_results
    if (is.null(bres) || nrow(bres$ms1_matches) == 0) {
      updateSelectizeInput(session, "tc_met_ids", choices = character(0))
      updateSelectizeInput(session, "mv_met_ids", choices = character(0))
      return()
    }
    met_info <- unique(bres$ms1_matches[, c("met_id", "met_name")])
    choices <- stats::setNames(met_info$met_id, paste0(met_info$met_name, " (", met_info$met_id, ")"))

    # Time Course default: top few metabolites by total signal, so the
    # trend plot isn't empty (nothing selected) or unreadable (every
    # metabolite at once) the first time this tab is opened.
    sig_col <- if ("area" %in% names(bres$ms1_matches) && any(!is.na(bres$ms1_matches$area))) "area" else "intensity"
    tot <- stats::aggregate(stats::as.formula(paste(sig_col, "~ met_id")),
                             data = bres$ms1_matches, FUN = sum, na.rm = TRUE)
    top_default <- tot$met_id[order(-tot[[sig_col]])][seq_len(min(6, nrow(tot)))]
    keep_tc <- intersect(isolate(input$tc_met_ids), met_info$met_id)
    updateSelectizeInput(session, "tc_met_ids", choices = choices,
                         selected = if (length(keep_tc) > 0) keep_tc else top_default)

    # Multivariate default: every metabolite (NULL met_ids in run_pca()/
    # run_hclust() below) -- narrowing is opt-in, unlike the trend view.
    keep_mv <- intersect(isolate(input$mv_met_ids), met_info$met_id)
    updateSelectizeInput(session, "mv_met_ids", choices = choices, selected = keep_mv)
  }, ignoreNULL = FALSE)

  .time_course_long <- reactive({
    bres <- rv$batch_ms_results
    if (is.null(bres) || nrow(bres$ms1_matches) == 0) return(data.frame())
    meta <- .study_sample_meta()
    if (!"timepoint" %in% names(meta)) return(data.frame())
    meta$timepoint <- suppressWarnings(as.numeric(meta$timepoint))
    meta <- meta[!is.na(meta$timepoint), c("sample", "timepoint")]
    if (nrow(meta) == 0) return(data.frame())
    abund <- build_abundance_matrix(bres$ms1_matches)
    abundance_long(abund, meta)
  })

  .time_course_summary <- reactive({
    long <- .time_course_long()
    if (nrow(long) == 0) return(data.frame())
    met_ids <- if (length(input$tc_met_ids) > 0) input$tc_met_ids else NULL
    multi_trend_summary(long, met_ids = met_ids, normalize = isTRUE(input$tc_normalize))
  })

  output$time_course_note <- renderUI({
    long <- .time_course_long()
    msg <- if (nrow(long) == 0) {
      "No Timepoint design detected (Batch Processing tab) after excluding calibration standards/QC/blanks."
    } else if (length(unique(long$timepoint)) < 2) {
      "Only one distinct timepoint found -- need at least two to show a trend."
    } else NA_character_
    req(!is.na(msg))
    tags$p(style = "font-size: 11px; color: #a3231b;", msg)
  })

  output$plot_time_course_trend <- renderPlot({
    plot_multi_trend(.time_course_summary(), normalize = isTRUE(input$tc_normalize))
  })

  output$time_course_table <- DT::renderDT({
    df <- .time_course_summary()
    req(nrow(df) > 0)
    DT::datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE)) |>
      DT::formatRound(c("mean_value", "sd", "sem"), digits = 3)
  })

  output$dl_time_course_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_time_course_trend.csv"),
    content = function(file) utils::write.csv(.time_course_summary(), file, row.names = FALSE)
  )

  ## ---- Multivariate tab: PCA + hierarchical clustering -------------------------
  # Same (log-transformed, study-samples-only) matrix underlies both --
  # run_pca()/run_hclust() (R/multivariate.R) each build it independently
  # via the shared .build_multivariate_matrix() helper, so the PCA plot and
  # the dendrogram below are always looking at the exact same filtered data.
  .pca_result <- reactive({
    bres <- rv$batch_ms_results
    if (is.null(bres) || nrow(bres$ms1_matches) == 0) {
      return(list(scores = data.frame(), loadings = data.frame(), var_explained = numeric(0),
                  dropped_zero_variance = character(0), note = "no batch results yet"))
    }
    met_ids <- if (length(input$mv_met_ids) > 0) input$mv_met_ids else NULL
    run_pca(bres$ms1_matches, sample_meta = rv$sample_meta, met_ids = met_ids,
            log_transform = isTRUE(input$mv_log), scale = isTRUE(input$mv_scale))
  })

  .hclust_result <- reactive({
    bres <- rv$batch_ms_results
    if (is.null(bres) || nrow(bres$ms1_matches) == 0) {
      return(list(hclust = NULL, clusters = data.frame(), k = 0L,
                  dropped_zero_variance = character(0), note = "no batch results yet"))
    }
    met_ids <- if (length(input$mv_met_ids) > 0) input$mv_met_ids else NULL
    k <- suppressWarnings(as.integer(input$mv_k))
    if (is.na(k) || k < 2) k <- 2
    run_hclust(bres$ms1_matches, sample_meta = rv$sample_meta, met_ids = met_ids,
               log_transform = isTRUE(input$mv_log), scale = isTRUE(input$mv_scale), k = k)
  })

  output$multivariate_note <- renderUI({
    note <- .pca_result()$note
    req(nzchar(note %||% ""))
    tags$p(style = "font-size: 11px; color: #a3231b;", note)
  })

  output$plot_pca_scores <- renderPlot({
    plot_pca_scores(.pca_result())
  })

  output$pca_loadings_table <- DT::renderDT({
    ld <- .pca_result()$loadings
    req(nrow(ld) > 0)
    ld <- ld[order(-abs(ld$PC1)), ]
    DT::datatable(ld, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE)) |>
      DT::formatRound(setdiff(names(ld), c("met_id", "met_name")), digits = 3)
  })

  output$dl_pca_scores_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_pca_scores.csv"),
    content = function(file) utils::write.csv(.pca_result()$scores, file, row.names = FALSE)
  )
  output$dl_pca_loadings_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_pca_loadings.csv"),
    content = function(file) utils::write.csv(.pca_result()$loadings, file, row.names = FALSE)
  )

  output$plot_dendrogram <- renderPlot({
    hc <- .hclust_result()
    plot_dendrogram(hc, k = hc$k)
  })

  output$hclust_clusters_table <- DT::renderDT({
    cl <- .hclust_result()$clusters
    req(nrow(cl) > 0)
    DT::datatable(cl, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

  output$dl_hclust_clusters_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_cluster_assignments.csv"),
    content = function(file) utils::write.csv(.hclust_result()$clusters, file, row.names = FALSE)
  )

  output$dl_batch_tsv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_batch_features.tsv"),
    content = function(file) {
      req(rv$batch_features)
      utils::write.table(rv$batch_features, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  ## ---- Statistical Analysis: data matrix downloads --------------------------
  # build_abundance_matrix()/abundance_long() already collapse to the
  # max-intensity match per (metabolite, sample) across charge states --
  # that's the SAME collapsing the batch pipeline's charge_group.py step
  # already does before matching, so this is one real matrix, not two
  # genuinely distinct raw-vs-grouped ones (see the Data Matrix tab's own
  # note; the Charge Grouping panel's "Re-run Aggregation" is a placeholder
  # for building a real second, post-hoc-grouped matrix).
  output$dl_data_matrix_wide <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_data_matrix_wide.csv"),
    content = function(file) {
      req(rv$batch_ms_results)
      utils::write.csv(build_abundance_matrix(rv$batch_ms_results$ms1_matches), file, row.names = FALSE)
    }
  )
  output$dl_data_matrix_long <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_data_matrix_long.csv"),
    content = function(file) {
      req(rv$batch_ms_results, rv$sample_meta)
      abund <- build_abundance_matrix(rv$batch_ms_results$ms1_matches)
      utils::write.csv(abundance_long(abund, rv$sample_meta), file, row.names = FALSE)
    }
  )
  output$dl_unmatched_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_unidentified_peaks.csv"),
    content = function(file) {
      req(rv$batch_ms_results)
      utils::write.csv(rv$batch_ms_results$unmatched, file, row.names = FALSE)
    }
  )
  output$dl_stats_csv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_statistics.csv"),
    content = function(file) {
      req(rv$stats_results)
      utils::write.csv(.stats_display_table(), file, row.names = FALSE)
    }
  )

  ## ---- Download all outputs as a single .zip -------------------------------
  # Re-derives the CSV/spectral-library files rather than reusing the
  # download handlers above (downloadHandlers aren't callable as plain
  # functions), writing everything into one scratch directory before zipping.
  # Uses the "zip" package when available (pure R, no external binary) and
  # falls back to utils::zip() (shells out to a system "zip") otherwise.
  output$dl_all <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_all_outputs.zip"),
    content = function(file) {
      req(rv$ready, rv$mets, rv$dict, rv$prm)
      withProgress(message = "Bundling all outputs...", value = 0, {
        prefix <- input$output_prefix
        bundle_dir <- tempfile("oligomet_bundle_")
        dir.create(bundle_dir)
        on.exit(unlink(bundle_dir, recursive = TRUE), add = TRUE)

        incProgress(0.05, detail = "Workbook")
        if (!is.null(rv$wb_path) && file.exists(rv$wb_path))
          file.copy(rv$wb_path, file.path(bundle_dir, paste0(prefix, "_library.xlsx")))

        incProgress(0.05, detail = "Report")
        if (!is.null(rv$report_path) && file.exists(rv$report_path))
          file.copy(rv$report_path, file.path(bundle_dir, paste0(prefix, "_report.html")))

        incProgress(0.05, detail = "PRM list")
        utils::write.csv(rv$prm, file.path(bundle_dir, paste0(prefix, "_prm_list.csv")),
                         row.names = FALSE)

        z_range <- input$z_min:input$z_max
        z2 <- min(input$ms2_z_min, input$ms2_z_max):max(input$ms2_z_min, input$ms2_z_max)

        incProgress(0.1, detail = "MS1 inclusion list")
        utils::write.csv(
          thermo_ms1_inclusion_list(rv$mets, rv$dict, z_range = z_range,
            h_offset = input$h_offset, max_oxid = input$max_oxid,
            rt_start = 0, rt_end = input$method_length,
            max_targets = input$ms1_target_cap),
          file.path(bundle_dir, paste0(prefix, "_MS1_inclusion_list.csv")),
          row.names = FALSE)

        incProgress(0.1, detail = "MS2 PRM target list")
        utils::write.csv(
          thermo_ms2_prm_target_list(rv$mets, rv$dict, z_range = z2,
            h_offset = input$h_offset, max_oxid = input$max_oxid,
            rt_start = 0, rt_end = input$method_length,
            nce = input$hcd_nce, max_targets = input$ms2_target_cap),
          file.path(bundle_dir, paste0(prefix, "_MS2_PRM_target_list.csv")),
          row.names = FALSE)

        incProgress(0.1, detail = "MS2 fragment reference")
        utils::write.csv(ms2_fragment_reference(rv$mets, rv$dict, z_range = z2),
                         file.path(bundle_dir, paste0(prefix, "_MS2_fragment_reference.csv")),
                         row.names = FALSE)

        incProgress(0.2, detail = "MS1 spectral library")
        ms1_lib <- .ms1_library()
        write_mgf(ms1_lib, file.path(bundle_dir, paste0(prefix, "_MS1_library.mgf")))
        write_msp(ms1_lib, file.path(bundle_dir, paste0(prefix, "_MS1_library.msp")))

        incProgress(0.2, detail = "MS2 spectral library")
        ms2_lib <- .ms2_library()
        write_mgf(ms2_lib, file.path(bundle_dir, paste0(prefix, "_MS2_library.mgf")))
        write_msp(ms2_lib, file.path(bundle_dir, paste0(prefix, "_MS2_library.msp")))

        if (!is.null(rv$batch_features)) {
          incProgress(0, detail = "Batch MS features")
          utils::write.table(rv$batch_features,
            file.path(bundle_dir, paste0(prefix, "_batch_features.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)
          utils::write.csv(rv$batch_ms_results$unmatched,
            file.path(bundle_dir, paste0(prefix, "_unidentified_peaks.csv")), row.names = FALSE)
          if (!is.null(rv$stats_results)) {
            utils::write.csv(.stats_display_table(),
              file.path(bundle_dir, paste0(prefix, "_statistics.csv")), row.names = FALSE)
          }
          if (!is.null(rv$kind_stats_results)) {
            utils::write.csv(.kind_stats_display_table(),
              file.path(bundle_dir, paste0(prefix, "_class_comparison.csv")), row.names = FALSE)
          }
          if (!is.null(rv$batch_ms_results$degradation) &&
              nrow(rv$batch_ms_results$degradation$per_sample) > 0) {
            utils::write.csv(rv$batch_ms_results$degradation$per_sample,
              file.path(bundle_dir, paste0(prefix, "_degradation_summary.csv")), row.names = FALSE)
          }
          ms2c <- rv$batch_ms_results$ms2_confirmations
          if (!is.null(ms2c) && nrow(ms2c) > 0) {
            incProgress(0, detail = "Batch MS2 mirror plots")
            batch_mirror_plots_pdf(rv$mets, rv$dict, ms2c, rv$batch_ms_results$ms2_spectra,
              file = file.path(bundle_dir, paste0(prefix, "_batch_MS2_mirror_plots.pdf")),
              tol_ppm = input$frag_tol_ppm, z_range = 1:input$frag_z_max,
              h_offset = input$h_offset)
            recs <- batch_annotated_msp_records(rv$mets, rv$dict, ms2c,
              rv$batch_ms_results$ms2_spectra, tol_ppm = input$frag_tol_ppm,
              z_range = 1:input$frag_z_max, h_offset = input$h_offset)
            write_msp(recs, file.path(bundle_dir, paste0(prefix, "_batch_MS2_annotated.msp")),
                      measured = TRUE)

            incProgress(0, detail = "Empirical MS2 library")
            emp <- build_empirical_ms2_library(rv$mets, rv$dict, ms2c,
              rv$batch_ms_results$ms2_spectra, frag_z_range = 1:input$frag_z_max,
              h_offset = input$h_offset, fragment_ppm = input$frag_tol_ppm)
            utils::write.csv(emp$summary,
              file.path(bundle_dir, paste0(prefix, "_empirical_MS2_library_summary.csv")),
              row.names = FALSE)
            if (length(emp$records) > 0) {
              write_msp(emp$records,
                file.path(bundle_dir, paste0(prefix, "_MS2_library_EMPIRICAL.msp")),
                measured = TRUE)
            }
          }
        }

        incProgress(0.15, detail = "Compressing")
        zipfile <- normalizePath(file, mustWork = FALSE)
        old_wd <- setwd(bundle_dir)
        on.exit(setwd(old_wd), add = TRUE)
        if (requireNamespace("zip", quietly = TRUE)) {
          zip::zip(zipfile, list.files("."))
        } else {
          utils::zip(zipfile, list.files("."))
        }
      })
    }
  )
}

## ---- Run app ---------------------------------------------------------------
shinyApp(ui = ui, server = server)
