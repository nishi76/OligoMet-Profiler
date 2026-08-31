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
               "ms_matching.R", "batch_ms_processing.R", "statistics.R",
               "build_workbook.R", "build_report.R",
               "export_acquisition.R", "export_spectral.R")) {
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
  "enable_ms", "ppm_tol", "adducts", "frag_tol_ppm", "frag_z_max",
  "enable_batch", "batch_run_ms2", "batch_n_workers", "batch_deconv_ppm",
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
  adducts        = function(s, v) updateCheckboxGroupInput(s, "adducts", selected = v),
  frag_tol_ppm   = function(s, v) updateNumericInput(s, "frag_tol_ppm", value = v),
  frag_z_max     = function(s, v) updateNumericInput(s, "frag_z_max", value = v),
  enable_batch   = function(s, v) updateCheckboxInput(s, "enable_batch", value = v),
  batch_run_ms2  = function(s, v) updateCheckboxInput(s, "batch_run_ms2", value = v),
  batch_n_workers  = function(s, v) updateNumericInput(s, "batch_n_workers", value = v),
  batch_deconv_ppm = function(s, v) updateNumericInput(s, "batch_deconv_ppm", value = v),
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
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

  titlePanel("OligoMet Profiler"),

  ## Research-use-only banner. Deliberately above the fold and not
  ## dismissible: every number this app produces is a prediction, and the
  ## outputs get shared as files that leave the app.
  tags$div(class = "ruo-banner",
    tags$span(class = "ruo-tag", "RESEARCH USE ONLY"),
    tags$span("Not for diagnostic, clinical, or regulatory submission use. ",
              "All values are computed predictions, not measurements -- ",
              "confirm every assignment experimentally. Provided without ",
              "warranty; the author accepts no liability for their use. "),
    tags$a(href = "#about-section", "Full disclaimer")
  ),

  sidebarLayout(

    ## ---- Sidebar: all parameters ------------------------------------------
    sidebarPanel(
      width = 4,
      tags$head(tags$style(HTML("
        .sidebar-section { margin-bottom: 18px; }
        .sidebar-section h5 { font-weight: 600; color: #2c3e50;
                              border-bottom: 1px solid #ecf0f1; padding-bottom: 4px; }
        .metric-card { background: #f8f9fa; border-radius: 6px; padding: 10px 14px;
                       text-align: center; border: 1px solid #dee2e6; }
        .metric-card .label { font-size: 11px; color: #6c757d; text-transform: uppercase;
                              letter-spacing: 0.5px; }
        .metric-card .value { font-size: 18px; font-weight: 600; color: #2c3e50; }
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
                      border-radius: 6px; padding: 8px 14px; margin-bottom: 14px;
                      font-size: 12px; color: #5c4813; line-height: 1.5; }
        .ruo-tag { display: inline-block; background: #a3231b; color: #fff;
                   font-weight: 700; font-size: 10.5px; letter-spacing: 0.6px;
                   border-radius: 3px; padding: 1px 7px; margin-right: 8px; }
        .about-block { border-top: 1px solid #dee2e6; margin-top: 24px;
                       padding-top: 12px; font-size: 12px; color: #6c757d; }
        .about-block h6 { font-weight: 600; color: #2c3e50; font-size: 13px;
                          margin-bottom: 4px; }
        .about-block p { margin-bottom: 8px; }
        .adv-panel { border: 1px solid #dee2e6; border-radius: 6px;
                     margin-bottom: 14px; background: #fff; }
        .adv-panel > summary { cursor: pointer; padding: 8px 12px;
                               font-weight: 600; color: #2c3e50;
                               list-style: revert; font-size: 14px; }
        .adv-panel[open] > summary { border-bottom: 1px solid #dee2e6; }
        .adv-body { padding: 10px 12px 4px; }
        .session-row .btn { font-size: 12px; }
        .chem-hint { font-size: 12px; color: #6c757d; }
        .summary-sticky { position: sticky; top: 0; z-index: 20;
                          background: #fff; padding: 6px 0 4px;
                          border-bottom: 1px solid #ecf0f1; }
        .landing-placeholder { border: 1px dashed #dee2e6; border-radius: 8px;
                               padding: 18px 20px; margin: 14px 0; background: #fbfbfc; }
        .landing-placeholder h5 { font-weight: 600; color: #2c3e50; }
        .landing-placeholder ul { padding-left: 20px; margin-bottom: 0; font-size: 13px; }
        .landing-placeholder li { margin-bottom: 6px; }
        .metric-card.placeholder .value { color: #adb5bd; }
      "))),

      ## -- Input section --
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
               "are also written directly to this folder (created if it ",
               "doesn't exist) when you click Run. \"Browse...\" browses the ",
               "filesystem of the machine R is running on -- if you're ",
               "running this app locally (RStudio / Rscript on your own ",
               "machine, as in the Quick Start), that's your own filesystem; ",
               "under a remote/hosted deployment it would be the server's."),
        selectInput("conj5", "5' conjugate", choices = .conj5_choices, selected = "none"),
        selectInput("conj3", "3' conjugate", choices = .conj3_choices, selected = "none")
      ),

      ## -- Metabolite generation section --
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
                    value = DEFAULT_PIPELINE_PARAMS$min_frag_len, min = 1, max = 20),
        tags$p(style = "font-size: 11px; color: #6c757d; margin-top: 6px;",
               "The Charge Envelopes sheet computes a full isotope pattern for ",
               "every metabolite x PS-oxidation level. Endonuclease fragments ",
               "and a high \u201cMax PS oxid.\u201d both multiply that count -- for ",
               "long sequences this can take several minutes.")
      ),

      ## -- Mass & isotope section --
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
                      "targets \u2014 each metabolite is modeled at every oxidation level ",
                      "from 0 up to this cap."),
                    value = DEFAULT_PIPELINE_PARAMS$max_oxid, min = 0, max = 30))
        ),
        numericInput("h_offset",
                    label = .with_info("Envelope offset (Da)",
                      "0 = standard [M-zH]^z- charge envelope. Nonzero only to ",
                      "reproduce a legacy or lab-specific mass convention."),
                    value = DEFAULT_PIPELINE_PARAMS$h_offset, step = 0.001),
        checkboxInput("use_envipat", "Use enviPat for isotopes",
                      value = DEFAULT_PIPELINE_PARAMS$use_envipat),
        tags$p(class = "chem-hint", style = "margin-top: -8px;",
               "Unchecked uses a built-in binomial-convolution approximation; ",
               "install enviPat (", tags$code("install.packages(\"enviPat\")"),
               ") for higher-accuracy, literature isotope-abundance patterns.")
      ),

      ## -- Orbitrap Exploris acquisition method section (advanced; collapsed
      ## by default -- most runs never touch these past their defaults) --
      tags$details(class = "adv-panel",
        tags$summary("Orbitrap Acquisition Method"),
        tags$div(class = "adv-body",
          numericInput("method_length",
                      label = .with_info("Method length (min)",
                        "Total LC-MS run time; sets the RT window end (Start/End Time) ",
                        "for every row in the MS1 inclusion and MS2 PRM target lists."),
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
                        "parameter \u2014 optimize per method and per instrument before ",
                        "relying on it."),
                      value = DEFAULT_PIPELINE_PARAMS$hcd_nce, min = 1, max = 200),
          fluidRow(
            column(6, numericInput("ms1_target_cap",
                      label = .with_info("Max MS1 targets",
                        "PRM/targeted duty cycle degrades quickly as target count grows; ",
                        "raising this cap can lengthen the instrument cycle time beyond ",
                        "what's practical for the method."),
                      value = DEFAULT_PIPELINE_PARAMS$ms1_target_cap, min = 1, max = 150000)),
            column(6, numericInput("ms2_target_cap", "Max MS2 targets",
                      value = DEFAULT_PIPELINE_PARAMS$ms2_target_cap, min = 1, max = 150000))
          ),
          tags$p(style = "font-size: 11px; color: #6c757d; margin-top: -6px;",
                 "MS1/MS2 lists match the Orbitrap Exploris Method Editor's ",
                 "Targeted Mass filter table (m/z & z, Start/End Time mode). ",
                 "PRM duty cycle degrades quickly with target count, so MS2 ",
                 "defaults to a narrower charge range than the full MS1 envelope.")
        )
      ),

      ## -- MS matching section (optional) --
      tags$div(class = "sidebar-section",
        tags$h5("MS Matching (optional)"),
        checkboxInput("enable_ms", "Enable MS matching",
                      value = DEFAULT_PIPELINE_PARAMS$enable_ms),
        conditionalPanel(
          condition = "input.enable_ms == true",
          fileInput("ms_file", "Upload MS file (.mzML, .mzXML, .csv)",
                    accept = c(".mzML", ".mzXML", ".mzml", ".mzxml", ".csv", ".txt")),
          numericInput("ppm_tol", "MS1 tolerance (ppm)",
                      value = DEFAULT_PIPELINE_PARAMS$ppm_tol, min = 1, max = 50),
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

      ## -- Batch MS processing section (optional, advanced; collapsed by
      ## default). Independent of the single-file "MS Matching" section
      ## above: uploads multiple raw files, runs them through the parallel
      ## Python charge-envelope deconvolution pipeline
      ## (inst/python/oligomet_deconv/), and matches/confirms/compares
      ## across samples. See R/batch_ms_processing.R and R/statistics.R.
      tags$details(class = "adv-panel",
        tags$summary("Batch MS Processing (optional)"),
        tags$div(class = "adv-body",
          checkboxInput("enable_batch", "Enable batch processing",
                        value = DEFAULT_PIPELINE_PARAMS$enable_batch),
          conditionalPanel(
            condition = "input.enable_batch == true",
            fileInput("batch_files", "Upload raw files (.mzML, .mzXML)", multiple = TRUE,
                      accept = c(".mzML", ".mzml", ".mzXML", ".mzxml")),
            DT::DTOutput("sample_meta_table"),
            tags$p(style = "font-size: 11px; color: #6c757d; margin-top: 4px;",
                   "One row per uploaded file. Fill in Group (2+ groups) or ",
                   "Timepoint (time series) before running -- leave both blank ",
                   "to only extract and match features, with no statistics."),
            checkboxInput("batch_run_ms2", "Confirm hits with MS2",
                          value = DEFAULT_PIPELINE_PARAMS$batch_run_ms2),
            fluidRow(
              column(6, numericInput("batch_n_workers", "Parallel workers",
                        value = DEFAULT_PIPELINE_PARAMS$batch_n_workers, min = 1, max = 64)),
              column(6, numericInput("batch_deconv_ppm", "Deconv mass tol (ppm)",
                        value = DEFAULT_PIPELINE_PARAMS$batch_deconv_ppm, min = 1, max = 100))
            )
          )
        )
      ),

      ## -- Run button --
      tags$hr(),
      actionButton("run", "Run Pipeline", class = "btn-primary btn-lg w-100"),

      ## -- Session save/load ---------------------------------------------
      ## Captures the sequence, every parameter above, and the Custom
      ## Chemistry table as one JSON file -- "what exactly produced this
      ## inclusion list" is a reproducibility question for a tool whose
      ## outputs feed instrument methods. Uploaded MS/batch files are not
      ## included (re-upload those after loading a session).
      tags$div(style = "margin-top: 10px;", class = "session-row",
        fluidRow(
          column(6, downloadButton("dl_session", "Save session",
                                   class = "btn-sm btn-outline-secondary w-100")),
          column(6, tags$div(
            fileInput("load_session_file", NULL, accept = ".json",
                      buttonLabel = "Load session...", placeholder = "")))
        ),
        tags$p(style = "font-size: 10.5px; color: #6c757d; margin-top: -8px;",
               "Saves/restores the sequence, custom chemistry, and all ",
               "parameters above as a .json file -- not uploaded MS/batch files.")
      )
    ),

    ## ---- Main panel: manual entry + status + summary + downloads ----------
    mainPanel(
      width = 8,

      ## -- Manual sequence entry ---------------------------------------------
      ## The other way in: instead of assembling triplet notation by hand in
      ## the sidebar, type the three lines a chemical analysis file (or a
      ## BioPharma Finder sequence entry) already gives you, and let the app
      ## assemble the triplet. See inst/help/SEQUENCE_GUIDE.md.
      tags$div(class = "manual-entry",
        tags$h5("Manual sequence entry"),
        tags$p(class = "hint",
               "Type the three lines from your chemical analysis file. ",
               "Bases and sugars need one code per position; linkages sit ",
               "between positions, so there is one fewer of them. A single ",
               "sugar or linkage code is applied to every position ",
               "(\"e\" = all MOE, \"s\" = all phosphorothioate). Separate ",
               "multi-character codes with commas or dashes ",
               "(\"MOE-MOE-d-d\"). Submit fills in the sequence box on the ",
               "left; then click Run Pipeline. New to this? Open ",
               tags$strong("Help & guides"), " below -- the sequence guide ",
               "walks through reading a chemical analysis file and filling ",
               "in these three fields."),
        fluidRow(
          column(4, textInput("man_bases", "Bases (5'->3')", value = "",
                              placeholder = "TSASTTTSATAATGSTGG")),
          column(4, textInput("man_sugars", "Sugars", value = "",
                              placeholder = "eeeeeeeeeeeeeeeeee")),
          column(4, textInput("man_linkages", "Linkages", value = "",
                              placeholder = "sssssssssssssssss"))
        ),
        fluidRow(
          column(3, actionButton("man_submit", "Submit",
                                 class = "btn-primary w-100")),
          column(3, actionButton("man_example", "Fill example",
                                 class = "btn-outline-secondary w-100")),
          column(3, actionButton("man_clear", "Clear",
                                 class = "btn-outline-secondary w-100")),
          column(3, downloadButton("dl_fasta", "BPF FASTA",
                                   class = "btn-outline-secondary w-100"))
        ),
        htmlOutput("man_feedback")
      ),

      ## -- Custom chemistry (advanced) -----------------------------------------
      ## Moved out of the sidebar and into the main panel: it's an advanced,
      ## occasional-use feature, and editing DT table cells in a narrow
      ## width-4 sidebar column was cramped. Full main-panel width gives the
      ## table room to breathe, and collapsing it by default keeps it out of
      ## the way of the primary Input -> Run flow for the common case where
      ## no custom chemistry is needed.
      tags$details(class = "help-panel",
        tags$summary("Custom Chemistry (advanced)"),
        tags$div(class = "help-body",
          tags$p(class = "chem-hint",
                 "Add custom base/sugar/linkage/conjugate entries not in the ",
                 "standard dictionary. Empty rows are ignored; a row with a ",
                 "code but no valid formula is flagged in the Status column ",
                 "and blocks Run until it's fixed or cleared."),
          DT::dataTableOutput("custom_chem"),
          tags$div(style = "height: 6px;"),
          fluidRow(
            column(3, actionButton("add_row", "Add Row", class = "btn-sm btn-outline-primary w-100")),
            column(3, actionButton("remove_row", "Remove Row", class = "btn-sm btn-outline-secondary w-100"))
          )
        )
      ),

      ## -- Help & guides -----------------------------------------------------
      ## The three guides render straight from inst/help/, so they are
      ## available to installed users, not just to people with a checkout.
      ## Collapsed by default: they are long, and the dashboard should not
      ## open on a wall of documentation.
      tags$details(class = "help-panel",
        tags$summary("Help & guides"),
        tags$div(class = "help-body",
          tabsetPanel(
            id = "help_tabs",
            tabPanel("Quick start",     uiOutput("help_quickstart")),
            tabPanel("Sequence guide",  uiOutput("help_sequence")),
            tabPanel("Modifications",   uiOutput("help_modifications"))
          )
        )
      ),

      ## Status / log
      verbatimTextOutput("status", placeholder = TRUE),

      ## -- Pre-run placeholder --------------------------------------------
      ## Before the first Run, everything below (metrics, plots, downloads)
      ## is conditionalPanel-hidden -- correct behavior, since there's
      ## nothing to show yet, but it leaves a lot of dead whitespace on
      ## first load. This teaches the output structure instead.
      conditionalPanel(
        condition = "!(input.run > 0 && output.status_ready == 'true')",
        tags$div(class = "landing-placeholder",
          tags$h5("What this produces"),
          tags$p(class = "chem-hint",
                 "Enter a sequence (left) and click Run Pipeline. A run computes:"),
          fluidRow(
            column(3, tags$div(class = "metric-card placeholder",
              tags$div(class = "label", "Formula"), tags$div(class = "value", "—"))),
            column(2, tags$div(class = "metric-card placeholder",
              tags$div(class = "label", "Mono Mass (Da)"), tags$div(class = "value", "—"))),
            column(2, tags$div(class = "metric-card placeholder",
              tags$div(class = "label", "Avg Mass (Da)"), tags$div(class = "value", "—"))),
            column(1, tags$div(class = "metric-card placeholder",
              tags$div(class = "label", "Length"), tags$div(class = "value", "—"))),
            column(2, tags$div(class = "metric-card placeholder",
              tags$div(class = "label", "Metabolites"), tags$div(class = "value", "—"))),
            column(2, tags$div(class = "metric-card placeholder",
              tags$div(class = "label", "PRM Entries"), tags$div(class = "value", "—")))
          ),
          tags$div(style = "height: 8px;"),
          tags$ul(
            tags$li(tags$strong("4 plots"), " -- charge envelope, truncation series, ",
                    "isotope pattern, and PS-oxidation series for the parent oligo, ",
                    "plus an interactive MS2 spectrum explorer."),
            tags$li(tags$strong("Downloads"), " -- an Excel workbook, an HTML report, ",
                    "a PRM inclusion list, Orbitrap Exploris MS1/MS2 acquisition method ",
                    "lists, and MGF/MSP spectral libraries (bundled as one .zip, or ",
                    "individually)."),
            tags$li(tags$strong("MS matching results"), " (if enabled) -- MS1 match ",
                    "counts and, for batch processing, per-sample feature tables and ",
                    "group/time-series statistics.")
          )
        )
      ),

      ## Summary dashboard (hidden until run completes)
      conditionalPanel(
        condition = "input.run > 0 && output.status_ready == 'true'",

        tags$hr(),

        ## Metrics cards: sticky to the top of the viewport, so they stay
        ## visible while scrolling down through the plot tabs and downloads
        ## below -- otherwise they scroll out of view exactly when you'd
        ## want to check them against a plot or a download.
        tags$div(class = "summary-sticky",
          tags$h4("Summary"),

          ## Metrics row
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

          ## MS metrics row (shown only if MS matching was enabled)
          conditionalPanel(
            condition = "input.enable_ms == true",
            fluidRow(
              column(4, tags$div(class = "metric-card",
                tags$div(class = "label", "MS1 Matches"), tags$div(class = "value", textOutput("m_ms1_matches")))),
              column(4, tags$div(class = "metric-card",
                tags$div(class = "label", "Annotated Mets"), tags$div(class = "value", textOutput("m_annotated")))),
              column(4, tags$div(class = "metric-card",
                tags$div(class = "label",
                  .with_info("Putative IDs (ppm match)",
                    "Matched within the configured MS1 ppm tolerance, with a ",
                    "consistent isotope pattern and charge state. Not confirmed by ",
                    "retention time. MS/MS-confirmed only when Batch MS Processing's ",
                    "\"Confirm hits with MS2\" was enabled and a spectrum was found ",
                    "near the precursor -- otherwise this counts m/z matches only, ",
                    "not orthogonally-confirmed identifications.")),
                tags$div(class = "value", textOutput("m_confident"))))
            ),
            tags$div(style = "height: 10px;")
          )
        ),

        ## Plot tabs
        tabsetPanel(
          tabPanel("Charge Envelope", plotOutput("plot_envelope", height = "400px")),
          tabPanel("Truncation Series", plotOutput("plot_truncation", height = "400px")),
          tabPanel("Isotope Pattern", plotOutput("plot_isotope", height = "400px")),
          tabPanel("Oxidation Series", plotOutput("plot_oxidation", height = "400px")),
          tabPanel("MS2 Explorer",
            tags$div(style = "padding-top: 12px;",
              fluidRow(
                column(8, tags$p(style = "font-size: 12px; color: #6c757d;",
                  "Builds the full MS2 spectral library (one spectrum per ",
                  "metabolite x precursor charge) for interactive browsing -- ",
                  "the same computation as the MS2 library download below, ",
                  "so it can take a while for long sequences or a wide charge ",
                  "range.")),
                column(4, actionButton("build_ms2_explorer",
                  "Build / Refresh Library",
                  class = "btn-sm btn-outline-primary w-100"))
              ),
              conditionalPanel(
                condition = "output.ms2_explorer_ready == 'true'",
                tags$p(style = "font-size: 12px; color: #6c757d;",
                  "Filter by precursor m/z with the range slider, or by ",
                  "charge with the dropdown, in the column headers below. ",
                  "Click a row to plot its predicted MS2 spectrum."),
                DT::dataTableOutput("ms2_explorer_table"),
                tags$hr(),
                uiOutput("ms2_explorer_title"),
                fluidRow(
                  column(7, plotOutput("ms2_explorer_plot", height = "320px")),
                  column(5, DT::dataTableOutput("ms2_explorer_peaks"))
                )
              )
            )
          ),
          tabPanel("Batch Results",
            conditionalPanel(
              condition = "output.batch_ready == 'true'",
              tags$div(style = "padding-top: 12px;",
                tags$p(style = "font-size: 12px; color: #6c757d;",
                  "One row per matched (sample, metabolite, charge, adduct) hit. ",
                  "Confirmation columns are populated only when MS2 confirmation ",
                  "was enabled and a spectrum was found near that precursor."),
                DT::DTOutput("batch_matches_table"),
                tags$div(style = "height: 8px;"),
                downloadButton("dl_batch_tsv", "Download combined features (.tsv)",
                               class = "btn-outline-primary")
              )
            ),
            conditionalPanel(
              condition = "output.batch_ready != 'true'",
              tags$p(style = "padding-top: 12px; color: #6c757d;",
                     "Enable batch processing, upload files, and run the pipeline to see results here.")
            )
          ),
          tabPanel("Unidentified Peaks",
            conditionalPanel(
              condition = "output.batch_ready == 'true'",
              tags$div(style = "padding-top: 12px;",
                tags$p(style = "font-size: 12px; color: #6c757d;",
                  "Deconvolved features that matched no theoretical metabolite ",
                  "within the MS1 tolerance -- retained for QC and follow-up, ",
                  "not discarded."),
                DT::DTOutput("unmatched_table"),
                tags$div(style = "height: 8px;"),
                downloadButton("dl_unmatched_csv", "Download unidentified peaks (.csv)",
                               class = "btn-outline-primary")
              )
            )
          ),
          tabPanel("Statistics",
            conditionalPanel(
              condition = "output.stats_ready == 'true'",
              tags$div(style = "padding-top: 12px;",
                uiOutput("stats_met_selector"),
                plotOutput("plot_stats_main", height = "360px"),
                DT::DTOutput("stats_table"),
                tags$div(style = "height: 8px;"),
                downloadButton("dl_stats_csv", "Download statistics table (.csv)",
                               class = "btn-outline-primary")
              )
            ),
            conditionalPanel(
              condition = "output.stats_ready != 'true'",
              tags$p(style = "padding-top: 12px; color: #6c757d;",
                     "Fill in Group or Timepoint in the batch sample table to see comparisons here.")
            )
          )
        ),

        tags$hr(),

        ## Download buttons
        tags$h5("Download"),
        fluidRow(
          column(12, downloadButton("dl_all", "Download All Outputs (.zip)",
                                    class = "btn-dark w-100"))
        ),
        tags$p(style = "font-size: 11px; color: #6c757d; margin: 4px 0 10px;",
               "Bundles every file below -- workbook, report, PRM list, ",
               "acquisition method lists, and spectral libraries -- into a ",
               "single .zip for local download. Individual files are also ",
               "available one at a time below."),
        fluidRow(
          column(4, downloadButton("dl_workbook", "Excel Workbook (.xlsx)",
                                   class = "btn-success w-100")),
          column(4, downloadButton("dl_report", "HTML Report (.html)",
                                   class = "btn-info w-100")),
          column(4, downloadButton("dl_prm", "PRM Inclusion List (.csv)",
                                   class = "btn-warning w-100"))
        ),
        tags$div(style = "height: 10px;"),
        tags$h5("Orbitrap Exploris Acquisition Method"),
        fluidRow(
          column(4, downloadButton("dl_ms1_inclusion", "MS1 Inclusion List (.csv)",
                                   class = "btn-outline-secondary w-100")),
          column(4, downloadButton("dl_ms2_prm", "MS2 PRM Target List (.csv)",
                                   class = "btn-outline-secondary w-100")),
          column(4, downloadButton("dl_frag_ref", "MS2 Fragment Reference (.csv)",
                                   class = "btn-outline-secondary w-100"))
        ),
        tags$p(style = "font-size: 11px; color: #6c757d; margin-top: 6px;",
               "The first two import directly into the Method Editor's Targeted ",
               "Mass filter. The fragment reference is for interpreting spectra ",
               "after acquisition (e.g. Skyline transitions) -- it isn't an ",
               "acquisition input and won't import into the Method Editor."),

        tags$div(style = "height: 10px;"),
        tags$h5("Spectral Libraries"),
        fluidRow(
          column(3, downloadButton("dl_ms1_mgf", "MS1 (.mgf)",
                                   class = "btn-outline-primary w-100")),
          column(3, downloadButton("dl_ms1_msp", "MS1 (.msp)",
                                   class = "btn-outline-primary w-100")),
          column(3, tagAppendAttributes(
            downloadButton("dl_ms2_mgf", "MS2 (.mgf) ⚠",
                          class = "btn-outline-primary w-100"),
            title = "MS2 intensities are placeholders (rule-based heuristic) -- match on m/z only, not for intensity-weighted scoring.")),
          column(3, tagAppendAttributes(
            downloadButton("dl_ms2_msp", "MS2 (.msp) ⚠",
                          class = "btn-outline-primary w-100"),
            title = "MS2 intensities are placeholders (rule-based heuristic) -- match on m/z only, not for intensity-weighted scoring."))
        ),
        tags$p(style = "font-size: 11px; color: #a3231b; font-weight: 600; margin-top: 6px; margin-bottom: 2px;",
               "⚠ MS2 intensities are placeholders, not measured or fitted values -- ",
               "match on m/z only. Intensity-weighted dot-product scoring (e.g. ",
               "MS-DIAL) against the MS2 library will produce meaningless confidence ",
               "values."),
        tags$p(style = "font-size: 11px; color: #6c757d; margin-top: 6px;",
               "Theoretical libraries for MS-DIAL, mzVault/Compound Discoverer, ",
               "MZmine, matchms and similar. MS1 spectra hold the isotope ",
               "cluster of each metabolite at each charge state, with real ",
               "relative abundances. MS2 spectra hold the McLuckey fragment ",
               "ions for each precursor charge state, with a rule-based ",
               "relative-intensity heuristic (PS-linkage lability, MOE-PS vs ",
               "DNA-PS dominant ion series, purine base-loss lability -- see ",
               "fragment_intensity_weight() in R/fragments.R) -- it is not ",
               "fit to any measured spectra, so match on m/z and treat it as ",
               "a coarse ranking, not a predicted abundance for quantitative ",
               "dot-product scoring."),
        tags$div(style = "height: 20px;")
      ),

      ## ---- About and full disclaimer --------------------------------------
      ## Always visible, not gated behind the run: the disclaimer has to be
      ## readable before anyone downloads anything.
      tags$div(id = "about-section", class = "about-block",
        tags$h6("About"),
        ## htmltools joins sibling arguments with a space, so anything that
        ## must sit flush against its neighbour (a colon, a bracket) is
        ## pasted into one string rather than passed as its own argument.
        tags$p(
          tags$strong("OligoMet Profiler"), HTML("&mdash;"),
          textOutput("about_version", inline = TRUE), tags$br(),
          paste0(OLIGOMET_AUTHOR_ROLE, ": ", OLIGOMET_AUTHOR),
          HTML(paste0("(<a href=\"mailto:", OLIGOMET_AUTHOR_EMAIL, "\">",
                      OLIGOMET_AUTHOR_EMAIL, "</a>)")), tags$br(),
          paste0(OLIGOMET_AUTHOR_TITLE, ", ", OLIGOMET_AUTHOR_AFFILIATION,
                 " -- an independent personal project, not a ",
                 OLIGOMET_AUTHOR_AFFILIATION, " product."), tags$br(),
          tags$a(href = OLIGOMET_URL, target = "_blank", rel = "noopener",
                 OLIGOMET_URL),
          HTML("&mdash; released under the MIT licence.")),
        tags$h6("Disclaimer"),
        lapply(OLIGOMET_DISCLAIMER, tags$p)
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
    plots = list(), ready = FALSE, status_text = "Enter a sequence and click Run Pipeline.\n",
    batch_features = NULL, batch_ms_results = NULL,
    sample_meta = NULL, stats_results = NULL
  )

  ## ---- Batch sample metadata table (group/timepoint assignment) -------------
  # Seeded from the uploaded batch_files' original filenames; edited in place
  # via DT's edit feature. Blank group/timepoint columns mean "extract and
  # match only, no statistics" -- see the Run handler below.
  batch_meta_data <- reactiveVal(data.frame(
    sample = character(0), group = character(0), timepoint = character(0),
    stringsAsFactors = FALSE))

  observeEvent(input$batch_files, {
    samples <- tools::file_path_sans_ext(input$batch_files$name)
    batch_meta_data(data.frame(sample = samples, group = "", timepoint = "",
                                stringsAsFactors = FALSE))
  })

  output$sample_meta_table <- DT::renderDT({
    DT::datatable(batch_meta_data(), editable = TRUE, rownames = FALSE,
                   options = list(dom = "t", paging = FALSE, scrollY = "180px"))
  })
  observeEvent(input$sample_meta_table_cell_edit, {
    df <- batch_meta_data()
    edit <- input$sample_meta_table_cell_edit
    df[edit$row, edit$col + 1] <- edit$value
    batch_meta_data(df)
  })

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
             sprintf("Built a %d-mer (%d linkages) and loaded it into the sequence box. Click Run Pipeline.",
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
    rv$status_text <- paste0(
      "Loaded session from '", input$load_session_file$name, "'. Re-upload any ",
      "MS/batch files (not saved in a session), review parameters, then click ",
      "Run Pipeline.\n")
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
  output$help_sequence      <- renderUI(.help_ui("SEQUENCE_GUIDE.md"))
  output$help_modifications <- renderUI(.help_ui("MODIFICATIONS.md"))

  ## ---- Status output -------------------------------------------------------
  output$status <- renderText({ rv$status_text })

  # Hidden flag for conditionalPanel
  output$status_ready <- reactive({ if (rv$ready) "true" else "false" })
  outputOptions(output, "status_ready", suspendWhenHidden = FALSE)

  ## ---- Run pipeline --------------------------------------------------------
  observeEvent(input$run, {
    # Reset
    rv$ready <- FALSE
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

    # Wrap the whole pipeline in a top-level tryCatch: any error anywhere
    # below (not just the points with their own tryCatch) gets reported
    # here instead of silently stopping the reactive with nothing visible
    # in the status panel -- previously an uncaught error from any single
    # step would make the app look like Run had done nothing at all.
    tryCatch({

    withProgress(message = "Running pipeline...", value = 0, {

      prog <- progress_tracker(c(
        "Building dictionary" = 1,
        "Parsing sequence" = 1,
        "Generating metabolites" = 2,
        "Computing masses" = 1,
        "Generating fragment ions" = 2,
        "Generating PRM list" = 1,
        "MS data import and matching" = 3,
        "Batch MS processing (parallel)" = 10,
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
        generate_metabolites(spec, opts = met_opts),
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

      # Step 7: Optional MS matching
      progress_next(prog)
      ms_results <- NULL
      if (input$enable_ms && !is.null(input$ms_file)) {
        incProgress(0.55, detail = "Importing MS data")
        ms_path <- input$ms_file$datapath
        ms_results <- tryCatch({
          if (grepl("\\.mzML$|\\.mzml$|\\.mzXML$|\\.mzxml$", ms_path, ignore.case = TRUE)) {
            ms_data <- parse_mzml(ms_path)
            ms1_features <- extract_ms1_features(ms_data$ms1, ppm = input$ppm_tol)
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
      if (input$enable_batch && !is.null(input$batch_files) && nrow(input$batch_files) > 0) {
        incProgress(0.58, detail = "Batch deconvolution (parallel)")
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

          deconv <- run_batch_deconvolution(
            input$batch_files$datapath, precursor_watchlist = watchlist_path,
            mass_tol_ppm = input$batch_deconv_ppm, n_workers = input$batch_n_workers,
            progress = function(msg) incProgress(0, detail = msg))

          # Shiny renames uploads to random tmp paths; recover the sample
          # name from the ORIGINAL filename, not the datapath basename.
          name_map <- stats::setNames(tools::file_path_sans_ext(input$batch_files$name),
                                       tools::file_path_sans_ext(basename(input$batch_files$datapath)))
          feats <- read_batch_features(deconv$features_path)
          feats$sample <- unname(name_map[feats$sample])
          ms2 <- if (!is.null(deconv$ms2_path)) read_batch_ms2(deconv$ms2_path) else NULL
          if (!is.null(ms2) && nrow(ms2) > 0) ms2$sample <- unname(name_map[ms2$sample])

          incProgress(0, detail = "Matching + MS2 confirmation")
          batch_results <- annotate_metabolites_batch(
            mets, feats, ms2, dict = dict, ppm_tol = input$ppm_tol, z_range = z_range,
            adducts = adducts, max_oxid = input$max_oxid, h_offset = input$h_offset,
            n_iso = input$n_iso, use_envipat = input$use_envipat,
            frag_tol_ppm = input$frag_tol_ppm, frag_z_range = 1:input$frag_z_max)

          list(features = feats, results = batch_results, meta = batch_meta_data())
        }, error = function(e) {
          rv$status_text <- paste0(rv$status_text,
            "WARNING: Batch MS processing failed: ", conditionMessage(e), "\n")
          NULL
        })

        if (!is.null(batch_out)) {
          rv$batch_features <- batch_out$features
          rv$batch_ms_results <- batch_out$results
          rv$sample_meta <- batch_out$meta

          meta <- batch_out$meta
          has_group <- !is.null(meta$group) && any(nzchar(meta$group))
          has_time <- !is.null(meta$timepoint) && any(nzchar(meta$timepoint))
          if (nrow(batch_out$results$ms1_matches) > 0 && (has_group || has_time)) {
            abund <- build_abundance_matrix(batch_out$results$ms1_matches)
            stats_res <- tryCatch({
              if (has_time) {
                sm <- meta[, c("sample", "timepoint")]
                sm$timepoint <- suppressWarnings(as.numeric(sm$timepoint))
                long <- abundance_long(abund, sm[!is.na(sm$timepoint), ])
                list(mode = "time_series", result = compare_time_series(long))
              } else {
                sm <- meta[nzchar(meta$group), c("sample", "group")]
                groups <- unique(sm$group)
                long <- abundance_long(abund, sm)
                if (length(groups) == 2) {
                  list(mode = "two_group", result = compare_two_groups(long, groups[1], groups[2]))
                } else if (length(groups) > 2) {
                  list(mode = "multi_group", result = compare_multi_groups(long, groups))
                } else NULL
              }
            }, error = function(e) {
              rv$status_text <- paste0(rv$status_text,
                "WARNING: statistical comparison failed: ", conditionMessage(e), "\n")
              NULL
            })
            rv$stats_results <- stats_res
          }
        }
      }

      # Step 8: Build workbook
      incProgress(0.65, detail = "Building Excel workbook")
      progress_next(prog)
      build_opts <- list(
        z_range = z_range, n_iso = input$n_iso,
        max_oxid = input$max_oxid, h_offset = input$h_offset,
        use_envipat = input$use_envipat, include_internal = TRUE,
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
                        stats_results = rv$stats_results)
      }, error = function(e) {
        rv$status_text <- paste0("ERROR building workbook: ", conditionMessage(e), "\n")
        NULL
      })
      if (is.null(wb_path)) return()
      rv$wb_path <- wb_path

      # Step 9: Build report
      incProgress(0.75, detail = "Building report")
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

      # Step 10: Generate plots
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

      # Optional: save copies directly to a user-chosen local folder,
      # for local RStudio/Rscript sessions where relying on the browser's
      # download flow isn't convenient.
      saved_to <- NULL
      out_dir <- trimws(input$output_dir %||% "")
      if (nzchar(out_dir)) {
        saved_to <- tryCatch({
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
      rv$status_text <- status
      rv$ready <- TRUE

    })  # end withProgress

    }, error = function(e) {
      rv$status_text <- paste0(
        "ERROR: unexpected failure -- ", conditionMessage(e), "\n",
        "Check the R console for a full traceback.\n")
      rv$ready <- FALSE
    })
  })  # end observeEvent(input$run)

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

  ## ---- Batch MS Processing / Unidentified Peaks / Statistics tabs -----------
  output$batch_ready <- reactive({
    if (!is.null(rv$batch_ms_results) && nrow(rv$batch_ms_results$ms1_matches) >= 0) "true" else "false"
  })
  outputOptions(output, "batch_ready", suspendWhenHidden = FALSE)

  output$stats_ready <- reactive({
    if (!is.null(rv$stats_results) && !is.null(rv$stats_results$result)) "true" else "false"
  })
  outputOptions(output, "stats_ready", suspendWhenHidden = FALSE)

  output$batch_matches_table <- DT::renderDT({
    req(rv$batch_ms_results)
    m <- rv$batch_ms_results$ms1_matches
    DT::datatable(m, filter = "top", rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$unmatched_table <- DT::renderDT({
    req(rv$batch_ms_results)
    DT::datatable(rv$batch_ms_results$unmatched, filter = "top", rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

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

  output$plot_stats_main <- renderPlot({
    sr <- rv$stats_results
    req(sr, input$stats_met_select)
    if (sr$mode == "two_group") {
      plot_volcano(sr$result)
    } else if (sr$mode == "time_series") {
      abund <- build_abundance_matrix(rv$batch_ms_results$ms1_matches)
      long <- abundance_long(abund, rv$sample_meta[, c("sample", "timepoint")])
      plot_trend(long, input$stats_met_select)
    } else {
      abund <- build_abundance_matrix(rv$batch_ms_results$ms1_matches)
      long <- abundance_long(abund, rv$sample_meta[, c("sample", "group")])
      plot_group_boxplot(long, input$stats_met_select)
    }
  })

  output$stats_table <- DT::renderDT({
    DT::datatable(.stats_display_table(), rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$dl_batch_tsv <- downloadHandler(
    filename = function() paste0(input$output_prefix, "_batch_features.tsv"),
    content = function(file) {
      req(rv$batch_features)
      utils::write.table(rv$batch_features, file, sep = "\t", row.names = FALSE, quote = FALSE)
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
