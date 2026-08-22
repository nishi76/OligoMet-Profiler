# =============================================================================
# app.R -- OligoMet Profiler
# Shiny app wrapper for the oligonucleotide metabolite identification
# pipeline. Runs on any oligonucleotide -- standard chemistry works out of
# the box, and the custom chemistry table below covers anything else.
#
# Run with:
#   Rscript -e 'shiny::runApp(".", launch.browser=TRUE)'
#   or in RStudio: open this file and click "Run App"
#
# The app sources all pipeline modules from the same directory and provides
# an interactive interface for:
#   - Sequence input (triplet or OligoDistiller notation, auto-detected)
#   - Custom chemistry overrides (editable table)
#   - Full parameter control (truncation depth, charge range, oxidation, etc.)
#   - Optional MS data upload for matching
#   - Summary dashboard with 4 plots + download buttons
# =============================================================================

## ---- Module sourcing -------------------------------------------------------
.module_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else "."
}, error = function(e) ".")

# shiny::runApp() sets the working directory to the app folder for the
# duration of the app, so "." resolves correctly in that case. This check
# only guards against the rare case where neither method finds the modules.
if (!file.exists(file.path(.module_dir, "R", "chemistry_dict.R"))) {
  stop("Cannot locate pipeline modules (chemistry_dict.R not found). ",
       "Run this app with shiny::runApp() pointed at the folder that ",
       "contains app.R, or set the working directory there first.")
}

# Install DT if missing
if (!requireNamespace("DT", quietly = TRUE)) {
  install.packages("DT", repos = "https://cloud.r-project.org", quiet = TRUE)
}

source(file.path(.module_dir, "R", "progress_utils.R"))
source(file.path(.module_dir, "R", "chemistry_dict.R"))
source(file.path(.module_dir, "R", "oligo_io.R"))
source(file.path(.module_dir, "R", "metabolites.R"))
source(file.path(.module_dir, "R", "mass_isotope.R"))
source(file.path(.module_dir, "R", "fragments.R"))
source(file.path(.module_dir, "R", "ms_matching.R"))
source(file.path(.module_dir, "R", "build_workbook.R"))
source(file.path(.module_dir, "R", "build_report.R"))
source(file.path(.module_dir, "R", "export_acquisition.R"))

library(shiny)
library(DT)

## ---- Example sequences (illustrative, not the assumed input) --------------
# Two worked examples: a generic gapmer built entirely from standard
# dictionary codes, and the ION337 reference sequence (validated to <1 ppm
# against a published mass -- see chemistry_dict.R). Neither is required;
# paste any sequence into the box below.
.EXAMPLE_SEQS <- list(
  "Generic 2'MOE/DNA gapmer" =
    "Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm",
  "ION337 reference (validated)" = ION337_TRIPLET
)

## ---- Conjugate options -----------------------------------------------------
.conj5_choices <- c("none", "5'-phosphate", "5'-thiophosphate", "biotin",
                    "cAG_cap", "cAU_cap", "ARCA_cap", "mCAP",
                    "GalNAc", "GalNAc3",
                    "cholesterol", "C6", "C12", "TEG", "FAM", "Cy3")
.conj3_choices <- c("none", "3'-phosphate", "3'-cyclophos", "3'-thiophosphate",
                    "GalNAc", "GalNAc3", "GalNAc3_triantennary",
                    "cholesterol", "C6", "C12", "TEG", "FAM", "Cy3")

## ---- Custom chemistry table initial data -----------------------------------
.custom_chem_init <- data.frame(
  Code    = c("", "", ""),
  Formula = c("", "", ""),
  Name    = c("", "", ""),
  Type    = c("base", "sugar", "linkage"),
  Attach  = c("add", "add", "add"),
  stringsAsFactors = FALSE
)

## =============================================================================
## UI
## =============================================================================
ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

  titlePanel("OligoMet Profiler"),

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
      "))),

      ## -- Input section --
      tags$div(class = "sidebar-section",
        tags$h5("Input"),
        textAreaInput("seq", "Sequence (triplet or OligoDistiller)",
                      value = .EXAMPLE_SEQS[[1]], rows = 3,
                      placeholder = "e.g. Ge-uAn-sGn-sSn-... or OH-Am*-Gm*-...-OH"),
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
          column(4, shinyFiles::shinyDirButton("browse_output_dir", "Browse...", "Choose a folder",
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

      ## -- Custom chemistry section --
      tags$div(class = "sidebar-section",
        tags$h5("Custom Chemistry"),
        tags$p(style = "font-size: 12px; color: #6c757d;",
               "Add custom base/sugar/linkage/conjugate entries. Empty rows are ignored."),
        DT::dataTableOutput("custom_chem", height = "180px"),
        fluidRow(
          column(6, actionButton("add_row", "Add Row", class = "btn-sm btn-outline-primary")),
          column(6, actionButton("remove_row", "Remove Row", class = "btn-sm btn-outline-secondary"))
        )
      ),

      ## -- Metabolite generation section --
      tags$div(class = "sidebar-section",
        tags$h5("Metabolite Generation"),
        fluidRow(
          column(6, numericInput("max_3p", "Max 3' trunc.", value = 10, min = 0, max = 50)),
          column(6, numericInput("max_5p", "Max 5' trunc.", value = 10, min = 0, max = 50))
        ),
        checkboxInput("endo", "Include endonuclease fragments", value = TRUE),
        radioButtons("endo_sites", "Endo cleavage sites",
                     choices = c("All positions" = "all", "DNA gap only" = "gap"),
                     selected = "all", inline = TRUE),
        numericInput("min_frag_len", "Min fragment length (nt)", value = 3, min = 1, max = 20),
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
          column(6, numericInput("z_min", "Min charge z", value = 3, min = 1, max = 20)),
          column(6, numericInput("z_max", "Max charge z", value = 12, min = 1, max = 30))
        ),
        fluidRow(
          column(6, numericInput("n_iso", "Isotope peaks", value = 5, min = 1, max = 20)),
          column(6, numericInput("max_oxid", "Max PS oxid.", value = 6, min = 0, max = 30))
        ),
        numericInput("h_offset", "Envelope offset (Da)", value = 0, step = 0.001),
        checkboxInput("use_envipat", "Use enviPat for isotopes", value = TRUE)
      ),

      ## -- Orbitrap Exploris acquisition method section --
      tags$div(class = "sidebar-section",
        tags$h5("Orbitrap Acquisition Method"),
        numericInput("method_length", "Method length (min)", value = 30, min = 1, max = 999),
        fluidRow(
          column(6, numericInput("ms2_z_min", "MS2 charge z min", value = 4, min = 1, max = 30)),
          column(6, numericInput("ms2_z_max", "MS2 charge z max", value = 7, min = 1, max = 30))
        ),
        numericInput("hcd_nce", "HCD NCE (%)", value = 20, min = 1, max = 200),
        fluidRow(
          column(6, numericInput("ms1_target_cap", "Max MS1 targets", value = 500, min = 1, max = 150000)),
          column(6, numericInput("ms2_target_cap", "Max MS2 targets", value = 300, min = 1, max = 150000))
        ),
        tags$p(style = "font-size: 11px; color: #6c757d; margin-top: -6px;",
               "MS1/MS2 lists match the Orbitrap Exploris Method Editor's ",
               "Targeted Mass filter table (m/z & z, Start/End Time mode). ",
               "HCD NCE is a starting value, not a validated instrument ",
               "parameter -- optimize per method. PRM duty cycle degrades ",
               "quickly with target count, so MS2 defaults to a narrower ",
               "charge range than the full MS1 envelope.")
      ),

      ## -- MS matching section (optional) --
      tags$div(class = "sidebar-section",
        tags$h5("MS Matching (optional)"),
        checkboxInput("enable_ms", "Enable MS matching", value = FALSE),
        conditionalPanel(
          condition = "input.enable_ms == true",
          fileInput("ms_file", "Upload MS file (.mzML, .mzXML, .csv)",
                    accept = c(".mzML", ".mzXML", ".mzml", ".mzxml", ".csv", ".txt")),
          numericInput("ppm_tol", "MS1 tolerance (ppm)", value = 10, min = 1, max = 50),
          checkboxGroupInput("adducts", "Adducts",
                             choices = c("H", "Na", "K", "NH4"),
                             selected = c("H", "Na", "K", "NH4"), inline = TRUE),
          fluidRow(
            column(6, numericInput("frag_tol_ppm", "Fragment tol (ppm)", value = 25, min = 5, max = 100)),
            column(6, numericInput("frag_z_max", "Fragment max z", value = 3, min = 1, max = 5))
          )
        )
      ),

      ## -- Run button --
      tags$hr(),
      actionButton("run", "Run Pipeline", class = "btn-primary btn-lg w-100")
    ),

    ## ---- Main panel: status + summary + downloads -------------------------
    mainPanel(
      width = 8,

      ## Status / log
      verbatimTextOutput("status", placeholder = TRUE),

      ## Summary dashboard (hidden until run completes)
      conditionalPanel(
        condition = "input.run > 0 && output.status_ready == 'true'",

        tags$hr(),
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
              tags$div(class = "label", "Confident IDs"), tags$div(class = "value", textOutput("m_confident"))))
          ),
          tags$div(style = "height: 10px;")
        ),

        ## Plot tabs
        tabsetPanel(
          tabPanel("Charge Envelope", plotOutput("plot_envelope", height = "400px")),
          tabPanel("Truncation Series", plotOutput("plot_truncation", height = "400px")),
          tabPanel("Isotope Pattern", plotOutput("plot_isotope", height = "400px")),
          tabPanel("Oxidation Series", plotOutput("plot_oxidation", height = "400px"))
        ),

        tags$hr(),

        ## Download buttons
        tags$h5("Download"),
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
        tags$div(style = "height: 20px;")
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
    plots = list(), ready = FALSE, status_text = "Enter a sequence and click Run Pipeline.\n"
  )

  ## ---- Custom chemistry table ----------------------------------------------
  custom_chem_data <- reactiveVal(.custom_chem_init)

  # Load a selected example sequence into the input box
  observeEvent(input$load_example, {
    sel <- input$example_seq
    if (nzchar(sel) && sel %in% names(.EXAMPLE_SEQS)) {
      updateTextAreaInput(session, "seq", value = .EXAMPLE_SEQS[[sel]])
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
  volumes <- c(Home = path.expand("~"), shinyFiles::getVolumes()())
  shinyFiles::shinyDirChoose(input, "browse_output_dir", roots = volumes, session = session)
  observeEvent(input$browse_output_dir, {
    sel <- input$browse_output_dir
    if (is.list(sel) && !is.null(sel$path)) {
      path <- shinyFiles::parseDirPath(volumes, sel)
      if (length(path) == 1 && nzchar(path)) {
        updateTextInput(session, "output_dir", value = path)
      }
    }
  })

  output$custom_chem <- DT::renderDataTable({
    DT::datatable(custom_chem_data(), editable = "cell", rownames = FALSE,
                  options = list(dom = "t", paging = FALSE, ordering = FALSE,
                                 autoWidth = TRUE,
                                 columnDefs = list(list(width = "60px", targets = 0),
                                                   list(width = "100px", targets = 1),
                                                   list(width = "80px", targets = 3))),
                  selection = "none")
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

  ## ---- Status output -------------------------------------------------------
  output$status <- renderText({ rv$status_text })

  # Hidden flag for conditionalPanel
  output$status_ready <- reactive({ rv$ready })
  outputOptions(output, "status_ready", suspendWhenHidden = FALSE)

  ## ---- Run pipeline --------------------------------------------------------
  observeEvent(input$run, {
    # Reset
    rv$ready <- FALSE
    rv$status_text <- "Running...\n"

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
        "Building Excel workbook" = 15,
        "Building report" = 4,
        "Generating plots" = 2
      ))

      # Step 1: Build dictionary with custom overrides
      incProgress(0.05, detail = "Building dictionary")
      progress_next(prog)
      custom_overrides <- list()
      cd <- custom_chem_data()
      for (i in seq_len(nrow(cd))) {
        code <- trimws(cd$Code[i])
        formula_str <- trimws(cd$Formula[i])
        if (nchar(code) == 0 || nchar(formula_str) == 0) next
        entry <- list(
          formula = formula_str,
          name = if (nchar(trimws(cd$Name[i])) > 0) trimws(cd$Name[i]) else code
        )
        type <- cd$Type[i]
        if (type == "conjugate") {
          entry$attach <- cd$Attach[i]
        }
        entry$kind <- type
        custom_overrides[[code]] <- entry
      }
      dict <- build_dictionary(overrides = custom_overrides)
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
                        console_tracker = prog)
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
          "  Confident IDs: ", n_conf, "\n")
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
}

## ---- Run app ---------------------------------------------------------------
shinyApp(ui = ui, server = server)
