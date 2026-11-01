# =============================================================================
# regulatory_report.R
# Regulatory-format metabolite-identification study report, populated from a
# live run. The scientific sections (test article, methods, results, structure
# maps, figures and tables) are filled from the metabolite library, the batch
# match table, quantification and the degradation summary; the administrative
# and GLP fields (sponsor, CRO, study number, dates, signatures) come from a
# study-details list and fall back to bracketed placeholders when not supplied.
#
# Two outputs from one assembled model (see .reg_model()):
#   build_regulatory_report_html()  -- self-contained HTML (inline SVG + PNG)
#   build_regulatory_report_docx()  -- native Word via officer (headers/footer,
#                                      page numbers, flextable tables, images)
# =============================================================================

## ---- Study-details fields and placeholder helper ---------------------------
.REG_DETAIL_FIELDS <- c(
  "sponsor", "sponsor_address", "cro", "cro_address", "study_title",
  "study_number", "project_code", "protocol", "test_article",
  "test_article_lot", "is_name", "study_director", "species", "route",
  "report_status", "report_date")

.reg_get <- function(details, key, default = NULL) {
  v <- details[[key]]
  if (is.null(v) || !nzchar(trimws(as.character(v)[1]))) {
    if (!is.null(default)) return(default)
    return(sprintf("[%s]", gsub("_", " ", key)))
  }
  trimws(as.character(v)[1])
}

## ---- Figure / image helpers ------------------------------------------------
.reg_fig_png <- function(plot, width = 7.2, height = 3.6, dpi = 150) {
  if (is.null(plot)) return(NULL)
  f <- tempfile(fileext = ".png")
  ok <- tryCatch({
    ggplot2::ggsave(f, plot = plot, width = width, height = height, dpi = dpi,
                    units = "in", bg = "white"); TRUE
  }, error = function(e) FALSE)
  if (ok && file.exists(f)) f else NULL
}

.reg_svg_png <- function(svg, width_px = 900) {
  if (is.null(svg) || !nzchar(svg)) return(NULL)
  f <- tempfile(fileext = ".png")
  ok <- tryCatch({
    rsvg::rsvg_png(charToRaw(svg), file = f, width = width_px); TRUE
  }, error = function(e) FALSE)
  if (ok && file.exists(f)) f else NULL
}

.reg_b64img <- function(png_path) {
  if (is.null(png_path) || !file.exists(png_path)) return(NULL)
  paste0("data:image/png;base64,", base64enc::base64encode(png_path))
}

.reg_xml_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

## ---- Assemble the report model ---------------------------------------------
# ctx fields (all optional; the report degrades gracefully):
#   details   named list of administrative/GLP fields (see .REG_DETAIL_FIELDS)
#   spec      parent oligo_spec (bases/sugars/linkages/conj5/conj3/n)
#   mets      metabolite library (list of metabolite objects)
#   dict      chemistry dictionary
#   matches   batch ms1_matches (met_id, met_name, kind, ppm_error, z, ...)
#   abundance parent_relative_abundance() table
#   degradation  degradation_summary() result (for degradation_by_condition)
#   quant     quantify_metabolites() result
#   params    list(ppm_tol, adducts, mods, z_min, z_max, max_oxid, pct_min)
#' Assemble the regulatory-report model from a run context
#'
#' Computes the structured content (front-matter fields, summary counts, test
#' article and methods tables, metabolite abundance table, degradation table,
#' structure maps and figures) shared by the HTML and Word renderers.
#'
#' @param ctx A context list; see the file header for the recognized fields.
#'   All are optional and the report fills gaps with bracketed placeholders.
#' @return A list consumed by [build_regulatory_report_html()] and
#'   [build_regulatory_report_docx()].
#' @seealso [build_regulatory_report_html()], [build_regulatory_report_docx()]
#' @export
.reg_model <- function(ctx) {
  details <- ctx$details %||% list()
  params <- ctx$params %||% list()
  dict <- ctx$dict %||% STANDARD_DICT
  threshold <- params$pct_min %||% 0.1

  # --- test article identity (from the parent spec / library) ---
  ta_rows <- list()
  parent_met <- NULL
  if (!is.null(ctx$mets)) parent_met <- Find(function(m) identical(m$kind, "parent"), ctx$mets)
  seq_str <- if (!is.null(ctx$spec)) .reg_sequence_string(ctx$spec) else NA
  mass_str <- NA; formula_str <- NA
  if (!is.null(parent_met)) {
    info <- tryCatch(metabolite_mass_info(parent_met, dict), error = function(e) NULL)
    if (!is.null(info)) { mass_str <- sprintf("%.3f", info$mono_mass); formula_str <- info$formula_str }
  } else if (!is.null(ctx$spec)) {
    f <- tryCatch(assemble_oligo_formula(ctx$spec$bases, ctx$spec$sugars, ctx$spec$linkages,
                    ctx$spec$conj5, ctx$spec$conj3, dict), error = function(e) NULL)
    if (!is.null(f)) { mass_str <- sprintf("%.3f", formula_mass(f, TRUE)); formula_str <- format_formula(f) }
  }
  n_mer <- if (!is.null(ctx$spec)) ctx$spec$n else if (!is.null(parent_met)) parent_met$n else NA
  test_article <- data.frame(
    Attribute = c("Name / code", "Modality", "Length", "Sequence (5'->3')",
                  "Molecular formula (free acid)", "Monoisotopic mass (Da)", "Lot",
                  "Internal standard"),
    Value = c(.reg_get(details, "test_article"),
              "oligonucleotide (per sequence)",
              if (is.na(n_mer)) "[length]" else sprintf("%d-mer", n_mer),
              if (is.na(seq_str)) "[sequence]" else seq_str,
              if (is.na(formula_str)) "[formula]" else formula_str,
              if (is.na(mass_str)) "[mass]" else mass_str,
              .reg_get(details, "test_article_lot"),
              .reg_get(details, "is_name")),
    stringsAsFactors = FALSE)

  # --- methods parameters ---
  adducts <- params$adducts %||% c("H", "Na", "K", "NH4")
  mods <- params$mods %||% character(0)
  methods <- data.frame(
    Parameter = c("Mass accuracy tolerance", "Charge states", "Adducts",
                  "Max PS->PO desulfurization", "In-place modifications searched",
                  "Metabolite reporting threshold"),
    Setting = c(sprintf("%s ppm", params$ppm_tol %||% 10),
                sprintf("%s to %s", params$z_min %||% 3, params$z_max %||% 12),
                paste(adducts, collapse = ", "),
                as.character(params$max_oxid %||% 6),
                if (length(mods) == 0) "none" else paste(mods, collapse = ", "),
                sprintf(">= %.3g%% of full-length parent", threshold)),
    stringsAsFactors = FALSE)

  # --- quantitation methodology (how signal is aggregated and normalized) ---
  xic_on <- isTRUE(params$xic_enabled)
  signal_basis <- params$signal_basis %||% "raw"
  basis_lab <- switch(signal_basis,
    raw = "raw peak response (area where available, else apex intensity)",
    bcorr = "blank-corrected peak response",
    isr = "internal-standard-normalized response ratio (analyte area / IS area)",
    bcorr_isr = "blank-corrected and internal-standard-normalized response",
    signal_basis)
  agg_lab <- if (xic_on) {
    sprintf("summed extracted-ion chromatogram: the top %s charge states and %s isotopes are co-integrated over one shared retention-time window per metabolite (summation of ions)",
            params$xic_n_charges %||% 5, params$xic_n_iso %||% 5)
  } else {
    "collapsed to one value per metabolite and sample, taken as the maximum response across charge states, adducts and desulfurization levels (the dominant ionization form)"
  }
  gis <- params$global_is %||% "none"
  is_lab <- if (!is.null(gis) && nzchar(gis) && !identical(gis, "none"))
    sprintf("response normalized to the internal standard (%s) in the same sample", gis)
  else "no internal-standard normalization applied"
  quant_methods <- data.frame(
    Aspect = c("Ionization", "Adducts considered", "Adduct handling",
               "Charge-state handling", "Signal aggregation", "Signal basis",
               "Internal standard", "Absolute vs relative", "Reporting threshold"),
    Approach = c(
      "negative-ion electrospray; analytes detected as multiply charged [M-zH]^z-",
      paste(adducts, collapse = ", "),
      "the deprotonated [M-zH] form is the dominant and quantified species; Na/K/NH4 adducts are the same metabolite with a cation for a proton and are used for confirmation, not counted as separate species",
      if (xic_on) sprintf("charge-state envelope (z %s to %s) co-summed per metabolite",
                          params$z_min %||% 3, params$z_max %||% 12)
      else sprintf("charge-state envelope (z %s to %s); the dominant charge state carries the reported response, with cross-charge-state mass agreement used as a confirmation check",
                   params$z_min %||% 3, params$z_max %||% 12),
      agg_lab,
      basis_lab,
      is_lab,
      "parent and selected metabolites quantified against a matrix-matched calibration curve; all other metabolites reported relative to the full-length parent",
      sprintf("metabolites at or above %.3g%% of the full-length parent in at least one sample", threshold)),
    stringsAsFactors = FALSE)
  quant_narrative <- paste0(
    "Each metabolite's response was ", agg_lab, ". ",
    "Adduct forms of a metabolite were not treated as separate species: the deprotonated [M-zH]^z- envelope is the dominant and quantified form, and sodium, potassium and ammonium adducts, when observed, were used to corroborate the identification. ",
    "Signal basis for the reported values was ", basis_lab, "; ", is_lab, ". ",
    "Parent and selected metabolites were quantified against a calibration curve; the remaining metabolites were reported relative to the full-length parent. ",
    sprintf("Only metabolites reaching %.3g%% of the full-length parent in at least one sample are reported.", threshold))

  # --- summary counts ---
  matches <- ctx$matches
  n_identified <- if (!is.null(matches) && nrow(matches) > 0) length(unique(matches$met_id)) else 0L
  ms2c <- ctx$ms2_confirmations
  n_ms2 <- if (!is.null(ms2c) && nrow(ms2c) > 0) length(unique(ms2c$met_id)) else 0L
  abund <- ctx$abundance
  n_floor <- if (!is.null(abund) && nrow(abund) > 0) sum(!is.na(abund$max_pct) & abund$max_pct >= threshold) else n_identified
  matrices <- if (!is.null(ctx$sample_meta) && "matrix" %in% names(ctx$sample_meta)) {
    sort(unique(trimws(ctx$sample_meta$matrix[nzchar(ctx$sample_meta$matrix)])))
  } else character(0)

  # --- metabolite abundance table (reported set) ---
  met_table <- NULL
  if (!is.null(abund) && nrow(abund) > 0) {
    keep <- abund[!is.na(abund$max_pct) & abund$max_pct >= threshold, , drop = FALSE]
    keep <- keep[order(-keep$max_pct), , drop = FALSE]
    ms2_ids <- if (!is.null(ms2c) && nrow(ms2c) > 0) unique(ms2c$met_id) else character(0)
    met_table <- data.frame(
      Metabolite = keep$met_name %||% keep$met_id,
      Class = keep$kind,
      `Max % of parent` = sprintf("%.3g", keep$max_pct),
      MS2 = ifelse(keep$met_id %in% ms2_ids, "confirmed", "-"),
      check.names = FALSE, stringsAsFactors = FALSE)
  }

  # --- degradation by condition (table + figure) ---
  deg_table <- NULL; deg_plot <- NULL; comp_plot <- NULL
  if (!is.null(ctx$degradation)) {
    bc <- tryCatch(degradation_by_condition(ctx$degradation), error = function(e) NULL)
    if (!is.null(bc) && bc$design != "none" && nrow(bc$per_condition) > 0) {
      deg_table <- bc$per_condition
      num <- vapply(deg_table, is.numeric, logical(1))
      deg_table[num] <- lapply(deg_table[num], function(x) round(x, 2))
      deg_plot <- tryCatch(plot_degradation_by_condition(ctx$degradation), error = function(e) NULL)
      comp_plot <- tryCatch(plot_degradation_composition(ctx$degradation, "condition"), error = function(e) NULL)
    }
  }

  # --- structure maps: parent + top reported metabolites ---
  structures <- list()
  if (!is.null(ctx$mets) && !is.null(ctx$spec)) {
    ids <- "parent-first"
    pick_ids <- character(0)
    if (!is.null(parent_met)) pick_ids <- parent_met$id
    if (!is.null(abund) && nrow(abund) > 0) {
      top <- abund[order(-abund$max_pct), , drop = FALSE]
      top <- top$met_id[!is.na(top$max_pct) & top$max_pct >= threshold]
      pick_ids <- unique(c(pick_ids, head(setdiff(top, pick_ids), 4)))
    }
    for (mid in pick_ids) {
      met <- Find(function(m) identical(m$id, mid), ctx$mets)
      if (is.null(met)) next
      svg <- tryCatch(annotate_metabolite_svg(met, parent = ctx$spec, dict = dict),
                      error = function(e) NULL)
      if (!is.null(svg)) structures[[met$name %||% met$id]] <- svg
    }
  }

  list(
    details = details, threshold = threshold,
    test_article = test_article, methods = methods,
    quant_methods = quant_methods, quant_narrative = quant_narrative,
    summary = list(n_identified = n_identified, n_ms2 = n_ms2, n_floor = n_floor,
                   matrices = matrices),
    met_table = met_table, deg_table = deg_table,
    deg_plot = deg_plot, comp_plot = comp_plot, structures = structures)
}

# Render a spec as a 5'->3' triplet string (sugar+base, linkage between).
.reg_sequence_string <- function(spec) {
  n <- spec$n
  res <- vapply(seq_len(n), function(i) paste0(spec$sugars[i], spec$bases[i]), "")
  lk <- spec$linkages
  out <- res[1]
  for (i in seq_len(n - 1)) out <- paste0(out, "-", if (!is.na(lk[i])) lk[i] else "", res[i + 1])
  c5 <- spec$conj5 %||% "none"; c3 <- spec$conj3 %||% "none"
  if (!identical(c5, "none")) out <- paste0(c5, "-", out)
  if (!identical(c3, "none")) out <- paste0(out, "-", c3)
  out
}

## ---- HTML renderer ---------------------------------------------------------
.reg_html_table <- function(df, caption = NULL) {
  if (is.null(df) || nrow(df) == 0) return("")
  th <- paste0("<th>", .reg_xml_escape(names(df)), "</th>", collapse = "")
  rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", .reg_xml_escape(r), "</td>", collapse = ""), "</tr>"))
  cap <- if (!is.null(caption)) paste0("<div class='tbl-cap'>", .reg_xml_escape(caption), "</div>") else ""
  paste0(cap, "<table><thead><tr>", th, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
}

#' Build the regulatory study report as self-contained HTML
#'
#' @param ctx A run context (see [.reg_model()]).
#' @return A single character string: a complete HTML document.
#' @seealso [build_regulatory_report_docx()]
#' @export
build_regulatory_report_html <- function(ctx) {
  m <- .reg_model(ctx)
  d <- m$details
  s <- m$summary
  img <- function(plot, w = 7.2, h = 3.6) {
    uri <- .reg_b64img(.reg_fig_png(plot, w, h))
    if (is.null(uri)) "" else sprintf("<img class='fig' src='%s'/>", uri)
  }
  struct_html <- paste(vapply(names(m$structures), function(nm)
    sprintf("<div class='struct'><div class='struct-name'>%s</div>%s</div>",
            .reg_xml_escape(nm), m$structures[[nm]]), ""), collapse = "")
  matrices_line <- if (length(s$matrices) > 0) paste(s$matrices, collapse = ", ") else "single / unspecified"

  css <- "
    body{font-family:'Segoe UI',system-ui,Helvetica,Arial,sans-serif;color:#1d2733;max-width:900px;margin:0 auto;padding:28px;line-height:1.5}
    h1{font-size:20px} h2{font-size:15px;border-bottom:1px solid #d9dee3;padding-bottom:4px;margin-top:28px}
    .conf{color:#a3231b;font-weight:700;letter-spacing:.05em;font-size:12px}
    .front td{padding:4px 10px;border-bottom:1px solid #eef1f4;font-size:13px}
    table{border-collapse:collapse;width:100%;margin:10px 0;font-size:12.5px}
    th,td{border:1px solid #e0e4e8;padding:5px 8px;text-align:left}
    th{background:#f5f7f9}
    .tbl-cap{font-size:12px;color:#6c757d;margin-top:10px;font-weight:600}
    .fig{max-width:100%;border:1px solid #eef1f4;border-radius:6px;margin:8px 0}
    .struct{margin:10px 0} .struct-name{font-size:12px;font-weight:700;margin-bottom:2px}
    .kpi{display:flex;gap:14px;margin:8px 0;flex-wrap:wrap}
    .kpi div{background:#f5f7f9;border:1px solid #e0e4e8;border-radius:8px;padding:8px 14px}
    .kpi b{font-size:20px;display:block}
    .muted{color:#6c757d;font-size:12px}"

  front <- sprintf("
    <div class='conf'>CONFIDENTIAL</div>
    <h1>%s</h1>
    <table class='front'>
      <tr><td>Sponsor</td><td>%s%s</td></tr>
      <tr><td>Test facility (CRO)</td><td>%s%s</td></tr>
      <tr><td>Study number</td><td>%s</td></tr>
      <tr><td>Sponsor project code</td><td>%s</td></tr>
      <tr><td>Protocol</td><td>%s</td></tr>
      <tr><td>Test article</td><td>%s (lot %s)</td></tr>
      <tr><td>Test system</td><td>%s; %s</td></tr>
      <tr><td>Study Director</td><td>%s</td></tr>
      <tr><td>Report status / date</td><td>%s / %s</td></tr>
    </table>",
    .reg_xml_escape(.reg_get(d, "study_title")),
    .reg_xml_escape(.reg_get(d, "sponsor")),
    { a <- .reg_get(d, "sponsor_address", ""); if (nzchar(a)) paste0(", ", .reg_xml_escape(a)) else "" },
    .reg_xml_escape(.reg_get(d, "cro")),
    { a <- .reg_get(d, "cro_address", ""); if (nzchar(a)) paste0(", ", .reg_xml_escape(a)) else "" },
    .reg_xml_escape(.reg_get(d, "study_number")),
    .reg_xml_escape(.reg_get(d, "project_code")),
    .reg_xml_escape(.reg_get(d, "protocol")),
    .reg_xml_escape(.reg_get(d, "test_article")), .reg_xml_escape(.reg_get(d, "test_article_lot")),
    .reg_xml_escape(.reg_get(d, "species")), .reg_xml_escape(.reg_get(d, "route")),
    .reg_xml_escape(.reg_get(d, "study_director")),
    .reg_xml_escape(.reg_get(d, "report_status", "Draft")),
    .reg_xml_escape(.reg_get(d, "report_date", format(Sys.Date(), "%d %b %Y"))))

  body <- paste0(
    front,
    "<h2>1. Summary</h2>",
    sprintf("<div class='kpi'><div><b>%d</b>metabolites identified</div><div><b>%d</b>MS2-confirmed</div><div><b>%d</b>&ge; %.3g%% of parent</div></div>",
            s$n_identified, s$n_ms2, s$n_floor, m$threshold),
    sprintf("<p>Biotransformation of %s was characterized across %s. Metabolites were identified by accurate mass, isotope pattern, charge-envelope consistency and MS/MS, and are reported relative to the full-length parent. %d metabolites cleared the %.3g%% reporting threshold. [Interpretation to be completed.]</p>",
            .reg_xml_escape(.reg_get(d, "test_article")), matrices_line, s$n_floor, m$threshold),
    "<h2>2. Objectives</h2><p>To identify and characterize the biotransformation products of the test article, determine the degradation kinetics per matrix, and report the abundance of each metabolite relative to the full-length parent. [Expand per protocol.]</p>",
    "<h2>3. Test Article</h2>", .reg_html_table(m$test_article),
    "<h2>5. Materials and Methods</h2>",
    "<p class='muted'>Ion-pairing reversed-phase LC-HRMS/MS, negative ionization. Deconvolution, theoretical metabolite library construction, accurate-mass matching, charge-envelope scoring and MS/MS confirmation as configured below. [Instrument, column, gradient and extraction per SOP.]</p>",
    .reg_html_table(m$methods, "Table. Data-processing parameters."),
    "<h3>5.1 Signal aggregation and quantitation</h3>",
    sprintf("<p>%s</p>", .reg_xml_escape(m$quant_narrative)),
    .reg_html_table(m$quant_methods, "Table. Quantitation methodology."),
    "<h2>6. Results</h2>",
    "<h3>6.1 Degradation by matrix / condition</h3>",
    if (!is.null(m$deg_table)) .reg_html_table(m$deg_table, "Table. Mean % degradation per condition.") else "<p class='muted'>No group/timepoint/matrix design in the sample metadata.</p>",
    img(m$deg_plot), img(m$comp_plot),
    "<h3>6.2 Metabolite profile (>= threshold)</h3>",
    if (!is.null(m$met_table)) .reg_html_table(m$met_table, "Table. Metabolites reported, as % of full-length parent.") else "<p class='muted'>No abundance data available.</p>",
    "<h3>6.3 Structural characterization</h3>",
    if (nzchar(struct_html)) struct_html else "<p class='muted'>Generate a library to embed structure maps.</p>",
    "<h2>7. Discussion</h2><p>[Interpretation of the matrix-specific biotransformation pathways and of any metabolite of toxicological interest.]</p>",
    "<h2>8. Conclusions</h2><p>[Summary conclusion.]</p>",
    "<h2>Appendices</h2><p class='muted'>A. Sequence and chemistry. B. LC gradient. C. Kinetic fits. D. Stability. E. Deviations. F. Annotated MS/MS. G. Full metabolite table. H. Chromatograms. I. Equipment/software. J. Protocol.</p>",
    "<p class='muted'>Generated by OligoMet Profiler. Scientific content is computed from the run; administrative, GLP and interpretive fields in brackets are completed by the study team.</p>")

  paste0("<!doctype html><html><head><meta charset='utf-8'><title>",
         .reg_xml_escape(.reg_get(d, "study_title")), "</title><style>", css, "</style></head><body>",
         body, "</body></html>")
}

## ---- Word (.docx) renderer -------------------------------------------------
.reg_ft <- function(df) {
  ft <- flextable::flextable(df)
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::fontsize(ft, size = 9, part = "all")
  ft <- flextable::autofit(ft)
  ft
}

#' Build the regulatory study report as a Word (.docx) file
#'
#' Native Word output with a running footer (study number and page numbering),
#' flextable tables and embedded figures and structure maps. Requires the
#' `officer`, `flextable` and `rsvg` packages.
#'
#' @param ctx A run context (see [.reg_model()]).
#' @param file Output path for the `.docx`.
#' @return `file`, invisibly.
#' @seealso [build_regulatory_report_html()]
#' @export
build_regulatory_report_docx <- function(ctx, file) {
  for (p in c("officer", "flextable", "rsvg")) {
    if (!requireNamespace(p, quietly = TRUE))
      stop("the .docx report needs the '", p, "' package")
  }
  m <- .reg_model(ctx)
  d <- m$details
  s <- m$summary
  H1 <- "heading 1"; H2 <- "heading 2"
  doc <- officer::read_docx()
  doc <- tryCatch(officer::set_doc_properties(doc, title = .reg_get(d, "study_title")),
                  error = function(e) doc)

  add_h <- function(doc, text, style = H2) officer::body_add_par(doc, text, style = style)
  add_p <- function(doc, text) officer::body_add_par(doc, text, style = "Normal")
  add_tbl <- function(doc, df, cap = NULL) {
    if (is.null(df) || nrow(df) == 0) return(doc)
    if (!is.null(cap)) doc <- officer::body_add_par(doc, cap, style = "Normal")
    flextable::body_add_flextable(doc, .reg_ft(df))
  }
  add_img <- function(doc, plot, w = 6.5, h = 3.2) {
    png <- .reg_fig_png(plot, w, h)
    if (is.null(png)) return(doc)
    officer::body_add_img(doc, png, width = w, height = h)
  }
  add_svg <- function(doc, svg, name) {
    png <- .reg_svg_png(svg, width_px = 1100)
    if (is.null(png)) return(doc)
    doc <- officer::body_add_par(doc, name, style = "Normal")
    info <- tryCatch(magick::image_info(magick::image_read(png)), error = function(e) NULL)
    w <- 6.5; h <- if (!is.null(info)) min(3.2, 6.5 * info$height / info$width) else 1.6
    officer::body_add_img(doc, png, width = w, height = h)
  }

  # Front page
  doc <- officer::body_add_par(doc, "CONFIDENTIAL", style = "Normal")
  doc <- officer::body_add_par(doc, .reg_get(d, "study_title"), style = "heading 1")
  front <- data.frame(Field = c("Sponsor", "Test facility (CRO)", "Study number",
                                 "Sponsor project code", "Protocol", "Test article",
                                 "Test system", "Study Director", "Report status / date"),
    Value = c(paste0(.reg_get(d, "sponsor"), { a <- .reg_get(d, "sponsor_address", ""); if (nzchar(a)) paste0(", ", a) else "" }),
              paste0(.reg_get(d, "cro"), { a <- .reg_get(d, "cro_address", ""); if (nzchar(a)) paste0(", ", a) else "" }),
              .reg_get(d, "study_number"), .reg_get(d, "project_code"), .reg_get(d, "protocol"),
              paste0(.reg_get(d, "test_article"), " (lot ", .reg_get(d, "test_article_lot"), ")"),
              paste0(.reg_get(d, "species"), "; ", .reg_get(d, "route")),
              .reg_get(d, "study_director"),
              paste0(.reg_get(d, "report_status", "Draft"), " / ", .reg_get(d, "report_date", format(Sys.Date(), "%d %b %Y")))),
    stringsAsFactors = FALSE)
  doc <- add_tbl(doc, front)
  doc <- officer::body_add_break(doc)

  doc <- add_h(doc, "1. Summary")
  doc <- add_p(doc, sprintf("Metabolites identified: %d. MS2-confirmed: %d. At or above %.3g%% of the full-length parent: %d.",
                            s$n_identified, s$n_ms2, m$threshold, s$n_floor))
  matrices_line <- if (length(s$matrices) > 0) paste(s$matrices, collapse = ", ") else "single / unspecified"
  doc <- add_p(doc, sprintf("Biotransformation of %s was characterized across %s. Metabolites were identified by accurate mass, isotope pattern, charge-envelope consistency and MS/MS, reported relative to the full-length parent. [Interpretation to be completed.]",
                            .reg_get(d, "test_article"), matrices_line))
  doc <- add_h(doc, "2. Objectives")
  doc <- add_p(doc, "To identify and characterize the biotransformation products of the test article, determine the degradation kinetics per matrix, and report metabolite abundance relative to the full-length parent. [Expand per protocol.]")
  doc <- add_h(doc, "3. Test Article")
  doc <- add_tbl(doc, m$test_article)
  doc <- add_h(doc, "5. Materials and Methods")
  doc <- add_p(doc, "Ion-pairing reversed-phase LC-HRMS/MS, negative ionization. [Instrument, column, gradient and extraction per SOP.]")
  doc <- add_tbl(doc, m$methods, "Table. Data-processing parameters.")
  doc <- add_h(doc, "5.1 Signal aggregation and quantitation")
  doc <- add_p(doc, m$quant_narrative)
  doc <- add_tbl(doc, m$quant_methods, "Table. Quantitation methodology.")
  doc <- add_h(doc, "6. Results")
  doc <- add_h(doc, "6.1 Degradation by matrix / condition")
  if (!is.null(m$deg_table)) doc <- add_tbl(doc, m$deg_table, "Table. Mean % degradation per condition.")
  doc <- add_img(doc, m$deg_plot); doc <- add_img(doc, m$comp_plot)
  doc <- add_h(doc, "6.2 Metabolite profile (>= threshold)")
  if (!is.null(m$met_table)) doc <- add_tbl(doc, m$met_table, "Table. Metabolites reported, as % of full-length parent.")
  doc <- add_h(doc, "6.3 Structural characterization")
  for (nm in names(m$structures)) doc <- add_svg(doc, m$structures[[nm]], nm)
  doc <- add_h(doc, "7. Discussion")
  doc <- add_p(doc, "[Interpretation of the matrix-specific biotransformation pathways and of any metabolite of toxicological interest.]")
  doc <- add_h(doc, "8. Conclusions")
  doc <- add_p(doc, "[Summary conclusion.]")
  doc <- add_p(doc, "Generated by OligoMet Profiler. Scientific content is computed from the run; administrative, GLP and interpretive fields in brackets are completed by the study team.")

  print(doc, target = file)
  invisible(file)
}
