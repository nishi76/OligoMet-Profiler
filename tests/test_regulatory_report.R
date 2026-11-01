# test_regulatory_report.R -- regulatory report generator (R/regulatory_report.R):
# HTML always builds from a run context; the .docx builds when officer/rsvg are
# present.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "structure_annotation.R", "statistics.R",
             "degradation.R", "fragments.R", "export_spectral.R", "mirror_plot.R", "batch_ms_processing.R", "xic_quant.R",
             "regulatory_report.R", "report_evidence.R")) {
  source(file.path(.pkg_root, "R", .f))
}
check <- function(label, ok) {
  cat(sprintf("%-78s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Regulatory report suite ====\n\n")

spec <- list(bases = strsplit("TCGTCATGACTTGCAG", "")[[1]],
             sugars = c(rep("e", 5), rep("d", 6), rep("e", 5)),
             linkages = c(rep("s", 15), NA), conj5 = "none", conj3 = "none", n = 16)
mets <- generate_metabolites(spec, list(oligo_name = "NWB-204", endo = TRUE,
                                        endo_sites = "gap", max_3p = 2, max_5p = 2))
pid <- mets[[1]]$id
e3 <- Filter(function(m) m$kind == "exo_3p", mets)[[1]]$id
mk <- function(s, mid, kind, a) data.frame(met_id = mid, met_name = mid, kind = kind,
  sample = s, area = a, intensity = a * 10, ppm_error = 0.8, z = 5, adduct = "H",
  stringsAsFactors = FALSE)
matches <- do.call(rbind, list(
  mk("p1", pid, "parent", 1000), mk("p2", pid, "parent", 400),
  mk("p1", e3, "exo_3p", 120), mk("p2", e3, "exo_3p", 260),
  mk("c1", pid, "parent", 990), mk("c2", pid, "parent", 900),
  mk("c1", e3, "exo_3p", 10), mk("c2", e3, "exo_3p", 40)))
smeta <- data.frame(sample = c("p1", "p2", "c1", "c2"),
  matrix = c("plasma", "plasma", "CSF", "CSF"),
  timepoint = c("1", "24", "1", "24"), sample_type = "unknown", stringsAsFactors = FALSE)

ctx <- list(
  details = list(sponsor = "Northwind Biosciences", cro = "Meridian Bioanalytical",
                 study_title = "Biotransformation of NWB-204", study_number = "MBS-2026-0417",
                 test_article = "NWB-204", species = "cynomolgus monkey", route = "intrathecal"),
  spec = spec, mets = mets, dict = STANDARD_DICT, matches = matches,
  abundance = parent_relative_abundance(matches, "area"),
  degradation = degradation_summary(matches, sample_meta = smeta),
  sample_meta = smeta,
  params = list(ppm_tol = 8, adducts = c("H", "Na", "K", "NH4"),
                mods = c("deamination"), z_min = 3, z_max = 12, max_oxid = 6, pct_min = 0.1))

html <- build_regulatory_report_html(ctx)
check("HTML report builds and closes", grepl("</html>", html, fixed = TRUE))
check("front page carries sponsor, CRO and study number",
      grepl("Northwind Biosciences", html) && grepl("Meridian Bioanalytical", html) &&
        grepl("MBS-2026-0417", html))
check("test article identity is filled from the library (formula + mass)",
      grepl("Monoisotopic mass", html) && grepl("Molecular formula", html))
check("methods table shows the 0.1% reporting threshold",
      grepl("0.1% of full-length", html) || grepl("0.1", html, fixed = TRUE))
check("results embed at least one figure (base64 PNG)",
      grepl("data:image/png;base64", html, fixed = TRUE))
check("results embed a structure map (inline SVG)",
      grepl("<svg", html, fixed = TRUE))
check("methods document the quantitation (adduct handling, aggregation)",
      grepl("Signal aggregation and quantitation", html) &&
        grepl("deprotonated [M-zH]", html, fixed = TRUE) &&
        grepl("Adduct handling", html, fixed = TRUE))
check("each reported metabolite gets an identification-evidence table",
      grepl("Identification evidence for", html, fixed = TRUE) &&
        grepl("Sequence (5'-&gt;3')", html, fixed = TRUE))
check("without raw m/z, the characterization falls back to theoretical values",
      grepl("Theoretical precursor m/z", html, fixed = TRUE) &&
        grepl("Theoretical isotope distribution", html, fixed = TRUE))
check("% of parent excludes non-study samples (blank with tiny parent signal)", local({
  m2 <- rbind(matches, mk("blank", pid, "parent", 1), mk("blank", e3, "exo_3p", 50))
  sm2 <- rbind(smeta, data.frame(sample = "blank", matrix = "plasma", timepoint = "", sample_type = "reagent_blank"))
  ab <- .reg_abundance(m2, .reg_sample_frame(unique(m2$sample), sm2), 0.1)
  ab$per_met$max_pct[ab$per_met$met_id == e3] < 100
}))
check("empty context still produces a valid HTML shell (placeholders)",
      grepl("</html>", build_regulatory_report_html(list()), fixed = TRUE))

if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("rsvg", quietly = TRUE)) {
  f <- tempfile(fileext = ".docx")
  build_regulatory_report_docx(ctx, f)
  check(".docx file is written and non-trivial", file.exists(f) && file.info(f)$size > 10000)
} else {
  cat("(.docx path skipped: officer/flextable/rsvg not installed)\n")
}

cat("\nAll regulatory report checks passed.\n")
