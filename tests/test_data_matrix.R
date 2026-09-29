# test_data_matrix.R -- validate pre-processed data-matrix import/export
# (R/data_matrix.R), calibration curves that fit with no Group/Timepoint
# design (R/statistics.R), and degradation relative to a reference
# timepoint (R/degradation.R), against the bundled example in
# inst/extdata/data_matrix_example/ and small hand-built tables.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R", "ms_matching.R",
             "batch_ms_processing.R", "degradation.R", "statistics.R",
             "blank_correction.R", "multivariate.R", "data_matrix.R")) {
  source(file.path(.pkg_root, "R", .f))
}

check <- function(label, ok) {
  cat(sprintf("%-78s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}

cat("==== Data matrix suite validation ====\n\n")
ex <- file.path(.pkg_root, "inst", "extdata", "data_matrix_example")

## ---- Import: wide and long, auto-detected --------------------------------
cat("--- read_data_matrix() ---\n")
w <- read_data_matrix(file.path(ex, "data_matrix_wide.csv"))
l <- read_data_matrix(file.path(ex, "data_matrix_long.csv"))
check("wide file auto-detected as wide", w$format == "wide")
check("long file auto-detected as long", l$format == "long")
check("wide and long give the same number of (metabolite, sample) values",
      nrow(w$matches) == nrow(l$matches))
check("blank cells are dropped as not detected (M13 missing in active t0)",
      !any(w$matches$met_id == "M13" & grepl("^active_t0_", w$matches$sample)))
check("long file carries its own sample information",
      !is.null(l$sample_meta) && all(c("group", "timepoint", "sample_type", "concentration") %in% names(l$sample_meta)))
check("match columns downstream functions need are present",
      all(c("met_id", "met_name", "kind", "sample", "intensity") %in% names(w$matches)))

# Unnamed identifier column, thousands separators, a raw-file extension on
# a sample header, and a stray text column.
df <- data.frame(check.names = FALSE, stringsAsFactors = FALSE,
                 Compound = c("PARENT", "N-1"), note = c("a", "b"),
                 `s1.mzML` = c("1,000", "200"), s2 = c(900, NA))
r <- read_data_matrix(df)
check("first text column becomes the identifier when no id header exists",
      setequal(unique(r$matches$met_id), c("PARENT", "N-1")))
check("thousands separators parse and file extensions are stripped from samples",
      r$matches$intensity[r$matches$sample == "s1" & r$matches$met_id == "PARENT"] == 1000)
check("a non-numeric extra column is reported, not treated as a sample",
      any(grepl("note", r$notes)) && !"note" %in% r$matches$sample)

## ---- Excel + TSV ------------------------------------------------------------
wide_df <- utils::read.csv(file.path(ex, "data_matrix_wide.csv"), check.names = FALSE)
xl <- tempfile(fileext = ".xlsx"); openxlsx::write.xlsx(wide_df, xl)
tsv <- tempfile(fileext = ".tsv")
utils::write.table(wide_df, tsv, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
check(".xlsx reads to the same values as .csv", nrow(read_data_matrix(xl)$matches) == nrow(w$matches))
check(".tsv reads to the same values as .csv", nrow(read_data_matrix(tsv)$matches) == nrow(w$matches))

## ---- Sample information -----------------------------------------------------
cat("\n--- read_sample_info() / normalize_sample_type() ---\n")
si <- read_sample_info(file.path(ex, "sample_info.csv"))
check("sample info has the five standard columns",
      all(c("sample", "group", "timepoint", "sample_type", "concentration") %in% names(si)))
check("sample-type aliases map onto the controlled vocabulary",
      identical(normalize_sample_type(c("Standard", "QC", "Quality Control", "Blank", "Matrix Blank", "cal", "foo")),
                c("standard", "quality_control", "quality_control", "reagent_blank", "matrix_blank", "standard", NA)))
alias_path <- tempfile(fileext = ".csv")
utils::write.csv(data.frame(`Sample Name` = "x.raw", Time = "4", Type = "Std", Conc = "5", check.names = FALSE),
                 alias_path, row.names = FALSE)
al <- read_sample_info(alias_path)
check("header aliases (Sample Name / Time / Type / Conc) are recognized",
      al$sample == "x" && al$timepoint == "4" && al$sample_type == "standard" && al$concentration == "5")

## ---- Annotation ---------------------------------------------------------------
m <- annotate_matrix_metabolites(w$matches)
check("kind column from the file is kept", all(m$kind[m$met_id == "M01"] == "parent"))
check("guess_parent_met_id() finds the parent", guess_parent_met_id(m) == "M01")
m2 <- annotate_matrix_metabolites(r$matches)
check("PARENT id is inferred as kind = parent without a library",
      all(m2$kind[m2$met_id == "PARENT"] == "parent") && all(m2$kind[m2$met_id == "N-1"] == "unknown"))
m3 <- annotate_matrix_metabolites(r$matches, parent_met_id = "N-1")
check("an explicit parent choice overrides inference (exactly one parent)",
      identical(unique(m3$met_id[m3$kind == "parent"]), "N-1"))
lib <- list(list(id = "PARENT", name = "my oligo", kind = "parent"),
            list(id = "N-1", name = "my oligo 3' N-1", kind = "exo_3p"))
m4 <- annotate_matrix_metabolites(r$matches, lib)
check("names and kinds come from a generated library when ids match",
      all(m4$kind[m4$met_id == "N-1"] == "exo_3p") && all(m4$met_name[m4$met_id == "N-1"] == "my oligo 3' N-1"))

## ---- Calibration on a matrix, incl. a design with no Group/Timepoint --------
cat("\n--- calibration curves from a data matrix ---\n")
q <- quantify_metabolites(m, si, absolute_met_ids = c("M01", "M02"),
                          mode = "time_series", weighting = "1/x2")
check("both calibration curves fit with r^2 > 0.99",
      all(vapply(q$calibration_curves, function(c) c$r_squared, numeric(1)) > 0.99))
qc <- q$absolute[q$absolute$sample_type == "quality_control", ]
check("QC back-calculated within +/-15% of nominal", nrow(qc) == 6 && all(abs(qc$percent_re) < 15))
check("every non-calibrated metabolite is relative-quantified",
      setequal(unique(q$relative$met_id), c("M03", "M04", "M12", "M13")))

std_only <- si[si$sample_type == "standard", ]
std_only$group <- ""; std_only$timepoint <- ""
q2 <- quantify_metabolites(m[m$sample %in% std_only$sample, ], std_only, absolute_met_ids = "M01",
                           mode = "group", control_group = NULL, weighting = "1/x2")
check("standards-only run (no control group) still fits the calibration curve",
      !is.null(q2$calibration_curves$M01$model) && q2$calibration_curves$M01$r_squared > 0.99)
check("...and says why relative quantification was skipped",
      length(q2$notes) == 1 && grepl("control group", q2$notes))

## ---- Degradation relative to the reference timepoint ------------------------
cat("\n--- degradation_vs_reference() ---\n")
deg <- degradation_summary(m, sample_meta = si)
check("standards/QC/blanks excluded from degradation", nrow(deg$per_sample) == 24)
dr <- degradation_vs_reference(deg$per_sample)
check("one row per (group, timepoint)", nrow(dr) == 8)
act <- dr[dr$group == "active", ]
check("reference timepoint is 100% parent remaining", act$pct_parent_remaining[act$timepoint == 0] == 100)
check("active arm loses parent over time (monotonic)", all(diff(act$pct_parent_remaining) < 0))
check("active arm 48 h: about 76% parent loss (known first-order k = 0.03/h)",
      abs(act$pct_parent_loss[act$timepoint == 48] - 100 * (1 - exp(-0.03 * 48))) < 5)
ctl <- dr[dr$group == "heat_inactivated", ]
check("heat-inactivated arm stays within 10% of reference", all(abs(ctl$pct_parent_remaining - 100) < 10))
dr24 <- degradation_vs_reference(deg$per_sample, reference_timepoint = 24)
check("explicit reference timepoint is honored",
      all(dr24$pct_parent_remaining[dr24$timepoint == 24] == 100))
check("no timepoint column -> empty result, not an error",
      nrow(degradation_vs_reference(deg$per_sample[, c("sample", "parent_signal")])) == 0)
invisible(plot_degradation_vs_reference(dr)); invisible(plot_degradation_vs_reference(data.frame()))

## ---- Export round trip ---------------------------------------------------------
cat("\n--- export_data_matrix() round trip ---\n")
wide <- export_data_matrix(m, format = "wide")
long <- export_data_matrix(m, si, format = "long")
rt_w <- read_data_matrix(wide); rt_l <- read_data_matrix(long)
check("wide export reads back as wide with identical values",
      rt_w$format == "wide" && nrow(rt_w$matches) == nrow(m) &&
      isTRUE(all.equal(sum(rt_w$matches$intensity), sum(m$intensity))))
check("long export reads back as long with sample information",
      rt_l$format == "long" && nrow(rt_l$matches) == nrow(m) && !is.null(rt_l$sample_meta))

cat("\nAll data matrix checks passed.\n")
