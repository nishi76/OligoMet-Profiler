# test_stress_example.R -- the bundled forced-degradation example data
# (inst/extdata/stress_example) imports through the data-matrix path and
# degrades as designed, and carries no supplied analyte name.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "internal_standard.R",
             "data_matrix.R", "statistics.R", "degradation.R"))
  source(file.path(.pkg_root, "R", .f))
check <- function(label, ok) {
  cat(sprintf("%-72s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Stress example-data suite ====\n\n")

dir <- file.path(.pkg_root, "inst", "extdata", "stress_example")
long <- file.path(dir, "data_matrix_long.csv")
info <- file.path(dir, "sample_info.csv")
check("bundled example files exist", file.exists(long) && file.exists(info))

raw <- readLines(long, warn = FALSE)
check("no supplied analyte name leaks into the example (zz001/ION582)",
      !any(grepl("zz001|ion582|zzz001", raw, ignore.case = TRUE)))

res <- read_data_matrix(long, format = "long")
check("imports as a long matrix", identical(res$format, "long"))
check("27 metabolites x 18 samples (3 arms x 6 timepoints)",
      length(unique(res$matches$met_id)) == 27 && length(res$samples) == 18)
check("embedded sample metadata carries group and timepoint",
      all(c("group", "timepoint") %in% names(res$sample_meta)))
check("three stress arms present",
      setequal(unique(res$sample_meta$group), c("acid", "base", "thermal")))

ann <- annotate_matrix_metabolites(res$matches, mets = NULL)
deg <- degradation_summary(ann, sample_meta = res$sample_meta, signal_col = "intensity")
check("degradation summary runs (one row per sample)", nrow(deg$per_sample) == 18)

vr <- degradation_vs_reference(deg$per_sample)
check("time 0 is the dosing reference, not a predose control", !any(vr$is_predose))
rem <- function(g) vr$pct_parent_remaining[vr$group == g & vr$timepoint == 168]
check("base degrades fastest, thermal slowest at 168 h (by design)",
      rem("base") < rem("acid") && rem("acid") < rem("thermal"))
check("parent is 100% at the dosing timepoint in every arm",
      all(abs(vr$pct_parent_remaining[vr$timepoint == 0] - 100) < 1e-6))

cat("\nAll stress example-data checks passed.\n")
