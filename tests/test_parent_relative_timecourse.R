# test_parent_relative_timecourse.R -- "% of initial parent" kinetic
# normalization (R/degradation.R): every species expressed vs the parent at
# the arm's reference timepoint, with a leading predose control excluded.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("palettes.R", "about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R", "ms_matching.R", "batch_ms_processing.R",
             "degradation.R")) {
  source(file.path(.pkg_root, "R", .f))
}
check <- function(label, ok) {
  cat(sprintf("%-72s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Parent-relative time-course suite ====\n\n")

# One arm; timepoints 0 (predose, no drug), 1, 24, 48 h; 2 reps each.
# parent: predose ~1% of t1; t1 = 10000 (reference), t24 = 8800, t48 = 8000.
# short-mer M1 (exo_3p): 0 up to t1, rises to 300 (t24), 500 (t48).
meta <- data.frame(sample = paste0("s", 1:8), group = "armA",
                   timepoint = rep(c(0, 1, 24, 48), each = 2),
                   sample_type = "", stringsAsFactors = FALSE)
mk <- function(id, name, kind, sig) data.frame(
  met_id = id, met_name = name, kind = kind, sample = meta$sample, area = sig,
  stringsAsFactors = FALSE)
matches <- rbind(
  mk("PARENT", "parent oligo", "parent", c(100,100, 10000,10000, 8800,8800, 8000,8000)),
  mk("M1", "3prime N-1", "exo_3p", c(0,0, 0,0, 300,300, 500,500)))

pr <- parent_relative_timecourse(matches, meta)
pp <- pr[pr$met_id == "PARENT", ]; mm <- pr[pr$met_id == "M1", ]

check("reference timepoint is the earliest non-predose (t = 1, not the predose t = 0)",
      unique(pr$reference_timepoint) == 1)
check("predose t = 0 is flagged as a control", isTRUE(pp$is_predose[pp$timepoint == 0]))
check("predose parent reads ~0% of initial (no drug), not 100%",
      pp$pct_of_initial_parent[pp$timepoint == 0] == 1)
check("parent is 100% at the reference timepoint",
      pp$pct_of_initial_parent[pp$timepoint == 1] == 100)
check("parent declines on the fixed denominator (88% at 24 h, 80% at 48 h)",
      pp$pct_of_initial_parent[pp$timepoint == 24] == 88 &&
      pp$pct_of_initial_parent[pp$timepoint == 48] == 80)
check("a constant-rising short-mer rises on the same denominator (3% -> 5%)",
      mm$pct_of_initial_parent[mm$timepoint == 24] == 3 &&
      mm$pct_of_initial_parent[mm$timepoint == 48] == 5)

# A genuine t = 0 DOSING point (full parent signal) must stay the reference,
# not be mistaken for a predose control.
meta2 <- data.frame(sample = paste0("d", 1:6), group = "armB",
                    timepoint = rep(c(0, 2, 6), each = 2),
                    sample_type = "", stringsAsFactors = FALSE)
matches2 <- rbind(
  data.frame(met_id = "PARENT", met_name = "p", kind = "parent", sample = meta2$sample,
             area = c(10000,10000, 7000,7000, 5000,5000), stringsAsFactors = FALSE))
pr2 <- parent_relative_timecourse(matches2, meta2)
check("a full-signal t = 0 dosing point stays the reference (not excluded)",
      unique(pr2$reference_timepoint) == 0 &&
      pr2$pct_of_initial_parent[pr2$timepoint == 0] == 100)

# Caller-supplied reference timepoint is honored exactly.
pr3 <- parent_relative_timecourse(matches, meta, reference_timepoint = 24)
check("an explicit reference timepoint is honored (t = 24 -> parent 100% there)",
      unique(pr3$reference_timepoint) == 24 &&
      pr3$pct_of_initial_parent[pr3$met_id == "PARENT" & pr3$timepoint == 24] == 100)

# Internal standards never enter the species set.
check("internal standards are excluded", !any(pr$kind == "internal_standard"))

# Empty / no-timepoint inputs degrade gracefully.
check("no timepoint column -> empty frame, not an error",
      nrow(parent_relative_timecourse(matches, meta[, c("sample","group")])) == 0)

cat("\nAll parent-relative time-course checks passed.\n")
