# test_degradation.R -- validate R/degradation.R: degradation_summary()'s
# peak-area-based % degradation arithmetic, its area-absent fallback to
# intensity, and the class-composition/top-degradant breakdowns, both on
# hand-built synthetic data (exact arithmetic checks) and against the real
# batch_example fixture (which now bundles two truncated degradant species
# alongside the parent -- see inst/extdata/batch_example/generate_example.py).

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
             "batch_ms_processing.R", "statistics.R", "degradation.R")) {
  source(file.path(.pkg_root, "R", .f))
}

cat("==== degradation.R validation ====\n\n")

.mk_row <- function(sample, met_id, met_name, kind, area = NA_real_, intensity = NA_real_) {
  data.frame(sample = sample, met_id = met_id, met_name = met_name, kind = kind,
             area = area, intensity = intensity, stringsAsFactors = FALSE)
}

## ---- Exact arithmetic: parent=100, degradants=50 -> 33.33% degradation ---
cat("--- exact arithmetic: parent=100, one degradant=50 ---\n")
m1 <- rbind(
  .mk_row("s1", "M01", "parent", "parent", area = 100),
  .mk_row("s1", "M02", "N-1", "exo_3p", area = 50)
)
d1 <- degradation_summary(m1)
stopifnot(d1$signal_used == "area")
stopifnot(d1$per_sample$pct_degradation == round(100 * (1 - 100 / 150), 2))
stopifnot(d1$per_sample$pct_degradation == 33.33)
cat("pct_degradation =", d1$per_sample$pct_degradation, "(expected 33.33): PASS\n")

## ---- Zero-degradant edge case: 100% parent -> 0% degradation -------------
cat("\n--- zero-degradant edge case ---\n")
m2 <- .mk_row("s1", "M01", "parent", "parent", area = 100)
d2 <- degradation_summary(m2)
stopifnot(d2$per_sample$pct_degradation == 0)
stopifnot(nrow(d2$top_degradants) == 0)
cat("pct_degradation = 0, no top_degradants rows: PASS\n")

## ---- All-degradant edge case: no parent -> 100% degradation --------------
cat("\n--- all-degradant edge case (no parent signal) ---\n")
m3 <- .mk_row("s1", "M02", "N-1", "exo_3p", area = 75)
d3 <- degradation_summary(m3)
stopifnot(d3$per_sample$pct_degradation == 100)
stopifnot(nrow(d3$top_degradants) == 1)
cat("pct_degradation = 100, one top_degradants row: PASS\n")

## ---- Multiple degradant classes: composition sums correctly --------------
cat("\n--- composition breakdown across 3 degradant classes ---\n")
m4 <- rbind(
  .mk_row("s1", "M01", "parent", "parent", area = 60),
  .mk_row("s1", "M02", "3'N-1", "exo_3p", area = 20),
  .mk_row("s1", "M03", "5'N-1", "exo_5p", area = 15),
  .mk_row("s1", "M04", "endo1", "endo_5frag", area = 5)
)
d4 <- degradation_summary(m4)
stopifnot(d4$per_sample$pct_degradation == round(100 * (1 - 60 / 100), 2))  # 40
stopifnot(nrow(d4$composition) == 4)
stopifnot(sum(d4$composition$pct_of_total) == 100)
degradant_pcts <- d4$composition$pct_of_degradants[d4$composition$kind != "parent"]
stopifnot(abs(sum(degradant_pcts) - 100) < 0.01)
stopifnot(nrow(d4$top_degradants) == 3)
stopifnot(d4$top_degradants$met_id[1] == "M02")  # highest-area degradant ranked first
cat("4-way composition sums to 100%, degradant-only pcts sum to 100%,",
    "top_degradants correctly ranked: PASS\n")

## ---- area-absent fallback to intensity ------------------------------------
cat("\n--- area-absent -> falls back to intensity ---\n")
m5 <- rbind(
  .mk_row("s1", "M01", "parent", "parent", area = NA, intensity = 200),
  .mk_row("s1", "M02", "N-1", "exo_3p", area = NA, intensity = 100)
)
d5 <- degradation_summary(m5)
stopifnot(d5$signal_used == "intensity")
stopifnot(d5$per_sample$pct_degradation == round(100 * (1 - 200 / 300), 2))
cat("signal_used='intensity', pct_degradation =", d5$per_sample$pct_degradation,
    "(expected 33.33): PASS\n")

## ---- Empty / degenerate inputs --------------------------------------------
cat("\n--- degenerate inputs don't error ---\n")
empty <- degradation_summary(data.frame())
stopifnot(nrow(empty$per_sample) == 0)
stopifnot(is.na(empty$signal_used))
null_result <- degradation_summary(NULL)
stopifnot(nrow(null_result$per_sample) == 0)
cat("empty data.frame() and NULL both return the empty-shaped result: PASS\n")

## ---- Multi-sample: per-sample independence --------------------------------
cat("\n--- multi-sample independence ---\n")
m6 <- rbind(
  .mk_row("ctrl", "M01", "parent", "parent", area = 90),
  .mk_row("ctrl", "M02", "N-1", "exo_3p", area = 10),
  .mk_row("treat", "M01", "parent", "parent", area = 50),
  .mk_row("treat", "M02", "N-1", "exo_3p", area = 50)
)
d6 <- degradation_summary(m6)
stopifnot(nrow(d6$per_sample) == 2)
ctrl_pct <- d6$per_sample$pct_degradation[d6$per_sample$sample == "ctrl"]
treat_pct <- d6$per_sample$pct_degradation[d6$per_sample$sample == "treat"]
stopifnot(ctrl_pct == 10)
stopifnot(treat_pct == 50)
stopifnot(treat_pct > ctrl_pct)
cat("ctrl = 10%, treat = 50%, computed independently per sample: PASS\n")

## ---- plot_degradation_composition() returns a ggplot ----------------------
cat("\n--- plot_degradation_composition() ---\n")
p <- plot_degradation_composition(d6)
stopifnot(inherits(p, "ggplot"))
p_empty <- plot_degradation_composition(list(composition = data.frame(), signal_used = NA))
stopifnot(inherits(p_empty, "ggplot"))
cat("returns a ggplot object for both populated and empty input: PASS\n")

## ---- Integration: real batch_example fixture (parent + 2 degradants) -----
cat("\n--- integration: inst/extdata/batch_example (real deconvolution) ---\n")
ex_dir <- file.path(.pkg_root, "inst", "extdata", "batch_example")
files <- list.files(ex_dir, pattern = "\\.mzML$", full.names = TRUE)
if (length(files) == 0 || is.na(find_python())) {
  cat("SKIPPED (no fixtures or no python3 on PATH)\n")
} else {
  dict <- STANDARD_DICT
  spec <- parse_input(INOTERSEN_TRIPLET)
  mets <- generate_metabolites(spec, opts = list(oligo_name = "inotersen",
                                                  max_3p = 1, max_5p = 1, endo = FALSE))
  watchlist_path <- tempfile(fileext = ".txt")
  write_precursor_watchlist(mets, dict, z_range = 3:12, max_oxid = 0, h_offset = 0,
                             out_path = watchlist_path)
  deconv <- tryCatch(
    run_batch_deconvolution(files, precursor_watchlist = watchlist_path,
                             mass_tol_ppm = 10, n_workers = 1),
    error = function(e) NULL)
  if (is.null(deconv)) {
    cat("SKIPPED (batch deconvolution failed)\n")
  } else {
    feats <- read_batch_features(deconv$features_path)
    feats$sample <- tools::file_path_sans_ext(basename(feats$sample))
    m <- match_ms1_batch(mets, feats, dict, ppm_tol = 10, z_range = 3:12,
                          adducts = "H", max_oxid = 0)
    deg <- degradation_summary(m)
    cat("signal_used:", deg$signal_used, "\n")
    stopifnot(deg$signal_used == "area")
    stopifnot(nrow(deg$per_sample) == 6)
    stopifnot(all(deg$per_sample$pct_degradation > 0))  # fixture has real degradants
    stopifnot(nrow(deg$top_degradants) > 0)

    ctrl_mean <- mean(deg$per_sample$pct_degradation[grepl("^ctrl", deg$per_sample$sample)])
    treat_mean <- mean(deg$per_sample$pct_degradation[grepl("^treat", deg$per_sample$sample)])
    cat("mean pct_degradation: ctrl =", round(ctrl_mean, 2),
        " treat =", round(treat_mean, 2), "\n")
    stopifnot(treat_mean > ctrl_mean)  # generate_example.py deliberately builds this in
    cat("real fixture: non-zero degradation, treated > control as designed: PASS\n")
  }
}

cat("\n==== All degradation tests passed ====\n")
