# test_quantification.R -- validate calibration-curve absolute quantification
# and baseline-relative quantification (R/statistics.R) against hand-built
# synthetic data with a known injected signal, same known-answer style as
# tests/test_statistics.R.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
source(file.path(.pkg_root, "R", "chemistry_dict.R"))  # %||%, .is_study_sample()
source(file.path(.pkg_root, "R", "degradation.R"))  # .best_signal_per_met()
source(file.path(.pkg_root, "R", "statistics.R"))

cat("==== Quantification suite validation ====\n\n")

## ---- fit_calibration_curve(): recovers a known linear relationship --------
cat("--- fit_calibration_curve() ---\n")
set.seed(42)
true_intercept <- 500
true_slope <- 2000
std_concs <- c(1, 5, 10, 50, 100, 500)
std_samples <- paste0("std_", seq_along(std_concs))

std_matches <- data.frame(
  met_id = "PARENT", met_name = "parent oligo", kind = "parent",
  sample = std_samples,
  area = true_intercept + true_slope * std_concs * (1 + rnorm(length(std_concs), 0, 0.02)),
  stringsAsFactors = FALSE
)
std_meta <- data.frame(sample = std_samples, group = "", timepoint = "",
                        sample_type = "standard", concentration = as.character(std_concs),
                        stringsAsFactors = FALSE)

curve <- fit_calibration_curve(std_matches, "PARENT", std_meta, weighting = "1/x2")
cat("Recovered slope =", round(curve$slope, 1), "(true", true_slope, ")\n")
cat("Recovered intercept =", round(curve$intercept, 1), "(true", true_intercept, ")\n")
cat("r_squared =", round(curve$r_squared, 4), "\n")
stopifnot(abs(curve$slope - true_slope) / true_slope < 0.05)
stopifnot(abs(curve$intercept - true_intercept) / true_intercept < 0.5)  # small n, intercept is noisier
stopifnot(curve$r_squared > 0.99)
stopifnot(curve$n_points == length(std_concs))
stopifnot(all(abs(curve$points$percent_re) < 15))  # standards back-calculate close to nominal
cat("Calibration curve recovers known slope/intercept with r_squared > 0.99: PASS\n")

## ---- fit_calibration_curve(): missing standards handled gracefully --------
cat("\n--- fit_calibration_curve() edge cases ---\n")
no_std_meta <- data.frame(sample = std_samples, group = "", timepoint = "",
                           sample_type = "unknown", concentration = "",
                           stringsAsFactors = FALSE)
empty_curve <- fit_calibration_curve(std_matches, "PARENT", no_std_meta, weighting = "1/x2")
stopifnot(is.null(empty_curve$model))
stopifnot(nzchar(empty_curve$note))
cat("No standard-type samples -> empty curve with a note, not an error: PASS\n")

too_few_meta <- std_meta[1:2, ]
few_curve <- fit_calibration_curve(std_matches[std_matches$sample %in% std_samples[1:2], ],
                                    "PARENT", too_few_meta, weighting = "1/x2", min_points = 3)
stopifnot(is.null(few_curve$model))
cat("Fewer than min_points standards -> empty curve, not a crash: PASS\n")

## ---- quantify_absolute(): back-calculates known unknown concentrations ----
cat("\n--- quantify_absolute() ---\n")
unk_concs <- c(unk_1 = 20, unk_2 = 80, unk_3 = 1000)  # unk_3 is above the standard range (extrapolated)
unk_samples <- names(unk_concs)
unk_matches <- data.frame(
  met_id = "PARENT", met_name = "parent oligo", kind = "parent",
  sample = unk_samples,
  area = true_intercept + true_slope * unk_concs,
  stringsAsFactors = FALSE
)
qc_matches <- data.frame(
  met_id = "PARENT", met_name = "parent oligo", kind = "parent",
  sample = "qc_1", area = true_intercept + true_slope * 30, stringsAsFactors = FALSE
)
all_matches <- do.call(rbind, list(std_matches, unk_matches, qc_matches))
all_meta <- rbind(
  std_meta,
  data.frame(sample = unk_samples, group = "control", timepoint = "",
             sample_type = "unknown", concentration = "", stringsAsFactors = FALSE),
  data.frame(sample = "qc_1", group = "", timepoint = "",
             sample_type = "quality_control", concentration = "30", stringsAsFactors = FALSE)
)

abs_res <- quantify_absolute(all_matches, "PARENT", all_meta, weighting = "1/x2")
quant <- abs_res$quant
cat("Rows returned:", nrow(quant), "(expected", length(std_samples) + length(unk_samples) + 1, ")\n")
stopifnot(nrow(quant) == length(std_samples) + length(unk_samples) + 1)

for (s in names(unk_concs)) {
  calc <- quant$concentration_calc[quant$sample == s]
  pct_err <- abs(calc - unk_concs[[s]]) / unk_concs[[s]] * 100
  cat(sprintf("  %s: nominal=%.0f calc=%.2f (%.1f%% error)\n", s, unk_concs[[s]], calc, pct_err))
  stopifnot(pct_err < 10)
}
cat("Unknown concentrations back-calculated within 10% of true value: PASS\n")

stopifnot(!quant$extrapolated[quant$sample == "unk_1"])
stopifnot(quant$extrapolated[quant$sample == "unk_3"])  # 1000 is above max standard (500)
cat("extrapolated flag correctly set only for the out-of-range unknown: PASS\n")

qc_row <- quant[quant$sample == "qc_1", ]
stopifnot(!is.na(qc_row$percent_re))
stopifnot(abs(qc_row$percent_re) < 15)
cat("QC sample gets a percent_re against its own nominal concentration: PASS\n")

unk_row <- quant[quant$sample == "unk_1", ]
stopifnot(is.na(unk_row$percent_re))  # unknowns have no nominal concentration to compare against
cat("Unknown (study) samples have no percent_re (no nominal to compare to): PASS\n")

## ---- quantify_relative(): time-series mode recovers a known trend --------
cat("\n--- quantify_relative(mode = 'time_series') ---\n")
ts_samples <- paste0("t", rep(c(0, 1, 2), each = 3), "_", rep(1:3, 3))
ts_timepoints <- rep(c(0, 1, 2), each = 3)
ts_fold <- c(1, 2, 4)[match(ts_timepoints, c(0, 1, 2))]  # doubles each timepoint
ts_matches <- data.frame(
  met_id = "DEG1", met_name = "degradant 1", kind = "exo_3p",
  sample = ts_samples, intensity = 1e5 * ts_fold, stringsAsFactors = FALSE
)
ts_meta <- data.frame(sample = ts_samples, group = "", timepoint = as.character(ts_timepoints),
                       sample_type = "unknown", concentration = "", stringsAsFactors = FALSE)

rel_ts <- quantify_relative(ts_matches, "DEG1", ts_meta, mode = "time_series", signal_col = "intensity")
t0 <- rel_ts$relative_signal[rel_ts$timepoint == "0"]
t2 <- rel_ts$relative_signal[rel_ts$timepoint == "2"]
cat("Relative signal at t0 (mean", round(mean(t0), 2), "), at t2 (mean", round(mean(t2), 2), ")\n")
stopifnot(abs(mean(t0) - 1) < 1e-6)
stopifnot(abs(mean(t2) - 4) < 1e-6)
cat("Baseline (time 0) normalizes to 1.0, later timepoint recovers 4x injected fold: PASS\n")

## ---- quantify_relative(): group mode relative to control ------------------
cat("\n--- quantify_relative(mode = 'group') ---\n")
grp_samples <- c("ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3")
grp_group <- c("control", "control", "control", "treated", "treated", "treated")
grp_fold <- ifelse(grp_group == "treated", 3, 1)
grp_matches <- data.frame(
  met_id = "DEG2", met_name = "degradant 2", kind = "exo_5p",
  sample = grp_samples, intensity = 1e4 * grp_fold, stringsAsFactors = FALSE
)
grp_meta <- data.frame(sample = grp_samples, group = grp_group, timepoint = "",
                        sample_type = "unknown", concentration = "", stringsAsFactors = FALSE)

rel_grp <- quantify_relative(grp_matches, "DEG2", grp_meta, mode = "group", control_group = "control",
                              signal_col = "intensity")
ctrl_rel <- rel_grp$relative_signal[rel_grp$group == "control"]
treat_rel <- rel_grp$relative_signal[rel_grp$group == "treated"]
cat("Control mean relative signal =", round(mean(ctrl_rel), 3),
    " Treated mean relative signal =", round(mean(treat_rel), 3), "\n")
stopifnot(abs(mean(ctrl_rel) - 1) < 1e-6)
stopifnot(abs(mean(treat_rel) - 3) < 1e-6)
cat("Control normalizes to 1.0, treated recovers 3x injected fold vs control: PASS\n")

## ---- quantify_relative(): calibration standards/QC never leak in ---------
cat("\n--- quantify_relative() excludes standards/QC/blanks ---\n")
mixed_meta <- rbind(
  grp_meta,
  data.frame(sample = "std_contaminant", group = "control", timepoint = "",
             sample_type = "standard", concentration = "100", stringsAsFactors = FALSE)
)
mixed_matches <- rbind(
  grp_matches,
  data.frame(met_id = "DEG2", met_name = "degradant 2", kind = "exo_5p",
             sample = "std_contaminant", intensity = 1e9, stringsAsFactors = FALSE)  # would badly skew the mean if included
)
rel_mixed <- quantify_relative(mixed_matches, "DEG2", mixed_meta, mode = "group", control_group = "control",
                                signal_col = "intensity")
stopifnot(!"std_contaminant" %in% rel_mixed$sample)
stopifnot(abs(mean(rel_mixed$relative_signal[rel_mixed$group == "control"]) - 1) < 1e-6)
cat("A standard-type row with a group value doesn't contaminate the control baseline: PASS\n")

## ---- quantify_metabolites(): splits absolute vs relative correctly -------
cat("\n--- quantify_metabolites() orchestrator ---\n")
# quantify_metabolites() picks ONE signal_col (via .auto_signal_col(), or
# the caller's own signal_col=) for the whole batch_matches table it's
# given -- realistic, since `area` is a property of how the batch/Python
# ROI pipeline ran (present for the whole run or not at all), never a
# per-metabolite choice. So this fixture uses `area` consistently across
# every metabolite, unlike the standalone quantify_relative() tests above
# which used `intensity` on their own (single-signal-column) fixtures.
combo_matches <- rbind(
  cbind(std_matches, group = "", timepoint = ""),
  cbind(unk_matches, group = "control", timepoint = ""),
  cbind(qc_matches, group = "", timepoint = ""),
  data.frame(met_id = grp_matches$met_id, met_name = grp_matches$met_name,
             kind = grp_matches$kind, sample = grp_matches$sample,
             area = grp_matches$intensity, group = grp_group, timepoint = "",
             stringsAsFactors = FALSE)
)
combo_meta <- rbind(all_meta, grp_meta)
combo_meta <- combo_meta[!duplicated(combo_meta$sample), ]

combo <- quantify_metabolites(combo_matches, combo_meta, absolute_met_ids = "PARENT",
                               mode = "group", control_group = "control", weighting = "1/x2")
stopifnot(all(combo$absolute$met_id == "PARENT"))
stopifnot(all(combo$relative$met_id == "DEG2"))
stopifnot("PARENT" %in% names(combo$calibration_curves))
cat("PARENT routed to absolute quant, DEG2 routed to relative quant: PASS\n")

cat("\n==== All quantification tests passed ====\n")
