# test_blank_correction.R -- validate reagent/matrix-blank detection and
# background subtraction (R/blank_correction.R, plus the RB/MB name-token
# detection helpers in R/chemistry_dict.R) against hand-built synthetic
# data with a known injected blank offset, the same known-answer style as
# tests/test_statistics.R.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
source(file.path(.pkg_root, "R", "chemistry_dict.R"))
source(file.path(.pkg_root, "R", "degradation.R"))
source(file.path(.pkg_root, "R", "statistics.R"))
source(file.path(.pkg_root, "R", "blank_correction.R"))

cat("==== Blank correction suite validation ====\n\n")

## ---- .detect_blank_type_from_name() ------------------------------------------
cat("--- .detect_blank_type_from_name() ---\n")
stopifnot(identical(.detect_blank_type_from_name("RB_01"), "reagent_blank"))
stopifnot(identical(.detect_blank_type_from_name("Sample-MB-2"), "matrix_blank"))
stopifnot(identical(.detect_blank_type_from_name("RB"), "reagent_blank"))
stopifnot(identical(.detect_blank_type_from_name("rb_1"), "reagent_blank"))  # case-insensitive
stopifnot(identical(.detect_blank_type_from_name("sample.MB.3"), "matrix_blank"))  # dot delimiter
stopifnot(is.na(.detect_blank_type_from_name("CARBON")))       # bare substring, not a token
stopifnot(is.na(.detect_blank_type_from_name("MBio_01")))      # bare substring, not a token
stopifnot(is.na(.detect_blank_type_from_name("Sample_1")))     # no RB/MB at all
cat("All name-detection cases: PASS\n")

## ---- .is_blank_sample() -------------------------------------------------------
cat("\n--- .is_blank_sample() ---\n")
stopifnot(.is_blank_sample("reagent_blank", "Sample_1"))          # explicit sample_type, name has no token
stopifnot(.is_blank_sample("unknown", "RB_01"))                   # name token alone, sample_type still "unknown"
stopifnot(!.is_blank_sample("unknown", "CARBON_dosed"))           # real sample, name happens to contain "CARBON"
stopifnot(!.is_blank_sample(NA_character_, "Sample_1"))
cat("All .is_blank_sample() cases: PASS\n")

## ---- Synthetic data with a known injected blank offset -----------------------
cat("\n--- compute_blank_signal() / apply_blank_correction() ---\n")

# M1: two blank samples (one reagent_blank by sample_type, one matrix_blank
# by RB/MB name token) with intensities 100 and 200 -> pooled mean 150.
# Study samples S1 (raw 1000, well above blank) and S2 (raw 50, BELOW the
# blank mean -- exercises the floor-at-zero case).
# M2: no blank observations at all -> blank_mean should be 0, corrected == raw.
batch_matches <- data.frame(
  met_id      = c("M1", "M1", "M1", "M1", "M2", "M2"),
  met_name    = c("met1", "met1", "met1", "met1", "met2", "met2"),
  kind        = "parent",
  sample      = c("RB_1", "MatrixBlank_A", "S1", "S2", "S1", "S2"),
  intensity   = c(100, 200, 1000, 50, 500, 700),
  stringsAsFactors = FALSE
)
sample_meta <- data.frame(
  sample      = c("RB_1", "MatrixBlank_A", "S1", "S2"),
  sample_type = c("unknown", "matrix_blank", "unknown", "unknown"),  # RB_1 detected by NAME, not sample_type
  stringsAsFactors = FALSE
)

blank <- compute_blank_signal(batch_matches, sample_meta)
stopifnot(nzchar(blank$note) == FALSE)
m1_row <- blank$table[blank$table$met_id == "M1", ]
m2_row <- blank$table[blank$table$met_id == "M2", ]
stopifnot(m1_row$n_blank == 2)
stopifnot(isTRUE(all.equal(m1_row$blank_mean, mean(c(100, 200)))))
stopifnot(isTRUE(all.equal(m1_row$blank_sd, stats::sd(c(100, 200)))))
stopifnot(m2_row$n_blank == 0)
stopifnot(m2_row$blank_mean == 0)
cat("compute_blank_signal() recovers the injected pooled blank mean/SD: PASS\n")

corrected <- apply_blank_correction(batch_matches, sample_meta)
stopifnot("intensity_bcorr" %in% names(corrected))
get_val <- function(met, samp) corrected$intensity_bcorr[corrected$met_id == met & corrected$sample == samp]
stopifnot(isTRUE(all.equal(get_val("M1", "S1"), max(0, 1000 - 150))))
stopifnot(isTRUE(all.equal(get_val("M1", "S2"), 0)))  # 50 - 150 < 0 -> floored at 0
stopifnot(isTRUE(all.equal(get_val("M2", "S1"), 500)))  # no blank -> unchanged
stopifnot(isTRUE(all.equal(get_val("M2", "S2"), 700)))
# Raw column must be untouched.
stopifnot(all(corrected$intensity == batch_matches$intensity))
cat("apply_blank_correction() matches max(0, raw - blank_mean) exactly, floors correctly: PASS\n")

## ---- No blanks present: empty note, no _bcorr columns added ------------------
cat("\n--- No blanks present ---\n")
no_blank_matches <- data.frame(
  met_id = c("M1", "M1"), met_name = c("met1", "met1"), kind = "parent",
  sample = c("S1", "S2"), intensity = c(1000, 700), stringsAsFactors = FALSE
)
no_blank_meta <- data.frame(sample = c("S1", "S2"), sample_type = c("unknown", "unknown"),
                             stringsAsFactors = FALSE)
blank_empty <- compute_blank_signal(no_blank_matches, no_blank_meta)
stopifnot(nzchar(blank_empty$note))
stopifnot(all(blank_empty$table$n_blank == 0))
no_blank_corrected <- apply_blank_correction(no_blank_matches, no_blank_meta)
stopifnot(!"intensity_bcorr" %in% names(no_blank_corrected))
cat("No blanks -> empty-with-note, no _bcorr column added: PASS\n")

## ---- build_abundance_matrix(signal_col = "area") regression guard ------------
cat("\n--- build_abundance_matrix(signal_col=) matches the old area wrapper ---\n")
area_matches <- data.frame(
  met_id = c("M1", "M1", "M2", "M2"), met_name = c("met1", "met1", "met2", "met2"),
  kind = "parent", sample = c("S1", "S2", "S1", "S2"),
  intensity = c(10, 20, 30, 40), area = c(100, 200, 300, 400),
  stringsAsFactors = FALSE
)
via_signal_col <- build_abundance_matrix(area_matches, signal_col = "area")
via_wrapper <- build_abundance_matrix_area(area_matches)
stopifnot(identical(via_signal_col, via_wrapper))
# Default (no signal_col given) must still behave exactly as before.
via_default <- build_abundance_matrix(area_matches)
stopifnot(all(via_default$S1 == c(10, 30)))
cat("build_abundance_matrix()/build_abundance_matrix_area() regression: PASS\n")

cat("\n==== All blank correction tests passed ====\n")
