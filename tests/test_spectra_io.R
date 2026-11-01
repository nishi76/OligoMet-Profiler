# test_spectra_io.R -- validate R/spectra_io.R: the Spectra/mzR reading
# layer, its automatic fallback to parse_mzml(), and the vendor-raw
# pre-conversion gate (resolve_ms_input_file()).
#
# read_ms_file()'s parity with parse_mzml() on a single fixture is already
# covered in tests/test_ms_matching.R; this file focuses on
# resolve_ms_input_file()'s dispatch logic and a broader round-trip across
# every bundled batch-example fixture.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R", "ms_matching.R", "spectra_io.R")) {
  source(file.path(.pkg_root, "R", .f))
}

cat("==== spectra_io.R validation ====\n\n")

## ---- resolve_ms_input_file(): pass-through for recognized formats --------
cat("--- resolve_ms_input_file() dispatch ---\n")
stopifnot(resolve_ms_input_file("/some/path/data.mzML") == "/some/path/data.mzML")
stopifnot(resolve_ms_input_file("/some/path/DATA.MZML") == "/some/path/DATA.MZML")
stopifnot(resolve_ms_input_file("/some/path/data.mzXML") == "/some/path/data.mzXML")
stopifnot(resolve_ms_input_file("/some/path/data.csv") == "/some/path/data.csv")
stopifnot(resolve_ms_input_file("/some/path/data.txt") == "/some/path/data.txt")
cat(".mzML/.mzXML/.csv/.txt pass through unchanged: PASS\n")

## ---- resolve_ms_input_file(): vendor formats need msconvert --------------
cat("\n--- resolve_ms_input_file() vendor formats ---\n")
have_msconvert <- !is.na(find_msconvert())
cat("msconvert available:", have_msconvert, "\n")
if (!have_msconvert) {
  err <- tryCatch({
    resolve_ms_input_file("/some/path/data.raw")
    NULL
  }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(err))
  stopifnot(grepl("msconvert", err, fixed = TRUE))
  cat("vendor .raw without msconvert fails with a clear error: PASS\n")
} else {
  cat("msconvert is available in this environment -- vendor-conversion",
      "success path is exercised implicitly by R/ms_matching.R's",
      "convert_vendor_raw(), not re-tested here (no real vendor file",
      "fixture is bundled).\n")
}

## ---- read_ms_file(): every bundled batch-example fixture -----------------
cat("\n--- read_ms_file() across all batch-example fixtures ---\n")
fixtures <- list.files(file.path(.pkg_root, "inst", "extdata", "batch_example"),
                        pattern = "\\.mzML$", full.names = TRUE)
stopifnot(length(fixtures) > 0)
cat("Fixtures:", paste(basename(fixtures), collapse = ", "), "\n")

for (fx in fixtures) {
  ref <- parse_mzml(fx)
  fb <- read_ms_file(fx, prefer_spectra = FALSE)
  stopifnot(identical(ref$ms1, fb$ms1))
  stopifnot(identical(ref$ms2, fb$ms2))
}
cat("Every fixture: read_ms_file(prefer_spectra=FALSE) == parse_mzml(): PASS\n")

if (.have_spectra()) {
  for (fx in fixtures) {
    ref <- parse_mzml(fx)
    sp <- read_ms_file(fx, prefer_spectra = TRUE)
    stopifnot(nrow(sp$ms1) == nrow(ref$ms1))
    stopifnot(nrow(sp$ms2) == nrow(ref$ms2))
    if (nrow(ref$ms1) > 0) {
      stopifnot(max(abs(sort(sp$ms1$mz) - sort(ref$ms1$mz))) < 1e-6)
    }
    if (nrow(ref$ms2) > 0) {
      stopifnot(all(sort(round(sp$ms2$mz, 4)) == sort(round(ref$ms2$mz, 4))))
    }
  }
  cat("Every fixture: read_ms_file(prefer_spectra=TRUE) agrees with parse_mzml(): PASS\n")
} else {
  cat("Spectra/mzR not installed -- Spectra-path round-trip SKIPPED\n",
      "(the fallback-path checks above still ran).\n")
}

cat("\n==== All spectra_io tests passed ====\n")
