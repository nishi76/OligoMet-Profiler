# tests/testthat.R -- entry point for the asserted unit-test suite.
#
# Run with:
#   Rscript tests/testthat.R
#
# Sources the pipeline modules straight from R/, the same approach the
# console-check scripts in tests/ use, so the suite exercises the current
# working tree rather than whatever version happens to be installed.
# Exits with a non-zero status if any test fails or errors, so it can gate
# a CI job (see .github/workflows/tests.yml) -- unlike the console-check
# scripts, which only stop() on the specific invariants each one bothered
# to check inline.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})

for (.f in c("about.R", "default_params.R", "progress_utils.R", "chemistry_dict.R",
             "oligo_io.R", "metabolites.R", "mass_isotope.R", "fragments.R",
             "ms_matching.R")) {
  source(file.path(.pkg_root, "R", .f))
}

library(testthat)

result <- test_dir(file.path(.pkg_root, "tests", "testthat"),
                    reporter = "summary", stop_on_failure = FALSE)
res_df <- as.data.frame(result)
n_bad <- sum(res_df$failed > 0 | res_df$error)
if (n_bad > 0) {
  cat("\n", n_bad, "test file(s) had failures or errors.\n")
  quit(status = 1, save = "no")
}
