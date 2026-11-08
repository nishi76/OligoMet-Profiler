# test_reference_oligos.R -- bundled reference oligos and report helpers
# Checks:
#   1. Givosiran sense strand uses the L96 triantennary GalNAc (P1, N11), not
#      the generic GalNAc3 entry (which is ~972.64 Da lighter).
#   2. Isotope-plot bar width stays well below the isotope spacing.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R", "mass_isotope.R")) {
  source(file.path(.pkg_root, "R", .f))
}
# .isotope_bar_width is a pure helper in build_report.R; extract it so the
# test does not need ggplot2.
.src <- parse(file.path(.pkg_root, "R", "build_report.R"))
for (.e in .src) if (is.call(.e) && identical(as.character(.e[[2]]), ".isotope_bar_width")) eval(.e)

check <- function(label, ok) {
  cat(sprintf("%-80s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Reference oligo suite ====\n\n")

giv <- REFERENCE_OLIGOS$givosiran
check("registry conj3 is GalNAc3_triantennary", giv$conj3 == "GalNAc3_triantennary")
check("structured spec conj3 is GalNAc3_triantennary", GIVOSIRAN_SENSE_SPEC$conj3 == "GalNAc3_triantennary")

mass_formula <- function(conj3) {
  sp <- GIVOSIRAN_SENSE_SPEC; sp$conj3 <- conj3
  parse_input(sp, dict = STANDARD_DICT, notation = "structured")
}
spec_new <- mass_formula("GalNAc3_triantennary")
spec_old <- mass_formula("GalNAc3")
mk <- function(sp) metabolite_mass_info(generate_metabolites(sp, list(max_3p = 0, max_5p = 0, endo = FALSE))[[1]])
inew <- mk(spec_new); iold <- mk(spec_old)
# 21 nt -> 20 linkage P + 1 P in the L96 linker; 79 nucleobase N + 11 in the linker
check("givosiran sense strand has P = 21", inew$formula_vec[["P"]] == 21)
check("givosiran sense strand has N = 90", inew$formula_vec[["N"]] == 90)
check("triantennary is 972.639 Da heavier than generic GalNAc3",
      abs((inew$mono_mass - iold$mono_mass) - 972.6388) < 0.01)

for (z in c(3, 6, 12)) {
  mz <- 900
  check(sprintf("isotope bar width < isotope spacing at z=%d", z),
        .isotope_bar_width(mz, z) < 1.00335 / z)
}
cat("\n==== All reference oligo tests passed ====\n")
