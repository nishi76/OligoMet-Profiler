# test_library_figures.R -- predicted-library characterization figures
# (R/library_figures.R): the mass-frame builder (base R) and, when ggplot2 is
# available, that each figure builds. No analyte name is ever hardcoded.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R", "library_figures.R"))
  source(file.path(.pkg_root, "R", .f))
check <- function(label, ok) {
  cat(sprintf("%-72s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Library-figures suite ====\n\n")

spec <- list(bases = strsplit("GCATGACTCC", "")[[1]],
             sugars = c(rep("e", 3), rep("d", 4), rep("e", 3)),
             linkages = c(rep("s", 9), NA), conj5 = "none", conj3 = "GalNAc3", n = 10)
mets <- generate_metabolites(spec, list(oligo_name = "OLM"))
mf <- library_mass_frame(mets)
prm <- prm_inclusion_list(mets, z_range = 5:10, max_oxid = 3)

## ---- mass frame -------------------------------------------------------------
check("mass frame has one row per metabolite", nrow(mf) == length(mets))
check("mass frame columns are present", all(c("id", "name", "kind", "n", "n_ps", "n_po",
      "formula", "mono_mass", "avg_mass", "mass_shift") %in% names(mf)))
check("parent mass shift is exactly 0", mf$mass_shift[mf$kind == "parent"][1] == 0)
check("truncation/fragment metabolites are lighter than parent", {
  frag <- mf$mass_shift[mf$kind %in% c("exo_3p", "exo_5p", "endo_5frag_p", "endo_3frag_p")]
  all(frag < 0)
})
check("every metabolite has a finite mass and a formula",
      all(is.finite(mf$mono_mass)) && all(nzchar(mf$formula)))
check("empty library yields an empty frame, not an error", nrow(library_mass_frame(list())) == 0)

## ---- Kendrick math ----------------------------------------------------------
check("PS->PO step and Kendrick base constants are the documented values",
      abs(.LIBFIG_PSPO_DA - 15.977156) < 1e-9 &&
      abs(.LIBFIG_KENDRICK_EXACT - 62.945833) < 1e-9 && .LIBFIG_KENDRICK_NOMINAL == 63)

## ---- figures build (needs ggplot2) ------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  is_gg <- function(p) inherits(p, "ggplot")
  NAME <- "MyOligo-7"   # a user-supplied name; must appear, never a hardcoded one
  check("mass-shift waterfall builds with the supplied name", {
    p <- plot_mass_shift_waterfall(mf, NAME)
    is_gg(p) && grepl(NAME, p$labels$title, fixed = TRUE)
  })
  check("mass-defect plot builds", is_gg(plot_mass_defect(mf, NAME)))
  check("Kendrick mass-defect plot builds", is_gg(plot_kendrick_defect(prm, mf, NAME)))
  check("charge-state heatmap builds", is_gg(plot_charge_heatmap(prm, NAME)))
  check("PRM m/z distribution builds", is_gg(plot_prm_mz_distribution(prm, NAME)))
  check("metabolite-class count builds", is_gg(plot_metabolite_categories(mf, NAME)))
  check("empty inputs give a placeholder plot, not an error",
        is_gg(plot_mass_defect(data.frame(), NAME)) && is_gg(plot_charge_heatmap(data.frame(), NAME)))
} else {
  cat("  (skipping figure-build checks: ggplot2 not installed)\n")
}

cat("\nAll library-figures checks passed.\n")
