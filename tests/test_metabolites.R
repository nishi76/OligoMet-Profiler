# test_metabolites.R -- validate metabolites.R against ION337
.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
source(file.path(.pkg_root, "R", "chemistry_dict.R"))
source(file.path(.pkg_root, "R", "oligo_io.R"))
source(file.path(.pkg_root, "R", "metabolites.R"))

sp <- parse_input(ION337_TRIPLET)
mets <- generate_metabolites(sp, list(oligo_name = "ION337", endo = FALSE))
tab <- metabolite_table(mets)

cat("=== Metabolite count ===\n")
cat("total:", length(mets), " (expect 21: parent + 10x3p + 10x5p)\n")
cat("kinds:", paste(names(table(tab$kind)), table(tab$kind), collapse = "  "), "\n")

cat("\n=== Parent mass check ===\n")
p <- mets[[1]]
f <- assemble_oligo_formula(p$bases, p$sugars, p$linkages, p$conj5, p$conj3)
cat("parent formula:", format_formula(f), "  mass:", sprintf("%.4f", formula_mass(f)), "\n")
cat("n_ps:", p$n_ps, " n_po:", p$n_po, "\n")

cat("\n=== 3' N-1 (remove last T+n+PS) ===\n")
m3 <- mets[[2]]
f3 <- assemble_oligo_formula(m3$bases, m3$sugars, m3$linkages, m3$conj5, m3$conj3)
cat("bases:", paste(m3$bases, collapse = ""), "  n:", m3$n, "  n_ps:", m3$n_ps, "\n")
cat("formula:", format_formula(f3), "  mass:", sprintf("%.4f", formula_mass(f3)), "\n")
cat("mass delta from parent:", sprintf("%.4f", formula_mass(f3) - formula_mass(f)), "\n")
cat("expected delta (T+n+PS-2H2O = C13H18N3O8PS):",
    sprintf("%.4f", formula_mass(parse_formula("C13H18N3O8PS"))), "\n")

cat("\n=== 5' N-1 (remove first G+e+PS) ===\n")
m5 <- mets[[12]]
f5 <- assemble_oligo_formula(m5$bases, m5$sugars, m5$linkages, m5$conj5, m5$conj3)
cat("bases:", paste(m5$bases, collapse = ""), "  n:", m5$n, "  n_ps:", m5$n_ps, "\n")
cat("formula:", format_formula(f5), "  mass:", sprintf("%.4f", formula_mass(f5)), "\n")
cat("mass delta from parent:", sprintf("%.4f", formula_mass(f5) - formula_mass(f)), "\n")

cat("\n=== Full table ===\n")
print(tab[, c("id", "name", "n", "n_ps", "bases")], right = FALSE)
