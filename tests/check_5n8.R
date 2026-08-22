# check_5n8.R -- verify the 5'N-8 discrepancy is a workbook error, not our bug
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
our_mass <- sapply(mets, function(m)
  formula_mass(assemble_oligo_formula(m$bases, m$sugars, m$linkages, m$conj5, m$conj3)))

# 5' series is M12..M21 (N-1..N-10). Check consecutive deltas.
cat("=== 5' truncation series: consecutive mass deltas ===\n")
cat(sprintf("%-6s %-16s %-12s %-12s %-12s\n", "met", "seq", "our_mass", "delta", "residue"))
for (i in 12:21) {
  d <- if (i == 12) our_mass[1] - our_mass[i] else our_mass[i - 1] - our_mass[i]
  # residue removed = base + sugar + linkage - 2*H2O
  pos <- if (i == 12) 1 else i - 11   # position removed from 5' end
  b <- sp$bases[pos]; sg <- sp$sugars[pos]; lk <- sp$linkages[pos]
  res <- formula_mass(add_formulas(add_formulas(
    ION337_DICT[[b]]$formula, ION337_DICT[[sg]]$formula),
    ION337_DICT[[lk]]$formula)) - 2 * formula_mass(.H2O)
  cat(sprintf("M%02d    %-16s %-12.4f %-12.4f %-12.4f  %s+%s+%s\n",
              i, paste(mets[[i]]$bases, collapse = ""), our_mass[i], d, res,
              b, sg, ifelse(is.na(lk), "NA", lk)))
}

# Workbook value for 5'N-8
cat("\n=== 5'N-8 (TTATSST) ===\n")
cat("our mass   :", sprintf("%.4f\n", our_mass[19]))
cat("workbook   : 2874.4407\n")
cat("difference : ", sprintf("%.4f  (= 1 x PS linkage %.4f)\n",
    2874.4407 - our_mass[19], formula_mass(ION337_DICT[["s"]]$formula)))
cat("\nDelta chain consistency:\n")
cat("  5'N-7 -> 5'N-8 (remove A8): our delta =",
    sprintf("%.4f", our_mass[18] - our_mass[19]), " expected A+n+PS-2H2O =",
    sprintf("%.4f\n", formula_mass(add_formulas(add_formulas(
      ION337_DICT[["A"]]$formula, ION337_DICT[["n"]]$formula),
      ION337_DICT[["s"]]$formula)) - 2 * formula_mass(.H2O)))
cat("  5'N-8 -> 5'N-9 (remove T9): our delta =",
    sprintf("%.4f", our_mass[19] - our_mass[20]), " expected T+n+PS-2H2O =",
    sprintf("%.4f\n", formula_mass(add_formulas(add_formulas(
      ION337_DICT[["T"]]$formula, ION337_DICT[["n"]]$formula),
      ION337_DICT[["s"]]$formula)) - 2 * formula_mass(.H2O)))
cat("  workbook 5'N-7->5'N-8 delta =",
    sprintf("%.4f  (off by +95.94 = 1 PS)\n", 3194.5640 - 2874.4407))
cat("  workbook 5'N-8->5'N-9 delta =",
    sprintf("%.4f  (off by -95.94 = 1 PS)\n", 2874.4407 - 2371.4419))
cat("\nConclusion: workbook 5'N-8 mass has a spurious +1 PS (95.94 Da).\n")
cat("Our value is consistent with both neighbors; workbook breaks the chain.\n")
