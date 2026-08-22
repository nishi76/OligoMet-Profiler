# validate_vs_workbook.R -- compare generated metabolite masses to ION337 workbook
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

library(openxlsx)
REFERENCE_WORKBOOK <- Sys.getenv("ION337_WORKBOOK_PATH", file.path(.pkg_root, "tests", "ION337_Workbook.xlsx"))
if (!file.exists(REFERENCE_WORKBOOK)) {
  stop("Reference workbook not found at ", REFERENCE_WORKBOOK, ". ",
       "Place ION337_Workbook.xlsx in tests/, or set the ",
       "ION337_WORKBOOK_PATH environment variable to its location. ",
       "This file is not distributed with the package.")
}
wb <- readWorkbook(REFERENCE_WORKBOOK, sheet = 1)

cat("=== Exact column names ===\n")
print(names(wb))

# Find columns by partial match (robust to bracket/space conversion)
nm <- names(wb)
col_name   <- nm[grep("^Name", nm)][1]
col_seq    <- nm[grep("^Seq", nm) & grepl("ence", nm)][1]
col_mass   <- nm[grep("Mono", nm)][1]
col_len    <- nm[grep("Length", nm)][1]
col_mod    <- nm[grep("Modif", nm)][1]
cat("mapped:", col_name, "|", col_seq, "|", col_mass, "|", col_len, "|", col_mod, "\n")

# Keep only metabolite header rows (Name + Seq.Length both present)
has_name <- !is.na(wb[[col_name]])
has_len  <- !is.na(wb[[col_len]])
hdr <- wb[has_name & has_len, c(col_name, col_seq, col_mass, col_len, col_mod)]
names(hdr) <- c("Name", "Sequence", "MonoMass", "SeqLen", "Mod")
hdr$MonoMass <- as.numeric(hdr$MonoMass)
# Drop the ION1884600 variant (different compound, same sequence)
hdr <- hdr[!grepl("1884600", hdr$Name), ]

cat("\n=== Workbook metabolite rows (ION337 only) ===\n")
print(hdr, right = FALSE)

# Build our library
sp <- parse_input(ION337_TRIPLET)
mets <- generate_metabolites(sp, list(oligo_name = "ION337", endo = FALSE))
our_mass <- sapply(mets, function(m) {
  formula_mass(assemble_oligo_formula(m$bases, m$sugars, m$linkages, m$conj5, m$conj3))
})

# Match by sequence string
cat("\n=== Mass comparison (by sequence) ===\n")
cat(sprintf("%-20s %-16s %-12s %-12s %-10s\n", "wb_name", "sequence", "wb_mass", "our_mass", "ppm"))
n_match <- 0; max_ppm <- 0
for (i in seq_len(nrow(hdr))) {
  seq_wb <- gsub(" ", "", hdr$Sequence[i])
  idx <- which(sapply(mets, function(m) paste(m$bases, collapse = "") == seq_wb))
  if (length(idx) == 1) {
    ppm <- abs(our_mass[idx] - hdr$MonoMass[i]) / hdr$MonoMass[i] * 1e6
    n_match <- n_match + 1; max_ppm <- max(max_ppm, ppm)
    cat(sprintf("%-20s %-16s %-12.4f %-12.4f %-10.4f %s\n",
                hdr$Name[i], seq_wb, hdr$MonoMass[i], our_mass[idx], ppm,
                ifelse(ppm < 1, "OK", "CHECK")))
  } else {
    cat(sprintf("%-20s %-16s %-12.4f %-12s\n", hdr$Name[i], seq_wb, hdr$MonoMass[i], "NO MATCH"))
  }
}
cat(sprintf("\nMatched %d/%d  |  max ppm = %.4f\n", n_match, nrow(hdr), max_ppm))
