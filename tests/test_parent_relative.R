# test_parent_relative.R -- abundance as a percentage of the full-length
# parent, and selection at or above a threshold (R/statistics.R).

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "statistics.R")) {
  source(file.path(.pkg_root, "R", .f))
}
check <- function(label, ok) {
  cat(sprintf("%-78s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Parent-relative abundance suite ====\n\n")

bm <- data.frame(
  met_id   = c("M01","M01","M02","M02","M03","M03","M09","M09"),
  met_name = c("parent","parent","N-1","N-1","endo","endo","IS","IS"),
  kind     = c("parent","parent","exo_3p","exo_3p","endo_5frag","endo_5frag",
               "internal_standard","internal_standard"),
  sample   = c("s1","s2","s1","s2","s1","s2","s1","s2"),
  area     = c(1000,1200, 100,90, 0.5,0.4, 800,820),
  intensity= c(1e6,1.1e6, 1e5,9e4, 500,400, 8e5,8.2e5),
  stringsAsFactors = FALSE)

rel <- parent_relative_abundance(bm, "area")
check("parent is 100% of itself", round(rel$max_pct[rel$met_id == "M01"]) == 100)
check("N-1 reads ~10% of parent", abs(rel$max_pct[rel$met_id == "M02"] - 10) < 0.5)
check("endo fragment reads ~0.05% of parent", rel$max_pct[rel$met_id == "M03"] < 0.1)

keep <- select_mets_by_parent_pct(bm, 0.1, "area")
check("0.1% floor keeps parent, N-1 and IS", setequal(keep, c("M01","M02","M09")))
check("0.1% floor drops the 0.05% endo fragment", !("M03" %in% keep))

# A lower floor keeps the endo fragment too.
keep2 <- select_mets_by_parent_pct(bm, 0.01, "area")
check("0.01% floor keeps the endo fragment", "M03" %in% keep2)

# IS is kept even when below the floor (it is a spiked reference).
bm_lowIS <- bm; bm_lowIS$area[bm_lowIS$met_id == "M09"] <- c(0.3, 0.3)
check("internal standard is always kept regardless of threshold",
      "M09" %in% select_mets_by_parent_pct(bm_lowIS, 0.1, "area"))

# No parent -> nothing to normalize against, keep everything.
bm_np <- bm[bm$kind != "parent", ]
check("no parent row -> all metabolites kept (nothing to filter on)",
      setequal(select_mets_by_parent_pct(bm_np, 0.1, "area"), unique(bm_np$met_id)))

# Falls back to intensity when area is absent.
bm_noarea <- bm; bm_noarea$area <- NULL
check("falls back to intensity when area column is absent",
      nrow(parent_relative_abundance(bm_noarea, "area")) == 4)

cat("\nAll parent-relative abundance checks passed.\n")
