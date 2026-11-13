# test_studio.R -- Studio tab data layer (R/studio.R): payload, detail, demo
.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R", "mass_isotope.R", "fragments.R",
             "ms_matching.R", "export_spectral.R", "studio.R")) source(file.path(.pkg_root, "R", .f))
check <- function(label, ok) {
  cat(sprintf("%-84s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Studio data layer ====\n\n")

## library-only payload (no batch results)
spec <- parse_input(GIVOSIRAN_SENSE_SPEC, dict = STANDARD_DICT, notation = "structured")
mets <- generate_metabolites(spec, list(max_3p = 3, max_5p = 3, endo = TRUE))
lib <- studio_payload(spec, mets, mode = "library")
check("library mode: no samples, one species per library entry", length(lib$samples) == 0 && length(lib$species) == length(mets))
check("oligo carries 21 bases/sugars and 20 linkage slots + terminal blank",
      length(lib$oligo$bases) == 21 && length(lib$oligo$sugars) == 21 && length(lib$oligo$linkages) == 21 && lib$oligo$linkages[21] == "")
check("parent species mass matches the engine", abs(lib$species[[1]]$mass - 8749.024) < 0.01)
sp <- Filter(function(s) s$id == "M06", lib$species)[[1]]
check("3' N-1 spans positions 1-20 and has lost the conjugate", sp$start == 1 && sp$end == 20 && sp$conj3 == "")
sp5 <- Filter(function(s) s$kind == "exo_5p", lib$species)[[1]]
check("5' N-1 spans positions 2-21 and keeps the L96 conjugate", sp5$start == 2 && sp5$end == 21 && grepl("L96", sp5$conj3))
j <- jsonlite::toJSON(lib, auto_unbox = TRUE, na = "null", digits = NA)
check("payload serialises to JSON", nchar(j) > 1000)

## demo inputs -> results payload
demo <- studio_demo_inputs()
pl <- studio_payload(demo$spec, demo$mets, demo$dict, demo$ms1_matches, demo$ms2_confirmations,
                     demo$sample_meta, mode = "demo")
check("demo: six samples ordered by timepoint", identical(vapply(pl$samples, function(s) s$time, numeric(1)), c(0, 4, 8, 24, 48, 72)))
ids <- vapply(pl$species, function(s) s$id, character(1))
check("demo: PS->PO variant is its own species row", "M01+ox1" %in% ids)
ox <- pl$species[[which(ids == "M01+ox1")]]
check("PS->PO variant mass is parent minus 15.9772 Da", abs(ox$mass - (8749.024 - 15.977157)) < 0.01)
tot <- Reduce(`+`, lapply(pl$species, function(s) ifelse(is.na(s$pct), 0, s$pct)))
check("percentages sum to 100 in every sample", all(abs(tot - 100) < 1e-6))
par_pct <- pl$species[[which(ids == "M01")]]$pct
check("parent share falls monotonically over time", all(diff(par_pct) < 0))
check("coverage is a fraction between 0 and 1", all(unlist(lapply(pl$species, function(s) s$cov)) >= 0 & unlist(lapply(pl$species, function(s) s$cov)) <= 1))

## detail
met1 <- demo$mets[[1]]
key <- names(demo$ms2_spectra)[grepl("24h.*\\|M01\\|0\\|", names(demo$ms2_spectra))][1]
det <- studio_detail(met1, demo$dict, demo$ms2_spectra[[key]], 0L, 9L)
check("detail: spectrum present, matched peaks carry fragment labels", det$has_spectrum && any(vapply(det$acquired, function(a) nzchar(a[[3]]), logical(1))))
check("detail: isotope pattern has 8 peaks, apex scaled to 1", length(det$iso) == 8 && max(vapply(det$iso, function(i) i[[3]], numeric(1))) == 1)
check("detail: monoisotopic peak is well below the apex for a 8.7 kDa strand", det$iso[[1]][[3]] < 0.5)
check("detail: charge envelope covers z = 3..12", length(det$env) == 10)
check("detail: score reports covered and total bond counts", det$score$total == 20 && det$score$covered > 0)
th <- studio_detail(met1, demo$dict, NULL, 0L)
check("detail without a spectrum returns the library only", !th$has_spectrum && length(th$theoretical) > 50 && length(th$acquired) == 0)
## data-matrix results: a minimal matrix (signal only, no ppm/rt/obs_mz/z)
## must not collapse the payload's per-sample vapply (regression).
ids_dm <- vapply(mets, function(m) m$id, "")
dm <- do.call(rbind, lapply(c("acid_0h", "acid_24h"), function(s) data.frame(
  met_id = ids_dm[1:3], met_name = "example", kind = c("parent", "exo_3p", "exo_3p"),
  sample = s, intensity = c(100, 30, 10), stringsAsFactors = FALSE)))
pl_dm <- tryCatch(studio_payload(spec, mets, ms1_matches = dm, mode = "results"),
                  error = function(e) conditionMessage(e))
check("data-matrix payload (signal-only schema) builds without error",
      is.list(pl_dm) && identical(pl_dm$mode, "results"))
check("data-matrix payload keeps per-sample signal; absent ppm/rt/z are NA", {
  sp <- pl_dm$species[[1]]
  "pct" %in% names(sp) && all(is.na(sp$ppm)) && all(is.na(sp$rt))
})

cat("\n==== All Studio data-layer tests passed ====\n")
