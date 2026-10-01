# test_internal_standard.R -- one analyte + internal standards
# (R/internal_standard.R): FASTA record roles, IS library entries, mass
# conflict checks, IS normalization (global + per-metabolite), IS response
# monitoring, and IS handling in data matrices, degradation, and class sums.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R", "mass_isotope.R",
             "fragments.R", "ms_matching.R", "batch_ms_processing.R", "degradation.R",
             "statistics.R", "blank_correction.R", "internal_standard.R", "data_matrix.R")) {
  source(file.path(.pkg_root, "R", .f))
}
check <- function(label, ok) {
  cat(sprintf("%-80s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Internal standard suite ====\n\n")

## ---- role detection ------------------------------------------------------------
check("IS / SIL / ISTD / IS2 / 'internal standard' tokens are internal standards",
      all(detect_is_role(c("ION582-IS", "ION582_SIL", "istd", "IS2", "Internal Standard A", "gapmer IS"))))
check("'ISIS 420915', 'analysis', 'This' stay analytes",
      !any(detect_is_role(c("ISIS 420915", "analysis", "This", "nusinersen"))))

## ---- sequence sets -------------------------------------------------------------
ino <- INOTERSEN_TRIPLET
gap <- "Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm"
one <- parse_sequence_set(ino, default_name = "my_oligo")
check("plain sequence (no '>') = one analyte, named by default_name",
      nrow(one$table) == 1 && one$table$role == "analyte" && one$analyte$name == "my_oligo" &&
        length(one$problems) == 0)
two <- parse_sequence_set(paste0(">inotersen 5'-none\n", ino, "\n>gapmer IS\n", gap))
check("two records: analyte + IS by name token",
      identical(two$table$role, c("analyte", "internal_standard")) &&
        two$table$role_source[2] == "IS token in name" && length(two$internal_standards) == 1)
check("each record keeps its own sequence (no concatenation)",
      two$analyte$spec$n == 20 && two$internal_standards[[1]]$spec$n == 16)
key <- parse_sequence_set(paste0(">ISIS 420915\n", ino, "\n>calibrant role=IS\n", gap))
check("role=IS header key makes an IS without any token; ISIS stays analyte",
      identical(key$table$role, c("analyte", "internal_standard")) &&
        key$table$role_source[2] == "header role= key" && key$table$name[2] == "calibrant")
dup <- parse_sequence_set(paste0(">A\n", ino, "\n>B\n", gap))
check("two analytes -> a problem asking to mark one role=IS", any(grepl("one analyte per run", dup$problems)))
bad <- parse_sequence_set(paste0(">A\n", ino, "\n>B IS\nXx-sYy"))
check("unparseable IS record -> reported problem, not a crash", any(grepl("record 2", bad$problems)))
conj <- parse_sequence_set(paste0(">A 3'-GalNAc3\n", ino))
check("header conjugate tag is applied", conj$analyte$spec$conj3 == "GalNAc3")

## ---- IS library entries and mass conflicts --------------------------------------
mets <- generate_metabolites(two$analyte$spec, opts = list(oligo_name = "inotersen", max_3p = 3,
                                                         max_5p = 3, endo = FALSE))
ism <- build_is_metabolites(two$internal_standards)
check("IS entry: id IS01, kind internal_standard, full length, no truncations",
      ism[[1]]$id == "IS01" && ism[[1]]$kind == "internal_standard" && ism[[1]]$n == 16 &&
        length(ism) == 1)
check("IS mass from its own sequence", abs(metabolite_mass_info(ism[[1]], STANDARD_DICT)$mono_mass -
                                             two$table$mono_mass[2]) < 1e-6)
check("unrelated analog IS -> no conflicts", nrow(check_is_interference(mets, ism)) == 0)
n1_is <- build_is_metabolites(list(list(name = "n-1 IS", spec = truncate_3p(two$analyte$spec, 1))))
cf <- check_is_interference(mets, n1_is)
check("IS identical to the 3' N-1 metabolite -> same-mass error",
      nrow(cf) >= 1 && any(cf$severity == "error" & cf$met_id == "M02"))

## ---- normalization -----------------------------------------------------------------
mk <- function(sample, met, kind, area) data.frame(met_id = met, met_name = met, kind = kind,
  sample = sample, intensity = area * 10, area = area, stringsAsFactors = FALSE)
m <- rbind(mk("s1", "M01", "parent", 1000), mk("s1", "M02", "exo_3p", 100),
           mk("s1", "IS01", "internal_standard", 500), mk("s1", "IS02", "internal_standard", 250),
           mk("s2", "M01", "parent", 2000), mk("s2", "M02", "exo_3p", 300),
           mk("s2", "IS01", "internal_standard", 1000), mk("s2", "IS02", "internal_standard", 100),
           mk("s3", "M01", "parent", 800))
n <- apply_is_normalization(m, "IS01")
check("global IS: area_isr = area / IS01 area in the same sample",
      isTRUE(all.equal(n$area_isr[n$met_id == "M01"], c(2, 2, NA))))
check("sample without IS -> NA ratio, listed as missing",
      is.na(n$area_isr[n$sample == "s3"]) && "s3" %in% attr(n, "is_missing_samples"))
check("IS rows themselves get NA", all(is.na(n$area_isr[n$kind == "internal_standard"])))
ov <- apply_is_normalization(m, "IS01", overrides = data.frame(met_id = "M02", is_id = "IS02"))
check("per-metabolite override: M02 / IS02, M01 still / IS01",
      isTRUE(all.equal(ov$area_isr[ov$met_id == "M02"], c(100 / 250, 300 / 100))) &&
        isTRUE(all.equal(ov$area_isr[ov$met_id == "M01" & ov$sample != "s3"], c(2, 2))))
check("assignment table records global vs per-metabolite",
      identical(attr(ov, "is_assignment")$assignment, c("global", "per-metabolite")))
mb <- m; mb$area_bcorr <- mb$area - 50
nb <- apply_is_normalization(mb, "IS01")
check("blank-corrected ratio = corrected analyte / RAW IS",
      isTRUE(all.equal(nb$area_bcorr_isr[nb$met_id == "M01" & nb$sample == "s1"], 950 / 500)))

## ---- IS response monitoring -----------------------------------------------------------
meta <- data.frame(sample = c("c1", "c2", "q1", "b1", "u1", "u2", "u3"),
                   sample_type = c("standard", "standard", "quality_control", "reagent_blank",
                                   "unknown", "unknown", "unknown"), stringsAsFactors = FALSE)
isd <- do.call(rbind, Map(function(s, a) mk(s, "IS01", "internal_standard", a),
                          c("c1", "c2", "q1", "b1", "u1", "u2"), c(100, 110, 90, 5, 40, 160)))
r <- is_response_summary(isd, meta, "IS01")
check("reference = mean of standards + QCs (100)", abs(attr(r, "reference_mean") - 100) < 1e-9)
check("statuses: low, high, not detected, blank not evaluated",
      identical(r$status[match(c("u1", "u2", "u3", "b1"), r$sample)],
                c("low", "high", "IS not detected", "blank (not evaluated)")))
check("window is editable", r$status[r$sample == "u1"] == "low" &&
        is_response_summary(isd, meta, "IS01", window = c(30, 170))$status[5] == "ok")
invisible(plot_is_response(r)); invisible(plot_is_response(data.frame()))

## ---- downstream exclusions + data matrix ------------------------------------------------
deg <- degradation_summary(m, top_n = 5)
check("degradation ignores the IS (parent/degradant sums)",
      !any(deg$composition$kind == "internal_standard") &&
        abs(deg$per_sample$total_signal[deg$per_sample$sample == "s1"] - 1100) < 1e-9)
kb <- build_kind_abundance_matrix(m, signal_col = "area")
check("class sums have no internal_standard class", !"internal_standard" %in% kb$met_id)
dm <- data.frame(met_id = c("M01", "N1", "IS-1", "calib"), kind = c("parent", "exo_3p", NA, "ISTD"),
                 s1 = c(1000, 100, 500, 400), check.names = FALSE, stringsAsFactors = FALSE)
am <- annotate_matrix_metabolites(read_data_matrix(dm)$matches)
check("matrix: IS token in id and kind = ISTD both become internal_standard",
      all(am$kind[am$met_id %in% c("IS-1", "calib")] == "internal_standard"))
check("matrix: parent guess never picks an IS row", guess_parent_met_id(am) == "M01")

cat("\nAll internal standard checks passed.\n")
