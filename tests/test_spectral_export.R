# test_spectral_export.R -- validate the MGF/MSP spectral library writers
# Checks:
#   1. MS1 spectra carry the isotope cluster and the correct precursor m/z
#   2. MS1 peak annotations are true isotopologue offsets (M+0 present)
#   3. MS2 spectra carry fragment ions at the right precursor m/z
#   4. MGF block structure: matched BEGIN/END IONS, parseable peak lines
#   5. MSP block structure: "Num Peaks:" agrees with the peak lines that follow
#   6. The same library written to both formats has the same peaks

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (m in c("about", "progress_utils", "chemistry_dict", "oligo_io", "metabolites",
            "mass_isotope", "fragments", "export_spectral")) {
  source(file.path(.pkg_root, "R", paste0(m, ".R")))
}

sp <- parse_input(NUSINERSEN_TRIPLET)
mets <- generate_metabolites(sp, list(oligo_name = "nusinersen", max_3p = 2,
                                      max_5p = 2, endo = FALSE))
info <- metabolite_mass_info(mets[[1]])
cat("=== 1. Library construction ===\n")
cat("metabolites:", length(mets), "  parent mono:",
    sprintf("%.4f Da\n", info$mono_mass))

z_ms1 <- 4:6
ms1 <- build_ms1_library(mets, z_range = z_ms1, n_iso = 5,
                         oligo_name = "nusinersen")
ms2 <- build_ms2_library(mets, precursor_z_range = 5:6, frag_z_range = 1:2,
                         oligo_name = "nusinersen")
cat("MS1 spectra:", length(ms1), " (expected", length(mets) * length(z_ms1),
    ")\n")
cat("MS2 spectra:", length(ms2), " (expected", length(mets) * 2, ")\n")
if (length(ms1) != length(mets) * length(z_ms1))
  stop("MS1 library should hold one spectrum per metabolite per charge state")
if (length(ms2) != length(mets) * 2)
  stop("MS2 library should hold one spectrum per metabolite per precursor z")

cat("\n=== 2. MS1 precursor m/z matches the charge envelope ===\n")
worst <- 0
for (r in ms1) {
  expect <- charge_envelope(r$mono_mass, r$z)$mz[1]
  worst <- max(worst, abs(r$precursor_mz - expect))
}
cat("max |precursor_mz - charge_envelope|:", sprintf("%.3e Da\n", worst))
if (worst > 1e-6)
  stop("MS1 precursor m/z does not agree with charge_envelope()")

# The monoisotopic peak must be in the cluster and labelled M+0, otherwise
# the library's precursor and its peak list disagree about what M is.
lab0 <- vapply(ms1, function(r) "M+0" %in% r$peaks$annotation, logical(1))
cat("spectra containing an M+0 peak:", sum(lab0), "/", length(ms1), "\n")
if (!all(lab0)) stop("Some MS1 spectra have no monoisotopic (M+0) peak")

cat("\n=== 3. MS2 precursor m/z and fragment count ===\n")
r <- ms2[[1]]
expect <- charge_envelope(r$mono_mass, r$z)$mz[1]
cat("first MS2 spectrum:", r$name, "\n")
cat("  precursor:", sprintf("%.4f", r$precursor_mz),
    " expected:", sprintf("%.4f", expect), "\n")
cat("  peaks:", nrow(r$peaks), " first annotations:",
    paste(head(r$peaks$annotation, 4), collapse = ", "), "\n")
if (abs(r$precursor_mz - expect) > 1e-3)
  stop("MS2 precursor m/z does not agree with charge_envelope()")
if (nrow(r$peaks) < 10) stop("MS2 spectrum has implausibly few fragment ions")
if (is.unsorted(r$peaks$mz)) stop("MS2 peak list is not sorted by m/z")

cat("\n=== 4. MGF structure ===\n")
out <- file.path(tempdir(), "spectral_test")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
paths <- export_spectral_libraries(mets, out_dir = out, prefix = "nusinersen",
                                   z_range = z_ms1, n_iso = 5,
                                   precursor_z_range = 5:6,
                                   frag_z_range = 1:2,
                                   oligo_name = "nusinersen")
for (nm in c("ms1_mgf", "ms2_mgf")) {
  ln <- readLines(paths[[nm]])
  nb <- sum(ln == "BEGIN IONS"); ne <- sum(ln == "END IONS")
  npep <- sum(grepl("^PEPMASS=", ln)); nchg <- sum(grepl("^CHARGE=", ln))
  cat(sprintf("  %-8s %d BEGIN / %d END / %d PEPMASS / %d CHARGE\n",
              nm, nb, ne, npep, nchg))
  if (nb != ne || nb != npep || nb != nchg)
    stop("MGF block structure is inconsistent in ", nm)
  # Every non-header line inside a block must parse as "m/z intensity".
  peak_lines <- ln[grepl("^[0-9]", ln)]
  bad <- peak_lines[!grepl("^[0-9.]+ [0-9.]+$", peak_lines)]
  if (length(bad) > 0)
    stop("Malformed MGF peak line in ", nm, ": ", bad[1])
  cat("    peak lines:", length(peak_lines), "all well-formed\n")
}

cat("\n=== 5. MSP structure (Num Peaks agrees with the block) ===\n")
for (nm in c("ms1_msp", "ms2_msp")) {
  ln <- readLines(paths[[nm]])
  starts <- grep("^NAME: ", ln)
  declared <- as.integer(sub("^Num Peaks: ", "",
                             ln[grep("^Num Peaks: ", ln)]))
  npk <- grep("^Num Peaks: ", ln)
  actual <- integer(length(npk))
  for (i in seq_along(npk)) {
    j <- npk[i] + 1L
    cnt <- 0L
    while (j <= length(ln) && nzchar(ln[j])) { cnt <- cnt + 1L; j <- j + 1L }
    actual[i] <- cnt
  }
  cat(sprintf("  %-8s %d spectra, Num Peaks matches block: %s\n",
              nm, length(starts), identical(declared, actual)))
  if (length(starts) != length(declared))
    stop("MSP has a NAME without a Num Peaks in ", nm)
  if (!identical(declared, actual))
    stop("MSP 'Num Peaks' disagrees with the peak lines in ", nm)
}

cat("\n=== 6. MGF and MSP carry the same peaks ===\n")
mgf_pk <- readLines(paths[["ms1_mgf"]])
mgf_pk <- as.numeric(sub(" .*$", "", mgf_pk[grepl("^[0-9]", mgf_pk)]))
msp_ln <- readLines(paths[["ms1_msp"]])
msp_pk <- as.numeric(sub("\t.*$", "", msp_ln[grepl("^[0-9]", msp_ln)]))
cat("  MGF peaks:", length(mgf_pk), " MSP peaks:", length(msp_pk),
    " max |diff|:", sprintf("%.1e\n", max(abs(sort(mgf_pk) - sort(msp_pk)))))
if (length(mgf_pk) != length(msp_pk) ||
    max(abs(sort(mgf_pk) - sort(msp_pk))) > 1e-6)
  stop("MGF and MSP exports of the same library disagree")

cat("\n=== 7. Files written ===\n")
for (p in paths) cat(sprintf("  %-34s %7d bytes\n", basename(p),
                             file.info(p)$size))

cat("\n=== 8. Empirical MS2 library: consensus denoising across replicates ===\n")
source(file.path(.pkg_root, "R", "batch_ms_processing.R"))

met <- mets[[1]]
frags <- generate_fragments(met, ion_types = c("aB", "w", "y", "b"), z_range = 1)
true_mz <- vapply(frags[1:3], function(f) f$mz_table$mz[1], numeric(1))

# Three "replicate" samples: each carries the same three real fragment
# peaks (small ppm jitter, as real scans would show) plus one noise peak
# unique to that sample -- a peak seen by only one of three spectra should
# not survive min_consensus_fraction = 0.5.
set.seed(7)
ms2_spectra <- list()
ms2_confirmations <- data.frame(
  sample = c("s1", "s2", "s3"), met_id = met$id, met_name = met$name,
  k_oxid = 0, z = 1, adduct = "H", stringsAsFactors = FALSE
)
for (i in 1:3) {
  jitter <- 1 + rnorm(3, 0, 2e-6)          # ~2 ppm jitter, well inside fragment_ppm=10
  noise_mz <- true_mz[1] * (1 + runif(1, 500, 900) / 1e6)  # unique per sample, well outside 10 ppm
  spec <- data.frame(mz = c(true_mz * jitter, noise_mz),
                     intensity = c(100, 80, 60, 40), stringsAsFactors = FALSE)
  ms2_spectra[[paste0("s", i, "|", met$id, "|0|1|H")]] <- spec
}

emp <- build_empirical_ms2_library(mets, STANDARD_DICT, ms2_confirmations, ms2_spectra,
                                    frag_z_range = 1, fragment_ppm = 10,
                                    min_consensus_fraction = 0.5,
                                    min_peaks_for_library = 2)
cat("groups written:", length(emp$records), " (expect 1)\n")
if (length(emp$records) != 1) stop("expected exactly one empirical MS2 library entry")
r <- emp$records[[1]]
cat("consensus peaks:", nrow(r$peaks), " (expect 3 -- the noise peaks must not recur)\n")
if (nrow(r$peaks) != 3) stop("consensus should keep exactly the 3 recurring real peaks")
if (emp$summary$n_source_spectra[1] != 3) stop("summary should report 3 contributing spectra")
if (emp$summary$n_consensus_peaks[1] != 3) stop("summary should report 3 consensus peaks")
cat("annotated peaks:", sum(nzchar(r$peaks$annotation)), "/", nrow(r$peaks),
    "(theoretical-fragment labels are cosmetic only)\n")

cat("\n--- A metabolite confirmed in only one sample keeps all its peaks ---\n")
met2 <- mets[[2]]
ms2_confirmations2 <- data.frame(sample = "s1", met_id = met2$id, met_name = met2$name,
                                  k_oxid = 0, z = 1, adduct = "H", stringsAsFactors = FALSE)
ms2_spectra2 <- list()
ms2_spectra2[[paste0("s1|", met2$id, "|0|1|H")]] <- data.frame(
  mz = c(300.1234, 450.5678), intensity = c(100, 50), stringsAsFactors = FALSE)
emp2 <- build_empirical_ms2_library(mets, STANDARD_DICT, ms2_confirmations2, ms2_spectra2,
                                     frag_z_range = 1, fragment_ppm = 10,
                                     min_peaks_for_library = 2)
cat("single-sample consensus peaks:", nrow(emp2$records[[1]]$peaks), "(expect 2, nothing to vote against)\n")
if (nrow(emp2$records[[1]]$peaks) != 2) stop("a lone contributing spectrum should keep all its peaks")

cat("\n--- Written MSP is well-formed ---\n")
emp_msp <- file.path(tempdir(), "empirical_ms2_test.msp")
write_msp(emp$records, emp_msp, measured = TRUE)
ln <- readLines(emp_msp)
if (!any(grepl("^Num Peaks: 3$", ln))) stop("empirical MSP should declare 3 peaks")
if (!any(grepl("^COMMENT: Acquired MS2 spectrum", ln))) stop("empirical MSP should be marked as acquired data")
cat("empirical MSP written and well-formed:", emp_msp, "\n")

cat("\n==== All spectral export tests passed ====\n")
