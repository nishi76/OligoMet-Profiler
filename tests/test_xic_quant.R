# test_xic_quant.R -- targeted summed-XIC quantitation (R/xic_quant.R).
# Pure-R checks always run (isotope clusters, charge parsing, target
# building, ion selection + summation on a hand-built ion table). The
# end-to-end check on inst/extdata/xic_example/ runs only when Python with
# inst/python/requirements.txt is available.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "progress_utils.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R", "ms_matching.R", "batch_ms_processing.R",
             "degradation.R", "statistics.R", "xic_quant.R")) {
  source(file.path(.pkg_root, "R", .f))
}

check <- function(label, ok) {
  cat(sprintf("%-80s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Targeted summed-XIC quantitation ====\n\n")

## ---- isotope clusters ------------------------------------------------------
sp <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(sp, opts = list(oligo_name = "ino", max_3p = 3, max_5p = 3,
                                             endo = FALSE, min_frag_len = 3))
info <- metabolite_mass_info(mets[[1]], STANDARD_DICT)
cl <- isotope_clusters(info$formula_vec, info$mono_mass, n = 12, use_envipat = FALSE)
check("clusters ranked by abundance, rank 1 first", identical(cl$iso_rank, seq_len(nrow(cl))) &&
        all(diff(cl$abundance) <= 0))
check("7.2 kDa PS oligo: most abundant cluster is M+3..M+5, not monoisotopic",
      cl$cluster[1] %in% 3:5 && cl$cluster[1] != 0)
check("cluster masses step by ~1.003 Da",
      all(abs(diff(sort(cl$mass)) - 1.0031) < 0.01))
check("top 5 clusters carry 65-75% of the pattern", sum(cl$abundance[1:5]) > 0.65 &&
        sum(cl$abundance[1:5]) < 0.75)

## ---- charge spec -------------------------------------------------------------
check("parse_charge_spec('3-8') -> 3:8", identical(parse_charge_spec("3-8"), 3:8))
check("parse_charge_spec('4, 6,5') -> 4:6", identical(parse_charge_spec("4, 6,5"), 4:6))
check("blank charge spec -> NULL (automatic ranking)", is.null(parse_charge_spec(" ")))

## ---- targets -----------------------------------------------------------------
mm <- data.frame(met_id = c("M01", "M01", "M01"), met_name = "ino", kind = "parent",
                 k_oxid = c(0, 0, 1), adduct = "H", z = c(6, 7, 7), rt = c(4.82, 4.84, 4.9),
                 intensity = c(1e6, 2e6, 1e4), area = c(1e5, 2e5, 1e3), sample = "s1",
                 stringsAsFactors = FALSE)
tg <- build_xic_targets(mets, mm, z_range = 3:12, n_iso_candidates = 8, use_envipat = FALSE)
check("one target per identified metabolite", identical(unique(tg$target_id), "M01"))
check("dominant variant chosen (k_oxid 0, not the weak oxidized hit)", all(tg$k_oxid == 0))
check("expected RT = median RT of that variant", abs(tg$rt_expected[1] - 4.83) < 1e-9)
check("8 isotopes x 10 charges inside the m/z range", nrow(tg) == 80)
z7 <- tg[tg$z == 7 & tg$iso_rank == 1, ]
check("z=7 m/z of the top cluster = (cluster mass - 7 x proton)/7",
      abs(z7$mz - (cl$mass[1] - 7 * .PROTON) / 7) < 1e-6)
check("monoisotopic peak (~2%) is not among the 8 candidates", !0 %in% tg$cluster)

## ---- selection + summation on a hand-built ion table ------------------------
# 3 charges x 3 isotopes, known areas; charge 5 is strongest in standards,
# charge 3 only matters in the unknown.
mk <- function(sample, z, iso, area) data.frame(
  target_id = "M01", ion_id = paste0("z", z, "_i", iso), sample = sample, area = area,
  height = area * 10, rt_apex = 4.8, rt_start = 4.7, rt_end = 4.9, sn_all_ions = 100,
  met_id = "M01", met_name = "ino", kind = "parent", k_oxid = 0, adduct = "H",
  rt_expected = 4.8, z = z, cluster = iso + 2, iso_rank = iso, rel_abundance = 0.1,
  stringsAsFactors = FALSE)
w_z <- c(`3` = 1, `4` = 3, `5` = 5)
ions <- do.call(rbind, c(
  lapply(1:3, function(i) do.call(rbind, lapply(3:5, function(z) mk("std", z, i, 100 * w_z[[as.character(z)]] / i)))),
  lapply(1:3, function(i) do.call(rbind, lapply(3:5, function(z) mk("unk", z, i, 100 * c(`3` = 6, `4` = 3, `5` = 1)[[as.character(z)]] / i))))))
meta <- data.frame(sample = c("std", "unk"), sample_type = c("standard", "unknown"),
                   stringsAsFactors = FALSE)
s <- summarize_xic_quant(ions, meta, n_charges = 2, n_isotopes = 2)
check("charges ranked on the standard (z5, z4), not the unknown",
      s$selection$charges_used == "4,5")
check("isotopes ranked by theoretical rank (ranks 1-2)", s$selection$isotopes_used == "M+3,M+4")
exp_std <- sum(100 * c(3, 5) %o% (1 / 1:2))
check("summed area = sum of the chosen ions", abs(s$quant$area[s$quant$sample == "std"] - exp_std) < 1e-6)
check("same ion list applied to every sample",
      s$quant$n_ions[s$quant$sample == "unk"] == 4 && s$quant$charges_used[2] == "4,5")
s2 <- summarize_xic_quant(ions, meta, charges = 3:5, n_isotopes = 3)
check("fixed charge list overrides the ranking", s2$selection$charges_used == "3,4,5")
check("all candidates chosen -> 100% of candidate signal", abs(s2$selection$pct_candidate_signal - 100) < 1e-9)
check("contribution map flags selected ions",
      sum(s$contributions$selected) == 4 && abs(sum(s$contributions$pct_of_candidate_signal) - 100) < 1e-9)
ions_nf <- ions; ions_nf$rt_apex[ions_nf$sample == "unk"] <- NA
check("no peak found -> sample left out (not detected), not zero",
      !"unk" %in% summarize_xic_quant(ions_nf, meta)$quant$sample)
invisible(plot_xic_contributions(s$contributions, "M01"))
invisible(plot_xic_contributions(data.frame(met_id = character(0)), "M01"))

## ---- end to end on the bundled example (needs Python) ------------------------
py <- find_python()
py_ok <- !is.na(py) && identical(suppressWarnings(system2(py, c("-c", shQuote("import pyteomics, scipy, pandas")),
                                                          stdout = FALSE, stderr = FALSE)), 0L)
if (!py_ok) {
  cat("\n(skipping end-to-end XIC check: Python with inst/python/requirements.txt not available)\n")
} else {
  cat("\n--- end to end: inst/extdata/xic_example ---\n")
  ex <- file.path(.pkg_root, "inst", "extdata", "xic_example")
  files <- list.files(ex, "\\.mzML$", full.names = TRUE)
  meta <- utils::read.csv(file.path(ex, "sample_meta.csv"), colClasses = "character")
  truth <- utils::read.csv(file.path(ex, "unknown_truth.csv"))
  mod <- file.path(.pkg_root, "inst", "python")
  dc <- run_batch_deconvolution(files, output_dir = tempfile(), z_range = 3:12,
                                min_intensity = 5000, n_workers = 2, module_dir = mod)
  res <- annotate_metabolites_batch(mets, read_batch_features(dc$features_path), z_range = 3:12,
                                    adducts = "H", max_oxid = 0, n_iso = 0, use_envipat = FALSE)
  tg <- build_xic_targets(mets, res$ms1_matches, z_range = 3:12, use_envipat = FALSE)
  x <- run_xic_quantitation(files, tg, n_workers = 2, module_dir = mod)
  check("every file gets a parent XIC area, incl. 1 ng/mL standards identification missed",
        all(c("std_L1_r1", "std_L1_r2") %in% x$ions$sample[x$ions$area > 0]) &&
          !"std_L1_r1" %in% res$ms1_matches$sample)
  bc <- function(n_z, n_i, charges = NULL) {
    s <- summarize_xic_quant(x$ions, meta, n_charges = n_z, n_isotopes = n_i, charges = charges)
    q <- quantify_metabolites(s$quant, meta, absolute_met_ids = "M01", mode = "group",
                              weighting = "1/x2")
    list(s = s, q = q, u = stats::setNames(q$absolute$concentration_calc, q$absolute$sample))
  }
  r55 <- bc(5, 5)
  check("5x5: all 10 standards on the curve, r^2 > 0.999",
        r55$q$calibration_curves$M01$n_points == 10 && r55$q$calibration_curves$M01$r_squared > 0.999)
  qc <- r55$q$absolute[r55$q$absolute$sample_type == "quality_control", ]
  check("5x5: QCs within +/-15%", all(abs(qc$percent_re) < 15))
  check("5x5: normal-CSD unknowns within +/-5%",
        all(abs(r55$u[c("U1", "U2")] / truth$true_concentration[1:2] - 1) < 0.05))
  r38 <- bc(0, 5, charges = 3:8)
  check("charges 3-8: shifted-CSD unknowns within +/-5% (broad set is CSD-robust)",
        all(abs(r38$u[c("U3", "U4")] / truth$true_concentration[3:4] - 1) < 0.05))
  check("auto top-5 charges under-reads the shifted-CSD unknowns (z=3 left out)",
        all(r55$u[c("U3", "U4")] / truth$true_concentration[3:4] < 0.92))
}
cat("\nAll XIC quantitation checks passed.\n")
