# =============================================================================
# batch_ms_processing.R
# R-side glue for the Python parallel charge-envelope deconvolution pipeline
# (inst/python/oligomet_deconv/) -- invoked as a subprocess via system2(),
# the same idiom already used for the ProteoWizard msconvert bridge in
# R/ms_matching.R. Reuses match_ms1()/envelope_consistency()/
# find_ms2_spectra()/confirm_metabolite()/prm_inclusion_list() unchanged;
# this file adds no new mass-spec math of its own, only per-sample looping
# and the R<->Python data handoff.
# =============================================================================

## ---- Locating Python -------------------------------------------------------
find_python <- function() {
  for (bin in c("python3", "python")) {
    path <- Sys.which(bin)
    if (nzchar(path)) return(bin)
  }
  NA_character_
}

# Resolve inst/python/ the same way ms_matching.R resolves py-adjacent
# assets: a source checkout first (so edits take effect without
# reinstalling), then the installed package, then a couple of interactive
# working-directory fallbacks. (The candidate order previously contradicted
# this comment -- the installed package was checked first, so an installed
# copy silently shadowed any local edits to inst/python/ until the package
# was reinstalled.)
.find_deconv_module_dir <- function() {
  candidates <- c(
    file.path(getwd(), "inst", "python"),
    file.path(getwd(), "..", "inst", "python"),
    tryCatch(system.file("python", package = "OligoMetProfiler"), error = function(e) "")
  )
  hit <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (length(hit) > 0) hit[1] else candidates[1]
}

## ---- Precursor watch-list (for targeted MS2 capture) ----------------------
# Thin wrapper around the EXISTING prm_inclusion_list() -- no new mass math.
# The watch-list only needs to be a broad net (Python's --ms2-watch-ppm
# defaults to 50 ppm), so using the H-adduct-only PRM list is sufficient:
# DDA fragmentation is triggered off the dominant charge/adduct envelope in
# practice, and the final confirmation match in R is exact regardless.
write_precursor_watchlist <- function(mets, dict = STANDARD_DICT, z_range = 3:12,
                                       max_oxid = 6, h_offset = 0, out_path) {
  prm <- prm_inclusion_list(mets, dict, z_range = z_range, h_offset = h_offset,
                             max_oxid = max_oxid)
  writeLines(as.character(unique(prm$precursor_mz)), out_path)
  out_path
}

## ---- Invoking the batch deconvolution CLI ----------------------------------
run_batch_deconvolution <- function(files, output_dir = tempdir(),
                                     output_file = "combined_features.tsv",
                                     ms2_output_file = NULL,
                                     precursor_watchlist = NULL, ms2_watch_ppm = 50,
                                     roi_ppm = 15, rt_tol = 0.15, mass_tol_ppm = 20,
                                     z_range = 3:20, min_intensity = 1e4, min_scans = 3,
                                     max_gap_scans = 2, min_charge_states = 2, n_workers = NULL,
                                     python_bin = find_python(),
                                     module_dir = .find_deconv_module_dir(),
                                     progress = NULL, console_tracker = NULL) {
  if (is.na(python_bin)) {
    stop("No python3/python interpreter found on PATH. Batch MS processing ",
         "requires Python 3.9+ with the packages in inst/python/requirements.txt installed.")
  }
  if (!dir.exists(module_dir)) {
    stop("Could not locate the oligomet_deconv Python module (looked in: ", module_dir, ").")
  }
  if (!is.null(precursor_watchlist) && is.null(ms2_output_file)) {
    ms2_output_file <- "combined_ms2.tsv"
  }

  args <- c(
    "-m", "oligomet_deconv.cli",
    "--input", files,
    "--output-dir", output_dir,
    "--output-file", output_file,
    "--roi-ppm", roi_ppm, "--rt-tol", rt_tol, "--mass-tol-ppm", mass_tol_ppm,
    "--z-min", min(z_range), "--z-max", max(z_range),
    "--min-intensity", min_intensity, "--min-scans", min_scans,
    "--max-gap-scans", max_gap_scans, "--min-charge-states", min_charge_states
  )
  if (!is.null(n_workers)) args <- c(args, "--n-workers", n_workers)
  if (!is.null(precursor_watchlist)) {
    args <- c(args, "--precursor-watchlist", precursor_watchlist,
              "--ms2-watch-ppm", ms2_watch_ppm, "--ms2-output-file", ms2_output_file)
  }

  if (!is.null(progress)) progress("Running parallel deconvolution (Python)...")
  # `python3 -m oligomet_deconv.cli` needs inst/python/ on PYTHONPATH to
  # find the module. An earlier version used setwd(module_dir) for this,
  # which silently resolved any RELATIVE output_dir/output_file against
  # inst/python/ instead of the caller's directory -- and, worse, setwd()
  # mutates the whole R process's working directory, which is unsafe if
  # this is ever called from a Shiny process serving concurrent sessions.
  # Prepending to PYTHONPATH instead leaves the working directory alone.
  old_pythonpath <- Sys.getenv("PYTHONPATH", unset = NA)
  new_pythonpath <- if (is.na(old_pythonpath) || !nzchar(old_pythonpath)) {
    module_dir
  } else {
    paste(module_dir, old_pythonpath, sep = .Platform$path.sep)
  }
  .reset_pythonpath <- function() {
    if (is.na(old_pythonpath)) Sys.unsetenv("PYTHONPATH") else Sys.setenv(PYTHONPATH = old_pythonpath)
  }
  on.exit(.reset_pythonpath(), add = TRUE)
  Sys.setenv(PYTHONPATH = new_pythonpath)

  status <- system2(python_bin, args = as.character(args), stdout = TRUE, stderr = TRUE)
  exit_code <- attr(status, "status")
  if (!is.null(console_tracker)) console_tracker(paste(status, collapse = "\n"))
  if (!is.null(exit_code) && exit_code != 0) {
    stop("Batch deconvolution failed:\n", paste(status, collapse = "\n"))
  }

  # ROI/charge-envelope detection is designed for centroided peaks; a
  # profile-mode input file still runs (just slowly, treating every raw
  # sample point as a candidate peak), so the Python side flags it via
  # this sidecar (present only when at least one file triggered it).
  profile_warn_path <- file.path(output_dir, "_profile_mode_warnings.tsv")
  profile_mode_files <- if (file.exists(profile_warn_path)) {
    utils::read.delim(profile_warn_path, stringsAsFactors = FALSE)
  } else NULL

  list(
    features_path = file.path(output_dir, output_file),
    ms2_path = if (!is.null(precursor_watchlist)) file.path(output_dir, ms2_output_file) else NULL,
    profile_mode_files = profile_mode_files,
    log = status
  )
}

## ---- Reading Python outputs back into R -------------------------------------
# `sample` values that look like bare integers (e.g. Shiny renames uploads
# to numeric temp names like "0.mzML"/"1.mzML", so Python's sample column
# ends up "0"/"1") get silently type-converted to numeric by read.delim()'s
# default type inference. Left as-is, that breaks any later name-based
# lookup like `name_map[feats$sample]` (R indexes a named vector by a
# NUMERIC vector positionally, not by name) -- coerce back to character
# right after reading so `sample` is always a stable join/lookup key.
read_batch_features <- function(tsv_path) {
  if (!file.exists(tsv_path)) stop("Feature table not found: ", tsv_path)
  df <- utils::read.delim(tsv_path, stringsAsFactors = FALSE)
  df$sample <- as.character(df$sample)
  df
}

read_batch_ms2 <- function(tsv_path) {
  empty <- data.frame(sample = character(), ms2_scan_id = character(), rt = numeric(),
                       precursor_mz = numeric(), precursor_z = integer(),
                       mz = numeric(), intensity = numeric(), stringsAsFactors = FALSE)
  if (is.null(tsv_path) || !file.exists(tsv_path)) return(empty)
  raw <- utils::read.delim(tsv_path, stringsAsFactors = FALSE)
  if (nrow(raw) == 0) return(empty)
  raw$sample <- as.character(raw$sample)  # see read_batch_features() for why

  rows <- lapply(seq_len(nrow(raw)), function(i) {
    mzs <- suppressWarnings(as.numeric(strsplit(raw$mz_list[i], ";", fixed = TRUE)[[1]]))
    ints <- suppressWarnings(as.numeric(strsplit(raw$intensity_list[i], ";", fixed = TRUE)[[1]]))
    n <- min(length(mzs), length(ints))
    if (n == 0) return(NULL)
    data.frame(sample = raw$sample[i], ms2_scan_id = as.character(raw$ms2_scan_id[i]),
               rt = raw$rt[i], precursor_mz = raw$precursor_mz[i],
               precursor_z = raw$precursor_z[i],
               mz = mzs[seq_len(n)], intensity = ints[seq_len(n)],
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

## ---- Multi-sample MS1 matching (reuses match_ms1() verbatim) ---------------
match_ms1_batch <- function(mets, features, dict = STANDARD_DICT,
                             ppm_tol = 10, z_range = 3:12,
                             adducts = c("H", "Na", "K", "NH4"),
                             max_oxid = 6, h_offset = 0,
                             n_iso = 5, use_envipat = TRUE) {
  if (is.null(features) || nrow(features) == 0) return(data.frame())
  samples <- unique(features$sample)
  out <- lapply(samples, function(s) {
    feat_s <- features[features$sample == s, c("mz", "rt", "max_intensity", "n_scans"), drop = FALSE]
    m <- match_ms1(mets, feat_s, dict, ppm_tol, z_range, adducts,
                    max_oxid, h_offset, n_iso, use_envipat)
    if (nrow(m) > 0) m$sample <- s
    m
  })
  out <- out[vapply(out, nrow, integer(1)) > 0]
  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

## ---- Retained unidentified peaks -------------------------------------------
# match_ms1() doesn't expose which feature row it picked as "best" for a
# given theoretical candidate, but it does copy `rt` and `intensity`
# (= max_intensity) through UNROUNDED from ms1_features -- so an exact-
# equality join on (sample, rt, intensity) reliably identifies which
# original features were used, without touching match_ms1() itself.
unmatched_features_batch <- function(features, ms1_matches) {
  if (is.null(features) || nrow(features) == 0) return(features)
  if (is.null(ms1_matches) || nrow(ms1_matches) == 0) return(features)
  matched_key <- paste(ms1_matches$sample, ms1_matches$rt, ms1_matches$intensity)
  feat_key <- paste(features$sample, features$rt, features$max_intensity)
  features[!feat_key %in% matched_key, , drop = FALSE]
}

## ---- MS2 confirmation of matched hits (reuses confirm_metabolite() verbatim)
confirm_ms2_batch <- function(mets, ms1_matches, ms2_by_sample, dict = STANDARD_DICT,
                               frag_tol_ppm = 25, frag_z_range = 1:2, h_offset = 0,
                               ms2_lookup_ppm_tol = 20, include_internal = FALSE) {
  if (is.null(ms1_matches) || nrow(ms1_matches) == 0) return(data.frame())
  if (is.null(ms2_by_sample) || nrow(ms2_by_sample) == 0) return(data.frame())

  unique_hits <- unique(ms1_matches[, c("sample", "met_id", "met_name", "k_oxid", "z", "adduct", "theo_mz")])
  results <- lapply(seq_len(nrow(unique_hits)), function(i) {
    hit <- unique_hits[i, ]
    ms2_s <- ms2_by_sample[ms2_by_sample$sample == hit$sample,
                            c("rt", "precursor_mz", "precursor_z", "mz", "intensity"), drop = FALSE]
    if (nrow(ms2_s) == 0) return(NULL)
    spectra <- find_ms2_spectra(ms2_s, hit$theo_mz, hit$z, ms2_lookup_ppm_tol)
    if (length(spectra) == 0) return(NULL)

    best_spec <- spectra[[which.max(vapply(spectra, nrow, integer(1)))]]
    met <- mets[[which(vapply(mets, function(m) m$id == hit$met_id, logical(1)))]]
    conf <- confirm_metabolite(met, best_spec, dict, tol_ppm = frag_tol_ppm,
                                z_range = frag_z_range, include_internal = include_internal,
                                h_offset = h_offset)

    data.frame(sample = hit$sample, met_id = hit$met_id, met_name = hit$met_name,
               k_oxid = hit$k_oxid, z = hit$z, adduct = hit$adduct,
               n_ms2_peaks = nrow(best_spec), n_frag_matches = conf$score$n_matches,
               coverage = conf$score$coverage, confirmation_score = conf$score$total_score,
               n_diagnostics = conf$score$n_diagnostics, confident = conf$score$confident,
               stringsAsFactors = FALSE)
  })
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0) return(data.frame())
  do.call(rbind, results)
}

## ---- Full batch annotation pipeline -----------------------------------------
annotate_metabolites_batch <- function(mets, features, ms2_by_sample = NULL,
                                        dict = STANDARD_DICT, ppm_tol = 10,
                                        z_range = 3:12, adducts = c("H", "Na", "K", "NH4"),
                                        max_oxid = 6, h_offset = 0, n_iso = 5,
                                        use_envipat = TRUE, frag_tol_ppm = 25,
                                        frag_z_range = 1:2, include_internal = FALSE) {
  ms1_matches <- match_ms1_batch(mets, features, dict, ppm_tol, z_range, adducts,
                                  max_oxid, h_offset, n_iso, use_envipat)

  env <- data.frame()
  if (nrow(ms1_matches) > 0) {
    env <- do.call(rbind, lapply(split(ms1_matches, ms1_matches$sample), function(g) {
      e <- envelope_consistency(g, h_offset = h_offset)
      if (nrow(e) > 0) e$sample <- g$sample[1]
      e
    }))
  }

  unmatched <- unmatched_features_batch(features, ms1_matches)

  ms2_conf <- data.frame()
  if (!is.null(ms2_by_sample) && nrow(ms1_matches) > 0) {
    ms2_conf <- confirm_ms2_batch(mets, ms1_matches, ms2_by_sample, dict,
                                   frag_tol_ppm = frag_tol_ppm, frag_z_range = frag_z_range,
                                   h_offset = h_offset, include_internal = include_internal)
  }

  list(ms1_matches = ms1_matches, envelope = env,
       unmatched = unmatched, ms2_confirmations = ms2_conf)
}
