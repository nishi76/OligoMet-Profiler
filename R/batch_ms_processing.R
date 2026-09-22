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
# Sys.which() only checks that SOMETHING on PATH resolves to that name -- on
# Windows that's not the same as it being runnable. Official python.org
# installers only ship python.exe, not python3.exe, so a `python3` entry on
# a Windows PATH is often a non-functional Microsoft Store app-execution-
# alias stub, or a stray non-Windows-executable shim (e.g. from Git Bash's
# usr/bin) -- Sys.which("python3") still finds it, but system2() then fails
# with Windows' generic "command not found" exit status 9009, which looks
# identical to Python genuinely being absent. Actually invoking each
# candidate with --version (not just resolving its path) catches this
# before run_batch_deconvolution() ever gets to the real work.
#' Locate a working Python 3 interpreter
#'
#' Tries `"python3"` then `"python"` on `PATH`, actually invoking each
#' with `--version` (not just resolving its path via `Sys.which()`) to
#' rule out a resolvable-but-non-functional entry -- see the comment
#' above for why that check matters on Windows in particular.
#'
#' @return `"python3"` or `"python"` (whichever works), or `NA_character_`
#'   if neither does. Used as the default `python_bin` for
#'   [run_batch_deconvolution()].
#' @seealso [run_batch_deconvolution()]
#' @export
find_python <- function() {
  for (bin in c("python3", "python")) {
    if (!nzchar(Sys.which(bin))) next
    ok <- tryCatch({
      out <- suppressWarnings(system2(bin, "--version", stdout = TRUE, stderr = TRUE))
      is.null(attr(out, "status")) || attr(out, "status") == 0
    }, error = function(e) FALSE)
    if (ok) return(bin)
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
#' Write a precursor m/z watch-list for targeted MS2 capture
#'
#' Writes the (H-adduct-only) theoretical precursor m/z list from
#' [prm_inclusion_list()] to a plain text file, one value per line, for
#' [run_batch_deconvolution()]'s `precursor_watchlist` argument. Only
#' needs to be a broad net -- the Python side's own `ms2_watch_ppm`
#' widens the match, and DDA fragmentation is triggered off the dominant
#' charge/adduct envelope in practice regardless.
#'
#' @param mets A list of metabolite objects (see [generate_metabolites()]).
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @param z_range,max_oxid,h_offset Passed to [prm_inclusion_list()].
#' @param out_path File path to write the watch-list to.
#' @return `out_path`, returned for chaining into
#'   [run_batch_deconvolution()]'s `precursor_watchlist` argument.
#' @seealso [prm_inclusion_list()], [run_batch_deconvolution()]
#' @export
write_precursor_watchlist <- function(mets, dict = STANDARD_DICT, z_range = 3:12,
                                       max_oxid = 6, h_offset = 0, out_path) {
  prm <- prm_inclusion_list(mets, dict, z_range = z_range, h_offset = h_offset,
                             max_oxid = max_oxid)
  writeLines(as.character(unique(prm$precursor_mz)), out_path)
  out_path
}

## ---- Invoking the batch deconvolution CLI ----------------------------------
#' Run the parallel Python charge-envelope deconvolution pipeline
#'
#' Shells out to `python -m oligomet_deconv.cli` (see
#' `inst/python/oligomet_deconv/`) to run ROI/charge-envelope
#' deconvolution across many mzML/mzXML files in parallel, then reads the
#' resulting per-file sidecars back into R.
#'
#' @param files Paths to mzML/mzXML files (already resolved -- see
#'   [resolve_ms_input_file()] for vendor-format conversion).
#' @param output_dir Directory the Python side writes its outputs into.
#' @param output_file,ms2_output_file Feature/MS2 table filenames within
#'   `output_dir`. `ms2_output_file` defaults to `"combined_ms2.tsv"`
#'   when `precursor_watchlist` is given.
#' @param precursor_watchlist Optional path to a precursor-mz watch-list
#'   (see [write_precursor_watchlist()]) enabling targeted MS2 capture.
#' @param ms2_watch_ppm,roi_ppm,rt_tol,mass_tol_ppm,z_range,min_scans,max_gap_scans,min_charge_states
#'   Deconvolution parameters passed straight through to the Python CLI
#'   -- see `inst/python/oligomet_deconv/cli.py`'s `--help` for what each
#'   one does.
#' @param min_intensity Fixed absolute intensity floor for a candidate
#'   ROI peak; used as a fallback when `sn_threshold` is `NULL` or a
#'   file's noise estimate comes back unavailable.
#' @param sn_threshold When given, OVERRIDES `min_intensity` with a
#'   PER-FILE threshold = (that file's own noise level) * `sn_threshold`
#'   -- see `_noise_thresholds.tsv` in the return value for what actually
#'   got applied per file.
#' @param n_workers Parallel worker processes; `NULL` lets the Python
#'   side pick (`os.cpu_count() - 1`).
#' @param python_bin,module_dir Interpreter and module location, normally
#'   left at their defaults ([find_python()]/`.find_deconv_module_dir()`).
#' @param progress,console_tracker Optional callbacks: `progress(msg)`
#'   for a one-line status update before the subprocess runs;
#'   `console_tracker(text)` for the subprocess's full captured
#'   stdout/stderr afterward.
#' @return A list: `features_path`, `ms2_path` (or `NULL`),
#'   `profile_mode_files` (data.frame, or `NULL` -- files that appear to
#'   be uncentroided profile-mode data), `noise_thresholds` (data.frame,
#'   or `NULL` -- only present in S/N-threshold mode, one row per file
#'   with its computed `noise_level`/`effective_min_intensity`), and
#'   `log` (the subprocess's captured output).
#' @seealso [read_batch_features()], [read_batch_ms2()],
#'   [match_ms1_batch()], [annotate_metabolites_batch()]
#' @export
run_batch_deconvolution <- function(files, output_dir = tempdir(),
                                     output_file = "combined_features.tsv",
                                     ms2_output_file = NULL,
                                     precursor_watchlist = NULL, ms2_watch_ppm = 50,
                                     roi_ppm = 15, rt_tol = 0.15, mass_tol_ppm = 20,
                                     z_range = 3:20, min_intensity = 1e4, sn_threshold = NULL,
                                     min_scans = 3,
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
  # sn_threshold OVERRIDES --min-intensity on the Python side (see
  # DeconvParams/ROIBuilder in inst/python/oligomet_deconv/) with a
  # threshold derived from each file's OWN noise level -- still passed
  # regardless, since a file where the noise estimate comes back NaN
  # (e.g. an empty file) falls back to --min-intensity there.
  if (!is.null(sn_threshold)) args <- c(args, "--sn-threshold", sn_threshold)
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

  # Only produced when sn_threshold was set -- records the per-file noise
  # level and the resulting absolute threshold actually applied, since the
  # whole point of sn_threshold is that the same multiple produces a
  # DIFFERENT number on every file depending on its own background level.
  noise_path <- file.path(output_dir, "_noise_thresholds.tsv")
  noise_thresholds <- if (file.exists(noise_path)) {
    utils::read.delim(noise_path, stringsAsFactors = FALSE)
  } else NULL

  list(
    features_path = file.path(output_dir, output_file),
    ms2_path = if (!is.null(precursor_watchlist)) file.path(output_dir, ms2_output_file) else NULL,
    profile_mode_files = profile_mode_files,
    noise_thresholds = noise_thresholds,
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
#' Read the combined feature table from a batch deconvolution run
#'
#' @param tsv_path Path to the feature TSV (see
#'   [run_batch_deconvolution()]'s `features_path`).
#' @return A data.frame with `sample` forced to character (see the
#'   comment above for why -- numeric-looking sample names, e.g. Shiny's
#'   renamed uploads, would otherwise silently become a numeric column).
#' @seealso [run_batch_deconvolution()], [read_batch_ms2()],
#'   [match_ms1_batch()]
#' @export
read_batch_features <- function(tsv_path) {
  if (!file.exists(tsv_path)) stop("Feature table not found: ", tsv_path)
  df <- utils::read.delim(tsv_path, stringsAsFactors = FALSE)
  df$sample <- as.character(df$sample)
  df
}

#' Read the combined MS2 table from a batch deconvolution run
#'
#' Unpacks the semicolon-delimited `mz_list`/`intensity_list` columns the
#' Python side writes (one row per MS2 scan) into one row per peak.
#'
#' @param tsv_path Path to the MS2 TSV (see [run_batch_deconvolution()]'s
#'   `ms2_path`), or `NULL`/nonexistent (returns the empty shape below).
#' @return A data.frame: `sample`, `ms2_scan_id`, `rt`, `precursor_mz`,
#'   `precursor_z`, `mz`, `intensity` -- one row per peak. Empty
#'   data.frame (same columns) if `tsv_path` is `NULL`, missing, or has
#'   no rows.
#' @seealso [run_batch_deconvolution()], [confirm_ms2_batch()]
#' @export
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
#' Match MS1 features against a metabolite library, per sample
#'
#' Thin per-sample loop over [match_ms1()]: splits a multi-sample feature
#' table by `sample`, runs the same single-sample matching logic on each
#' (theoretical metabolite x oxidation x charge x adduct enumeration
#' against that sample's features), and stacks the results back together
#' with a `sample` column. No new matching math of its own.
#'
#' @param mets A list of metabolite objects (see [generate_metabolites()]).
#' @param features A multi-sample feature table (see
#'   [read_batch_features()]) with `sample`, `mz`, `rt`, `max_intensity`,
#'   `n_scans`, and optionally `area` columns.
#' @param dict,ppm_tol,z_range,adducts,max_oxid,h_offset,n_iso,use_envipat
#'   Passed straight through to [match_ms1()] for each sample.
#' @return A data.frame, one row per (sample, metabolite, oxidation,
#'   charge, adduct) match within tolerance, with the same columns
#'   [match_ms1()] produces plus `sample`. Empty data.frame if `features`
#'   is empty or nothing matches.
#' @seealso [match_ms1()], [annotate_metabolites_batch()],
#'   [unmatched_features_batch()]
#' @export
match_ms1_batch <- function(mets, features, dict = STANDARD_DICT,
                             ppm_tol = 10, z_range = 3:12,
                             adducts = c("H", "Na", "K", "NH4"),
                             max_oxid = 6, h_offset = 0,
                             n_iso = 5, use_envipat = TRUE) {
  if (is.null(features) || nrow(features) == 0) return(data.frame())
  samples <- unique(features$sample)
  # `area` (trapezoidal AUC, computed by the Python ROI pipeline -- see
  # roi.py) is only present when features came from run_batch_deconvolution();
  # a peak-list-derived feature table won't have it. Pass it through when
  # available so match_ms1() can carry it onto ms1_matches.
  feat_cols <- c("mz", "rt", "max_intensity", "n_scans",
                 if ("area" %in% names(features)) "area")
  out <- lapply(samples, function(s) {
    feat_s <- features[features$sample == s, feat_cols, drop = FALSE]
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
#' Retain the batch features that matched no theoretical metabolite
#'
#' [match_ms1()] doesn't expose which feature row it picked as "best"
#' for a given theoretical candidate, but it does copy `rt` and
#' `intensity` (= `max_intensity`) through UNROUNDED from the original
#' features -- so an exact-equality join on (`sample`, `rt`, `intensity`)
#' reliably identifies which original features were used, without
#' needing to touch [match_ms1()] itself.
#'
#' @param features A multi-sample feature table (see
#'   [read_batch_features()]).
#' @param ms1_matches Output of [match_ms1_batch()]/
#'   [annotate_metabolites_batch()].
#' @return `features`, filtered down to rows with no corresponding row
#'   in `ms1_matches` -- retained for QC/follow-up rather than discarded.
#' @seealso [annotate_metabolites_batch()]
#' @export
unmatched_features_batch <- function(features, ms1_matches) {
  if (is.null(features) || nrow(features) == 0) return(features)
  if (is.null(ms1_matches) || nrow(ms1_matches) == 0) return(features)
  matched_key <- paste(ms1_matches$sample, ms1_matches$rt, ms1_matches$intensity)
  feat_key <- paste(features$sample, features$rt, features$max_intensity)
  features[!feat_key %in% matched_key, , drop = FALSE]
}

## ---- MS2 confirmation of matched hits (reuses confirm_metabolite() verbatim)
#' Confirm batch MS1 hits with their acquired MS2 spectra
#'
#' For each unique (sample, metabolite, oxidation, charge, adduct) MS1
#' hit, finds its acquired MS2 spectrum ([find_ms2_spectra()]) and scores
#' fragment confirmation against it ([confirm_metabolite()]) -- unlike
#' [annotate_metabolites()]'s single-file path, this searches the exact
#' adduct the MS1 hit matched (since [match_ms1_batch()] tracks it per
#' hit), not the full requested adduct set.
#'
#' @param mets A list of metabolite objects (see [generate_metabolites()]).
#' @param ms1_matches Output of [match_ms1_batch()].
#' @param ms2_by_sample A multi-sample MS2 table (see [read_batch_ms2()]).
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @param frag_tol_ppm,frag_z_range,include_internal Fragment matching
#'   parameters -- see [confirm_metabolite()].
#' @param h_offset Charge-envelope offset -- see [match_ms1()].
#' @param ms2_lookup_ppm_tol Precursor m/z tolerance for finding the MS2
#'   spectrum to confirm against -- see [find_ms2_spectra()].
#' @return A data.frame, one row per confirmed (sample, metabolite,
#'   oxidation, charge, adduct) hit: `sample`, `met_id`, `met_name`,
#'   `k_oxid`, `z`, `adduct`, `n_ms2_peaks`, `n_frag_matches`,
#'   `coverage`, `confirmation_score`, `n_diagnostics`, `confident`.
#'   Also carries a `"spectra"` attribute: a named list (keyed
#'   `"sample|met_id|k_oxid|z|adduct"`) of the acquired spectrum used for
#'   each row's confirmation, for callers (the mirror-plot UI) that want
#'   to re-render the comparison without re-running
#'   [find_ms2_spectra()] -- this rides along as an attribute rather than
#'   changing the return shape, so `nrow()`/`$sample`/column-selection
#'   callers keep working unmodified.
#' @seealso [annotate_metabolites_batch()], [confirm_metabolite()],
#'   [plot_mirror_spectrum()]
#' @export
confirm_ms2_batch <- function(mets, ms1_matches, ms2_by_sample, dict = STANDARD_DICT,
                               frag_tol_ppm = 25, frag_z_range = 1:2, h_offset = 0,
                               ms2_lookup_ppm_tol = 20, include_internal = FALSE) {
  if (is.null(ms1_matches) || nrow(ms1_matches) == 0) return(data.frame())
  if (is.null(ms2_by_sample) || nrow(ms2_by_sample) == 0) return(data.frame())

  unique_hits <- unique(ms1_matches[, c("sample", "met_id", "met_name", "k_oxid", "z", "adduct", "theo_mz")])
  spectra <- list()
  results <- lapply(seq_len(nrow(unique_hits)), function(i) {
    hit <- unique_hits[i, ]
    ms2_s <- ms2_by_sample[ms2_by_sample$sample == hit$sample,
                            c("rt", "precursor_mz", "precursor_z", "mz", "intensity"), drop = FALSE]
    if (nrow(ms2_s) == 0) return(NULL)
    spec <- find_ms2_spectra(ms2_s, hit$theo_mz, hit$z, ms2_lookup_ppm_tol)
    if (length(spec) == 0) return(NULL)

    best_spec <- spec[[which.max(vapply(spec, nrow, integer(1)))]]
    met <- mets[[which(vapply(mets, function(m) m$id == hit$met_id, logical(1)))]]
    # Unlike annotate_metabolites()'s single-file path, unique_hits here
    # already carries the exact adduct this MS1 hit matched (from
    # match_ms1_batch()), so MS2 fragment matching can search that one
    # adduct precisely rather than the full requested set.
    conf <- confirm_metabolite(met, best_spec, dict, tol_ppm = frag_tol_ppm,
                                z_range = frag_z_range, include_internal = include_internal,
                                h_offset = h_offset, adducts = hit$adduct)

    key <- paste(hit$sample, hit$met_id, hit$k_oxid, hit$z, hit$adduct, sep = "|")
    spectra[[key]] <<- best_spec

    data.frame(sample = hit$sample, met_id = hit$met_id, met_name = hit$met_name,
               k_oxid = hit$k_oxid, z = hit$z, adduct = hit$adduct,
               n_ms2_peaks = nrow(best_spec), n_frag_matches = conf$score$n_matches,
               coverage = conf$score$coverage, confirmation_score = conf$score$total_score,
               n_diagnostics = conf$score$n_diagnostics, confident = conf$score$confident,
               stringsAsFactors = FALSE)
  })
  results <- results[!vapply(results, is.null, logical(1))]
  out <- if (length(results) == 0) data.frame() else do.call(rbind, results)
  attr(out, "spectra") <- spectra
  out
}

## ---- Full batch annotation pipeline -----------------------------------------
#' Run the full batch annotation pipeline: MS1 match, envelope, MS2, degradation
#'
#' The end-to-end batch analysis step: MS1 matching ([match_ms1_batch()]),
#' per-sample charge-envelope consistency, unmatched-peak retention
#' ([unmatched_features_batch()]), optional MS2 confirmation
#' ([confirm_ms2_batch()]), and an optional degradation summary
#' ([degradation_summary()]) -- everything [run_batch_deconvolution()]'s
#' output needs before it's ready for statistics/quantification or
#' export.
#'
#' @param mets A list of metabolite objects (see [generate_metabolites()]).
#' @param features A multi-sample feature table (see
#'   [read_batch_features()]).
#' @param ms2_by_sample Optional multi-sample MS2 table (see
#'   [read_batch_ms2()]) for MS2 confirmation; `NULL` skips that step.
#' @param dict,ppm_tol,z_range,adducts,max_oxid,h_offset,n_iso,use_envipat
#'   MS1 matching parameters -- see [match_ms1()].
#' @param frag_tol_ppm,frag_z_range,include_internal MS2 confirmation
#'   parameters -- see [confirm_metabolite()].
#' @param compute_degradation Whether to also compute
#'   [degradation_summary()] from the MS1 matches.
#' @param sample_meta Optional sample metadata (`sample`, `sample_type`,
#'   `group`/`timepoint`) forwarded to [degradation_summary()] so
#'   calibration standards/QC/blanks are excluded from it and its
#'   per-sample/composition tables carry group/timepoint context. `NULL`
#'   (the default) keeps every sample, same as calling
#'   [degradation_summary()] with no `sample_meta` of its own.
#' @return A list: `ms1_matches`, `envelope` (charge-envelope consistency,
#'   per sample), `unmatched` (retained unmatched features),
#'   `ms2_confirmations`, `ms2_spectra` (named list keyed
#'   `"sample|met_id|k_oxid|z|adduct"`, for mirror-plot rendering), and
#'   `degradation` (`NULL` if `compute_degradation = FALSE` or there are
#'   no MS1 matches).
#' @seealso [match_ms1_batch()], [run_batch_deconvolution()],
#'   [quantify_metabolites()]
#' @export
annotate_metabolites_batch <- function(mets, features, ms2_by_sample = NULL,
                                        dict = STANDARD_DICT, ppm_tol = 10,
                                        z_range = 3:12, adducts = c("H", "Na", "K", "NH4"),
                                        max_oxid = 6, h_offset = 0, n_iso = 5,
                                        use_envipat = TRUE, frag_tol_ppm = 25,
                                        frag_z_range = 1:2, include_internal = FALSE,
                                        compute_degradation = TRUE, sample_meta = NULL) {
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
  ms2_spectra <- list()
  if (!is.null(ms2_by_sample) && nrow(ms1_matches) > 0) {
    ms2_conf <- confirm_ms2_batch(mets, ms1_matches, ms2_by_sample, dict,
                                   frag_tol_ppm = frag_tol_ppm, frag_z_range = frag_z_range,
                                   h_offset = h_offset, include_internal = include_internal)
    ms2_spectra <- attr(ms2_conf, "spectra") %||% list()
  }

  degradation <- if (compute_degradation) degradation_summary(ms1_matches, sample_meta = sample_meta) else NULL

  list(ms1_matches = ms1_matches, envelope = env,
       unmatched = unmatched, ms2_confirmations = ms2_conf,
       ms2_spectra = ms2_spectra, degradation = degradation)
}
