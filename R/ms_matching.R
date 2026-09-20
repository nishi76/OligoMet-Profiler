# =============================================================================
# ms_matching.R
# MS data import and targeted matching for oligonucleotide metabolite ID.
#
# Imports:
#   - mzML files (parsed with xml2 + native base64/zlib decoding of binary data)
#   - Simple text/CSV peak lists (mz, intensity columns)
#   - Vendor raw via msconvert bridge (if msconvert is on PATH)
#
# Matching:
#   - MS1 targeted matching: metabolite m/z vs observed features within ppm tol.
#     Nearest-feature (which.min()) matching reports the single closest
#     target as if it were unambiguous; match_ms1() also flags how many
#     DISTINCT metabolites' theoretical m/z fall within tolerance of the
#     same observed feature (n_candidates/ambiguous columns) -- at the
#     library sizes and tolerances this app targets, most matches have
#     more than one candidate, and that isn't otherwise visible in the
#     results table.
#   - Isotope pattern fit: observed vs theoretical isotope envelope (cosine sim)
#   - Charge envelope consistency: multiple charge states form consistent envelope
#   - MS2 fragment matching: McLuckey fragments vs MS/MS peaks
#   - PS diagnostic ions: m/z 94.9452, 192.9746
#
# Scoring (Standard mode, per user spec):
#   - 10 ppm mass tolerance
#   - Charge states z = 3-12
#   - Adducts: H, Na, K, NH4
#   - Isotope match (cosine similarity)
#   - Charge-envelope consistency
#   - Scored by mass error + isotope fit + envelope consistency
#
# Grounded in:
#   - Eluforsen (Kim et al., 2019): LC-MS/MS method, PS diagnostic ions
#   - FMVS (Ye et al., 2025): total confirmation score + sequence coverage
#   - OligoDistiller (Liu et al., 2025): isotope envelope goodness-of-fit
# =============================================================================

## ---- Vendor raw bridge ----------------------------------------------------
#' Locate ProteoWizard's msconvert on this machine
#'
#' Checks `PATH` first, then a few common Windows/macOS/Linux install
#' locations. Used by [convert_vendor_raw()]/[resolve_ms_input_file()] to
#' decide whether vendor raw files (Thermo `.raw`, Sciex `.wiff`,
#' Bruker `.baf`/`.yep`) can be converted automatically.
#'
#' @return The path to `msconvert`, or `NA` if it isn't found anywhere
#'   checked (install from proteowizard.org).
#' @seealso [convert_vendor_raw()]
#' @export
find_msconvert <- function() {
  path <- Sys.which("msconvert")
  if (nzchar(path)) return(path)
  # Check common Windows/MSI install paths (won't exist in sandbox, but for user)
  candidates <- c(
    "/usr/local/bin/msconvert",
    "C:/Program Files/ProteoWizard/ProteoWizard 3.0/msconvert.exe",
    "/Applications/ProteoWizard/msconvert"
  )
  for (c in candidates) if (file.exists(c)) return(c)
  NA
}

#' Convert a vendor raw file to mzML via ProteoWizard msconvert
#'
#' Shells out to `msconvert` (found via [find_msconvert()]) to convert a
#' vendor raw file to `.mzML`, optionally centroiding it. If `msconvert`
#' isn't found, prints the equivalent command to run externally and
#' returns `NA` rather than erroring, so a caller (e.g.
#' [resolve_ms_input_file()]) can decide how to handle it.
#'
#' @param raw_file Path to the vendor raw file.
#' @param output_dir Directory to write the converted `.mzML` into.
#' @param centroid Whether to pass `--centroid` to msconvert. ROI/
#'   charge-envelope detection ([extract_ms1_features()],
#'   the Python batch pipeline) is designed for centroided peaks.
#' @param peak_picking Currently unused (reserved).
#' @return The path to the converted `.mzML` file, or `NA` if msconvert
#'   isn't available or the conversion didn't produce the expected
#'   output file.
#' @seealso [find_msconvert()], [resolve_ms_input_file()]
#' @export
convert_vendor_raw <- function(raw_file, output_dir = tempdir(),
                                centroid = TRUE, peak_picking = TRUE) {
  msconvert <- find_msconvert()
  if (is.na(msconvert)) {
    message("msconvert (ProteoWizard) not found. ",
            "Please convert vendor raw files to mzML externally:\n",
            "  msconvert <raw_file> --mzML --centroid --filter 'peakPicking true 1-'")
    return(NA)
  }
  out_file <- file.path(output_dir, sub("\\.[^.]+$", ".mzML", basename(raw_file)))
  cmd <- sprintf("%s \"%s\" --mzML -o \"%s\"", msconvert, raw_file, output_dir)
  if (centroid) cmd <- paste0(cmd, " --centroid")
  system(cmd, intern = FALSE)
  if (file.exists(out_file)) out_file else NA
}

## ---- mzML parser (lightweight, xml2-based) --------------------------------
#' Parse an mzML/mzXML file with a hand-rolled xml2 parser
#'
#' Reads every spectrum's m/z/intensity binary data arrays (base64 +
#' zlib, decoded natively -- no mzR/compiled dependency) and separates
#' MS1 peaks from MS2 spectra. This is the automatic fallback
#' [read_ms_file()] uses when the (much faster, indexed) Spectra/mzR
#' backend isn't available or fails -- loads the whole file into memory
#' as an xml2 DOM, so it doesn't scale to very large files the way
#' [read_ms_file()]'s primary path or the Python batch pipeline's
#' streaming reader do.
#'
#' @param file Path to the mzML or mzXML file.
#' @return A list: `ms1` (data.frame of `rt`, `mz`, `intensity`, one row
#'   per peak across every MS1 spectrum), `ms2` (data.frame of `rt`,
#'   `precursor_mz`, `precursor_z`, `mz`, `intensity`), and `info` (list:
#'   `file`, `n_spectra`, `n_ms1`, `n_ms2`, `instrument`, and
#'   `profile_mode` -- `TRUE`/`FALSE` once the first MS1 spectrum's
#'   profile/centroid cvParam settles it, else `NA`).
#' @seealso [read_ms_file()], [extract_ms1_features()]
#' @export
parse_mzml <- function(file) {
  if (!file.exists(file)) stop("mzML file not found: ", file)
  doc <- xml2::read_xml(file)

  # mzML/mzXML normally declare a default namespace (e.g.
  # xmlns="http://psi.hupo.org/ms/mzml"), which an unprefixed xpath tag
  # test (e.g. ".//spectrum") never matches -- XPath 1.0 treats an
  # unprefixed name test as "no namespace" regardless of any default
  # namespace in scope. Matching on local-name() instead makes every
  # lookup below namespace-agnostic, so it works whether or not the
  # source file declares one.
  find_all <- function(node, tag)
    xml2::xml_find_all(node, sprintf(".//*[local-name()='%s']", tag))
  find_first <- function(node, tag)
    xml2::xml_find_first(node, sprintf(".//*[local-name()='%s']", tag))

  # Helper: get cvParam value by accession
  get_cv <- function(node, accession) {
    cv <- xml2::xml_find_first(node,
      sprintf(".//*[local-name()='cvParam' and @accession='%s']", accession))
    if (!is.na(cv)) {
      v <- xml2::xml_attr(cv, "value")
      if (!is.na(v) && nzchar(v)) return(v)
    }
    NA
  }

  # Helper: check if cvParam exists
  has_cv <- function(node, accession) {
    !is.na(xml2::xml_find_first(node,
      sprintf(".//*[local-name()='cvParam' and @accession='%s']", accession)))
  }

  # Find all spectra
  spectra <- find_all(doc, "spectrum")
  n_spectra <- length(spectra)
  if (n_spectra == 0) {
    # Try mzXML format
    spectra <- find_all(doc, "scan")
    if (length(spectra) == 0) stop("No spectra found in mzML/mzXML file")
  }

  ms1_peaks <- list()
  ms2_spectra <- list()
  profile_mode_detected <- NA  # NA = undetermined; TRUE/FALSE once the first MS1 spectrum settles it

  for (i in seq_along(spectra)) {
    sp <- spectra[[i]]
    # Determine MS level
    ms_level <- get_cv(sp, "MS:1000511")  # MS:1000511 = ms level
    if (is.na(ms_level)) {
      # Try defaultMSLevel or scan attribute
      ms_level <- get_cv(sp, "MS:1000511")
      if (is.na(ms_level)) ms_level <- "1"
    }
    ms_level <- as.integer(ms_level)

    if (ms_level == 1 && is.na(profile_mode_detected)) {
      # MS:1000128/MS:1000127 (profile/centroid spectrum) are presence-only
      # flag cvParams -- has_cv() checks for the node existing, not its
      # (nonexistent) value, so this works the same way as everywhere else
      # has_cv() is used above.
      if (has_cv(sp, "MS:1000128")) profile_mode_detected <- TRUE
      else if (has_cv(sp, "MS:1000127")) profile_mode_detected <- FALSE
    }

    # Get retention time (scan start time, MS:1000016)
    rt <- get_cv(sp, "MS:1000016")
    if (is.na(rt)) {
      # Try in scanWindow or parent scan
      rt <- get_cv(find_first(sp, "scan"), "MS:1000016")
    }
    rt <- as.numeric(rt) %||% NA_real_

    # Extract binary data arrays
    bdas <- find_all(sp, "binaryDataArray")
    mz_arr <- NULL
    int_arr <- NULL
    for (bda in bdas) {
      # Check if this is m/z or intensity array
      is_mz <- has_cv(bda, "MS:1000514")  # m/z array
      is_int <- has_cv(bda, "MS:1000515")  # intensity array
      # Check precision
      is_32 <- has_cv(bda, "MS:1000521")  # 32-bit float
      is_64 <- has_cv(bda, "MS:1000523")  # 64-bit float
      dtype <- if (is_64) "64" else "32"
      # Check compression
      is_zlib <- has_cv(bda, "MS:1000574")  # zlib compression

      # Get binary data
      bin_node <- find_first(bda, "binary")
      if (is.na(bin_node)) next
      b64 <- xml2::xml_text(bin_node)
      if (!nzchar(b64)) next

      vals <- .decode_binary(b64, dtype, is_zlib)
      if (is_mz) mz_arr <- vals
      else if (is_int) int_arr <- vals
    }

    if (is.null(mz_arr) || is.null(int_arr) || length(mz_arr) == 0) next

    if (ms_level == 1) {
      ms1_peaks[[length(ms1_peaks) + 1]] <- data.frame(
        rt = rt, mz = mz_arr, intensity = int_arr,
        stringsAsFactors = FALSE)
    } else if (ms_level == 2) {
      # Get precursor info
      prec <- find_first(sp, "precursor")
      prec_mz <- NA_real_; prec_z <- NA_integer_
      if (!is.na(prec)) {
        prec_mz <- as.numeric(get_cv(prec, "MS:1000744"))  # selected ion m/z
        if (is.na(prec_mz)) prec_mz <- as.numeric(get_cv(prec, "MS:1000040"))  # isolation m/z
        prec_z <- as.integer(get_cv(prec, "MS:1000041"))  # charge state
      }
      ms2_spectra[[length(ms2_spectra) + 1]] <- data.frame(
        rt = rt, precursor_mz = prec_mz, precursor_z = prec_z,
        mz = mz_arr, intensity = int_arr,
        stringsAsFactors = FALSE)
    }
  }

  # Combine MS1 peaks
  ms1_df <- if (length(ms1_peaks) > 0) do.call(rbind, ms1_peaks) else data.frame()
  # Combine MS2 spectra
  ms2_df <- if (length(ms2_spectra) > 0) do.call(rbind, ms2_spectra) else data.frame()

  # File info
  instrument <- get_cv(doc, "MS:1000031")  # instrument model
  info <- list(
    file = file, n_spectra = n_spectra,
    n_ms1 = length(ms1_peaks), n_ms2 = length(ms2_spectra),
    instrument = if (is.na(instrument)) "unknown" else instrument,
    profile_mode = profile_mode_detected
  )

  list(ms1 = ms1_df, ms2 = ms2_df, info = info)
}

## ---- Binary data decoder ----------------------------------------------------
# mzML binary data arrays are base64-encoded and, per MS:1000574, zlib-
# compressed (RFC 1950 -- distinct from gzip/RFC 1952). Decoded natively in R:
# base64 via xfun, decompression via memDecompress() (its "unknown" type
# auto-detects and correctly inflates zlib streams, unlike gzcon() which
# expects a gzip header and silently passes zlib data through unchanged).
.decode_binary <- function(b64, dtype = "64", is_zlib = TRUE) {
  raw <- xfun::base64_decode(b64)
  decomp <- if (is_zlib) memDecompress(raw, type = "unknown") else raw
  if (dtype == "32") {
    readBin(decomp, "numeric", n = length(decomp) / 4, size = 4)
  } else {
    readBin(decomp, "numeric", n = length(decomp) / 8, size = 8)
  }
}

## ---- Peak list import (text/CSV) ------------------------------------------
#' Import a simple peak-list file (mz, intensity columns)
#'
#' For a caller with a plain text/CSV export instead of raw mzML/mzXML
#' data -- e.g. a peak-picked table from vendor software. Column names
#' are normalized (lowercased, spaces/dots to underscores) and missing
#' optional columns are filled with `NA`/`1` as appropriate, so the
#' minimum required input is just an `mz` column.
#'
#' @param file Path to the peak-list file.
#' @param type `"ms1"` or `"ms2"`; `"ms2"` additionally ensures
#'   `precursor_mz`/`precursor_z` columns exist.
#' @param sep Field separator passed to `utils::read.table()`.
#' @return A data.frame: `rt`, `mz`, `intensity` for `type = "ms1"`
#'   (`rt` is `NA` if not present in the file); adds `precursor_mz`,
#'   `precursor_z` for `type = "ms2"`.
#' @seealso [extract_ms1_features()], [annotate_metabolites()]
#' @export
import_peak_list <- function(file, type = c("ms1", "ms2"), sep = ",") {
  type <- match.arg(type)
  df <- utils::read.table(file, header = TRUE, sep = sep,
                          stringsAsFactors = FALSE, check.names = FALSE)
  # Normalize column names
  names(df) <- tolower(gsub("[. ]", "_", names(df)))
  need_mz <- "mz" %in% names(df)
  need_int <- "intensity" %in% names(df) || "int" %in% names(df)
  if (!need_mz) stop("Peak list must have an 'mz' column")
  if (!need_int) df$intensity <- df$int %||% 1
  if (!("rt" %in% names(df))) df$rt <- NA_real_
  if (type == "ms2") {
    if (!("precursor_mz" %in% names(df))) df$precursor_mz <- NA_real_
    if (!("precursor_z" %in% names(df))) df$precursor_z <- NA_integer_
  }
  df
}

## ---- MS1 feature extraction ------------------------------------------------
# Extract features from raw MS1 peaks (group nearby m/z values that also
# co-elute in time). Two-stage grouping -- cluster by RT proximity first,
# then chain by ppm within each RT cluster -- mirroring the RT-then-mass
# clustering the Python batch pipeline already uses (_rt_clusters() +
# mass chaining in inst/python/oligomet_deconv/charge_group.py).
#
# Grouping by m/z alone (the previous behavior here) is RT-blind: any two
# peaks sharing an m/z within tolerance get merged into one "feature"
# regardless of how far apart in time they actually eluted, silently
# averaging together unrelated chromatographic events for any real run
# with more than a handful of scans. rt_tol is in minutes, matching
# MS:1000016 "scan start time"'s conventional unit (mzML files declare it
# via unitAccession UO:0000031, but this parser -- like the rest of this
# module -- doesn't convert units, so a file using different units would
# need a correspondingly different rt_tol).
# Median of ALL observed MS1 peak intensities -- most individual centroided
# peaks in a real LC-MS run are background/chemical noise rather than real
# analyte signal, so this is a simple, defensible proxy for the noise floor
# (not a windowed local-RMS baseline the way dedicated vendor software
# computes it, which needs per-region statistics this doesn't attempt).
# Mirrors _NoiseReservoir.noise_level() in the Python batch pipeline
# (inst/python/oligomet_deconv/roi.py), same idea applied to the R-native
# single-file path.
#' Estimate a file's MS1 background/noise level
#'
#' Median of every observed MS1 peak intensity. Most individual
#' centroided peaks in a real LC-MS run are background/chemical noise
#' rather than real analyte signal, so this is a simple, defensible proxy
#' for the noise floor -- not a windowed local-RMS baseline the way
#' dedicated vendor software computes it, which needs per-region
#' statistics this doesn't attempt.
#'
#' @param ms1_peaks A data.frame with an `intensity` column (e.g.
#'   `read_ms_file(file)$ms1`).
#' @return A single numeric noise level (the median of all positive,
#'   non-`NA` intensities), or `NA_real_` if there are none.
#' @seealso [extract_ms1_features()]'s `sn_threshold` argument, which
#'   uses this to compute `noise_level * sn_threshold` as the effective
#'   intensity floor.
#' @export
estimate_ms_noise_level <- function(ms1_peaks) {
  ints <- ms1_peaks$intensity
  ints <- ints[!is.na(ints) & ints > 0]
  if (length(ints) == 0) return(NA_real_)
  stats::median(ints)
}

# `sn_threshold`, when given, OVERRIDES min_intensity with a threshold
# derived from THIS file's own noise level: effective_threshold =
# estimate_ms_noise_level(ms1_peaks) * sn_threshold ("Absolute MS Signal
# Threshold = MS Noise Level x S/N Threshold"). A fixed min_intensity
# picked for one file's background level is rarely right for another file
# with a different one -- sn_threshold corrects for that, at the cost of
# the applied number no longer being a literal, predictable value up front.
#' Extract MS1 features from raw peaks
#'
#' Groups nearby raw MS1 (rt, mz, intensity) peaks into features:
#' clusters by retention-time proximity first, then chains by ppm
#' tolerance within each RT cluster (mirroring the RT-then-mass
#' clustering the Python batch pipeline uses -- see
#' `inst/python/oligomet_deconv/charge_group.py`), so peaks sharing an
#' m/z but eluting at unrelated times are never merged into one feature.
#'
#' @param ms1_peaks A data.frame with `rt`, `mz`, `intensity` columns
#'   (e.g. `read_ms_file(file)$ms1`).
#' @param ppm m/z chaining tolerance within an RT cluster.
#' @param min_intensity Fixed absolute intensity floor; ignored (falls
#'   back to it only if the noise estimate is `NA`) when `sn_threshold`
#'   is given.
#' @param rt_tol Retention-time clustering tolerance, in the same units
#'   `ms1_peaks$rt` uses (conventionally minutes; this parser doesn't
#'   convert units).
#' @param sn_threshold When given, OVERRIDES `min_intensity` with a
#'   threshold derived from THIS file's own noise level:
#'   `effective_min_intensity = estimate_ms_noise_level(ms1_peaks) *
#'   sn_threshold` ("Absolute MS Signal Threshold = MS Noise Level x S/N
#'   Threshold"). Corrects for files/instruments with different
#'   background levels, at the cost of the applied number no longer
#'   being a fixed, predictable value picked up front.
#' @return A data.frame of extracted features (columns include `mz`,
#'   `rt`, `max_intensity`, `n_scans`), one row per chained feature.
#'   Empty data.frame if `ms1_peaks` is empty or nothing clears the
#'   threshold.
#' @seealso [estimate_ms_noise_level()], [read_ms_file()]
#' @export
extract_ms1_features <- function(ms1_peaks, ppm = 10, min_intensity = 100,
                                  rt_tol = 0.15, sn_threshold = NULL) {
  if (nrow(ms1_peaks) == 0) return(data.frame())
  effective_min_intensity <- min_intensity
  if (!is.null(sn_threshold)) {
    noise <- estimate_ms_noise_level(ms1_peaks)
    if (!is.na(noise)) effective_min_intensity <- noise * sn_threshold
  }
  df <- ms1_peaks[ms1_peaks$intensity >= effective_min_intensity, ]
  if (nrow(df) == 0) return(data.frame())

  # Stage 1: cluster by RT proximity (adjacent-gap chaining).
  df <- df[order(df$rt), ]
  rt_groups <- integer(nrow(df))
  cur_rt <- 1L
  rt_groups[1] <- cur_rt
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      gap <- df$rt[i] - df$rt[i - 1]
      if (is.na(gap) || gap > rt_tol) cur_rt <- cur_rt + 1L
      rt_groups[i] <- cur_rt
    }
  }

  # Stage 2: within each RT cluster, chain by m/z proximity (ppm).
  out <- lapply(split(df, rt_groups), function(rt_grp) {
    rt_grp <- rt_grp[order(rt_grp$mz), ]
    n <- nrow(rt_grp)
    mz_groups <- integer(n)
    cur_mz <- 1L
    mz_groups[1] <- cur_mz
    if (n > 1) {
      for (i in 2:n) {
        tol <- rt_grp$mz[i] * ppm / 1e6
        if (rt_grp$mz[i] - rt_grp$mz[i - 1] > tol) cur_mz <- cur_mz + 1L
        mz_groups[i] <- cur_mz
      }
    }
    do.call(rbind, lapply(split(rt_grp, mz_groups), function(g) {
      data.frame(
        mz = weighted.mean(g$mz, g$intensity),
        rt = mean(g$rt, na.rm = TRUE),
        max_intensity = max(g$intensity),
        n_scans = nrow(g),
        stringsAsFactors = FALSE)
    }))
  })
  do.call(rbind, out)
}

## ---- Targeted MS1 matching ------------------------------------------------
#' Match a metabolite library against observed MS1 features
#'
#' For each metabolite, enumerates every (oxidation level x charge state
#' x adduct) theoretical target and reports the nearest observed feature
#' within `ppm_tol`, along with an isotope-pattern fit score. Nearest-
#' feature matching reports the single closest target as if it were
#' unambiguous; `n_candidates`/`ambiguous` flag when more than one
#' DISTINCT metabolite's theoretical m/z falls within tolerance of the
#' same observed feature, which is common at the tolerances this package
#' targets and otherwise invisible in the results table.
#'
#' @param mets A list of metabolite objects (see [generate_metabolites()]).
#' @param ms1_features Observed MS1 features (see
#'   [extract_ms1_features()]): a data.frame with `mz`, `rt`,
#'   `max_intensity`, and optionally `area`.
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @param ppm_tol Mass tolerance for a match.
#' @param z_range Charge states to consider.
#' @param adducts Adducts to consider (`"H"` is the unmodified `[M-zH]^z-`
#'   envelope; others add `adduct_shift()`).
#' @param max_oxid Max PS -> PO oxidation events to model per metabolite
#'   (capped at that metabolite's own phosphorothioate count).
#' @param h_offset Charge-envelope offset: `0` for the standard
#'   `[M-zH]^z-` convention, non-zero (e.g. `3.0046`) only to match a
#'   legacy workbook's convention.
#' @param n_iso Number of isotope peaks to check for the fit score; `0`
#'   skips isotope fitting entirely (`iso_fit` stays `NA`).
#' @param use_envipat Whether to use enviPat for isotope patterns
#'   (`FALSE` uses a built-in convolution instead).
#' @return A data.frame, one row per (metabolite, oxidation, charge,
#'   adduct) combination that matched within `ppm_tol`: `met_id`,
#'   `met_name`, `kind`, `n`, `k_oxid`, `z`, `adduct`, `theo_mz`,
#'   `obs_mz`, `ppm_error`, `rt`, `intensity`, `area` (`NA` unless
#'   `ms1_features` carries it), `iso_fit`, `formula`, `n_candidates`,
#'   `ambiguous`. Empty data.frame if nothing matches.
#' @seealso [match_ms1_fast()] for a faster path when the caller already
#'   has a pre-expanded target list, [match_ms1_batch()] for multi-sample
#'   matching, [annotate_metabolites()]
#' @export
match_ms1 <- function(mets, ms1_features, dict = STANDARD_DICT,
                       ppm_tol = 10, z_range = 3:12,
                       adducts = c("H", "Na", "K", "NH4"),
                       max_oxid = 6, h_offset = 0,
                       n_iso = 5, use_envipat = TRUE) {
  if (nrow(ms1_features) == 0) return(data.frame())

  matches <- list()
  for (met in mets) {
    info <- metabolite_mass_info(met, dict)
    kmax <- min(met$n_ps, max_oxid)

    for (k in 0:kmax) {
      mass <- ps_oxid_mass(info$mono_mass, k)
      fv <- ps_oxid_formula(info$formula_vec, k)

      for (z in z_range) {
        for (ad in adducts) {
          # Calculate theoretical m/z for this adduct
          ad_shift <- if (ad == "H") 0 else adduct_shift(ad)
          theo_mass <- mass + ad_shift
          theo_mz <- (theo_mass + h_offset - z * .PROTON) / z
          if (theo_mz < 100 || theo_mz > 3000) next

          # Find matching features
          dmz <- abs(ms1_features$mz - theo_mz)
          best <- which.min(dmz)
          ppm_err <- dmz[best] / theo_mz * 1e6

          if (ppm_err <= ppm_tol) {
            # Isotope pattern fit
            iso_fit <- NA
            if (n_iso > 0) {
              theo_iso <- isotope_mz_cluster(fv, z, n_top = n_iso,
                                              h_offset = h_offset,
                                              use_envipat = use_envipat)
              if (!is.null(theo_iso) && nrow(theo_iso) > 1) {
                iso_fit <- .isotope_fit(theo_iso, ms1_features, ppm_tol)
              }
            }

            matches[[length(matches) + 1]] <- data.frame(
              met_id = met$id, met_name = met$name, kind = met$kind,
              n = met$n, k_oxid = k, z = z, adduct = ad,
              theo_mz = round(theo_mz, 4),
              obs_mz = round(ms1_features$mz[best], 4),
              ppm_error = round(ppm_err, 3),
              rt = ms1_features$rt[best],
              intensity = ms1_features$max_intensity[best],
              # Peak area (trapezoidal AUC) is only ever present when
              # ms1_features came from the batch/Python ROI pipeline (see
              # match_ms1_batch()) -- always emit the column, NA when not,
              # so every caller (rbind across metabolites/samples, Excel
              # export, degradation_summary()) sees a stable column set.
              area = if ("area" %in% names(ms1_features)) ms1_features$area[best] else NA_real_,
              iso_fit = round(iso_fit, 3),
              formula = format_formula(fv),
              .feat_idx = best,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  if (length(matches) == 0) return(data.frame())
  out <- do.call(rbind, matches)
  # Ambiguity flag: at the tolerances this function defaults to, most
  # observed features sit within ppm_tol of more than one theoretical
  # target (metabolite x oxidation x charge x adduct) -- which.min()-style
  # nearest-match logic above silently reports only the single closest
  # one as if it were an unambiguous identification. n_candidates counts
  # how many DISTINCT library metabolites (not oxidation/charge/adduct
  # variants of the same metabolite) picked this exact feature as their
  # best match; ambiguous flags any feature more than one metabolite
  # could explain, so a caller can require MS2/RT confirmation before
  # trusting a match that would otherwise look identical in this table.
  n_candidates <- tapply(out$met_id, out$.feat_idx, function(x) length(unique(x)))
  out$n_candidates <- as.integer(n_candidates[as.character(out$.feat_idx)])
  out$ambiguous <- out$n_candidates > 1L
  out$.feat_idx <- NULL
  out
}

## ---- Vectorised MS1 matching over a pre-expanded target list ---------------
# match_ms1() above enumerates metabolite x oxidation x charge x adduct with
# an R loop and does a full which.min() scan over every feature inside it --
# for a 55-metabolite library at the app's Standard-mode defaults that is
# 55 x 7 x 10 x 4 = 15,400 full-vector scans per file. This is an
# alternative entry point for callers that already have (or can cheaply
# build) the expanded theoretical target list as a flat data.frame: it
# replaces the nested scan with one sorted merge (features sorted once,
# findInterval() locates each target's tolerance window in it), so cost
# scales with n_targets*log(n_features) instead of n_targets*n_features.
# It does not compute isotope fit or adduct/oxidation bookkeeping itself
# (the caller's `targets` frame already carries kind/k_oxid/z/adduct) --
# match_ms1() remains the entry point when those are needed.
#
#' Match a pre-expanded target list against MS1 features (fast path)
#'
#' A faster alternative to [match_ms1()] for a caller that already has
#' (or can cheaply build) the expanded theoretical target list as a flat
#' data.frame: sorts `features` once and uses `findInterval()` to locate
#' each target's tolerance window, so cost scales with
#' `n_targets * log(n_features)` instead of `match_ms1()`'s
#' `n_targets * n_features`. Does NOT compute isotope fit or adduct/
#' oxidation bookkeeping itself -- the caller's `targets` frame is
#' assumed to already carry `kind`/`k_oxid`/`z`/`adduct`; use
#' [match_ms1()] when those need to be derived.
#'
#' @param targets A data.frame with (at least) `met_id`, `met_name`,
#'   `kind`, `k_oxid`, `z`, `adduct`, `theo_mz` -- one row per
#'   theoretical target.
#' @param features Observed MS1 features (see [extract_ms1_features()]):
#'   a data.frame with `mz`, `rt`, `max_intensity`.
#' @param ppm_tol Mass tolerance for a match.
#' @return `targets`, filtered to rows with a match within `ppm_tol`,
#'   with `obs_mz`, `rt`, `intensity`, `ppm_error` appended, plus the
#'   same `n_candidates`/`ambiguous` ambiguity flag as [match_ms1()]
#'   (computed by distinct `met_id` competing for the same observed
#'   feature).
#' @seealso [match_ms1()]
#' @export
match_ms1_fast <- function(targets, features, ppm_tol = 10) {
  stopifnot("theo_mz" %in% names(targets), "mz" %in% names(features))
  if (nrow(targets) == 0 || nrow(features) == 0) return(targets[0, ])
  f <- features[order(features$mz), ]
  lo <- targets$theo_mz * (1 - ppm_tol / 1e6)
  hi <- targets$theo_mz * (1 + ppm_tol / 1e6)
  i_lo <- findInterval(lo, f$mz) + 1L
  i_hi <- findInterval(hi, f$mz)
  hit <- i_hi >= i_lo
  if (!any(hit)) return(targets[0, ])
  idx <- mapply(function(a, b, t) {
    j <- a:b
    j[which.min(abs(f$mz[j] - t))]
  }, i_lo[hit], i_hi[hit], targets$theo_mz[hit])
  out <- targets[hit, ]
  out$obs_mz <- f$mz[idx]
  out$rt <- f$rt[idx]
  out$intensity <- f$max_intensity[idx]
  out$ppm_error <- (out$obs_mz - out$theo_mz) / out$theo_mz * 1e6
  # Ambiguity: distinct metabolites competing for the same observed feature.
  n_candidates <- tapply(out$met_id, idx, function(x) length(unique(x)))
  out$n_candidates <- as.integer(n_candidates[as.character(idx)])
  out$ambiguous <- out$n_candidates > 1L
  out
}

## ---- Isotope pattern fit --------------------------------------------------
# Compare theoretical isotope pattern to observed features.
# Uses cosine similarity of relative abundances.
.isotope_fit <- function(theo_iso, ms1_features, ppm_tol) {
  if (is.null(theo_iso) || nrow(theo_iso) < 2) return(NA)
  # For each theoretical isotope peak, find the nearest observed feature
  obs_intensities <- numeric(nrow(theo_iso))
  for (i in seq_len(nrow(theo_iso))) {
    dmz <- abs(ms1_features$mz - theo_iso$mz[i])
    best <- which.min(dmz)
    ppm_err <- dmz[best] / theo_iso$mz[i] * 1e6
    obs_intensities[i] <- if (ppm_err <= ppm_tol) ms1_features$max_intensity[best] else 0
  }
  # Cosine similarity between theoretical and observed abundance vectors
  theo_abund <- theo_iso$abundance
  obs_abund <- obs_intensities
  if (sum(obs_abund) == 0) return(0)
  cos_sim <- sum(theo_abund * obs_abund) /
    (sqrt(sum(theo_abund^2)) * sqrt(sum(obs_abund^2)))
  # Penalize partial matches: cosine similarity alone can stay high even
  # when most theoretical isotope peaks went unmatched (e.g. 2/5 peaks
  # matched with correct ratios can still score ~0.86), so scale by the
  # fraction of theoretical peaks that were actually matched.
  cos_sim <- cos_sim * (sum(obs_abund > 0) / length(obs_abund))
  cos_sim
}

## ---- Charge envelope consistency ------------------------------------------
# Check if matched charge states form a consistent envelope.
# Consistency = the observed m/z values across charge states should all
# correspond to the same neutral mass within tolerance.
#
# h_offset must match whatever was used to generate the theoretical m/z
# these matches were scored against (0 = standard [M-zH]^z-, 3.0046 =
# legacy workbook offset -- see charge_envelope() in R/mass_isotope.R).
# It was hardcoded to 0 here previously; with a non-zero h_offset every
# back-calculated mass was off by z * h_offset (a different amount per
# charge state), which corrupts the cross-charge-state agreement this
# function exists to measure.
#' Check whether matched charge states agree on one neutral mass
#'
#' Groups [match_ms1()]'s output by (metabolite, oxidation, adduct) and
#' back-calculates the neutral mass implied by each matched charge
#' state's observed m/z. A real multiply-charged species should
#' back-calculate to the SAME neutral mass at every charge state it was
#' observed at; this reports how tightly they agree.
#'
#' @param matches Output of [match_ms1()]/[match_ms1_batch()] (must have
#'   `met_id`, `k_oxid`, `adduct`, `z`, `obs_mz`).
#' @param h_offset Must match whatever `h_offset` [match_ms1()] used to
#'   generate the theoretical m/z these matches were scored against --
#'   see [match_ms1()]. A mismatch here silently corrupts the agreement
#'   this function measures, since each charge state's back-calculated
#'   mass would be off by a DIFFERENT amount (`z * h_offset`).
#' @return A data.frame, one row per (metabolite, oxidation, adduct)
#'   group: `group` (the grouping key), `n_z` (number of charge states
#'   observed), `consistency` (0-1 scale, `NA` when `n_z < 2` since a
#'   single charge state can't confirm consistency), `mass_cv_ppm`
#'   (coefficient of variation of the back-calculated masses, in ppm),
#'   `mean_mass`.
#' @seealso [match_ms1()]
#' @export
envelope_consistency <- function(matches, h_offset = 0) {
  if (nrow(matches) == 0) return(data.frame())
  # Group by metabolite (met_id + k_oxid + adduct)
  matches$group <- paste(matches$met_id, matches$k_oxid, matches$adduct, sep = "_")
  # Every branch must return the SAME columns -- these frames get rbind()ed
  # together below (and again by callers grouping across samples), and a
  # data.frame with fewer/differently-named columns silently made rbind()
  # error out as soon as a real match set had a mix of single- and multi-
  # charge-state groups (previously untested: the n_z<2 branch returned
  # `mass_cv`/no `mean_mass`, the n_z>=2 branch returned `mass_cv_ppm` and
  # `mean_mass`).
  do.call(rbind, lapply(split(matches, matches$group), function(g) {
    # Back-calculate neutral mass from each charge state
    # M = z * mz + z * proton - h_offset - adduct_shift
    ad_shift <- if (g$adduct[1] == "H") 0 else adduct_shift(g$adduct[1])
    masses <- g$z * g$obs_mz + g$z * .PROTON - h_offset - ad_shift
    if (nrow(g) < 2) {
      # a single charge state can't confirm envelope consistency, but its
      # own back-calculated mass is still a useful point estimate
      return(data.frame(group = g$group[1], n_z = nrow(g),
                        consistency = NA_real_, mass_cv_ppm = NA_real_,
                        mean_mass = round(masses[1], 4),
                        stringsAsFactors = FALSE))
    }
    cv <- sd(masses) / mean(masses) * 1e6  # coefficient of variation in ppm
    data.frame(group = g$group[1], n_z = nrow(g),
               consistency = 1 / (1 + cv / 100),  # 0-1 scale
               mass_cv_ppm = round(cv, 2),
               mean_mass = round(mean(masses), 4),
               stringsAsFactors = FALSE)
  }))
}

## ---- MS2 spectrum lookup --------------------------------------------------
#' Find MS2 spectra matching a given precursor m/z and charge
#'
#' Groups `ms2_data` into scans (by unique `rt`/`precursor_mz`) and
#' returns the peaks from every scan whose precursor matches within
#' tolerance.
#'
#' @param ms2_data A data.frame of MS2 peaks with `rt`, `precursor_mz`,
#'   `precursor_z`, `mz`, `intensity` (see [read_ms_file()]/
#'   [import_peak_list()]/[read_batch_ms2()]).
#' @param precursor_mz Target precursor m/z to match against.
#' @param precursor_z Target precursor charge; `NA` (the default) skips
#'   the charge check (matches on m/z alone). When both the target and a
#'   scan's own recorded charge are non-`NA`, they must agree exactly.
#' @param ppm_tol Precursor m/z matching tolerance.
#' @return A list of data.frames (`mz`, `intensity`), one per matching
#'   MS2 scan. Empty list if `ms2_data` is empty or nothing matches.
#' @seealso [confirm_metabolite()], [annotate_metabolites()]
#' @export
find_ms2_spectra <- function(ms2_data, precursor_mz, precursor_z = NA,
                              ppm_tol = 20) {
  if (nrow(ms2_data) == 0) return(list())
  # Group MS2 data by scan (unique rt + precursor_mz)
  scans <- unique(ms2_data[, c("rt", "precursor_mz", "precursor_z")])
  matched <- list()
  for (i in seq_len(nrow(scans))) {
    if (is.na(scans$precursor_mz[i])) next
    dmz <- abs(scans$precursor_mz[i] - precursor_mz)
    ppm_err <- dmz / precursor_mz * 1e6
    if (ppm_err <= ppm_tol) {
      if (!is.na(precursor_z) && !is.na(scans$precursor_z[i])) {
        if (scans$precursor_z[i] != precursor_z) next
      }
      idx <- ms2_data$rt == scans$rt[i] &
        ms2_data$precursor_mz == scans$precursor_mz[i]
      matched[[length(matched) + 1]] <- ms2_data[idx, c("mz", "intensity")]
    }
  }
  matched
}

## ---- Combined metabolite annotation ----------------------------------------
#' Run the single-file annotation pipeline: MS1 match + MS2 confirmation
#'
#' The single-file counterpart to [annotate_metabolites_batch()]: MS1
#' matching ([match_ms1()]), charge-envelope consistency
#' ([envelope_consistency()]), and -- when `ms2_data` is supplied -- MS2
#' fragment confirmation ([find_ms2_spectra()], [match_fragments()],
#' [confirmation_score()]) for every matched metabolite, rolled up into
#' one summary table.
#'
#' @param mets A list of metabolite objects (see [generate_metabolites()]).
#' @param ms1_features Observed MS1 features (see
#'   [extract_ms1_features()]/[import_peak_list()]).
#' @param ms2_data Optional MS2 peaks (see [read_ms_file()]/
#'   [import_peak_list()]); `NULL` skips MS2 confirmation entirely.
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @param ppm_tol,z_range,adducts,max_oxid,h_offset,n_iso,use_envipat
#'   MS1 matching parameters -- see [match_ms1()].
#' @param frag_tol_ppm,frag_z_range MS2 fragment matching tolerance/
#'   charge range -- see [match_fragments()].
#' @param include_internal Whether to also generate and match internal
#'   fragments (see [generate_internal_fragments()]).
#' @return A list: `ms1_matches` (see [match_ms1()]), `ms2_results`
#'   (named list by `met_id`, each with `met_name`, `n_peaks`,
#'   `best_spec`, `matched_frags`, `diagnostics`, `score`), and
#'   `summary` (one row per matched metabolite: best ppm error, best
#'   isotope fit, charge states/adducts found, and -- when MS2 was run
#'   -- confirmation score/coverage/fragment count).
#' @seealso [annotate_metabolites_batch()] for the multi-sample version
#' @export
annotate_metabolites <- function(mets, ms1_features, ms2_data = NULL,
                                  dict = STANDARD_DICT,
                                  ppm_tol = 10, z_range = 3:12,
                                  adducts = c("H", "Na", "K", "NH4"),
                                  max_oxid = 6, h_offset = 0,
                                  frag_tol_ppm = 25, frag_z_range = 1:2,
                                  n_iso = 5, use_envipat = TRUE,
                                  include_internal = FALSE) {
  # MS1 matching
  ms1_matches <- match_ms1(mets, ms1_features, dict, ppm_tol, z_range,
                            adducts, max_oxid, h_offset, n_iso = n_iso,
                            use_envipat = use_envipat)
  if (nrow(ms1_matches) == 0) {
    return(list(ms1_matches = data.frame(), ms2_results = list(),
                summary = data.frame()))
  }

  # Envelope consistency
  env_cons <- envelope_consistency(ms1_matches, h_offset = h_offset)

  # MS2 fragment confirmation (if MS2 data provided)
  ms2_results <- list()
  if (!is.null(ms2_data) && nrow(ms2_data) > 0) {
    # For each unique matched metabolite, find MS2 spectra and confirm
    unique_mets <- unique(ms1_matches[, c("met_id", "met_name", "theo_mz", "z")])
    for (i in seq_len(nrow(unique_mets))) {
      met_id <- unique_mets$met_id[i]
      met <- mets[[which(sapply(mets, function(m) m$id == met_id))]]
      spectra <- find_ms2_spectra(ms2_data, unique_mets$theo_mz[i],
                                   unique_mets$z[i], ppm_tol * 2)
      if (length(spectra) > 0) {
        # Use the spectrum with most peaks
        best_spec <- spectra[[which.max(sapply(spectra, nrow))]]
        # Generate fragments and match
        frags <- generate_fragments(met, dict, z_range = frag_z_range)
        if (include_internal) {
          frags <- c(frags, generate_internal_fragments(met, dict, z_range = frag_z_range))
        }
        # Reuses the same adducts requested for MS1 matching (a superset
        # search, not narrowed to the specific adduct this metabolite's MS1
        # hit matched -- unique_mets above doesn't carry per-hit adduct
        # identity). confirm_ms2_batch() in batch_ms_processing.R does the
        # more precise per-hit version since match_ms1_batch() does track it.
        matched_frags <- match_fragments(frags, best_spec, frag_tol_ppm, frag_z_range,
                                          adducts = adducts)
        diags <- check_ps_diagnostic(best_spec, max(frag_tol_ppm, 50))
        score <- confirmation_score(matched_frags, met$n, diags)
        ms2_results[[met_id]] <- list(
          met_name = unique_mets$met_name[i],
          n_peaks = nrow(best_spec),
          best_spec = best_spec,
          matched_frags = matched_frags,
          diagnostics = diags,
          score = score
        )
      }
    }
  }

  # Summary table
  summary <- do.call(rbind, lapply(unique(ms1_matches$met_id), function(mid) {
    sub <- ms1_matches[ms1_matches$met_id == mid, ]
    met <- mets[[which(sapply(mets, function(m) m$id == mid))]]
    ms2 <- ms2_results[[mid]]
    data.frame(
      met_id = mid, met_name = sub$met_name[1], kind = sub$kind[1],
      n = sub$n[1], n_ms1_matches = nrow(sub),
      best_ppm = round(min(sub$ppm_error), 3),
      best_iso_fit = {
        v <- sub$iso_fit[!is.na(sub$iso_fit)]
        if (length(v) > 0) round(max(v), 3) else NA
      },
      n_charge_states = length(unique(sub$z)),
      adducts_found = paste(unique(sub$adduct), collapse = ","),
      has_ms2 = !is.null(ms2),
      ms2_score = if (!is.null(ms2)) ms2$score$total_score else NA,
      ms2_coverage = if (!is.null(ms2)) ms2$score$coverage else NA,
      n_ms2_frags = if (!is.null(ms2)) ms2$score$n_matches else NA,
      n_diagnostics = if (!is.null(ms2)) ms2$score$n_diagnostics else NA,
      confident = if (!is.null(ms2)) ms2$score$confident else FALSE,
      stringsAsFactors = FALSE
    )
  }))

  list(ms1_matches = ms1_matches, envelope = env_cons,
       ms2_results = ms2_results, summary = summary)
}
