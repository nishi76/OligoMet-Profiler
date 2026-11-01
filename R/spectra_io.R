# =============================================================================
# spectra_io.R
# Robust mzML/mzXML reading via Bioconductor's Spectra + MsBackendMzR
# (rformassspectrometry.org), with the hand-rolled xml2 parser in
# ms_matching.R (parse_mzml()) kept as an automatic fallback.
#
# Why: parse_mzml() loads the whole file into an in-memory xml2 DOM (no
# streaming) and materializes MS1 peaks via a rbind()-in-a-loop that scales
# poorly on large/profile-mode files -- it also has no numpress support and
# no indexed-mzML fast seeking. Spectra::MsBackendMzR() wraps mzR's compiled
# backend, which handles all of that (as of recent Spectra releases the
# MsBackendMzR class ships inside the Spectra package itself, not a separate
# "MsBackendMzR" package -- mzR remains the actual compiled-code dependency).
# Spectra/mzR are optional (Suggests, not Imports) -- mzR pulls in
# Bioconductor and a compiled, platform-sensitive dependency, so a user who
# can't get it installed must still be able to run the app; read_ms_file()
# silently falls back to parse_mzml() whenever Spectra isn't available or
# fails to read a file.
#
# Batch (multi-file) mode keeps reading via the Python/pyteomics pipeline
# (inst/python/oligomet_deconv/io_mzml.py) unchanged -- Spectra is not
# threaded through that subprocess boundary. See run_batch_deconvolution()
# in R/batch_ms_processing.R for where that boundary is.
# =============================================================================

## ---- Availability check ----------------------------------------------------
# Spectra::MsBackendMzR() needs the compiled mzR package at runtime even
# though it's part of the Spectra package's own R code -- check both.
.have_spectra <- function() {
  requireNamespace("Spectra", quietly = TRUE) &&
    requireNamespace("mzR", quietly = TRUE)
}

## ---- Primary entry point ----------------------------------------------------
# Drop-in replacement for parse_mzml()'s caller-facing contract: same return
# shape -- list(ms1 = data.frame(rt, mz, intensity), ms2 = data.frame(rt,
# precursor_mz, precursor_z, mz, intensity), info = list(file, n_spectra,
# n_ms1, n_ms2, instrument, profile_mode)) -- so call sites need a one-line
# swap (parse_mzml(path) -> read_ms_file(path)), not a rewrite.
read_ms_file <- function(file, prefer_spectra = TRUE) {
  if (prefer_spectra && .have_spectra()) {
    result <- tryCatch(.read_ms_file_spectra(file), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    message("Spectra/mzR read failed (", conditionMessage(result),
            "); falling back to the built-in mzML/mzXML parser.")
  }
  parse_mzml(file)
}

## ---- Spectra/MsBackendMzR reader --------------------------------------------
.read_ms_file_spectra <- function(file) {
  sp <- Spectra::Spectra(file, source = Spectra::MsBackendMzR())

  levels <- Spectra::msLevel(sp)

  # MS1
  sp1 <- sp[levels == 1L]
  ms1_df <- if (length(sp1) == 0) {
    data.frame(rt = numeric(0), mz = numeric(0), intensity = numeric(0))
  } else {
    .stack_peaks(Spectra::peaksData(sp1), Spectra::rtime(sp1))
  }

  # MS2
  sp2 <- sp[levels == 2L]
  ms2_df <- if (length(sp2) == 0) {
    data.frame(rt = numeric(0), precursor_mz = numeric(0),
               precursor_z = integer(0), mz = numeric(0), intensity = numeric(0))
  } else {
    .stack_ms2_peaks(Spectra::peaksData(sp2), Spectra::rtime(sp2),
                      Spectra::precursorMz(sp2), Spectra::precursorCharge(sp2))
  }

  # Profile-vs-centroid: mirror parse_mzml()'s convention of deciding this
  # from the first MS1 spectrum only. Spectra::centroided() is a per-spectrum
  # logical (NA if the file didn't declare it) -- profile_mode is the
  # inverse, consistent with parse_mzml()'s info$profile_mode semantics
  # (TRUE = fires the "appears to be PROFILE mode" warning in inst/app/app.R).
  profile_mode <- if (length(sp1) == 0) NA else {
    c1 <- Spectra::centroided(sp1)[1]
    if (is.na(c1)) NA else !isTRUE(c1)
  }

  list(
    ms1 = ms1_df, ms2 = ms2_df,
    info = list(file = file, n_spectra = length(sp),
                n_ms1 = length(sp1), n_ms2 = length(sp2),
                instrument = "unknown",  # not exposed uniformly across mzR backends
                profile_mode = profile_mode)
  )
}

# peaksData() returns a list of 2-column (mz, intensity) numeric matrices,
# one per spectrum -- stack them into one long data.frame with the matching
# per-spectrum retention time attached to every peak row.
.stack_peaks <- function(peaks_list, rt) {
  rows <- lapply(seq_along(peaks_list), function(i) {
    m <- peaks_list[[i]]
    if (nrow(m) == 0) return(NULL)
    data.frame(rt = rt[i], mz = m[, 1], intensity = m[, 2], stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(data.frame(rt = numeric(0), mz = numeric(0), intensity = numeric(0)))
  do.call(rbind, rows)
}

.stack_ms2_peaks <- function(peaks_list, rt, precursor_mz, precursor_z) {
  rows <- lapply(seq_along(peaks_list), function(i) {
    m <- peaks_list[[i]]
    if (nrow(m) == 0) return(NULL)
    data.frame(rt = rt[i], precursor_mz = precursor_mz[i], precursor_z = precursor_z[i],
               mz = m[, 1], intensity = m[, 2], stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(data.frame(rt = numeric(0), precursor_mz = numeric(0),
                       precursor_z = integer(0), mz = numeric(0), intensity = numeric(0)))
  }
  do.call(rbind, rows)
}

## ---- Vendor-raw pre-conversion gate ------------------------------------------
# Given an uploaded/local file path, return a path read_ms_file() (or the
# Python batch pipeline) can consume directly: passes .mzML/.mzXML through
# unchanged, converts vendor formats via the ProteoWizard msconvert bridge
# (find_msconvert()/convert_vendor_raw(), R/ms_matching.R) since neither
# reader here nor pyteomics on the Python side understands vendor formats
# natively. Standardizing on msconvert -> .mzML for both readers avoids
# needing two different vendor-handling code paths.
#
# .d (Agilent/Bruker) is a directory, not a single file, so it cannot arrive
# via a Shiny fileInput at all -- only reachable through the local-folder
# batch_dir path, where it's still routed through this same function.
resolve_ms_input_file <- function(path, output_dir = tempdir()) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("mzml", "mzxml")) return(path)
  if (ext %in% c("raw", "wiff", "baf", "yep", "d")) {
    converted <- convert_vendor_raw(path, output_dir = output_dir, centroid = TRUE)
    if (is.na(converted)) {
      stop("Vendor file '", basename(path), "' requires ProteoWizard msconvert, ",
           "which was not found on PATH. Install ProteoWizard (proteowizard.org) ",
           "or convert externally: msconvert <file> --mzML --centroid")
    }
    return(converted)
  }
  path  # not a recognized MS format extension -- let the caller handle it
        # (e.g. import_peak_list() for .csv/.txt)
}
