# =============================================================================
# data_matrix.R
# Import/export of a pre-processed metabolite x sample data matrix, so every
# downstream step (statistics, calibration curve, back-calculated
# concentration, relative quantification, % degradation) can run on the
# matrix alone -- no peak picking, deconvolution, or MS1/MS2 matching.
#
# The matrix is converted into the same long shape annotate_metabolites_batch()
# produces for `ms1_matches` (one row per sample x metabolite, with `met_id`,
# `met_name`, `kind`, `sample`, `intensity`), so build_abundance_matrix(),
# quantify_metabolites(), degradation_summary(), apply_blank_correction(),
# run_pca()/run_hclust() all consume it unchanged.
#
# Accepted layouts (auto-detected, or forced with `format`):
#   wide -- one row per metabolite, one numeric column per sample, plus
#           optional `met_id`/`met_name`/`kind` identifier columns. This is
#           exactly what the app's "Data Matrix -- wide" download writes.
#   long -- one row per (metabolite, sample): a sample column, a metabolite
#           column, and a value column (intensity/area/signal/...). Any
#           group/timepoint/sample_type/concentration columns are carried
#           back out as sample information. This is exactly what the app's
#           "Data Matrix -- long" download writes.
# File types: .csv, .tsv, .txt (tab or comma), .xlsx (openxlsx), and .xls
# (only when the readxl package is installed).
# =============================================================================

.DM_ID_ALIASES <- c("met_id", "metabolite_id", "id", "feature_id", "feature",
                    "compound_id", "compound", "metabolite", "analyte")
.DM_NAME_ALIASES <- c("met_name", "metabolite_name", "name", "compound_name",
                      "analyte_name", "description")
.DM_KIND_ALIASES <- c("kind", "class", "metabolite_class", "type")
.DM_SAMPLE_ALIASES <- c("sample", "sample_name", "sample_id", "file", "filename",
                        "file_name", "raw_file", "injection")
.DM_VALUE_ALIASES <- c("intensity", "area", "signal", "abundance", "value",
                       "response", "peak_area", "peak_intensity", "height")
.DM_META_ALIASES <- list(
  group = c("group", "condition", "treatment", "arm"),
  timepoint = c("timepoint", "time", "time_point", "time_h", "hours", "day"),
  sample_type = c("sample_type", "sampletype"),
  concentration = c("concentration", "conc", "nominal_concentration",
                    "nominal_conc", "spiked_concentration", "amount")
)

## ---- Low-level readers -------------------------------------------------------

# Normalizes a header for alias matching only -- the original header text
# is what gets kept as the sample name.
.dm_norm_names <- function(x) {
  gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(x))))
}

.dm_first_match <- function(norm_names, aliases) {
  hit <- match(aliases, norm_names)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) NA_integer_ else hit[1]
}

# Strips a raw-file extension and surrounding whitespace, so "ctrl_1.mzML"
# in a sample information sheet still matches the "ctrl_1" column of a
# data matrix (and vice versa).
.dm_clean_sample <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.(mzml|mzxml|raw|wiff|d|baf|yep|csv|txt)$", "", x, ignore.case = TRUE)
}

#' Read a delimited text or Excel table
#'
#' Reads `.csv`, `.tsv`, `.txt` (tab- or comma-delimited, sniffed from the
#' first line), `.xlsx` (via openxlsx), or `.xls` (via readxl, when
#' installed). Column headers are kept verbatim (`check.names = FALSE`) so
#' sample names survive the round trip.
#'
#' @param path File path.
#' @param sheet Worksheet (name or 1-based index) for Excel files.
#' @param all_character If `TRUE`, every column is read as character.
#' @param name Optional original file name, used for the extension when
#'   `path` is a temporary upload path without one.
#' @return A data.frame.
#' @export
read_table_file <- function(path, sheet = 1, all_character = FALSE, name = NULL) {
  ext <- tolower(tools::file_ext(if (!is.null(name) && nzchar(name)) name else path))
  df <- if (ext %in% c("xlsx", "xlsm")) {
    openxlsx::read.xlsx(path, sheet = sheet, check.names = FALSE, sep.names = " ")
  } else if (ext == "xls") {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Reading .xls needs the 'readxl' package -- install.packages(\"readxl\"), ",
           "or save the sheet as .xlsx/.csv.")
    }
    as.data.frame(readxl::read_excel(path, sheet = sheet), check.names = FALSE)
  } else {
    first <- readLines(path, n = 1, warn = FALSE)
    sep <- if (ext == "tsv" || (length(first) > 0 && grepl("\t", first))) "\t"
           else if (length(first) > 0 && !grepl(",", first) && grepl(";", first)) ";"
           else ","
    utils::read.table(path, sep = sep, header = TRUE, quote = "\"", comment.char = "",
                      check.names = FALSE, stringsAsFactors = FALSE,
                      colClasses = if (all_character) "character" else NA,
                      na.strings = c("", "NA", "NaN", "N/A", "#N/A"),
                      fill = TRUE, strip.white = TRUE)
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- trimws(names(df))
  df <- df[, nzchar(names(df)) | vapply(df, function(v) any(!is.na(v) & nzchar(as.character(v))), logical(1)),
           drop = FALSE]
  if (all_character) df[] <- lapply(df, function(v) ifelse(is.na(v), "", as.character(v)))
  df
}

# as.numeric() that tolerates thousands separators and stray whitespace,
# e.g. "1,234,567" or " 5.2e6 ".
.dm_as_numeric <- function(v) {
  if (is.numeric(v)) return(as.numeric(v))
  v <- gsub("[ ,]", "", trimws(as.character(v)))
  suppressWarnings(as.numeric(v))
}

# A column counts as numeric (a sample column, in a wide matrix) when every
# non-empty value parses as a number.
.dm_is_numeric_col <- function(v) {
  if (is.numeric(v)) return(TRUE)
  s <- trimws(as.character(v))
  s <- s[!is.na(s) & nzchar(s) & !toupper(s) %in% c("NA", "NAN", "N/A", "#N/A", "ND")]
  length(s) > 0 && all(!is.na(.dm_as_numeric(s)))
}

## ---- Data matrix import --------------------------------------------------------

#' Read a pre-processed data matrix (wide or long) into match format
#'
#' Converts an externally processed (or previously exported) metabolite x
#' sample table into the long `ms1_matches` shape every downstream function
#' in this package consumes: one row per (sample, metabolite) with
#' `met_id`, `met_name`, `kind`, `sample`, and `intensity` (plus `area` when
#' a long file carries both). Empty, `NA`, and non-numeric cells are
#' treated as "not detected" and dropped.
#'
#' @param x A file path, or an already-loaded data.frame.
#' @param format `"auto"` (default), `"wide"`, or `"long"`. Auto picks long
#'   when the table has both a sample column and a value column
#'   (intensity/area/signal/...), else wide.
#' @param sheet Worksheet for Excel files.
#' @param name Original file name, when `x` is a temporary upload path.
#' @return A list: `matches` (the long data.frame), `format` (the layout
#'   actually used), `samples` (sample names, in file order), `sample_meta`
#'   (group/timepoint/sample_type/concentration found in a long file, else
#'   `NULL`), and `notes` (character vector of anything worth telling the
#'   user, e.g. dropped non-numeric cells).
#' @seealso [export_data_matrix()] for the inverse, [annotate_matrix_metabolites()]
#'   to fill in `met_name`/`kind` from a generated library.
#' @export
read_data_matrix <- function(x, format = c("auto", "wide", "long"), sheet = 1, name = NULL) {
  format <- match.arg(format)
  df <- if (is.data.frame(x)) x else read_table_file(x, sheet = sheet, name = name)
  if (nrow(df) == 0 || ncol(df) < 2) stop("the data matrix is empty or has fewer than 2 columns")
  nn <- .dm_norm_names(names(df))

  i_sample <- .dm_first_match(nn, .DM_SAMPLE_ALIASES)
  i_value <- .dm_first_match(nn, .DM_VALUE_ALIASES)
  if (format == "auto") format <- if (!is.na(i_sample) && !is.na(i_value)) "long" else "wide"

  notes <- character(0)
  sample_meta <- NULL

  if (format == "long") {
    if (is.na(i_sample)) stop("long format needs a sample column (one of: ",
                              paste(.DM_SAMPLE_ALIASES, collapse = ", "), ")")
    if (is.na(i_value)) stop("long format needs a value column (one of: ",
                             paste(.DM_VALUE_ALIASES, collapse = ", "), ")")
    i_id <- .dm_first_match(nn, .DM_ID_ALIASES)
    i_name <- .dm_first_match(nn, .DM_NAME_ALIASES)
    i_kind <- .dm_first_match(nn, .DM_KIND_ALIASES)
    if (is.na(i_id)) i_id <- i_name
    if (is.na(i_id)) stop("long format needs a metabolite column (one of: ",
                          paste(c(.DM_ID_ALIASES, .DM_NAME_ALIASES), collapse = ", "), ")")

    out <- data.frame(
      met_id = trimws(as.character(df[[i_id]])),
      met_name = if (!is.na(i_name)) trimws(as.character(df[[i_name]])) else NA_character_,
      kind = if (!is.na(i_kind)) trimws(as.character(df[[i_kind]])) else NA_character_,
      sample = .dm_clean_sample(df[[i_sample]]),
      intensity = .dm_as_numeric(df[[i_value]]),
      stringsAsFactors = FALSE)
    # Keep area alongside intensity when both are present, so the
    # area-vs-intensity auto-selection downstream still has a choice.
    i_area <- match("area", nn)
    if (!is.na(i_area) && i_area != i_value) out$area <- .dm_as_numeric(df[[i_area]])
    if (nn[i_value] == "area" && !is.na(match("intensity", nn))) {
      out$area <- out$intensity
      out$intensity <- .dm_as_numeric(df[[match("intensity", nn)]])
    }
    samples <- unique(out$sample)

    meta_idx <- vapply(.DM_META_ALIASES, function(a) .dm_first_match(nn, a), integer(1))
    meta_idx <- meta_idx[!is.na(meta_idx) & !(meta_idx %in% c(i_id, i_name, i_kind, i_value))]
    if (length(meta_idx) > 0) {
      sm <- data.frame(sample = .dm_clean_sample(df[[i_sample]]), stringsAsFactors = FALSE)
      for (col in names(meta_idx)) {
        v <- as.character(df[[meta_idx[[col]]]])
        sm[[col]] <- ifelse(is.na(v), "", trimws(v))
      }
      sample_meta <- sm[!duplicated(sm$sample), , drop = FALSE]
      rownames(sample_meta) <- NULL
    }
  } else {
    i_id <- .dm_first_match(nn, .DM_ID_ALIASES)
    i_name <- .dm_first_match(nn, .DM_NAME_ALIASES)
    i_kind <- .dm_first_match(nn, .DM_KIND_ALIASES)
    is_num <- vapply(df, .dm_is_numeric_col, logical(1))
    id_cols <- c(i_id, i_name, i_kind)
    id_cols <- id_cols[!is.na(id_cols)]
    if (is.na(i_id) && is.na(i_name)) {
      # No recognized identifier header: the first non-numeric column is
      # the metabolite identifier (e.g. an unnamed first column).
      first_text <- which(!is_num)[1]
      if (is.na(first_text)) stop("wide format needs a metabolite identifier column ",
                                  "(e.g. met_id) -- every column looks numeric")
      i_id <- first_text
      id_cols <- c(id_cols, i_id)
    }
    if (is.na(i_id)) i_id <- i_name
    sample_idx <- setdiff(which(is_num), id_cols)
    skipped <- setdiff(seq_along(df), c(sample_idx, id_cols))
    if (length(skipped) > 0) {
      notes <- c(notes, paste0("Ignored non-numeric column(s): ",
                               paste(names(df)[skipped], collapse = ", ")))
    }
    if (length(sample_idx) == 0) stop("wide format needs at least one numeric sample column")

    ids <- trimws(as.character(df[[i_id]]))
    nms <- if (!is.na(i_name)) trimws(as.character(df[[i_name]])) else rep(NA_character_, nrow(df))
    kinds <- if (!is.na(i_kind)) trimws(as.character(df[[i_kind]])) else rep(NA_character_, nrow(df))
    samples <- .dm_clean_sample(names(df)[sample_idx])
    out <- do.call(rbind, lapply(seq_along(sample_idx), function(j) {
      data.frame(met_id = ids, met_name = nms, kind = kinds, sample = samples[j],
                 intensity = .dm_as_numeric(df[[sample_idx[j]]]), stringsAsFactors = FALSE)
    }))
  }

  blank_id <- is.na(out$met_id) | !nzchar(out$met_id)
  if (any(blank_id)) {
    notes <- c(notes, sprintf("Dropped %d row(s) with no metabolite identifier.", sum(blank_id)))
    out <- out[!blank_id, , drop = FALSE]
  }
  out$met_name <- ifelse(is.na(out$met_name) | !nzchar(out$met_name), out$met_id, out$met_name)
  out$kind[!is.na(out$kind) & !nzchar(out$kind)] <- NA_character_

  sig_cols <- intersect(c("intensity", "area"), names(out))
  keep <- Reduce(`|`, lapply(sig_cols, function(cc) !is.na(out[[cc]])))
  out <- out[keep, , drop = FALSE]
  if (nrow(out) == 0) stop("no numeric signal values found in the data matrix")

  dup <- duplicated(out[, c("met_id", "sample")])
  if (any(dup)) {
    notes <- c(notes, sprintf(paste0("%d duplicate (metabolite, sample) row(s) found -- ",
                                     "kept all; downstream steps use the max per sample."), sum(dup)))
  }
  rownames(out) <- NULL
  list(matches = out, format = format, samples = samples,
       sample_meta = sample_meta, notes = notes)
}

#' Fill in metabolite names and classes for an imported data matrix
#'
#' Looks each `met_id` up in a generated metabolite library (from
#' [generate_metabolites()]; matched by id first, then by name) to fill in
#' `met_name` and `kind`. Anything still without a `kind` is inferred:
#' `"parent"` for an id/name of PARENT, intact, full-length, or FLP, else
#' `"unknown"` (counted as a degradant by [degradation_summary()]).
#'
#' @param matches `read_data_matrix()$matches`.
#' @param mets Optional metabolite library (list of lists with `id`,
#'   `name`, `kind`).
#' @param parent_met_id Optional metabolite id to force as the intact
#'   parent (`kind = "parent"`); any other row previously inferred as
#'   parent is demoted to `"unknown"` so the degradation ratio has exactly
#'   one parent.
#' @return `matches` with `met_name` and `kind` filled in.
#' @export
annotate_matrix_metabolites <- function(matches, mets = NULL, parent_met_id = NULL) {
  if (is.null(matches) || nrow(matches) == 0) return(matches)
  if (!"kind" %in% names(matches)) matches$kind <- NA_character_
  if (!"met_name" %in% names(matches)) matches$met_name <- matches$met_id

  if (!is.null(mets) && length(mets) > 0) {
    lib_id <- vapply(mets, function(m) as.character(m$id), character(1))
    lib_name <- vapply(mets, function(m) as.character(m$name %||% m$id), character(1))
    lib_kind <- vapply(mets, function(m) as.character(m$kind %||% NA_character_), character(1))
    hit <- match(matches$met_id, lib_id)
    by_name <- is.na(hit)
    hit[by_name] <- match(matches$met_id[by_name], lib_name)
    by_name2 <- is.na(hit)
    hit[by_name2] <- match(matches$met_name[by_name2], lib_name)
    found <- !is.na(hit)
    rename <- found & (matches$met_name == matches$met_id)
    matches$met_name[rename] <- lib_name[hit[rename]]
    need_kind <- found & is.na(matches$kind)
    matches$kind[need_kind] <- lib_kind[hit[need_kind]]
  }

  # Internal standards: an IS-like kind value ("IS", "ISTD", "internal
  # standard") or, with no kind, an IS token in the id/name (see
  # detect_is_role()). Checked before parent inference.
  k_norm <- gsub("[^a-z]", "", tolower(matches$kind))
  matches$kind[!is.na(k_norm) & k_norm %in% c("is", "istd", "internalstandard", "sil")] <- "internal_standard"
  no_kind <- is.na(matches$kind) | !nzchar(matches$kind)
  is_like <- no_kind & (detect_is_role(matches$met_id) | detect_is_role(matches$met_name))
  matches$kind[is_like] <- "internal_standard"
  still <- is.na(matches$kind) | !nzchar(matches$kind)
  looks_parent <- grepl("^(parent|intact|flp|full[ _-]?length)", matches$met_id, ignore.case = TRUE) |
    grepl("^(parent|intact|flp|full[ _-]?length)", matches$met_name, ignore.case = TRUE)
  matches$kind[still] <- ifelse(looks_parent[still], "parent", "unknown")

  if (!is.null(parent_met_id) && length(parent_met_id) == 1 && nzchar(parent_met_id) &&
      parent_met_id %in% matches$met_id) {
    matches$kind[matches$kind == "parent" & matches$met_id != parent_met_id] <- "unknown"
    matches$kind[matches$met_id == parent_met_id] <- "parent"
  }
  matches
}

#' Guess which metabolite in a data matrix is the intact parent
#'
#' @param matches `read_data_matrix()$matches`, optionally annotated.
#' @return A single `met_id`: the one with `kind == "parent"` if present,
#'   else the non-IS metabolite with the highest total signal.
#' @export
guess_parent_met_id <- function(matches) {
  if (is.null(matches) || nrow(matches) == 0) return(NA_character_)
  if ("kind" %in% names(matches)) {
    p <- unique(matches$met_id[!is.na(matches$kind) & matches$kind == "parent"])
    if (length(p) > 0) return(p[1])
  }
  if ("kind" %in% names(matches)) {
    keep <- is.na(matches$kind) | matches$kind != "internal_standard"
    if (any(keep)) matches <- matches[keep, , drop = FALSE]
  }
  sig <- if ("intensity" %in% names(matches)) "intensity" else "area"
  tot <- tapply(matches[[sig]], matches$met_id, sum, na.rm = TRUE)
  names(tot)[which.max(tot)]
}

## ---- Sample information --------------------------------------------------------

#' Map free-text sample types onto the controlled vocabulary
#'
#' Case-, spacing- and punctuation-insensitive: "Standard", "QC",
#' "Quality Control", "blank", "Matrix Blank" all map. Unmappable values
#' come back as `NA`.
#'
#' @param x Character vector.
#' @return Character vector over `unknown`, `standard`, `quality_control`,
#'   `reagent_blank`, `matrix_blank`, or `NA`.
#' @export
normalize_sample_type <- function(x) {
  norm <- .dm_norm_names(x)
  canon <- c(unknown = "unknown", sample = "unknown", study = "unknown",
             study_sample = "unknown", unk = "unknown",
             standard = "standard", std = "standard", calibrator = "standard",
             calibration = "standard", cal = "standard", calibration_standard = "standard",
             quality_control = "quality_control", qc = "quality_control",
             reagent_blank = "reagent_blank", blank = "reagent_blank", rb = "reagent_blank",
             matrix_blank = "matrix_blank", mb = "matrix_blank")
  unname(canon[norm])
}

#' Read a sample information table
#'
#' Accepts the same file types as [read_table_file()]. Header aliases are
#' mapped onto `sample`, `group`, `timepoint`, `sample_type`, and
#' `concentration` (e.g. "Sample Name", "Time (h)", "Conc", "Type");
#' `sample_type` values are normalized with [normalize_sample_type()]
#' (unmappable values are kept as typed so the caller can flag them; a
#' blank stays blank, which downstream counts as a study sample).
#'
#' @param path File path.
#' @param sheet Worksheet for Excel files.
#' @param name Original file name, when `path` is a temporary upload path.
#' @return A character data.frame with at least `sample`, `group`,
#'   `timepoint`, `sample_type`, `concentration`.
#' @export
read_sample_info <- function(path, sheet = 1, name = NULL) {
  df <- read_table_file(path, sheet = sheet, all_character = TRUE, name = name)
  nn <- .dm_norm_names(names(df))
  i_sample <- .dm_first_match(nn, .DM_SAMPLE_ALIASES)
  if (is.na(i_sample)) stop("sample information needs a 'sample' column")
  out <- data.frame(sample = .dm_clean_sample(df[[i_sample]]), stringsAsFactors = FALSE)
  for (col in names(.DM_META_ALIASES)) {
    aliases <- .DM_META_ALIASES[[col]]
    # "type" is ambiguous in a sample sheet only when a proper sample_type
    # header is absent -- accept it as a fallback.
    if (col == "sample_type") aliases <- c(aliases, "type")
    i <- .dm_first_match(nn, aliases)
    out[[col]] <- if (is.na(i)) "" else trimws(df[[i]])
  }
  mapped <- normalize_sample_type(out$sample_type)
  ok <- !is.na(mapped)
  out$sample_type[ok] <- mapped[ok]
  out <- out[nzchar(out$sample), , drop = FALSE]
  rownames(out) <- NULL
  out
}

## ---- Export ------------------------------------------------------------------

#' Export a data matrix (wide or long) from match results
#'
#' The inverse of [read_data_matrix()]: collapses match results to one
#' value per (metabolite, sample) -- the max across charge states/adducts/
#' oxidation levels, same rule as [build_abundance_matrix()] -- and writes
#' it in a layout [read_data_matrix()] reads straight back.
#'
#' @param matches `ms1_matches` from [annotate_metabolites_batch()] or
#'   `read_data_matrix()$matches`.
#' @param sample_meta Optional sample information, appended to the long
#'   layout.
#' @param format `"wide"` (metabolites x samples) or `"long"`.
#' @param signal_col Signal column; `NULL` auto-detects (area, else
#'   intensity).
#' @return A data.frame.
#' @export
export_data_matrix <- function(matches, sample_meta = NULL, format = c("wide", "long"),
                               signal_col = NULL) {
  format <- match.arg(format)
  if (is.null(matches) || nrow(matches) == 0) return(data.frame())
  if (is.null(signal_col)) signal_col <- .auto_signal_col(matches)
  wide <- build_abundance_matrix(matches, signal_col = signal_col)
  if (format == "wide" || nrow(wide) == 0) return(wide)
  samples <- setdiff(names(wide), c("met_id", "met_name", "kind"))
  meta <- if (!is.null(sample_meta) && nrow(sample_meta) > 0) sample_meta
          else data.frame(sample = samples, stringsAsFactors = FALSE)
  missing <- setdiff(samples, meta$sample)
  if (length(missing) > 0) {
    pad <- meta[rep(1, length(missing)), , drop = FALSE]
    pad[] <- lapply(pad, function(v) if (is.character(v)) "" else NA)
    pad$sample <- missing
    meta <- rbind(meta, pad)
  }
  long <- abundance_long(wide, meta[meta$sample %in% samples, , drop = FALSE])
  long[!is.na(long$intensity), , drop = FALSE]
}
