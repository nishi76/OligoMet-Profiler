# =============================================================================
# internal_standard.R
# One analyte plus one or more internal standards (IS), entered together as
# a multi-record FASTA. The analyte gets the full metabolite library; each IS
# is a single library entry (kind = "internal_standard", no truncations, no
# endonuclease fragments) whose mass comes from its own sequence -- the
# usual analog-IS case (a related but different oligo).
#
# A record is an IS when its header carries role=IS, or, without a role key,
# when its name contains a standalone IS token: IS, SIL, ISTD, IntStd, IS1,
# or the phrase "internal standard". Tokens are matched on word boundaries,
# so Ionis-style "ISIS 420915" stays an analyte.
#
# Quantitation: every metabolite's signal is divided by its assigned IS's
# signal in the same sample (the response ratio). One global IS serves every
# metabolite unless a per-metabolite assignment overrides it. The IS itself
# is never blank-corrected (a zero blank contains IS by design), so a
# blank-corrected ratio is corrected analyte / raw IS.
# =============================================================================

## ---- Role detection -------------------------------------------------------------

#' Does a sequence name mark an internal standard?
#'
#' @param name Character vector of record names/headers.
#' @return Logical vector: `TRUE` when the name contains a standalone `IS`,
#'   `SIL`, `ISTD`, `IntStd`, `IS<n>` token or the phrase "internal
#'   standard" (case-insensitive). "ISIS 420915" or "analysis" do not match.
#' @export
detect_is_role <- function(name) {
  vapply(tolower(as.character(name)), function(n) {
    if (is.na(n)) return(FALSE)
    if (grepl("internal[ _-]*standard", n)) return(TRUE)
    tok <- strsplit(n, "[^a-z0-9]+")[[1]]
    any(tok %in% c("is", "sil", "istd", "intstd") | grepl("^(is|istd)[0-9]+$", tok))
  }, logical(1), USE.NAMES = FALSE)
}

## ---- Multi-record sequence input -----------------------------------------------

#' Parse one analyte plus internal standards from FASTA-style text
#'
#' Each `>` header starts a record; the lines after it are that record's
#' sequence in any notation [parse_input()] accepts (triplet or
#' OligoDistiller). Text with no `>` line is a single analyte, as before.
#' Header keys: `role=IS` / `role=analyte` (overrides name detection),
#' `5'-<conj>` / `3'-<conj>` (terminal conjugates, same as [parse_fasta()]).
#'
#' @param text The sequence text (pasted or read from a `.fasta`/`.txt`).
#' @param dict Chemistry dictionary.
#' @param default_name Analyte name when the text has no header.
#' @return A list: `analyte` (list with `name`, `spec`), `internal_standards`
#'   (list of the same), `table` (one row per record: `record`, `name`,
#'   `role`, `role_source`, `length`, `mono_mass`, `problem`), and
#'   `problems` (character; empty when the set is usable: exactly one
#'   analyte, every record parsed).
#' @export
parse_sequence_set <- function(text, dict = STANDARD_DICT, default_name = "analyte") {
  lines <- trimws(strsplit(text %||% "", "\r?\n")[[1]])
  lines <- lines[nzchar(lines)]
  if (!any(startsWith(lines, ">"))) {
    recs <- list(list(header = default_name, seq = paste(lines, collapse = "")))
  } else {
    start <- which(startsWith(lines, ">"))
    if (start[1] > 1) start <- c(1L, start)  # sequence before the first header
    ends <- c(start[-1] - 1L, length(lines))
    recs <- lapply(seq_along(start), function(i) {
      blk <- lines[start[i]:ends[i]]
      hdr <- if (startsWith(blk[1], ">")) sub("^>\\s*", "", blk[1]) else default_name
      body <- blk[!startsWith(blk, ">")]
      list(header = hdr, seq = paste(body, collapse = ""))
    })
  }

  rows <- list(); parsed <- list()
  for (i in seq_along(recs)) {
    h <- recs[[i]]$header
    role_key <- regmatches(h, regexpr("(?i)\\brole\\s*=\\s*[A-Za-z_ ]+?(?=\\s|$)", h, perl = TRUE))
    m5 <- regmatches(h, regexpr("5'-[A-Za-z0-9_]+", h))
    m3 <- regmatches(h, regexpr("3'-[A-Za-z0-9_]+", h))
    name <- trimws(gsub("\\s+", " ", gsub("(?i)\\brole\\s*=\\s*\\S+|5'-[A-Za-z0-9_]+|3'-[A-Za-z0-9_]+",
                                          "", h, perl = TRUE)))
    if (!nzchar(name)) name <- sprintf("record_%d", i)
    if (length(role_key) > 0) {
      v <- tolower(trimws(sub("(?i)^role\\s*=\\s*", "", role_key, perl = TRUE)))
      role <- if (v %in% c("is", "istd", "internal_standard", "internal standard", "sil")) "internal_standard" else "analyte"
      src <- "header role= key"
    } else if (detect_is_role(name)) {
      role <- "internal_standard"; src <- "IS token in name"
    } else {
      role <- "analyte"; src <- "default"
    }
    spec <- NULL; problem <- ""
    if (!nzchar(recs[[i]]$seq)) {
      problem <- "no sequence lines"
    } else {
      spec <- tryCatch(parse_input(recs[[i]]$seq, dict = dict), error = function(e) {
        problem <<- conditionMessage(e); NULL
      })
    }
    mono <- NA_real_
    if (!is.null(spec)) {
      if (length(m5) > 0) spec$conj5 <- sub("^5'-", "", m5)
      if (length(m3) > 0) spec$conj3 <- sub("^3'-", "", m3)
      mono <- tryCatch(metabolite_mass_info(.make_met("tmp", name, role, role, NA, spec, dict = dict),
                                            dict)$mono_mass, error = function(e) NA_real_)
    }
    parsed[[i]] <- list(name = name, role = role, spec = spec)
    rows[[i]] <- data.frame(record = i, name = name, role = role, role_source = src,
                            length = if (is.null(spec)) NA_integer_ else spec$n,
                            mono_mass = mono, problem = problem, stringsAsFactors = FALSE)
  }
  tbl <- do.call(rbind, rows)
  problems <- sprintf("record %d (%s): %s", tbl$record[nzchar(tbl$problem)],
                      tbl$name[nzchar(tbl$problem)], tbl$problem[nzchar(tbl$problem)])
  n_an <- sum(tbl$role == "analyte")
  if (n_an == 0) problems <- c(problems, "no analyte record -- every record was read as an internal standard")
  if (n_an > 1) problems <- c(problems, sprintf(paste0(
    "%d analyte records (%s); one analyte per run is supported -- mark the others ",
    "with role=IS in the header, or remove them"), n_an, paste(tbl$name[tbl$role == "analyte"], collapse = ", ")))
  ok <- vapply(parsed, function(p) !is.null(p$spec), logical(1))
  list(analyte = if (n_an >= 1) parsed[[which(tbl$role == "analyte")[1]]] else NULL,
       internal_standards = parsed[tbl$role == "internal_standard" & ok],
       table = tbl, problems = problems)
}

#' Library entries for internal standards
#'
#' @param is_records `parse_sequence_set()$internal_standards`.
#' @param dict Chemistry dictionary.
#' @return A list of metabolite objects with ids `IS01`, `IS02`, ...,
#'   `kind = "internal_standard"`, no truncations or fragments.
#' @export
build_is_metabolites <- function(is_records, dict = STANDARD_DICT) {
  lapply(seq_along(is_records), function(i) {
    r <- is_records[[i]]
    .make_met(sprintf("IS%02d", i), r$name, "internal_standard", "internal standard",
              NA, r$spec, dict = dict)
  })
}

#' Check internal standards against the analyte's metabolites
#'
#' An IS whose mass equals a metabolite's cannot be told apart from it by
#' MS1; one within an isotope-envelope width shares ions with it in the
#' summed XIC unless the two separate chromatographically.
#'
#' @param mets Analyte metabolite library.
#' @param is_mets [build_is_metabolites()] output.
#' @param dict Chemistry dictionary.
#' @param ppm Same-mass tolerance.
#' @param envelope_da Mass difference below which isotope envelopes overlap.
#' @return data.frame of conflicts (`is_id`, `is_name`, `met_id`,
#'   `met_name`, `delta_da`, `severity`, `message`); empty when clean.
#' @export
check_is_interference <- function(mets, is_mets, dict = STANDARD_DICT, ppm = 10,
                                  envelope_da = 12) {
  if (length(is_mets) == 0) return(data.frame())
  mm <- vapply(mets, function(m) metabolite_mass_info(m, dict)$mono_mass, numeric(1))
  out <- list()
  for (ism in is_mets) {
    im <- metabolite_mass_info(ism, dict)$mono_mass
    d <- mm - im
    for (j in which(abs(d) < envelope_da)) {
      same <- abs(d[j]) / im * 1e6 <= ppm
      out[[length(out) + 1]] <- data.frame(
        is_id = ism$id, is_name = ism$name, met_id = mets[[j]]$id, met_name = mets[[j]]$name,
        delta_da = round(d[j], 4), severity = if (same) "error" else "warning",
        message = if (same) "same monoisotopic mass: indistinguishable by MS1"
                  else "isotope envelopes overlap: needs chromatographic separation",
        stringsAsFactors = FALSE)
    }
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

## ---- IS normalization ---------------------------------------------------------------

#' Ids of internal-standard rows in a match table
#' @param matches `ms1_matches`-shaped table with a `kind` column.
#' @return Character vector of met_ids with `kind == "internal_standard"`.
#' @export
is_met_ids <- function(matches) {
  if (is.null(matches) || nrow(matches) == 0 || !"kind" %in% names(matches)) return(character(0))
  unique(matches$met_id[!is.na(matches$kind) & matches$kind == "internal_standard"])
}

#' Which IS each metabolite is normalized to
#'
#' @param met_ids Metabolites to assign.
#' @param global_is The IS used for every metabolite without an override.
#' @param overrides Optional data.frame(`met_id`, `is_id`).
#' @return data.frame(`met_id`, `is_id`, `assignment`).
#' @export
is_assignment <- function(met_ids, global_is, overrides = NULL) {
  out <- data.frame(met_id = met_ids, is_id = global_is %||% NA_character_,
                    assignment = "global", stringsAsFactors = FALSE)
  if (!is.null(overrides) && nrow(overrides) > 0) {
    hit <- match(out$met_id, overrides$met_id)
    ov <- !is.na(hit) & nzchar(overrides$is_id[hit])
    out$is_id[ov] <- overrides$is_id[hit[ov]]
    out$assignment[ov] <- "per-metabolite"
  }
  out
}

#' Divide every metabolite's signal by its internal standard's
#'
#' Adds `<col>_isr` (response ratio) next to each signal column present:
#' `intensity`, `area`, and their `_bcorr` blank-corrected versions. The
#' denominator is always the RAW IS signal (max per sample across the IS's
#' match rows), never a blank-corrected one. A sample with no IS signal
#' gets `NA`.
#'
#' @param matches `ms1_matches`-shaped table.
#' @param global_is IS met_id used for every metabolite without an override.
#' @param overrides Optional data.frame(`met_id`, `is_id`).
#' @return `matches` with `_isr` columns added; IS rows get `NA`. The
#'   assignment used is attached as `attr(, "is_assignment")`, and samples
#'   lacking IS signal as `attr(, "is_missing_samples")`.
#' @export
apply_is_normalization <- function(matches, global_is, overrides = NULL) {
  if (is.null(matches) || nrow(matches) == 0 || is.null(global_is) || !nzchar(global_is)) return(matches)
  is_ids <- unique(c(global_is, if (!is.null(overrides)) overrides$is_id, is_met_ids(matches)))
  is_ids <- is_ids[!is.na(is_ids) & nzchar(is_ids)]
  analytes <- setdiff(unique(matches$met_id), is_ids)
  asg <- is_assignment(analytes, global_is, overrides)
  cols <- intersect(c("intensity", "area", "intensity_bcorr", "area_bcorr"), names(matches))
  cols <- cols[vapply(cols, function(cc) any(!is.na(matches[[cc]])), logical(1))]
  row_is <- asg$is_id[match(matches$met_id, asg$met_id)]
  missing <- character(0)
  for (cc in cols) {
    raw <- sub("_bcorr$", "", cc)
    isd <- matches[matches$met_id %in% is_ids & !is.na(matches[[raw]]), c("sample", "met_id", raw)]
    denom <- if (nrow(isd) > 0) {
      agg <- stats::aggregate(isd[[raw]], list(sample = isd$sample, met_id = isd$met_id), max)
      agg$x[match(paste(matches$sample, row_is), paste(agg$sample, agg$met_id))]
    } else rep(NA_real_, nrow(matches))
    denom[!is.na(denom) & denom <= 0] <- NA_real_
    v <- matches[[cc]] / denom
    v[matches$met_id %in% is_ids] <- NA_real_
    matches[[paste0(cc, "_isr")]] <- v
    missing <- union(missing, unique(matches$sample[!(matches$met_id %in% is_ids) & is.na(denom)]))
  }
  attr(matches, "is_assignment") <- asg
  attr(matches, "is_missing_samples") <- missing
  matches
}

## ---- IS response monitoring ---------------------------------------------------------

#' Internal-standard response per sample, against an acceptance window
#'
#' The reference is the mean IS signal of the calibration standards and
#' QCs (every non-blank sample when there are none). Each sample's IS
#' signal is expressed as % of that mean and flagged outside `window`.
#' Reagent/matrix blanks are listed but not evaluated.
#'
#' @param matches `ms1_matches`-shaped table.
#' @param sample_meta Sample table (`sample`, `sample_type`); its row order
#'   is taken as injection order.
#' @param is_id The IS met_id.
#' @param window Acceptance window, % of the reference mean.
#' @param signal_col Signal column (`NULL`: area if present, else intensity).
#' @return data.frame: `run_order`, `sample`, `sample_type`, `is_signal`,
#'   `pct_of_reference`, `status`; the reference mean and window as
#'   attributes.
#' @export
is_response_summary <- function(matches, sample_meta, is_id, window = c(50, 150),
                                signal_col = NULL) {
  if (is.null(matches) || nrow(matches) == 0 || is.null(is_id) || !nzchar(is_id)) return(data.frame())
  if (is.null(signal_col)) signal_col <- .auto_signal_col(matches)
  d <- matches[matches$met_id == is_id & !is.na(matches[[signal_col]]), , drop = FALSE]
  samples <- if (!is.null(sample_meta) && nrow(sample_meta) > 0) sample_meta$sample else unique(matches$sample)
  st <- if (!is.null(sample_meta) && "sample_type" %in% names(sample_meta))
    sample_meta$sample_type[match(samples, sample_meta$sample)] else rep("unknown", length(samples))
  st[is.na(st) | !nzchar(st)] <- "unknown"
  sig <- if (nrow(d) > 0) tapply(d[[signal_col]], d$sample, max)[samples] else rep(NA_real_, length(samples))
  sig <- as.numeric(sig)
  is_blank <- st %in% c("reagent_blank", "matrix_blank")
  ref_pool <- st %in% c("standard", "quality_control") & !is.na(sig)
  if (!any(ref_pool)) ref_pool <- !is_blank & !is.na(sig)
  ref <- if (any(ref_pool)) mean(sig[ref_pool]) else NA_real_
  pct <- 100 * sig / ref
  status <- ifelse(is_blank, "blank (not evaluated)",
            ifelse(is.na(sig), "IS not detected",
            ifelse(pct < window[1], "low", ifelse(pct > window[2], "high", "ok"))))
  out <- data.frame(run_order = seq_along(samples), sample = samples, sample_type = st,
                    is_signal = sig, pct_of_reference = round(pct, 1), status = status,
                    stringsAsFactors = FALSE)
  attr(out, "reference_mean") <- ref
  attr(out, "window") <- window
  out
}

#' Plot IS response across the run
#'
#' @param summary [is_response_summary()] output.
#' @return A ggplot: % of reference IS response by injection order, coloured
#'   by sample type, with the acceptance window shaded.
#' @export
plot_is_response <- function(summary) {
  if (is.null(summary) || nrow(summary) == 0 || all(is.na(summary$pct_of_reference))) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "No internal-standard signal to plot"))
  }
  w <- attr(summary, "window") %||% c(50, 150)
  d <- summary[!is.na(summary$pct_of_reference), , drop = FALSE]
  d$flagged <- d$status %in% c("low", "high")
  types <- sort(unique(d$sample_type))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$run_order, y = .data$pct_of_reference)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = w[1], ymax = w[2],
                      fill = "#75A025", alpha = 0.08) +
    ggplot2::geom_hline(yintercept = c(w[1], 100, w[2]), linetype = c("dotted", "solid", "dotted"),
                        color = "#B7B2A7") +
    ggplot2::geom_point(ggplot2::aes(color = .data$sample_type, shape = .data$flagged), size = 2.8) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
    ggplot2::scale_color_manual(values = stats::setNames(grDevices::colorRampPalette(
      c("#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C"))(length(types)), types), name = NULL) +
    ggplot2::labs(x = "Injection order (sample table order)", y = "IS response (% of reference mean)",
                  title = "Internal standard response",
                  subtitle = sprintf("Reference: mean of standards/QCs. Window %g-%g%%; x = outside window.",
                                     w[1], w[2])) +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
}
