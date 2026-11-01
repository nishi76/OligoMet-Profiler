# =============================================================================
# regulatory_report.R
# Regulatory-format metabolite-identification study report, populated from a
# live run. The scientific sections (test article, study design, methods,
# results, kinetics, structure maps, figures and tables) are computed from the
# metabolite library, the batch match table, the sample information and the
# degradation summary; the administrative and GLP fields (sponsor, CRO, study
# number, dates, signatures) come from a study-details list and fall back to
# bracketed placeholders when not supplied.
#
# One assembled model (see .reg_model()) holds the cover-page fields and an
# ordered list of content blocks (headings, paragraphs, tables, figures,
# structure maps). Two renderers walk the same blocks:
#   build_regulatory_report_html()  -- self-contained HTML (inline SVG + PNG)
#   build_regulatory_report_docx()  -- native Word via officer (cover page,
#                                      approval page, table of contents,
#                                      running header/footer with page X of Y,
#                                      numbered headings, flextable tables)
#
# All quantities are computed on study samples only: calibration standards,
# QCs and blanks (see .is_study_sample()) are listed in the sample inventory
# but never enter the "% of parent", degradation or kinetics figures.
# =============================================================================

## ---- Study-details fields and placeholder helper ---------------------------
.REG_DETAIL_FIELDS <- c(
  "sponsor", "sponsor_address", "cro", "cro_address", "study_title",
  "study_number", "project_code", "protocol", "test_article",
  "test_article_lot", "is_name", "study_director", "species", "route",
  "matrix", "report_status", "report_date")

.reg_get <- function(details, key, default = NULL) {
  v <- details[[key]]
  if (is.null(v) || !nzchar(trimws(as.character(v)[1]))) {
    if (!is.null(default)) return(default)
    return(sprintf("[%s]", gsub("_", " ", key)))
  }
  trimws(as.character(v)[1])
}

# Text width of the Word body (A4 portrait, 1 inch margins), in inches.
.REG_BODY_W <- 6.27

# Non-ASCII characters built from code points, so they stay UTF-8 even when
# R runs in a C/POSIX locale (officer otherwise writes their escaped bytes).
.REG_BULLET <- intToUtf8(0x2022)
.REG_SUP2 <- intToUtf8(0x00b2)

.REG_KIND_LABELS <- c(
  parent = "Parent (full-length)",
  exo_3p = "3'-exonuclease product",
  exo_5p = "5'-exonuclease product",
  endo_5frag = "Endonuclease fragment (5')",
  endo_3frag = "Endonuclease fragment (3')",
  endo_5frag_p = "Endonuclease fragment (5', phosphorylated)",
  endo_3frag_p = "Endonuclease fragment (3', phosphorylated)",
  conj_catab = "Conjugate catabolite",
  internal_standard = "Internal standard",
  unknown = "Unassigned")

.reg_kind_label <- function(k) {
  k <- as.character(k)
  out <- unname(.REG_KIND_LABELS[k])
  out[is.na(out)] <- k[is.na(out)]
  out
}

.reg_dict_name <- function(dict, code) {
  nm <- dict[[code]]$name
  if (is.null(nm)) return(code)
  sub("\\s+--.*$", "", nm)
}

# Significant-figure formatting that never drops into scientific notation.
.reg_sig <- function(x, digits = 3) {
  out <- formatC(signif(x, digits), digits = digits, format = "fg", flag = "#")
  out <- sub("\\.$", "", trimws(out))
  out[is.na(x)] <- "-"
  out
}
.reg_fix <- function(x, digits = 1) {
  out <- formatC(x, digits = digits, format = "f", big.mark = ",")
  out[is.na(x)] <- "-"
  out
}

.reg_version <- function() {
  v <- tryCatch(as.character(utils::packageVersion("OligoMetProfiler")), error = function(e) NULL)
  if (!is.null(v)) return(v)
  for (p in c("DESCRIPTION", file.path("..", "DESCRIPTION"), file.path("..", "..", "DESCRIPTION"))) {
    d <- tryCatch(read.dcf(p, fields = c("Package", "Version")), error = function(e) NULL)
    if (!is.null(d) && identical(unname(d[1, "Package"]), "OligoMetProfiler")) return(unname(d[1, "Version"]))
  }
  "unknown"
}

## ---- Figure / image helpers ------------------------------------------------
.reg_fig_png <- function(plot, width = 6.2, height = 3.4, dpi = 200) {
  if (is.null(plot)) return(NULL)
  f <- tempfile(fileext = ".png")
  ok <- tryCatch({
    ggplot2::ggsave(f, plot = plot, width = width, height = height, dpi = dpi,
                    units = "in", bg = "white"); TRUE
  }, error = function(e) FALSE)
  if (ok && file.exists(f)) f else NULL
}

.reg_svg_png <- function(svg, width_px = 1400) {
  if (is.null(svg) || !nzchar(svg)) return(NULL)
  f <- tempfile(fileext = ".png")
  ok <- tryCatch({
    rsvg::rsvg_png(charToRaw(svg), file = f, width = width_px); TRUE
  }, error = function(e) FALSE)
  if (ok && file.exists(f)) f else NULL
}

# Width and height in pixels, read from the PNG IHDR chunk.
.reg_png_dims <- function(png_path) {
  b <- tryCatch(readBin(png_path, "raw", 24), error = function(e) raw(0))
  if (length(b) < 24) return(NULL)
  be <- function(x) sum(as.numeric(x) * 256^(3:0))
  c(be(b[17:20]), be(b[21:24]))
}

.reg_b64img <- function(png_path) {
  if (is.null(png_path) || !file.exists(png_path)) return(NULL)
  paste0("data:image/png;base64,", base64enc::base64encode(png_path))
}

.reg_xml_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

## ---- Test-article chemistry ------------------------------------------------
# Display letter per nucleobase code: canonical bases as-is, 5-methyl
# pyrimidines as the conventional "mC"/"mU", anything else bracketed.
.reg_base_letter <- function(code) {
  map <- c(A = "A", G = "G", C = "C", T = "T", U = "U", S = "mC", m5C = "mC",
           hm5C = "hmC", m5U = "mU", psi = intToUtf8(0x03a8), I = "I", D = "D")
  out <- unname(map[code])
  out[is.na(out)] <- paste0("[", code[is.na(out)], "]")
  out
}

.reg_runs <- function(x) {
  r <- rle(as.character(x))
  ends <- cumsum(r$lengths)
  data.frame(value = r$values, start = ends - r$lengths + 1L, end = ends,
             len = r$lengths, stringsAsFactors = FALSE)
}

.reg_chemistry <- function(spec, dict) {
  if (is.null(spec) || is.null(spec$n)) return(NULL)
  n <- spec$n
  bases <- .reg_base_letter(spec$bases)
  seq_txt <- paste(bases, collapse = "")
  sr <- .reg_runs(spec$sugars)
  sugar_txt <- paste(sprintf("%s %s",
    ifelse(sr$start == sr$end, sprintf("position %d", sr$start),
           sprintf("positions %d-%d", sr$start, sr$end)),
    vapply(sr$value, function(s) .reg_dict_name(dict, s), "")), collapse = "; ")
  design <- if (nrow(sr) == 3 && sr$value[2] == "d" && sr$value[1] != "d" && sr$value[3] != "d") {
    sprintf("%d-%d-%d gapmer (%s wings, DNA gap)", sr$len[1], sr$len[2], sr$len[3],
            .reg_dict_name(dict, sr$value[1]))
  } else if (nrow(sr) == 1) {
    sprintf("uniformly modified (%s)", .reg_dict_name(dict, sr$value[1]))
  } else "mixed-chemistry"
  lk <- spec$linkages[seq_len(max(0, n - 1))]
  lk <- lk[!is.na(lk)]
  lt <- table(lk)
  backbone <- if (length(lt) == 0) "-" else paste(sprintf("%d %s", as.integer(lt),
    vapply(names(lt), function(l) .reg_dict_name(dict, l), "")), collapse = ", ")
  mod_b <- spec$bases[!spec$bases %in% c("A", "G", "C", "T", "U")]
  mod_txt <- if (length(mod_b) == 0) "none" else {
    tb <- table(mod_b)
    paste(sprintf("%d x %s", as.integer(tb),
                  vapply(names(tb), function(b) .reg_dict_name(dict, b), "")), collapse = ", ")
  }
  c5 <- spec$conj5 %||% "none"; c3 <- spec$conj3 %||% "none"
  conj_txt <- paste(c(if (!identical(c5, "none")) paste0("5': ", .reg_dict_name(dict, c5)),
                      if (!identical(c3, "none")) paste0("3': ", .reg_dict_name(dict, c3))),
                    collapse = "; ")
  if (!nzchar(conj_txt)) conj_txt <- "none"
  per_pos <- data.frame(
    Position = seq_len(n),
    Nucleobase = vapply(spec$bases, function(b) .reg_dict_name(dict, b), ""),
    Sugar = vapply(spec$sugars, function(s) .reg_dict_name(dict, s), ""),
    `3' linkage` = c(vapply(lk, function(l) .reg_dict_name(dict, l), "")[seq_len(n - 1)], "3'-terminus"),
    check.names = FALSE, stringsAsFactors = FALSE)
  list(n = n, seq_txt = seq_txt, sugar_txt = sugar_txt, design = design,
       backbone = backbone, mod_txt = mod_txt, conj_txt = conj_txt, per_pos = per_pos)
}

## ---- Sample bookkeeping ----------------------------------------------------
.reg_signal_col <- function(matches) {
  if (!is.null(matches) && "area" %in% names(matches) && any(!is.na(matches$area))) "area" else "intensity"
}

.reg_sample_frame <- function(samples, meta) {
  df <- data.frame(sample = samples, stringsAsFactors = FALSE)
  for (col in c("group", "timepoint", "matrix", "sample_type")) {
    v <- if (!is.null(meta) && col %in% names(meta)) as.character(meta[[col]][match(samples, meta$sample)]) else rep(NA_character_, length(samples))
    v[is.na(v)] <- ""
    df[[col]] <- trimws(v)
  }
  df$study <- .is_study_sample(df$sample_type)
  df$condition <- vapply(seq_len(nrow(df)), function(i) {
    parts <- c(df$matrix[i], df$group[i]); parts <- parts[nzchar(parts)]
    if (length(parts) == 0) "all samples" else gsub("_", " ", paste(parts, collapse = " / "), fixed = TRUE)
  }, "")
  df$time <- suppressWarnings(as.numeric(df$timepoint))
  df
}

.reg_type_label <- function(t) {
  map <- c(standard = "Calibration standard", quality_control = "Quality control",
           reagent_blank = "Reagent blank", matrix_blank = "Matrix blank")
  out <- unname(map[t]); out[is.na(out)] <- "Study sample"
  out
}

# Per-sample metabolite signal (wide), restricted to study samples, plus the
# per-metabolite parent-relative and total-drug-related metrics.
.reg_abundance <- function(matches, sframe, threshold) {
  if (is.null(matches) || nrow(matches) == 0) return(NULL)
  sig <- .reg_signal_col(matches)
  mat <- build_abundance_matrix(matches, signal_col = sig)
  if (nrow(mat) == 0) return(NULL)
  study <- intersect(sframe$sample[sframe$study], names(mat))
  if (length(study) == 0) return(NULL)
  p_row <- which(mat$kind == "parent")[1]
  if (is.na(p_row)) return(NULL)
  X <- as.matrix(mat[, study, drop = FALSE]); storage.mode(X) <- "double"
  X[is.na(X)] <- 0
  parent <- X[p_row, ]
  analyte <- mat$kind != "internal_standard"
  total <- colSums(X[analyte, , drop = FALSE])
  pct_parent <- sweep(X, 2, ifelse(parent > 0, parent, NA), "/") * 100
  pct_total <- sweep(X, 2, ifelse(total > 0, total, NA), "/") * 100
  cond <- sframe$condition[match(study, sframe$sample)]
  tm <- sframe$timepoint[match(study, sframe$sample)]
  where <- ifelse(nzchar(tm), paste0(cond, ", t = ", tm), cond)
  rows <- lapply(seq_len(nrow(mat)), function(i) {
    pp <- pct_parent[i, ]; pt <- pct_total[i, ]
    k <- if (all(is.na(pp))) NA_integer_ else which.max(pp)
    data.frame(met_id = mat$met_id[i], met_name = mat$met_name[i], kind = mat$kind[i],
               n_detected = sum(X[i, ] > 0), n_study = length(study),
               max_pct = if (is.na(k)) NA_real_ else pp[k],
               max_pct_total = if (all(is.na(pt))) NA_real_ else max(pt, na.rm = TRUE),
               at = if (is.na(k)) "-" else where[k], stringsAsFactors = FALSE)
  })
  per_met <- do.call(rbind, rows)
  per_met$reported <- per_met$kind == "parent" |
    (per_met$kind != "internal_standard" & !is.na(per_met$max_pct) & per_met$max_pct >= threshold)
  list(signal = sig, mat = mat, study = study, pct_parent = pct_parent,
       pct_total = pct_total, per_met = per_met, parent_row = p_row)
}

# First-order loss of the parent per condition: ln(% remaining) vs time.
.reg_kinetics <- function(per_cond_samples) {
  out <- lapply(split(per_cond_samples, per_cond_samples$condition), function(d) {
    d <- d[!is.na(d$time) & !is.na(d$remaining) & d$remaining > 0, , drop = FALSE]
    res <- data.frame(Condition = d$condition[1] %||% NA, k = NA_real_, t_half = NA_real_,
                      r2 = NA_real_, n = nrow(d), t_max = if (nrow(d)) max(d$time) else NA_real_,
                      stringsAsFactors = FALSE)
    if (nrow(d) < 3 || length(unique(d$time)) < 2) return(res)
    fit <- tryCatch(stats::lm(log(remaining) ~ time, data = d), error = function(e) NULL)
    if (is.null(fit)) return(res)
    res$k <- -unname(stats::coef(fit)[2])
    res$r2 <- summary(fit)$r.squared
    res$t_half <- if (res$k > 0) log(2) / res$k else Inf
    res
  })
  do.call(rbind, out)
}

## ---- Assemble the report model ---------------------------------------------
# ctx fields (all optional; the report degrades gracefully):
#   details   named list of administrative/GLP fields (see .REG_DETAIL_FIELDS)
#   spec      parent oligo_spec (bases/sugars/linkages/conj5/conj3/n)
#   mets      metabolite library (list of metabolite objects)
#   dict      chemistry dictionary
#   matches   batch ms1_matches (met_id, met_name, kind, sample, signal, ...)
#   ms2_confirmations  MS2 confirmation table (met_id, ...), or NULL
#   abundance parent_relative_abundance() table (unused; recomputed on study samples)
#   degradation  degradation_summary() result
#   sample_meta  sample information (sample, group, timepoint, matrix, sample_type)
#   params    list(ppm_tol, adducts, mods, z_min, z_max, max_oxid, pct_min, ...)
#' Assemble the regulatory-report model from a run context
#'
#' Computes the structured content shared by the HTML and Word renderers:
#' cover-page fields, the approval table, and an ordered list of content
#' blocks (headings, paragraphs, tables, figures and structure maps) covering
#' the summary, test article, study design, methods, results, kinetics,
#' discussion, conclusions and appendices. Quantities are computed on study
#' samples only; standards, QCs and blanks are inventoried but excluded.
#'
#' @param ctx A context list; see the file header for the recognized fields.
#'   All are optional and the report fills gaps with bracketed placeholders.
#' @return A list consumed by [build_regulatory_report_html()] and
#'   [build_regulatory_report_docx()].
#' @seealso [build_regulatory_report_html()], [build_regulatory_report_docx()]
#' @export
.reg_model <- function(ctx) {
  details <- ctx$details %||% list()
  params <- ctx$params %||% list()
  dict <- ctx$dict %||% STANDARD_DICT
  threshold <- params$pct_min %||% 0.1
  thr <- format(threshold, drop0trailing = TRUE, scientific = FALSE)
  ta <- .reg_get(details, "test_article")
  title <- .reg_get(details, "study_title")
  report_date <- .reg_get(details, "report_date", format(Sys.Date(), "%d %B %Y"))
  status <- .reg_get(details, "report_status", "Draft")

  # ---- block builders (numbered tables and figures) ----
  blocks <- list(); n_tab <- 0L; n_fig <- 0L
  add <- function(b) blocks[[length(blocks) + 1L]] <<- b
  H1 <- function(t) add(list(type = "h1", text = t))
  H2 <- function(t) add(list(type = "h2", text = t))
  H3 <- function(t) add(list(type = "h3", text = t))
  P <- function(...) add(list(type = "p", text = paste0(...)))
  NOTE <- function(...) add(list(type = "note", text = paste0(...)))
  BUL <- function(items) add(list(type = "bullets", items = items))
  TAB <- function(df, caption, widths = NULL, header = TRUE) {
    if (is.null(df) || nrow(df) == 0) return(invisible(NA_integer_))
    n_tab <<- n_tab + 1L
    add(list(type = "table", df = df, caption = sprintf("Table %d. %s", n_tab, caption),
             widths = widths, header = header))
    invisible(n_tab)
  }
  FIG <- function(plot, caption, w = 6.2, h = 3.4) {
    if (is.null(plot)) return(invisible(NA_integer_))
    n_fig <<- n_fig + 1L
    add(list(type = "figure", plot = plot, caption = sprintf("Figure %d. %s", n_fig, caption), w = w, h = h))
    invisible(n_fig)
  }
  SVG <- function(svg, caption) {
    n_fig <<- n_fig + 1L
    add(list(type = "svg", svg = svg, caption = sprintf("Figure %d. %s", n_fig, caption)))
    invisible(n_fig)
  }
  BREAK <- function() add(list(type = "pagebreak"))

  # ---- parent identity ----
  parent_met <- if (!is.null(ctx$mets)) Find(function(m) identical(m$kind, "parent"), ctx$mets) else NULL
  chem <- .reg_chemistry(ctx$spec, dict)
  mass_info <- NULL
  if (!is.null(parent_met)) {
    mass_info <- tryCatch(metabolite_mass_info(parent_met, dict), error = function(e) NULL)
  } else if (!is.null(ctx$spec)) {
    f <- tryCatch(assemble_oligo_formula(ctx$spec$bases, ctx$spec$sugars, ctx$spec$linkages,
                    ctx$spec$conj5, ctx$spec$conj3, dict), error = function(e) NULL)
    if (!is.null(f)) mass_info <- list(formula_str = format_formula(f), mono_mass = formula_mass(f, TRUE),
                                       avg_mass = tryCatch(formula_mass(f, FALSE), error = function(e) NA_real_))
  }
  lib_prefix <- if (!is.null(parent_met)) parent_met$name %||% "" else ""

  # ---- samples, abundance, degradation ----
  matches <- ctx$matches
  meta <- ctx$sample_meta
  samples <- if (!is.null(matches) && nrow(matches) > 0) unique(matches$sample) else character(0)
  sframe <- .reg_sample_frame(samples, meta)
  ab <- .reg_abundance(matches, sframe, threshold)
  ms2c <- ctx$ms2_confirmations
  ms2_ran <- !is.null(ms2c)
  ms2_ids <- if (ms2_ran && nrow(ms2c) > 0) unique(ms2c$met_id) else character(0)

  st <- sframe[sframe$study, , drop = FALSE]
  has_time <- nrow(st) > 0 && any(!is.na(st$time))
  conds <- unique(st$condition)
  times <- sort(unique(st$time[!is.na(st$time)]))
  tunit <- params$time_unit %||% "h"

  per_sample <- NULL
  if (!is.null(ab)) {
    sig <- ab$mat[ab$parent_row, ab$study]
    per_sample <- data.frame(sample = ab$study, parent = as.numeric(sig),
      total = colSums(as.matrix(ab$mat[ab$mat$kind != "internal_standard", ab$study, drop = FALSE]), na.rm = TRUE),
      stringsAsFactors = FALSE)
    per_sample <- merge(per_sample, sframe[, c("sample", "condition", "group", "matrix", "timepoint", "time")], by = "sample")
    per_sample$pct_deg <- ifelse(per_sample$total > 0, 100 * (1 - per_sample$parent / per_sample$total), NA)
    base <- vapply(per_sample$condition, function(cn) {
      d <- per_sample[per_sample$condition == cn & !is.na(per_sample$time), , drop = FALSE]
      if (nrow(d) == 0) return(NA_real_)
      mean(d$parent[d$time == min(d$time)])
    }, numeric(1))
    per_sample$remaining <- ifelse(!is.na(base) & base > 0, 100 * per_sample$parent / base, NA)
    per_sample <- per_sample[order(per_sample$condition, per_sample$time, per_sample$sample), , drop = FALSE]
  }

  cond_tab <- NULL; kin <- NULL
  if (!is.null(per_sample) && nrow(per_sample) > 0) {
    key <- if (has_time) paste(per_sample$condition, per_sample$time, sep = "\r") else per_sample$condition
    cond_tab <- do.call(rbind, lapply(split(per_sample, factor(key, levels = unique(key))), function(d) {
      data.frame(condition = d$condition[1], time = d$time[1], n = nrow(d),
                 parent = mean(d$parent), rem = mean(d$remaining), rem_sd = stats::sd(d$remaining),
                 deg = mean(d$pct_deg), deg_sd = stats::sd(d$pct_deg), stringsAsFactors = FALSE)
    }))
    if (has_time) kin <- .reg_kinetics(per_sample)
  }

  pm <- if (!is.null(ab)) ab$per_met else NULL
  rep_mets <- if (!is.null(pm)) pm[pm$reported, , drop = FALSE] else NULL
  if (!is.null(rep_mets)) rep_mets <- rep_mets[order(rep_mets$kind != "parent", -rep_mets$max_pct), , drop = FALSE]
  degr <- if (!is.null(rep_mets)) rep_mets[rep_mets$kind != "parent", , drop = FALSE] else NULL
  is_rows <- if (!is.null(pm)) pm[pm$kind == "internal_standard", , drop = FALSE] else NULL

  # library lookup for a reported metabolite: by name suffix after the
  # oligo prefix (ids in an imported matrix need not match library ids),
  # then by id.
  lib_lookup <- function(met_id, met_name, kind) {
    if (is.null(ctx$mets)) return(NULL)
    if (identical(kind, "parent")) return(parent_met)
    suffix <- function(x) trimws(sub("^\\S+", "", x))
    want <- suffix(met_name)
    hit <- integer(0)
    if (nzchar(want)) {
      lib_suf <- vapply(ctx$mets, function(m) if (identical(m$kind, "parent")) "" else suffix(m$name %||% ""), "")
      hit <- which(lib_suf == want)
    }
    if (length(hit) == 0) hit <- which(vapply(ctx$mets, function(m) identical(m$id, met_id), TRUE))
    if (length(hit) == 0) NULL else ctx$mets[[hit[1]]]
  }

  # class composition of the degradant signal per condition/time
  comp <- NULL
  if (!is.null(ab) && !is.null(degr) && nrow(degr) > 0) {
    dk <- ab$mat$kind
    deg_rows <- which(!dk %in% c("parent", "internal_standard"))
    X <- as.matrix(ab$mat[deg_rows, ab$study, drop = FALSE]); X[is.na(X)] <- 0
    kinds <- unique(dk[deg_rows])
    per_s <- sapply(kinds, function(k) colSums(X[dk[deg_rows] == k, , drop = FALSE]))
    if (is.null(dim(per_s))) per_s <- matrix(per_s, nrow = 1, dimnames = list(ab$study, kinds))
    tot <- rowSums(per_s)
    share <- sweep(per_s, 1, ifelse(tot > 0, tot, NA), "/") * 100
    sf <- sframe[match(rownames(share), sframe$sample), ]
    key <- if (has_time) paste(sf$condition, sf$time, sep = "\r") else sf$condition
    comp <- do.call(rbind, lapply(split(seq_len(nrow(share)), factor(key, levels = unique(key[order(sf$condition, sf$time)]))), function(ix) {
      data.frame(condition = sf$condition[ix[1]], time = sf$time[ix[1]],
                 t(colMeans(share[ix, , drop = FALSE], na.rm = TRUE)), check.names = FALSE)
    }))
  }

  # ---- numbers used in the narrative ----
  n_lib <- if (!is.null(ctx$mets)) length(ctx$mets) else 0L
  n_identified <- if (!is.null(pm)) sum(pm$kind != "internal_standard" & pm$n_detected > 0) else 0L
  n_reported_deg <- if (!is.null(degr)) nrow(degr) else 0L
  n_ms2 <- length(intersect(ms2_ids, rep_mets$met_id %||% character(0)))
  last_t <- if (has_time) max(times) else NA
  fmt_t <- function(t) paste0(format(t, trim = TRUE), " ", tunit)
  kin_line <- function(r) {
    if (is.na(r$k)) return(sprintf("%s: not determined", r$Condition))
    if (!is.finite(r$t_half) || r$t_half > 10 * r$t_max)
      return(sprintf("%s: no measurable first-order loss over %s (stable)", r$Condition, fmt_t(r$t_max)))
    sprintf("%s: k = %s /%s, t1/2 = %s %s (R%s = %.3f)", r$Condition, .reg_sig(r$k, 3), tunit,
            .reg_sig(r$t_half, 3), tunit, .REG_SUP2, r$r2)
  }
  end_state <- if (!is.null(cond_tab) && has_time) {
    e <- cond_tab[cond_tab$time == last_t, , drop = FALSE]
    paste(sprintf("%s retained %s%% of the initial parent at %s", e$condition, .reg_fix(e$rem, 1), fmt_t(last_t)), collapse = "; ")
  } else NULL
  major <- if (!is.null(degr) && nrow(degr) > 0) head(degr, 3) else NULL
  major_txt <- if (!is.null(major)) paste(sprintf("%s (%s, up to %s%% of parent)", major$met_name,
    tolower(.reg_kind_label(major$kind)), .reg_sig(major$max_pct, 3)), collapse = "; ") else NULL
  pathway_txt <- NULL
  if (!is.null(comp) && nrow(comp) > 0) {
    cls <- setdiff(names(comp), c("condition", "time"))
    lc <- if (has_time) comp[comp$time == last_t, , drop = FALSE] else comp
    pathway_txt <- paste(vapply(seq_len(nrow(lc)), function(i) {
      v <- unlist(lc[i, cls]); v <- sort(v[!is.na(v)], decreasing = TRUE)
      sprintf("%s: %s", lc$condition[i], paste(sprintf("%s %s%%", tolower(.reg_kind_label(names(v))), .reg_fix(v, 1)), collapse = ", "))
    }, ""), collapse = "; ")
  }
  stable_conds <- if (!is.null(kin)) kin$Condition[!is.na(kin$k) & (!is.finite(kin$t_half) | kin$t_half > 10 * kin$t_max)] else character(0)
  labile_conds <- if (!is.null(kin)) kin$Condition[!is.na(kin$k) & is.finite(kin$t_half) & kin$t_half <= 10 * kin$t_max] else character(0)
  draft_flag <- "[Draft text generated from the processed data; to be reviewed and finalized by the Study Director.]"

  # ===================== 1. Summary =====================
  H1("Summary")
  P(sprintf("This report describes the identification and relative profiling of the metabolites and degradation products of %s%s. ",
            ta, if (nrow(st) > 0) sprintf(" in %d study samples across %d condition%s%s", nrow(st), length(conds),
                                            if (length(conds) == 1) "" else "s",
                                            if (has_time) sprintf(" and %d timepoints (%s)", length(times), paste(format(times, trim = TRUE), collapse = ", ")) else "") else ""),
    sprintf("Metabolites were screened against a theoretical library of %d candidate species and reported when they reached %s%% of the full-length parent in at least one study sample.",
            n_lib, thr))
  key_rows <- data.frame(
    Finding = c("Candidate species in the theoretical library", "Species detected in study samples",
                sprintf("Metabolites reported (>= %s%% of parent)", thr),
                "MS/MS-confirmed reported metabolites", "Study samples analyzed"),
    Result = c(as.character(n_lib), as.character(n_identified), as.character(n_reported_deg),
               if (ms2_ran) as.character(n_ms2) else "MS/MS not acquired", as.character(nrow(st))),
    stringsAsFactors = FALSE)
  if (!is.null(kin)) key_rows <- rbind(key_rows, data.frame(Finding = paste0("Parent half-life, ", kin$Condition),
    Result = vapply(seq_len(nrow(kin)), function(i) {
      r <- kin[i, ]
      if (is.na(r$k)) "not determined"
      else if (!is.finite(r$t_half) || r$t_half > 10 * r$t_max) sprintf("stable over %s", fmt_t(r$t_max))
      else sprintf("%s %s (R%s %.2f)", .reg_sig(r$t_half, 3), tunit, .REG_SUP2, r$r2)
    }, ""), stringsAsFactors = FALSE))
  TAB(key_rows, "Key findings.", widths = c(3.9, 2.37))
  if (!is.null(end_state)) P("Parent depletion: ", end_state, ".")
  if (!is.null(major_txt)) P("Principal metabolites: ", major_txt, ".")
  if (!is.null(pathway_txt)) P(sprintf("Composition of the degradant signal%s: ", if (has_time) paste0(" at ", fmt_t(last_t)) else ""), pathway_txt, ".")
  NOTE(draft_flag)

  # ===================== 2. Introduction and objectives =====================
  H1("Introduction and Objectives")
  P(if (!is.null(chem)) sprintf("%s is a %d-mer %s oligonucleotide. ", ta, chem$n, chem$design) else sprintf("%s is an oligonucleotide test article. ", ta),
    "Oligonucleotide therapeutics are cleared predominantly by nuclease-mediated chain shortening from the 3' and 5' termini and by endonucleolytic cleavage; ",
    "characterizing these products, their abundance relative to the parent and the rate of parent loss supports the safety assessment of circulating and tissue metabolites.")
  P("The objectives of this study were to:")
  BUL(c(sprintf("identify the metabolites and degradation products of %s by liquid chromatography-high resolution mass spectrometry (LC-HRMS);", ta),
        "report each metabolite's abundance relative to the full-length parent and its share of the total drug-related signal;",
        if (has_time) "describe the time course of parent loss and estimate first-order degradation kinetics per condition;" else NULL,
        "provide structural assignments (truncation position and chemistry) for the reported metabolites."))
  P("[Expand per protocol.]")

  # ===================== 3. Test article =====================
  H1("Test Article")
  ta_tab <- data.frame(
    Attribute = c("Name / code", "Lot", "Length", "Sequence (5'->3')", "Chemistry design",
                  "Sugar pattern", "Backbone", "Modified nucleobases", "Conjugates",
                  "Molecular formula (free acid)", "Monoisotopic mass (Da)", "Average mass (Da)",
                  "Internal standard"),
    Value = c(ta, .reg_get(details, "test_article_lot"),
              if (is.null(chem)) "[length]" else sprintf("%d nucleotides", chem$n),
              chem$seq_txt %||% "[sequence]", chem$design %||% "[design]",
              chem$sugar_txt %||% "[sugars]", chem$backbone %||% "[backbone]",
              chem$mod_txt %||% "[modified bases]", chem$conj_txt %||% "[conjugates]",
              mass_info$formula_str %||% "[formula]",
              if (is.null(mass_info$mono_mass)) "[mass]" else sprintf("%.4f", mass_info$mono_mass),
              if (is.null(mass_info$avg_mass) || is.na(mass_info$avg_mass)) "[mass]" else sprintf("%.2f", mass_info$avg_mass),
              .reg_get(details, "is_name", if (!is.null(is_rows) && nrow(is_rows) > 0) is_rows$met_name[1] else "none")),
    stringsAsFactors = FALSE)
  TAB(ta_tab, "Test article identity and chemistry.", widths = c(2.0, 4.27))
  P("In the sequence, mC denotes 5-methylcytidine. The per-position chemistry is listed in Appendix 1. ",
    "Molecular formula and masses are calculated for the free acid from the per-residue chemistry dictionary of OligoMet Profiler.")

  # ===================== 4. Study design =====================
  H1("Study Design")
  sys_tab <- data.frame(
    Item = c("Test system", "Route / incubation", "Matrix", "Conditions", "Timepoints", "Replicates"),
    Description = c(.reg_get(details, "species"), .reg_get(details, "route"),
                    .reg_get(details, "matrix", if (any(nzchar(st$matrix))) paste(unique(st$matrix[nzchar(st$matrix)]), collapse = ", ") else "[matrix]"),
                    if (length(conds)) paste(conds, collapse = ", ") else "-",
                    if (has_time) paste(paste(format(times, trim = TRUE), collapse = ", "), tunit) else "not applicable",
                    if (nrow(st) > 0) {
                      reps <- table(paste(st$condition, st$timepoint))
                      if (length(unique(reps)) == 1) sprintf("n = %d per condition%s", reps[1], if (has_time) " and timepoint" else "")
                      else sprintf("n = %d-%d per condition", min(reps), max(reps))
                    } else "-"),
    stringsAsFactors = FALSE)
  TAB(sys_tab, "Test system and design.", widths = c(1.8, 4.47))
  if (nrow(sframe) > 0) {
    tl <- .reg_type_label(ifelse(sframe$study, "", sframe$sample_type))
    inv <- as.data.frame(table(factor(tl, levels = c("Study sample", "Calibration standard", "Quality control", "Reagent blank", "Matrix blank"))),
                         stringsAsFactors = FALSE)
    inv <- inv[inv$Freq > 0, , drop = FALSE]
    names(inv) <- c("Sample type", "Number of injections")
    inv$Use <- ifelse(inv$`Sample type` == "Study sample", "metabolite profiling, degradation and kinetics",
                      ifelse(inv$`Sample type` == "Calibration standard", "calibration (where absolute quantitation applied)",
                             ifelse(inv$`Sample type` == "Quality control", "assay performance", "carry-over and interference check")))
    TAB(inv, "Sample inventory. Only study samples contribute to the results in Section 6.", widths = c(1.9, 1.5, 2.87))
  }
  P("[Dosing, sample collection, storage and shipment per protocol.]")

  # ===================== 5. Materials and methods =====================
  adducts <- params$adducts %||% c("H", "Na", "K", "NH4")
  mods <- params$mods %||% character(0)
  H1("Materials and Methods")
  H2("Sample Preparation and LC-HRMS Analysis")
  P("Samples were analyzed by ion-pairing reversed-phase liquid chromatography coupled to high-resolution mass spectrometry in negative electrospray ionization, with data-dependent MS/MS where acquired. ",
    "[Extraction procedure, LC column, mobile phases, gradient, mass spectrometer and acquisition settings per SOP.]")
  H2("Theoretical Metabolite Library")
  lib_counts <- if (!is.null(ctx$mets) && length(ctx$mets) > 0) {
    tb <- table(vapply(ctx$mets, function(m) m$kind %||% "unknown", ""))
    paste(sprintf("%d %s", as.integer(tb), tolower(.reg_kind_label(names(tb)))), collapse = ", ")
  } else NULL
  P("A theoretical library of candidate metabolites was generated from the parent sequence and chemistry, covering stepwise 3'- and 5'-exonuclease truncation, endonucleolytic cleavage and, where applicable, conjugate catabolism. ",
    "Each candidate carries an exact elemental formula, from which isotope distributions and multiply charged ion series were computed. ",
    if (!is.null(lib_counts)) sprintf("The library for this study contained %d species: %s.", n_lib, lib_counts) else "")
  H2("Identification Criteria")
  methods <- data.frame(
    Parameter = c("Mass accuracy tolerance", "Charge states searched", "Adducts considered",
                  "Phosphorothioate desulfurization (PS->PO) searched", "In-place modifications searched",
                  "Isotope pattern", "MS/MS confirmation", "Metabolite reporting threshold"),
    Setting = c(sprintf("+/- %s ppm", params$ppm_tol %||% 10),
                sprintf("z = %s to %s", params$z_min %||% 3, params$z_max %||% 12),
                paste(adducts, collapse = ", "),
                sprintf("up to %s per molecule", params$max_oxid %||% 6),
                if (length(mods) == 0) "none" else paste(mods, collapse = ", "),
                "theoretical isotope envelope compared with the observed charge-state envelope",
                if (ms2_ran) "fragment-ion matching against the predicted a-B/w/c/y ladders" else "not acquired for this dataset",
                sprintf(">= %s%% of the full-length parent in at least one study sample", thr)),
    stringsAsFactors = FALSE)
  TAB(methods, "Data-processing and identification parameters.", widths = c(2.4, 3.87))

  # quantitation methodology
  xic_on <- isTRUE(params$xic_enabled)
  signal_basis <- params$signal_basis %||% "raw"
  basis_lab <- switch(signal_basis,
    raw = "raw peak response (area where available, else apex intensity)",
    bcorr = "blank-corrected peak response",
    isr = "internal-standard-normalized response ratio (analyte area / IS area)",
    bcorr_isr = "blank-corrected and internal-standard-normalized response",
    signal_basis)
  agg_lab <- if (xic_on) {
    sprintf("summed extracted-ion chromatogram: the top %s charge states and %s isotopes are co-integrated over one shared retention-time window per metabolite (summation of ions)",
            params$xic_n_charges %||% 5, params$xic_n_iso %||% 5)
  } else {
    "collapsed to one value per metabolite and sample, taken as the maximum response across charge states, adducts and desulfurization levels (the dominant ionization form)"
  }
  gis <- params$global_is %||% "none"
  is_lab <- if (!is.null(gis) && nzchar(gis) && !identical(gis, "none"))
    sprintf("response normalized to the internal standard (%s) in the same sample", gis)
  else "no internal-standard normalization applied"
  quant_methods <- data.frame(
    Aspect = c("Ionization", "Adduct handling", "Charge-state handling", "Signal aggregation",
               "Signal basis", "Internal standard", "Relative abundance", "Degradation and kinetics"),
    Approach = c(
      "negative-ion electrospray; analytes detected as multiply charged [M-zH]z- ions",
      sprintf("the deprotonated [M-zH] form is the dominant and quantified species; %s adducts are the same metabolite with a cation in place of a proton and are used for confirmation, not counted as separate species",
              paste(setdiff(adducts, "H"), collapse = "/")),
      if (xic_on) sprintf("charge-state envelope (z %s to %s) co-summed per metabolite", params$z_min %||% 3, params$z_max %||% 12)
      else sprintf("charge-state envelope (z %s to %s); the dominant charge state carries the reported response, with cross-charge-state mass agreement used as a confirmation check",
                   params$z_min %||% 3, params$z_max %||% 12),
      agg_lab, basis_lab, is_lab,
      "metabolite response / parent response in the same sample x 100 (% of parent); metabolite response / summed response of all drug-related species x 100 (% of total drug-related signal)",
      if (has_time) sprintf("percent parent remaining relative to the mean of the earliest timepoint of each condition; first-order rate constant k from linear regression of ln(%% remaining) on time; t1/2 = ln 2 / k") else "percent degradation = summed degradant response / total drug-related response x 100"),
    stringsAsFactors = FALSE)
  quant_narrative <- paste0(
    "Each metabolite's response was ", agg_lab, ". ",
    "Adduct forms of a metabolite were not treated as separate species: the deprotonated [M-zH]^z- envelope is the dominant and quantified form, and sodium, potassium and ammonium adducts, when observed, were used to corroborate the identification. ",
    "Signal basis for the reported values was ", basis_lab, "; ", is_lab, ". ",
    "Abundances are relative (equal ionization response assumed for parent and metabolites) unless a calibration curve is reported. ",
    "Calibration standards, quality controls and blanks were excluded from all relative-abundance, degradation and kinetic calculations.")
  H2("Signal aggregation and quantitation")
  P(quant_narrative)
  TAB(quant_methods, "Quantitation and calculation methodology.", widths = c(1.7, 4.57))
  H2("Software")
  P(sprintf("Data processing, library generation, matching and report generation: OligoMet Profiler version %s, running on %s. ",
            .reg_version(), R.version.string),
    "[Instrument acquisition and vendor processing software per SOP.]")

  # ===================== 6. Results =====================
  H1("Results")
  H2("Metabolite Identification")
  if (!is.null(rep_mets) && nrow(rep_mets) > 0) {
    has_ppm <- !is.null(matches) && "ppm_error" %in% names(matches) && any(!is.na(matches$ppm_error))
    has_z <- !is.null(matches) && "z" %in% names(matches) && any(!is.na(matches$z))
    stud <- if (!is.null(matches)) matches[matches$sample %in% ab$study, , drop = FALSE] else NULL
    id_tab <- data.frame(
      ID = rep_mets$met_id,
      Metabolite = rep_mets$met_name,
      Class = .reg_kind_label(rep_mets$kind),
      Length = vapply(seq_len(nrow(rep_mets)), function(i) {
        lm <- lib_lookup(rep_mets$met_id[i], rep_mets$met_name[i], rep_mets$kind[i]); if (is.null(lm)) "-" else paste0(lm$n, "-mer") }, ""),
      Formula = vapply(seq_len(nrow(rep_mets)), function(i) {
        lm <- lib_lookup(rep_mets$met_id[i], rep_mets$met_name[i], rep_mets$kind[i]); if (is.null(lm)) return("-")
        tryCatch(metabolite_mass_info(lm, dict)$formula_str, error = function(e) "-") }, ""),
      `Mono. mass (Da)` = vapply(seq_len(nrow(rep_mets)), function(i) {
        lm <- lib_lookup(rep_mets$met_id[i], rep_mets$met_name[i], rep_mets$kind[i]); if (is.null(lm)) return("-")
        tryCatch(sprintf("%.3f", metabolite_mass_info(lm, dict)$mono_mass), error = function(e) "-") }, ""),
      check.names = FALSE, stringsAsFactors = FALSE)
    if (has_ppm) id_tab$`Median error (ppm)` <- vapply(rep_mets$met_id, function(id) {
      v <- stud$ppm_error[stud$met_id == id]; if (length(v) == 0 || all(is.na(v))) "-" else sprintf("%.1f", stats::median(v, na.rm = TRUE)) }, "")
    if (has_z) id_tab$`Charge states` <- vapply(rep_mets$met_id, function(id) {
      z <- sort(unique(abs(stud$z[stud$met_id == id & !is.na(stud$z)]))); if (length(z) == 0) "-" else paste(z, collapse = ", ") }, "")
    id_tab$MS2 <- if (!ms2_ran) "n/a" else ifelse(rep_mets$met_id %in% ms2_ids, "confirmed", "not confirmed")
    TAB(id_tab, sprintf("Reported metabolites and their structural assignments (parent plus species >= %s%% of parent).", thr))
    if (!has_ppm) NOTE("Accurate-mass error and charge-state evidence are not available for results imported from a pre-processed data matrix; identities are as assigned in the source matrix and matched to the theoretical library by name.")
    if (!ms2_ran) NOTE("MS/MS confirmation was not acquired or not processed for this dataset; assignments rest on accurate mass and isotope/charge-envelope agreement.")
  } else {
    P("No metabolite abundance data were available for this run.")
  }

  H2("Relative Abundance")
  if (!is.null(rep_mets) && nrow(rep_mets) > 0) {
    ab_tab <- data.frame(
      Metabolite = rep_mets$met_name,
      `Detected (study samples)` = sprintf("%d / %d", rep_mets$n_detected, rep_mets$n_study),
      `Max % of parent` = ifelse(rep_mets$kind == "parent", "100 (reference)", .reg_sig(rep_mets$max_pct, 3)),
      `Max % of total drug-related` = .reg_sig(rep_mets$max_pct_total, 3),
      `Observed at` = rep_mets$at,
      check.names = FALSE, stringsAsFactors = FALSE)
    TAB(ab_tab, "Maximum abundance of each reported species in study samples.", widths = c(1.6, 1.0, 1.0, 1.15, 1.52))
    excl <- if (!is.null(pm)) pm[!pm$reported & pm$kind != "internal_standard" & pm$n_detected > 0, , drop = FALSE] else NULL
    if (!is.null(excl) && nrow(excl) > 0)
      P(sprintf("%d further species were detected below the %s%% reporting threshold (%s) and are not reported.",
                nrow(excl), thr, paste(excl$met_name, collapse = ", ")))
    if (has_time && !is.null(degr) && nrow(degr) > 0) {
      for (cn in conds) {
        cols <- ab$study[sframe$condition[match(ab$study, sframe$sample)] == cn]
        ct <- sframe$time[match(cols, sframe$sample)]
        tpts <- sort(unique(ct[!is.na(ct)]))
        if (length(tpts) == 0) next
        ft <- data.frame(Metabolite = degr$met_name, stringsAsFactors = FALSE)
        for (t in tpts) {
          sc <- cols[ct == t & !is.na(ct)]
          v <- rowMeans(ab$pct_parent[match(degr$met_id, ab$mat$met_id), sc, drop = FALSE], na.rm = TRUE)
          ft[[fmt_t(t)]] <- .reg_sig(v, 3)
        }
        TAB(ft, sprintf("Formation of metabolites over time, %s (mean %% of parent, n = %d per timepoint).", cn,
                        round(length(cols) / length(tpts))))
      }
    }
  }

  if (!is.null(cond_tab)) {
    H2(if (has_time) "Parent Depletion and Degradation Kinetics" else "Degradation by Condition")
    dt <- data.frame(Condition = cond_tab$condition, stringsAsFactors = FALSE)
    if (has_time) dt[[sprintf("Time (%s)", tunit)]] <- format(cond_tab$time, trim = TRUE)
    dt$n <- as.character(cond_tab$n)
    dt$`Mean parent response` <- .reg_fix(cond_tab$parent, 0)
    msd <- function(m, s) ifelse(is.na(s), .reg_fix(m, 1), sprintf("%s +/- %s", .reg_fix(m, 1), .reg_fix(s, 1)))
    if (has_time) dt$`Parent remaining (% of t0, mean +/- SD)` <- msd(cond_tab$rem, cond_tab$rem_sd)
    dt$`Degradation (% of total, mean +/- SD)` <- msd(cond_tab$deg, cond_tab$deg_sd)
    TAB(dt, if (has_time) "Parent remaining and percent degradation by condition and timepoint." else "Percent degradation by condition.")
    if (!is.null(kin) && nrow(kin) > 0) {
      kt <- data.frame(Condition = kin$Condition,
        `k (1/h)` = ifelse(is.na(kin$k), "-", .reg_sig(kin$k, 3)),
        `t1/2` = vapply(seq_len(nrow(kin)), function(i) {
          r <- kin[i, ]
          if (is.na(r$k)) "not determined" else if (!is.finite(r$t_half) || r$t_half > 10 * r$t_max) sprintf("> %s (stable)", fmt_t(10 * r$t_max))
          else fmt_t(.reg_sig(r$t_half, 3)) }, ""),
        `R squared` = ifelse(is.na(kin$r2), "-", sprintf("%.3f", kin$r2)),
        `Points fitted` = as.character(kin$n), check.names = FALSE, stringsAsFactors = FALSE)
      names(kt)[2] <- sprintf("k (1/%s)", tunit)
      TAB(kt, "First-order kinetics of parent loss.", widths = c(2.0, 1.1, 1.4, 0.85, 0.92))
      P(paste(vapply(seq_len(nrow(kin)), function(i) kin_line(kin[i, ]), ""), collapse = ". "), ".")
    }
    if (!is.null(ctx$degradation)) {
      FIG(tryCatch(plot_degradation_by_condition(ctx$degradation), error = function(e) NULL),
          "Percent degradation over time by condition (line: mean; error bars: SD; points: replicates).")
      FIG(tryCatch(plot_degradation_composition(ctx$degradation, "condition"), error = function(e) NULL),
          "Composition of the drug-related signal by metabolite class, per condition and timepoint (mean of replicates).", h = 3.6)
    }
  }

  if (!is.null(comp) && nrow(comp) > 0) {
    H2("Biotransformation Pathways")
    cls <- setdiff(names(comp), c("condition", "time"))
    ct <- data.frame(Condition = comp$condition, stringsAsFactors = FALSE)
    if (has_time) ct[[sprintf("Time (%s)", tunit)]] <- format(comp$time, trim = TRUE)
    for (k in cls) ct[[.reg_kind_label(k)]] <- .reg_fix(comp[[k]], 1)
    TAB(ct, "Share of the degradant signal by pathway (mean % of summed degradant response).")
    if (!is.null(pathway_txt)) P(sprintf("At %s the degradant signal was distributed as follows. ", if (has_time) fmt_t(last_t) else "each condition"), pathway_txt, ".")
  }

  H2("Metabolite Characterization")
  char_ids <- if (!is.null(rep_mets)) head(rep_mets, params$max_characterized %||% 15) else NULL
  raw_ms <- !is.null(matches) && "obs_mz" %in% names(matches)
  ev <- if (raw_ms && !is.null(char_ids) && nrow(char_ids) > 0) {
    tryCatch(collect_metabolite_evidence(char_ids$met_id, matches, ctx$features, ctx$mets, dict,
                                         study_samples = ab$study, ms2_confirmations = ms2c,
                                         ms2_spectra = ctx$ms2_spectra, h_offset = params$h_offset %||% 0,
                                         use_envipat = params$use_envipat %||% TRUE,
                                         ppm = params$ppm_tol %||% 10, extract = !isFALSE(ctx$extract_raw)),
             error = function(e) list())
  } else list()
  n_raw <- sum(vapply(ev, function(e) !is.null(e$spectrum), logical(1)))
  if (!is.null(char_ids) && nrow(char_ids) > 0) {
    P("Each reported species is documented below: its sequence and position in the parent, molecular formula and masses, the precursor ion and charge states used for identification, the originating sample (the study sample with the strongest signal), ",
      if (raw_ms) "the extracted-ion chromatogram and apex MS1 isotope pattern re-read from that sample's raw data file, the charge-state envelope, " else "the theoretical isotope distribution, ",
      "the MS/MS spectrum where acquired, and a structure map against the parent (retained residues drawn with their sugar chemistry, removed residues greyed out; solid linkage: phosphorothioate, dotted: phosphodiester).")
    if (raw_ms && n_raw == 0 && length(ev) > 0)
      NOTE("The raw data files were not available when this report was generated (or Python could not read them), so chromatograms and acquired MS1 spectra are omitted; the matched m/z values and theoretical patterns are shown.")
    if (!raw_ms)
      NOTE("Results were imported from a pre-processed data matrix, so observed m/z values, charge states and spectra are not available; theoretical values are shown for each species.")
    parent_bases <- ctx$spec$bases
    z_all <- seq(params$z_min %||% 3, params$z_max %||% 12)
    for (i in seq_len(nrow(char_ids))) {
      r <- char_ids[i, ]
      lm <- lib_lookup(r$met_id, r$met_name, r$kind)
      e <- ev[[r$met_id]]
      info <- if (!is.null(lm)) tryCatch(metabolite_mass_info(lm, dict), error = function(x) NULL) else NULL
      if (i > 1) BREAK()
      H3(r$met_name)
      pos_txt <- ""
      if (!is.null(lm) && !is.null(parent_bases)) {
        L <- length(lm$bases)
        hit <- which(vapply(seq_len(max(0, length(parent_bases) - L + 1)), function(st)
          identical(parent_bases[st:(st + L - 1)], lm$bases), logical(1)))
        if (length(hit) == 1) pos_txt <- sprintf(" (positions %d-%d of the parent)", hit, hit + L - 1)
      }
      # originating sample: from the match evidence, else the abundance matrix
      orig <- if (!is.null(e)) e$row$sample else if (!is.null(ab)) {
        rr <- match(r$met_id, ab$mat$met_id)
        v <- as.numeric(ab$mat[rr, ab$study]); if (all(is.na(v))) NA else ab$study[which.max(v)]
      } else NA
      orig_txt <- if (is.na(orig)) "-" else {
        sf <- sframe[sframe$sample == orig, ]
        ctxt <- paste(c(if (nrow(sf) && nzchar(sf$condition)) sf$condition, if (nrow(sf) && nzchar(sf$timepoint)) paste0("t = ", fmt_t(sf$timepoint))), collapse = ", ")
        sig_v <- if (!is.null(ab)) as.numeric(ab$mat[match(r$met_id, ab$mat$met_id), orig]) else NA
        sprintf("%s%s%s", orig, if (nzchar(ctxt)) paste0(" (", ctxt, ")") else "",
                if (!is.na(sig_v)) paste0("; response ", .reg_fix(sig_v, 0)) else "")
      }
      kv <- list(
        Class = .reg_kind_label(r$kind),
        `Sequence (5'->3')` = if (!is.null(lm)) paste0(paste(.reg_base_letter(lm$bases), collapse = ""), pos_txt) else "-",
        Length = if (!is.null(lm)) sprintf("%d nucleotides", lm$n) else "-",
        `Molecular formula (free acid)` = info$formula_str %||% "-",
        `Monoisotopic / average mass (Da)` = if (!is.null(info)) sprintf("%.4f / %.2f", info$mono_mass, info$avg_mass) else "-")
      if (!is.null(e)) {
        rw <- e$row
        if (identical(rw$adduct, "H") && rw$k_oxid == 0 && !is.null(info)) {
          meas <- rw$obs_mz * rw$z + rw$z * .PROTON - (params$h_offset %||% 0)
          kv[["Measured monoisotopic mass (Da)"]] <- sprintf("%.4f (%+.1f ppm)", meas, (meas - info$mono_mass) / info$mono_mass * 1e6)
        }
        kv[["Precursor ion"]] <- sprintf("%s: theoretical m/z %.4f, observed m/z %.4f (%+.1f ppm)%s",
          .ev_ion_label(rw$z, rw$adduct), rw$theo_mz, rw$obs_mz, rw$ppm_error,
          if (rw$k_oxid > 0) sprintf("; %d PS->PO desulfurization", rw$k_oxid) else "")
        kv[["Charge states observed"]] <- paste(sprintf("%d- (m/z %.4f)", e$envelope$z, e$envelope$obs_mz), collapse = ", ")
        if (!is.na(rw$rt)) kv[["Retention time (apex)"]] <- sprintf("%.2f min", rw$rt)
        if (!is.null(rw$iso_fit) && !is.na(rw$iso_fit)) kv[["Isotope pattern fit"]] <- sprintf("%.2f (0-1; cosine similarity x fraction of isotopes matched)", rw$iso_fit)
        if (nrow(e$envelope) >= 2) {
          ec <- tryCatch(envelope_consistency(matches[matches$met_id == r$met_id & matches$sample == rw$sample &
                                                      matches$adduct == rw$adduct & matches$k_oxid == rw$k_oxid, , drop = FALSE],
                                              h_offset = params$h_offset %||% 0), error = function(x) NULL)
          if (!is.null(ec) && nrow(ec) > 0 && !is.na(ec$mass_cv_ppm[1]))
            kv[["Charge-state mass agreement"]] <- sprintf("%d charge states, neutral-mass CV %.2f ppm", ec$n_z[1], ec$mass_cv_ppm[1])
        }
      } else if (!is.null(info)) {
        zz <- z_all[((info$mono_mass - z_all * .PROTON) / z_all) >= 400 & ((info$mono_mass - z_all * .PROTON) / z_all) <= 2000]
        if (length(zz)) kv[["Theoretical precursor m/z"]] <- paste(sprintf("%d- %.4f", zz, (info$mono_mass - zz * .PROTON) / zz), collapse = "; ")
      }
      kv[["Originating sample"]] <- orig_txt
      kv[["Detected in"]] <- sprintf("%d of %d study samples", r$n_detected, r$n_study)
      kv[["MS/MS"]] <- if (!is.null(e) && !is.null(e$ms2)) {
        cf <- e$ms2$conf
        sprintf("%s: score %.1f, %d fragment ion(s) matched (%.0f%% sequence coverage), %d PS diagnostic ion(s); %s, sample %s",
                if (isTRUE(cf$confident)) "confirmed" else "acquired, not confirmed",
                cf$confirmation_score, cf$n_frag_matches, 100 * cf$coverage, cf$n_diagnostics %||% 0L,
                .ev_ion_label(cf$z, cf$adduct), cf$sample)
      } else if (ms2_ran) "no MS/MS spectrum acquired for this precursor" else "not acquired"
      kvdf <- data.frame(Item = names(kv), Value = unlist(kv, use.names = FALSE), stringsAsFactors = FALSE)
      TAB(kvdf, sprintf("Identification evidence for %s.", r$met_name), widths = c(1.9, 4.37), header = FALSE)

      if (!is.null(e)) {
        FIG(plot_evidence_xic(e), sprintf("Extracted-ion chromatogram of %s, %s, in %s.", r$met_name, .ev_ion_label(e$row$z, e$row$adduct), e$row$sample), h = 2.2)
        FIG(plot_evidence_isotopes(e, params$ppm_tol %||% 10), sprintf("MS1 isotope pattern of %s at the chromatographic apex, observed against theoretical.", r$met_name), h = 2.7)
        FIG(plot_evidence_envelope(e), sprintf("Charge-state envelope of %s in %s.", r$met_name, e$row$sample), h = 2.2)
        if (!is.null(e$ms2))
          FIG(tryCatch(plot_mirror_spectrum(e$ms2$data, title = sprintf("%s, %s precursor, sample %s", r$met_name,
                                                                        .ev_ion_label(e$ms2$conf$z, e$ms2$conf$adduct), e$ms2$conf$sample)),
                       error = function(x) NULL),
              sprintf("MS/MS spectrum of %s: acquired (top) against predicted fragment ions (bottom).", r$met_name), h = 3.2)
      } else if (!is.null(info)) {
        zz <- z_all[((info$mono_mass - z_all * .PROTON) / z_all) >= 400]
        zt <- if (length(zz)) zz[which.min(abs((info$mono_mass - zz * .PROTON) / zz - 1100))] else NA
        if (!is.na(zt)) {
          te <- list(row = data.frame(z = zt, adduct = "H", theo_mz = (info$mono_mass - zt * .PROTON) / zt, sample = "-"),
                     theo_iso = .ev_nominal_isotopes(info$formula_str, zt, (info$mono_mass - zt * .PROTON) / zt,
                                                     params$h_offset %||% 0, params$use_envipat %||% TRUE))
          FIG(plot_evidence_isotopes(te), sprintf("Theoretical isotope distribution of %s, %s.", r$met_name, .ev_ion_label(zt)), h = 2.5)
        }
      }
      if (!is.null(lm) && !is.null(ctx$spec)) {
        svg <- tryCatch(annotate_metabolite_svg(lm, parent = ctx$spec, dict = dict), error = function(x) NULL)
        if (!is.null(svg)) {
          if (nzchar(lib_prefix) && !identical(lib_prefix, ta) && !grepl("^\\[", ta)) svg <- gsub(lib_prefix, ta, svg, fixed = TRUE)
          SVG(svg, sprintf("Structure map: %s (%s).", r$met_name, tolower(.reg_kind_label(r$kind))))
        }
      }
    }
  } else {
    P("No reported metabolites to characterize.")
  }

  # ===================== 7. Discussion =====================
  H1("Discussion")
  disc <- character(0)
  if (length(labile_conds) > 0)
    disc <- c(disc, sprintf("%s was degraded in %s, with first-order loss of the parent over the incubation (%s).", ta,
                            paste(labile_conds, collapse = " and "),
                            paste(vapply(which(kin$Condition %in% labile_conds), function(i) kin_line(kin[i, ]), ""), collapse = "; ")))
  if (length(stable_conds) > 0)
    disc <- c(disc, sprintf("No meaningful loss of parent occurred in %s, consistent with enzyme-mediated (nuclease) rather than chemical degradation.", paste(stable_conds, collapse = " and ")))
  if (!is.null(major_txt))
    disc <- c(disc, sprintf("The principal products were %s.", major_txt))
  if (!is.null(comp) && nrow(comp) > 0) {
    cls <- setdiff(names(comp), c("condition", "time"))
    lc <- if (has_time) comp[comp$time == last_t, , drop = FALSE] else comp
    if (length(labile_conds) > 0) lc <- lc[lc$condition %in% labile_conds, , drop = FALSE]
    if (nrow(lc) > 0) {
      v <- colMeans(lc[, cls, drop = FALSE], na.rm = TRUE); top <- names(which.max(v))
      disc <- c(disc, sprintf("Chain shortening was dominated by %s activity (%s%% of the degradant signal), the pattern expected for a phosphorothioate oligonucleotide with nuclease-resistant modified termini being trimmed residue by residue.",
                              sub(" product$", "", tolower(.reg_kind_label(top))), .reg_fix(v[top], 0)))
    }
  }
  if (!ms2_ran) disc <- c(disc, "Assignments should be confirmed by MS/MS before the report is finalized.")
  if (length(disc) == 0) disc <- "[Interpretation of the matrix-specific biotransformation pathways and of any metabolite of toxicological interest.]"
  for (d in disc) P(d)
  P("[Assess any metabolite exceeding 10% of total drug-related material against regulatory guidance on metabolites in safety testing.]")
  NOTE(draft_flag)

  # ===================== 8. Conclusions =====================
  H1("Conclusions")
  concl <- c(
    if (n_reported_deg > 0) sprintf("%d metabolite%s of %s reached the %s%% reporting threshold%s.", n_reported_deg,
                                     if (n_reported_deg == 1) "" else "s", ta, thr,
                                     if (!is.null(major)) sprintf("; the most abundant was %s at up to %s%% of parent", major$met_name[1], .reg_sig(major$max_pct[1], 3)) else "") else NULL,
    if (!is.null(end_state)) paste0(end_state, ".") else NULL,
    if (!is.null(kin)) paste0(paste(vapply(seq_len(nrow(kin)), function(i) kin_line(kin[i, ]), ""), collapse = "; "), ".") else NULL)
  if (length(concl) == 0) concl <- "[Summary conclusion.]"
  BUL(concl)
  NOTE(draft_flag)

  # ===================== 9. Appendices =====================
  H1("Appendices")
  if (!is.null(chem)) {
    H2("Sequence and Chemistry by Position")
    TAB(chem$per_pos, sprintf("Per-position chemistry of %s.", ta), widths = c(0.8, 1.9, 2.0, 1.57))
  }
  if (!is.null(per_sample) && nrow(per_sample) > 0) {
    H2("Individual Study-Sample Results")
    ist <- data.frame(Sample = per_sample$sample, Condition = per_sample$condition, stringsAsFactors = FALSE)
    if (has_time) ist[[sprintf("Time (%s)", tunit)]] <- format(per_sample$time, trim = TRUE)
    ist$`Parent response` <- .reg_fix(per_sample$parent, 0)
    ist$`Total drug-related response` <- .reg_fix(per_sample$total, 0)
    if (has_time) ist$`Parent remaining (%)` <- .reg_fix(per_sample$remaining, 1)
    ist$`Degradation (%)` <- .reg_fix(per_sample$pct_deg, 1)
    TAB(ist, sprintf("Per-sample parent and drug-related response (%s).", if (!is.null(ab)) ab$signal else "signal"),
        widths = if (has_time) c(1.5, 1.1, 0.55, 0.8, 0.85, 0.75, 0.72) else c(1.7, 1.3, 1.0, 1.2, 1.07))
  }
  if (nrow(sframe) > 0 && any(!sframe$study) && !is.null(matches)) {
    H2("Calibration, QC and Blank Samples")
    ns <- sframe[!sframe$study, , drop = FALSE]
    sig <- .reg_signal_col(matches)
    mat <- build_abundance_matrix(matches, signal_col = sig)
    pr <- which(mat$kind == "parent")[1]
    nt <- data.frame(Sample = ns$sample, Type = .reg_type_label(ns$sample_type),
                     `Parent response` = if (is.na(pr)) "-" else .reg_fix(as.numeric(mat[pr, ns$sample]), 0),
                     check.names = FALSE, stringsAsFactors = FALSE)
    TAB(nt, "Non-study injections (excluded from all results in Section 6).", widths = c(2.0, 2.0, 2.27))
  }
  H2("Deviations")
  P("[None / list protocol and SOP deviations and their impact.]")

  list(
    details = details, threshold = threshold, title = title, test_article = ta,
    status = status, report_date = report_date,
    cover = c(
      `Study number` = .reg_get(details, "study_number"),
      Sponsor = paste0(.reg_get(details, "sponsor"), { a <- .reg_get(details, "sponsor_address", ""); if (nzchar(a)) paste0(", ", a) else "" }),
      `Sponsor project code` = .reg_get(details, "project_code"),
      `Test facility` = paste0(.reg_get(details, "cro"), { a <- .reg_get(details, "cro_address", ""); if (nzchar(a)) paste0(", ", a) else "" }),
      `Test article` = paste0(ta, " (lot ", .reg_get(details, "test_article_lot"), ")"),
      `Test system` = paste0(.reg_get(details, "species"), "; ", .reg_get(details, "route")),
      Protocol = .reg_get(details, "protocol"),
      `Study Director` = .reg_get(details, "study_director"),
      `Report status` = status,
      `Report date` = report_date),
    signatures = data.frame(
      Role = c("Study Director", "Principal Investigator (bioanalysis)", "Quality Assurance", "Sponsor representative"),
      Name = c(.reg_get(details, "study_director"), "[name]", "[name]", "[name]"),
      Signature = "", Date = "", stringsAsFactors = FALSE),
    abbreviations = data.frame(
      Abbreviation = c("DNA", "HRMS", "IS", "LC", "m/z", "mC", "MOE", "MS/MS", "PO", "ppm", "PS", "SD", "t1/2"),
      Definition = c("2'-deoxyribonucleotide", "high-resolution mass spectrometry", "internal standard",
                     "liquid chromatography", "mass-to-charge ratio", "5-methylcytidine",
                     "2'-O-(2-methoxyethyl)", "tandem mass spectrometry", "phosphodiester", "parts per million",
                     "phosphorothioate", "standard deviation", "half-life"),
      stringsAsFactors = FALSE),
    blocks = blocks,
    # kept for callers that read the model directly
    test_article_table = ta_tab, methods = methods, quant_methods = quant_methods,
    quant_narrative = quant_narrative,
    summary = list(n_identified = n_identified, n_ms2 = n_ms2, n_floor = n_reported_deg,
                   matrices = unique(st$matrix[nzchar(st$matrix)])))
}

# Render a spec as a 5'->3' triplet string (sugar+base, linkage between).
.reg_sequence_string <- function(spec) {
  n <- spec$n
  res <- vapply(seq_len(n), function(i) paste0(spec$sugars[i], spec$bases[i]), "")
  lk <- spec$linkages
  out <- res[1]
  for (i in seq_len(n - 1)) out <- paste0(out, "-", if (!is.na(lk[i])) lk[i] else "", res[i + 1])
  c5 <- spec$conj5 %||% "none"; c3 <- spec$conj3 %||% "none"
  if (!identical(c5, "none")) out <- paste0(c5, "-", out)
  if (!identical(c3, "none")) out <- paste0(out, "-", c3)
  out
}

## ---- HTML renderer ---------------------------------------------------------
.reg_html_table <- function(df, caption = NULL, header = TRUE, cls = "") {
  if (is.null(df) || nrow(df) == 0) return("")
  th <- if (header) paste0("<thead><tr>", paste0("<th>", .reg_xml_escape(names(df)), "</th>", collapse = ""), "</tr></thead>") else ""
  rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", .reg_xml_escape(r), "</td>", collapse = ""), "</tr>"))
  cap <- if (!is.null(caption)) paste0("<div class='tbl-cap'>", .reg_xml_escape(caption), "</div>") else ""
  paste0(cap, "<table class='", cls, "'>", th, "<tbody>", paste(rows, collapse = ""), "</tbody></table>")
}

#' Build the regulatory study report as self-contained HTML
#'
#' @param ctx A run context (see [.reg_model()]).
#' @return A single character string: a complete HTML document.
#' @seealso [build_regulatory_report_docx()]
#' @export
build_regulatory_report_html <- function(ctx) {
  m <- .reg_model(ctx)
  h1 <- 0L; h2 <- 0L; h3 <- 0L
  toc <- character(0)
  parts <- vapply(m$blocks, function(b) {
    switch(b$type,
      h1 = { h1 <<- h1 + 1L; h2 <<- 0L; h3 <<- 0L
             id <- sprintf("s%d", h1)
             toc <<- c(toc, sprintf("<li><a href='#%s'>%d. %s</a></li>", id, h1, .reg_xml_escape(b$text)))
             sprintf("<h2 id='%s'>%d. %s</h2>", id, h1, .reg_xml_escape(b$text)) },
      h2 = { h2 <<- h2 + 1L; h3 <<- 0L
             sprintf("<h3>%d.%d %s</h3>", h1, h2, .reg_xml_escape(b$text)) },
      h3 = { h3 <<- h3 + 1L
             sprintf("<h4>%d.%d.%d %s</h4>", h1, h2, h3, .reg_xml_escape(b$text)) },
      p = sprintf("<p>%s</p>", .reg_xml_escape(b$text)),
      note = sprintf("<p class='note'>%s</p>", .reg_xml_escape(b$text)),
      bullets = paste0("<ul>", paste0("<li>", .reg_xml_escape(b$items), "</li>", collapse = ""), "</ul>"),
      table = .reg_html_table(b$df, b$caption, b$header),
      figure = {
        uri <- .reg_b64img(.reg_fig_png(b$plot, b$w, b$h))
        if (is.null(uri)) "" else sprintf("<figure><img class='fig' src='%s'/><figcaption>%s</figcaption></figure>", uri, .reg_xml_escape(b$caption))
      },
      svg = sprintf("<figure class='struct'>%s<figcaption>%s</figcaption></figure>", b$svg, .reg_xml_escape(b$caption)),
      pagebreak = "",
      "")
  }, "")

  css <- "
    body{font-family:'Segoe UI',system-ui,Helvetica,Arial,sans-serif;color:#1d2733;max-width:920px;margin:0 auto;padding:28px;line-height:1.5;font-size:14px}
    .cover{border-bottom:2px solid #1f3864;padding-bottom:18px;margin-bottom:24px}
    .cover h1{font-size:26px;color:#1f3864;margin:6px 0 4px}
    .cover .sub{font-size:15px;color:#44546a;margin-bottom:14px}
    h2{font-size:18px;color:#1f3864;border-bottom:1px solid #c9d3df;padding-bottom:4px;margin-top:34px}
    h3{font-size:15px;color:#1f3864;margin-top:22px} h4{font-size:14px}
    .conf{color:#a3231b;font-weight:700;letter-spacing:.08em;font-size:12px}
    table{border-collapse:collapse;width:100%;margin:6px 0 14px;font-size:12.5px}
    th,td{border:1px solid #d5dbe2;padding:5px 8px;text-align:left;vertical-align:top}
    th{background:#e8edf3;font-weight:600}
    table.kv td:first-child{width:30%;background:#f5f7f9;font-weight:600}
    .tbl-cap{font-size:12.5px;color:#1f3864;margin-top:14px;font-weight:600}
    figure{margin:12px 0} figcaption{font-size:12.5px;color:#1f3864;font-weight:600;margin-top:4px}
    .fig{max-width:100%;border:1px solid #eef1f4}
    .struct svg{max-width:100%;height:auto}
    .note{color:#6c757d;font-style:italic;font-size:12.5px}
    .toc ol{list-style:none;padding-left:0} .toc a{color:#1f3864;text-decoration:none}
    .sig td{height:34px}
    .muted{color:#6c757d;font-size:12px}"

  cover <- paste0(
    "<div class='cover'><div class='conf'>CONFIDENTIAL</div>",
    sprintf("<div class='muted'>%s REPORT</div>", .reg_xml_escape(toupper(m$status))),
    sprintf("<h1>%s</h1>", .reg_xml_escape(m$title)),
    sprintf("<div class='sub'>Metabolite identification and profiling of %s</div>", .reg_xml_escape(m$test_article)),
    .reg_html_table(data.frame(k = names(m$cover), v = unname(m$cover), stringsAsFactors = FALSE), header = FALSE, cls = "kv"),
    "</div>")
  approval <- paste0("<h2>Approval</h2>",
    "<p>The undersigned confirm that this report accurately reflects the raw data and the methods used.</p>",
    gsub("<table class=''>", "<table class='sig'>", .reg_html_table(m$signatures), fixed = TRUE),
    "<h2>Compliance and Quality Assurance Statements</h2>",
    "<p>[GLP compliance statement, or statement of non-GLP status, per the test facility's procedures.]</p>",
    "<p>[Quality Assurance statement listing the phases inspected and the dates reported to the Study Director and management.]</p>")
  toc_html <- paste0("<div class='toc'><h2>Contents</h2><ol>", paste(toc, collapse = ""), "</ol></div>")
  abbr <- paste0("<h2>Abbreviations</h2>", .reg_html_table(m$abbreviations))

  paste0("<!doctype html><html><head><meta charset='utf-8'><title>",
         .reg_xml_escape(m$title), "</title><style>", css, "</style></head><body>",
         cover, approval, toc_html, abbr, paste(parts, collapse = "\n"),
         "<p class='muted'>Generated by OligoMet Profiler. Scientific content is computed from the run; administrative, GLP and interpretive fields in brackets are completed by the study team.</p>",
         "</body></html>")
}

## ---- Word (.docx) renderer -------------------------------------------------
.reg_col_widths <- function(df, total = .REG_BODY_W) {
  chars <- vapply(seq_along(df), function(j) {
    head_w <- max(nchar(strsplit(names(df)[j], " ", fixed = TRUE)[[1]]))
    body_w <- if (nrow(df)) stats::quantile(nchar(as.character(df[[j]])), 0.9, names = FALSE) else 0
    max(head_w, min(body_w, 42), 4)
  }, numeric(1))
  w <- sqrt(chars); w / sum(w) * total
}

.reg_ft <- function(df, widths = NULL, header = TRUE, fontsize = 8.5) {
  df[] <- lapply(df, as.character)
  ft <- flextable::flextable(df)
  if (!header) ft <- flextable::delete_part(ft, part = "header")
  ft <- flextable::font(ft, fontname = "Calibri", part = "all")
  ft <- flextable::fontsize(ft, size = fontsize, part = "all")
  bd <- officer::fp_border(color = "#9aa8b8", width = 0.5)
  ft <- flextable::border_remove(ft)
  ft <- flextable::border_outer(ft, border = officer::fp_border(color = "#1f3864", width = 1))
  ft <- flextable::border_inner_h(ft, border = bd, part = "all")
  ft <- flextable::border_inner_v(ft, border = bd, part = "all")
  if (header) {
    ft <- flextable::bg(ft, bg = "#dfe6ee", part = "header")
    ft <- flextable::bold(ft, part = "header")
  } else {
    ft <- flextable::bg(ft, j = 1, bg = "#f2f5f8", part = "body")
    ft <- flextable::bold(ft, j = 1, part = "body")
  }
  ft <- flextable::valign(ft, valign = "top", part = "all")
  ft <- flextable::padding(ft, padding.top = 2, padding.bottom = 2, padding.left = 4, padding.right = 4, part = "all")
  if (is.null(widths)) widths <- .reg_col_widths(df)
  widths <- widths / sum(widths) * min(sum(widths), .REG_BODY_W)
  ft <- flextable::width(ft, j = seq_along(df), width = widths)
  if (nrow(df) <= 30) ft <- tryCatch(flextable::paginate(ft, init = TRUE, hdr_ftr = TRUE), error = function(e) ft)
  flextable::set_table_properties(ft, layout = "fixed")
}

# Ask Word to refresh fields (table of contents, page counts) on open.
.reg_docx_update_fields <- function(file) {
  tryCatch({
    td <- tempfile("regdocx"); dir.create(td)
    utils::unzip(file, exdir = td)
    sp <- file.path(td, "word", "settings.xml")
    s <- paste(readLines(sp, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    if (!grepl("updateFields", s, fixed = TRUE)) {
      # schema order: updateFields precedes these elements, when present
      after <- c("hdrShapeDefaults", "footnotePr", "endnotePr", "compat", "docVars", "rsids",
                 "mathPr", "attachedSchema", "themeFontLang", "clrSchemeMapping",
                 "doNotIncludeSubdocsInStats", "doNotAutoCompressPictures", "forceUpgrade",
                 "captions", "readModeInkLockDown", "smartTagType", "schemaLibrary",
                 "shapeDefaults", "doNotEmbedSmartTags", "decimalSymbol", "listSeparator")
      pos <- regexpr(sprintf("<(w|m|sl|o):(%s)[ />]", paste(after, collapse = "|")), s)
      if (pos < 0) pos <- regexpr("</w:settings>", s, fixed = TRUE)
      s <- paste0(substr(s, 1, pos - 1), "<w:updateFields w:val=\"true\"/>", substr(s, pos, nchar(s)))
      writeLines(s, sp, useBytes = TRUE)
      files <- list.files(td, recursive = TRUE, all.files = TRUE, no.. = TRUE)
      out <- tempfile(fileext = ".docx")
      zip::zip(out, files = files, root = td, mode = "mirror")
      file.copy(out, file, overwrite = TRUE)
    }
    TRUE
  }, error = function(e) FALSE)
}

#' Build the regulatory study report as a Word (.docx) file
#'
#' Native Word output with a cover page, approval and compliance page,
#' table of contents, abbreviations, numbered sections, numbered table and
#' figure captions, a running header and a "Page X of Y" footer, flextable
#' tables sized to the page, embedded figures and structure maps. Requires
#' the `officer`, `flextable` and `rsvg` packages.
#'
#' @param ctx A run context (see [.reg_model()]).
#' @param file Output path for the `.docx`.
#' @return `file`, invisibly.
#' @seealso [build_regulatory_report_html()]
#' @export
build_regulatory_report_docx <- function(ctx, file) {
  for (p in c("officer", "flextable", "rsvg")) {
    if (!requireNamespace(p, quietly = TRUE))
      stop("the .docx report needs the '", p, "' package")
  }
  m <- .reg_model(ctx)
  d <- m$details
  navy <- "#1f3864"
  tx <- function(size = 10, bold = FALSE, color = "#1d2733", italic = FALSE)
    officer::fp_text(font.size = size, bold = bold, color = color, italic = italic, font.family = "Calibri")
  par <- function(..., align = "left", before = 0, after = 6, keep = FALSE)
    officer::fp_par(text.align = align, padding.top = before, padding.bottom = after, keep_with_next = keep)

  doc <- officer::read_docx()
  doc <- tryCatch(officer::set_doc_properties(doc, title = m$title, subject = "Metabolite identification report",
                                              creator = "OligoMet Profiler"), error = function(e) doc)

  add_fpar <- function(doc, text, fpt = tx(), fpp = par()) officer::body_add_fpar(doc, officer::fpar(officer::ftext(text, fpt), fp_p = fpp))
  add_cap <- function(doc, text) add_fpar(doc, text, tx(9, TRUE, navy), par(before = 8, after = 3, keep = TRUE))
  add_tbl <- function(doc, df, cap = NULL, widths = NULL, header = TRUE) {
    if (is.null(df) || nrow(df) == 0) return(doc)
    if (!is.null(cap)) doc <- add_cap(doc, cap)
    doc <- flextable::body_add_flextable(doc, .reg_ft(df, widths, header), align = "left")
    add_fpar(doc, "", tx(4), par(after = 2))
  }

  # ---- cover page ----
  doc <- add_fpar(doc, "CONFIDENTIAL", tx(10, TRUE, "#a3231b"), par(after = 60))
  doc <- add_fpar(doc, sprintf("%s REPORT", toupper(m$status)), tx(11, TRUE, "#44546a"), par(after = 6))
  doc <- add_fpar(doc, m$title, tx(24, TRUE, navy), par(after = 6))
  doc <- add_fpar(doc, sprintf("Metabolite identification and profiling of %s", m$test_article), tx(14, FALSE, "#44546a"), par(after = 36))
  cover <- data.frame(Field = names(m$cover), Value = unname(m$cover), stringsAsFactors = FALSE)
  doc <- flextable::body_add_flextable(doc, .reg_ft(cover, c(2.0, 4.27), header = FALSE, fontsize = 10), align = "left")
  doc <- add_fpar(doc, "For research use only. Generated by OligoMet Profiler; administrative, GLP and interpretive fields in brackets are completed by the study team.",
                  tx(8, FALSE, "#6c757d", TRUE), par(before = 36))
  doc <- officer::body_add_break(doc)

  # ---- approval, compliance, QA ----
  doc <- add_fpar(doc, "Approval", tx(14, TRUE, navy), par(after = 8))
  doc <- add_fpar(doc, "The undersigned confirm that this report accurately reflects the raw data and the methods used in the study.", tx(10), par(after = 10))
  sig <- .reg_ft(m$signatures, c(1.9, 1.7, 1.67, 1.0), fontsize = 9.5)
  sig <- flextable::height(sig, height = 0.5, part = "body")
  sig <- flextable::hrule(sig, rule = "atleast", part = "body")
  doc <- flextable::body_add_flextable(doc, sig, align = "left")
  doc <- add_fpar(doc, "Compliance Statement", tx(12, TRUE, navy), par(before = 18, after = 4))
  doc <- add_fpar(doc, "[GLP compliance statement, or statement of non-GLP status, per the test facility's procedures.]", tx(10))
  doc <- add_fpar(doc, "Quality Assurance Statement", tx(12, TRUE, navy), par(before = 12, after = 4))
  doc <- add_fpar(doc, "[Phases inspected, dates of inspection and dates reported to the Study Director and management.]", tx(10))
  doc <- officer::body_add_break(doc)

  # ---- contents and abbreviations ----
  doc <- add_fpar(doc, "Contents", tx(14, TRUE, navy), par(after = 8))
  doc <- officer::body_add_toc(doc, level = 2)
  doc <- add_fpar(doc, "Abbreviations", tx(14, TRUE, navy), par(before = 18, after = 8))
  doc <- flextable::body_add_flextable(doc, .reg_ft(m$abbreviations, c(1.4, 4.87), fontsize = 9), align = "left")
  doc <- officer::body_add_break(doc)

  # ---- body ----
  first_h1 <- TRUE
  for (b in m$blocks) {
    doc <- switch(b$type,
      h1 = {
        if (!first_h1 && b$text %in% c("Results", "Appendices")) doc <- officer::body_add_break(doc)
        first_h1 <- FALSE
        officer::body_add_par(doc, b$text, style = "heading 1")
      },
      h2 = officer::body_add_par(doc, b$text, style = "heading 2"),
      h3 = officer::body_add_par(doc, b$text, style = "heading 3"),
      p = add_fpar(doc, b$text, tx(10), par(after = 6)),
      note = add_fpar(doc, b$text, tx(9, FALSE, "#6c757d", TRUE), par(after = 6)),
      bullets = {
        for (it in b$items) doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(paste0(.REG_BULLET, "  "), tx(10)), officer::ftext(it, tx(10)),
                                                                         fp_p = officer::fp_par(padding.left = 14, padding.bottom = 3)))
        doc
      },
      table = add_tbl(doc, b$df, b$caption, b$widths, b$header),
      figure = {
        png <- .reg_fig_png(b$plot, b$w, b$h)
        if (is.null(png)) doc else {
          doc <- add_cap(doc, b$caption)
          officer::body_add_img(doc, png, width = b$w, height = b$h)
        }
      },
      svg = {
        png <- .reg_svg_png(b$svg)
        if (is.null(png)) doc else {
          dims <- .reg_png_dims(png)
          w <- .REG_BODY_W
          h <- if (!is.null(dims)) w * dims[2] / dims[1] else 1.2
          doc <- add_cap(doc, b$caption)
          officer::body_add_img(doc, png, width = w, height = h)
        }
      },
      pagebreak = officer::body_add_break(doc),
      doc)
  }

  # ---- page setup: running header and "Page X of Y" footer ----
  hf <- tx(8, FALSE, "#44546a")
  header <- officer::block_list(officer::fpar(
    officer::ftext(sprintf("Study %s  |  %s", .reg_get(d, "study_number"), m$title), hf),
    fp_p = officer::fp_par(text.align = "right")))
  footer <- officer::block_list(officer::fpar(
    officer::ftext("CONFIDENTIAL    |    Page ", hf), officer::run_word_field("PAGE", prop = hf),
    officer::ftext(" of ", hf), officer::run_word_field("NUMPAGES", prop = hf),
    fp_p = officer::fp_par(text.align = "center")))
  doc <- officer::body_set_default_section(doc, officer::prop_section(
    page_size = officer::page_size(orient = "portrait", width = 8.27, height = 11.69),
    page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1, header = 0.45, footer = 0.45, gutter = 0),
    type = "continuous", header_default = header, footer_default = footer))

  print(doc, target = file)
  .reg_docx_update_fields(file)
  invisible(file)
}
