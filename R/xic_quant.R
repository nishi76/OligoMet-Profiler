# =============================================================================
# xic_quant.R
# Targeted summed-XIC quantitation, the same "summation of ions" approach
# Chromeleon CDS uses when several ions are assigned to one quantitation
# component: each ion's extracted-ion chromatogram (XIC) is taken within a
# ppm tolerance, the XICs are summed point by point into one chromatogram,
# and that single chromatogram is integrated once over one window with a
# baseline.
#
# Workflow:
#   1. build_xic_targets()      -- for every metabolite identified by the MS1
#      matching step, the candidate ions: every charge state in z_range x the
#      n most abundant isotope clusters (theoretical, from the formula).
#   2. run_xic_quantitation()   -- one streaming pass per raw file (Python,
#      inst/python/oligomet_deconv/xic.py): per-ion areas over a window shared
#      by all ions of the target, each with its own linear baseline.
#   3. summarize_xic_quant()    -- picks WHICH ions to sum (top-N charge
#      states x top-N isotopes) and adds their areas. With linear baselines
#      the sum of per-ion areas equals the area of the summed chromatogram,
#      so changing N here is exact and needs no re-extraction.
#
# Ion selection is fixed per metabolite across all samples, the way a
# Chromeleon quantitation method fixes its ion list:
#   * isotopes: ranked by THEORETICAL abundance (stable, instrument-free);
#   * charge states: ranked by observed signal in the reference samples
#     (calibration standards + QCs when present, else every sample), since
#     the charge-state distribution depends on source, mobile phase and
#     ion-pairing conditions, not on the formula. A fixed list (e.g. "3-8")
#     overrides the ranking.
# =============================================================================

.C13_SPACING <- 1.0033548

#' Theoretical isotope clusters of a formula
#'
#' Collapses the isotopic fine structure (e.g. 13C2 vs 34S at M+2) into the
#' nominal clusters an Orbitrap resolves for a multiply charged oligo, with
#' an abundance-weighted cluster mass.
#'
#' @param formula Formula vector or string.
#' @param mono_mass Monoisotopic mass (cluster 0 anchor).
#' @param n Number of clusters to return, most abundant first.
#' @param use_envipat Use enviPat when installed (else the built-in engine).
#' @return data.frame: `cluster` (0 = monoisotopic, 1 = M+1, ...), `mass`,
#'   `abundance` (fraction of the whole pattern), `iso_rank` (1 = most
#'   abundant).
#' @export
isotope_clusters <- function(formula, mono_mass, n = 8, use_envipat = TRUE) {
  pat <- isotope_pattern(formula, threshold = 1e-6, n_top = 400, use_envipat = use_envipat)
  if (is.null(pat) || nrow(pat) == 0) return(data.frame())
  k <- round((pat$mass - mono_mass) / .C13_SPACING)
  ab <- tapply(pat$abundance, k, sum)
  wm <- tapply(pat$mass * pat$abundance, k, sum) / ab
  out <- data.frame(cluster = as.integer(names(ab)), mass = as.numeric(wm),
                    abundance = as.numeric(ab) / sum(ab), stringsAsFactors = FALSE)
  out <- out[order(-out$abundance), ]
  out$iso_rank <- seq_len(nrow(out))
  out <- utils::head(out, n)
  rownames(out) <- NULL
  out
}

#' Parse a charge-state list such as "3-8" or "4,5,6,7"
#'
#' @param x Character string; blank returns `NULL` (use automatic ranking).
#' @return Sorted unique integer vector, or `NULL`.
#' @export
parse_charge_spec <- function(x) {
  if (is.null(x) || length(x) == 0 || !nzchar(trimws(x))) return(NULL)
  parts <- trimws(strsplit(x, "[,; ]+")[[1]])
  parts <- parts[nzchar(parts)]
  z <- unlist(lapply(parts, function(p) {
    if (grepl("^[0-9]+\\s*-\\s*[0-9]+$", p)) {
      r <- as.integer(strsplit(p, "\\s*-\\s*")[[1]])
      seq(min(r), max(r))
    } else suppressWarnings(as.integer(p))
  }))
  z <- sort(unique(z[!is.na(z) & z > 0]))
  if (length(z) == 0) NULL else z
}

#' Build candidate XIC ions for every identified metabolite
#'
#' One target per metabolite id found in `ms1_matches`. The (PS-oxidation,
#' adduct) variant used is the one with the most matched signal across
#' samples; the expected retention time is the median matched RT of that
#' variant.
#'
#' @param mets Metabolite library ([generate_metabolites()]).
#' @param ms1_matches Identification results
#'   ([annotate_metabolites_batch()]`$ms1_matches`).
#' @param dict Chemistry dictionary.
#' @param z_range Candidate charge states.
#' @param n_iso_candidates Isotope clusters extracted per charge state (the
#'   final quantitation may use fewer; see [summarize_xic_quant()]).
#' @param h_offset Envelope offset, as elsewhere in the package.
#' @param use_envipat Passed to [isotope_clusters()].
#' @param mz_range Only ions inside this m/z range (the scan range) are kept.
#' @return data.frame, one row per ion: `target_id`, `ion_id`, `met_id`,
#'   `met_name`, `kind`, `k_oxid`, `adduct`, `rt_expected`, `z`, `cluster`,
#'   `iso_rank`, `rel_abundance`, `mz`.
#' @export
build_xic_targets <- function(mets, ms1_matches, dict = STANDARD_DICT, z_range = 3:12,
                              n_iso_candidates = 8, h_offset = 0, use_envipat = TRUE,
                              mz_range = c(200, 4000)) {
  if (is.null(ms1_matches) || nrow(ms1_matches) == 0) return(data.frame())
  sig <- if ("area" %in% names(ms1_matches) && any(!is.na(ms1_matches$area))) "area" else "intensity"
  lib_ids <- vapply(mets, function(m) m$id, character(1))
  rows <- lapply(unique(ms1_matches$met_id), function(mid) {
    met <- mets[[match(mid, lib_ids)]]
    if (is.null(met)) return(NULL)
    mm <- ms1_matches[ms1_matches$met_id == mid, , drop = FALSE]
    k_ox <- if ("k_oxid" %in% names(mm)) mm$k_oxid else 0L
    ad <- if ("adduct" %in% names(mm)) mm$adduct else "H"
    key <- paste(k_ox, ad, sep = "|")
    tot <- tapply(mm[[sig]], key, sum, na.rm = TRUE)
    best <- names(tot)[which.max(tot)]
    k <- as.integer(strsplit(best, "|", fixed = TRUE)[[1]][1])
    adduct <- strsplit(best, "|", fixed = TRUE)[[1]][2]
    rt_exp <- stats::median(mm$rt[key == best], na.rm = TRUE)

    info <- metabolite_mass_info(met, dict)
    fv <- ps_oxid_formula(info$formula_vec, k)
    mono <- ps_oxid_mass(info$mono_mass, k)
    cl <- isotope_clusters(fv, mono, n = n_iso_candidates, use_envipat = use_envipat)
    if (nrow(cl) == 0) return(NULL)
    shift <- if (identical(adduct, "H")) 0 else adduct_shift(adduct)
    grid <- expand.grid(i = seq_len(nrow(cl)), z = z_range)
    mz <- (cl$mass[grid$i] + shift + h_offset - grid$z * .PROTON) / grid$z
    out <- data.frame(
      target_id = mid, ion_id = sprintf("z%d_M+%d", grid$z, cl$cluster[grid$i]),
      met_id = mid, met_name = met$name, kind = met$kind, k_oxid = k, adduct = adduct,
      rt_expected = rt_exp, z = grid$z, cluster = cl$cluster[grid$i],
      iso_rank = cl$iso_rank[grid$i], rel_abundance = cl$abundance[grid$i],
      mz = mz, stringsAsFactors = FALSE)
    out[out$mz >= mz_range[1] & out$mz <= mz_range[2], , drop = FALSE]
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

## ---- Python invocation -------------------------------------------------------
.run_python_module <- function(args, python_bin, module_dir, what) {
  if (is.na(python_bin)) {
    stop("No python3/python interpreter found on PATH. ", what, " requires Python 3.9+ ",
         "with the packages in inst/python/requirements.txt installed.")
  }
  if (!dir.exists(module_dir)) {
    stop("Could not locate the oligomet_deconv Python module (looked in: ", module_dir, ").")
  }
  old <- Sys.getenv("PYTHONPATH", unset = NA)
  new <- if (is.na(old) || !nzchar(old)) module_dir else paste(module_dir, old, sep = .Platform$path.sep)
  on.exit(if (is.na(old)) Sys.unsetenv("PYTHONPATH") else Sys.setenv(PYTHONPATH = old), add = TRUE)
  Sys.setenv(PYTHONPATH = new)
  status <- system2(python_bin, args = as.character(args), stdout = TRUE, stderr = TRUE)
  code <- attr(status, "status")
  if (!is.null(code) && code != 0) stop(what, " failed:\n", paste(status, collapse = "\n"))
  status
}

#' Run targeted summed-XIC extraction over raw files
#'
#' @param files mzML/mzXML paths.
#' @param targets [build_xic_targets()] output.
#' @param output_dir Where the Python step writes its tables.
#' @param ppm XIC extraction half-width (ppm).
#' @param rt_window Apex search window around each target's expected RT (min).
#' @param max_half_width Maximum peak half-width on each side of the apex (min).
#' @param edge_frac Peak edges where the smoothed summed trace falls within
#'   this fraction of the apex height above local background.
#' @param n_workers Parallel worker processes (one file each).
#' @param python_bin,module_dir Python interpreter and `inst/python` path.
#' @param progress Optional `function(msg)` progress callback.
#' @return A list: `ions` (one row per sample x ion, joined to the target
#'   metadata: `area`, `height`, shared `rt_apex`/`rt_start`/`rt_end`,
#'   `sn_all_ions`), `traces` (summed candidate-ion chromatogram per sample x
#'   target around the search window), `params`, and `log`.
#' @export
run_xic_quantitation <- function(files, targets, output_dir = tempfile("xic_"),
                                 ppm = 10, rt_window = 0.5, max_half_width = 0.6,
                                 edge_frac = 0.01, n_workers = NULL,
                                 python_bin = find_python(),
                                 module_dir = .find_deconv_module_dir(),
                                 progress = NULL) {
  if (is.null(targets) || nrow(targets) == 0) stop("no XIC targets to extract")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  tpath <- file.path(output_dir, "xic_targets.tsv")
  utils::write.table(targets[, c("target_id", "ion_id", "mz", "rt_expected")], tpath,
                     sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  args <- c("-m", "oligomet_deconv.xic_cli", "--input", files, "--targets", tpath,
            "--output-dir", output_dir, "--ppm", ppm, "--rt-window", rt_window,
            "--max-half-width", max_half_width, "--edge-frac", edge_frac)
  if (!is.null(n_workers)) args <- c(args, "--n-workers", n_workers)
  if (!is.null(progress)) progress("Extracting summed XICs (Python)...")
  log <- .run_python_module(args, python_bin, module_dir, "Targeted XIC quantitation")

  ions <- utils::read.delim(file.path(output_dir, "xic_ions.tsv"), stringsAsFactors = FALSE)
  traces <- utils::read.delim(file.path(output_dir, "xic_traces.tsv"), stringsAsFactors = FALSE)
  ions$sample <- as.character(ions$sample)
  traces$sample <- as.character(traces$sample)
  ions$target_id <- as.character(ions$target_id)
  traces$target_id <- as.character(traces$target_id)
  meta_cols <- setdiff(names(targets), c("mz", "rt_expected"))
  ions <- merge(ions, targets[, c(meta_cols, "mz", "rt_expected")], by = c("target_id", "ion_id"),
                sort = FALSE)
  list(ions = ions, traces = traces,
       params = list(ppm = ppm, rt_window = rt_window, max_half_width = max_half_width,
                     edge_frac = edge_frac),
       log = log)
}

## ---- Ion selection + summation -------------------------------------------------

# Reference samples for ranking charge states: calibration standards and QCs
# when the design has them (clean, high-signal, the samples a curve is built
# from), otherwise every sample.
.xic_reference_samples <- function(samples, sample_meta) {
  if (!is.null(sample_meta) && all(c("sample", "sample_type") %in% names(sample_meta))) {
    ref <- sample_meta$sample[!is.na(sample_meta$sample_type) &
                                sample_meta$sample_type %in% c("standard", "quality_control")]
    ref <- intersect(ref, samples)
    if (length(ref) > 0) return(ref)
  }
  samples
}

#' Choose ions and sum them into one quantitative value per sample
#'
#' @param ions `run_xic_quantitation()$ions`.
#' @param sample_meta Optional sample table (`sample`, `sample_type`), used to
#'   rank charge states on standards/QCs.
#' @param n_charges Number of charge states to sum (ranked by signal in the
#'   reference samples). Ignored when `charges` is given.
#' @param n_isotopes Number of isotope clusters to sum (ranked by
#'   theoretical abundance).
#' @param charges Optional fixed charge-state vector (see
#'   [parse_charge_spec()]), applied to every metabolite.
#' @return A list:
#'   * `quant`: one row per (metabolite, sample) in the `ms1_matches` shape --
#'     `met_id`, `met_name`, `kind`, `sample`, `area` (summed XIC area, NA
#'     when no peak), `intensity` (summed apex height), `rt`, `rt_start`,
#'     `rt_end`, `n_ions`, `n_ions_detected`, `sn`, `charges_used`,
#'     `isotopes_used`, `k_oxid`, `adduct`;
#'   * `selection`: one row per metabolite -- the ions used and
#'     `pct_candidate_signal`, the share of all extracted candidate signal
#'     (reference samples) the chosen ions capture;
#'   * `contributions`: per metabolite x charge x isotope, % of the
#'     reference-sample candidate signal, with a `selected` flag.
#' @export
summarize_xic_quant <- function(ions, sample_meta = NULL, n_charges = 5, n_isotopes = 5,
                                charges = NULL) {
  empty <- list(quant = data.frame(), selection = data.frame(), contributions = data.frame())
  if (is.null(ions) || nrow(ions) == 0) return(empty)
  ions$area_pos <- pmax(ions$area, 0)
  ions$area_pos[is.na(ions$rt_apex)] <- 0
  ref <- .xic_reference_samples(unique(ions$sample), sample_meta)

  quant <- list(); sel <- list(); contrib <- list()
  for (tid in unique(ions$target_id)) {
    d <- ions[ions$target_id == tid, , drop = FALSE]
    iso_keep <- sort(unique(d$iso_rank))[seq_len(min(n_isotopes, length(unique(d$iso_rank))))]
    dref <- d[d$sample %in% ref, , drop = FALSE]
    z_sig <- tapply(dref$area_pos[dref$iso_rank %in% iso_keep], dref$z[dref$iso_rank %in% iso_keep], sum)
    z_sig <- z_sig[order(-z_sig)]
    z_keep <- if (!is.null(charges)) intersect(charges, unique(d$z))
              else as.integer(names(z_sig))[seq_len(min(n_charges, length(z_sig)))]
    chosen <- d$z %in% z_keep & d$iso_rank %in% iso_keep

    tot_ref <- sum(dref$area_pos)
    cz <- stats::aggregate(area_pos ~ z + cluster + iso_rank, data = dref, FUN = sum)
    cz$pct <- if (tot_ref > 0) 100 * cz$area_pos / tot_ref else NA_real_
    cz$selected <- cz$z %in% z_keep & cz$iso_rank %in% iso_keep
    contrib[[tid]] <- data.frame(met_id = d$met_id[1], z = cz$z, cluster = cz$cluster,
                                 iso_rank = cz$iso_rank, pct_of_candidate_signal = cz$pct,
                                 selected = cz$selected, stringsAsFactors = FALSE)
    clusters_used <- sort(unique(d$cluster[d$iso_rank %in% iso_keep]))
    charges_txt <- paste(sort(z_keep), collapse = ",")
    iso_txt <- paste0("M+", clusters_used, collapse = ",")
    sel[[tid]] <- data.frame(
      met_id = d$met_id[1], met_name = d$met_name[1], kind = d$kind[1],
      k_oxid = d$k_oxid[1], adduct = d$adduct[1], rt_expected = d$rt_expected[1],
      charges_used = charges_txt, isotopes_used = iso_txt,
      n_ions = sum(d$z %in% z_keep & d$iso_rank %in% iso_keep & d$sample == d$sample[1]),
      pct_candidate_signal = if (tot_ref > 0) 100 * sum(dref$area_pos[dref$z %in% z_keep & dref$iso_rank %in% iso_keep]) / tot_ref else NA_real_,
      n_reference_samples = length(unique(dref$sample)),
      stringsAsFactors = FALSE)

    dc <- d[chosen, , drop = FALSE]
    for (s in unique(d$sample)) {
      x <- dc[dc$sample == s, , drop = FALSE]
      w <- d[d$sample == s, , drop = FALSE][1, ]
      found <- nrow(x) > 0 && !is.na(x$rt_apex[1])
      a <- if (found) sum(x$area) else NA_real_
      if (!is.na(a) && a <= 0) a <- NA_real_
      quant[[length(quant) + 1]] <- data.frame(
        met_id = w$met_id, met_name = w$met_name, kind = w$kind, sample = s,
        area = a, intensity = if (is.na(a)) NA_real_ else sum(pmax(x$height, 0)),
        rt = w$rt_apex, rt_start = w$rt_start, rt_end = w$rt_end,
        n_ions = nrow(x), n_ions_detected = sum(x$area > 0, na.rm = TRUE),
        sn = w$sn_all_ions, charges_used = charges_txt, isotopes_used = iso_txt,
        k_oxid = w$k_oxid, adduct = w$adduct, stringsAsFactors = FALSE)
    }
  }
  q <- do.call(rbind, quant)
  rownames(q) <- NULL
  list(quant = q[!is.na(q$area), , drop = FALSE],
       selection = do.call(rbind, unname(sel)),
       contributions = do.call(rbind, unname(contrib)))
}

## ---- Plots -------------------------------------------------------------------------

# Colour-blind-safe categorical order (Okabe-Ito, black dropped), checked
# with a CVD validator: worst adjacent deutan Delta E 9.6. Assigned in this
# fixed order, never cycled -- a 7th series is not given a recycled colour.
.OKABE_ITO <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7")

# 7-point quadratic Savitzky-Golay smoothing -- the same filter the XIC
# integrator uses to find apex and peak edges (areas use the raw trace).
.sg7 <- function(y) {
  if (length(y) < 7) return(y)
  k <- c(-2, 3, 6, 7, 6, 3, -2) / 21
  out <- stats::filter(y, k, sides = 2)
  out[is.na(out)] <- y[is.na(out)]
  as.numeric(out)
}

#' Plot summed XIC traces with their integration windows
#'
#' @param traces `run_xic_quantitation()$traces`.
#' @param quant `summarize_xic_quant()$quant` (for the integration windows).
#' @param met_id Metabolite to show.
#' @param samples Optional subset of samples.
#' @param sample_meta Optional sample table (`sample`, `sample_type`,
#'   `group`) for colouring by sample type or group.
#' @param color_by `"sample_type"` (default), `"group"`, or `"sample"`. A
#'   colour-blind-safe palette of six colours is used in a fixed order;
#'   colouring by sample shows at most six samples, so no colour repeats.
#' @param smoothed Draw the 7-point Savitzky-Golay trace the integrator uses
#'   for peak finding instead of the raw summed trace.
#' @return A ggplot object.
#' @export
plot_xic_traces <- function(traces, quant, met_id, samples = NULL, sample_meta = NULL,
                            color_by = c("sample_type", "group", "sample"), smoothed = FALSE) {
  color_by <- match.arg(color_by)
  tr <- traces[traces$target_id == met_id, , drop = FALSE]
  if (!is.null(samples) && length(samples) > 0) tr <- tr[tr$sample %in% samples, , drop = FALSE]
  if (nrow(tr) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "No XIC trace for this metabolite"))
  }
  note <- NULL
  smp <- unique(tr$sample)
  if (color_by == "sample" && length(smp) > length(.OKABE_ITO)) {
    note <- sprintf("Showing the first %d of %d samples -- pick samples to compare, or colour by sample type.",
                    length(.OKABE_ITO), length(smp))
    smp <- smp[seq_along(.OKABE_ITO)]
    tr <- tr[tr$sample %in% smp, , drop = FALSE]
  }
  tr <- tr[order(tr$sample, tr$rt), , drop = FALSE]
  if (smoothed) tr$intensity <- stats::ave(tr$intensity, tr$sample, FUN = .sg7)

  lookup <- function(col) {
    if (is.null(sample_meta) || !col %in% names(sample_meta)) return(rep("", nrow(tr)))
    v <- sample_meta[[col]][match(tr$sample, sample_meta$sample)]
    ifelse(is.na(v) | !nzchar(v), "", v)
  }
  tr$.col <- switch(color_by,
    sample = tr$sample,
    sample_type = { v <- lookup("sample_type"); ifelse(nzchar(v), gsub("_", " ", v), "unassigned") },
    group = { v <- lookup("group"); ifelse(nzchar(v), v, "no group") })
  lv <- unique(tr$.col)
  if (color_by == "sample_type") {
    pref <- c("standard", "quality control", "unknown", "reagent blank", "matrix blank", "unassigned")
    lv <- c(intersect(pref, lv), setdiff(lv, pref))
  }
  if (length(lv) > length(.OKABE_ITO)) {
    keep <- lv[seq_len(length(.OKABE_ITO) - 1)]
    tr$.col[!tr$.col %in% keep] <- "other"
    lv <- c(keep, "other")
    note <- c(note, "More categories than distinct colours: the rest are shown as 'other'.")
  }
  pal <- stats::setNames(.OKABE_ITO[seq_along(lv)], lv)
  if ("other" %in% lv) pal["other"] <- "#8A8A8A"
  tr$.col <- factor(tr$.col, levels = lv)

  win <- quant[quant$met_id == met_id & quant$sample %in% unique(tr$sample), , drop = FALSE]
  p <- ggplot2::ggplot(tr, ggplot2::aes(x = .data$rt, y = .data$intensity,
                                        color = .data$.col, group = .data$sample))
  if (nrow(win) > 0 && any(!is.na(win$rt_start))) {
    # One light band for the integration window (median start/end across
    # the shown samples) rather than one overlapping band per sample.
    p <- p + ggplot2::annotate("rect", xmin = stats::median(win$rt_start, na.rm = TRUE),
                               xmax = stats::median(win$rt_end, na.rm = TRUE),
                               ymin = -Inf, ymax = Inf, fill = "#1f2430", alpha = 0.05)
  }
  p + ggplot2::geom_hline(yintercept = 0, color = "#9AA0A6", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 1.1, lineend = "round", linejoin = "round", alpha = 0.95) +
    ggplot2::scale_color_manual(values = pal, name = NULL) +
    ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                                expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::labs(x = "Retention time (min)", y = "Summed XIC intensity",
                  title = paste0("Summed XIC -- ", met_id),
                  subtitle = paste(c(
                    paste0(if (smoothed) "Savitzky-Golay smoothed (7 pt, as used for peak finding)" else "Raw summed trace",
                           "; shaded = integration window (median across shown samples)"),
                    note), collapse = "\n")) +
    ggplot2::theme_minimal(base_size = 12, base_family = "Liberation Sans") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#EEEEEE", linewidth = 0.3),
      axis.line.x = ggplot2::element_line(color = "#5F6368", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(color = "#5F6368", linewidth = 0.3),
      axis.text = ggplot2::element_text(color = "#3C4043"),
      plot.subtitle = ggplot2::element_text(color = "#5F6368", size = 9.5),
      legend.position = "right", legend.key.width = ggplot2::unit(18, "pt"))
}

#' Plot the charge-state x isotope contribution map
#'
#' @param contributions `summarize_xic_quant()$contributions`.
#' @param met_id Metabolite to show.
#' @return A ggplot tile map: % of candidate signal per charge state and
#'   isotope cluster in the reference samples, selected ions outlined.
#' @export
plot_xic_contributions <- function(contributions, met_id) {
  d <- contributions[contributions$met_id == met_id, , drop = FALSE]
  if (nrow(d) == 0 || all(is.na(d$pct_of_candidate_signal))) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "No XIC signal for this metabolite"))
  }
  d$iso_lab <- factor(paste0("M+", d$cluster), levels = paste0("M+", sort(unique(d$cluster))))
  d$z_lab <- factor(paste0("z=", d$z), levels = paste0("z=", sort(unique(d$z))))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$iso_lab, y = .data$z_lab)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$pct_of_candidate_signal),
                       color = "white", linewidth = 0.6) +
    ggplot2::geom_tile(data = d[d$selected, , drop = FALSE], fill = NA, color = "#1f2430",
                       linewidth = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", .data$pct_of_candidate_signal)),
                       size = 3, color = "#1f2430") +
    ggplot2::scale_fill_gradient(low = "#F4F1EA", high = "#C8912E", name = "% signal") +
    ggplot2::labs(x = "Isotope cluster", y = "Charge state",
                  title = paste0("Where the signal is -- ", met_id),
                  subtitle = "% of all extracted candidate signal (reference samples); outlined = summed") +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans") +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
