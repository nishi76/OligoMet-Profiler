# =============================================================================
# studio.R
# "Studio" tab: a three-pane workbench (constructs/samples | sequence map,
# species table, MS/MS mirror | species card, isotope pattern, composition).
#
# Data flow
#   rv$spec / rv$mets / rv$batch_ms_results  -->  studio_payload()  -->  JSON
#   click on a species or sample (JS)        -->  studio_detail()   -->  JSON
# The browser side (inst/app/www/studio.js) only draws what it is sent; every
# mass, formula, charge state, isotope abundance and fragment m/z comes from
# the same engine functions the other tabs use.
#
# Species are keyed by (met_id, k_oxid), so phosphorothioate -> phosphodiester
# variants of one library entry are separate rows here instead of being
# collapsed to the tallest peak. Signal for a species in a sample is the
# largest charge-state/adduct signal for that species (the same rule the
# Degradation Summary uses) and percentages are shares of the summed signal of
# all species shown, with no response-factor correction.
#
# Demo mode builds the same inputs from a simulated acquisition of the bundled
# Givosiran sense strand (masses/fragments are real engine output; signals,
# mass errors and spectra are simulated) so the screen can be explored before
# any batch has been run.
# =============================================================================

## ---- Helpers ---------------------------------------------------------------
.studio_with_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

.studio_hash <- function(s) {
  v <- utf8ToInt(paste(s, collapse = "|"))
  as.integer((sum(v * seq_along(v)) %% 2147483000) + 1)
}

.studio_conj_label <- function(code, dict) {
  if (is.null(code) || length(code) == 0 || is.na(code) || code == "none" || !nzchar(code)) return("")
  known <- c(GalNAc3_triantennary = "L96 · 3 GalNAc", GalNAc3t_m1 = "L96 · 2 GalNAc",
             GalNAc3t_m2 = "L96 · 1 GalNAc", GalNAc3t_m3 = "L96 scaffold")
  if (code %in% names(known)) return(unname(known[code]))
  nm <- dict[[code]]$name
  if (is.null(nm)) code else sub(" cap$", "", nm)
}

.studio_class <- function(kind, k_oxid = 0L) {
  if (isTRUE(k_oxid > 0)) return("ps")
  switch(kind,
         parent = "parent", conj_catab = "gal", exo_3p = "3p", exo_5p = "5p",
         endo_5frag = "endo", endo_5frag_p = "endo", endo_3frag = "endo", endo_3frag_p = "endo",
         "other")
}

.studio_range <- function(kind, n, N) {
  if (kind %in% c("exo_5p", "endo_3frag", "endo_3frag_p")) c(N - n + 1L, N) else c(1L, n)
}

.studio_key <- function(sample, met_id, k_oxid, z, adduct) {
  paste(sample, met_id, k_oxid, z, adduct, sep = "|")
}

.studio_num <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x[!is.finite(x)] <- NA; x
}

## One record per library entry (unmodified species).
.studio_records <- function(mets, spec, dict) {
  N <- spec$n
  recs <- list()
  for (met in mets) {
    if (identical(met$kind, "internal_standard")) next
    info <- tryCatch(metabolite_mass_info(met, dict), error = function(e) NULL)
    if (is.null(info)) next
    rg <- .studio_range(met$kind, met$n, N)
    recs[[met$id]] <- list(
      id = met$id, met_id = met$id, k_oxid = 0L,
      name = if (identical(met$name, "OLIGO") || !nzchar(met$name %||% "")) "Intact strand" else sub("^OLIGO\\s*", "", met$name),
      kind = met$kind, cls = .studio_class(met$kind), n = met$n,
      mass = round(info$mono_mass, 4), formula = format_formula(info$formula_vec),
      start = rg[1], end = rg[2],
      conj5 = .studio_conj_label(met$conj5, dict), conj3 = .studio_conj_label(met$conj3, dict))
  }
  recs
}

.studio_oxid_record <- function(base, k, mets, dict) {
  met <- mets[[which(vapply(mets, function(m) m$id == base$met_id, logical(1)))[1]]]
  info <- metabolite_mass_info(met, dict)
  r <- base
  r$id <- paste0(base$met_id, "+ox", k); r$k_oxid <- as.integer(k)
  r$name <- sprintf("%s (PS→PO ×%d)", base$name, k)
  r$cls <- "ps"
  r$mass <- round(ps_oxid_mass(info$mono_mass, k), 4)
  r$formula <- format_formula(ps_oxid_formula(info$formula_vec, k))
  r
}

.studio_sample_order <- function(samples, sample_meta) {
  tp <- rep(NA_real_, length(samples))
  if (!is.null(sample_meta) && "sample" %in% names(sample_meta) && "timepoint" %in% names(sample_meta)) {
    t <- sample_meta$timepoint[match(samples, sample_meta$sample)]
    tp <- .studio_num(gsub("[^0-9.]", "", as.character(t)))
  }
  if (all(!is.na(tp)) && length(samples) > 0) {
    o <- order(tp, samples)
  } else {
    o <- order(samples); tp <- rep(NA_real_, length(samples))
  }
  list(samples = samples[o], time = tp[o])
}

## Best (highest-signal) hit per (sample, species) -----------------------------
.studio_best_rows <- function(m, sig) {
  m$.sp <- paste(m$met_id, m$k_oxid, sep = "\r")
  m$.s <- m[[sig]]; m$.s[is.na(m$.s)] <- -Inf
  grp <- paste(m$sample, m$.sp, sep = "\t")
  idx <- tapply(seq_len(nrow(m)), grp, function(i) i[which.max(m$.s[i])])
  m[unlist(idx), , drop = FALSE]
}

## ---- Payload ---------------------------------------------------------------
#' Build the JSON-ready payload behind the Studio tab
#'
#' @param spec A parsed oligo spec (see [parse_input()]).
#' @param mets Library metabolites (see [generate_metabolites()]).
#' @param dict A chemistry dictionary.
#' @param ms1_matches Optional batch MS1 match table (see [match_ms1_batch()]).
#' @param ms2_confirmations Optional MS2 confirmation table
#'   (see [confirm_ms2_batch()]), used for the coverage column.
#' @param sample_meta Optional sample metadata (`sample`, `timepoint`).
#' @param mode `"results"`, `"library"` or `"demo"`; only changes the notice
#'   the page shows.
#' @param title Display name of the oligo.
#' @return A list; serialise with `jsonlite::toJSON(auto_unbox = TRUE)`.
#' @export
studio_payload <- function(spec, mets, dict = STANDARD_DICT, ms1_matches = NULL,
                           ms2_confirmations = NULL, sample_meta = NULL,
                           mode = "library", title = NULL) {
  N <- spec$n
  recs <- .studio_records(mets, spec, dict)
  parent <- Filter(function(m) identical(m$kind, "parent"), mets)
  pinfo <- if (length(parent) > 0) metabolite_mass_info(parent[[1]], dict) else NULL
  oligo <- list(
    name = title %||% (if (!is.null(spec$raw) && nchar(spec$raw[1]) < 40) spec$raw[1] else "Active oligo"),
    n = N, bases = spec$bases, sugars = spec$sugars,
    linkages = ifelse(is.na(spec$linkages), "", spec$linkages),
    conj5 = .studio_conj_label(spec$conj5, dict), conj3 = .studio_conj_label(spec$conj3, dict),
    formula = if (!is.null(pinfo)) format_formula(pinfo$formula_vec) else "",
    mass = if (!is.null(pinfo)) round(pinfo$mono_mass, 4) else NA_real_)

  have <- !is.null(ms1_matches) && is.data.frame(ms1_matches) && nrow(ms1_matches) > 0
  samples <- list(); species <- unname(recs); signal_col <- NA_character_

  if (have) {
    if (!"k_oxid" %in% names(ms1_matches)) ms1_matches$k_oxid <- 0L
    ms1_matches$k_oxid <- as.integer(ms1_matches$k_oxid)
    signal_col <- if ("area" %in% names(ms1_matches) && any(!is.na(ms1_matches$area))) "area" else "intensity"
    m <- ms1_matches[ms1_matches$met_id %in% names(recs) & !is.na(ms1_matches[[signal_col]]), , drop = FALSE]
    if (nrow(m) == 0) have <- FALSE
  }
  if (have) {
    ord <- .studio_sample_order(unique(m$sample), sample_meta)
    sid <- ord$samples
    samples <- lapply(seq_along(sid), function(i) list(id = sid[i], label = sid[i],
                                                       time = if (is.na(ord$time[i])) NULL else ord$time[i]))
    best <- .studio_best_rows(m, signal_col)
    cov_tbl <- NULL
    if (!is.null(ms2_confirmations) && is.data.frame(ms2_confirmations) && nrow(ms2_confirmations) > 0 &&
        "coverage" %in% names(ms2_confirmations)) {
      if (!"k_oxid" %in% names(ms2_confirmations)) ms2_confirmations$k_oxid <- 0L
      cov_tbl <- ms2_confirmations
      cov_tbl$.k <- paste(cov_tbl$sample, cov_tbl$met_id, cov_tbl$k_oxid, cov_tbl$z, cov_tbl$adduct, sep = "|")
    }
    sp_keys <- unique(best[, c("met_id", "k_oxid")])
    sp_keys <- sp_keys[order(sp_keys$met_id, sp_keys$k_oxid), , drop = FALSE]
    species <- list(); sums <- numeric(length(sid))
    for (i in seq_len(nrow(sp_keys))) {
      mid <- sp_keys$met_id[i]; k <- sp_keys$k_oxid[i]
      rec <- if (k == 0L) recs[[mid]] else .studio_oxid_record(recs[[mid]], k, mets, dict)
      pick <- function(col) vapply(sid, function(s) {
        r <- best[best$sample == s & best$met_id == mid & best$k_oxid == k, , drop = FALSE]
        if (nrow(r) == 0) NA_real_ else .studio_num(r[[col]][1])
      }, numeric(1))
      sg <- pick(signal_col)
      rec$signal <- unname(sg); rec$ppm <- unname(pick("ppm_error")); rec$rt <- unname(pick("rt"))
      rec$obs <- unname(pick("obs_mz")); rec$z0 <- unname(pick("z"))
      rec$cov <- rep(NA_real_, length(sid))
      if (!is.null(cov_tbl)) {
        rec$cov <- vapply(seq_along(sid), function(j) {
          r <- best[best$sample == sid[j] & best$met_id == mid & best$k_oxid == k, , drop = FALSE]
          if (nrow(r) == 0) return(NA_real_)
          cr <- cov_tbl[cov_tbl$.k == .studio_key(sid[j], mid, k, r$z[1], r$adduct[1]), , drop = FALSE]
          if (nrow(cr) == 0) NA_real_ else .studio_num(cr$coverage[1])
        }, numeric(1))
      }
      sums <- sums + ifelse(is.na(sg), 0, sg)
      species[[length(species) + 1]] <- rec
    }
    species <- lapply(species, function(r) {
      r$pct <- unname(ifelse(is.na(r$signal) | sums <= 0, NA_real_, 100 * r$signal / sums))
      r$signal <- NULL; r
    })
  } else {
    species <- lapply(species, function(r) { r$z0 <- NULL; r })
  }
  list(mode = mode, signal_col = signal_col, oligo = oligo, samples = samples, species = species)
}

## ---- Detail for one (species, sample) --------------------------------------
.studio_theory_only <- function(met, dict, z_range = 1:2) {
  frags <- generate_fragments(met, dict, c("a", "aB", "b", "bB", "w", "y"), z_range, 0, FALSE)
  rows <- list()
  for (f in frags) {
    w <- fragment_intensity_weight(f, met, dict)
    for (i in seq_len(nrow(f$mz_table)))
      rows[[length(rows) + 1]] <- list(round(f$mz_table$mz[i], 4), w, .fragment_label(f, f$mz_table$z[i]), 0L)
  }
  rows
}

#' Detail behind the Studio mirror plot, card and isotope pattern
#'
#' @param met A library metabolite.
#' @param dict A chemistry dictionary.
#' @param spectrum Optional acquired MS2 spectrum (`mz`, `intensity`); when
#'   `NULL` only the theoretical library is returned.
#' @param k_oxid PS -> PO count of the species (isotope pattern only).
#' @param z Precursor charge for the isotope pattern; default is the charge
#'   nearest 1000 Da per charge.
#' @return A list ready for JSON.
#' @export
studio_detail <- function(met, dict = STANDARD_DICT, spectrum = NULL, k_oxid = 0L,
                          z = NULL, tol_ppm = 25) {
  info <- metabolite_mass_info(met, dict)
  fv <- if (k_oxid > 0) ps_oxid_formula(info$formula_vec, k_oxid) else info$formula_vec
  mass <- if (k_oxid > 0) ps_oxid_mass(info$mono_mass, k_oxid) else info$mono_mass
  z0 <- if (is.null(z) || is.na(z)) max(3L, as.integer(round(mass / 1000))) else as.integer(z)
  env <- charge_envelope(mass, z_range = 3:12)
  iso <- tryCatch(isotope_mz_cluster(fv, z0, n_top = 8), error = function(e) NULL)
  out <- list(met = met$id, k_oxid = as.integer(k_oxid), z0 = z0,
              env = lapply(seq_len(nrow(env)), function(i) list(env$z[i], round(env$mz[i], 4))),
              iso = if (is.null(iso)) list() else lapply(seq_len(nrow(iso)), function(i)
                list(iso$iso[i], round(iso$mz[i], 4), round(iso$abundance[i] / max(iso$abundance), 5))),
              acquired = list(), theoretical = list(), score = NULL, has_spectrum = FALSE)

  has <- !is.null(spectrum) && is.data.frame(spectrum) && nrow(spectrum) > 0
  if (!has) {
    out$theoretical <- .studio_theory_only(met, dict)
    return(out)
  }
  conf <- confirm_metabolite(met, spectrum, dict, tol_ppm = tol_ppm, z_range = 1:2,
                             ion_types = c("a", "aB", "b", "bB", "w", "y"), include_internal = FALSE)
  matched <- conf$matched
  theo <- list(); keys <- character(0)
  for (f in conf$fragments) {
    w <- fragment_intensity_weight(f, met, dict)
    for (i in seq_len(nrow(f$mz_table))) {
      z <- f$mz_table$z[i]
      key <- paste(f$ion_type, f$direction, f$cleavage_site, f$frag_length, z)
      theo[[length(theo) + 1]] <- list(mz = round(f$mz_table$mz[i], 4), w = w, label = .fragment_label(f, z), key = key)
    }
  }
  wmax <- max(vapply(theo, function(t) t$w, numeric(1)), 0)
  mkeys <- if (nrow(matched) > 0) paste(matched$ion_type, matched$direction, matched$cleavage_site,
                                         matched$frag_length, matched$z) else character(0)
  out$theoretical <- lapply(theo, function(t) list(t$mz, if (wmax > 0) round(100 * t$w / wmax, 2) else 0,
                                                    t$label, as.integer(t$key %in% mkeys)))
  mx <- max(spectrum$intensity, na.rm = TRUE)
  acq <- lapply(seq_len(nrow(spectrum)), function(i) list(round(spectrum$mz[i], 4), round(100 * spectrum$intensity[i] / mx, 3), "", NULL))
  if (nrow(matched) > 0) for (i in seq_len(nrow(matched))) {
    idx <- which(spectrum$mz == matched$obs_mz[i])[1]
    if (is.na(idx)) next
    f <- Filter(function(t) t$key == mkeys[i], theo)
    lbl <- if (length(f) > 0) f[[1]]$label else ""
    acq[[idx]][[3]] <- if (nzchar(acq[[idx]][[3]])) paste(acq[[idx]][[3]], lbl, sep = " / ") else lbl
    acq[[idx]][[4]] <- round(matched$ppm_error[i], 2)
  }
  d <- conf$diagnostics
  if (is.data.frame(d) && nrow(d) > 0) for (i in seq_len(nrow(d))) {
    idx <- which(spectrum$mz == d$obs_mz[i])[1]
    if (is.na(idx)) next
    acq[[idx]][[3]] <- "PS"; acq[[idx]][[4]] <- round(d$ppm_error[i], 2)
  }
  out$acquired <- acq; out$has_spectrum <- TRUE
  sc <- conf$score
  out$score <- list(coverage = sc$coverage, covered = length(sc$covered_sites), total = sc$total_sites,
                    n_matches = sc$n_matches, median_ppm = sc$median_ppm, total_score = sc$total_score,
                    confident = isTRUE(sc$confident))
  out
}

## ---- Demo inputs -----------------------------------------------------------
.STUDIO_DEMO_TIMES <- c(0, 4, 8, 24, 48, 72)
.STUDIO_DEMO_AMP <- c(M02 = .30, M03 = .17, M04 = .09, M05 = .13, M06 = .12, M07 = .07, M08 = .04,
                      M09 = .035, M10 = .02, M11 = .015, M13 = .03, M17 = .02, M20 = .02, M21 = .015,
                      M22 = .03, M40 = .015, M41 = .01, M42 = .02)

.studio_demo_signal <- function(id, k, t) {
  if (id == "M01" && k == 0) return(exp(-log(2) * t / 22))
  if (id == "M01" && k == 1) return(0.06 * (0.015 + (1 - exp(-t / 14))) * exp(-t / 90))
  keys <- names(.STUDIO_DEMO_AMP); ord <- match(id, keys) - 1
  .STUDIO_DEMO_AMP[[id]] * (0.015 + (1 - exp(-t / (9 + 2.2 * ord)))) * exp(-t / 90)
}

#' Simulated acquisition of the bundled Givosiran sense strand
#'
#' Masses, formulas, charge states and fragment m/z are engine output; signal
#' kinetics, mass errors, retention times and MS/MS spectra are simulated. Used
#' by the Studio tab's demo mode and its tests.
#' @param dict A chemistry dictionary.
#' @return A list with `spec`, `mets`, `dict`, `ms1_matches`,
#'   `ms2_confirmations`, `ms2_spectra` and `sample_meta`.
#' @export
studio_demo_inputs <- function(dict = STANDARD_DICT) {
  spec <- parse_input(GIVOSIRAN_SENSE_SPEC, dict = dict, notation = "structured")
  all_mets <- generate_metabolites(spec, list(max_3p = 3, max_5p = 3, endo = TRUE))
  ids <- c("M01", names(.STUDIO_DEMO_AMP))
  mets <- Filter(function(m) m$id %in% ids, all_mets)
  byid <- setNames(mets, vapply(mets, function(m) m$id, character(1)))
  ion_types <- c("a", "aB", "b", "bB", "w", "y")
  frag_cache <- lapply(byid, function(m) generate_fragments(m, dict, ion_types, 1:2, 0, FALSE))
  pmatch <- c(w = .78, y = .66, `a-B` = .58, b = .34, a = .30, `b-B` = .25)
  basei <- c(w = 1, y = .8, `a-B` = .7, b = .4, a = .35, `b-B` = .3)

  times <- .STUDIO_DEMO_TIMES
  sample_ids <- sprintf("RatLiver_%02dh_A", times)
  species <- rbind(data.frame(id = names(byid), k = 0L, stringsAsFactors = FALSE),
                   data.frame(id = "M01", k = 1L, stringsAsFactors = FALSE))
  rows <- list(); conf_rows <- list(); spectra <- list()
  for (si in seq_along(times)) {
    raw <- mapply(function(id, k) .studio_demo_signal(id, k, times[si]), species$id, species$k)
    pct <- 100 * raw / sum(raw)
    for (j in seq_len(nrow(species))) {
      id <- species$id[j]; k <- species$k[j]; met <- byid[[id]]
      info <- metabolite_mass_info(met, dict)
      mass <- if (k > 0) ps_oxid_mass(info$mono_mass, k) else info$mono_mass
      z0 <- max(3L, as.integer(round(mass / 1000)))
      theo_mz <- (mass - z0 * .PROTON) / z0
      seed <- .studio_hash(c(sample_ids[si], id, k))
      sim <- .studio_with_seed(seed, {
        ppm <- max(-5, min(5, rnorm(1, 0, 1.4)))
        level <- min(1, 0.35 + pct[j] / 14)
        peaks_mz <- numeric(0); peaks_int <- numeric(0); sites <- integer(0)
        for (f in frag_cache[[id]]) {
          ty <- f$ion_type
          for (i in seq_len(nrow(f$mz_table))) {
            zz <- f$mz_table$z[i]; mz <- f$mz_table$mz[i]
            p <- pmatch[[ty]] * level * (if (zz == 1) 1 else .8) * (if (mz > 1500) .7 else 1)
            if (runif(1) < p) {
              e <- max(-5.5, min(5.5, rnorm(1, 0, 1.9)))
              peaks_mz <- c(peaks_mz, mz * (1 + e * 1e-6))
              peaks_int <- c(peaks_int, basei[[ty]] * (0.25 + 0.75 * runif(1)) * (if (zz == 1) 1 else .6))
              sites <- c(sites, f$cleavage_site)
            }
          }
        }
        if (isTRUE(met$n_ps > 0) && !(met$kind %in% c("exo_5p", "endo_3frag", "endo_3frag_p")) && runif(1) < level) {
          peaks_mz <- c(peaks_mz, 94.9452 * (1 + rnorm(1, 0, 1e-6)), 192.9746 * (1 + rnorm(1, 0, 1e-6)))
          peaks_int <- c(peaks_int, .55 + .3 * runif(1), .45 + .3 * runif(1))
        }
        nn <- round(55 + 40 * (1 - level))
        peaks_mz <- c(peaks_mz, 120 + runif(nn) * 2000); peaks_int <- c(peaks_int, runif(nn)^3 * .22)
        list(ppm = ppm, spec = data.frame(mz = peaks_mz, intensity = peaks_int * 1e5),
             cov = length(unique(sites)) / max(1, met$n - 1), rt = 14.28 - (21 - met$n) * 0.07 - if (met$kind == "conj_catab") .15 else 0)
      })
      area <- 1e8 * raw[j] / sum(raw[species$id == "M01" & species$k == 0])
      rows[[length(rows) + 1]] <- data.frame(
        sample = sample_ids[si], met_id = id, met_name = met$name, kind = met$kind, k_oxid = k,
        z = z0, adduct = "H", theo_mz = theo_mz, obs_mz = theo_mz * (1 + sim$ppm * 1e-6),
        ppm_error = sim$ppm, rt = sim$rt, intensity = area / 40, area = area, stringsAsFactors = FALSE)
      conf_rows[[length(conf_rows) + 1]] <- data.frame(
        sample = sample_ids[si], met_id = id, met_name = met$name, k_oxid = k, z = z0, adduct = "H",
        n_ms2_peaks = nrow(sim$spec), n_frag_matches = NA_integer_, coverage = sim$cov,
        confirmation_score = NA_real_, n_diagnostics = NA_integer_, confident = NA, stringsAsFactors = FALSE)
      spectra[[.studio_key(sample_ids[si], id, k, z0, "H")]] <- sim$spec
    }
  }
  list(spec = spec, mets = mets, dict = dict,
       ms1_matches = do.call(rbind, rows), ms2_confirmations = do.call(rbind, conf_rows),
       ms2_spectra = spectra,
       sample_meta = data.frame(sample = sample_ids, timepoint = as.character(times), stringsAsFactors = FALSE))
}

## ---- UI --------------------------------------------------------------------
.studio_asset <- function(name, module_dir = NULL) {
  p <- if (!is.null(module_dir)) file.path(module_dir, "inst", "app", "www", name) else ""
  if (!file.exists(p)) p <- system.file(file.path("app", "www", name), package = "OligoMetProfiler")
  if (!nzchar(p) || !file.exists(p)) NULL else p
}

#' UI for the Studio tab
#'
#' @param module_dir Repository checkout root, when running from source.
#' @return A `shiny::tagList`.
#' @export
studio_ui <- function(module_dir = NULL) {
  css <- .studio_asset("studio.css", module_dir); js <- .studio_asset("studio.js", module_dir)
  shiny::tagList(
    shiny::tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    shiny::tags$link(rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Instrument+Sans:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap"),
    if (!is.null(css)) shiny::tags$style(shiny::HTML(paste(readLines(css, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))),
    if (!is.null(js)) shiny::tags$script(shiny::HTML(paste(readLines(js, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))),
    shiny::tags$div(id = "oms-studio", class = "oms", `data-theme` = "light",
                    shiny::tags$p(style = "padding:16px;color:#6c757d;", "Loading Studio…")))
}

## ---- Server ----------------------------------------------------------------
#' Server logic for the Studio tab
#'
#' Call once from the app's `server()`; `rv` is the app's reactiveValues
#' (needs `spec`, `mets`, `dict`, `batch_ms_results`, `sample_meta`).
#' @param input,output,session Shiny server arguments.
#' @param rv The app's reactiveValues.
#' @param tab Name of the Studio tab in `main_nav`.
#' @export
studio_server <- function(input, output, session, rv, tab = "Studio") {
  st <- shiny::reactiveValues(demo = NULL, use_demo = FALSE)
  .send <- function(type, obj) {
    session$sendCustomMessage(type, jsonlite::toJSON(obj, auto_unbox = TRUE, na = "null", null = "null", digits = NA))
  }
  .inputs <- shiny::reactive({
    if (isTRUE(st$use_demo) && !is.null(st$demo)) return(c(st$demo, list(mode = "demo")))
    if (is.null(rv$spec) || is.null(rv$mets)) return(NULL)
    br <- rv$batch_ms_results
    has <- !is.null(br) && !is.null(br$ms1_matches) && nrow(br$ms1_matches) > 0
    list(mode = if (has) "results" else "library", spec = rv$spec, mets = rv$mets, dict = rv$dict %||% STANDARD_DICT,
         ms1_matches = if (has) br$ms1_matches else NULL, ms2_confirmations = br$ms2_confirmations,
         ms2_spectra = br$ms2_spectra, sample_meta = rv$sample_meta)
  })

  shiny::observe({
    shiny::req(identical(input$main_nav, tab))
    inp <- .inputs()
    if (is.null(inp)) { .send("oms-studio-empty", list(has_library = FALSE)); return() }
    ttl <- if (identical(inp$mode, "demo")) "Givosiran sense strand" else {
      nm <- input$oligo_name; if (is.null(nm) || !nzchar(nm)) NULL else nm }
    pl <- tryCatch(studio_payload(inp$spec, inp$mets, inp$dict, inp$ms1_matches, inp$ms2_confirmations,
                                  inp$sample_meta, mode = inp$mode, title = ttl),
                   error = function(e) { message("Studio payload failed: ", conditionMessage(e)); NULL })
    if (is.null(pl)) .send("oms-studio-empty", list(has_library = TRUE, error = TRUE)) else .send("oms-studio-load", pl)
  })

  shiny::observeEvent(input$studio_pick, {
    pk <- input$studio_pick; inp <- .inputs(); shiny::req(pk$id, inp)
    parts <- strsplit(pk$id, "\\+ox")[[1]]; mid <- parts[1]; k <- if (length(parts) > 1) as.integer(parts[2]) else 0L
    mi <- which(vapply(inp$mets, function(m) m$id == mid, logical(1)))
    shiny::req(length(mi) == 1)
    met <- inp$mets[[mi]]; spectrum <- NULL; z <- NULL
    m <- inp$ms1_matches
    if (!is.null(m) && !is.null(pk$sample) && nzchar(pk$sample)) {
      if (!"k_oxid" %in% names(m)) m$k_oxid <- 0L
      sig <- if ("area" %in% names(m) && any(!is.na(m$area))) "area" else "intensity"
      d <- m[m$sample == pk$sample & m$met_id == mid & m$k_oxid == k, , drop = FALSE]
      if (nrow(d) > 0) {
        r <- d[which.max(replace(d[[sig]], is.na(d[[sig]]), -Inf)), ]
        z <- r$z; spectrum <- inp$ms2_spectra[[.studio_key(pk$sample, mid, k, r$z, r$adduct)]]
      }
    }
    det <- tryCatch(studio_detail(met, inp$dict, spectrum, k, z),
                    error = function(e) { message("Studio detail failed: ", conditionMessage(e)); NULL })
    if (!is.null(det)) { det$sample <- pk$sample %||% ""; det$id <- pk$id; .send("oms-studio-detail", det) }
  })

  shiny::observeEvent(input$studio_demo, {
    shiny::withProgress(message = "Building demo study", value = 0.2, {
      if (is.null(st$demo)) st$demo <- studio_demo_inputs(rv$dict %||% STANDARD_DICT)
      shiny::incProgress(0.8)
    })
    st$use_demo <- TRUE
  })
  shiny::observeEvent(input$studio_mine, { st$use_demo <- FALSE })
  shiny::observeEvent(input$studio_nav, {
    target <- input$studio_nav
    if (target %in% c("Library Generation", "Batch Processing", "Empirical MS2 Library"))
      shiny::updateTabsetPanel(session, "main_nav", selected = target)
  })
  invisible(NULL)
}
