# =============================================================================
# library_figures.R
# Characterization figures for the predicted metabolite library, driven by the
# live library (not a CSV round-trip): a mass-shift waterfall, a mass-defect
# plot, a phosphorothioate-base Kendrick mass-defect plot (which aligns the
# PS->PO desulfurization series horizontally), a charge-state m/z heatmap, the
# PRM target m/z distribution, and a metabolite-class count. These complement
# the app's existing truncation, oxidation and charge-envelope plots rather
# than repeat them.
#
# Every title and label takes the oligo name from the caller (the app passes
# input$oligo_name); nothing here hardcodes an analyte name.
# =============================================================================

# Generic metabolite-class colors (Okabe-Ito, colorblind-safe) and labels. The
# phosphate-terminus endo variants share their family's hue; conjugate
# catabolites get their own.
.LIBFIG_COL <- c(parent = "#000000", exo_3p = "#0072B2", exo_5p = "#E69F00",
                 endo_5frag = "#009E73", endo_5frag_p = "#009E73",
                 endo_3frag = "#CC79A7", endo_3frag_p = "#CC79A7",
                 conj_catab = "#D55E00")
.LIBFIG_LAB <- c(parent = "Parent (full-length)", exo_3p = "3' Exonuclease",
                 exo_5p = "5' Exonuclease", endo_5frag = "Endo 5' fragment",
                 endo_5frag_p = "Endo 5' fragment (p)", endo_3frag = "Endo 3' fragment",
                 endo_3frag_p = "Endo 3' fragment (p)", conj_catab = "Conjugate catabolite")
# Exact S->O (desulfurization) mass loss per PS->PO conversion, and the PS
# backbone Kendrick base (nominal 63) used to align the oxidation series.
.LIBFIG_PSPO_DA <- 15.977156
.LIBFIG_KENDRICK_EXACT <- 62.945833
.LIBFIG_KENDRICK_NOMINAL <- 63

.libfig_theme <- function() {
  ggplot2::theme_bw(base_size = 11, base_family = "Liberation Sans") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", size = 12),
                   plot.subtitle = ggplot2::element_text(size = 9, color = "grey40"),
                   legend.position = "right",
                   strip.background = ggplot2::element_rect(fill = "grey92", colour = NA))
}
.libfig_kinds_present <- function(kinds) {
  k <- intersect(names(.LIBFIG_COL), unique(kinds))
  list(col = .LIBFIG_COL[k], lab = .LIBFIG_LAB[k])
}
.libfig_empty <- function(msg) ggplot2::ggplot() + ggplot2::theme_void() +
  ggplot2::labs(title = msg)

#' Assemble the metabolite mass table for the library figures
#'
#' One row per metabolite: `id`, `name`, `kind`, `n`, `n_ps`, `n_po`,
#' `formula`, `mono_mass`, `avg_mass`, and `mass_shift` (monoisotopic mass
#' minus the full-length parent's). Built from the live library, so the
#' figures need no CSV export.
#'
#' @param mets A metabolite library (see [generate_metabolites()]).
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A data.frame, or an empty data.frame when `mets` is empty.
#' @seealso [plot_mass_defect()], [plot_kendrick_defect()]
#' @export
library_mass_frame <- function(mets, dict = STANDARD_DICT) {
  if (length(mets) == 0) return(data.frame())
  rows <- lapply(mets, function(m) {
    mi <- tryCatch(metabolite_mass_info(m, dict), error = function(e) NULL)
    if (is.null(mi)) return(NULL)
    data.frame(id = m$id, name = m$name %||% m$id,
               kind = (m$kind_all %||% m$kind)[1], n = m$n,
               n_ps = m$n_ps %||% NA_integer_, n_po = m$n_po %||% NA_integer_,
               formula = mi$formula_str %||% NA_character_,
               mono_mass = mi$mono_mass, avg_mass = mi$avg_mass %||% NA_real_,
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(data.frame())
  df <- do.call(rbind, rows)
  pm <- df$mono_mass[df$kind == "parent"]
  df$mass_shift <- if (length(pm) >= 1 && is.finite(pm[1])) df$mono_mass - pm[1] else NA_real_
  df
}

#' Metabolite mass-shift waterfall
#' @param mass_frame A [library_mass_frame()] result.
#' @param oligo_name Analyte name for the title (from user input).
#' @return A ggplot.
#' @export
plot_mass_shift_waterfall <- function(mass_frame, oligo_name = "oligo") {
  if (is.null(mass_frame) || nrow(mass_frame) == 0) return(.libfig_empty("No metabolites to plot"))
  pal <- .libfig_kinds_present(mass_frame$kind)
  ggplot2::ggplot(mass_frame, ggplot2::aes(x = stats::reorder(.data$id, .data$mono_mass),
                                           y = .data$mass_shift, fill = .data$kind)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = pal$col, labels = pal$lab, name = NULL) +
    ggplot2::labs(title = sprintf("%s predicted metabolite mass shifts", oligo_name),
                  subtitle = sprintf("%d metabolites from the predicted library", nrow(mass_frame)),
                  x = "Metabolite ID", y = "Mass shift from parent (Da)") +
    .libfig_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7))
}

#' Mass-defect plot (fractional vs monoisotopic mass)
#' @inheritParams plot_mass_shift_waterfall
#' @return A ggplot.
#' @export
plot_mass_defect <- function(mass_frame, oligo_name = "oligo") {
  if (is.null(mass_frame) || nrow(mass_frame) == 0) return(.libfig_empty("No metabolites to plot"))
  d <- mass_frame
  d$mass_defect <- d$mono_mass - floor(d$mono_mass)
  pal <- .libfig_kinds_present(d$kind)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$mono_mass, y = .data$mass_defect, color = .data$kind)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::scale_color_manual(values = pal$col, labels = pal$lab, name = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = sprintf("Mass defect plot — %s metabolites", oligo_name),
                  subtitle = "Fractional mass vs monoisotopic mass; clusters track elemental composition",
                  x = "Monoisotopic mass (Da)", y = "Mass defect (Da)") +
    .libfig_theme()
}

#' Phosphorothioate-base Kendrick mass-defect plot
#'
#' Rescales mass onto a PS backbone Kendrick axis so the PS->PO
#' desulfurization series lines up horizontally, one panel per oxidation
#' level. Needs the PRM table (which carries the per-oxidation `k_oxid`).
#'
#' @param prm A [prm_inclusion_list()] result (`met_id`, `kind`, `k_oxid`, ...).
#' @param mass_frame A [library_mass_frame()] result (for `mono_mass`).
#' @param oligo_name Analyte name for the title (from user input).
#' @param max_oxid Highest oxidation level to panel (default 3).
#' @return A ggplot.
#' @export
plot_kendrick_defect <- function(prm, mass_frame, oligo_name = "oligo", max_oxid = 3) {
  if (is.null(prm) || nrow(prm) == 0 || is.null(mass_frame) || nrow(mass_frame) == 0)
    return(.libfig_empty("No PRM targets to plot"))
  d <- merge(prm, mass_frame[, c("id", "mono_mass")], by.x = "met_id", by.y = "id", all.x = TRUE)
  d <- d[is.finite(d$mono_mass) & d$k_oxid <= max_oxid, , drop = FALSE]
  if (nrow(d) == 0) return(.libfig_empty("No PRM targets to plot"))
  d$eff_mass <- d$mono_mass - d$k_oxid * .LIBFIG_PSPO_DA
  d$kendrick_mass <- d$eff_mass * (.LIBFIG_KENDRICK_NOMINAL / .LIBFIG_KENDRICK_EXACT)
  d$kendrick_defect <- d$kendrick_mass - round(d$kendrick_mass)
  pal <- .libfig_kinds_present(d$kind)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$eff_mass, y = .data$kendrick_defect, color = .data$kind)) +
    ggplot2::geom_point(size = 1.5, alpha = 0.6) +
    ggplot2::scale_color_manual(values = pal$col, labels = pal$lab, name = NULL) +
    ggplot2::facet_wrap(~k_oxid, labeller = ggplot2::labeller(
      k_oxid = function(x) paste0("PS→PO: ", x))) +
    ggplot2::labs(title = sprintf("Kendrick mass-defect plot — %s (PS base)", oligo_name),
                  subtitle = "Horizontal alignment marks the PS→PO desulfurization series; each panel = oxidation level",
                  x = "Effective monoisotopic mass (Da)", y = "Kendrick mass defect") +
    .libfig_theme()
}

#' Charge-state m/z heatmap (unoxidized metabolites)
#' @param prm A [prm_inclusion_list()] result.
#' @param oligo_name Analyte name for the title (from user input).
#' @param z_lo,z_hi Charge-state range to show (default 5-10).
#' @return A ggplot.
#' @export
plot_charge_heatmap <- function(prm, oligo_name = "oligo", z_lo = 5, z_hi = 10) {
  if (is.null(prm) || nrow(prm) == 0) return(.libfig_empty("No PRM targets to plot"))
  d <- prm[prm$k_oxid == 0 & prm$z >= z_lo & prm$z <= z_hi, , drop = FALSE]
  if (nrow(d) == 0) return(.libfig_empty("No PRM targets in the selected charge range"))
  d$met_label <- paste0(d$met_id, " (n=", d$n, ")")
  ggplot2::ggplot(d, ggplot2::aes(x = factor(.data$z),
                                  y = stats::reorder(.data$met_label, .data$precursor_mz),
                                  fill = .data$precursor_mz)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", .data$precursor_mz)), size = 2.2, color = "white") +
    ggplot2::scale_fill_viridis_c(name = "m/z", option = "D") +
    ggplot2::labs(title = sprintf("Charge-state m/z heatmap — %s (unoxidized)", oligo_name),
                  subtitle = sprintf("Precursor m/z per metabolite at z = %d–%d; a reference for EIC extraction", z_lo, z_hi),
                  x = "Charge state (z)", y = "Metabolite") +
    .libfig_theme() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
}

#' PRM target m/z distribution
#' @inheritParams plot_charge_heatmap
#' @return A ggplot.
#' @export
plot_prm_mz_distribution <- function(prm, oligo_name = "oligo") {
  if (is.null(prm) || nrow(prm) == 0) return(.libfig_empty("No PRM targets to plot"))
  pal <- .libfig_kinds_present(prm$kind)
  ggplot2::ggplot(prm, ggplot2::aes(x = .data$precursor_mz, fill = .data$kind)) +
    ggplot2::geom_histogram(bins = 80, alpha = 0.75, position = "identity", color = NA) +
    ggplot2::scale_fill_manual(values = pal$col, labels = pal$lab, name = NULL) +
    ggplot2::labs(title = sprintf("PRM target m/z distribution — %s", oligo_name),
                  subtitle = sprintf("%d targeted m/z across metabolites, oxidation levels and charge states", nrow(prm)),
                  x = "Precursor m/z", y = "Count") +
    .libfig_theme()
}

#' Metabolite-class count
#' @inheritParams plot_mass_shift_waterfall
#' @return A ggplot.
#' @export
plot_metabolite_categories <- function(mass_frame, oligo_name = "oligo") {
  if (is.null(mass_frame) || nrow(mass_frame) == 0) return(.libfig_empty("No metabolites to plot"))
  cc <- as.data.frame(table(kind = mass_frame$kind), stringsAsFactors = FALSE)
  cc$label <- .LIBFIG_LAB[cc$kind]; cc$label[is.na(cc$label)] <- cc$kind[is.na(cc$label)]
  pal <- .libfig_kinds_present(mass_frame$kind)
  ggplot2::ggplot(cc, ggplot2::aes(x = stats::reorder(.data$label, .data$Freq),
                                   y = .data$Freq, fill = .data$kind)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Freq), hjust = -0.2, size = 3.5) +
    ggplot2::scale_fill_manual(values = pal$col, guide = "none") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = sprintf("Metabolite classes — %s", oligo_name),
                  subtitle = "Predicted degradation products by cleavage mechanism",
                  x = NULL, y = "Number of metabolites") +
    .libfig_theme()
}
