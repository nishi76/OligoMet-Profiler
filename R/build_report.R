# =============================================================================
# build_report.R
# Full publication-ready report generation for oligonucleotide metabolite ID.
#
# Produces an HTML or PDF report via rmarkdown with:
#   - Input summary and chemistry details
#   - Metabolite library table
#   - Mass/formula validation
#   - Charge envelope plot
#   - Isotope pattern plot
#   - Truncation series plot
#   - PS oxidation series plot
#   - Fragment ion map
#   - MS matching results (if data provided)
#   - Methods section
#   - References
#
# Plots use ggplot2; saved as PNG per user preference.
# =============================================================================

## ---- Plot 1: Charge envelope ----------------------------------------------
plot_charge_envelope <- function(met, dict = STANDARD_DICT, z_range = 3:12,
                                  h_offset = 0, max_oxid = 3) {
  info <- metabolite_mass_info(met, dict)
  kmax <- min(met$n_ps, max_oxid)
  rows <- list()
  for (k in 0:kmax) {
    mass <- ps_oxid_mass(info$mono_mass, k)
    for (z in z_range) {
      rows[[length(rows) + 1]] <- data.frame(
        z = z, mz = (mass + h_offset - z * .PROTON) / z,
        k_oxid = k, stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  df$k_oxid <- factor(df$k_oxid, levels = 0:kmax,
                       labels = paste0("k=", 0:kmax))
  ggplot2::ggplot(df, ggplot2::aes(x = z, y = mz, color = k_oxid)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = grDevices::colorRampPalette(
      c("#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C"))(kmax + 1)) +
    ggplot2::labs(x = "Charge State (z)", y = "m/z",
                  color = "PS→PO\nOxidation") +
    ggplot2::theme_minimal(base_size = 11,
                           base_family = "Liberation Sans") +
    ggplot2::theme(legend.position = "right")
}

## ---- Plot 2: Isotope pattern ----------------------------------------------
plot_isotope_pattern <- function(met, dict = STANDARD_DICT, z = 6,
                                  n_iso = 8, h_offset = 0,
                                  use_envipat = TRUE) {
  info <- metabolite_mass_info(met, dict)
  cl <- isotope_mz_cluster(info$formula_vec, z, n_top = n_iso,
                            h_offset = h_offset, use_envipat = use_envipat)
  if (is.null(cl)) return(NULL)
  cl$iso <- factor(cl$iso)
  ggplot2::ggplot(cl, ggplot2::aes(x = mz, y = abundance * 100)) +
    ggplot2::geom_col(fill = "#0279EE", width = 0.003 * max(cl$mz)) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", abundance * 100)),
                       vjust = -0.3, size = 3) +
    ggplot2::labs(x = sprintf("m/z (z=%d)", z), y = "Relative Abundance (%)",
                  title = sprintf("Isotope Pattern: %s", met$name)) +
    ggplot2::theme_minimal(base_size = 11,
                           base_family = "Liberation Sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
}

## ---- Plot 3: Truncation series --------------------------------------------
plot_truncation_series <- function(mets, dict = STANDARD_DICT) {
  rows <- list()
  for (met in mets) {
    if (met$kind %in% c("parent", "exo_3p", "exo_5p")) {
      info <- metabolite_mass_info(met, dict)
      rows[[length(rows) + 1]] <- data.frame(
        name = met$name, kind = met$kind, n = met$n,
        mass = info$mono_mass, stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  if (nrow(df) == 0) return(NULL)
  df$series <- ifelse(df$kind == "exo_3p", "3' Exonuclease",
                      ifelse(df$kind == "exo_5p", "5' Exonuclease", "Parent"))
  ggplot2::ggplot(df, ggplot2::aes(x = n, y = mass, color = series)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(values = c("Parent" = "#333333",
      "3' Exonuclease" = "#0279EE", "5' Exonuclease" = "#FF9400")) +
    ggplot2::labs(x = "Length (nt)", y = "Monoisotopic Mass (Da)",
                  color = "Series") +
    ggplot2::theme_minimal(base_size = 11,
                           base_family = "Liberation Sans")
}

## ---- Plot 4: PS oxidation mass shift --------------------------------------
plot_oxidation_series <- function(met, dict = STANDARD_DICT, max_oxid = 6) {
  series <- ps_oxidation_series(met, max_oxid = max_oxid, dict = dict)
  series$delta <- series$mono_mass - series$mono_mass[1]
  ggplot2::ggplot(series, ggplot2::aes(x = k, y = delta)) +
    ggplot2::geom_col(fill = "#0279EE", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", delta)),
                       vjust = -0.3, size = 3.5) +
    ggplot2::labs(x = "PS→PO Oxidation Events (k)",
                  y = "Mass Shift (Da)",
                  title = sprintf("PS Oxidation Series: %s", met$name)) +
    ggplot2::theme_minimal(base_size = 11,
                           base_family = "Liberation Sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
}

## ---- Plot 5: Fragment ion map ---------------------------------------------
plot_fragment_map <- function(met, dict = STANDARD_DICT, z_range = 1:2) {
  frags <- generate_fragments(met, dict, z_range = z_range)
  if (length(frags) == 0) return(NULL)
  rows <- do.call(rbind, lapply(frags, function(f) {
    data.frame(ion_type = f$ion_type, direction = f$direction,
               cleavage_site = f$cleavage_site,
               frag_length = f$frag_length,
               mass = f$mono_mass, stringsAsFactors = FALSE)
  }))
  rows$ion_type <- factor(rows$ion_type,
                           levels = c("a", "a-B", "b", "b-B", "c", "d",
                                      "w", "x", "y", "z"))
  ggplot2::ggplot(rows, ggplot2::aes(x = cleavage_site, y = ion_type,
                                      fill = mass)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient(low = "#ECE9E2", high = "#0279EE",
                                  name = "Mass (Da)") +
    ggplot2::labs(x = "Cleavage Site", y = "Ion Type",
                  title = sprintf("Fragment Ion Map: %s", met$name)) +
    ggplot2::theme_minimal(base_size = 11,
                           base_family = "Liberation Sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
}

## ---- Save plots to files ---------------------------------------------------
save_plots <- function(spec, mets, dict, opts, output_dir = tempdir()) {
  plots <- list()
  z_range <- opts$z_range %||% 3:12
  h_offset <- opts$h_offset %||% 0
  max_oxid <- opts$max_oxid %||% 6
  use_envipat <- opts$use_envipat %||% TRUE
  parent <- mets[[1]]

  # Charge envelope
  p <- plot_charge_envelope(parent, dict, z_range, h_offset, max_oxid)
  f <- file.path(output_dir, "plot_charge_envelope.png")
  ggplot2::ggsave(f, p, width = 7, height = 4, dpi = 150, bg = "white")
  plots$charge_envelope <- f

  # Truncation series
  p <- plot_truncation_series(mets, dict)
  if (!is.null(p)) {
    f <- file.path(output_dir, "plot_truncation_series.png")
    ggplot2::ggsave(f, p, width = 7, height = 4, dpi = 150, bg = "white")
    plots$truncation <- f
  }

  # PS oxidation
  p <- plot_oxidation_series(parent, dict, max_oxid)
  f <- file.path(output_dir, "plot_oxidation_series.png")
  ggplot2::ggsave(f, p, width = 6, height = 4, dpi = 150, bg = "white")
  plots$oxidation <- f

  # Fragment map
  p <- plot_fragment_map(parent, dict)
  if (!is.null(p)) {
    f <- file.path(output_dir, "plot_fragment_map.png")
    ggplot2::ggsave(f, p, width = 7, height = 4, dpi = 150, bg = "white")
    plots$fragment_map <- f
  }

  # Isotope pattern (for a representative charge state)
  z_rep <- z_range[which.min(abs(z_range - median(z_range)))]
  p <- plot_isotope_pattern(parent, dict, z = z_rep, h_offset = h_offset,
                             use_envipat = use_envipat)
  if (!is.null(p)) {
    f <- file.path(output_dir, "plot_isotope_pattern.png")
    ggplot2::ggsave(f, p, width = 7, height = 3.5, dpi = 150, bg = "white")
    plots$isotope <- f
  }

  plots
}

## ---- Build the full report -------------------------------------------------
# Generate an HTML or PDF report via rmarkdown.
# output_format: "html" or "pdf" (pdf requires LaTeX)
build_report <- function(spec, mets, dict, validation, ms_results = NULL,
                          ms_info = NULL, opts = list(),
                          output_file = "oligo_metabolite_report",
                          output_format = c("html", "pdf"),
                          plot_dir = tempdir()) {
  output_format <- match.arg(output_format)
  plots <- save_plots(spec, mets, dict, opts, plot_dir)

  # Build report content
  parent <- mets[[1]]
  parent_info <- metabolite_mass_info(parent, dict)
  mtab <- metabolite_table(mets)

  # Create rmarkdown content
  md <- character()

  # Title
  md <- c(md, "---",
    "title: 'Oligonucleotide Metabolite Identification Report'",
    paste0("date: '", format(Sys.Date(), "%B %d, %Y"), "'"),
    "output:",
    if (output_format == "pdf") {
      c("  pdf_document:", "    toc: true", "    toc_depth: 3",
        "    fig_caption: yes")
    } else {
      c("  html_document:", "    toc: true", "    toc_depth: 3",
        "    toc_float: true", "    theme: flatly", "    code_folding: hide")
    },
    "---", "")

  # 1. Summary
  md <- c(md, "## 1. Summary", "",
    paste0("**Oligonucleotide:** ", spec$raw %||% "N/A"), "",
    paste0("**Notation:** ", spec$notation), "",
    paste0("**Length:** ", spec$n, " nucleotides"), "",
    paste0("**Bases (5'→3'):** ", paste(spec$bases, collapse = "")), "",
    paste0("**Sugars:** ", paste(spec$sugars, collapse = "")), "",
    paste0("**Linkages:** ", paste(ifelse(is.na(spec$linkages), ".", spec$linkages), collapse = "")), "",
    paste0("**5' Conjugate:** ", spec$conj5, "  |  **3' Conjugate:** ", spec$conj3), "",
    paste0("**Parent Formula:** ", parent_info$formula_str), "",
    paste0("**Parent Monoisotopic Mass:** ", sprintf("%.4f Da", parent_info$mono_mass)), "",
    paste0("**Total Metabolites Generated:** ", length(mets)), "")

  # 2. Validation
  if (!is.null(validation)) {
    md <- c(md, "## 2. Formula & Mass Validation", "",
      paste0("| Check | Result |",
             "|-------|--------|",
             paste0("| Formula | ", validation$formula, " |"),
             paste0("| Monoisotopic Mass | ", sprintf("%.6f Da", validation$mass), " |"),
             paste0("| Mass Error | ", sprintf("%.4f ppm", validation$ppm), " |"),
             paste0("| Status | ", ifelse(validation$ok, "PASS (<1 ppm)", "FAIL"), " |")),
      "")
  }

  # 3. Metabolite Library
  md <- c(md, "## 3. Metabolite Library", "",
    "The following metabolites were generated from the parent oligonucleotide:", "")
  # Table header
  md <- c(md, "| ID | Name | Kind | Length | n_PS | Formula | Mono Mass |",
    "|----|------|------|--------|------|---------|-----------|")
  for (met in mets) {
    info <- metabolite_mass_info(met, dict)
    md <- c(md, sprintf("| %s | %s | %s | %d | %d | %s | %.2f |",
      met$id, met$name, met$kind, met$n, met$n_ps,
      info$formula_str, info$mono_mass))
  }
  md <- c(md, "")

  # 4. Charge Envelope
  md <- c(md, "## 4. Charge State Envelope", "",
    "Theoretical m/z values for each charge state and PS oxidation level:", "")
  if (!is.null(plots$charge_envelope))
    md <- c(md, paste0("![](", basename(plots$charge_envelope), ")"), "")
  if (!is.null(plots$isotope))
    md <- c(md, "### Isotope Pattern", "",
      paste0("![](", basename(plots$isotope), ")"), "")

  # 5. Truncation Series
  md <- c(md, "## 5. Exonuclease Truncation Series", "",
    "3' and 5' exonuclease truncation series with monoisotopic masses:", "")
  if (!is.null(plots$truncation))
    md <- c(md, paste0("![](", basename(plots$truncation), ")"), "")

  # 6. PS Oxidation
  md <- c(md, "## 6. PS→PO Oxidation Series", "",
    "Phosphorothioate to phosphodiester oxidation (desulfurization) results in a mass shift of -15.977 Da per event:", "")
  if (!is.null(plots$oxidation))
    md <- c(md, paste0("![](", basename(plots$oxidation), ")"), "")

  # 7. Fragment Ions
  md <- c(md, "## 7. MS/MS Fragment Ions (McLuckey)", "",
    "Theoretical McLuckey fragment ions for sequence confirmation:", "")
  if (!is.null(plots$fragment_map))
    md <- c(md, paste0("![](", basename(plots$fragment_map), ")"), "")

  # 8. MS Matching Results
  if (!is.null(ms_results) && !is.null(ms_results$summary) &&
      nrow(ms_results$summary) > 0) {
    md <- c(md, "## 8. MS Data Matching Results", "",
      "### MS1 Matches", "")
    ms1 <- ms_results$ms1_matches
    if (nrow(ms1) > 0) {
      md <- c(md, "| Met Name | z | Theo m/z | Obs m/z | ppm Error |",
        "|----------|---|----------|---------|-----------|")
      for (i in seq_len(min(20, nrow(ms1)))) {
        md <- c(md, sprintf("| %s | %d | %.4f | %.4f | %.2f |",
          ms1$met_name[i], ms1$z[i], ms1$theo_mz[i],
          ms1$obs_mz[i], ms1$ppm_error[i]))
      }
      md <- c(md, "")
    }
    md <- c(md, "### Annotation Summary", "")
    s <- ms_results$summary
    md <- c(md, "| Met Name | MS1 Matches | Best ppm | MS2 Score | Coverage | Confident |",
      "|----------|------------|----------|-----------|----------|------------|")
    for (i in seq_len(nrow(s))) {
      md <- c(md, sprintf("| %s | %d | %.2f | %s | %s | %s |",
        s$met_name[i], s$n_ms1_matches[i], s$best_ppm[i],
        ifelse(is.na(s$ms2_score[i]), "-", as.character(s$ms2_score[i])),
        ifelse(is.na(s$ms2_coverage[i]), "-",
               sprintf("%.0f%%", s$ms2_coverage[i] * 100)),
        ifelse(s$confident[i], "Yes", "No")))
    }
    md <- c(md, "")
  }

  # 9. Methods
  md <- c(md, "## 9. Methods", "",
    "### Chemistry Dictionary", "",
    "Molecular formulas for nucleobases, sugars, linkages, and conjugates were defined from standard atomic compositions and validated against the published molecular formulas of approved oligonucleotide therapeutics (nusinersen, inotersen) to <1 ppm mass accuracy.", "",
    "### Metabolite Generation", "",
    "The theoretical metabolite library was generated by simulating:",
    "- 3' exonuclease truncation (n-1 through n-k)",
    "- 5' exonuclease truncation (n-1 through n-k)",
    "- Endonuclease internal cleavage (single-cut fragments)",
    "- PS→PO oxidation series (0 through k events, -15.977 Da each)", "",
    "### Mass Calculation", "",
    "Monoisotopic masses were computed from atomic masses (IUPAC). Charge state envelopes were calculated for negative ESI mode: [M - zH]^z- with m/z = (M - z * proton) / z. Isotope patterns were computed via enviPat or built-in convolution.", "",
    "### MS/MS Fragmentation", "",
    "McLuckey fragment ions (a, a-B, b, c, w, x, y) and internal fragments (w-a, w-b) were generated for each metabolite. PS diagnostic ions at m/z 94.9452 and 192.9746 were checked. Sequence coverage was computed as the fraction of backbone cleavage sites confirmed by matched fragments. The confirmation score combines coverage (0-50 pts), match count (0-25 pts), mass accuracy (0-15 pts), and diagnostic ion presence (0-10 pts).", "",
    "### References", "",
    "1. McLuckey SA, Van Berkel GJ, Glish GL. Tandem mass spectrometry of small, multiply charged oligonucleotides. J Am Soc Mass Spectrom. 1992;3(1):60-72.",
    "2. Kim J et al. Metabolite Profiling of the Antisense Oligonucleotide Eluforsen Using LC-MS. Mol Ther Nucleic Acids. 2019;17:714-725.",
    "3. Liu Y et al. OligoDistiller: A Platform Agnostic Software Tool for MS and MS2 Data Analysis. Anal Chem. 2025.",
    "4. Ye et al. Automatic identification of oligonucleotide metabolites using FMVS. J Chromatogr B. 2025.",
    "")

  # Write rmd file
  rmd_file <- file.path(plot_dir, paste0(output_file, ".Rmd"))
  writeLines(md, rmd_file)

  # Render
  out_format <- if (output_format == "pdf") {
    rmarkdown::pdf_document(toc = TRUE, toc_depth = 3)
  } else {
    rmarkdown::html_document(toc = TRUE, toc_depth = 3,
                             toc_float = TRUE, theme = "flatly")
  }

  tryCatch({
    rmarkdown::render(rmd_file, output_format = out_format,
                      output_dir = plot_dir, quiet = TRUE)
    ext <- if (output_format == "pdf") ".pdf" else ".html"
    out_path <- file.path(plot_dir, paste0(output_file, ext))
    cat("Report saved to", out_path, "\n")
    out_path
  }, error = function(e) {
    cat("Report rendering failed:", conditionMessage(e), "\n")
    cat("Rmd source saved at", rmd_file, "\n")
    rmd_file
  })
}
