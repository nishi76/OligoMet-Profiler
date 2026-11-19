# =============================================================================
# palettes.R
# One colour policy for every figure in the package:
#   * categorical data  -> the Okabe-Ito colourblind-safe qualitative palette
#   * continuous data    -> viridis (ggplot2's built-in scale_*_viridis_c/d)
# Every plot function draws its categorical colours from okabe_ito()/
# kind_colors() and its continuous fills from scale_*_viridis_c, so colour is
# consistent and colourblind-safe across the app, the report, and exports.
# =============================================================================

# Okabe & Ito (2008) qualitative palette. High-contrast colours first; yellow
# (low contrast on white) and black are placed last.
OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
               "#D55E00", "#CC79A7", "#F0E442", "#000000")

#' First `n` Okabe-Ito colours (categorical palette)
#'
#' Returns up to eight distinct colourblind-safe colours. For `n > 8` the
#' palette is interpolated (distinctness degrades; prefer folding rare
#' categories into "Other" past eight).
#'
#' @param n Number of colours.
#' @return A character vector of hex colours of length `max(1, n)`.
#' @export
okabe_ito <- function(n) {
  n <- max(1L, as.integer(n))
  if (n <= length(OKABE_ITO)) return(OKABE_ITO[seq_len(n)])
  grDevices::colorRampPalette(OKABE_ITO)(n)
}

# Canonical metabolite-class -> Okabe-Ito colour, so a given kind reads the
# same colour in every figure (library, degradation, report). The phosphate-
# terminus endo variants share their family's hue.
KIND_COLORS <- c(parent = "#000000", exo_3p = "#0072B2", exo_5p = "#E69F00",
                 endo_5frag = "#009E73", endo_5frag_p = "#009E73",
                 endo_3frag = "#CC79A7", endo_3frag_p = "#CC79A7",
                 conj_catab = "#D55E00", internal_standard = "#56B4E9")

#' Okabe-Ito colours for a set of metabolite kinds
#'
#' Maps each kind to its canonical colour (see `KIND_COLORS`); any kind not
#' in the canonical set is assigned a remaining Okabe-Ito colour.
#'
#' @param kinds Character vector of metabolite kinds.
#' @return A named character vector, one colour per unique kind.
#' @export
kind_colors <- function(kinds) {
  kinds <- unique(as.character(kinds))
  out <- unname(KIND_COLORS[kinds])
  miss <- which(is.na(out))
  if (length(miss) > 0) {
    used <- setdiff(OKABE_ITO, out[!is.na(out)])
    fill <- if (length(used) >= length(miss)) used[seq_along(miss)] else okabe_ito(length(miss))
    out[miss] <- fill[seq_along(miss)]
  }
  stats::setNames(out, kinds)
}
