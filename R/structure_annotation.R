# =============================================================================
# structure_annotation.R
# Sequence-level structure annotation for an oligonucleotide metabolite --
# the "cleavage map" view (as in Sciex Molecule Profiler / Thermo BioPharma
# Finder): the parent sequence drawn as a residue chain, with the
# metabolite's retained span highlighted, the lost portion faded, the
# cleavage site marked, the 5'/3' termini labelled, the backbone drawn as
# phosphorothioate vs phosphodiester, each residue coloured by its 2'-sugar
# chemistry, and any in-place modification (depurination/abasic, +O
# oxidation, deamination, 2',3'-cyclic phosphate) badged at its position.
#
# Two layers:
#   metabolite_annotation_model()  -- the tidy model (testable, no drawing)
#   annotate_metabolite_svg()      -- renders the model to an inline SVG
#
# The SVG is self-contained (no external assets) so it drops straight into a
# Shiny htmlOutput or an HTML report.
# =============================================================================

## ---- Sugar-class colouring (colourblind-safe, Okabe-Ito) -------------------
# Border colour per 2'-sugar class; the box fill is a light tint of it.
.SUGAR_CLASS <- function(code) {
  switch(as.character(code),
    d = "DNA", r = "RNA",
    m = "2'-OMe", OMe = "2'-OMe",
    f = "2'-F",
    e = "2'-MOE", MOE = "2'-MOE",
    cEt = "locked", LNA = "locked", ENA = "locked",
    "other")
}
.SUGAR_COLORS <- c(
  DNA      = "#E69F00",  # orange
  RNA      = "#56B4E9",  # sky blue
  `2'-OMe` = "#009E73",  # green
  `2'-F`   = "#D55E00",  # vermilion
  `2'-MOE` = "#0072B2",  # blue
  locked   = "#CC79A7",  # purple
  other    = "#999999"   # grey
)

## ---- Terminus labels -------------------------------------------------------
.term5_label <- function(conj5) {
  conj5 <- conj5 %||% "none"
  if (conj5 == "none") "5'-OH"
  else if (conj5 == "5'-phosphate") "5'-p"
  else if (conj5 == "5'-thiophosphate") "5'-ps"
  else paste0("5'-", conj5)
}
.term3_label <- function(conj3, cyclic = FALSE) {
  if (isTRUE(cyclic)) return("2',3'-cyclic-p")
  conj3 <- conj3 %||% "none"
  if (conj3 == "none") "3'-OH"
  else if (conj3 == "3'-phosphate") "3'-p"
  else if (conj3 == "3'-thiophosphate") "3'-ps"
  else paste0("3'-", conj3)
}

## ---- Fragment span within the parent ---------------------------------------
# Locate the metabolite's residues within the parent, from its kind + site.
.fragment_span <- function(met, N) {
  kind <- met$kind %||% "parent"
  k <- met$site
  if (is.null(k) || is.na(k) || kind == "parent") return(c(1L, N))
  k <- as.integer(k)
  if (kind == "exo_3p")                         return(c(1L, N - k))
  if (kind == "exo_5p")                         return(c(k + 1L, N))
  if (kind %in% c("endo_5frag", "endo_5frag_p")) return(c(1L, k))
  if (kind %in% c("endo_3frag", "endo_3frag_p")) return(c(k + 1L, N))
  c(1L, N)
}

## ---- Is the linkage at position i a phosphorothioate? ----------------------
.is_ps_code <- function(code, dict) {
  if (is.na(code) || !nzchar(code)) return(NA)
  f <- dict[[code]]$formula
  !is.null(f) && isTRUE(unname(.as_formula(f)[["S"]]) > 0)
}

## ---- Candidate modification sites (metabolite-local 1..n) ------------------
.mod_candidate_sites <- function(met, mod) {
  switch(mod,
    depurination     = which(met$bases %in% c("A", "G")),
    deamination      = which(met$bases %in% .DEAMINABLE_BASES),
    cyclic_phosphate = met$n,
    integer(0))
}

#' Build the sequence-annotation model for a metabolite
#'
#' Produces the tidy data behind the cleavage-map view: one row per parent
#' residue (base, 2'-sugar class and colour, backbone type to the next
#' residue, whether it is retained in this metabolite), the retained span,
#' the 5'/3' terminus labels, and any modification badges. Rendering is
#' separate (see [annotate_metabolite_svg()]).
#'
#' @param met A metabolite object (see [generate_metabolites()]).
#' @param parent Optional parent `oligo_spec` (see [parse_input()]) so the
#'   lost portion can be drawn faded; without it only the metabolite's own
#'   residues are modelled.
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @param mod In-place modification to badge: `"none"`, `"depurination"`,
#'   `"oxidation"`, `"deamination"`, or `"cyclic_phosphate"` (see
#'   [metabolite_mod_variants()]).
#' @param mod_site 1-indexed metabolite residue the modification sits on,
#'   or `NA` when MS1 could not localize it (then every candidate site is
#'   badged lightly).
#' @param k_oxid Number of PS -> PO desulfurizations to note (a count, not
#'   a site).
#' @return A list of class `oligomet_annotation`: `title` (name/kind/
#'   formula/mono_mass), `N`, `span`, `residues` (data.frame), `term5`,
#'   `term3`, `mods` (data.frame: `position`, `type`, `label`,
#'   `candidate`), and `notes` (character).
#' @seealso [annotate_metabolite_svg()]
#' @export
metabolite_annotation_model <- function(met, parent = NULL, dict = STANDARD_DICT,
                                        mod = "none", mod_site = NA_integer_,
                                        k_oxid = 0) {
  mod <- mod %||% "none"
  # Parent chain (full) if given, else just the metabolite's own residues.
  if (!is.null(parent)) {
    N <- parent$n
    bases <- parent$bases; sugars <- parent$sugars; linkages <- parent$linkages
    span <- .fragment_span(met, N)
  } else {
    N <- met$n
    bases <- met$bases; sugars <- met$sugars; linkages <- met$linkages
    span <- c(1L, N)
  }
  in_frag <- seq_len(N) >= span[1] & seq_len(N) <= span[2]

  sclass <- vapply(sugars, .SUGAR_CLASS, "")
  residues <- data.frame(
    position = seq_len(N),
    base = as.character(bases),
    sugar = as.character(sugars),
    sugar_class = sclass,
    color = unname(.SUGAR_COLORS[sclass]),
    in_frag = in_frag,
    is_ps_next = vapply(seq_len(N), function(i) {
      if (i >= N) return(NA)
      isTRUE(.is_ps_code(linkages[i], dict))
    }, logical(1)),
    stringsAsFactors = FALSE
  )

  info <- metabolite_mass_info(met, dict)
  title <- list(name = met$name %||% met$id, kind = met$kind %||% "parent",
                formula = info$formula_str, mono_mass = info$mono_mass, n = met$n)

  cyclic <- identical(mod, "cyclic_phosphate") && can_form_cyclic_phosphate(met)
  term5 <- .term5_label(met$conj5)
  term3 <- .term3_label(met$conj3, cyclic = cyclic)

  # Modification badges, mapped from metabolite-local to parent coordinates.
  mods_df <- data.frame(position = integer(0), type = character(0),
                        label = character(0), candidate = logical(0))
  notes <- character(0)
  to_parent <- function(local) span[1] + as.integer(local) - 1L
  if (!identical(mod, "none")) {
    if (identical(mod, "cyclic_phosphate") && !cyclic) {
      notes <- c(notes, "2',3'-cyclic phosphate not applicable (needs a 3'-phosphate terminus on a ribonucleotide)")
    } else if (identical(mod, "oxidation")) {
      notes <- c(notes, "+O oxidation (position not localized by MS1)")
    } else if (!is.na(mod_site)) {
      lab <- switch(mod,
        depurination = "abasic",
        deamination  = paste0(met$bases[mod_site], ">",
                              .DEAMINATION_PRODUCT[[met$bases[mod_site]]]),
        cyclic_phosphate = "cyclic-p", mod)
      mods_df <- rbind(mods_df, data.frame(position = to_parent(mod_site),
                       type = mod, label = lab, candidate = FALSE,
                       stringsAsFactors = FALSE))
    } else {
      sites <- .mod_candidate_sites(met, mod)
      for (s in sites) {
        lab <- switch(mod,
          depurination = "abasic?",
          deamination  = paste0(met$bases[s], ">",
                               .DEAMINATION_PRODUCT[[met$bases[s]]], "?"),
          cyclic_phosphate = "cyclic-p", mod)
        mods_df <- rbind(mods_df, data.frame(position = to_parent(s),
                         type = mod, label = lab, candidate = TRUE,
                         stringsAsFactors = FALSE))
      }
      if (nrow(mods_df) > 0)
        notes <- c(notes, sprintf("%s site not localized by MS1 (%d candidate position%s)",
                                  mod, nrow(mods_df), if (nrow(mods_df) == 1) "" else "s"))
      else
        notes <- c(notes, sprintf("no applicable %s site in this metabolite", mod))
    }
  }
  if (!is.na(k_oxid) && k_oxid > 0)
    notes <- c(notes, sprintf("-%d x S (PS->PO desulfurization; positions not localized)", k_oxid))

  structure(list(title = title, N = N, span = span, residues = residues,
                 term5 = term5, term3 = term3, mods = mods_df, notes = notes),
            class = "oligomet_annotation")
}

## ---- SVG rendering ---------------------------------------------------------
.xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

#' Render a metabolite's sequence annotation as an inline SVG
#'
#' Draws the model from [metabolite_annotation_model()] as a cleavage map:
#' the parent sequence as a chain of residue boxes coloured by 2'-sugar
#' chemistry, the retained metabolite span solid and the lost portion
#' faded, the backbone drawn thick/dark for phosphorothioate and thin/grey
#' for phosphodiester, a scissors at the cleavage boundary, the 5'/3'
#' terminus labels, and modification badges above their residue. The
#' returned string is a complete `<svg>` element with no external assets.
#'
#' @param x A metabolite object (see [generate_metabolites()]) or a model
#'   from [metabolite_annotation_model()].
#' @param parent,dict,mod,mod_site,k_oxid Passed to
#'   [metabolite_annotation_model()] when `x` is a metabolite.
#' @param per_row Residues per row before wrapping (default 25).
#' @return A single character string: an `<svg>` element.
#' @seealso [metabolite_annotation_model()]
#' @export
annotate_metabolite_svg <- function(x, parent = NULL, dict = STANDARD_DICT,
                                    mod = "none", mod_site = NA_integer_,
                                    k_oxid = 0, per_row = 25) {
  m <- if (inherits(x, "oligomet_annotation")) x else
    metabolite_annotation_model(x, parent, dict, mod, mod_site, k_oxid)

  cw <- 30; ch <- 30; gap <- 8; padx <- 16
  row_h <- ch + 46            # residue box + badge space + ruler
  n_rows <- ceiling(m$N / per_row)
  width <- padx * 2 + per_row * (cw + gap)
  head_h <- 54
  legend_h <- 26
  height <- head_h + n_rows * row_h + legend_h + 16

  pos_xy <- function(p) {
    r <- (p - 1L) %/% per_row
    c_in <- (p - 1L) %% per_row
    x0 <- padx + c_in * (cw + gap)
    y0 <- head_h + r * row_h + 20
    c(x = x0, y = y0)
  }

  el <- character(0)
  add <- function(...) el[[length(el) + 1L]] <<- paste0(...)

  # Header
  t <- m$title
  add(sprintf('<text x="%d" y="22" font-size="15" font-weight="700" fill="#111">%s</text>',
              padx, .xml_escape(t$name)))
  add(sprintf('<text x="%d" y="40" font-size="11" fill="#555">%s &#183; %s &#183; %.3f Da &#183; n=%d</text>',
              padx, .xml_escape(t$kind), .xml_escape(t$formula), t$mono_mass, t$n))

  # Backbone connectors first (under the boxes)
  for (i in seq_len(m$N - 1L)) {
    a <- pos_xy(i); b <- pos_xy(i + 1L)
    if (a[["y"]] != b[["y"]]) next   # skip wrap boundaries
    both_in <- m$residues$in_frag[i] && m$residues$in_frag[i + 1L]
    ps <- isTRUE(m$residues$is_ps_next[i])
    col <- if (!both_in) "#d0d0d0" else if (ps) "#333333" else "#9aa0a6"
    wdt <- if (ps) 3.2 else 1.6
    y <- a[["y"]] + ch / 2
    add(sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%.1f"%s/>',
                a[["x"]] + cw, y, b[["x"]], y, col, wdt,
                if (!ps && both_in) ' stroke-dasharray="2,2"' else ""))
  }

  # Residue boxes
  for (i in seq_len(m$N)) {
    xy <- pos_xy(i); x0 <- xy[["x"]]; y0 <- xy[["y"]]
    rr <- m$residues[i, ]
    col <- rr$color
    if (rr$in_frag) {
      add(sprintf('<rect x="%.1f" y="%.1f" width="%d" height="%d" rx="6" fill="%s" fill-opacity="0.18" stroke="%s" stroke-width="2"/>',
                  x0, y0, cw, ch, col, col))
      add(sprintf('<text x="%.1f" y="%.1f" font-size="15" font-weight="700" text-anchor="middle" fill="#111">%s</text>',
                  x0 + cw / 2, y0 + ch / 2 + 5, .xml_escape(rr$base)))
    } else {
      add(sprintf('<rect x="%.1f" y="%.1f" width="%d" height="%d" rx="6" fill="none" stroke="#c8c8c8" stroke-width="1.3" stroke-dasharray="3,2"/>',
                  x0, y0, cw, ch))
      add(sprintf('<text x="%.1f" y="%.1f" font-size="14" text-anchor="middle" fill="#bbb">%s</text>',
                  x0 + cw / 2, y0 + ch / 2 + 5, .xml_escape(rr$base)))
    }
    # position ruler every 5th and at the ends
    if (i %% 5 == 0 || i == 1 || i == m$N) {
      add(sprintf('<text x="%.1f" y="%.1f" font-size="9" text-anchor="middle" fill="#999">%d</text>',
                  x0 + cw / 2, y0 + ch + 12, i))
    }
  }

  # Modification badges
  if (nrow(m$mods) > 0) {
    for (j in seq_len(nrow(m$mods))) {
      md <- m$mods[j, ]
      xy <- pos_xy(md$position); x0 <- xy[["x"]]; y0 <- xy[["y"]]
      fill <- if (isTRUE(md$candidate)) "#bbbbbb" else "#CC79A7"
      by <- y0 - 8
      add(sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2"/>',
                  x0 + cw / 2, y0, x0 + cw / 2, by, fill))
      add(sprintf('<circle cx="%.1f" cy="%.1f" r="3.2" fill="%s"/>', x0 + cw / 2, by, fill))
      add(sprintf('<text x="%.1f" y="%.1f" font-size="9" text-anchor="middle" fill="%s" font-weight="600">%s</text>',
                  x0 + cw / 2, by - 5, fill, .xml_escape(md$label)))
    }
  }

  # 5' / 3' terminus labels at the fragment ends
  s5 <- pos_xy(m$span[1]); s3 <- pos_xy(m$span[2])
  add(sprintf('<text x="%.1f" y="%.1f" font-size="10" text-anchor="end" fill="#0072B2" font-weight="700">%s</text>',
              s5[["x"]] - 3, s5[["y"]] + ch / 2 + 4, .xml_escape(m$term5)))
  add(sprintf('<text x="%.1f" y="%.1f" font-size="10" text-anchor="start" fill="#0072B2" font-weight="700">%s</text>',
              s3[["x"]] + cw + 3, s3[["y"]] + ch / 2 + 4, .xml_escape(m$term3)))

  # Scissors at cleavage boundaries (where the retained span meets the lost part)
  scissor <- function(after_pos) {
    a <- pos_xy(after_pos); b <- pos_xy(after_pos + 1L)
    if (a[["y"]] != b[["y"]]) return()
    xm <- (a[["x"]] + cw + b[["x"]]) / 2
    y <- a[["y"]] + ch / 2
    add(sprintf('<text x="%.1f" y="%.1f" font-size="15" text-anchor="middle" fill="#D55E00">&#9986;</text>',
                xm, y + 5))
  }
  if (m$span[1] > 1) scissor(m$span[1] - 1L)
  if (m$span[2] < m$N) scissor(m$span[2])

  # Legend
  ly <- head_h + n_rows * row_h + 2
  lx <- padx
  classes <- names(.SUGAR_COLORS)
  present <- unique(m$residues$sugar_class)
  classes <- classes[classes %in% present]
  for (cl in classes) {
    add(sprintf('<rect x="%.1f" y="%.1f" width="12" height="12" rx="3" fill="%s" fill-opacity="0.18" stroke="%s" stroke-width="1.6"/>',
                lx, ly, .SUGAR_COLORS[[cl]], .SUGAR_COLORS[[cl]]))
    add(sprintf('<text x="%.1f" y="%.1f" font-size="10" fill="#555">%s</text>', lx + 16, ly + 10, .xml_escape(cl)))
    lx <- lx + 20 + nchar(cl) * 6.6 + 14
  }
  add(sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#333" stroke-width="3.2"/>', lx, ly + 6, lx + 16, ly + 6))
  add(sprintf('<text x="%.1f" y="%.1f" font-size="10" fill="#555">PS</text>', lx + 20, ly + 10))
  lx <- lx + 44
  add(sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#9aa0a6" stroke-width="1.6" stroke-dasharray="2,2"/>', lx, ly + 6, lx + 16, ly + 6))
  add(sprintf('<text x="%.1f" y="%.1f" font-size="10" fill="#555">PO</text>', lx + 20, ly + 10))

  # Notes
  if (length(m$notes) > 0) {
    add(sprintf('<text x="%d" y="%.1f" font-size="10" fill="#777">%s</text>',
                padx, height - 4, .xml_escape(paste(m$notes, collapse = "   |   "))))
  }

  paste0(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%" font-family="system-ui,Segoe UI,Helvetica,Arial,sans-serif">',
            width, as.integer(height)),
    sprintf('<rect x="0" y="0" width="%d" height="%d" fill="#ffffff"/>', width, as.integer(height)),
    paste(el, collapse = ""),
    "</svg>")
}
