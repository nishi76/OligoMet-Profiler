# test_degradation_map.R -- position-resolved degradation map model
# (R/degradation_map.R): class/position assignment from the match table,
# coverage, parent-remaining, timepoint ordering, and the HTML widget build.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R", "degradation_map.R"))
  source(file.path(.pkg_root, "R", .f))
check <- function(label, ok) {
  cat(sprintf("%-74s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Degradation-map model suite ====\n\n")

## ---- helpers ----------------------------------------------------------------
check("timepoint parse: predose/hours/days/weeks", {
  h <- c(.dm_parse_hours("predose"), .dm_parse_hours("4h"), .dm_parse_hours("30 days"),
         .dm_parse_hours("1 week"), .dm_parse_hours("12"))
  all(abs(h - c(0, 4, 720, 168, 12)) < 1e-6)
})
check("timepoint ordering is numeric, not lexical", identical(
  .dm_order_levels(c("168h", "pre", "24h", "4h")), c("pre", "4h", "24h", "168h")))
check("span: 3'-exo N-2 on a 10-mer covers 1..8", identical(.dm_met_span("exo_3p", 2L, 10L), c(1L, 8L)))
check("span: 5'-exo N-3 covers 4..10", identical(.dm_met_span("exo_5p", 3L, 10L), c(4L, 10L)))
check("span: endo 5' frag @6 covers 1..6", identical(.dm_met_span("endo_5frag", 6L, 10L), c(1L, 6L)))
check("span: parent covers the whole oligo", identical(.dm_met_span("parent", NA, 10L), c(1L, 10L)))

## ---- build a synthetic matched batch ----------------------------------------
spec <- list(bases = strsplit("GCATGACTCC", "")[[1]],
             sugars = c(rep("e", 3), rep("d", 4), rep("e", 3)),
             linkages = c(rep("s", 9), NA), conj5 = "none", conj3 = "GalNAc3", n = 10)
mets <- generate_metabolites(spec, list(oligo_name = "OLM"))
ids <- vapply(mets, function(m) m$id, ""); kinds <- vapply(mets, function(m) m$kind, "")
parent_id <- ids[kinds == "parent"][1]
exo3 <- mets[kinds == "exo_3p"][[1]]
conj <- mets[kinds == "conj_catab"]

rows <- list()
addrow <- function(mid, kind, sample, area, mod = "none", mod_site = NA, ms2 = NA, nn = 10) {
  rows[[length(rows) + 1L]] <<- data.frame(
    met_id = mid, met_name = paste0(mid, "_", kind), kind = kind, n = nn,
    mod = mod, mod_site = mod_site, ppm_error = 1.1, ms2_score = ms2,
    area = area, intensity = area, sample = sample, stringsAsFactors = FALSE)
}
meta <- list()
for (mx in c("plasma", "csf")) for (tp in c("pre", "24h", "168h")) {
  s <- paste(mx, tp, sep = "_")
  meta[[length(meta) + 1L]] <- data.frame(sample = s, timepoint = tp, matrix = mx, stringsAsFactors = FALSE)
  pbase <- if (mx == "plasma") c(pre = 1e7, `24h` = 5e6, `168h` = 3e5)[[tp]]
           else c(pre = 1e7, `24h` = 9e6, `168h` = 6e6)[[tp]]
  addrow(parent_id, "parent", s, pbase)
  addrow(exo3$id, "exo_3p", s, pbase * 0.3, ms2 = 0.8)
  addrow(conj[[1]]$id, "conj_catab", s, pbase * 0.6)
  addrow(conj[[2]]$id, "conj_catab", s, pbase * 0.2)
  dm <- c(pre = 0, `24h` = 0.03, `168h` = 0.05)[[tp]]
  if (dm > 0) addrow(parent_id, "parent", s, pbase * dm, mod = "deamination", mod_site = 5, ms2 = 0.7)
  ox <- c(pre = 0, `24h` = 0.02, `168h` = 0.06)[[tp]]
  if (ox > 0) addrow(parent_id, "parent", s, pbase * ox, mod = "oxidation", mod_site = 2)
}
matches <- do.call(rbind, rows); sm <- do.call(rbind, meta)
model <- degradation_position_model(matches, mets, spec, sm, oligo_name = "OLM-10", min_pct = 0.1)

## ---- model assertions -------------------------------------------------------
check("times resolve and order pre -> 24h -> 168h", identical(model$times, c("pre", "24h", "168h")))
check("both matrices present", setequal(model$matrices, c("plasma", "csf")))
check("DNA gap detected (residues 4..7)", identical(as.integer(model$gap), c(4L, 7L)))
check("base cells carry oxidation@2 and deamination@5", {
  b <- model$base
  any(b$pos == 2 & b$cls == "ox") && any(b$pos == 5 & b$cls == "dA")
})
check("base-mod abundance grows pre<24h<168h (deamination@5, plasma)", {
  b <- model$base[model$base$matrix == "plasma" & model$base$pos == 5, ]
  p24 <- b$pct[b$time == "24h"]; p168 <- b$pct[b$time == "168h"]
  length(p24) == 1 && length(p168) == 1 && p168 > p24
})
check("deamination@5 is MS2-confirmed", all(model$base$ms2[model$base$cls == "dA"]))
check("cleavage lane holds the 3'-exo product", {
  cl <- model$cleave
  nrow(cl) > 0 && all(cl$cls[cl$cls != "e3"] %in% c("e3", "e5", "en")) && any(cl$cls == "e3")
})
check("conjugate track holds GalNAc catabolites", nrow(model$conj) > 0 && all(model$conj$cls == "gn"))
check("parent remaining: plasma collapses, CSF holds", {
  pl <- model$parent$pct[model$parent$matrix == "plasma" & model$parent$time == "168h"]
  cs <- model$parent$pct[model$parent$matrix == "csf" & model$parent$time == "168h"]
  pl < 10 && cs > 40
})
check("parent remaining starts at 100 at the earliest timepoint", {
  all(abs(model$parent$pct[model$parent$time == "pre"] - 100) < 1e-6)
})
check("min_pct filter drops sub-threshold base cells", {
  m2 <- degradation_position_model(matches, mets, spec, sm, min_pct = 4.5)
  all(m2$base$pct >= 4.5)
})
check("empty matches returns an empty model, no error", {
  m0 <- degradation_position_model(data.frame(), mets, spec, sm)
  nrow(m0$base) == 0 && nrow(m0$cleave) == 0 && length(m0$matrices) == 0
})
check("no sample metadata collapses to one timepoint/matrix", {
  m1 <- degradation_position_model(matches, mets, spec, sample_meta = NULL)
  identical(m1$matrices, "all") && identical(m1$times, "all")
})

## ---- HTML widget (needs jsonlite) -------------------------------------------
if (requireNamespace("jsonlite", quietly = TRUE)) {
  html <- degradation_map_html(model, id = "dmapT")
  check("widget HTML is self-contained and scoped under .dmap",
        grepl('class="dmap"', html, fixed = TRUE) && grepl("<style>", html, fixed = TRUE))
  check("widget injects the model JSON (sequence + cells)",
        grepl('"seq"', html, fixed = TRUE) && grepl('"base"', html, fixed = TRUE))
  check("widget builds one matrix tab per matrix",
        lengths(regmatches(html, gregexpr('class="tab"', html)))[1] == length(model$matrices))
  check("empty model yields the empty-state message, not a crash",
        grepl("No positioned degradation", degradation_map_html(
          degradation_position_model(data.frame(), mets, spec, sm))))
} else {
  cat("  (skipping HTML-widget checks: jsonlite not installed)\n")
}

cat("\nAll degradation-map checks passed.\n")
