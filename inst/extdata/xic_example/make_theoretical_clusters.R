# Writes theoretical_clusters.json for generate_xic_example.py: isotope
# clusters of the inotersen parent (M01) and its 3' N-1 metabolite (M02),
# computed with this package's own formula and isotope engine.
# Run from the repository root: Rscript inst/extdata/xic_example/make_theoretical_clusters.R
for (f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R", "mass_isotope.R",
            "fragments.R", "ms_matching.R", "batch_ms_processing.R", "xic_quant.R")) {
  source(file.path("R", f))
}
sp <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(sp, opts = list(oligo_name = "inotersen_example", max_3p = 3,
                                             max_5p = 3, endo = FALSE, min_frag_len = 3))
out <- lapply(mets[1:2], function(m) {
  info <- metabolite_mass_info(m, STANDARD_DICT)
  cl <- isotope_clusters(info$formula_vec, info$mono_mass, n = 12, use_envipat = FALSE)
  cl <- cl[order(cl$cluster), ]
  list(id = m$id, name = m$name, mono_mass = info$mono_mass,
       cluster = cl$cluster, mass = cl$mass, abundance = cl$abundance)
})
jsonlite::write_json(out, "inst/extdata/xic_example/theoretical_clusters.json",
                     auto_unbox = TRUE, digits = 10, pretty = TRUE)
