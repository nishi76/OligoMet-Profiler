# =============================================================================
# run_library_figures.R
# Generate the predicted-library characterization figures for an oligonucleotide
# as publication-quality PNGs, using OligoMetProfiler's library-figure
# functions. A standalone counterpart to the "Library Figures" tab.
#
# The analyte name is set once, here, by you. Nothing downstream hardcodes it.
# =============================================================================

suppressPackageStartupMessages(library(OligoMetProfiler))

# ---- configure -------------------------------------------------------------
oligo_name <- "my_oligo"          # <- your analyte name; used in every title
output_dir <- "library_figures"   # where the PNGs are written

# Describe the oligo as a base/sugar/linkage triplet (see parse_input()), or
# load your own CSV exports and skip to the plotting section below. This
# example uses a 5-10-5 MOE/DNA gapmer with a phosphorothioate backbone.
triplet <- paste0(
  "[A,e,s][G,e,s][C,e,s][A,e,s][T,e,s]",   # 2'-MOE wing
  "[G,d,s][A,d,s][C,d,s][T,d,s][C,d,s][A,d,s][T,d,s][G,d,s][A,d,s][C,d,s]",  # DNA gap
  "[T,e,s][T,e,s][G,e,s][C,e,s][A,e,o]")   # 2'-MOE wing

# ---- build the library -----------------------------------------------------
dict <- build_dictionary()
spec <- parse_input(triplet, dict = dict)
mets <- generate_metabolites(spec, opts = list(oligo_name = oligo_name), dict = dict)
prm  <- prm_inclusion_list(mets, dict, z_range = 4:12, max_oxid = 6)
mf   <- library_mass_frame(mets, dict)

# ---- figures ---------------------------------------------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
figs <- list(
  mass_shift_waterfall = plot_mass_shift_waterfall(mf, oligo_name),
  mass_defect          = plot_mass_defect(mf, oligo_name),
  kendrick_mass_defect = plot_kendrick_defect(prm, mf, oligo_name),
  charge_heatmap       = plot_charge_heatmap(prm, oligo_name),
  prm_mz_distribution  = plot_prm_mz_distribution(prm, oligo_name),
  metabolite_classes   = plot_metabolite_categories(mf, oligo_name)
)
sizes <- list(
  mass_shift_waterfall = c(10, 5), mass_defect = c(8, 5), kendrick_mass_defect = c(10, 7),
  charge_heatmap = c(8, 10), prm_mz_distribution = c(9, 5), metabolite_classes = c(7, 4)
)
for (nm in names(figs)) {
  f <- file.path(output_dir, sprintf("%s_%s.png", oligo_name, nm))
  ggplot2::ggsave(f, figs[[nm]], width = sizes[[nm]][1], height = sizes[[nm]][2], dpi = 300)
  cat("  wrote", f, "\n")
}

# ---- reference table -------------------------------------------------------
utils::write.csv(mf, file.path(output_dir, sprintf("%s_metabolite_reference.csv", oligo_name)),
                 row.names = FALSE)
cat("\nDone. Metabolites:", nrow(mf), " | PRM targets:", nrow(prm), "\n")
