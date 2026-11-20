# Abundance normalization in OligoMet-Profiler

How the app expresses metabolite abundance, why it offers two bases, and the
published basis for each. This documents the methods behind the Degradation
Summary, the Degradation Map, and the time-course views.

## The shared input

Every metabolite is first reduced to one signal per sample: peak **area (AUC)**
when available, otherwise **max intensity**, taken as the **maximum** across the
species' charge states, adducts and PS→PO oxidation levels
(`.best_signal_per_met`, `R/degradation.R`). Calibration standards, QC and
blank samples are excluded from the biological summaries; internal standards
are neither parent nor degradant and are excluded from the species set.

## Two normalization bases, for two different questions

### 1. Within-sample composition (snapshot)

Each species as a percentage of the parent, or of the total drug-related
signal, **in the same injection**:

```
% of parent (same sample) = 100 × species(sample) / parent(sample)
% of total                = 100 × species(sample) / total(sample)
```

This answers "of the material present now, what fraction is each species." It
is the basis used in oligonucleotide **purity and impurity profiling**, where
each impurity is reported as a percent of the target-plus-impurities total in
one run:

- **Wong & Rye, Agilent BioConfirm Target Plus Impurities** application note.
  "% Quant" = species area / total of target + all impurities per injection
  (their 21-mer target reads 58.49%, impurities a few percent each, summing to
  ~100%).
- **Zhao et al., J. Pharm. Biomed. Anal. 270 (2026)** (inclisiran). "Relative
  abundance (%)" reported per timepoint and summing to ~100% **within each
  strand family** (sense and antisense normalized separately).

It is correct for composition and purity. It is **not** a reliable kinetic
trend, because the denominator (the same-sample parent or total) is itself
changing over time: a species that is constant in absolute terms appears to
rise simply because the parent shrinks.

### 2. Parent-relative, reference timepoint (kinetic)

Every species divided by the intact parent at the **reference timepoint** of
its arm, a fixed denominator (`parent_relative_timecourse`, `R/degradation.R`):

```
% of initial parent = 100 × species(t) / parent(t_ref)
```

On this basis the parent declines from ~100% and degradants rise from ~0% on
one common axis, so each trajectory reflects its own change. This is the basis
for degradation and metabolism **kinetics**:

- **Ye et al., J. Chromatogr. B 1254 (2025)** (fomivirsen metabolite ID).
  Parent degradation is computed as the parent peak area at each timepoint
  relative to the **T0 sample** (17.65 / 32.85 / 37.83% degraded at 24 h in
  rat / human / dog serum). Same-analyte, reference-timepoint denominator.
- **Gorrochategui et al., Trends Anal. Chem. 82 (2016)** (LC-MS metabolomics
  data analysis). Normalizing to the total signal "is not appropriate and can
  cause serious data distortions" when the analytes' own concentrations change
  across the experiment, which is exactly a degradation time course.
  Internal-standard / reference normalization is preferred for quantitative
  cross-sample comparison.

### Pre-dose / negative controls

A leading pre-dose sample with no drug has near-zero parent signal. Using it as
the 100% reference would make every later timepoint read as fully degraded. The
app detects it by **signal, not the clock**: the leading run of timepoints
whose mean parent signal is below 10% of the arm's maximum is flagged as a
control and excluded from the reference. A genuine t=0 **dosing** point with
full parent signal is kept as the reference. Interior low points are not
flagged. A caller- or user-supplied reference timepoint is honored exactly.
(Same rule in `degradation_vs_reference`, `parent_relative_timecourse`, and the
degradation map's control detection.)

## Quantitative caveats

- **Signal is not amount.** `parent(t)/parent(t_ref)` is a clean within-analyte
  ratio: ionization efficiency cancels. A degradant divided by the parent mixes
  two molecules with different response factors, so a degradant value is a
  **response ratio**, not a molar percentage. Zhao et al. apply intensity
  scaling factors of 5–10 between species for this reason.
- **Use internal-standard-normalized signal when available.** Across injections,
  `parent(t)/parent(t_ref)` still drifts with recovery, injection volume and
  matrix. The robust form is
  `(species(t)/IS(t)) / (parent(t_ref)/IS(t_ref))`. Internal standards are the
  recommended control for matrix effects and drift in oligonucleotide
  bioanalysis (**Kumar et al., J. Chromatogr. B 1270 (2026)**), which also
  defines recovery against a post-extraction spike.
- **"100%" is the first measured point, not absolute recovery.** The parent at
  the reference timepoint already includes extraction and matrix losses, so the
  anchor is a relative reference, not percent of the spiked dose.
- **Duplexes.** For an siRNA, compute the parent-relative basis per strand so
  each strand's parent anchors its own family (as Zhao et al. do).
- **True amounts** require a calibration curve and, ideally, a stable-isotope
  internal standard, which the app provides separately in the calibration /
  absolute-quantification path.

## Where each is used in the app

- **Degradation Map** — selectable abundance basis for cell opacity
  (`initial_parent` default, or `sample_parent`). The parent track is always
  shown as % remaining vs the reference timepoint.
- **Degradation Summary** — `% of total` composition by class (within-sample),
  the parent-remaining-vs-reference metric, and the parent-relative
  (% of initial parent) time course for every species.
- **Report (Relative Abundance section)** — both bases documented, with the
  references above.

## Functions

- `parent_relative_timecourse()` / `plot_parent_relative_timecourse()` —
  `R/degradation.R`
- `degradation_vs_reference()` / `plot_degradation_vs_reference()` —
  `R/degradation.R`
- `degradation_summary()` (composition, `% of total`) — `R/degradation.R`
- `degradation_position_model(..., norm = c("sample_parent","initial_parent"))`
  — `R/degradation_map.R`
