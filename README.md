# OligoMet Profiler

> **FOR RESEARCH USE ONLY.** Not for diagnostic, clinical, or regulatory
> submission use. Every value is a computed prediction, not a
> measurement, and must be confirmed experimentally. No warranty; see
> [DISCLAIMER.md](DISCLAIMER.md).

An R pipeline that generates theoretical metabolite libraries, charge
envelopes, isotope patterns, McLuckey MS/MS fragment ions, and PRM
inclusion lists for any therapeutic oligonucleotide — DNA, RNA, 2'-OMe,
2'-F, MOE, cEt, LNA, mixed PO/PS backbones, and custom chemistry via
dictionary overrides — with optional matching against uploaded MS data.
Ships as an installable package (`OligoMetProfiler`) with a Shiny
dashboard.

Method background comes from SynONIM (Lippens et al. 2024), the
eluforsen metabolite profiling study (Kim et al. 2019), OligoDistiller
(Liu et al. 2025), and FMVS automatic metabolite ID (Ye et al. 2025);
see the [Bibliography](#bibliography). Four approved drugs are bundled
as worked examples, one per modality class in Takakusa et al. (2023) —
they are illustrations, not default targets.

**New here?** Start with [QUICKSTART.md](inst/help/QUICKSTART.md), which
covers reading a Certificate of Analysis and building a valid input
string. Full workflow documentation is in the
[vignette](vignettes/OligoMetProfiler.Rmd).

## Installation

```r
# install.packages(c("remotes", "shiny", "DT", "shinyFiles"))
remotes::install_github("nishi76/OligoMet-Profiler")
library(OligoMetProfiler)   # engine functions
OligoMetProfiler::run_app() # launches the dashboard
```

Or from a clone, for the dashboard and CLI drivers:

```r
install.packages("shiny")
source("install_packages.R")  # installs the rest
shiny::runApp(".")
```

The dashboard carries its own documentation in a **Help & guides** panel
(quick start, sequence guide, modification reference).

Rather than assembling triplet notation by hand, the **Manual sequence
entry** panel takes the three lines off your analysis document — bases,
sugars, linkages — and builds the sequence:

| Field | Example |
|---|---|
| Bases (5'→3') | `TSASTTTSATAATGSTGG` |
| Sugars | `eeeeeeeeeeeeeeeeee` (or just `e` for all positions) |
| Linkages | `sssssssssssssssss` (or just `s`) |

[SEQUENCE_GUIDE.md](inst/help/SEQUENCE_GUIDE.md) walks a first-time user
through reading those three lines and exporting a BioPharma Finder FASTA.

## Bundled reference examples

Selectable from the app's "Load example" dropdown, and available
programmatically as constants and via the `REFERENCE_OLIGOS` registry:

| Modality | Drug (brand) | nt | Chemistry | Constant |
|---|---|---:|---|---|
| Antisense (SSO) | nusinersen (SPINRAZA) | 18 | PS, uniform 2'-MOE, 5-methyl pyrimidines | `NUSINERSEN_TRIPLET` |
| Antisense (gapmer) | inotersen (TEGSEDI) | 20 | PS, 5-10-5 2'-MOE/DNA gapmer, 5-methyl-C | `INOTERSEN_TRIPLET` |
| siRNA (LNP) | patisiran (ONPATTRO) | 21 | 2'-OMe/ribose, all-PO backbone, dTdT overhang | `PATISIRAN_SENSE_TRIPLET` |
| siRNA (GalNAc) | givosiran (GIVLAARI) | 21 | 2'-OMe/2'-F, partial PS, 3' trivalent GalNAc | `GIVOSIRAN_SENSE_TRIPLET` |

```r
spec <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(spec, opts = list(oligo_name = "inotersen",
                                               endo = TRUE, endo_sites = "gap"))
```

siRNA entries are the **sense strand** — the pipeline profiles one
strand at a time, so run each strand separately. Givosiran's GalNAc is
carried in the registry's `conj3` field, since triplet notation has no
conjugate field.

Each modality stresses a different part of the pipeline: gapmers cleave
endonucleolytically in the DNA gap, sugar-modified SSOs degrade from the
termini, the all-PO siRNA has no PS oxidation series, and the GalNAc
conjugate is retained on 5' truncations but lost from 3' ones cutting
past it.

## BioPharma Finder compatibility

Triplet sequences from BPF's sequence editor (`Ad-pTd-pCd-pAd`) parse
directly — `p` is recognized as phosphodiester alongside this pipeline's
`o`, and `s` already matched. BPF's terminal modification codes (biotin,
cap analogs, triantennary GalNAc) appear in the app's conjugate
dropdowns under descriptive names, since BPF's single letters (`a`, `u`,
`r`, `c`) already mean specific sugars or linkages here. Formulas come
from the BioPharma Finder 5.2 user guide.

Going the other way, `format_biopharma_fasta()` — the app's **BPF
FASTA** button — writes the sequence as a FASTA record for BPF's
Sequence Manager:

```
>nusinersen | 18-mer
Te-sSe-sAe-sSe-sTe-sTe-sTe-sSe-sAe-sTe-sAe-sAe-sTe-sGe-sSe-sTe-sGe-sGe
```

Two things depend on your BPF installation: **sugar codes** resolve
against its Building Block editor, so a locally-defined code needs
renaming or adding there; and **terminal conjugates** have nowhere to
live in triplet notation, so the export writes them into the header as a
reminder and you set them yourself. See
[SEQUENCE_GUIDE.md](inst/help/SEQUENCE_GUIDE.md) section 7.

## Modification coverage

| Class | Codes |
|---|---|
| Sugars | `d` `r` `m` `f` `FANA` `e`/`MOE` **`NMA`** `allyl` `AP` `NH2` `UNA` `GNA` `LNA` `ENA` `cEt` `NP` |
| Backbones | `o`/`p` `s`/`u` `mp` `prp` `ibu` `chx` `mop` `pace` `tpace` **`msp`** `pgo` `tmg` |
| Conjugates | GalNAc (mono/tri/triantennary), cholesterol, C6/C12, TEG, FAM, Cy3, phosphate and thiophosphate caps, cap analogs, fatty acids |

The mass differences people most often need:

- **2'-O-NMA** is **+12.9953 Da per position vs MOE** — an amide where
  MOE has an ether.
- **Mesyl phosphoramidate** (msPA) is **+61.0164 Da per bond vs PS**.
- **Phosphoryl guanidine** has two variants differing by 2.0157 Da per
  bond (`pgo` cyclic, `tmg` tetramethyl) — check which yours is.
- **FANA** is exactly isobaric with 2'-F ribo; MS cannot distinguish
  them.
- **UNA** and **GNA** are siRNA destabilizers at +2.0157 and −58.0055 Da
  vs RNA.

[MODIFICATIONS.md](inst/help/MODIFICATIONS.md) is the full reference:
every code, residue formula, mass difference, and how to add one the
dictionary lacks. `tests/test_modifications.R` recomputes every number.

Morpholino (PMO), PNA, and boranophosphate are deliberately absent — the
first two do not decompose into base/sugar/linkage residues, and
boranophosphate needs boron added to the element table with its isotope
abundances. The guide explains what to do about each.

## Orbitrap Exploris acquisition export

Targeted mass lists follow the Orbitrap Exploris Software Manual's
"Targeted Inclusion — Targeted Mass filter" column layout exactly
(`Compound, Formula, Adduct, m/z, z, Intensity Threshold, t start, t
stop, HCD Collision Energies (%), Maximum Injection Time`). Import with
Mass List Type "m/z & z" and Time Mode "Start/End Time" (or
"Unscheduled").

- **MS1 Inclusion List** — every (metabolite, PS-oxidation level, charge
  state) across the configured envelope, for a Full Scan experiment.
- **MS2 PRM Target List** — the same, narrowed to a smaller charge range
  (PRM duty cycle degrades with target count) and tagged with an HCD
  collision energy, for a tMS2/PRM scan.
- **MS2 Fragment Reference** — theoretical McLuckey fragment ions
  (a/a-B/b/b-B/c/w/x/y plus internal). **Not** a Method Editor import;
  PRM targets precursors, so this is for interpreting the resulting
  spectra afterward or building a Skyline-style transition list.

The default HCD NCE (20%) is a starting point for PS backbones, not a
validated parameter — optimize per method. Both list functions cap rows
(`max_targets`), since PRM cycle time degrades well before the Method
Editor's own 150,000-row limit.

## Spectral library export

The library also exports as MGF and MSP, read by MS-DIAL,
mzVault/Compound Discoverer, MZmine, matchms, SIRIUS, and GNPS — from
the app's **Spectral Libraries** download row or from R:

```r
export_spectral_libraries(mets, dict, out_dir = "results", prefix = "my_oligo")
#> results/my_oligo_MS1_library.mgf   my_oligo_MS1_library.msp
#> results/my_oligo_MS2_library.mgf   my_oligo_MS2_library.msp
```

- **MS1 library** — one spectrum per (metabolite, oxidation level,
  charge state). Peaks are the theoretical isotope cluster, so
  intensities are real relative abundances matchable against an acquired
  pattern, annotated with isotopologue offset (`M+0`, `M+1`, …).
- **MS2 library** — one spectrum per (metabolite, precursor charge),
  holding McLuckey terminal and internal ions at the requested fragment
  charges, annotated in the MSP (`w4^2-`, `a-B5^1-`, `w-a(3,8)^1-`).
  Intensities are flat placeholders — match on *m/z* only.

## Performance

**Charge Envelopes** is the heaviest step: one isotope pattern per
metabolite × oxidation level. Patterns are memoized per formula, so
charge states and adducts add no repeated cost, and the workbook is
written in bulk. A full default run on a 16-mer (endonuclease fragments
on, `max_oxid = 6`, z = 3–12) builds and saves in well under a minute on
a typical laptop.

To speed up a long run: uncheck "Include endonuclease fragments", lower
"Max PS oxid.", or narrow the charge range — each directly reduces
isotope-pattern calls.

If the app seems unresponsive after clicking the folder picker, you are
on an old version: `tcltk`/`utils::choose.dir()` used blocking native
dialogs that could freeze the single-threaded Shiny process. That was
replaced with an in-app `shinyFiles` picker. Separately, the Run handler
is wrapped in a top-level `tryCatch`, so unexpected errors are reported
in the status panel instead of silently doing nothing.

## Batch processing and statistical comparison

**Try it now, no data of your own needed:** a bundled example
(`inst/extdata/batch_example/`) ships 6 synthetic mzML files — 3
"control" + 3 "treated" replicates of the inotersen reference sequence,
including a contaminant trace that matches nothing (exercises the
retained "unidentified peaks" path) and a real MS2 spectrum per file for
confirmation. Run it end to end with:

```bash
Rscript inst/examples/run_batch_example.R
```

which matches, confirms, and statistically compares (Welch t-test) all
six files in about a minute, writing a workbook, an unmatched-peaks CSV,
and a statistics CSV to `results_batch_example/`. `inst/extdata/batch_example/generate_example.py`
shows exactly how the synthetic data was built, with fixed seeds, if you
want to regenerate or adapt it.

For many raw files at once — a multi-sample experiment rather than one
sequence against one file — a parallel Python pipeline
(`inst/python/oligomet_deconv/`) streams each mzML/mzXML file, detects
chromatographic features, and groups them into charge-state envelopes by
neutral-mass agreement across observed charge states (the same
`M = z*(mz - proton)` relationship the R side already uses, not a
peptide-style isotope/averagine model — oligonucleotides don't fit that).
One worker process per file, so batches scale with core count; a bad file
is isolated and logged rather than aborting the run.

R matches the resulting features against the theoretical library
(reusing `match_ms1()` unchanged), retains anything that didn't match as
its own "unidentified peaks" table for QC, and — when MS2 confirmation is
enabled — confirms matched hits against the theoretical fragment library
using a *targeted* precursor watch-list, so only scans near a real
candidate are captured (cheap on top of the deconvolution pass, and
reuses `confirm_metabolite()` unchanged). A small statistics suite then
compares each confirmed metabolite's abundance across 2+ experimental
groups (Welch t-test or one-way ANOVA + Tukey HSD, both BH-adjusted) or
across a time series (linear trend).

```r
source("R/batch_ms_processing.R"); source("R/statistics.R")

deconv <- run_batch_deconvolution(list.files("raw_data", pattern = "\\.mzML$", full.names = TRUE))
features <- read_batch_features(deconv$features_path)
batch <- annotate_metabolites_batch(mets, features, dict = dict)

abund <- build_abundance_matrix(batch$ms1_matches)
long  <- abundance_long(abund, data.frame(sample = c("ctrl_1","ctrl_2","treat_1","treat_2"),
                                           group = c("control","control","treated","treated")))
compare_two_groups(long, "control", "treated")
```

Or from the command line: copy `run_batch_ms.R`, edit its `CONFIG` block
(file glob, sample metadata, analysis mode), and run
`Rscript run_batch_ms.R`. The Shiny app has the same workflow under
**Batch MS Processing (optional)** — upload multiple files, fill in a
Group or Timepoint column in the sample table, and the **Batch Results**,
**Unidentified Peaks**, and **Statistics** tabs populate after Run.

Requires Python 3.9+ with `inst/python/requirements.txt` installed
(`pip install -r inst/python/requirements.txt`); everything else in the
package works without it.

## Dependencies

Required: `openxlsx`, `ggplot2`, `xml2`, `xfun` (installed with the
package). Dashboard: `shiny`, `DT`, `bslib`, `shinyFiles`. Optional:
`enviPat` (higher-accuracy isotope patterns; a built-in convolution is
used if absent) and `rmarkdown` (HTML/PDF reports). `install_packages.R`
installs the required set and offers the optional one.

mzML import is native R (no external dependency); vendor raw conversion
uses ProteoWizard `msconvert` if on PATH, optional — the app runs
without it, with reduced MS-import coverage. Batch/parallel processing of
many raw files (see above) is a separate, optional Python component —
`pip install -r inst/python/requirements.txt` — needed only for that
workflow.

## Known limitations

- No `testthat` suite — the `tests/` scripts are console-output checks,
  not asserted unit tests.
- Default parameters (charge range, oxidation cap, ppm tolerance) are
  set independently in `inst/app/app.R` and the CLI driver, so a change
  in one does not propagate to the other.

## Bibliography

### Methodology

1. Lippens J.L. *et al.* SynONIM: a Synthetic Oligonucleotide Nomenclature
   for Impurities and Modifications. *J Am Soc Mass Spectrom* (2024).
2. Kim J. *et al.* Metabolite Profiling of the Antisense Oligonucleotide
   Eluforsen Using Liquid Chromatography-Mass Spectrometry. *Mol Ther
   Nucleic Acids* (2019).
3. Liu R. *et al.* OligoDistiller: untargeted profiling of therapeutic
   oligonucleotide drugs and metabolites. *Anal Chem* (2025).
4. Ye X. *et al.* Automatic metabolite identification of antisense
   oligonucleotides by full MS variable scanning. *J Chromatogr B* (2025).
5. McLuckey S.A., Van Berkel G.J., Glish G.L. Tandem mass spectrometry of
   small, multiply charged oligonucleotides. *J Am Soc Mass Spectrom*
   3:60-70 (1992). — fragment ion nomenclature (a/a-B/b/c/w/x/y/z).

### Reference drugs and modality classification

6. Takakusa H. *et al.* Drug Metabolism and Pharmacokinetics of Antisense
   Oligonucleotide Therapeutics. *Nucleic Acid Therapeutics* 33(2):83-94
   (2023). [doi:10.1089/nat.2022.0054](https://doi.org/10.1089/nat.2022.0054)
   — modality classification (Table 1) for the four bundled examples, and
   the class-specific metabolism routes.
7. Egli M., Manoharan M. Chemistry, structure and function of approved
   oligonucleotide therapeutics. *Nucleic Acids Research* 51(6):2529-2573
   (2023). [doi:10.1093/nar/gkad067](https://doi.org/10.1093/nar/gkad067)
8. SPINRAZA (nusinersen) prescribing information, Biogen — 18-mer sequence,
   uniform 2'-MOE/PS chemistry, free-acid formula C234H340N61O128P17S17.
9. TEGSEDI (inotersen) prescribing information, Akcea/Ionis — 5-10-5 MOE
   gapmer, TTR 3'UTR target (bases 618-637), free-acid formula
   C230H318N69O121P19S19 (average MW 7183.08).
10. ONPATTRO (patisiran) prescribing information, Alnylam — duplex
    sequences, 2'-OMe pattern, dTdT overhangs, LNP formulation.
11. GIVLAARI (givosiran) prescribing information, Alnylam — ESC-GalNAc
    sequences, 2'-F/2'-OMe alternation, terminal PS placement, trivalent
    GalNAc (L96) 3' conjugate.

### Instrument and software documentation

12. Thermo Fisher Scientific. *BioPharma Finder 5.2 Oligonucleotide Analysis
    User Guide* — triplet notation and terminal modification codes.
13. Thermo Fisher Scientific. *Orbitrap Exploris 120 Software Manual*,
    "Targeted Inclusion — Targeted Mass filter" — mass-list column layout.

## Author

**Nishikant Wase, PhD** — author and developer
Research Scientist, Thermo Fisher Scientific
<nishikant.wase@gmail.com>

If this software contributes to work you publish, please cite it as:

> Wase, N. *OligoMetProfiler: theoretical metabolite libraries, charge
> envelopes, and MS/MS fragment ions for therapeutic oligonucleotides.*
> R package. <https://github.com/nishi76/OligoMet-Profiler>

## Disclosure and disclaimer

**Research use only.** Not a medical device; not validated for
diagnostic use, clinical decision-making, quality-control release
testing, or regulatory submission. **Everything it reports is a
prediction** computed from a chemistry dictionary — confirm every
assignment experimentally. **No warranty, no liability:** provided "as
is" under the MIT licence; the author is not liable for any damages,
loss of data, wasted instrument time, or erroneous conclusions, and the
user is solely responsible for validating every result.

**Affiliation.** The author is a Research Scientist at Thermo Fisher
Scientific, disclosed because the package interoperates with Thermo
Fisher products. OligoMetProfiler is an independent personal project:
not a Thermo Fisher Scientific product, and not supplied, reviewed,
endorsed, or approved by Thermo Fisher Scientific or any other company.
All interoperability uses publicly documented formats only. Product and
drug names are used for identification only and remain their owners'
trademarks; the bundled reference drugs are worked examples from
published literature and public filings. All views and outputs are the
author's own.

Full statement: [DISCLAIMER.md](DISCLAIMER.md), or `oligomet_about()` in R.

## License

MIT — see [LICENSE.md](LICENSE.md). Copyright (c) 2026 Nishikant Wase.
