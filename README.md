# OligoMet Profiler

A general-purpose R pipeline for generating theoretical metabolite
libraries, charge envelopes, isotope patterns, McLuckey MS/MS fragment ions,
and PRM inclusion lists for any therapeutic oligonucleotide -- DNA, RNA,
2'OMe, 2'F, MOE, cEt, LNA, mixed PO/PS backbones, and custom chemistry via
dictionary overrides -- with optional matching against uploaded MS data.
Ships as an installable R package (`OligoMetProfiler`) with a Shiny
dashboard and a CLI driver on top.

Grounded in SynONIM (Lippens et al., JASMS 2024), the Eluforsen metabolite
profiling study (Kim et al., Mol Ther Nucleic Acids 2019), OligoDistiller
(Liu et al., Anal Chem 2025), and the FMVS automatic metabolite ID method
(Ye et al., J Chromatogr B 2025). Four approved oligonucleotide
therapeutics are bundled as worked examples, one per modality class in
Takakusa et al. (2023) -- **nusinersen** (antisense SSO), **inotersen**
(antisense gapmer), **patisiran** (siRNA in LNP) and **givosiran**
(GalNAc-siRNA). The two antisense examples double as the formula
engine's regression anchors: the pipeline reproduces their published
free-acid molecular formulas exactly (`validate_reference()`). They are
illustrations, not default targets -- run the pipeline on any sequence.
See the [Bibliography](#bibliography) for sources.

**New here?** Start with [QUICKSTART.md](inst/help/QUICKSTART.md) -- it walks
through reading your Certificate of Analysis / synthesis documentation and
constructing a valid input string. The full workflow documentation is the
package vignette, [vignettes/OligoMetProfiler.Rmd](vignettes/OligoMetProfiler.Rmd).

## Installation

As an R package, straight from GitHub -- the dashboard ships with it:

```r
# install.packages(c("remotes", "shiny", "DT", "shinyFiles"))
remotes::install_github("nishi76/OligoMet-Profiler")
library(OligoMetProfiler)   # engine functions
OligoMetProfiler::run_app() # launches the dashboard
```

`run_app()` passes any extra arguments to `shiny::runApp()`, so
`run_app(launch.browser = TRUE, port = 8080)` works as expected.

Or clone the repository to use the Shiny dashboard and CLI drivers:

```r
install.packages("shiny")   # if not already installed
source("install_packages.R") # installs the rest of the dependencies
shiny::runApp(".")          # launches the dashboard
```

Or from the shell, once dependencies are installed:

```bash
Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'
```

No path editing is required in either case. From a clone, `shiny::runApp()`
sets the working directory to this folder for the duration of the app and
every module locates its neighbours relative to that; from an installed
package, the app uses the functions in the package namespace. The app opens with a generic 2'MOE/DNA gapmer
example loaded in the sequence box; paste in your own sequence, or use the
"Load example" dropdown to load any of the four bundled reference drugs.
Terminal conjugates travel with the example -- loading the GalNAc-siRNA
sets its 3' conjugate for you.

The dashboard carries its own documentation: a **Help & guides** panel
on the front page renders the quick start, the sequence guide and the
modification reference in three tabs, so an installed user never has to
go looking for them.

If you would rather not assemble triplet notation by hand, the **Manual
sequence entry** panel in the middle of the dashboard takes the three
lines a chemical analysis file gives you directly -- bases, sugars and
linkages -- and builds the sequence for you:

| Field | Example |
|---|---|
| Bases (5'->3') | `TSASTTTSATAATGSTGG` |
| Sugars | `eeeeeeeeeeeeeeeeee` (or just `e` for all positions) |
| Linkages | `sssssssssssssssss` (or just `s`) |

[SEQUENCE_GUIDE.md](inst/help/SEQUENCE_GUIDE.md) walks a first-time user
through reading those three lines off a chemical analysis file, and
through exporting the same sequence as a BioPharma Finder FASTA.

## What's in here

| Path | Role |
|---|---|
| `DESCRIPTION`, `NAMESPACE` | R package metadata -- makes the whole engine installable via `remotes::install_github()`. |
| `inst/app/app.R` | Shiny dashboard -- sequence input, custom chemistry table, full parameter control, optional MS upload, 4-plot summary, Excel/report/PRM downloads. Ships inside the package; launch it with `OligoMetProfiler::run_app()`. |
| `app.R` | Repository-root launcher, so `shiny::runApp(".")`, `shiny::runGitHub()`, and RStudio's "Run App" button work from a clone. Sources `inst/app/app.R`. |
| `R/run_app.R` | `run_app()` -- launches the bundled dashboard from an installed package. |
| `R/chemistry_dict.R` | Element table, formula arithmetic, and `STANDARD_DICT` -- the default chemistry dictionary (DNA/RNA/2'OMe/2'F/MOE/NMA/UNA/GNA/LNA/ENA/cEt sugars; PO/PS, alkyl phosphonate, PACE, mesyl phosphoramidate and phosphoryl guanidine linkages; GalNAc, lipid and fatty acid conjugates) every module falls back to.  Also holds the four bundled reference drugs and the `validate_reference()` self-test. |
| `R/oligo_io.R` | Sequence parsing (triplet / OligoDistiller / structured / three-line bases-sugars-linkages notation), plus BioPharma Finder FASTA export. |
| `R/metabolites.R` | Theoretical metabolite library generation (truncations, endo fragments). |
| `R/mass_isotope.R` | Mass, charge envelope, isotope pattern, PS oxidation series. |
| `R/fragments.R` | McLuckey MS/MS fragment ions, matching, confirmation scoring. |
| `R/ms_matching.R` | mzML/mzXML/peak-list import and MS1/MS2 matching. |
| `R/build_workbook.R` | 7-sheet Excel workbook export. |
| `R/build_report.R` | Plots and HTML/PDF report export. |
| `R/export_acquisition.R` | Thermo Orbitrap Exploris MS1 inclusion / MS2 PRM target list export, plus a fragment-ion reference table (see "Orbitrap Exploris acquisition method export" below). |
| `R/export_spectral.R` | Theoretical MS1 and MS2 spectral libraries as MGF and MSP (see "Spectral library export" below). |
| `R/progress_utils.R` | Console progress bar with elapsed time and adaptive ETA (see "Console progress reporting" below). |
| `inst/py_decode.py` | Python helper for mzML base64/zlib binary decoding, called via `system2()`. Falls back to a pure-R decoder automatically if python3 is not on PATH. |
| `run_custom_oligo.R` | **Primary CLI entry point.** Template driver for any sequence -- copy it, edit the CONFIG block (sequence, chemistry overrides, parameters), and run. Works with standard chemistry out of the box. |
| `tests/` | Validation scripts (print-based, not testthat -- see note below). |
| `inst/help/QUICKSTART.md` | From CoA/synthesis documentation to a valid input string, plus the fastest path to a first run. |
| `inst/help/SEQUENCE_GUIDE.md` | **Start here if you are new to oligo work.** Reading a chemical analysis file, writing the bases/sugars/linkages lines, and preparing a BioPharma Finder input, with a worked example. |
| `inst/help/MODIFICATIONS.md` | Every sugar, backbone and conjugate chemistry the dictionary knows, with per-position mass differences and how to add one it does not. |
| `inst/help/` | These three guides ship inside the package, so the dashboard's Help panel can render them for installed users too. |
| `vignettes/OligoMetProfiler.Rmd` | Full package vignette: input notations, chemistry dictionary and overrides, library generation, masses/envelopes/isotopes, fragment ions, workbook/report/acquisition exports, MS matching. |

## Bundled reference examples

Four approved oligonucleotide therapeutics ship with the package, one per
modality class in Takakusa et al. (2023) Table 1. They are selectable from
the app's "Load example" dropdown and available programmatically as
constants and through the `REFERENCE_OLIGOS` registry:

| Modality | Drug (brand) | nt | Chemistry | Constant |
|---|---|---:|---|---|
| Antisense (SSO) | nusinersen (SPINRAZA) | 18 | PS, uniform 2'-MOE, 5-methyl pyrimidines | `NUSINERSEN_TRIPLET` |
| Antisense (gapmer) | inotersen (TEGSEDI) | 20 | PS, 5-10-5 2'-MOE/DNA gapmer, 5-methyl-C | `INOTERSEN_TRIPLET` |
| siRNA (LNP) | patisiran (ONPATTRO) | 21 | 2'-OMe/ribose, all-PO backbone, dTdT overhang | `PATISIRAN_SENSE_TRIPLET` |
| siRNA (GalNAc) | givosiran (GIVLAARI) | 21 | 2'-OMe/2'-F, partial PS, 3' trivalent GalNAc | `GIVOSIRAN_SENSE_TRIPLET` |

siRNA entries are the **sense strand** -- this pipeline profiles one strand
at a time, so run each strand of a duplex separately. Givosiran's GalNAc is
carried in the registry's `conj3` field, since triplet notation has no
conjugate field (`GIVOSIRAN_SENSE_SPEC` gives the equivalent structured
form).

```r
library(OligoMetProfiler)
spec <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(spec, opts = list(oligo_name = "inotersen",
                                               endo = TRUE, endo_sites = "gap"))
```

Each modality stresses a different part of the pipeline: gapmers cleave
endonucleolytically in the DNA gap, sugar-modified SSOs degrade from the
termini by exonuclease, the all-PO siRNA has no PS oxidation series to
model, and the GalNAc conjugate is retained on 5' truncations but lost from
3' truncations that cut past it.

## Running the pipeline on your own sequence

Copy `run_custom_oligo.R`, edit the `CONFIG` block at the top (sequence,
optional custom chemistry overrides, parameters), and run:

```bash
Rscript run_custom_oligo.R
```

Or use the Shiny app (`run_app()`) for an interactive session. The app's "Save
to folder (optional)" field writes the workbook, report, PRM list, and
acquisition method lists directly to a folder you choose, on top of the
usual download buttons -- useful when running locally in RStudio/Rscript
rather than through a browser-hosted deployment. "Browse..." opens an
in-app folder picker (via the `shinyFiles` package) that browses the
filesystem of the machine R is running on -- your own machine, in the
local setup above.

Or drive the engine directly from R after installing the package:

```r
library(OligoMetProfiler)
spec <- parse_input("Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm")
mets <- generate_metabolites(spec, opts = list(max_3p = 10, max_5p = 10, endo = TRUE))
build_workbook(spec, mets, opts = list(z_range = 3:12, max_oxid = 6),
               output_dir = ".")
```

## Validating the formula engine

The chemistry engine is anchored to two approved drugs whose free-acid
molecular formulas appear on their product labelling. It reproduces both
exactly (0 ppm):

| Drug | Chemistry exercised | Published formula | Avg MW |
|---|---|---|---|
| nusinersen | uniform 2'-MOE, full PS, 5-methyl pyrimidines | C234H340N61O128P17S17 | 7127.2 |
| inotersen | 5-10-5 MOE/DNA gapmer, full PS, 5-methyl-C | C230H318N69O121P19S19 | 7183.08 |

To re-run that regression check:

```r
library(OligoMetProfiler)   # or source the R/ modules from a clone
validate_reference()
```

And from the shell, the validation scripts:

```bash
Rscript tests/test_metabolites.R
Rscript tests/test_mass_isotope.R
Rscript tests/test_fragments.R
Rscript tests/test_ms_matching.R
Rscript tests/test_outputs.R
```

Most of these print computed values to the console rather than asserting
pass/fail -- read the output rather than the exit code. `test_mass_isotope.R`
is the exception: it stops with an error if the formula engine no longer
reproduces the published molecular formulas, or if the charge envelope
stops being self-consistent across charge states.

## BioPharma Finder compatibility

Triplet sequences copied from Thermo BioPharma Finder's sequence editor
(e.g. `Ad-pTd-pCd-pAd`) are accepted directly -- `p` is recognized as the
phosphodiester linkage, matching BPF's own notation, alongside this
pipeline's original `o`. `s` (phosphorothioate) already matched BPF's
convention. BPF's 5'/3' terminal modification codes (biotin, cAG/cAU/ARCA/
mCAP cap analogs, triantennary GalNAc) are available from the conjugate
dropdowns in the app (or `conj5`/`conj3` in the CLI drivers) under
descriptive names rather than BPF's own single-letter codes, since those
letters (`a`, `u`, `r`, `c`) already mean specific sugars or linkages in
this dictionary. Their formulas come straight from the BioPharma Finder
5.2 Oligonucleotide Analysis User Guide's "Modification notation" topic
and are flagged for verification (`verify = TRUE` in `chemistry_dict.R`)
until checked against a real BPF-computed mass.

Going the other way, `format_biopharma_fasta()` -- the app's **BPF
FASTA** button -- writes the current sequence as a FASTA record for
BioPharma Finder's Sequence Manager ("Import FASTA File"):

```
>nusinersen | 18-mer
Te-sSe-sAe-sSe-sTe-sTe-sTe-sSe-sAe-sTe-sAe-sAe-sTe-sGe-sSe-sTe-sGe-sGe
```

Two things to check on the BPF side, since they depend on your
installation rather than on this package: **sugar codes** resolve against
BPF's own Building Block editor, so a modification your site defined
under a different code needs renaming or adding there; and **terminal
conjugates** have nowhere to live in triplet notation, so the export
writes them into the header as a reminder and you set them as terminus
modifications yourself. See
[SEQUENCE_GUIDE.md](inst/help/SEQUENCE_GUIDE.md) section 7.

## Modification coverage

The dictionary covers the chemistries used in approved and clinical
oligonucleotide therapeutics, plus the next-generation ones now in
development:

| Class | Codes |
|---|---|
| Sugars | `d` `r` `m` `f` `FANA` `e`/`MOE` **`NMA`** `allyl` `AP` `NH2` `UNA` `GNA` `LNA` `ENA` `cEt` `NP` |
| Backbones | `o`/`p` `s`/`u` `mp` `prp` `ibu` `chx` `mop` `pace` `tpace` **`msp`** `pgo` `tmg` |
| Conjugates | GalNAc (mono/trivalent/triantennary), cholesterol, C6/C12, TEG, FAM, Cy3, phosphate and thiophosphate caps, cap analogs, fatty acids (`myristoyl` `palmitoyl` `stearoyl` `docosanoyl`) |

A few worth calling out, because their mass differences are the ones
people most often need:

- **2'-O-NMA** (2'-O-[2-(methylamino)-2-oxoethyl]) is **+12.9953 Da per
  position vs MOE** -- an amide where MOE has an ether.
- **Mesyl phosphoramidate** (msPA) is **+61.0164 Da per bond vs PS**.
- **Phosphoryl guanidine** comes in two variants that differ by
  2.0157 Da per bond (`pgo` cyclic, `tmg` tetramethyl) -- check which
  one your synthesis used before assigning anything.
- **FANA** is exactly isobaric with 2'-F ribo; MS cannot distinguish
  them, and the separate code exists only to keep the sequence record
  correct.
- **UNA** and **GNA** are duplex destabilizers used positionally in
  siRNA, at +2.0157 and −58.0055 Da vs RNA.

[inst/help/MODIFICATIONS.md](inst/help/MODIFICATIONS.md) is the full
reference: every code, its residue formula, its mass difference from the
parent chemistry, how to enter it, and how to add one the dictionary
does not have. `tests/test_modifications.R` recomputes every number in
it and anchors each sugar against an independently known nucleoside
formula.

Morpholino (PMO), PNA and boranophosphate are deliberately absent --
the first two do not decompose into base/sugar/linkage residues, and
boranophosphate needs boron added to the element table (with its
isotope abundances, not just a mass). The guide explains what to do
about each.

## Orbitrap Exploris acquisition method export

The app and `run_custom_oligo.R` both export targeted mass lists formatted
for the Orbitrap Exploris Method Editor's Targeted Mass filter (Mass List
Type "m/z & z", Time Mode "Start/End Time" -- select that when importing,
or "Unscheduled" if you'd rather ignore the time columns). Column layout
follows the Orbitrap Exploris 120 Software Manual's "Targeted Inclusion --
Targeted Mass filter" topic exactly: `Compound, Formula, Adduct, m/z, z,
Intensity Threshold, t start (min), t stop (min), HCD Collision Energies
(%), Maximum Injection Time (ms)`.

- **MS1 Inclusion List** -- every (metabolite, PS-oxidation level, charge
  state) combination across the configured charge envelope. Import into a
  Full Scan experiment's Targeted Mass filter to prioritize these masses.
- **MS2 PRM Target List** -- the same style of list, narrowed to a smaller
  charge range (PRM duty cycle degrades quickly with target count) and
  tagged with an HCD collision energy. Import as the target list for a
  Targeted MS2 (tMS2)/PRM scan.
- **MS2 Fragment Reference** -- theoretical McLuckey fragment ions (a/a-B/
  b/b-B/c/w/x/y, plus internal fragments) per metabolite. This is **not**
  a Method Editor import -- PRM targets a precursor and records the full
  fragment spectrum, it doesn't take individual fragment ions as
  acquisition input. This table is for interpreting the resulting spectra
  afterward (manual annotation, or building a Skyline-style transition
  list).

The default HCD NCE (20%) is a starting point for HCD of phosphorothioate
oligonucleotide backbones, not a validated instrument parameter -- optimize
it per method. Both list functions cap the number of rows (`max_targets`,
configurable) since the Method Editor's own limit is 150,000 rows per file
but PRM cycle time degrades well before that.

## Spectral library export

The theoretical library also exports as MGF and MSP -- the two formats
that spectral-library tools read (MS-DIAL, mzVault/Compound Discoverer,
MZmine, matchms, SIRIUS, GNPS). Four files, from the app's **Spectral
Libraries** download row, from the "Save to folder" field, or from
`export_spectral_libraries()`:

```r
export_spectral_libraries(mets, dict, out_dir = "results", prefix = "my_oligo")
#> results/my_oligo_MS1_library.mgf   my_oligo_MS1_library.msp
#> results/my_oligo_MS2_library.mgf   my_oligo_MS2_library.msp
```

- **MS1 library** -- one spectrum per (metabolite, PS-oxidation level,
  charge state). The peaks are the theoretical isotope cluster at that
  charge, so their intensities are real relative abundances and can be
  matched against an acquired isotope pattern. Peaks are annotated with
  their true isotopologue offset (`M+0`, `M+1`, ...), and the cluster is
  anchored on the m/z computed directly from the formula.
- **MS2 library** -- one spectrum per (metabolite, precursor charge
  state), holding the McLuckey fragment ions (terminal a/a-B/b/b-B/c/w/x/y
  plus w-a/w-b/w-d internal ions) at the requested fragment charges, each
  annotated in the MSP (`w4^2-`, `a-B5^1-`, `w-a(3,8)^1-`).

**One caveat, and it matters.** This pipeline has no fragment-intensity
model, so every MS2 peak is written at a flat intensity of 100. Match on
*m/z*; do not run intensity-weighted dot-product or cosine scoring
against the MS2 libraries, because with flat intensities that score
reduces to a function of peak count. The MS1 libraries carry genuine
isotope abundances and are not subject to this caveat.

Both builders cap their output (`max_spectra`) since a full metabolite
library across a wide charge envelope runs to thousands of spectra: at
default settings a 51-metabolite library produces ~3,400 MS1 spectra and
~200 MS2 spectra (about 20 seconds each to build, a few MB on disk).

## Console progress reporting

Every entry point (`run_custom_oligo.R` and the Shiny app's R console
output) prints a weighted, multi-step progress bar as it runs:

```
[#########---------------] 38%  Step 5/11: Generating McLuckey MS/MS fragment ions
  elapsed 12s | ETA (remaining) 19s
```

Plus a live in-place updating sub-bar during the Charge Envelopes sheet
specifically, since that step dominates runtime:

```
  [############----] 75%  Charge envelopes (44/58)  step elapsed 31s  step ETA 10s
```

The ETA is adaptive, not a fixed prediction -- each step is given a rough
starting weight (`R/progress_utils.R`; the workbook step is weighted
heaviest since it's the actual bottleneck), and after every completed step
the tracker recomputes "time per unit of weight" from what's actually
elapsed so far and applies that to the remaining weight. So the estimate
tightens up as the run progresses rather than trusting the initial guesses
for the whole run. A final line reports total elapsed time once the
pipeline finishes.

## Performance

The **Charge Envelopes** computation is the heaviest step: one isotope
pattern per metabolite × PS-oxidation level. Isotope patterns are memoized
per formula (both the enviPat engine and the built-in fallback), so charge
states and adducts add no repeated cost, and the Excel workbook itself is
written in bulk -- one write per sheet, styling restricted to headers and a
few accent cells -- so `saveWorkbook()` completes in seconds even for
libraries with tens of thousands of envelope rows. A full default run on a
16-mer (endonuclease fragments on, `max_oxid = 6`, z = 3–12) builds and
saves the workbook in well under a minute on a typical laptop.

To speed a long run up further: uncheck "Include endonuclease fragments",
lower "Max PS oxid.", or narrow the charge range -- all directly reduce the
number of isotope-pattern calls.

If the app appears completely unresponsive (Run does nothing, no error, no
progress) after clicking a folder picker button, an earlier version of this
app used `tcltk`/`utils::choose.dir()` for the "Browse..." button --
blocking native OS dialogs, which on a single-threaded Shiny process can
freeze the entire app (not just the picker) if the dialog fails to
initialize properly, which happened intermittently on Windows depending on
how the R process was spawned. This was replaced with an in-app
`shinyFiles` picker, which can't block the R process the same way. If
you're on a version with the old picker, update to this one. Separately,
the Run handler is wrapped in a top-level `tryCatch`, so any unexpected
error anywhere in the pipeline gets reported in the status panel instead
of silently doing nothing.

## Dependencies

Required: `openxlsx`, `ggplot2`, `xml2`, `xfun` (installed automatically
with the package). For the Shiny app: `shiny`, `DT`, `bslib`,
`shinyFiles`. Optional: `enviPat` (higher-accuracy isotope patterns; a
built-in convolution method is used if absent), `rmarkdown` (needed to
render HTML/PDF reports). `install_packages.R` installs the required set
and offers to install the optional set.

mzML import additionally uses `python3` if it's on PATH, and vendor raw
conversion uses `msconvert` (ProteoWizard) if it's on PATH. Both are
optional -- the app runs without them, with reduced MS-import coverage.

## Known limitations

- No `testthat` suite -- the `tests/` scripts are console-output checks,
  not asserted unit tests.
- Default parameters (charge range, oxidation cap, ppm tolerance, etc.) are
  set independently in `inst/app/app.R` and in the CLI driver, so changing a
  default in one place doesn't propagate to the other.

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
   3:60-70 (1992). -- fragment ion nomenclature (a/a-B/b/c/w/x/y/z).

### Reference drugs and modality classification

6. Takakusa H., Iwazaki N., Nishikawa M., Yoshida T., Obika S., Inoue T.
   Drug Metabolism and Pharmacokinetics of Antisense Oligonucleotide
   Therapeutics: Typical Profiles, Evaluation Approaches, and Points to
   Consider Compared with Small Molecule Drugs. *Nucleic Acid Therapeutics*
   33(2):83-94 (2023). [doi:10.1089/nat.2022.0054](https://doi.org/10.1089/nat.2022.0054)
   -- source of the modality classification (Table 1) used for the four
   bundled examples, and of the class-specific metabolism routes.
7. Egli M., Manoharan M. Chemistry, structure and function of approved
   oligonucleotide therapeutics. *Nucleic Acids Research* 51(6):2529-2573
   (2023). [doi:10.1093/nar/gkad067](https://doi.org/10.1093/nar/gkad067)
   -- chemistry and structural review of the approved oligonucleotide drugs.
8. SPINRAZA (nusinersen) prescribing information, Biogen -- 18-mer sequence,
   uniform 2'-MOE/PS chemistry with 5-methyl pyrimidines, and free-acid
   molecular formula C234H340N61O128P17S17.
9. TEGSEDI (inotersen) prescribing information, Akcea Therapeutics / Ionis
   Pharmaceuticals -- 5-10-5 MOE gapmer design, TTR 3'UTR target sequence
   (bases 618-637), and free-acid molecular formula C230H318N69O121P19S19
   (average MW 7183.08).
10. ONPATTRO (patisiran) prescribing information, Alnylam Pharmaceuticals --
    siRNA duplex sequences, per-position 2'-OMe pattern, dTdT overhangs, and
    lipid-nanoparticle formulation.
11. GIVLAARI (givosiran) prescribing information, Alnylam Pharmaceuticals --
    ESC-GalNAc siRNA sequences, 2'-F/2'-OMe alternation, terminal
    phosphorothioate placement, and the trivalent GalNAc (L96) 3' conjugate.

### Instrument and software documentation

12. Thermo Fisher Scientific. *BioPharma Finder 5.2 Oligonucleotide Analysis
    User Guide* -- triplet sequence notation and terminal modification codes.
13. Thermo Fisher Scientific. *Orbitrap Exploris 120 Software Manual*,
    "Targeted Inclusion -- Targeted Mass filter" -- the mass-list column
    layout used by the acquisition exports.

## License

MIT -- see [LICENSE.md](LICENSE.md).
