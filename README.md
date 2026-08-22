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
(Ye et al., J Chromatogr B 2025). The published ION337 gapmer is bundled as
a known-mass reference for the formula-engine self-test
(`validate_ION337()`, cross-checked against a reference workbook in
`tests/`) -- it is not a default target or an assumed input; run the
pipeline on any sequence.

**New here?** Start with [QUICKSTART.md](QUICKSTART.md) -- it walks
through reading your Certificate of Analysis / synthesis documentation and
constructing a valid input string. The full workflow documentation is the
package vignette, [vignettes/OligoMetProfiler.Rmd](vignettes/OligoMetProfiler.Rmd).

## Installation

As an R package, straight from GitHub:

```r
# install.packages("remotes")
remotes::install_github("nishi76/OligoMet-Profiler")
library(OligoMetProfiler)
```

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

No path editing is required -- `shiny::runApp()` sets the working directory
to this folder for the duration of the app, and every module locates its
neighbours relative to that. The app opens with a generic 2'MOE/DNA gapmer
example loaded in the sequence box; paste in your own sequence, or use the
"Load example" dropdown to switch between that generic example and the
ION337 reference case.

## What's in here

| Path | Role |
|---|---|
| `DESCRIPTION`, `NAMESPACE` | R package metadata -- makes the whole engine installable via `remotes::install_github()`. |
| `app.R` | Shiny dashboard -- sequence input, custom chemistry table, full parameter control, optional MS upload, 4-plot summary, Excel/report/PRM downloads. |
| `R/chemistry_dict.R` | Element table, formula arithmetic, and `STANDARD_DICT` -- the default chemistry dictionary (DNA/RNA/2'OMe/2'F/MOE/cEt/LNA sugars, PO/PS linkages, common conjugates) every module falls back to. Also holds the ION337 reference example and `validate_ION337()` self-test. |
| `R/oligo_io.R` | Sequence parsing (triplet / OligoDistiller / structured notation). |
| `R/metabolites.R` | Theoretical metabolite library generation (truncations, endo fragments). |
| `R/mass_isotope.R` | Mass, charge envelope, isotope pattern, PS oxidation series. |
| `R/fragments.R` | McLuckey MS/MS fragment ions, matching, confirmation scoring. |
| `R/ms_matching.R` | mzML/mzXML/peak-list import and MS1/MS2 matching. |
| `R/build_workbook.R` | 7-sheet Excel workbook export. |
| `R/build_report.R` | Plots and HTML/PDF report export. |
| `R/export_acquisition.R` | Thermo Orbitrap Exploris MS1 inclusion / MS2 PRM target list export, plus a fragment-ion reference table (see "Orbitrap Exploris acquisition method export" below). |
| `R/progress_utils.R` | Console progress bar with elapsed time and adaptive ETA (see "Console progress reporting" below). |
| `inst/py_decode.py` | Python helper for mzML base64/zlib binary decoding, called via `system2()`. Falls back to a pure-R decoder automatically if python3 is not on PATH. |
| `run_custom_oligo.R` | **Primary CLI entry point.** Template driver for any sequence -- copy it, edit the CONFIG block (sequence, chemistry overrides, parameters), and run. Works with standard chemistry out of the box. |
| `tests/` | Validation scripts (print-based, not testthat -- see note below). |
| `QUICKSTART.md` | From CoA/synthesis documentation to a valid input string, plus the fastest path to a first run. |
| `vignettes/OligoMetProfiler.Rmd` | Full package vignette: input notations, chemistry dictionary and overrides, library generation, masses/envelopes/isotopes, fragment ions, workbook/report/acquisition exports, MS matching. |

## Running the pipeline on your own sequence

Copy `run_custom_oligo.R`, edit the `CONFIG` block at the top (sequence,
optional custom chemistry overrides, parameters), and run:

```bash
Rscript run_custom_oligo.R
```

Or use the Shiny app (`app.R`) for an interactive session. The app's "Save
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

The chemistry engine is anchored to a published, known-mass reference
case (the ION337 gapmer, < 1 ppm). To re-run that regression check:

```r
library(OligoMetProfiler)   # or source the R/ modules from a clone
validate_ION337()
```

And from the shell, the validation scripts:

```bash
Rscript tests/test_metabolites.R
Rscript tests/test_mass_isotope.R
Rscript tests/test_fragments.R
Rscript tests/test_ms_matching.R
Rscript tests/test_outputs.R
```

These print computed values and comparisons to the console rather than
asserting pass/fail -- read the output rather than the exit code. Three
scripts (`test_mass_isotope.R`, `validate_vs_workbook.R`,
`inspect_workbook.R`) additionally cross-check against a reference
workbook, `ION337_Workbook.xlsx`, which is **not** distributed with this
package. Drop a copy into `tests/` or point to it with:

```bash
ION337_WORKBOOK_PATH=/path/to/ION337_Workbook.xlsx Rscript tests/validate_vs_workbook.R
```

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
  set independently in `app.R` and in the CLI driver, so changing a
  default in one place doesn't propagate to the other.

## License

MIT -- see [LICENSE.md](LICENSE.md).
