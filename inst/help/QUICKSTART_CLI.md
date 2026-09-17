# Quick Start: running without the Shiny app

> **FOR RESEARCH USE ONLY.** All outputs are computed predictions, not
> measurements — confirm every assignment experimentally. No warranty;
> see DISCLAIMER.md.

This covers the same pipeline as [QUICKSTART.md](QUICKSTART.md), but
entirely from R scripts — no `shiny`/`DT`/`shinyFiles` needed, no
browser, nothing to click. Useful for a headless server, an HPC/cluster
job, or just scripting a batch of runs. Three entry points, one per use
case:

| Script | Use case |
|---|---|
| `run_custom_oligo.R` | One sequence, optionally one MS file |
| `run_batch_ms.R` | One sequence, many MS files in parallel, with groups/timepoints/quantification |
| `run_app()` | The Shiny dashboard (see QUICKSTART.md instead) |

All three live in the repository root (or, once the package is
installed, `run_app()` is the only one that ships as an actual package
function — copy the other two from the repository or
`system.file(package = "OligoMetProfiler")` if you installed via
`remotes::install_github()`).

## 1. Install (no Shiny needed)

```r
remotes::install_github("nishi76/OligoMet-Profiler")
```

or from a clone:

```r
install.packages(c("openxlsx", "ggplot2", "xml2", "xfun", "jsonlite"))
```

Batch (multi-file) processing additionally needs Python 3.9+ on PATH
with `inst/python/requirements.txt` installed — see that file, or the
"Batch processing" section of README.md. Single-file mode works without
Python at all.

## 2. Single sequence, single (optional) MS file: `run_custom_oligo.R`

Copy `run_custom_oligo.R`, edit its `CONFIG` block (your sequence, in
any of the notations QUICKSTART.md/SEQUENCE_GUIDE.md cover, plus the
usual mass/matching parameters), then:

```
Rscript run_custom_oligo.R
```

This builds the theoretical library, writes the Excel workbook/HTML
report/PRM lists/spectral libraries, and — if you pointed `MS_FILE` at
an `.mzML`/`.mzXML`/peak-list file in the CONFIG block — matches it
against the library. No group/timepoint comparison or quantification at
this scale; that's what `run_batch_ms.R` is for.

## 3. Many files in parallel, with groups/timepoints/quantification: `run_batch_ms.R`

Copy `run_batch_ms.R`, edit its `CONFIG` block, then:

```
Rscript run_batch_ms.R
```

The parts of the CONFIG block worth knowing about beyond the sequence
itself:

**Sample metadata.** `SAMPLE_META` is a `sample, group, timepoint,
sample_type, concentration` data.frame — `sample` must match each
file's basename (no extension). Fill in `group` (two or more groups) OR
`timepoint` (a time course), not both. `sample_type` is one of
`unknown` / `standard` / `quality_control` / `reagent_blank` /
`matrix_blank`, defaulting to `unknown`; `concentration` is the nominal
value for a `standard`/`quality_control` row, blank otherwise. This is
the same CSV shape the Shiny app's batch sample table accepts — see
`inst/extdata/batch_example/sample_meta.csv` for a minimal example, or
`utils::read.csv()` a filled-in copy of the app's own "Download CSV
template" straight into `SAMPLE_META`.

**Background/noise threshold.** `PARAMS$sn_threshold` (default `3`):
absolute intensity floor = (each file's own noise level) × this number,
computed separately per file (see `estimate_ms_noise_level()` /
`?run_batch_deconvolution`) rather than one fixed number applied to
every file regardless of its background level. Set `PARAMS$sn_threshold
<- NULL` to use `PARAMS$min_intensity` (a fixed absolute value) instead.

**Quantification.** `ABSOLUTE_QUANT_MET_IDS` names which metabolite IDs
(`met$id` — print `mets` after a run, or check the workbook's metabolite
sheet, to find them) get absolute quantification from their own
calibration curve, built from `SAMPLE_META` rows with `sample_type ==
"standard"` and a `concentration`. Every other metabolite present in the
data is relative-quantified instead: fold-change vs. the earliest
timepoint for a time-course design, or vs. `CONTROL_GROUP` for a
group-comparison design. Leave `ABSOLUTE_QUANT_MET_IDS` empty
(`character(0)`, the default) for relative quantification across the
board — see `?quantify_metabolites` for the full contract, and
`?fit_calibration_curve` for `CALIBRATION_WEIGHTING`'s options.

Outputs land in `PARAMS$results_dir` (`results_batch/` by default):
the combined feature table, MS1 matches, unmatched peaks, statistics
table, `*_absolute_quantification.csv` / `*_relative_quantification.csv`
(when applicable), and the Excel workbook. Console output reports
per-file noise thresholds when S/N mode is active, and prints a
one-line reason for any file that failed rather than just its name.

## 4. Function-level help

Every function mentioned above (and the rest of the quantification/
statistics/matching suite) has a normal R help page once the package is
loaded — `?quantify_metabolites`, `?run_batch_deconvolution`,
`?compare_two_groups`, `?extract_ms1_features`, and so on. `??quantify`
or `??degradation` to search by keyword.

---

*OligoMetProfiler — Nishikant Wase, PhD. MIT licence. Research use
only; see DISCLAIMER.md.*
