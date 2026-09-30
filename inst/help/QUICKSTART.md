# Quick Start Guide

> **FOR RESEARCH USE ONLY.** All outputs are computed predictions, not
> measurements — confirm every assignment experimentally. No warranty;
> see DISCLAIMER.md.

## 1. Install and launch

```r
# install.packages(c("remotes", "shiny", "DT", "shinyFiles"))
remotes::install_github("nishi76/OligoMet-Profiler")
OligoMetProfiler::run_app()
```

Or from a clone: `source("install_packages.R")` then `shiny::runApp(".")`.
For scripted runs, copy `run_custom_oligo.R`, edit its CONFIG block, and
`Rscript run_custom_oligo.R`.

## 2. Read your documentation

From your Certificate of Analysis or synthesis report you need, per
position: the **base sequence** (watch for 5-methylcytosine — `mC`,
`5-Me-C`, or a footnote), the **sugar** (e.g. "5-10-5 MOE gapmer"), the
**backbone linkages** (*n* bases = *n−1* linkages), and any **terminal
conjugates**.

Codes:

- **Bases:** `A` `G` `C` `T` `U`, `S` = 5-methyl-C, `D` = 2,6-diaminopurine, `I` = inosine
- **Sugars:** `d` DNA, `r` RNA, `m` 2'-OMe, `f` 2'-F, `e` MOE, `cEt`, `LNA`
- **Linkages:** `s` = PS, `o` or `p` = PO
- **Conjugates:** `GalNAc`, `cholesterol`, `C6`, `TEG`, `FAM`, … (see the app dropdown)

## 3. Enter the sequence

**Easiest: the three-line entry panel.** Type the bases, sugars, and
linkages as separate lines; a single code is applied to every position.
Nusinersen:

| Field | Value |
|---|---|
| Bases (5'→3') | `TSASTTTSATAATGSTGG` |
| Sugars | `e` |
| Linkages | `s` |

Click Submit and the app builds the triplet string and shows the
computed formula and mass. From R: `parse_three_line("TSAS...", "e", "s")`.
New to this? [SEQUENCE_GUIDE.md](SEQUENCE_GUIDE.md) walks through it slowly.

**Explicit: triplet tokens.** One dash-separated token per position,
5'→3': `[linkage][BASE][sugar]`, with no linkage on the first token.
Inotersen (5-10-5 MOE gapmer, full PS, all C are 5-Me-C):

```
Te-sSe-sTe-sTe-sGe-sGd-sTd-sTd-sAd-sSd-sAd-sTd-sGd-sAd-sAd-sAe-sTe-sSe-sSe-sSe
```

**Check yourself:** the computed formula should match the published one
(inotersen: `C230H318N69O121P19S19`, avg MW 7183.08 — reproduced exactly;
see `validate_reference()`).

Also accepted, auto-detected: **BioPharma Finder triplets**
(`Ad-pTd-pCd`), **OligoDistiller** notation, and a structured R list
(needed for terminal conjugates).

**Bundled examples** (app's "Load example" dropdown): nusinersen,
inotersen, patisiran (sense), givosiran (sense). siRNA duplexes are run
one strand at a time.

### Common pitfalls

- **5-Me-C** must be `S`, not `C` — 14.016 Da per site.
- **5-Me-U is thymine** — use `T`.
- ***n* bases need *n−1* linkage codes** (on the *following* token).
- **Unknown modification?** Add it as a custom override with its formula
  (app: custom chemistry table) — see [MODIFICATIONS.md](MODIFICATIONS.md).

## 4. Run and collect outputs

Click **1. Generate Library**. You get the 7-sheet Excel workbook, an
HTML report, Orbitrap Exploris MS1 inclusion and MS2 PRM target lists
(CSV), and MGF/MSP spectral libraries — no MS data needed yet. Save the
session (below) if you want to come back and import acquired data later.

## 5. (Optional) match against LC-MS data

Upload an mzML/mzXML/raw file or two-column peak list in the MS panel,
then click **2. Import & Process MS Data** (enabled once Step 1 has run).
Vendor raw files (Thermo `.raw`, Sciex `.wiff`, Bruker `.baf`/`.yep`)
convert automatically to `.mzML` when ProteoWizard `msconvert` is on the
PATH (install from proteowizard.org) — without it, upload an `.mzML`/
`.mzXML` export instead. Agilent/Bruker `.d` folders can only be reached
through the batch mode's local-folder input, not file upload. The MS
Matching sheet reports matches with ppm error, isotope-fit,
envelope-consistency, and MS2 confirmation scores.

## 6. Internal standard

Add the IS as a second FASTA record with IS in its name, e.g.
`>ION582 analog IS`, or `role=IS` in the header; load both from a `.fasta`
file if you prefer. The IS gets no metabolites. After processing, the
global IS in **Calibration & Quantification** defaults to it and Signal basis
switches to IS-normalized, so curves are fitted on area ratios. Check the
**Internal Standard** tab for samples outside the 50-150% response window,
and assign a different IS to individual metabolites there if needed.

## 7. Quantitation signal: summed XICs

With raw files, **Targeted summed-XIC quantitation** (Batch Processing >
Processing Options, on by default) sums the XICs of the chosen charge
states x isotopes and integrates them once per sample, the way Chromeleon
sums a component's ions. Defaults: 5 charge states (ranked on the
standards) x 5 isotopes (ranked by theoretical abundance), +/-10 ppm,
+/-0.5 min search window. Open the **XIC Quantitation** tab to see where
the signal is. If the charge-state distribution differs between standards
and samples, type a fixed list covering it (e.g. `3-8`) in **Fixed charge
states**. The counts and the list re-sum instantly.

## 8. Start from a data matrix instead of raw files

Once identification has run, download the data matrix from **Batch
Processing > Batch Results**. Next time, open **Batch Processing**, set
**Data Source** to *Pre-processed data matrix*, and upload:

1. the matrix (`.csv`, `.tsv`, `.txt`, or `.xlsx`, wide or long), and
2. the sample information sheet (`sample`, `group`, `timepoint`,
   `sample_type`, `concentration`). Either file can go first.

Mark calibrators as `sample_type = standard` with a concentration, QCs as
`quality_control`, blanks as `reagent_blank`/`matrix_blank`. Select the
metabolites to calibrate under **Calibration & Quantification**. The
**Calibration Curves** tab shows each curve, its fit, and back-calculated
concentrations for QCs and unknowns. **Degradation Summary** shows %
parent remaining and % degradation against the reference timepoint (set on
the Statistical Analysis tab, default earliest). No peak picking runs on
this path. Try it with `inst/extdata/data_matrix_example/`.

---

*OligoMetProfiler — Nishikant Wase, PhD. MIT licence. Research use
only; see DISCLAIMER.md.*
