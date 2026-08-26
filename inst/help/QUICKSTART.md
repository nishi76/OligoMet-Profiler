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

Click **Run**. You get the 7-sheet Excel workbook, an HTML report,
Orbitrap Exploris MS1 inclusion and MS2 PRM target lists (CSV), and
MGF/MSP spectral libraries.

## 5. (Optional) match against LC-MS data

Upload an mzML/mzXML file or two-column peak list in the MS panel.
Vendor raw files convert automatically when ProteoWizard `msconvert` is
on the PATH. The MS Matching sheet reports matches with ppm error,
isotope-fit, envelope-consistency, and MS2 confirmation scores.

---

*OligoMetProfiler — Nishikant Wase, PhD. MIT licence. Research use
only; see DISCLAIMER.md.*
