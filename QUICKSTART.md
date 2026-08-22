# Quick Start Guide

From zero to a full theoretical metabolite library in about ten minutes.

## 1. Install and launch

Install the package from GitHub and launch the dashboard -- no clone
needed:

```r
# install.packages(c("remotes", "shiny", "DT", "shinyFiles"))
remotes::install_github("nishi76/OligoMet-Profiler")
OligoMetProfiler::run_app()
```

Or, from a clone of this repository:

```r
install.packages("shiny")
source("install_packages.R")     # installs remaining dependencies
shiny::runApp(".")
```

To drive the engine from R instead of the dashboard:

```r
library(OligoMetProfiler)
```

Or run the command-line driver on your own sequence: copy
`run_custom_oligo.R`, edit the CONFIG block at the top, and
`Rscript run_custom_oligo.R`.

## 2. Read your chemical analysis file

Before you can type the sequence in, you need three pieces of
information **per nucleotide position**, plus anything attached to the
ends. All of it is normally in the documentation that came with your
oligonucleotide — a Certificate of Analysis (CoA), a synthesis report,
a patent sequence listing, or a BioPharma Finder sequence-editor
export. Look for:

| What you need | Where it usually appears on the document |
|---|---|
| **Base sequence** (A/G/C/T/U, 5'→3') | "Sequence", often in IUPAC letters. Watch for 5-methylcytosine — CoAs write it as `mC`, `5-Me-C`, `C*`, or a footnote ("all cytosines are 5-methyl"). |
| **Sugar at each position** | "Chemistry", "Modification pattern", or a design shorthand like *"5-10-5 MOE gapmer"* (5 MOE wings, 10 DNA gap, 5 MOE) or *"2'-OMe/2'-F alternating"*. |
| **Backbone linkages** | "Backbone: full phosphorothioate (PS)" is most common for ASOs. Mixed backbones are listed per position or as a pattern ("PS wings, PO gap"). One oligo of *n* bases has *n−1* linkages. |
| **Terminal conjugates** | "5'-conjugate / 3'-conjugate": GalNAc, cholesterol, C6 amino linker, phosphate cap, dyes. If nothing is stated, there is none. |

Translate each item with these code tables:

**Bases:** `A` `G` `C` `T` `U`, plus `S` = 5-methylcytosine,
`D` = 2,6-diaminopurine, `I` = inosine.

**Sugars:** `d` = DNA (2'-deoxy), `r` = RNA, `m` = 2'-OMe,
`f` = 2'-F, `e` = MOE (2'-O-methoxyethyl), `cEt`, `LNA`.

**Linkages:** `s` = phosphorothioate (PS), `o` or `p` = phosphodiester
(PO — `p` matches BioPharma Finder's notation).

**Conjugates:** `GalNAc`, `GalNAc3`, `cholesterol`, `C6`, `C12`,
`TEG`, `FAM`, `Cy3`, `5'-phosphate`, `3'-phosphate` (see the app's
dropdown or `chemistry_dict.R` for the full list).

## 3. Construct the input string

Write one dash-separated token per nucleotide, 5'→3'. Each token is:

```
[linkage-coming-into-this-position][BASE][sugar]
```

The first (5') token has **no** linkage prefix — there is no incoming
bond.

### Worked example: inotersen

Take a real one. Inotersen's published description reads:

> 20-mer antisense oligonucleotide targeting the TTR 3'UTR (bases
> 618-637). Sequence: 5'-TCTTGGTTACATGAAATCCC-3'. Design: 5-10-5
> gapmer; positions 1-5 and 16-20 are 2'-MOE, positions 6-15 are
> 2'-deoxy. Backbone: full phosphorothioate. All cytosines are
> 5-methylcytosine. No terminal conjugates.

Build it position by position:

| Pos | Base | Base code | Sugar | Incoming linkage | Token |
|----:|------|-----------|-------|------------------|-------|
| 1 | T | `T` | 2'-MOE → `e` | (none, 5' end) | `Te` |
| 2 | 5-Me-C | `S` | 2'-MOE → `e` | PS → `s` | `sSe` |
| 3–5 | T, T, G | `T` `T` `G` | 2'-MOE → `e` | PS → `s` | `sTe-sTe-sGe` |
| 6 | G | `G` | DNA → `d` (gap starts) | PS → `s` | `sGd` |
| 7–15 | ...gap... | | DNA → `d` | PS → `s` | `sTd-sTd-sAd-sSd-...` |
| 16 | A | `A` | 2'-MOE → `e` (wing resumes) | PS → `s` | `sAe` |
| 17–20 | T, C, C, C | `T` `S` `S` `S` | 2'-MOE → `e` | PS → `s` | `sTe-sSe-sSe-sSe` |

Joined with dashes:

```
Te-sSe-sTe-sTe-sGe-sGd-sTd-sTd-sAd-sSd-sAd-sTd-sGd-sAd-sAd-sAe-sTe-sSe-sSe-sSe
```

Paste that into the app's sequence box, or into the CONFIG block of the
CLI driver. **How you know it's right:** the computed molecular formula
should match the published one. For inotersen that is
`C230H318N69O121P19S19`, average MW 7183.08 — and this pipeline
reproduces it exactly. Run `validate_reference()` in R to see that check
run against both inotersen and nusinersen.

Sanity checks the parser enforces: token count = oligo length; every
token after the first starts with a linkage code; every code must exist
in the dictionary (you get a clear error naming any unknown code).

### The four bundled examples

Rather than typing them out, four approved drugs ship with the package —
one per modality class — and are loadable from the app's "Load example"
dropdown:

| Modality | Drug | Constant | Chemistry |
|---|---|---|---|
| Antisense (SSO) | nusinersen | `NUSINERSEN_TRIPLET` | PS, uniform 2'-MOE, 5-methyl pyrimidines |
| Antisense (gapmer) | inotersen | `INOTERSEN_TRIPLET` | PS, 5-10-5 2'-MOE/DNA gapmer |
| siRNA (LNP) | patisiran (sense) | `PATISIRAN_SENSE_TRIPLET` | 2'-OMe/ribose, all-PO, dTdT overhang |
| siRNA (GalNAc) | givosiran (sense) | `GIVOSIRAN_SENSE_TRIPLET` | 2'-OMe/2'-F, partial PS, 3' GalNAc |

siRNA drugs are given as their **sense strand** — the pipeline profiles
one strand at a time, so run each strand of a duplex separately. The
GalNAc conjugate on givosiran is set through the `conj3` field (the app
sets it for you when you load that example), because triplet notation
has no conjugate field.

### Alternative notations (also accepted, auto-detected)

- **BioPharma Finder triplets** — copy the sequence straight out of
  BPF's sequence editor (`Ad-pTd-pCd-pAd` style); `p` and `s` parse
  directly.
- **OligoDistiller** — `OH-Gm*-Tm*-...-Gm-OH`, where `*` = PS outgoing
  bond.
- **Structured R list** — explicit `bases` / `sugars` / `linkages`
  vectors plus `conj5` / `conj3`; required when you have terminal
  conjugates. See the vignette for a GalNAc siRNA example.

### Common pitfalls

- **5-methylcytosine**: 2'-MOE antisense drugs almost always use 5-Me-C
  at every cytosine. If the CoA says so, use base code `S`, not `C` —
  the mass differs by CH₂ (14.016 Da) per site. In inotersen that is
  five sites, so getting it wrong shifts the parent mass by 70 Da.
- **5-methyluracil is thymine**: CoAs for MOE drugs often write `MeU` or
  `5-Me-U`. That is the same nucleobase as thymine — use base code `T`.
- **Linkage count**: *n* bases means *n−1* linkage codes. The parser
  puts the linkage on the *following* token.
- **A modification you don't have a code for**: add it as a custom
  override with its elemental formula (app: custom chemistry table;
  CLI: `custom_overrides` block) — see "Custom chemistry overrides" in
  the vignette.

## 4. Run and collect outputs

Click **Run** in the app (or run the CLI driver). You get:

- `*_library.xlsx` — 7-sheet workbook (library, charge envelopes,
  PS-oxidation series, fragment ions, PRM list, summary),
- `*_report.html` — plots and summary report,
- Orbitrap Exploris **MS1 inclusion** and **MS2 PRM** target lists
  (CSV, Method Editor import-ready).

## 5. (Optional) match against your LC-MS data

Upload an `mzML`/`mzXML` file (or a two-column m/z–intensity peak
list) in the app's MS panel, or set `MS_FILE` in the CLI driver.
Vendor raw files convert automatically when ProteoWizard `msconvert`
is on the PATH. The MS Matching sheet of the workbook then reports
matched features with ppm error, isotope-fit, envelope-consistency,
and MS2 confirmation scores.
