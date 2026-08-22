# Quick Start Guide

From zero to a full theoretical metabolite library in about ten minutes.

## 1. Install and launch

```r
# One-time setup
install.packages("shiny")
# from a clone of this repository:
source("install_packages.R")     # installs remaining dependencies

# Launch the dashboard
shiny::runApp(".")
```

Or install as a package and drive it from R directly:

```r
remotes::install_github("nishi76/OligoMet-Profiler")
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

### Worked example

Your CoA says:

> 16-mer antisense oligonucleotide. Sequence: 5'-GTCTCTCTCTTCTCTG-3'.
> Design: 5-6-5 gapmer; positions 1–5 and 12–16 are 2'-O-methyl,
> positions 6–11 are 2'-deoxy. Backbone: full phosphorothioate.
> No terminal conjugates.

Build it position by position:

| Pos | Base | Sugar | Incoming linkage | Token |
|----:|------|-------|------------------|-------|
| 1 | G | 2'-OMe → `m` | (none, 5' end) | `Gm` |
| 2 | T | 2'-OMe → `m` | PS → `s` | `sTm` |
| ... | | | | |
| 6 | T | DNA → `d` | PS → `s` | `sTd` |
| ... | | | | |
| 16 | G | 2'-OMe → `m` | PS → `s` | `sGm` |

Join with dashes:

```
Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm
```

Paste that into the app's sequence box, or into the `--seq` flag /
CONFIG block of the CLI drivers. Sanity checks the parser enforces:
token count = oligo length; every token after the first starts with a
linkage code; every code must exist in the dictionary (you get a clear
error naming any unknown code).

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

- **5-methylcytosine**: gapmer ASOs almost always use 5-Me-C in the
  wings. If the CoA says so, use base code `S`, not `C` — the mass
  differs by CH₂ (14.016 Da) per site, which is very noticeable.
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
