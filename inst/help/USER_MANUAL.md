# OligoMet Profiler — User Manual

A complete walk-through, from typing a sequence to a quantified, annotated
metabolite-identification experiment. It follows one storyline: define the
drug, predict its metabolites, run a forced-degradation study to find which
ones are real, build MS1 and empirical MS2 libraries, process the whole set
in batch, annotate the hits, quantify the ones you care about, and report
only the species that rise above 0.1% of the parent.

The app has these top-level tabs, used roughly in this order:

1. **Library Generation** — define the oligo, predict the metabolite library, preview structures.
2. **Empirical MS2 Library** — build a spectral library from real data.
3. **Batch Processing** — run many raw files, match, confirm, quantify.
4. **Statistical Analysis** — data matrix, univariate, class, time course, multivariate.
5. **Ask OligoMet** — the built-in assistant.
6. **Help & Quick Start** — the guides, including this manual.

---

## 1. Launching the app

**Hosted.** Open the public instance at the link in the project README. Nothing to install. Your uploaded data stays in that session.

**Local (R).** From a checkout:

```r
shiny::runApp("inst/app")
```

You need R 4.2+, plus `shiny`, `bs4Dash`, `DT`, `ggplot2`, and for in-app
guide rendering `markdown`. Isotope modelling uses `enviPat` when present and
falls back to a built-in convolution when it is not. Raw-file deconvolution
calls a bundled Python module; the app finds `python3` on the PATH.

---

## 2. Building the sequence

This is the one input everything else depends on. Get it right and the whole
library is correct.

### 2.1 Three ways to enter an oligo

**Plain bases.** Type the base sequence 5' to 3', e.g. `ATGCATGC`. Pick a
global sugar and linkage from the dropdowns. Good for a quick uniform-chemistry
oligo.

**Triplet / per-residue notation.** Each residue is written as
`<sugar><base>` and each linkage as a one-letter code between residues. This is
the precise form, and it is what you want for a gapmer or a mixed-chemistry
siRNA strand. The grammar:

- **Base**: `A C G T U`, plus modified bases such as `S` (5-methyl-C), `m5C`, `hm5C`, `D` (2,6-diaminopurine), `I` (hypoxanthine/inosine), `psi` (pseudouridine).
- **Sugar** (2' position): `d` = DNA (2'-deoxy), `r` = RNA (ribose), `m` = 2'-O-methyl, `f` = 2'-F, `e` or `MOE` = 2'-MOE, `cEt`, `LNA`.
- **Linkage**: `s` = phosphorothioate (PS), `o` = phosphodiester (PO).

A 5-10-5 MOE gapmer, all-PS, reads like:

```
Ge-sTe-sCe-sAe-sTe-sGd-sAd-sCd-sTd-sGd-sAd-sCd-sTd-sGd-sAd-sGe-sTe-sCe-sAe-sTe
```

The **Sequence Guide** tab reproduces a full nusinersen worked example and the
complete code tables. Open it side by side the first time.

### 2.2 Terminal conjugates and caps

Set the 5' and 3' conjugate from their dropdowns, or add them inline in a FASTA
header. The dictionary carries GalNAc, trivalent GalNAc, cholesterol, fatty
acids, and amino-linkers. A terminal phosphate cap (`5'-phosphate`,
`3'-phosphate`) is also a conjugate choice. On truncation, a conjugate is
carried or lost with its own end, so a 3' exonuclease product drops a 3'
conjugate automatically.

### 2.3 FASTA and the internal standard

Paste one or more records in FASTA form, or load a `.fasta`/`.txt` file:

```
>Inotersen
Ge-sTe-sCe-s...
>Inotersen analog IS
Ae-sGe-s...
```

A record becomes the **internal standard** when its name holds a standalone
`IS`, `SIL`, `ISTD`, `IntStd`, or `IS1` token, the phrase "internal standard",
or `role=IS` in the header. The analyte gets the full metabolite library; the
IS is a single library entry with no metabolites of its own, and its mass comes
from its own sequence. One analyte per run. See the **Internal Standard**
section below for how its signal normalizes everything else.

### 2.4 Custom chemistry

If your residue is not in the dictionary, add it in the Custom Chemistry table:
a code, a molecular formula, and a type (base, sugar, linkage, conjugate). The
formula is validated before it reaches any mass. A sulfur-bearing custom
linkage is treated as a phosphorothioate for the desulfurization series.

---

## 3. Predicting the metabolite library

Open **Library Generation** and run Phase 1. The app builds the theoretical
metabolite library from the sequence alone. No data needed yet.

### 3.1 What gets predicted

- **Parent** (full length).
- **3' exonuclease shortmers** (n-1, n-2, ... up to your 3' cap).
- **5' exonuclease shortmers** (n-1, n-2, ... up to your 5' cap).
- **Endonuclease internal fragments** at each cleavage site, both mass-balanced termini (the bridging phosphate lands on one side, so each cut gives a free-OH product and a terminal-phosphate product, with the correct phosphorothioate cap where the cut bond was PS).
- **Gapmer-targeted cleavage**: set `endo_sites = "gap"` to cut only inside the DNA gap, the dominant route for a MOE/DNA/MOE gapmer.

Two mass dimensions ride on top of every species during matching:

- **PS to PO desulfurization** (0 to k events, each -16 Da).
- **Adducts**: H, Na, K, NH4.

### 3.2 In-place chemical modifications (optional)

Beyond chain shortening, four in-place modifications can be added to the search
from the sidebar (all off by default):

- **Depurination (abasic)** — loss of an A or G base.
- **Oxidation (+O)** — one added oxygen.
- **Deamination (+0.984 Da)** — C to U, 5-methyl-C to T, A to hypoxanthine, and the other aminated bases.
- **2',3'-cyclic phosphate (-H2O)** — on a 3'-phosphate terminus that sits on a ribonucleotide.

Each is modeled as a single event layered on the oxidation, charge, and adduct
grid. MS1 cannot say which of several equivalent bases carried the change, so a
modification whose sites share a mass is reported with the position
unresolved; the structure view shows the candidate positions.

### 3.3 Previewing structures

The **Structure Annotation** sub-tab draws any metabolite as a cleavage map:
the parent sequence as a chain of residues colored by 2'-sugar chemistry, the
retained part solid and the lost part faded, the backbone drawn thick for
phosphorothioate and thin-dashed for phosphodiester, a scissors at the
cleavage, the 5' and 3' termini labeled, and a badge on any in-place
modification. Pick a metabolite and a modification from the two dropdowns. Use
**Download SVG** to pull a figure for a report or slide. The other sub-tabs
(Charge Envelope, Truncation Series, Isotope Pattern, Oxidation Series) show
the mass and isotope detail for the selected species.

At the end of Phase 1 you have the library, a workbook, and a PRM/method list
you can export.

---

## 4. From predicted to identified: the forced-degradation study

Prediction gives you every species that could form. A forced-degradation (or
in-vitro metabolism) study tells you which ones actually do. The design that
feeds this app cleanly:

- **Time-zero control (t0).** The drug, unstressed. Establishes the parent and any synthesis impurities so you do not call them metabolites.
- **Stressed samples.** Nuclease digest (S1, phosphodiesterase, or a tissue/serum/liver homogenate), or chemical stress, sampled at several timepoints.
- **Matrix blank.** The matrix with no drug. Anything here is background, subtracted later.
- **Replicates.** At least duplicate per timepoint so the statistics have a footing.

Acquire each sample by LC-HRMS in negative mode, with data-dependent or
targeted MS2 on the dominant precursors. Keep the gradient and column identical
across the set so retention times line up.

The logic you are testing: a predicted species is **identified** when its
accurate mass matches within tolerance, its charge states agree on one neutral
mass, it rises over time in the stressed samples and is absent (or low) at t0
and in the blank, and, when a spectrum is available, its MS2 fragments confirm
the sequence.

---

## 5. MS1 features and the empirical MS2 library

### 5.1 MS1 feature detection

Raw files are deconvolved by the bundled Python pipeline: it finds regions of
interest, builds extracted-ion chromatograms, and reports features as
(m/z, retention time, intensity, area). You choose a noise model (a signal-to-
noise multiple of the local noise, or a fixed intensity floor) and a matching
tolerance in ppm. For a single file you can run this from the MS panel; for
many files use Batch Processing (next section).

### 5.2 Empirical MS2 library

The **Empirical MS2 Library** tab builds a spectral library from your real MS2
scans rather than from theory. For each confirmed precursor it keeps the best
acquired spectrum, so later matches are scored against spectra you measured on
your own instrument. This is what separates a confident sequence call from a
bare mass match. Compare an acquired spectrum to the predicted fragment ladder
in the MS2 mirror plot; a good match shows b/y (and a-B) ions walking along the
sequence.

The predicted-fragment browser lives under Library Generation (**MS2 Explorer**)
and is useful for method setup; the empirical library is what you trust for
identification.

### 5.3 Method export

Phase 1 produces a PRM / targeted inclusion list (precursor m/z by charge for
every library species, H-adduct) you can load onto the instrument to trigger
MS2 on the metabolites, not just the parent.

---

## 6. Running the whole experiment in batch mode

Open **Batch Processing**, enable it, and upload the raw files for the entire
study (t0, every stressed timepoint, blanks, and, for absolute quant, the
calibration standards and QCs).

### 6.1 The sample-information sheet

Supply a sample metadata table (CSV/TSV/Excel). The columns the app reads:

| Column | Meaning |
|---|---|
| `sample` | File/sample name, matching the raw file basename. |
| `sample_type` | `standard`, `quality_control`, `unknown`, or `blank`. |
| `concentration` | Nominal amount for `standard` rows (drives the calibration curve). |
| `group` | Experimental group / condition (e.g. stressed vs control). |
| `timepoint` | Numeric time, for degradation kinetics and trend tests. |

`sample_type` is not just a label. Calibration curves are fit from `standard`
rows with a positive `concentration`; QCs are checked against their nominal
value; blanks drive background subtraction and are never used as calibrators.

### 6.2 Parameters

Set the deconvolution ppm, the noise model, charge range, maximum
desulfurization level, adducts, and (if you want them) the in-place
modification classes. Turn on MS2 confirmation so each hit is scored against
the empirical library. Set the number of parallel workers to match the machine.

Turn on **Align retention times across samples** when column or gradient drift
moves a peak's apex from run to run, across a long sequence or across matrices.
The alignment groups features across samples by mass and retention time, picks
stable groups as anchors, and fits a smooth per-sample correction from them
(the peakGroups method XCMS uses by default). Spiked internal standards are
seeded as anchors, since they are present in every injection. A sample with too
few anchors keeps its raw retention time and is flagged, and the raw value is
kept alongside the aligned one. Identity-based matching still runs on accurate
mass, so alignment sharpens the retention-time steps rather than driving the
identification.

### 6.3 Run

The batch step deconvolves every file, matches the library against each
sample's features (metabolite × modification × desulfurization × charge ×
adduct), checks charge-envelope consistency, confirms MS2 where spectra exist,
and computes degradation. Results land in the sub-tabs: **Batch Results**, **XIC
Quantitation**, **Unidentified Peaks**, **Calibration Curves**, **Internal
Standard**, **Degradation Summary**.

---

## 7. Annotating the results

Read the match table top down:

- **ppm error** and **isotope fit** — mass and pattern agreement.
- **charge-envelope consistency** — do the charge states back-calculate to one neutral mass. A real multiply-charged species does; a coincidental single-charge match does not.
- **mod / mod_site** — which in-place modification (if any), and the residue when MS1 can localize it.
- **MS2 confirmation score** — fragment coverage against the empirical library, with phosphorothioate diagnostic ions where relevant.
- **ambiguity flag** — set when more than one distinct metabolite can explain the same feature at this tolerance. Treat an ambiguous hit as unconfirmed until MS2 or retention time settles it.

Open the **Structure Annotation** view for any identified metabolite to see
exactly where it was cut and what changed. The **Unidentified Peaks** tab
retains everything that matched nothing; scan it for metabolites your library
did not predict (a non-nuclease product, a conjugate loss), which you can then
add and re-run.

---

## 8. Quantification: absolute for the few, relative for the rest

The usual plan: put a calibration curve behind the handful of molecules that
need a real concentration (the parent, a regulated metabolite), and report
everything else relative to the parent or to an earlier timepoint.

### 8.1 Internal standard normalization

If you entered an IS, it drives normalization. In **Calibration &
Quantification** the global IS defaults to your IS record and the signal basis
switches to IS-normalized, so curves and ratios are built on area ratios
(analyte area / IS area in the same sample) rather than raw area. Check the
**Internal Standard** tab for any sample whose IS response falls outside the
50-150% window, and assign a different IS to individual metabolites there if a
single global IS does not fit them all. The IS is never blank-corrected.

### 8.2 Absolute quantification (selected molecules)

Choose the molecules that get a calibration curve. For each, the app fits the
curve from the `standard` samples (optionally weighted), back-calculates
concentrations for the unknowns, and reports percent relative error against the
nominal value for standards and QCs. The **Calibration Curves** tab shows the
line, the calibrated range, and the standards and QCs on it. QCs outside their
window flag a quant you should not trust.

### 8.3 Relative quantification (everything else)

The remaining metabolites are quantified relative to a reference: as a fraction
of the parent, or as a fold-change over an earlier timepoint. This is where you
see degradation build. The summed-XIC engine (Chromeleon-style summation of the
top charge states and isotopes) gives each species one integrated area over a
shared window, so relative numbers are stable and re-summing a subset needs no
re-extraction.

---

## 9. Selecting molecules at or above 0.1% of the full-length compound

A metabolite-ID report should not list noise. The common cutoff is to keep only
species whose signal reaches some fraction of the parent (full-length) in at
least one sample. 0.1% is a typical floor.

How to apply it with the current outputs:

1. For each sample, take the parent's summed area as 100%.
2. For every other metabolite, compute its summed area as a percentage of that parent area in the same sample.
3. Keep a metabolite if it reaches 0.1% in any stressed sample. Drop the rest from the reported set.

The **Degradation Summary** and the XIC area table give you the per-sample areas
to do this; the relative-abundance view expresses each species against the
parent. Apply the 0.1% floor before you quantify, so the calibration and
relative-quant tables carry only real species.

> A one-click "keep metabolites >= X% of full-length" filter is a natural
> addition to the batch results. It is not yet a dedicated control; for now the
> step above reproduces it from the area tables. Ask if you want it wired in as
> a slider that prunes the results and every downstream table in one move.

---

## 10. Degradation and statistics

- **Degradation Summary** reports parent vs degradant composition, grouped by metabolite class and by experimental group/timepoint, with mean, SD, SEM, and replicate count. A class absent in a replicate counts as 0%, not missing, so averages are honest.
- **Univariate Statistics** and **Class Comparison** test which metabolites or classes differ between groups, with multiple-testing correction.
- **Time Course** fits a linear trend across timepoints and normalizes each metabolite to its own signal at a reference time, which is how you show a shortmer rising as the parent falls.
- **Multivariate** gives PCA and clustering across the sample set.
- The **Data Matrix** sub-tab exports the whole metabolite × sample matrix (long or wide) and re-imports a preprocessed matrix, so all of the above can run without re-processing raw files.

---

## 11. Export and reproducibility

- **Workbook** (xlsx): the library, charge envelopes, matches, quant, and degradation.
- **Report** (HTML): a shareable summary.
- **PRM / method list** (CSV): the targeted inclusion list.
- **Structure SVGs**: per-metabolite figures.
- **Session save/restore**: the app saves every parameter so a run is reproducible. The data matrix round-trip lets a reviewer rebuild the analysis from the matrix alone.

---

## 12. Troubleshooting

- **No calibration curve shows.** You need `standard` rows with positive `concentration` in the sample sheet. Without them there is nothing to fit.
- **A metabolite matches but looks wrong.** Check the ambiguity flag and the charge-envelope consistency, then confirm with MS2 before reporting it.
- **The IS flags many samples outside 50-150%.** Recovery or injection drifted; review the IS response plot before trusting any IS-normalized number.
- **A real metabolite is missing from the library.** It is probably a class you did not enable (turn on the relevant in-place modification) or a non-nuclease product; look in Unidentified Peaks, add the chemistry, and re-run.
- **Deconvolution does not start.** The Python module needs `python3` on the PATH. The hosted app has it; a local run may not.

---

*For the exact notation tables see the Sequence Guide; for the chemistry
dictionary see Modifications; for the command-line path see No-Shiny (CLI).*
