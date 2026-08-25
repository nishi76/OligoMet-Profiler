# OligoMet Profiler

> **FOR RESEARCH USE ONLY.** Not for diagnostic, clinical, or regulatory
> submission use. Every value this software produces is a computed
> prediction, not a measurement, and must be confirmed experimentally.
> Provided without warranty; the author accepts no liability for its use.
> See [DISCLAIMER.md](DISCLAIMER.md).

A general-purpose R pipeline for generating theoretical metabolite
libraries, charge envelopes, isotope patterns, McLuckey MS/MS fragment ions,
and PRM inclusion lists for any therapeutic oligonucleotide -- DNA, RNA,
2'OMe, 2'F, MOE, cEt, LNA, mixed PO/PS backbones, and custom chemistry via
dictionary overrides -- with optional matching against uploaded MS data.
Ships as an installable R package (`OligoMetProfiler`) with a Shiny
dashboard.

Background information is obtained from published work SynONIM (Lippens et al., JASMS 2024), the Eluforsen metabolite
profiling study (Kim et al., Mol Ther Nucleic Acids 2019), OligoDistiller
(Liu et al., Anal Chem 2025), and the FMVS automatic metabolite ID method
(Ye et al., J Chromatogr B 2025). Four approved oligonucleotide
therapeutics are bundled as working examples, one per modality class in
Takakusa et al. (2023) -- **nusinersen** (antisense SSO), **inotersen**
(antisense gapmer), **patisiran** (siRNA in LNP) and **givosiran**
(GalNAc-siRNA). They are illustrations, not default targets -- run the pipeline on any sequence.
See the [Bibliography](#bibliography) for sources.

**New here?** Start with [QUICKSTART.md](inst/help/QUICKSTART.md) -- This document will helps user 
through reading Certificate of Analysis / synthesis documentation and
constructing a valid input string. The full workflow documentation is the
package vignette, [vignettes/OligoMetProfiler.Rmd](vignettes/OligoMetProfiler.Rmd).

## Installation

As an R package, straight from GitHub:

```r
# install.packages(c("remotes", "shiny", "DT", "shinyFiles"))
remotes::install_github("nishi76/OligoMet-Profiler")
library(OligoMetProfiler)   # engine functions
OligoMetProfiler::run_app() # launches the dashboard
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

The dashboard carries its own documentation: a **Help & guides** panel
on the front page quick start, the sequence guide and the
modification reference.

If the user do not want to assemble triplet notation by hand, the **Manual
sequence entry** panel in the middle of the dashboard takes information provided by the chemical analysis file  -- bases, sugars and
linkages -- and it will build the sequence for you:

| Field | Example |
|---|---|
| Bases (5'->3') | `TSASTTTSATAATGSTGG` |
| Sugars | `eeeeeeeeeeeeeeeeee` (or just `e` for all positions) |
| Linkages | `sssssssssssssssss` (or just `s`) |

[SEQUENCE_GUIDE.md](inst/help/SEQUENCE_GUIDE.md) walks a first-time user
through reading those three lines off a chemical analysis file, and
through exporting the same sequence as a BioPharma Finder FASTA.

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

## BioPharma Finder compatibility

Triplet sequences copied from Thermo BioPharma Finder's sequence editor
(e.g. `Ad-pTd-pCd-pAd`) are accepted directly -- `p` is recognized as the
phosphodiester linkage, matching BPF's own notation, alongside this
pipeline's original `o`. `s` (phosphorothioate) already matched BPF's
convention. BPF's 5'/3' terminal modification codes (biotin, cAG/cAU/ARCA/
mCAP cap analogs, triantennary GalNAc) are available from the conjugate
dropdowns in the app under
descriptive names rather than BPF's own single-letter codes, since those
letters (`a`, `u`, `r`, `c`) already mean specific sugars or linkages in
this dictionary. Their formulas come straight from the BioPharma Finder
5.2 Oligonucleotide Analysis User Guide's "Modification notation" topic.

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

The app export targeted mass lists formatted
for the Orbitrap Exploris Method Editor's Targeted Mass filter (Mass List
Type "m/z & z", Time Mode "Start/End Time" -- select that when importing,
or "Unscheduled" if you'd rather ignore the time columns). Column layout
follows the Orbitrap Exploris 240 Software Manual's "Targeted Inclusion --
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

## Author

**Nishikant Wase, PhD** -- author and developer
Research Scientist, Thermo Fisher Scientific
<nishikant.wase@gmail.com>

Designed, written and maintained by the author. If this software
contributes to work you publish, please cite it as:

> Wase, N. *OligoMetProfiler: theoretical metabolite libraries, charge
> envelopes, and MS/MS fragment ions for therapeutic oligonucleotides.*
> R package. <https://github.com/nishi76/OligoMet-Profiler>

## Disclosure and disclaimer

**Research use only.** OligoMetProfiler is a research tool. It is not a
medical device and is not intended or validated for diagnostic use,
clinical decision-making, patient care, quality control release testing,
or inclusion in a regulatory submission. Any such use is outside the
scope of this software and is at the user's own risk.

**Everything it reports is a prediction.** Metabolite libraries,
formulas, masses, envelopes, isotope patterns, fragment ions, target
lists and spectral libraries are computed from a chemistry dictionary.
Confirm every assignment experimentally.

**No warranty, no liability.** This software is provided "as is",
without warranty of any kind, express or implied, under the MIT licence.
In no event shall the author be liable for any claim, damages or other
liability including any loss of data, wasted instrument time or
materials, or erroneous scientific conclusion arising from this software or its use. The user is solely responsible
for verifying that it is fit for their purpose and for validating every
result it produces.

**Affiliation and disclosure.** The author is employed as a Research
Scientist at Thermo Fisher Scientific. This is disclosed because the
package interoperates with Thermo Fisher products BioPharma Finder
sequence notation and Orbitrap Exploris method-editor mass lists.
OligoMetProfiler is an independent personal project: it is not a Thermo
Fisher Scientific product and has not been supplied, reviewed,
supported, endorsed or approved by Thermo Fisher Scientific or any other
company. All interoperability is implemented from publicly documented
formats only.

**No endorsement, and trademarks.** Nothing here implies endorsement by
any instrument vendor, software vendor or pharmaceutical company named
in the package. Product and drug names are used for identification and
interoperability only and remain the trademarks of their owners. The
bundled reference drugs are working examples from published literature
and public regulatory filings. All views and outputs are the author's
own and do not represent Thermo Fisher Scientific or any other employer
or institution.

The full statement is in [DISCLAIMER.md](DISCLAIMER.md). From R: `oligomet_about()`.

## License

MIT -- see [LICENSE.md](LICENSE.md). Copyright (c) 2026 Nishikant Wase.
