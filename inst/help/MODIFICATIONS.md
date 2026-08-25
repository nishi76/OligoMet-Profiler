# Modification reference

Every sugar and backbone chemistry the dictionary knows, what it is for,
what it does to the mass, and how to enter it. If your modification is
not here, section 5 shows how to add it yourself in about a minute --
you do not need to modify the package.

New to oligonucleotide notation? Read
[SEQUENCE_GUIDE.md](SEQUENCE_GUIDE.md) first; this page assumes you know
what the bases/sugars/linkages lines are.

---

## 1. How a mass difference is computed

The formula engine is additive:

```
oligo = sum(bases) + sum(sugars) + sum(linkages) - (2n-1) H2O + conjugates
```

Two consequences that make the tables below directly usable:

- **A sugar swap costs exactly the residue difference, per position.**
  Changing one position from DNA to MOE adds 74.0368 Da. Changing all 18
  positions of an 18-mer adds 18 x 74.0368 = 1332.66 Da.
- **A linkage swap costs the residue difference, per bond.** An n-mer has
  n-1 bonds, so converting an 18-mer from PO to PS adds
  17 x 15.9772 = 271.61 Da.

All numbers here are **monoisotopic**, in daltons, and come from
`tests/test_modifications.R`, which recomputes them from the dictionary
on every run. Average masses differ slightly; the app reports both.

---

## 2. Sugar modifications

Residue mass is the free sugar; a nucleoside is `base + sugar - H2O`.

| Code | Modification | Formula | Residue mass | vs DNA (`d`) | vs RNA (`r`) |
|---|---|---:|---:|---:|---:|
| `d` | 2'-deoxyribose (DNA) | C5H10O4 | 134.0579 | 0 | −15.9949 |
| `r` | ribose (RNA) | C5H10O5 | 150.0528 | +15.9949 | 0 |
| `m` | 2'-O-methyl | C6H12O5 | 164.0685 | +30.0106 | +14.0157 |
| `f` | 2'-fluoro (ribo) | C5H9O4F | 152.0485 | +17.9906 | +1.9957 |
| `FANA` | 2'-fluoro-arabino | C5H9O4F | 152.0485 | +17.9906 | +1.9957 |
| `e` / `MOE` | 2'-O-methoxyethyl | C8H16O6 | 208.0947 | +74.0368 | +58.0419 |
| **`NMA`** | **2'-O-NMA** | **C8H15NO6** | **221.0899** | **+87.0320** | **+71.0371** |
| `allyl` | 2'-O-allyl | C8H14O5 | 190.0841 | +56.0262 | +40.0313 |
| `AP` | 2'-O-aminopropyl | C8H17NO5 | 207.1107 | +73.0528 | +57.0578 |
| `UNA` | unlocked (2',3'-seco) | C5H12O5 | 152.0685 | +18.0106 | +2.0157 |
| `GNA` | glycol nucleic acid | C3H8O3 | 92.0473 | −42.0106 | −58.0055 |
| `LNA` | locked nucleic acid | C6H10O5 | 162.0528 | +27.9949 | +12.0000 |
| `ENA` | 2'-O,4'-C-ethylene | C7H12O5 | 176.0685 | +42.0106 | +26.0157 |
| `cEt` | constrained ethyl | C7H12O4 | 160.0736 | +26.0157 | +10.0207 |
| `NH2` | 2'-amino-2'-deoxy | C5H11NO4 | 149.0688 | +15.0109 | −0.9840 |
| `NP` | 3'-amino (N3'→P5') | C5H11NO3 | 133.0739 | −0.9840 | −16.9789 |

### 2'-O-NMA — the one most people are looking for

**What it is.** 2'-O-[2-(methylamino)-2-oxoethyl]: the 2'-OH hydrogen is
replaced by `-CH2-C(=O)-NH-CH3`. An amide where MOE has an ether.

**Why it exists.** MOE-like binding affinity, better metabolic
stability, and in splice-switching oligonucleotides better potency than
the MOE equivalent. It is a next-generation chemistry rather than an
established one -- there is no approved NMA drug yet.

**The number to remember.** NMA is **+12.9953 Da heavier than MOE** per
position, because its substituent (C3H6NO) differs from MOE's (C3H7O) by
exactly N minus H. On a uniformly modified 18-mer that is +233.9 Da
overall -- large, unambiguous, and easy to confirm on a deconvoluted
mass.

```r
parse_three_line("TSASTTTSATAATGSTGG", "NMA", "s")   # uniform NMA
```

### UNA — unlocked nucleic acid

The C2'-C3' bond of ribose is absent, so the residue is acyclic and
flexible. Formally ribose + H2, i.e. **+2.0157 Da vs RNA**. UNA
*lowers* duplex melting temperature -- it is used deliberately, usually
at one or two positions in an siRNA seed region, to suppress off-target
silencing rather than to add stability. Watch the mass: 2 Da is close
enough to an isotope spacing that a single UNA in a large oligo is easy
to misassign; check the full isotope pattern, not the monoisotopic peak
alone.

### GNA — glycol nucleic acid

The sugar is replaced outright by a three-carbon glycerol unit, making
this by far the largest per-position mass change in the table:
**−58.0055 Da vs RNA**. Like UNA, it is a destabilizer used positionally
in siRNA.

### ENA, LNA, cEt — the bridged (BNA) family

All three lock the sugar pucker for high affinity and nuclease
resistance. LNA is ribose plus a methylene bridge (which, because
forming two bonds costs 2 H, works out to exactly +12.0000 Da, one
carbon); ENA has an ethylene bridge, **+14.0157 Da more than LNA**; cEt
is a constrained ethyl on a deoxy scaffold. LNA, ENA and cEt are flagged
`verify = TRUE` in the dictionary -- confirm against a known mass before
relying on them for a formal assignment.

### FANA — isobaric, and worth knowing about

2'-deoxy-2'-fluoro-**arabino**nucleic acid has the same atoms as 2'-F
ribo, in a different 2' configuration. It is **exactly isobaric**: mass
spectrometry cannot distinguish `f` from `FANA`. The separate code
exists so your sequence record is chemically correct; it will never
change a computed mass.

### N3'→P5' phosphoramidate — a sugar change, not a linkage change

The 3'-OH becomes a 3'-NH2 and the phosphate bonds to nitrogen instead
of oxygen. In the additive model that lives entirely in the sugar
(`NP`, **−0.9840 Da vs DNA**, an O swapped for NH); keep the linkage as
an ordinary `o` or `s`.

---

## 3. Backbone (linkage) modifications

Every linkage is the phosphodiester bridge with its single non-bridging
`-OH` swapped for something else. An n-mer has n−1 of them.

| Code | Modification | Formula | Residue mass | vs PO | vs PS |
|---|---|---:|---:|---:|---:|
| `o` / `p` | phosphodiester (PO) | HO3P | 79.9663 | 0 | −15.9772 |
| `s` / `u` | phosphorothioate (PS) | HO2PS | 95.9435 | +15.9772 | 0 |
| `mp` | methylphosphonate | CH3O2P | 77.9871 | −1.9793 | −17.9564 |
| `prp` | propyl phosphonate | C3H7O2P | 106.0184 | +26.0520 | +10.0749 |
| `ibu` | isobutyl phosphonate (iBu) | C4H9O2P | 120.0340 | +40.0677 | +24.0905 |
| `chx` | cyclohexyl phosphonate (cHex) | C6H11O2P | 146.0497 | +66.0833 | +50.1062 |
| `mop` | methoxypropyl phosphonate (MOP) | C4H9O3P | 136.0289 | +56.0626 | +40.0854 |
| `pace` | phosphonoacetate (PACE) | C2H3O4P | 121.9769 | +42.0106 | +26.0334 |
| `tpace` | thiophosphonoacetate | C2H3O3PS | 137.9541 | +57.9877 | +42.0106 |
| **`msp`** | **mesyl phosphoramidate** | **CH4NO4PS** | **156.9599** | **+76.9935** | **+61.0164** |
| `pgo` | phosphoryl guanidine (cyclic) | C5H10N3O2P | 175.0511 | +95.0847 | +79.1076 |
| `tmg` | phosphoryl guanidine (Tmg) | C5H12N3O2P | 177.0667 | +97.1004 | +81.1232 |

### Mesyl phosphoramidate (msPA, "µ")

`-O-P(=O)(NH-SO2-CH3)-O-`. The non-bridging OH is replaced by a
methanesulfonamide. It gives PS-like nuclease resistance without PS's
promiscuous protein binding, and unlike PS it is achiral at phosphorus,
so there is no Rp/Sp diastereomer mixture. **+61.0164 Da per bond
relative to PS** -- on a fully modified 18-mer (17 bonds) that is
+1037.3 Da, so a partially converted batch is very visible in a
deconvoluted spectrum.

### PACE and thioPACE

Phosphonoacetate replaces the OH with `-CH2-COOH` (**+42.0106 Da vs
PO**); thioPACE is PACE with the P=O replaced by P=S, a further
**+15.9772 Da**. Both improve cell uptake and nuclease resistance and
are used in siRNA.

### Phosphoryl guanidine (PGO) — check which variant you have

A guanidine group replaces the OH, giving a **neutral** backbone with a
chiral phosphorus (so synthesis produces diastereomers). Two variants
are in common use and **they are not the same mass**:

| Code | Guanidine | vs PO |
|---|---|---:|
| `pgo` | 1,3-dimethylimidazolidin-2-imine (cyclic) | +95.0847 |
| `tmg` | tetramethylguanidine, Tmg (acyclic) | +97.1004 |

They differ by H2 — 2.0157 Da per bond, which on a heavily modified
oligo compounds into a large and very confusing discrepancy. Confirm
which one your synthesis used before assigning anything.

A neutral backbone also changes what you should expect in the mass
spectrometer: it ionizes far less readily in negative ESI, so the charge
envelope shifts to lower charge states than the `z_min`/`z_max` defaults
assume. Widen the range downward if a PGO or alkyl phosphonate oligo
gives no matches.

### Alkyl phosphonates — methyl, propyl, isobutyl, cyclohexyl, MOP

The same swap with a plain alkyl group in place of the OH, giving
neutral linkages. Placed at siRNA internucleotide positions 6–7 from the
5′ end they suppress off-target activity; in ASOs they steer RNase H1
cleavage.

| Code | Alkyl group | vs PO |
|---|---|---:|
| `mp` | methyl | −1.9793 |
| `prp` | propyl | +26.0520 |
| `ibu` | isobutyl | +40.0677 |
| `chx` | cyclohexyl | +66.0833 |
| `mop` | methoxypropyl | +56.0626 |

Methylphosphonate is the awkward one: at **−1.9793 Da vs PO** it is a
small shift, and close enough to other common differences that intact
mass alone will not settle it. Confirm those with MS/MS.

### A note on chiral phosphorothioates (Rp/Sp)

Stereopure PS oligonucleotides are chemically identical to the racemic
mixture -- same atoms, same mass. MS cannot distinguish Rp from Sp. The
dictionary's `u` code exists to let you *record* a stereochemistry
variant; it carries the same formula as `s` and will never change a
computed mass.

---

## 4. Entering a modification in the app

Three routes, all equivalent.

**Manual sequence entry** (middle of the dashboard). Multi-character
codes need a separator, and a single code applies to every position:

| Design | Bases | Sugars | Linkages |
|---|---|---|---|
| Uniform NMA, all msPA | `TSASTTTSATAATGSTGG` | `NMA` | `msp` |
| 5-10-5 NMA gapmer, all PS | `TSASTTTSATAATGSTGG` | `NMA,NMA,NMA,NMA,NMA,d,d,d,d,d,d,d,d,NMA,NMA,NMA,NMA,NMA` | `s` |
| LNA/DNA mixmer | `AGSTAGST` | `LNA-d-LNA-d-LNA-d-LNA-d` | `s` |

A field holding one code is recycled to every position, which is why
`NMA` and `msp` above are not repeated 18 times.

**Triplet notation** (sequence box). `[linkage][base][sugar]` per token;
the parser matches multi-character linkage codes longest-first, so
`mspANMA` reads as msPA + adenine + NMA:

```
TNMA-mspSNMA-mspANMA-mspSNMA-mspTNMA
```

**From R:**

```r
spec <- parse_three_line(
  bases    = "TSASTTTSATAATGSTGG",
  sugars   = c(rep("NMA", 5), rep("d", 8), rep("NMA", 5)),
  linkages = c(rep("msp", 5), rep("s", 7), rep("msp", 5)))
metabolite_mass_info(spec)$formula_str
#> "C220H312N81O132P17S17"
```

---

## 5. Adding a modification the dictionary does not have

You do not need to edit the package. In the app, use the **Custom
Chemistry** table in the sidebar: one row per code.

| Column | What to put |
|---|---|
| Code | Whatever you want to type in your sequence, e.g. `tcDNA` |
| Formula | The **residue** formula (see below), e.g. `C8H14O4` |
| Name | A human-readable label |
| Type | `base`, `sugar`, `linkage`, or `conjugate` |
| Attach | Conjugates only: `add`, `replace_H`, or `replace_OH` |

From R, pass the same thing as an override:

```r
dict <- build_dictionary(overrides = list(
  tcDNA = list(formula = "C8H14O4", name = "tricyclo-DNA", kind = "sugar")))
spec <- parse_three_line("AGST", "tcDNA", "s", dict = dict)
```

### Working out the residue formula

This is the only part that takes thought. The residue is what the
molecule contributes *before* condensation, and the engine subtracts the
water for you.

- **A 2'-O-R sugar** is ribose with the 2'-OH hydrogen replaced by R:
  `C5H9O5 + R`. For 2'-O-NMA, R = `-CH2C(O)NHCH3` = C3H6NO, giving
  C8H15NO6.
- **A bridged sugar** is the parent plus the bridge, minus 2 H for the
  two new bonds: LNA = ribose + CH2 − H2 = C6H10O5.
- **A linkage** is `HPO3` with the `-OH` replaced by your substituent X:
  `PO2 + X`. For PACE, X = `-CH2COOH` = C2H3O2, giving C2H3O4P.
- **A base** is the free nucleobase, not the nucleoside.

**Then check it.** Build a single nucleoside by hand and compare against
a known value -- `base + sugar - H2O` should reproduce the published
nucleoside formula. That one check catches nearly every mistake:

```r
f <- add_formulas(STANDARD_DICT[["A"]]$formula, dict[["tcDNA"]]$formula)
format_formula(f - c(H = 2, O = 1))
```

`tests/test_modifications.R` does exactly this for every built-in code;
copy a row from its `nucleoside_anchors` list to add your own.

---

## 6. What is deliberately not in the dictionary

Three chemistries do not fit the additive residue model, and a wrong
number is worse than a missing one:

- **Morpholino (PMO/PPMO)** replaces both the sugar (a morpholine ring)
  and the linkage (a phosphorodiamidate). The model can express it in
  principle, as a custom sugar plus a custom linkage, but the subunit
  formula must come from your own reference -- do not guess it.
- **Peptide nucleic acid (PNA)** has no sugar and no phosphate at all;
  its backbone is N-(2-aminoethyl)glycine with the base on an acetyl
  linker. The base-plus-sugar-plus-linkage decomposition does not
  describe it.
- **Boranophosphate** replaces a non-bridging oxygen with `BH3`, and
  boron is not in the engine's element table (`.ELEMENTS` in
  `chemistry_dict.R`). Adding it means adding boron's atomic masses
  *and* its isotope abundances (10B/11B is roughly 20/80, which visibly
  changes isotope patterns), not just a residue formula.
- **Tricyclo-DNA (tcDNA)** is a straightforward custom sugar -- it is
  left out only because the residue formula should come from a source
  you trust rather than from this package's guess. Add it per section 5.

---

## 7. Which formulas are provisional

Codes flagged `verify = TRUE` in `chemistry_dict.R` are best estimates
that have not been checked against an independently measured mass:
`cEt`, `LNA`, `ENA`, `AP`, `msp`, `pace`, `tpace`, `pgo`, `mp`, and the
BioPharma Finder terminal-modification conjugates. They are internally
consistent -- each is derived from the substitution rules above and
verified by `tests/test_modifications.R` -- but "internally consistent"
is not "confirmed against a real spectrum". Check one before you rely on
it for a regulatory assignment, and override it if your value differs.

The unflagged entries reproduce published nucleoside formulas, and the
whole engine reproduces the published molecular formulas of nusinersen
and inotersen exactly (`validate_reference()`).

### Sources

- Prakash T.P. et al., *Comparing in vitro and in vivo activity of
  2′-O-[2-(methylamino)-2-oxoethyl]- and 2′-O-methoxyethyl-modified
  antisense oligonucleotides*, J Med Chem 2008,
  [doi:10.1021/jm701537z](https://pubs.acs.org/doi/10.1021/jm701537z)
  (2'-O-NMA structure and properties)
- *Enhanced splicing modulation by NMA-modified antisense
  oligonucleotides*, Nucleic Acids Research 2026,
  [gkag484](https://academic.oup.com/nar/article/54/10/gkag484/8688739)
- Bio-Synthesis technology briefs:
  [2'-O-NMA](https://www.biosyn.com/nma-oligonucleotide-synthesis.aspx),
  [unlocked/flexible nucleic acids](https://www.biosyn.com/flexible-nucleic-acid-overview.aspx),
  [bridged nucleic acids](https://www.biosyn.com/bridged-nucleic-acid.aspx),
  [phosphonoacetate (PACE)](https://www.biosyn.com/pace-phosphonoacetate-oligo-modification.aspx),
  [phosphoramidate backbones](https://www.biosyn.com/phosphoramidate-pn-backbone-linkage-reengineering.aspx),
  [chiral phosphorothioates](https://www.biosyn.com/custom-chiral-oligo-synthesis.aspx),
  [morpholinos](https://www.biosyn.com/custom-morpholino-synthesis.aspx),
  [PNA](https://www.biosyn.com/custom-pna-synthesis.aspx)
