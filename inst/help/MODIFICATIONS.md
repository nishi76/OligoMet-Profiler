# Modification reference

> **FOR RESEARCH USE ONLY.** All outputs are computed predictions, not
> measurements — confirm every assignment experimentally. No warranty;
> see DISCLAIMER.md.

Every sugar and backbone chemistry in the dictionary, its mass effect,
and how to enter it. Not listed? Section 5 shows how to add it without
modifying the package. New to the notation? Read
[SEQUENCE_GUIDE.md](SEQUENCE_GUIDE.md) first.

## 1. How mass differences work

The formula engine is additive:

```
oligo = sum(bases) + sum(sugars) + sum(linkages) - (2n-1) H2O + conjugates
```

So a **sugar swap costs the residue difference per position** (DNA→MOE
= +74.0368 Da; ×18 positions = +1332.66 Da) and a **linkage swap costs
it per bond** (an *n*-mer has *n−1* bonds; PO→PS on an 18-mer =
17 × 15.9772 = +271.61 Da).

All values below are **monoisotopic** daltons, recomputed from the
dictionary by `tests/test_modifications.R` on every run. The app also
reports average masses.

## 2. Sugar modifications

Residue mass is the free sugar; a nucleoside is `base + sugar − H2O`.

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

Notes on the ones that cause trouble:

- **`NMA`** (2'-O-[2-(methylamino)-2-oxoethyl]) is an amide where MOE
  has an ether: **+12.9953 Da vs MOE** per position, or +233.9 Da across
  a uniform 18-mer. A next-generation chemistry; no approved drug yet.
- **`UNA`** is ribose + H2 (**+2.0157 Da**), a duplex destabilizer used
  at one or two siRNA seed positions. 2 Da is close to isotope spacing —
  check the full pattern, not the monoisotopic peak alone.
- **`GNA`** replaces the sugar with a glycerol unit: **−58.0055 Da vs
  RNA**, the largest per-position change here. Also a destabilizer.
- **`LNA`/`ENA`/`cEt`** (bridged/BNA family) lock the sugar pucker. LNA
  is ribose + CH2 − H2 = exactly +12.0000 Da; ENA is +14.0157 Da beyond
  LNA. All three are flagged `verify = TRUE`.
- **`FANA`** is **exactly isobaric** with 2'-F ribo — MS cannot
  distinguish them. The code exists only to keep the record correct.
- **`NP`** (N3'→P5' phosphoramidate) is modelled as a sugar change
  (**−0.9840 Da vs DNA**, O→NH), not a linkage change; keep the linkage
  as ordinary `o` or `s`.

## 3. Backbone (linkage) modifications

Each linkage is the phosphodiester bridge with its non-bridging `-OH`
swapped for something else. An *n*-mer has *n−1*.

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

- **`msp`** (mesyl phosphoramidate, µ) gives PS-like nuclease resistance
  without PS's protein binding, and is achiral at phosphorus (no Rp/Sp
  mixture). **+61.0164 Da per bond vs PS** — +1037.3 Da on a fully
  modified 18-mer, so partial conversion is very visible.
- **`pace`/`tpace`** replace the OH with `-CH2-COOH` (+42.0106 vs PO);
  thioPACE adds P=S for a further +15.9772. Used in siRNA.
- **`pgo` vs `tmg`** — phosphoryl guanidines differ by H2, **2.0157 Da
  per bond**, which compounds into a badly confusing discrepancy.
  Confirm which variant your synthesis used before assigning anything.
- **Alkyl phosphonates** (`mp` `prp` `ibu` `chx` `mop`) are neutral
  linkages that suppress siRNA off-target activity and steer RNase H1
  cleavage. `mp` is the awkward one: at −1.9793 Da vs PO, intact mass
  alone will not settle it — confirm with MS/MS.
- **Neutral backbones ionize poorly** in negative ESI, shifting the
  charge envelope below the `z_min`/`z_max` defaults. Widen downward if
  a PGO or alkyl-phosphonate oligo gives no matches.
- **Stereopure PS** (`u`) carries the same formula as `s` — MS cannot
  distinguish Rp from Sp; the code exists only to record it.

## 4. Entering a modification in the app

**Three-line entry.** Multi-character codes need a separator; a single
code is recycled to every position:

| Design | Bases | Sugars | Linkages |
|---|---|---|---|
| Uniform NMA, all msPA | `TSASTTTSATAATGSTGG` | `NMA` | `msp` |
| 5-10-5 NMA gapmer, all PS | `TSASTTTSATAATGSTGG` | `NMA,NMA,NMA,NMA,NMA,d,d,d,d,d,d,d,d,NMA,NMA,NMA,NMA,NMA` | `s` |
| LNA/DNA mixmer | `AGSTAGST` | `LNA-d-LNA-d-LNA-d-LNA-d` | `s` |

**Triplet notation.** `[linkage][base][sugar]` per token; multi-character
linkage codes match longest-first, so `mspANMA` reads as msPA + A + NMA:

```
TNMA-mspSNMA-mspANMA-mspSNMA-mspTNMA
```

**From R:**

```r
spec <- parse_three_line(
  bases    = "TSASTTTSATAATGSTGG",
  sugars   = c(rep("NMA", 5), rep("d", 8), rep("NMA", 5)),
  linkages = c(rep("msp", 5), rep("s", 7), rep("msp", 5)))
```

## 5. Adding a modification the dictionary lacks

No need to edit the package. In the app use the **Custom Chemistry**
table (one row per code: Code, Formula, Name, Type — `base`/`sugar`/
`linkage`/`conjugate` — and Attach for conjugates). From R:

```r
dict <- build_dictionary(overrides = list(
  tcDNA = list(formula = "C8H14O4", name = "tricyclo-DNA", kind = "sugar")))
spec <- parse_three_line("AGST", "tcDNA", "s", dict = dict)
```

**Working out the residue formula** — what the molecule contributes
before condensation; the engine subtracts the water:

- **2'-O-R sugar** = `C5H9O5 + R`. NMA: R = C3H6NO → C8H15NO6.
- **Bridged sugar** = parent + bridge − 2 H. LNA = ribose + CH2 − H2.
- **Linkage** = `PO2 + X`, X being the substituent. PACE: X = C2H3O2 →
  C2H3O4P.
- **Base** = the free nucleobase, not the nucleoside.

Then check it: `base + sugar − H2O` should reproduce a published
nucleoside formula. That one check catches nearly every mistake.

```r
f <- add_formulas(STANDARD_DICT[["A"]]$formula, dict[["tcDNA"]]$formula)
format_formula(f - c(H = 2, O = 1))
```

`tests/test_modifications.R` does this for every built-in code; copy a
row from its `nucleoside_anchors` list to add your own.

## 6. Deliberately absent

These do not fit the additive residue model, and a wrong number is worse
than a missing one:

- **Morpholino (PMO/PPMO)** — replaces both sugar and linkage; can be
  expressed as custom entries, but take the subunit formula from your
  own reference.
- **PNA** — no sugar and no phosphate; the decomposition doesn't apply.
- **Boranophosphate** — needs boron added to `.ELEMENTS` with its
  isotope abundances (10B/11B ≈ 20/80 visibly changes patterns), not
  just a residue formula.
- **Tricyclo-DNA** — a straightforward custom sugar; left out only so
  the formula comes from a source you trust. Add it per section 5.

## 7. Which formulas are provisional

Codes flagged `verify = TRUE` in `chemistry_dict.R` are best estimates
not checked against an independently measured mass: `cEt`, `LNA`, `ENA`,
`AP`, `msp`, `pace`, `tpace`, `pgo`, `mp`, and the BioPharma Finder
terminal-modification conjugates. They are internally consistent and
verified by `tests/test_modifications.R`, but that is not the same as
confirmed against a real spectrum. Check one before relying on it for a
formal assignment, and override it if your value differs.

Unflagged entries reproduce published nucleoside formulas, and the
engine reproduces the published molecular formulas of nusinersen and
inotersen exactly (`validate_reference()`).

### Sources

- Prakash T.P. et al., *Comparing in vitro and in vivo activity of
  2′-O-[2-(methylamino)-2-oxoethyl]- and 2′-O-methoxyethyl-modified
  antisense oligonucleotides*, J Med Chem 2008,
  [doi:10.1021/jm701537z](https://pubs.acs.org/doi/10.1021/jm701537z)
- *Enhanced splicing modulation by NMA-modified antisense
  oligonucleotides*, Nucleic Acids Research 2026,
  [gkag484](https://academic.oup.com/nar/article/54/10/gkag484/8688739)
- Bio-Synthesis technology briefs:
  [2'-O-NMA](https://www.biosyn.com/nma-oligonucleotide-synthesis.aspx),
  [flexible nucleic acids](https://www.biosyn.com/flexible-nucleic-acid-overview.aspx),
  [bridged nucleic acids](https://www.biosyn.com/bridged-nucleic-acid.aspx),
  [PACE](https://www.biosyn.com/pace-phosphonoacetate-oligo-modification.aspx),
  [phosphoramidate backbones](https://www.biosyn.com/phosphoramidate-pn-backbone-linkage-reengineering.aspx),
  [chiral phosphorothioates](https://www.biosyn.com/custom-chiral-oligo-synthesis.aspx),
  [morpholinos](https://www.biosyn.com/custom-morpholino-synthesis.aspx),
  [PNA](https://www.biosyn.com/custom-pna-synthesis.aspx)

---

*OligoMetProfiler — Nishikant Wase, PhD. MIT licence. Research use
only; see DISCLAIMER.md.*
