# From a chemical analysis file to a sequence you can run

A guide for people who have never typed an oligonucleotide sequence into
software before. It takes about fifteen minutes, and by the end you will
have the same sequence in three usable forms: three plain lines you can
type into OligoMet Profiler, a triplet string, and a FASTA file for
BioPharma Finder.

You do not need to understand oligonucleotide chemistry to do this. You
need to be able to read a table.

---

## 1. Why a sequence is not just letters

For ordinary DNA, a sequence is a string of letters: `ATGCATGC`. A
therapeutic oligonucleotide needs more than that, because two drugs can
have the identical letters and completely different masses.

Every position in a therapeutic oligo has **three** pieces of information:

| Piece | What it is | Everyday analogy |
|---|---|---|
| **Base** | the letter -- A, G, C, T, U | which character |
| **Sugar** | the ring the base is attached to (DNA? RNA? a modified sugar?) | which font |
| **Linkage** | the chemical bond joining this position to the next | the glue between characters |

The base tells you the letter. The sugar and the linkage are the
modifications that make it a drug: they are what stops the molecule being
chewed up in blood, and they are what changes its mass.

One detail that trips up everybody once, so learn it now:

> **Bases and sugars have one entry per position. Linkages have one
> fewer.** An 18-mer has 18 bases, 18 sugars, and **17** linkages, because
> linkages sit *between* positions. Eighteen fence posts, seventeen
> stretches of fence.

---

## 2. What to look for in the chemical analysis file

The document goes by several names -- certificate of analysis, chemical
analysis file, structure sheet, product specification. Whatever it is
called, you are looking for the part that describes the molecule
position by position. It usually looks like one of these three shapes.

**Shape A -- a table with one row per position.** The easiest case.

| Position | Base | Sugar | 3' linkage |
|---|---|---|---|
| 1 | T | 2'-MOE | PS |
| 2 | 5-Me-C | 2'-MOE | PS |
| 3 | A | 2'-MOE | PS |
| ... | ... | ... | ... |

**Shape B -- a written description plus a plain sequence.** Very common:

> 5'-TCACTTTCATAATGCTGG-3'
> All 2'-O-(2-methoxyethyl) ribonucleosides; all internucleoside linkages
> are phosphorothioate; all cytosines are 5-methylcytosine.

**Shape C -- a structural formula or a drawn diagram.** If this is all
you have, find the accompanying text. Do not try to count atoms off a
picture.

In every shape, you are hunting for exactly three answers:

1. **What are the letters, 5' to 3'?**
2. **What sugar is at each position?**
3. **What linkage is at each position?**

"5' to 3'" is just the reading direction, and it is the standard one. If
the document does not say, assume 5'→3'; if it gives a sequence labelled
3'→5', reverse it before continuing.

---

## 3. The code tables

Write each answer using these codes. Both OligoMet Profiler and (with the
caveat in section 7) BioPharma Finder read them.

### Bases

| Code | Means | Look for this wording |
|---|---|---|
| `A` | adenine | A, dA, Ade |
| `G` | guanine | G, dG, Gua |
| `C` | cytosine | C, dC, Cyt |
| `T` | thymine | T, dT, Thy |
| `U` | uracil | U, rU, Ura |
| `S` | **5-methylcytosine** | 5-Me-C, 5mC, m5C, "methylated cytosine" |
| `D` | 2,6-diaminopurine | 2,6-DAP |
| `I` | hypoxanthine (inosine) | I, dI, Ino |

`S` is the one people miss. Most 2'-MOE antisense drugs replace **every**
cytosine with 5-methylcytosine. If the document says so, every `C`
becomes `S`. Getting this wrong shifts the mass by 14 Da per cytosine,
which is enough to make nothing match.

`T` and 5-methyluracil are the same nucleobase, so use `T` for both.

### Sugars

| Code | Means | Look for this wording |
|---|---|---|
| `d` | 2'-deoxyribose (ordinary DNA) | DNA, deoxy, dA/dC/dG/dT |
| `r` | ribose (ordinary RNA) | RNA, ribo |
| `m` | 2'-O-methyl | 2'-OMe, 2'-O-Me, 2'-methoxy |
| `f` | 2'-fluoro | 2'-F |
| `e` | 2'-O-methoxyethyl (MOE) | 2'-MOE, 2'-O-(2-methoxyethyl) |
| `cEt` | constrained ethyl | cEt, (S)-cEt |
| `LNA` | locked nucleic acid | LNA, BNA |

### Linkages

| Code | Means | Look for this wording |
|---|---|---|
| `o` or `p` | phosphodiester -- the natural bond | PO, phosphodiester, "natural backbone" |
| `s` | phosphorothioate | PS, phosphorothioate, "sulfur-modified" |
| `mp` | methylphosphonate | MP, methylphosphonate |

`o` and `p` mean the same thing and weigh the same. `p` exists because
that is what BioPharma Finder writes; use whichever you prefer.

A **gapmer** -- an extremely common antisense design -- has modified
sugars at both ends and DNA in the middle. "5-10-5 MOE gapmer" means 5
MOE, then 10 DNA, then 5 MOE, so its sugar line is:

```
eeeeeddddddddddeeeee
```

---

## 4. Worked example

Take the 18-mer used as the example throughout the app. Its chemical
analysis text reads, in Shape B form:

> 5'-TCACTTTCATAATGCTGG-3'; uniform 2'-O-methoxyethyl (MOE) sugars;
> all internucleoside linkages phosphorothioate; all cytosines are
> 5-methylcytosine.

Work through the three questions.

**Bases.** Start from `TCACTTTCATAATGCTGG`. The text says every cytosine
is 5-methylcytosine, so replace every `C` with `S`:

```
TSASTTTSATAATGSTGG
```

Count them: 18.

**Sugars.** "Uniform 2'-MOE" means every position is `e`. Eighteen of
them:

```
eeeeeeeeeeeeeeeeee
```

**Linkages.** "All phosphorothioate" means every bond is `s`. An 18-mer
has 17 bonds:

```
sssssssssssssssss
```

Count that one carefully -- 17, not 18. (The app will tell you if you get
it wrong, and it accepts a single `s` meaning "all of them", which is
easier and safer.)

Those three lines are the answer. That is the whole job.

**Check yourself.** Run those three lines through the app and it reports
the molecular formula `C234H340N61O128P17S17` and an average mass of
7127.17 Da. This sequence is nusinersen, whose published formula and
average mass (7127.2 Da) are exactly that. If your own sequence's
computed mass lands within a fraction of a dalton of the value on your
chemical analysis file, you typed it correctly. If it is off by a
multiple of 14, suspect a `C` that should be an `S`; off by ~16 per
position, suspect a PO/PS mixup.

---

## 5. Entering it in OligoMet Profiler

Launch the app (`OligoMetProfiler::run_app()`), and use the **Manual
sequence entry** panel in the middle of the page:

| Field | What you type |
|---|---|
| Bases (5'→3') | `TSASTTTSATAATGSTGG` |
| Sugars | `eeeeeeeeeeeeeeeeee` |
| Linkages | `sssssssssssssssss` |

Click **Submit**. The app assembles the sequence, loads it into the
sequence box on the left, and shows you the formula and mass it
computed. Then click **Run Pipeline**.

Three conveniences worth knowing:

- **A single code means "all positions".** Typing `e` in Sugars and `s`
  in Linkages gives exactly the same result as typing them out. For a
  uniform drug this is much less error-prone.
- **Multi-character codes need separators.** `cEt` and `LNA` are more
  than one letter, so separate every code with a comma or dash:
  `MOE-MOE-d-d-d`, `A,G,S,T`.
- **Conjugates are set separately.** A GalNAc, cholesterol or 5'-phosphate
  goes in the 5'/3' conjugate dropdowns in the sidebar, not in these
  three fields.

The **Fill example** button loads the worked example above if you just
want to see it run.

---

## 6. Triplet notation (what the app builds for you)

The sequence box holds the same information in one line, called triplet
notation. For the example:

```
Te-sSe-sAe-sSe-sTe-sTe-sTe-sSe-sAe-sTe-sAe-sAe-sTe-sGe-sSe-sTe-sGe-sGe
```

Read one token at a time: `sSe` is **s** = the phosphorothioate bond
coming *into* this position, **S** = 5-methylcytosine, **e** = MOE
sugar. The first token has no linkage prefix, because nothing comes
before the 5' end -- which is why there are 18 tokens but only 17 `s`
prefixes.

You can type this form directly into the sequence box if you prefer. The
three-field panel exists so that you do not have to.

---

## 7. Making a BioPharma Finder input

BioPharma Finder annotates oligonucleotides in the same 5'-1'-2'
building-block style: a linkage prefix, then the base, then the sugar,
with `p` for phosphodiester and `s` for phosphorothioate. An unmodified
DNA `CAG` is `Cd-pAd-pGd` there; RNA is `Cr-pAr-pGr`; and with
phosphorothioate linkages it becomes `Cd-sAd-sGd`. That is the same
notation the app produces.

**To get a file:** with a sequence loaded, click **BPF FASTA** in the
manual entry panel. You get a `.fasta` file like this:

```
>nusinersen | 18-mer
Te-sSe-sAe-sSe-sTe-sTe-sTe-sSe-sAe-sTe-sAe-sAe-sTe-sGe-sSe-sTe-sGe-sGe
```

Import it through BioPharma Finder's Sequence Manager → **Import FASTA
File**. The file must have a `.fasta` extension and be under 1 MB;
BioPharma Finder rejects the file with an "invalid oligo building blocks
or an invalid format" message if it cannot resolve a code.

**Two things to verify on the BioPharma Finder side**, because they
depend on your installation rather than on this app:

1. **Sugar codes.** The exported file uses this package's codes
   (`d`, `r`, `m`, `f`, `e`, or `MOE`/`cEt`/`LNA` spelled out).
   BioPharma Finder resolves sugar codes against its own **Building
   Block and Variable Modification Editor**. If your site has defined a
   modification under a different code, either rename it in the exported
   file or add a matching building block in the editor. If an import is
   rejected, this is the first thing to check.
2. **Terminal conjugates.** Triplet notation has nowhere to put a GalNAc
   or a 5'-phosphate. The export writes them into the header line as a
   reminder, but you must set them as terminus modifications in
   BioPharma Finder yourself.

If you would rather type it in by hand, BioPharma Finder's manual
sequence entry takes the same three lines you prepared in section 4.

---

## 8. Mistakes people actually make

| Symptom | Likely cause |
|---|---|
| "Got 18 bases but 18 linkages" | You typed one linkage per position. There is one fewer -- or just type a single `s`. |
| Mass is ~14 Da per cytosine too low | Cytosines should be `S` (5-methylcytosine), not `C`. |
| Mass is ~16 Da per bond too low | Linkages should be `s` (PS), not `o`/`p` (PO). |
| Mass is roughly double what you expect | You entered a duplex. siRNA has two strands -- run the sense and antisense strands separately. |
| "Unknown sugar code 'E'" | Codes are case-sensitive. MOE is lowercase `e`; `LNA` and `cEt` are capitalized as shown. |
| "Unknown base code 'X'" | An unusual chemistry. Add it in the Custom Chemistry table in the sidebar with its formula, then use your own code. |
| Everything is right but nothing matches the MS data | Check the sequence direction. A sequence entered 3'→5' parses fine and gives a wrong answer silently. |

---

## 9. Where to go next

- [QUICKSTART.md](QUICKSTART.md) -- running the pipeline once you have a
  sequence
- [MODIFICATIONS.md](MODIFICATIONS.md) -- every modification the dictionary
  knows, its mass difference, and how to add one it does not
- [README.md](../../README.md) -- what the outputs contain, and the CLI drivers
- The vignette (`vignettes/OligoMetProfiler.Rmd` in the repository) -- the chemistry and the
  literature the mass calculations are grounded in

### Sources

- Thermo Fisher, [BioPharma Finder 5.2 Oligonucleotide Analysis User
  Guide -- Import a FASTA file containing an oligonucleotide
  sequence](https://docs.thermofisher.com/r/BioPharma-Finder-5.2-Oligonucleotide-Analysis-User-Guide/1288844299v2)
  (import procedure, `.fasta` extension, 1 MB limit, error behaviour)
- Thermo Fisher, [Oligonucleotide mapping using BioPharma Finder
  software, Application Note
  73789](https://documents.thermofisher.com/TFS-Assets/CMD/Application-Notes/an-73789-lc-ms-identification-mapping-oligonucleotides-an73789-en.pdf)
  (5'-1'-2' building-block notation; `Cd-pAd-pGd` / `Cr-pAr-pGr` /
  `Cd-sAd-sGd` examples)
