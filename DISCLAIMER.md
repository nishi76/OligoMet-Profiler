# Disclosure and Disclaimer

## Author

**Nishikant Wase, PhD** — author and developer
Research Scientist, Thermo Fisher Scientific
<nishikant.wase@gmail.com>
<https://github.com/nishi76/OligoMet-Profiler>

OligoMetProfiler was designed, written and is maintained by the author.
It is distributed under the [MIT licence](LICENSE.md), which includes the
warranty and liability terms restated below.

---

## FOR RESEARCH USE ONLY

**OligoMetProfiler is a research tool.** It is not a medical device, and
it is not intended or validated for:

- diagnostic use of any kind;
- clinical decision-making or patient care;
- quality control, batch release or stability testing;
- inclusion in a regulatory submission to any health authority;
- any application governed by GxP, CLIA, IVDR, or equivalent regulation.

Any such use falls outside the scope of this software and is undertaken
entirely at the user's own risk.

---

## Everything this software reports is a prediction

Metabolite libraries, molecular formulas, monoisotopic and average
masses, charge envelopes, isotope patterns, MS/MS fragment ions,
phosphorothioate oxidation series, acquisition target lists and spectral
libraries are **computed** from a chemistry dictionary and from
assumptions about how oligonucleotides are metabolised and how they
fragment.

They are not measurements. They are not evidence that any species exists
in any sample. A match between a predicted mass and an observed peak is
a hypothesis, not an identification. **Every assignment must be
confirmed experimentally before it is acted on.**

Specific limitations the user should be aware of:

- **Some formulas are best estimates.** Dictionary entries flagged
  `verify = TRUE` in `R/chemistry_dict.R` have not been checked against
  an independently measured mass. They are internally consistent, but
  confirming them against a known standard is the user's responsibility.
  See [inst/help/MODIFICATIONS.md](inst/help/MODIFICATIONS.md) section 7
  for the current list.
- **MS2 library intensities are placeholders.** There is no
  fragment-intensity model in this pipeline. Every MS2 peak in the
  exported MGF/MSP libraries is written at a flat value of 100. Match on
  *m/z*; intensity-weighted scoring against these libraries is
  meaningless.
- **Metabolite generation is a model, not a prediction of biology.** The
  truncation and endonuclease series reflect published degradation
  patterns for a few well-studied chemistries. Your compound, matrix,
  species or timepoint may behave differently, and real metabolites not
  in the library will simply be absent from it.
- **Isotope patterns depend on the engine available.** Results computed
  with the built-in convolution fallback differ in the last decimal
  places from those computed with `enviPat`.
- **Collision energies and instrument parameters are starting points.**
  The default HCD NCE and injection times in the acquisition exports are
  not validated instrument settings and must be optimised per method.

---

## No warranty and no liability

This software is provided **"as is"**, without warranty of any kind,
express or implied, including but not limited to the warranties of
merchantability, fitness for a particular purpose and non-infringement,
as set out in the MIT licence under which it is distributed.

**In no event shall the author be liable** for any claim, damages or
other liability, whether in an action of contract, tort or otherwise,
arising from, out of or in connection with this software or its use.
This includes, without limitation, any direct, indirect, incidental,
consequential or special loss; any loss or corruption of data; any
wasted instrument time, sample or materials; any missed or misassigned
metabolite; and any erroneous scientific, analytical or business
conclusion drawn from its output.

**The user is solely responsible** for determining that the software is
appropriate for their intended purpose, for validating it within their
own quality system where one applies, and for verifying every result it
produces before relying on it.

---

## Affiliation and conflict-of-interest disclosure

The author is employed as a **Research Scientist at Thermo Fisher
Scientific**. This is disclosed because the package interoperates with
Thermo Fisher products — it reads and writes BioPharma Finder sequence
notation and exports mass lists for the Orbitrap Exploris Method Editor.

**OligoMetProfiler is an independent personal project.** It is not a
Thermo Fisher Scientific product. It has not been supplied, reviewed,
supported, endorsed or approved by Thermo Fisher Scientific, and no
company has reviewed or approved its content. All interoperability with
third-party software and instruments is implemented from publicly
documented formats only; no proprietary or internal information is used.

## No endorsement, and trademarks

Nothing in this package implies endorsement or sponsorship by any
instrument vendor, software vendor or pharmaceutical company named in
it. Product, instrument, software and drug names are used solely for
identification and interoperability, and remain the trademarks or
registered trademarks of their respective owners.

The bundled reference drugs (nusinersen, inotersen, patisiran,
givosiran) are worked examples compiled from published literature and
public regulatory filings; they are not supplied, reviewed or approved
by their manufacturers.

All views and outputs are the author's own and do not represent Thermo
Fisher Scientific or any other employer or institution.

---

## Citation

If this software contributes to work you publish, please cite it as:

> Wase, N. *OligoMetProfiler: theoretical metabolite libraries, charge
> envelopes, and MS/MS fragment ions for therapeutic oligonucleotides.*
> R package. <https://github.com/nishi76/OligoMet-Profiler>

---

*This statement is also shown in the dashboard (the About section at the
bottom of the page), and travels inside the generated HTML/PDF report,
the Excel workbook and the exported spectral libraries. It is available
in R with `oligomet_about()`.*
