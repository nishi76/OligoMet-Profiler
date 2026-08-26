# Disclosure and Disclaimer

**Nishikant Wase, PhD** — author and developer
<nishikant.wase@gmail.com> · <https://github.com/nishi76/OligoMet-Profiler>
Distributed under the [MIT licence](LICENSE.md).

## For research use only

OligoMetProfiler is a research tool. It is not a medical device and is
not intended or validated for diagnostic use, clinical decision-making,
quality control or release testing, or inclusion in any regulatory
submission. Any such use is outside its scope and at the user's own risk.

## Everything it reports is a prediction

All outputs metabolite libraries, formulas, masses, charge envelopes,
isotope patterns, fragment ions, target lists, and spectral libraries —
are computed from a chemistry dictionary and modelling assumptions. They
are not measurements: a match between a predicted mass and an observed
peak is a hypothesis, and **every assignment must be confirmed
experimentally**. In particular:

- Dictionary entries flagged `verify = TRUE` are best estimates not yet
  checked against a measured mass (see
  [MODIFICATIONS.md](inst/help/MODIFICATIONS.md) section 7).
- MS2 library intensities are flat placeholders — match on *m/z* only.
- The metabolite series model published degradation patterns; your
  compound or matrix may behave differently.
- Default collision energies and instrument parameters are starting
  points, not validated settings.

## No warranty, no liability

This software is provided **"as is"**, without warranty of any kind,
express or implied, under the MIT licence. In no event shall the author
be liable for any claim, damages, or other liability arising from this
software or its use — including loss of data, wasted instrument time or
materials, or erroneous conclusions drawn from its output. The user is
solely responsible for verifying fitness for purpose and validating
every result.

## Affiliation disclosure

The author is a Research Scientist at Thermo Fisher Scientific,
disclosed because the package interoperates with Thermo Fisher products
(BioPharma Finder notation, Orbitrap Exploris mass lists).
**OligoMetProfiler is an independent personal project**  not a Thermo
Fisher Scientific product, and not supplied, reviewed, supported,
endorsed, or approved by Thermo Fisher Scientific or any other company.
All interoperability is implemented from publicly documented formats
only. Product, instrument, and drug names are used solely for
identification and remain their owners' trademarks; the bundled
reference drugs are worked examples from published literature and public
regulatory filings. All views and outputs are the author's own.

## Citation

> Wase, N. *OligoMetProfiler: theoretical metabolite libraries, charge
> envelopes, and MS/MS fragment ions for therapeutic oligonucleotides.*
> R package. <https://github.com/nishi76/OligoMet-Profiler>

*Also shown in the dashboard's About section, in generated reports and
workbooks, and in R via `oligomet_about()`.*
