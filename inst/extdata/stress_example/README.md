# Forced-degradation stress study (example data)

A simulated forced-degradation dataset for trying every analysis tab without
any raw files, a generated library, or the Python backend. Load it from the
app: **Batch Processing - "Try an example dataset" - Forced-degradation stress
study - Load**. It imports as a pre-processed data matrix.

The analyte is a generic placeholder (`example_gapmer`, a 20-mer 5-10-5
MOE/DNA phosphorothioate gapmer). Any real analyte name comes from your own
input, never from this file.

## Design

- Three stress arms (the `group` column): `acid`, `base`, `thermal`.
- Six timepoints (`timepoint`, hours): 0, 6, 24, 48, 72, 168.
- 27 species: the full-length parent plus its exonuclease and gap
  endonuclease products.
- `intensity` is a simulated LC-MS peak area, percent of the t=0 parent. The
  parent's oxidation (PS to PO) series is summed into total parent signal.

Time 0 is the dosing point (full intact signal), so it is the 100% reference,
not a predose control.

## Files

- `data_matrix_long.csv` - the importable long matrix: `met_id`, `met_name`,
  `kind`, `sample`, `intensity`, `group`, `timepoint`, `sample_type`,
  `concentration`.
- `sample_info.csv` - the sample table (sample, group, timepoint), the same
  information embedded in the long matrix.

The numbers are simulated for demonstration and are not measured data.
