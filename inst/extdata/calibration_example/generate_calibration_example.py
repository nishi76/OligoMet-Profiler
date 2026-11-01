#!/usr/bin/env python3
"""Generates the bundled calibration-curve example: 10 small synthetic mzML
files -- a 5-level standard curve, n=2 replicate injections per level --
for the SAME inotersen reference sequence as inst/extdata/batch_example/,
plus sample_meta.csv. Run from the repo root:

    python3 inst/extdata/calibration_example/generate_calibration_example.py

Regenerate only if you want different synthetic data -- the checked-in
files are already reproducible from this script (fixed seeds).

Concentration levels (ng/mL): 1, 5, 25, 100, 500 -- a ~500x dynamic range,
geometric-ish spacing typical of a real LC-MS/MS bioanalytical curve.
Each standard's parent charge envelope (z=5,6,7,8, same theoretical m/z as
batch_example's inotersen) scales LINEARLY with concentration
(signal = SLOPE * concentration, no intercept), with the same per-scan
Gaussian-ish elution shape, per-charge-state weighting, and multiplicative
CV/ppm noise convention as generate_example.py -- constant RELATIVE noise
at every level, which is exactly the heteroscedastic case 1/x^2-weighted
regression (fit_calibration_curve()'s recommended default) is meant for.

Deliberately minimal compared to batch_example: no degradants, no
contaminant trace, no MS2 scan -- these are pure calibration standards
(a clean spike of parent reference material), not a degraded biological
sample, so there is nothing else to model here. Sample Type is "standard"
for every file; Group/Timepoint are blank (calibration standards are not
part of a group/timepoint study design -- see .is_study_sample() in
R/chemistry_dict.R, which excludes them from every biological comparison
for exactly this reason).
"""
from __future__ import annotations

import base64
import json
import os
import random
import struct
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))

# (level, concentration in ng/mL)
LEVELS = [(1, 1.0), (2, 5.0), (3, 25.0), (4, 100.0), (5, 500.0)]
N_REPLICATES = 2
SLOPE = 50000.0  # arbitrary intensity units per ng/mL -- see module docstring


def _b64(vals):
    raw = struct.pack(f"<{len(vals)}d", *vals)
    return base64.b64encode(zlib.compress(raw)).decode()


def _binarr(kind, vals):
    acc = "MS:1000514" if kind == "mz" else "MS:1000515"
    name = "m/z array" if kind == "mz" else "intensity array"
    b64 = _b64(vals)
    return (
        f'<binaryDataArray encodedLength="{len(b64)}">'
        f'<cvParam accession="{acc}" name="{name}"/>'
        f'<cvParam accession="MS:1000523" name="64-bit float"/>'
        f'<cvParam accession="MS:1000574" name="zlib compression"/>'
        f'<binary>{b64}</binary></binaryDataArray>'
    )


def _ms1(idx, rt, mzs, ints):
    return (
        f'<spectrum id="scan={idx}" index="{idx - 1}" defaultArrayLength="{len(mzs)}">'
        f'<cvParam accession="MS:1000511" name="ms level" value="1"/>'
        f'<scanList count="1"><scan><cvParam accession="MS:1000016" name="scan start time" '
        f'value="{rt}" unitAccession="UO:0000031" unitName="minute"/></scan></scanList>'
        f'<binaryDataArrayList count="2">{_binarr("mz", mzs)}{_binarr("int", ints)}</binaryDataArrayList>'
        f'</spectrum>'
    )


def build_file(path, concentration, rep, theo, seed):
    rng = random.Random(seed)
    z_mz = {int(k[1:]): v for k, v in theo["mz"].items()}
    base_intensity = SLOPE * concentration
    profile = [0.16, 0.40, 0.72, 1.00, 0.72, 0.40, 0.16]  # Gaussian-ish elution shape
    rts = [4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3]
    per_charge_weight = {5: 0.6, 6: 1.0, 7: 0.9, 8: 0.5}

    spectra, idx = [], 1
    for rt, shape in zip(rts, profile):
        mzs, ints = [], []
        for z, mz in sorted(z_mz.items()):
            ppm_noise = rng.uniform(-2.5, 2.5)
            mzs.append(mz * (1 + ppm_noise / 1e6))
            cv_noise = rng.uniform(0.9, 1.1)  # constant RELATIVE noise at every level
            ints.append(base_intensity * shape * per_charge_weight[z] * cv_noise)
        order = sorted(range(len(mzs)), key=lambda i: mzs[i])
        spectra.append(_ms1(idx, rt, [mzs[i] for i in order], [ints[i] for i in order]))
        idx += 1

    xml = (
        '<?xml version="1.0"?>\n<mzML xmlns="http://psi.hupo.org/ms/mzml">\n<run>\n'
        f'<spectrumList count="{idx - 1}">\n' + "\n".join(spectra) + "\n</spectrumList>\n</run>\n</mzML>\n"
    )
    with open(path, "w") as f:
        f.write(xml)


def main():
    with open(os.path.join(HERE, "theoretical_values.json")) as f:
        theo = json.load(f)

    rows = ["sample,group,timepoint,sample_type,concentration"]
    seed = 2000
    for level, conc in LEVELS:
        for rep in range(1, N_REPLICATES + 1):
            sample = f"std_L{level}_r{rep}"
            path = os.path.join(HERE, f"{sample}.mzML")
            build_file(path, conc, rep, theo, seed)
            rows.append(f"{sample},,,standard,{conc}")
            seed += 1

    with open(os.path.join(HERE, "sample_meta.csv"), "w") as f:
        f.write("\n".join(rows) + "\n")

    print(f"Wrote {len(LEVELS) * N_REPLICATES} mzML files + sample_meta.csv to {HERE}")


if __name__ == "__main__":
    main()
