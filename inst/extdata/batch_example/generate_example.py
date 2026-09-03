#!/usr/bin/env python3
"""Generates the bundled batch-processing example: 6 small synthetic mzML
files (3 "control" + 3 "treated" replicates) for the inotersen reference
sequence, plus sample_meta.csv. Run from the repo root:

    python3 inst/extdata/batch_example/generate_example.py

Regenerate only if you want different synthetic data -- the checked-in
files are already reproducible from this script (fixed seeds).

Each file has:
  - a 7-scan MS1 charge envelope (z=5,6,7,8) for the inotersen parent,
    with realistic per-replicate mass (ppm) and intensity (CV) noise,
    "treated" replicates at ~3x "control" intensity;
  - a 7-scan MS1 charge envelope (z=5,6,7,8) for each of two truncated
    degradant forms (3' N-1 and 5' N-1, from theoretical_values.json's
    "degradants" list), at a lower level than the parent in both groups
    but noticeably higher under "treated" (~2.3x "control") -- gives the
    degradation-percentage calculation (R/degradation.R) something
    non-trivial and non-zero to compute, and gives the kind-level group
    comparison (R/statistics.R's build_kind_abundance_matrix()) a real,
    detectable control-vs-treated difference to find. Every charge
    state's m/z sits tens of Da clear of the parent's and the
    contaminant's, so ROI extraction keeps all three as separate traces;
  - a 7-scan "contaminant" trace at an unrelated m/z that matches no
    theoretical metabolite at any charge/adduct -- persists across
    enough scans to form a real ROI, so it lands in "unidentified
    peaks" downstream instead of being filtered out as noise;
  - one MS2 scan on the z=6 precursor with real McLuckey fragment m/z
    (from R/fragments.R, see theoretical_values.json) plus the two PS
    diagnostic ions, so MS2 confirmation has something real to score.
"""
from __future__ import annotations

import base64
import json
import os
import random
import struct
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
CONTAMINANT_MZ = 743.61  # arbitrary, matches no theoretical charge state


def _b64(vals):
    raw = struct.pack(f"<{len(vals)}d", *vals)
    return base64.b64encode(zlib.compress(raw)).decode()


def _binarr(kind, vals):
    acc = "MS:1000514" if kind == "mz" else "MS:1000515"
    name = "m/z array" if kind == "mz" else "intensity array"
    b64 = _b64(vals)
    # encodedLength (length of the base64 string, in bytes) is a REQUIRED
    # attribute of binaryDataArray per the mzML spec. Our own xml2-based
    # parser (parse_mzml() in R/ms_matching.R) never checks it, but mzR's
    # stricter parser does and rejects a binaryDataArray missing it -- so it
    # must be set for these fixtures to be readable by mzR/Spectra as well
    # as the built-in fallback parser.
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


def _ms2(idx, rt, prec_mz, prec_z, mzs, ints):
    return (
        f'<spectrum id="scan={idx}" index="{idx - 1}" defaultArrayLength="{len(mzs)}">'
        f'<cvParam accession="MS:1000511" name="ms level" value="2"/>'
        f'<scanList count="1"><scan><cvParam accession="MS:1000016" name="scan start time" '
        f'value="{rt}" unitAccession="UO:0000031" unitName="minute"/></scan></scanList>'
        f'<precursorList count="1"><precursor spectrumRef="scan=1">'
        f'<selectedIonList count="1"><selectedIon>'
        f'<cvParam accession="MS:1000744" name="selected ion m/z" value="{prec_mz}"/>'
        f'<cvParam accession="MS:1000041" name="charge state" value="{prec_z}"/>'
        f'</selectedIon></selectedIonList></precursor></precursorList>'
        f'<binaryDataArrayList count="2">{_binarr("mz", mzs)}{_binarr("int", ints)}</binaryDataArrayList>'
        f'</spectrum>'
    )


def build_file(path, group, replicate, theo, seed):
    rng = random.Random(seed)
    z_mz = {int(k[1:]): v for k, v in theo["mz"].items()}
    base_intensity = 90000 if group == "control" else 270000  # ~3x fold change
    profile = [0.16, 0.40, 0.72, 1.00, 0.72, 0.40, 0.16]  # Gaussian-ish elution shape
    rts = [4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3]

    per_charge_weight = {5: 0.6, 6: 1.0, 7: 0.9, 8: 0.5}
    # Degradant level relative to that group's OWN parent base_intensity --
    # low in both groups (so the parent still dominates), clearly higher
    # under "treated". At the dominant charge states (z=6,7; weight >= 0.9)
    # this clears the ROI builder's default min_intensity=1e4 floor in both
    # groups; the weaker z=5,8 charge states may not, same as the parent.
    degradant_frac = {"control": 0.15, "treated": 0.35}[group]
    degradants = theo.get("degradants", [])

    spectra, idx = [], 1
    for rt, shape in zip(rts, profile):
        mzs, ints = [], []
        for z, mz in sorted(z_mz.items()):
            ppm_noise = rng.uniform(-2.5, 2.5)
            mzs.append(mz * (1 + ppm_noise / 1e6))
            cv_noise = rng.uniform(0.9, 1.1)
            ints.append(base_intensity * shape * per_charge_weight[z] * cv_noise)
        for degr in degradants:
            degr_z_mz = {int(k[1:]): v for k, v in degr["mz"].items()}
            for z, mz in sorted(degr_z_mz.items()):
                ppm_noise = rng.uniform(-2.5, 2.5)
                mzs.append(mz * (1 + ppm_noise / 1e6))
                cv_noise = rng.uniform(0.9, 1.1)
                ints.append(base_intensity * degradant_frac * shape * per_charge_weight[z] * cv_noise)
        # a contaminant trace present in every scan (own elution shape),
        # unrelated to any theoretical metabolite -- exercises the
        # "unidentified peaks" (retained, not discarded) path downstream
        ppm_noise = rng.uniform(-3, 3)
        mzs.append(CONTAMINANT_MZ * (1 + ppm_noise / 1e6))
        ints.append(4e4 * (0.5 + 0.5 * shape) * rng.uniform(0.9, 1.1))
        order = sorted(range(len(mzs)), key=lambda i: mzs[i])
        spectra.append(_ms1(idx, rt, [mzs[i] for i in order], [ints[i] for i in order]))
        idx += 1

    frag_mzs = list(theo["frag_mz"]) + list(theo["ps_diag"])
    frag_ints = [rng.uniform(2000, 9000) for _ in theo["frag_mz"]] + \
                [rng.uniform(6000, 9000) for _ in theo["ps_diag"]]
    prec_mz = z_mz[6] * (1 + rng.uniform(-2, 2) / 1e6)
    spectra.append(_ms2(idx, 5.0, prec_mz, 6, frag_mzs, frag_ints))
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

    rows = ["sample,group"]
    seed = 1000
    for group in ("control", "treated"):
        for rep in (1, 2, 3):
            sample = f"{'ctrl' if group == 'control' else 'treat'}_{rep}"
            path = os.path.join(HERE, f"{sample}.mzML")
            build_file(path, group, rep, theo, seed)
            rows.append(f"{sample},{group}")
            seed += 1

    with open(os.path.join(HERE, "sample_meta.csv"), "w") as f:
        f.write("\n".join(rows) + "\n")

    print(f"Wrote 6 mzML files + sample_meta.csv to {HERE}")


if __name__ == "__main__":
    main()
