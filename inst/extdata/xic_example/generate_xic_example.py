#!/usr/bin/env python3
"""Generates the bundled targeted summed-XIC example: synthetic centroided
mzML files for the inotersen reference sequence with FULL isotope envelopes
and a realistic charge-state distribution (CSD), so that summing charge
states x isotopes can be compared against single-channel quantitation.

    python3 inst/extdata/xic_example/generate_xic_example.py   (from the repo root)

theoretical_clusters.json (isotope cluster masses/abundances, from
make_theoretical_clusters.R) supplies the chemistry; this script only adds
chromatography, CSD, and noise. Fixed seeds, so the output is reproducible.

Design (sample_meta.csv):
  * 5-level calibration, n=2 (1, 5, 25, 100, 500 ng/mL), sample_type standard
  * QCs at 3, 50, 400 ng/mL
  * one zero blank (RB_01): matrix + internal standard, no analyte
  * four unknowns (U1..U4 at 20, 80, 200, 300 ng/mL) carrying 10% of their
    parent amount as the 3' N-1 metabolite. U3 and U4 have a SHIFTED CSD
    (z=3 up, z=7/8 down), the way a different matrix or ion-pairing
    condition shifts it -- the case where a fixed, broad charge-state set
    matters.
CSD modelled on a real 7 kDa PS-oligo full scan: z3..z9 relative apex
0.42, 0.55, 0.64, 0.88, 1.00, 0.58, 0.13.

Internal standard: an analog 16-mer gapmer (IS01 in theoretical_clusters.json,
second record of sequences.fasta) spiked at a constant amount into every
sample except U2, where it is under-spiked to 35% (an IS pipetting error the
IS response monitor should flag). Every injection also gets a random volume
factor (0.8-1.2) that scales analyte and IS alike -- the variability an IS
corrects for.
"""
from __future__ import annotations

import base64
import json
import math
import os
import random
import struct
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
PROTON = 1.007276466
RESPONSE = 20000.0            # summed apex intensity (all ions) per ng/mL
APEX_RT, SIGMA = 4.83, 0.03   # min
RTS = [round(4.60 + 0.01 * i, 3) for i in range(46)]

CSD = {3: 0.42, 4: 0.55, 5: 0.64, 6: 0.88, 7: 1.00, 8: 0.58, 9: 0.13}
CSD_SHIFTED = {3: 1.00, 4: 0.66, 5: 0.68, 6: 0.98, 7: 0.80, 8: 0.40, 9: 0.10}


def _b64(vals, fmt):
    return base64.b64encode(zlib.compress(struct.pack(f"<{len(vals)}{fmt}", *vals))).decode()


def _binarr(kind, vals):
    acc = "MS:1000514" if kind == "mz" else "MS:1000515"
    name = "m/z array" if kind == "mz" else "intensity array"
    # m/z as 64-bit (ppm accuracy), intensity as 32-bit (keeps files small)
    fmt, bits, facc = ("d", "64", "MS:1000523") if kind == "mz" else ("f", "32", "MS:1000521")
    b64 = _b64(vals, fmt)
    return (f'<binaryDataArray encodedLength="{len(b64)}"><cvParam accession="{acc}" name="{name}"/>'
            f'<cvParam accession="{facc}" name="{bits}-bit float"/>'
            f'<cvParam accession="MS:1000574" name="zlib compression"/><binary>{b64}</binary></binaryDataArray>')


def _ms1(idx, rt, mzs, ints):
    return (f'<spectrum id="scan={idx}" index="{idx - 1}" defaultArrayLength="{len(mzs)}">'
            f'<cvParam accession="MS:1000511" name="ms level" value="1"/>'
            f'<cvParam accession="MS:1000127" name="centroid spectrum"/>'
            f'<scanList count="1"><scan><cvParam accession="MS:1000016" name="scan start time" '
            f'value="{rt}" unitAccession="UO:0000031" unitName="minute"/></scan></scanList>'
            f'<binaryDataArrayList count="2">{_binarr("mz", mzs)}{_binarr("int", ints)}</binaryDataArrayList>'
            f'</spectrum>')


def ions_for(species, amount, csd):
    """(mz, apex intensity) for every charge x isotope cluster of one species."""
    csd_tot = sum(csd.values())
    out = []
    for z, wz in csd.items():
        for m, ab in zip(species["mass"], species["abundance"]):
            out.append(((m - z * PROTON) / z, RESPONSE * amount * (wz / csd_tot) * ab))
    return out


IS_AMOUNT = 150.0


def build_file(path, parent_amt, n1_amt, csd, seed, is_frac=1.0):
    rng = random.Random(seed)
    parent, n1, istd = THEO[0], THEO[1], THEO[2]
    vol = rng.uniform(0.8, 1.2)  # injection-volume factor, common to every species
    ions = ions_for(parent, parent_amt * vol, csd) + (ions_for(n1, n1_amt * vol, csd) if n1_amt > 0 else [])
    ions += ions_for(istd, IS_AMOUNT * is_frac * vol, csd)
    # fixed chemical-background ions for this file (random m/z, flat in time)
    bg = [(rng.uniform(600, 2450), rng.lognormvariate(math.log(2500), 0.6)) for _ in range(25)]
    spectra = []
    for i, rt in enumerate(RTS, start=1):
        g = math.exp(-0.5 * ((rt - APEX_RT) / SIGMA) ** 2)
        mzs, ints = [], []
        for mz, h in ions:
            val = h * g * rng.uniform(0.93, 1.07) + abs(rng.gauss(0, 120))
            if val > 150:  # centroiding threshold
                mzs.append(mz * (1 + rng.uniform(-2.0, 2.0) / 1e6))
                ints.append(val)
        for mz, h in bg:
            mzs.append(mz)
            ints.append(h * rng.uniform(0.8, 1.2))
        order = sorted(range(len(mzs)), key=lambda k: mzs[k])
        spectra.append(_ms1(i, rt, [mzs[k] for k in order], [ints[k] for k in order]))
    xml = ('<?xml version="1.0"?>\n<mzML xmlns="http://psi.hupo.org/ms/mzml">\n<run>\n'
           f'<spectrumList count="{len(spectra)}">\n' + "\n".join(spectra) + "\n</spectrumList>\n</run>\n</mzML>\n")
    with open(path, "w") as f:
        f.write(xml)


def main():
    rows = ["sample,group,timepoint,sample_type,concentration"]
    seed = 5000
    plan = []
    for lvl, conc in enumerate([1, 5, 25, 100, 500], start=1):
        for rep in (1, 2):
            plan.append((f"std_L{lvl}_r{rep}", conc, 0.0, CSD, "standard", conc))
    for name, conc in (("QC_low", 3), ("QC_mid", 50), ("QC_high", 400)):
        plan.append((name, conc, 0.0, CSD, "quality_control", conc))
    plan.append(("RB_01", 0.0, 0.0, CSD, "reagent_blank", ""))
    for name, conc, csd in (("U1", 20, CSD), ("U2", 80, CSD), ("U3", 200, CSD_SHIFTED),
                            ("U4", 300, CSD_SHIFTED)):
        plan.append((name, conc, 0.1 * conc, csd, "unknown", ""))
    for name, amt, n1, csd, stype, conc in plan:
        build_file(os.path.join(HERE, f"{name}.mzML"), amt, n1, csd, seed,
                   is_frac=0.35 if name == "U2" else 1.0)
        rows.append(f"{name},,,{stype},{conc}")
        seed += 1
    with open(os.path.join(HERE, "sample_meta.csv"), "w") as f:
        f.write("\n".join(rows) + "\n")
    # true concentrations of the unknowns, for checking back-calculation
    with open(os.path.join(HERE, "unknown_truth.csv"), "w") as f:
        f.write("sample,true_concentration,csd,is_note\nU1,20,normal,\nU2,80,normal,IS under-spiked to 35%\n"
                "U3,200,shifted,\nU4,300,shifted,\n")


with open(os.path.join(HERE, "theoretical_clusters.json")) as _f:
    THEO = json.load(_f)

if __name__ == "__main__":
    main()
