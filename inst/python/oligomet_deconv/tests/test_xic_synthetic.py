"""Console PASS/FAIL check for targeted summed-XIC quantitation (xic.py)
against a hand-built mzML with known Gaussian peaks, a constant baseline, and
an unrelated later-eluting peak on the same ions.

Run with:  python -m oligomet_deconv.tests.test_xic_synthetic
(run from inst/python/, or with inst/python/ on PYTHONPATH)
"""

from __future__ import annotations

import math
import os
import sys
import tempfile

import numpy as np

from oligomet_deconv.tests.test_deconv_synthetic import _ms1_spectrum, check, _FAILED

PROTON = 1.007276466
M = 7051.2
CHARGES = (5, 6, 7)
ISO = (2, 3, 4, 5)           # isotope offsets (Da, x 1.00335)
Z_W = {5: 0.6, 6: 0.9, 7: 1.0}
ISO_W = {2: 0.7, 3: 0.95, 4: 1.0, 5: 0.9}
APEX, SIGMA, HEIGHT = 4.83, 0.03, 1.0e6
BASELINE = 2000.0            # constant offset on every ion, every scan
INTERFERENT_RT = 5.60        # same m/z, outside the search window


def ion_mz(z, k):
    return (M + k * 1.003355 - z * PROTON) / z


def build(path):
    rts = np.round(np.arange(4.40, 5.90, 0.005), 4)
    rng = np.random.default_rng(11)
    spectra = []
    for i, rt in enumerate(rts, start=1):
        mzs, ints = [], []
        g = math.exp(-0.5 * ((rt - APEX) / SIGMA) ** 2)
        gi = math.exp(-0.5 * ((rt - INTERFERENT_RT) / SIGMA) ** 2)
        for z in CHARGES:
            for k in ISO:
                mzs.append(ion_mz(z, k))
                ints.append(BASELINE + HEIGHT * Z_W[z] * ISO_W[k] * g + 3 * HEIGHT * gi
                            + rng.normal(0, 200.0))
        # an unrelated ion 50 ppm away must not leak into the +/-10 ppm XIC
        mzs.append(ion_mz(7, 4) * (1 + 50e-6))
        ints.append(5 * HEIGHT * g)
        order = np.argsort(mzs)
        spectra.append(_ms1_spectrum(i, float(rt), list(np.array(mzs)[order]), list(np.array(ints)[order])))
    xml = ('<?xml version="1.0"?>\n<mzML xmlns="http://psi.hupo.org/ms/mzml">\n<run>\n'
           f'<spectrumList count="{len(spectra)}">\n' + "\n".join(spectra)
           + "\n</spectrumList>\n</run>\n</mzML>\n")
    with open(path, "w") as f:
        f.write(xml)


def main():
    print("=== targeted summed-XIC synthetic test ===")
    from oligomet_deconv.xic import XICParams, process_file_xic
    from oligomet_deconv.xic_cli import run_xic_batch

    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "std_1.mzML")
    build(path)
    targets = [{"target_id": "T1", "ion_id": f"z{z}_i{k}", "mz": ion_mz(z, k), "rt_expected": 4.85}
               for z in CHARGES for k in ISO]
    ions, traces, meta = process_file_xic(path, targets, XICParams(ppm=10))
    area = {r["ion_id"]: r["area"] for r in ions}
    gauss = HEIGHT * SIGMA * math.sqrt(2 * math.pi)
    expected = {f"z{z}_i{k}": gauss * Z_W[z] * ISO_W[k] for z in CHARGES for k in ISO}

    total = sum(area.values())
    total_exp = sum(expected.values())
    check(f"summed area within 3% of analytic ({total:.4g} vs {total_exp:.4g})",
          abs(total / total_exp - 1) < 0.03)
    worst = max(abs(area[i] / expected[i] - 1) for i in expected)
    check(f"every ion within 3% of analytic (worst {100 * worst:.2f}%)", worst < 0.03)
    windows = {(r["rt_start"], r["rt_end"]) for r in ions}
    check("all ions share one integration window", len(windows) == 1)
    r0 = ions[0]
    check(f"apex found at {r0['rt_apex']:.3f} min (expected {APEX})", abs(r0["rt_apex"] - APEX) < 0.01)
    check("later interferent at 5.60 min not integrated", r0["rt_end"] < INTERFERENT_RT - 0.2)
    check("constant baseline removed (areas not inflated by 2000 x width)",
          total < total_exp * 1.03)
    check("neighbouring ion 50 ppm away excluded at 10 ppm",
          abs(area["z7_i4"] / expected["z7_i4"] - 1) < 0.03)
    sub = ["z6_i3", "z6_i4", "z7_i3", "z7_i4"]
    check("subset sum = sum of subset areas (post-hoc selection is exact)",
          abs(sum(area[i] for i in sub) - sum(expected[i] for i in sub)) / sum(expected[i] for i in sub) < 0.03)
    check("S/N is reported and large", r0["sn_all_ions"] > 100)
    check("summed trace written for review", len(traces) > 50)

    # Search window placed on the interferent picks the interferent instead.
    t2 = [dict(t, rt_expected=5.6) for t in targets]
    ions2, _, _ = process_file_xic(path, t2, XICParams(ppm=10))
    check("rt_expected steers which peak is integrated", abs(ions2[0]["rt_apex"] - INTERFERENT_RT) < 0.01)

    out = os.path.join(tmp, "out")
    ions_path, traces_path, failures = run_xic_batch([path, path + ".missing"], targets,
                                                     XICParams(ppm=10), out, n_workers=1)
    check("batch writes ion table", os.path.exists(ions_path))
    check("unreadable file is isolated as a failure, not a crash", len(failures) == 1)

    if _FAILED:
        print(f"\n{len(_FAILED)} check(s) FAILED")
        return 1
    print("\nAll XIC checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
