"""Console PASS/FAIL check for cross-sample RT alignment (align.py).

Builds a synthetic multi-sample feature table: a set of analytes each with a
true consensus RT, replicated across samples with a per-sample drift (a
smooth gradient shift plus jitter). Alignment should recover the consensus
and shrink the cross-sample RT spread, use internal-standard groups as
anchors, and leave an already-aligned table essentially untouched.

Run with:  python -m oligomet_deconv.tests.test_align_synthetic
(run from inst/python/, or with inst/python/ on PYTHONPATH)
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from oligomet_deconv.align import align_retention_times, _lowess_predict

_FAILED = []


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        _FAILED.append(name)


def _make_table(drift_fn, seed=0, jitter=0.002):
    """One feature per (analyte, sample). true_rt is the consensus; observed
    rt = true_rt + drift_fn(sample, true_rt) + jitter."""
    rng = np.random.default_rng(seed)
    analytes = [
        # (mz, true_rt, is_internal_standard)
        (6597.9890, 2.10, True),    # IS 1 (early)
        (5123.4000, 6.80, True),    # IS 2 (late)
        (7051.2000, 3.50, False),
        (6848.1200, 4.25, False),
        (6645.0400, 5.10, False),
        (4001.0000, 7.40, False),
        (8800.0000, 8.90, False),
    ]
    samples = [f"S{i:02d}" for i in range(8)]
    rows = []
    for si, s in enumerate(samples):
        for (mz, trt, _is) in analytes:
            obs = trt + drift_fn(si, trt) + rng.normal(0, jitter)
            rows.append({"sample": s, "mz": mz, "rt": obs,
                         "max_intensity": 1e6 * (1.0 + 0.1 * rng.random()),
                         "true_rt": trt})
    return pd.DataFrame(rows), samples, analytes


def _spread_by_analyte(df, rt_col):
    """Mean cross-sample RT range per analyte (grouped by true_rt)."""
    g = df.groupby("true_rt")[rt_col]
    return float((g.max() - g.min()).mean())


def main():
    print("==== RT alignment synthetic suite ====")

    # ---- loess sanity: recover a known linear deviation --------------------
    x = np.linspace(2, 9, 20)
    y = 0.05 * x - 0.1                       # deviation grows with RT
    pred = _lowess_predict(x, y, x, span=0.6)
    check("loess recovers a linear deviation (max resid < 5e-3)",
          np.max(np.abs(pred - y)) < 5e-3)
    # flat extrapolation beyond the anchor range
    pe = _lowess_predict(x, y, np.array([0.0, 12.0]), span=0.6)
    check("loess holds the deviation flat outside the anchor range",
          abs(pe[0] - y[0]) < 1e-6 and abs(pe[1] - y[-1]) < 1e-6)

    # ---- a smooth per-sample gradient drift --------------------------------
    # each sample shifts RT by a sample-specific slope times RT (column drift)
    def drift(si, trt):
        return (si - 3.5) * 0.01 * trt

    df, samples, analytes = _make_table(drift, seed=1)
    is_mz = [mz for (mz, _t, isx) in analytes if isx]

    before = _spread_by_analyte(df, "rt")
    out, summary, anchors = align_retention_times(
        df, is_mz=is_mz, ppm=20.0, rt_tol=0.20, min_fraction=0.5, loess_span=0.6)
    after = _spread_by_analyte(out, "rt_aligned")

    check("rt_aligned and rt_raw columns are added", "rt_aligned" in out.columns and "rt_raw" in out.columns)
    check("no features dropped or added by alignment", len(out) == len(df))
    check("cross-sample RT spread shrinks after alignment", after < 0.4 * before)
    # residual is floored by the per-point jitter (random, not drift); a
    # working aligner lands within a few jitter widths of the true RT, far
    # below the ~0.19 half-drift a broken one would leave
    check("aligned RT lands near the true consensus (mean abs err < 8e-3)",
          float(np.mean(np.abs(out["rt_aligned"] - out["true_rt"]))) < 8e-3)

    # internal standards were seeded as anchors
    if not anchors.empty:
        is_anchor_mz = anchors.loc[anchors["is_internal_standard"], "mz"].tolist()
        got_both_is = all(any(abs(a - m) / m < 20e-6 for a in is_anchor_mz) for m in is_mz)
        check("both internal standards are used as anchors", got_both_is)
    else:
        check("both internal standards are used as anchors", False)
    check("every sample was aligned (enough anchors)", bool(summary["aligned"].all()))

    # ---- a no-drift table is left essentially untouched --------------------
    df0, _s, _a = _make_table(lambda si, trt: 0.0, seed=2, jitter=0.0)
    out0, _sum0, _anc0 = align_retention_times(df0, is_mz=is_mz, ppm=20.0, rt_tol=0.20)
    check("a drift-free table is left essentially unchanged (max |rt_aligned - rt| < 5e-3)",
          float(np.max(np.abs(out0["rt_aligned"] - out0["rt_raw"]))) < 5e-3)

    # ---- empty input is handled -------------------------------------------
    empty = pd.DataFrame(columns=["sample", "mz", "rt", "max_intensity"])
    oe, se, ae = align_retention_times(empty)
    check("empty feature table returns cleanly with rt_aligned present",
          "rt_aligned" in oe.columns and oe.empty)

    print()
    if _FAILED:
        print(f"==== {len(_FAILED)} check(s) FAILED ====")
        return 1
    print("==== all RT-alignment checks passed ====")
    return 0


if __name__ == "__main__":
    sys.exit(main())
