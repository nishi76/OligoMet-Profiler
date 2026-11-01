"""Command-line entry point for report evidence extraction:

    python -m oligomet_deconv.evidence_cli --requests requests.tsv --output-dir out/

requests.tsv: one row per request with columns req_id, file, rt, mz_lo, mz_hi,
rt_lo, rt_hi, xic_mz (semicolon-separated m/z list), ppm. For each request the
file is streamed once (requests on the same file share one pass) and two
tables are written:

  evidence_spectra.tsv  req_id, scan_id, scan_rt, mz, intensity
      the MS1 peaks inside [mz_lo, mz_hi] from the MS1 scan closest to rt
      (the apex spectrum shown as the isotope pattern in the report);
  evidence_xic.tsv      req_id, rt, intensity
      the summed extracted-ion chromatogram of the xic_mz ions (+/- ppm)
      over every MS1 scan in [rt_lo, rt_hi].

Only the requested windows are kept in memory, so multi-GB files stream with
a small footprint. A file that cannot be read is reported on stderr and its
requests are skipped; the other files still produce output.
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback

import numpy as np
import pandas as pd

from .io_mzml import iter_scans


def _process_file(path: str, reqs: pd.DataFrame):
    spectra, xic = [], []
    best = {}  # req_id -> (|drt|, scan_id, scan_rt, mz, intensity)
    rows = list(reqs.itertuples(index=False))
    targets = {r.req_id: np.array([float(x) for x in str(r.xic_mz).split(";") if x], dtype=np.float64)
               for r in rows}
    for scan in iter_scans(path):
        if scan["ms_level"] != 1 or scan["rt"] is None:
            continue
        rt = scan["rt"]
        mz = scan["mz"]
        inten = scan["intensity"]
        order = None
        for r in rows:
            if r.rt_lo <= rt <= r.rt_hi:
                if order is None:
                    order = np.argsort(mz)
                    mz_s, in_s = mz[order], inten[order]
                tg = targets[r.req_id]
                total = 0.0
                for t in tg:
                    tol = t * r.ppm * 1e-6
                    lo, hi = np.searchsorted(mz_s, [t - tol, t + tol])
                    total += float(in_s[lo:hi].sum())
                xic.append((r.req_id, rt, total))
            d = abs(rt - r.rt)
            if r.req_id not in best or d < best[r.req_id][0]:
                keep = (mz >= r.mz_lo) & (mz <= r.mz_hi)
                best[r.req_id] = (d, scan["scan_id"], rt, mz[keep].copy(), inten[keep].copy())
    for req_id, (_, scan_id, scan_rt, m, i) in best.items():
        for a, b in zip(m, i):
            spectra.append((req_id, scan_id, scan_rt, float(a), float(b)))
    return spectra, xic


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--requests", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args(argv)

    os.makedirs(args.output_dir, exist_ok=True)
    reqs = pd.read_csv(args.requests, sep="\t", dtype={"req_id": str, "file": str, "xic_mz": str})
    all_spec, all_xic, failed = [], [], 0
    for path, grp in reqs.groupby("file", sort=False):
        try:
            s, x = _process_file(path, grp)
            all_spec.extend(s)
            all_xic.extend(x)
        except Exception as exc:  # noqa: BLE001 -- isolate a single-file failure
            failed += 1
            sys.stderr.write(f"FAILED {path}: {type(exc).__name__}: {exc}\n{traceback.format_exc()}\n")
    pd.DataFrame(all_spec, columns=["req_id", "scan_id", "scan_rt", "mz", "intensity"]).to_csv(
        os.path.join(args.output_dir, "evidence_spectra.tsv"), sep="\t", index=False)
    pd.DataFrame(all_xic, columns=["req_id", "rt", "intensity"]).to_csv(
        os.path.join(args.output_dir, "evidence_xic.tsv"), sep="\t", index=False)
    print(f"evidence: {len(reqs)} request(s), {len(all_spec)} spectrum peak(s), "
          f"{len(all_xic)} XIC point(s), {failed} file(s) failed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
