"""Command-line entry point for targeted summed-XIC quantitation:

    python -m oligomet_deconv.xic_cli --input a.mzML b.mzML --targets targets.tsv \
        --output-dir out/ [--ppm 10] [--rt-window 0.5] [--n-workers 4]

targets.tsv: one row per ion with columns target_id, ion_id, mz, rt_expected
(rt_expected may be empty). Writes xic_ions.tsv (one row per sample x ion:
area, height, shared integration window) and xic_traces.tsv (summed trace of
every target's ions around its search window). One worker process per file,
same per-file failure isolation as the deconvolution batch (batch.py).
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback

import pandas as pd

from .cli import _expand_inputs
from .xic import XICParams, process_file_xic

ION_COLUMNS = ["sample", "target_id", "ion_id", "area", "height",
               "rt_apex", "rt_start", "rt_end", "sn_all_ions"]
TRACE_COLUMNS = ["sample", "target_id", "rt", "intensity"]


def _one(args):
    path, targets, params = args
    try:
        ions, traces, meta = process_file_xic(path, targets, params)
        return path, ions, traces, meta, None
    except Exception as exc:  # noqa: BLE001 -- isolate a single-file failure
        return path, [], [], {}, f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"


def run_xic_batch(paths, targets, params: XICParams, output_dir: str, n_workers=None):
    import concurrent.futures as cf

    os.makedirs(output_dir, exist_ok=True)
    n_workers = n_workers or max(1, (os.cpu_count() or 2) - 1)
    all_ions, all_traces, failures = [], [], []
    tasks = [(p, targets, params) for p in paths]
    if n_workers == 1 or len(paths) == 1:
        results = map(_one, tasks)
        for path, ions, traces, meta, err in results:
            if err:
                failures.append({"file": path, "error": err})
            else:
                all_ions.extend(ions)
                all_traces.extend(traces)
    else:
        try:
            with cf.ProcessPoolExecutor(max_workers=n_workers) as ex:
                for path, ions, traces, meta, err in ex.map(_one, tasks):
                    if err:
                        failures.append({"file": path, "error": err})
                    else:
                        all_ions.extend(ions)
                        all_traces.extend(traces)
        except Exception as exc:
            raise RuntimeError(
                f"XIC worker pool failed with {n_workers} worker(s): {type(exc).__name__}: {exc}\n"
                "Retry with --n-workers 1 to process the files sequentially.") from exc

    ions_path = os.path.join(output_dir, "xic_ions.tsv")
    traces_path = os.path.join(output_dir, "xic_traces.tsv")
    pd.DataFrame(all_ions, columns=ION_COLUMNS).to_csv(ions_path, sep="\t", index=False)
    pd.DataFrame(all_traces, columns=TRACE_COLUMNS).to_csv(traces_path, sep="\t", index=False)
    if failures:
        pd.DataFrame(failures).to_csv(os.path.join(output_dir, "_xic_failed_files.tsv"),
                                      sep="\t", index=False)
    return ions_path, traces_path, failures


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Targeted summed-XIC quantitation")
    p.add_argument("--input", nargs="+", required=True)
    p.add_argument("--targets", required=True, help="TSV: target_id, ion_id, mz, rt_expected")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--ppm", type=float, default=10.0, help="XIC extraction half-width (ppm)")
    p.add_argument("--rt-window", type=float, default=0.5,
                   help="apex search window around rt_expected (min)")
    p.add_argument("--max-half-width", type=float, default=0.6,
                   help="maximum peak half-width on each side of the apex (min)")
    p.add_argument("--edge-frac", type=float, default=0.01,
                   help="peak start/end where the smoothed trace falls within this fraction "
                        "of the apex height above local background")
    p.add_argument("--smooth-points", type=int, default=7,
                   help="Savitzky-Golay window for apex/edge finding (areas use raw data)")
    p.add_argument("--n-workers", type=int, default=None)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    paths = _expand_inputs(args.input)
    if not paths:
        print("No input files matched.", file=sys.stderr)
        return 1
    tdf = pd.read_csv(args.targets, sep="\t", dtype={"target_id": str, "ion_id": str})
    missing = {"target_id", "ion_id", "mz"} - set(tdf.columns)
    if missing:
        print(f"targets file missing column(s): {', '.join(sorted(missing))}", file=sys.stderr)
        return 1
    if "rt_expected" not in tdf.columns:
        tdf["rt_expected"] = float("nan")
    targets = tdf[["target_id", "ion_id", "mz", "rt_expected"]].to_dict("records")
    params = XICParams(ppm=args.ppm, rt_window=args.rt_window, max_half_width=args.max_half_width,
                       edge_frac=args.edge_frac, smooth_points=args.smooth_points)
    ions_path, traces_path, failures = run_xic_batch(paths, targets, params, args.output_dir,
                                                     n_workers=args.n_workers)
    print(f"XIC quantitation: {len(paths)} file(s), {tdf['target_id'].nunique()} target(s), "
          f"{len(tdf)} ion(s); {len(failures)} failed.")
    print(f"Ion table: {ions_path}")
    for f in failures:
        print(f"  FAILED: {f['file']}: {f['error'].splitlines()[0]}", file=sys.stderr)
    return 1 if failures and len(failures) == len(paths) else 0


if __name__ == "__main__":
    sys.exit(main())
