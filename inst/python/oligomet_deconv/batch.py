"""Parallel batch driver: one worker process per input file, with per-file
error isolation so a corrupt/unreadable file logs to `_failed_files.tsv`
rather than aborting the whole batch.
"""

from __future__ import annotations

import os
import traceback
from typing import Optional

import pandas as pd

from .deconvolve import DeconvParams, process_file

FEATURE_COLUMNS = ["sample", "source_file", "feature_id", "mz", "rt", "max_intensity",
                    "n_scans", "charge", "neutral_mass", "n_charge_states", "mass_cv_ppm",
                    "rt_start", "rt_end", "area"]
MS2_COLUMNS = ["sample", "source_file", "ms2_scan_id", "precursor_mz", "precursor_z",
               "rt", "mz_list", "intensity_list"]


def _process_one(args):
    path, params, watchlist_path = args
    try:
        features, ms2 = process_file(path, params, watchlist_path)
        return path, features, ms2, None
    except Exception as exc:  # noqa: BLE001 -- isolate any single-file failure
        return path, [], [], f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"


def run_batch(paths, params: DeconvParams, output_dir: str,
              output_file: str = "combined_features.tsv",
              ms2_output_file: Optional[str] = None,
              precursor_watchlist_path: Optional[str] = None,
              n_workers: Optional[int] = None):
    import concurrent.futures as cf

    os.makedirs(output_dir, exist_ok=True)
    n_workers = n_workers or max(1, (os.cpu_count() or 2) - 1)

    all_features, all_ms2, failures = [], [], []
    tasks = [(p, params, precursor_watchlist_path) for p in paths]

    with cf.ProcessPoolExecutor(max_workers=n_workers) as executor:
        for path, features, ms2, error in executor.map(_process_one, tasks):
            if error:
                failures.append({"file": path, "error": error})
            else:
                all_features.extend(features)
                all_ms2.extend(ms2)

    feat_df = pd.DataFrame(all_features, columns=FEATURE_COLUMNS)
    feat_path = os.path.join(output_dir, output_file)
    feat_df.to_csv(feat_path, sep="\t", index=False)

    ms2_path = None
    if ms2_output_file:
        ms2_df = pd.DataFrame(all_ms2, columns=MS2_COLUMNS)
        ms2_path = os.path.join(output_dir, ms2_output_file)
        ms2_df.to_csv(ms2_path, sep="\t", index=False)

    if failures:
        pd.DataFrame(failures).to_csv(os.path.join(output_dir, "_failed_files.tsv"), sep="\t", index=False)

    return feat_path, ms2_path, failures
