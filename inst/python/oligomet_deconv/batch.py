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
        features, ms2, meta = process_file(path, params, watchlist_path)
        return path, features, ms2, meta, None
    except Exception as exc:  # noqa: BLE001 -- isolate any single-file failure
        return path, [], [], {}, f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"


def run_batch(paths, params: DeconvParams, output_dir: str,
              output_file: str = "combined_features.tsv",
              ms2_output_file: Optional[str] = None,
              precursor_watchlist_path: Optional[str] = None,
              n_workers: Optional[int] = None):
    import concurrent.futures as cf

    os.makedirs(output_dir, exist_ok=True)
    n_workers = n_workers or max(1, (os.cpu_count() or 2) - 1)

    all_features, all_ms2, failures, profile_mode_files, noise_thresholds = [], [], [], [], []
    tasks = [(p, params, precursor_watchlist_path) for p in paths]

    # _process_one() isolates any per-file failure and returns it as data
    # (the `error` string below) -- but a WORKER POOL failure (a spawned
    # worker process dying, or on Windows -- where multiprocessing always
    # uses "spawn", never "fork" -- a pickling/import failure re-creating
    # the worker) escapes that isolation entirely: it surfaces as an
    # exception from executor.map() itself, outside the per-task try/except,
    # and previously crashed run_batch() with whatever raw
    # concurrent.futures traceback happened to result. That's opaque to
    # diagnose from the R side (system2() just reports "non-zero exit"),
    # so turn it into one clear, actionable error instead.
    try:
        with cf.ProcessPoolExecutor(max_workers=n_workers) as executor:
            for path, features, ms2, meta, error in executor.map(_process_one, tasks):
                if error:
                    failures.append({"file": path, "error": error})
                else:
                    all_features.extend(features)
                    all_ms2.extend(ms2)
                    if meta.get("profile_mode_detected"):
                        profile_mode_files.append({"sample": meta["sample"], "source_file": meta["source_file"]})
                    if params.sn_threshold is not None:
                        noise_thresholds.append({
                            "sample": meta["sample"], "source_file": meta["source_file"],
                            "noise_level": meta.get("noise_level"),
                            "sn_threshold": params.sn_threshold,
                            "effective_min_intensity": meta.get("effective_min_intensity"),
                        })
    except Exception as exc:
        raise RuntimeError(
            f"Batch worker pool failed with {n_workers} worker(s): {type(exc).__name__}: {exc}\n"
            "This is a process-pool-level failure (a worker process crashed or "
            "could not be started/pickled), not a single bad input file -- "
            "retrying with --n-workers 1 will confirm whether it's parallelism-"
            "related and still process the files sequentially."
        ) from exc

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

    # ROI/charge-envelope detection is designed for centroided peaks; a
    # profile-mode file still runs (just slowly, treating every raw sample
    # point as a candidate peak), so this is surfaced as its own sidecar
    # rather than silently produced or bundled into `failures` (it's not
    # a failure -- the run still completes and returns a features table).
    if profile_mode_files:
        pd.DataFrame(profile_mode_files).to_csv(
            os.path.join(output_dir, "_profile_mode_warnings.tsv"), sep="\t", index=False)

    # Only produced in S/N-threshold mode (sn_threshold is None otherwise,
    # so this list stays empty) -- the whole point of deriving the
    # threshold from each file's own noise level is that it comes out
    # DIFFERENT per file, so it's worth a visible record of what number
    # actually got applied where, not just the one sn_threshold input.
    if noise_thresholds:
        pd.DataFrame(noise_thresholds).to_csv(
            os.path.join(output_dir, "_noise_thresholds.tsv"), sep="\t", index=False)

    return feat_path, ms2_path, failures, profile_mode_files
