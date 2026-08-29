"""Per-file pipeline: stream a file once, running MS1 ROI/charge-envelope
deconvolution (Stages A+B) and, when a precursor watch-list is supplied,
targeted MS2 spectrum capture (Stage C) in the same pass -- so MS2 capture
costs almost nothing extra on top of the MS1 pass.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .charge_group import group_charge_states
from .io_mzml import iter_scans
from .roi import ROIBuilder


@dataclass
class DeconvParams:
    roi_ppm: float = 15.0
    rt_tol: float = 0.15
    mass_tol_ppm: float = 20.0
    z_min: int = 3
    z_max: int = 20
    min_intensity: float = 1e4
    min_scans: int = 3
    max_gap_scans: int = 2
    min_mass: float = 200.0
    max_mass: float = 50000.0
    ms2_watch_ppm: float = 50.0


def _load_watchlist(path: Optional[str]) -> Optional[np.ndarray]:
    if not path:
        return None
    with open(path) as f:
        vals = [float(line.strip()) for line in f if line.strip()]
    return np.array(sorted(vals))


def _on_watchlist(precursor_mz, watchlist: Optional[np.ndarray], ppm: float) -> bool:
    if precursor_mz is None or watchlist is None or len(watchlist) == 0:
        return False
    tol = precursor_mz * ppm / 1e6
    lo = np.searchsorted(watchlist, precursor_mz - tol, side="left")
    hi = np.searchsorted(watchlist, precursor_mz + tol, side="right")
    return hi > lo


def process_file(path: str, params: DeconvParams, precursor_watchlist_path: Optional[str] = None):
    """Returns (features: list[dict], ms2: list[dict]) for one file."""
    watchlist = _load_watchlist(precursor_watchlist_path)
    roi_builder = ROIBuilder(roi_ppm=params.roi_ppm, max_gap_scans=params.max_gap_scans,
                              min_intensity=params.min_intensity, min_scans=params.min_scans)
    ms2_rows = []
    sample = os.path.splitext(os.path.basename(path))[0]

    for scan in iter_scans(path):
        if scan["ms_level"] == 1:
            roi_builder.add_scan(scan["rt"], scan["mz"], scan["intensity"])
        elif scan["ms_level"] == 2 and watchlist is not None:
            if _on_watchlist(scan["precursor_mz"], watchlist, params.ms2_watch_ppm):
                ms2_rows.append({
                    "sample": sample, "source_file": path, "ms2_scan_id": scan["scan_id"],
                    "precursor_mz": scan["precursor_mz"], "precursor_z": scan["precursor_z"],
                    "rt": scan["rt"],
                    "mz_list": ";".join(f"{v:.6f}" for v in scan["mz"]),
                    "intensity_list": ";".join(f"{v:.2f}" for v in scan["intensity"]),
                })

    roi_peaks = roi_builder.finalize()
    groups = group_charge_states(
        roi_peaks, z_min=params.z_min, z_max=params.z_max, rt_tol=params.rt_tol,
        mass_tol_ppm=params.mass_tol_ppm, min_mass=params.min_mass, max_mass=params.max_mass,
    )

    features = []
    for i, g in enumerate(groups):
        features.append({
            "sample": sample, "source_file": path, "feature_id": f"{sample}_F{i + 1}",
            "mz": g.mz, "rt": g.rt, "max_intensity": g.apex_intensity, "n_scans": g.n_scans,
            "charge": g.charge, "neutral_mass": g.neutral_mass,
            "n_charge_states": g.n_charge_states, "mass_cv_ppm": g.mass_cv_ppm,
            "rt_start": g.rt_start, "rt_end": g.rt_end, "area": g.area,
        })
    return features, ms2_rows
