"""Targeted summed-XIC quantitation (Chromeleon-style "summation of ions").

For every target (one identified metabolite), a fixed list of ions -- charge
states x isotope clusters, with theoretical m/z computed on the R side (see
build_xic_targets() in R/xic_quant.R) -- is extracted from every MS1 scan as
an extracted-ion chromatogram (XIC): the sum of all centroid (or profile)
intensities within +/- ppm of the ion's m/z. The XICs of a target are then
summed point by point into one chromatogram, which is integrated ONCE:

  1. apex = maximum of the smoothed summed trace inside rt_expected +/- rt_window
  2. start/end = walk out from the apex until the smoothed trace drops below
     background + edge_frac x (apex - background) or turns back up (a
     valley), capped at max_half_width
  3. baseline = straight line between the trace values at start and end
  4. area = trapezoidal integral of (trace - baseline) from start to end

Every ion is integrated over that SAME window with its own linear baseline,
and the per-ion areas are written out. Because the baseline is linear, the
sum of per-ion areas over any ion subset equals the area of that subset's
summed chromatogram over the same window -- so the choice of how many charge
states and isotopes to sum can be made afterwards, on the R side, without
re-reading any file.

One streaming pass per file, O(n_points log n_points + n_ions) per scan via
a cumulative-sum lookup, so hundreds of ions cost almost nothing extra.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass

import numpy as np

from .io_mzml import iter_scans

try:
    _trapz = np.trapezoid  # numpy >= 2.0
except AttributeError:  # pragma: no cover - older numpy
    _trapz = np.trapz


@dataclass
class XICParams:
    ppm: float = 10.0
    rt_window: float = 0.5        # min, apex search window around rt_expected
    max_half_width: float = 0.6   # min, cap on each side of the apex
    edge_frac: float = 0.01       # peak edge: smoothed trace < background + edge_frac x (apex - background)
    smooth_points: int = 7        # Savitzky-Golay window (odd); <5 disables smoothing
    trace_margin: float = 0.5     # min, extra trace written out around the search window


def extract_xics(path: str, ion_mz: np.ndarray, ppm: float):
    """Stream MS1 scans once; return (rts[n_scans], xic[n_scans, n_ions]).

    xic[i, j] = sum of intensities in scan i within ion_mz[j] +/- ppm.
    """
    ion_mz = np.asarray(ion_mz, dtype=np.float64)
    lo_mz = ion_mz * (1.0 - ppm / 1e6)
    hi_mz = ion_mz * (1.0 + ppm / 1e6)
    rts, rows = [], []
    for scan in iter_scans(path):
        if scan["ms_level"] != 1 or scan["rt"] is None:
            continue
        mz, it = scan["mz"], scan["intensity"]
        if len(mz) > 1 and np.any(mz[1:] < mz[:-1]):
            order = np.argsort(mz)
            mz, it = mz[order], it[order]
        csum = np.concatenate(([0.0], np.cumsum(it)))
        lo = np.searchsorted(mz, lo_mz, side="left")
        hi = np.searchsorted(mz, hi_mz, side="right")
        rts.append(float(scan["rt"]))
        rows.append((csum[hi] - csum[lo]).astype(np.float32))
    if not rts:
        return np.zeros(0), np.zeros((0, len(ion_mz)), dtype=np.float32)
    rts = np.asarray(rts)
    order = np.argsort(rts, kind="stable")
    return rts[order], np.vstack(rows)[order]


def _smooth(y: np.ndarray, points: int) -> np.ndarray:
    if points < 5 or len(y) < points:
        return y.astype(np.float64)
    try:
        from scipy.signal import savgol_filter
        return savgol_filter(y.astype(np.float64), points if points % 2 else points + 1, 2)
    except Exception:  # pragma: no cover - scipy missing
        k = np.ones(3) / 3.0
        return np.convolve(y, k, mode="same")


def integrate_target(rts: np.ndarray, xic: np.ndarray, rt_expected: float, params: XICParams):
    """Integrate one target's ions over one shared window.

    `xic` is [n_scans, n_ions] for this target only. Returns a dict with the
    window (rt_apex, rt_start, rt_end, i_start, i_end), the summed trace's
    `height` and `sn`, and per-ion `areas`/`heights` arrays.
    """
    n_ions = xic.shape[1]
    empty = dict(found=False, rt_apex=math.nan, rt_start=math.nan, rt_end=math.nan,
                 i_start=-1, i_end=-1, height=0.0, sn=math.nan,
                 areas=np.zeros(n_ions), heights=np.zeros(n_ions))
    if len(rts) < 3 or n_ions == 0:
        return empty
    y = xic.sum(axis=1).astype(np.float64)
    s = _smooth(y, params.smooth_points)

    if rt_expected is None or not np.isfinite(rt_expected):
        search = np.ones(len(rts), dtype=bool)
    else:
        search = np.abs(rts - rt_expected) <= params.rt_window
    idx = np.flatnonzero(search)
    if len(idx) == 0:
        return empty
    apex = int(idx[np.argmax(s[idx])])
    if s[apex] <= 0:
        return empty

    # Edge threshold measured from the local background (5th percentile of
    # the smoothed trace in a region three search windows wide), so a
    # constant chemical background doesn't hold the edges up near the apex.
    region = np.abs(rts - rts[apex]) <= 3 * max(params.rt_window, params.max_half_width)
    bg = float(np.percentile(s[region], 5)) if region.any() else 0.0
    floor = bg + params.edge_frac * (s[apex] - bg)
    left = apex
    while left > 0 and rts[apex] - rts[left - 1] <= params.max_half_width:
        if s[left - 1] <= floor:
            left -= 1
            break
        if s[left - 1] > s[left]:  # trace turns back up: valley
            break
        left -= 1
    right = apex
    n = len(rts)
    while right < n - 1 and rts[right + 1] - rts[apex] <= params.max_half_width:
        if s[right + 1] <= floor:
            right += 1
            break
        if s[right + 1] > s[right]:
            break
        right += 1
    if right - left < 2:
        return empty

    t = rts[left:right + 1]
    seg = xic[left:right + 1, :].astype(np.float64)
    frac = (t - t[0]) / (t[-1] - t[0]) if t[-1] > t[0] else np.zeros_like(t)
    base = seg[0, :][None, :] + frac[:, None] * (seg[-1, :] - seg[0, :])[None, :]
    corr = seg - base
    areas = _trapz(corr, t, axis=0)
    heights = corr[apex - left, :]

    # Noise from the summed trace outside the integration window, within a
    # region three windows wide -- robust MAD of first differences / sqrt(2),
    # so a sloping baseline doesn't inflate it.
    outside = region.copy()
    outside[left:right + 1] = False
    noise = math.nan
    if outside.sum() >= 5:
        d = np.diff(y[outside])
        mad = np.median(np.abs(d - np.median(d))) * 1.4826 / math.sqrt(2)
        noise = float(mad) if mad > 0 else math.nan
    height = float(heights.sum())
    sn = height / noise if noise and np.isfinite(noise) and noise > 0 else math.nan

    return dict(found=True, rt_apex=float(rts[apex]), rt_start=float(t[0]), rt_end=float(t[-1]),
                i_start=left, i_end=right, height=height, sn=sn,
                areas=np.asarray(areas, dtype=np.float64),
                heights=np.asarray(heights, dtype=np.float64))


def process_file_xic(path: str, targets: list[dict], params: XICParams):
    """targets: dicts with target_id, ion_id, mz, rt_expected (one row per ion).

    Returns (ion_rows, trace_rows, meta) for one file.
    """
    sample = os.path.splitext(os.path.basename(path))[0]
    mz = np.array([t["mz"] for t in targets], dtype=np.float64)
    rts, xic = extract_xics(path, mz, params.ppm)

    by_target: dict[str, list[int]] = {}
    for j, t in enumerate(targets):
        by_target.setdefault(t["target_id"], []).append(j)

    ion_rows, trace_rows = [], []
    for tid, cols in by_target.items():
        rt_exp = targets[cols[0]].get("rt_expected")
        rt_exp = float(rt_exp) if rt_exp is not None and str(rt_exp) not in ("", "nan", "NA") else math.nan
        res = integrate_target(rts, xic[:, cols], rt_exp, params)
        for k, j in enumerate(cols):
            ion_rows.append({
                "sample": sample, "target_id": tid, "ion_id": targets[j]["ion_id"],
                "area": float(res["areas"][k]), "height": float(res["heights"][k]),
                "rt_apex": res["rt_apex"], "rt_start": res["rt_start"], "rt_end": res["rt_end"],
                "sn_all_ions": res["sn"],
            })
        # Summed trace of ALL candidate ions around the search window, for
        # reviewing the integration in the app.
        if len(rts):
            centre = rt_exp if np.isfinite(rt_exp) else (res["rt_apex"] if res["found"] else float(np.median(rts)))
            keep = np.abs(rts - centre) <= params.rt_window + params.trace_margin
            ysum = xic[:, cols].sum(axis=1)
            for i in np.flatnonzero(keep):
                trace_rows.append({"sample": sample, "target_id": tid,
                                   "rt": float(rts[i]), "intensity": float(ysum[i])})
    return ion_rows, trace_rows, {"sample": sample, "source_file": path, "n_ms1_scans": int(len(rts))}
