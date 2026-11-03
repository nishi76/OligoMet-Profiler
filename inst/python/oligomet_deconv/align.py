"""Cross-sample retention-time alignment (XCMS peakGroups-style).

The batch deconvolver extracts features one file at a time, so nothing ties
the same analyte's retention time together across injections. Over a long
sequence, or across matrices, column and gradient drift move a peak's apex
from run to run. Identity-based matching (accurate mass + charge envelope +
isotope fit) does not depend on RT, but every RT-based step downstream does:
the apex search window in xic_quant, the per-sample consistency QC, and any
match-between-runs rescue of a feature seen in some samples but not others.

This module corrects that drift with the method XCMS uses by default
(groupChromPeaks + adjustRtime "peakGroups"), implemented here in numpy so
the deconv backend keeps a single dependency surface:

  1. Group features across samples by m/z (ppm) and RT (tol) into consensus
     groups, one member per sample (the most intense when a sample has
     several).
  2. Pick anchor groups: those present in a large enough fraction of samples
     with a tight RT spread. Groups whose m/z matches a spiked internal
     standard are seeded as anchors regardless of spread, since an IS is in
     every injection by construction.
  3. For each sample, fit a smooth local regression (loess) of the RT
     deviation (consensus - observed) against observed RT over its anchor
     pairs, then add the smoothed deviation back to every feature's RT.
     Beyond the anchor range the nearest anchor's deviation is held flat,
     as peakGroups does.

The correction is a deviation model, not a hard remap, so a run with no
drift is left essentially untouched and the raw RT is always kept alongside
the aligned one.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd


def _intensity_column(features: pd.DataFrame) -> str:
    for col in ("max_intensity", "area", "apex_intensity", "intensity"):
        if col in features.columns:
            return col
    raise KeyError("features needs one of: max_intensity, area, apex_intensity, intensity")


def _lowess_predict(x: np.ndarray, y: np.ndarray, xq: np.ndarray,
                    span: float = 0.5) -> np.ndarray:
    """Local linear regression (loess) with tricube weights.

    Fits y ~ x locally at each query point in xq. `span` is the fraction of
    anchor points in each local window. Query points outside [min(x),
    max(x)] take the fitted value at the nearest end (flat extrapolation of
    the deviation, matching XCMS peakGroups). With fewer than two anchors it
    returns a constant (the mean deviation, or zero)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xq = np.asarray(xq, dtype=float)
    n = x.size
    if n == 0:
        return np.zeros_like(xq)
    if n == 1:
        return np.full_like(xq, y[0])

    order = np.argsort(x)
    x, y = x[order], y[order]
    # collapse exact-duplicate x (same observed RT) to their mean deviation,
    # so the local fit is not dominated by a repeated abscissa
    ux, inv = np.unique(x, return_inverse=True)
    if ux.size < x.size:
        uy = np.zeros_like(ux)
        np.add.at(uy, inv, y)
        counts = np.bincount(inv, minlength=ux.size)
        uy = uy / counts
        x, y = ux, uy
        n = x.size
        if n == 1:
            return np.full_like(xq, y[0])

    k = max(2, int(np.ceil(span * n)))
    k = min(k, n)
    xr = np.clip(xq, x[0], x[-1])
    out = np.empty_like(xr)
    for i, q in enumerate(xr):
        d = np.abs(x - q)
        idx = np.argsort(d)[:k]
        xd, yd, dd = x[idx], y[idx], d[idx]
        dmax = dd.max()
        if dmax <= 0:
            out[i] = yd.mean()
            continue
        w = (1.0 - (dd / dmax) ** 3) ** 3
        w = np.clip(w, 0.0, None)
        sw = w.sum()
        if sw <= 0:
            out[i] = yd.mean()
            continue
        # weighted linear fit y = a + b*x
        mx = np.sum(w * xd) / sw
        my = np.sum(w * yd) / sw
        vxx = np.sum(w * (xd - mx) ** 2)
        if vxx <= 0:
            out[i] = my
            continue
        b = np.sum(w * (xd - mx) * (yd - my)) / vxx
        a = my - b * mx
        out[i] = a + b * q
    return out


def _consensus_groups(features: pd.DataFrame, ppm: float, rt_tol: float,
                      int_col: str) -> list:
    """Cluster features across samples by m/z (ppm) then RT (tol).

    Returns a list of groups, each a dict with the per-sample member rows
    (one per sample, the most intense), the consensus RT, the m/z, and the
    RT spread."""
    if features.empty:
        return []
    f = features[["sample", "mz", "rt", int_col]].copy()
    f = f.dropna(subset=["mz", "rt"])
    f = f.sort_values("mz").reset_index(drop=True)
    mz = f["mz"].to_numpy()
    # split into m/z clusters where a gap exceeds the ppm tolerance
    rel_gap = np.diff(mz) / mz[1:]
    cut = np.where(rel_gap > ppm / 1e6)[0] + 1
    mz_blocks = np.split(np.arange(len(f)), cut) if len(f) > 1 else [np.arange(len(f))]

    groups = []
    for block in mz_blocks:
        if len(block) == 0:
            continue
        sub = f.iloc[block].sort_values("rt")
        rts = sub["rt"].to_numpy()
        if len(rts) > 1:
            rt_cut = np.where(np.diff(rts) > rt_tol)[0] + 1
            rt_blocks = np.split(np.arange(len(sub)), rt_cut)
        else:
            rt_blocks = [np.arange(len(sub))]
        for rb in rt_blocks:
            members = sub.iloc[rb]
            # one feature per sample: keep the most intense
            best = members.sort_values(int_col, ascending=False).drop_duplicates("sample")
            if best.empty:
                continue
            member_rts = best["rt"].to_numpy()
            groups.append({
                "mz": float(best["mz"].median()),
                "samples": best["sample"].tolist(),
                "rt_by_sample": dict(zip(best["sample"], best["rt"])),
                "consensus_rt": float(np.median(member_rts)),
                "rt_spread": float(member_rts.max() - member_rts.min()),
                "n_samples": int(best["sample"].nunique()),
            })
    return groups


def _is_match(mz: float, is_mz: Optional[np.ndarray], ppm: float) -> bool:
    if is_mz is None or is_mz.size == 0:
        return False
    tol = mz * ppm / 1e6
    lo = np.searchsorted(is_mz, mz - tol, side="left")
    hi = np.searchsorted(is_mz, mz + tol, side="right")
    return hi > lo


def align_retention_times(features: pd.DataFrame,
                          is_mz: Optional[Sequence[float]] = None,
                          ppm: float = 15.0,
                          rt_tol: float = 0.15,
                          min_fraction: float = 0.5,
                          min_anchors: int = 3,
                          loess_span: float = 0.5,
                          max_group_rt_spread: Optional[float] = None):
    """Align retention times across samples (peakGroups-style).

    Parameters
    ----------
    features : DataFrame with at least `sample`, `mz`, `rt`, and an intensity
        column (`max_intensity`/`area`/`apex_intensity`/`intensity`).
    is_mz : optional internal-standard theoretical m/z values. Groups within
        `ppm` of one are seeded as anchors regardless of RT spread.
    ppm, rt_tol : grouping tolerances (m/z in ppm, RT in minutes).
    min_fraction : a group must appear in at least this fraction of samples to
        be a (non-IS) anchor.
    min_anchors : a sample needs at least this many anchor pairs to be warped;
        otherwise its RT is left unaligned and flagged.
    loess_span : fraction of anchors in each local regression window.
    max_group_rt_spread : optional cap on a non-IS anchor's RT spread
        (minutes). Defaults to 3 * rt_tol.

    Returns
    -------
    (features_out, summary, anchors) : `features_out` is `features` with a new
        `rt_aligned` column (and `rt_raw` preserving the input RT); `summary`
        is one row per sample; `anchors` is one row per anchor group used.
    """
    out = features.copy()
    if out.empty or "sample" not in out.columns or "rt" not in out.columns:
        out["rt_aligned"] = out.get("rt", pd.Series(dtype=float))
        out["rt_raw"] = out.get("rt", pd.Series(dtype=float))
        return out, pd.DataFrame(), pd.DataFrame()

    int_col = _intensity_column(out)
    samples = list(dict.fromkeys(out["sample"].tolist()))
    n_samples = len(samples)
    if max_group_rt_spread is None:
        max_group_rt_spread = 3.0 * rt_tol
    is_arr = None
    if is_mz is not None and len(is_mz) > 0:
        is_arr = np.sort(np.asarray([float(m) for m in is_mz], dtype=float))

    groups = _consensus_groups(out, ppm=ppm, rt_tol=rt_tol, int_col=int_col)
    min_present = max(2, int(np.ceil(min_fraction * n_samples)))

    anchors = []
    for g in groups:
        is_is = _is_match(g["mz"], is_arr, ppm)
        present_ok = g["n_samples"] >= (2 if is_is else min_present)
        spread_ok = is_is or g["rt_spread"] <= max_group_rt_spread
        if present_ok and spread_ok:
            g["is_is"] = is_is
            anchors.append(g)

    out["rt_raw"] = out["rt"].astype(float)
    out["rt_aligned"] = out["rt"].astype(float)

    summary_rows = []
    for s in samples:
        pairs_obs, pairs_delta = [], []
        for g in anchors:
            if s in g["rt_by_sample"]:
                obs = g["rt_by_sample"][s]
                pairs_obs.append(obs)
                pairs_delta.append(g["consensus_rt"] - obs)
        n_anchor = len(pairs_obs)
        mask = out["sample"] == s
        rt_before = out.loc[mask, "rt"].to_numpy()
        if n_anchor >= min_anchors:
            delta = _lowess_predict(np.array(pairs_obs), np.array(pairs_delta),
                                    rt_before, span=loess_span)
            out.loc[mask, "rt_aligned"] = rt_before + delta
            applied = True
            med_shift = float(np.median(delta))
        else:
            applied = False
            med_shift = 0.0
        summary_rows.append({
            "sample": s,
            "n_anchors": n_anchor,
            "aligned": applied,
            "median_shift": med_shift,
            "max_abs_shift": float(np.max(np.abs(out.loc[mask, "rt_aligned"].to_numpy()
                                                 - rt_before))) if applied else 0.0,
        })

    anchor_rows = [{
        "mz": g["mz"], "consensus_rt": g["consensus_rt"], "rt_spread": g["rt_spread"],
        "n_samples": g["n_samples"], "is_internal_standard": bool(g.get("is_is", False)),
    } for g in anchors]

    return out, pd.DataFrame(summary_rows), pd.DataFrame(anchor_rows)
