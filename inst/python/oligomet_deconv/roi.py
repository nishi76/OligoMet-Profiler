"""Streaming ROI (chromatographic-trace) detection.

A simplified, online centWave-style algorithm: peaks are chained across
scans into regions-of-interest (ROIs) by m/z proximity, the same "chain
adjacent values within a ppm tolerance" idea the R package already uses for
single-scan feature grouping (`extract_ms1_features()` in R/ms_matching.R),
extended here across the retention-time axis. ROIs are processed one scan
at a time and only currently-open ROIs are held in memory, so peak memory
is bounded by the number of co-eluting traces, not by file size.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

try:
    _trapz = np.trapezoid  # numpy >= 2.0
except AttributeError:  # pragma: no cover - older numpy
    _trapz = np.trapz


@dataclass
class _OpenROI:
    mz_sum: float
    weight_sum: float
    points: list = field(default_factory=list)  # (rt, mz, intensity)
    gap: int = 0

    @property
    def mz_mean(self) -> float:
        return self.mz_sum / self.weight_sum


@dataclass
class ROIPeak:
    mz: float
    rt_apex: float
    rt_start: float
    rt_end: float
    apex_intensity: float
    area: float
    n_scans: int


class ROIBuilder:
    """Feed scans one at a time via `add_scan()`, then call `finalize()`."""

    def __init__(self, roi_ppm: float = 15.0, max_gap_scans: int = 2,
                 min_intensity: float = 1e4, min_scans: int = 3):
        self.roi_ppm = roi_ppm
        self.max_gap_scans = max_gap_scans
        self.min_intensity = min_intensity
        self.min_scans = min_scans
        self._open: list[_OpenROI] = []
        self._closed_points: list[list] = []

    def add_scan(self, rt: float, mz_array: np.ndarray, intensity_array: np.ndarray) -> None:
        if rt is None or len(mz_array) == 0:
            return
        order = np.argsort(mz_array)
        mzs = mz_array[order]
        ints = intensity_array[order]
        self._open.sort(key=lambda r: r.mz_mean)

        matched = np.zeros(len(mzs), dtype=bool)
        n_open = len(self._open)
        if n_open > 0:
            # A prior version scanned every scan point for every open ROI
            # (a full `for j in range(len(mzs))`, restarting from index 0
            # each time) -- O(n_open_rois * n_points_per_scan) per scan.
            # That's fine for a handful of centroided peaks, but profile-
            # mode HRMS scans routinely carry thousands to tens of
            # thousands of points, and a real acquisition has thousands of
            # scans -- the combination made this pipeline appear to hang
            # indefinitely on real data (a 150-scan x 3000-point/scan
            # synthetic stress test alone took 86s, 88% of it in this
            # loop). `mzs` is already sorted, so each ROI's ppm window can
            # be bounded in O(log n) via searchsorted instead of a full
            # rescan -- and calling searchsorted once on arrays covering
            # ALL open ROIs (rather than once per ROI in a Python loop)
            # avoids per-call numpy dispatch overhead, which otherwise
            # dominates once there are thousands of ROIs per scan. Only
            # the (typically tiny) window itself is then linearly scanned
            # for the closest unmatched point.
            roi_means = np.fromiter((r.mz_mean for r in self._open), dtype=np.float64, count=n_open)
            tols = roi_means * (self.roi_ppm / 1e6)
            los = np.searchsorted(mzs, roi_means - tols, side="left")
            his = np.searchsorted(mzs, roi_means + tols, side="right")

            for idx, roi in enumerate(self._open):
                mean, tol, lo, hi = roi_means[idx], tols[idx], los[idx], his[idx]
                best_j, best_d = -1, tol
                for j in range(lo, hi):
                    if matched[j]:
                        continue
                    d = abs(mzs[j] - mean)
                    if d < best_d:
                        best_d, best_j = d, j
                if best_j >= 0:
                    matched[best_j] = True
                    roi.mz_sum += mzs[best_j] * ints[best_j]
                    roi.weight_sum += ints[best_j]
                    roi.points.append((rt, mzs[best_j], ints[best_j]))
                    roi.gap = 0
                else:
                    roi.gap += 1

        still_open = []
        for roi in self._open:
            if roi.gap > self.max_gap_scans:
                self._closed_points.append(roi.points)
            else:
                still_open.append(roi)
        self._open = still_open

        for j in range(len(mzs)):
            if not matched[j]:
                self._open.append(_OpenROI(mz_sum=mzs[j] * ints[j], weight_sum=ints[j],
                                            points=[(rt, mzs[j], ints[j])]))

    def finalize(self) -> list[ROIPeak]:
        for roi in self._open:
            self._closed_points.append(roi.points)
        self._open = []

        from scipy.signal import find_peaks

        peaks: list[ROIPeak] = []
        for points in self._closed_points:
            if len(points) < self.min_scans:
                continue
            points = sorted(points, key=lambda p: p[0])
            rts = np.array([p[0] for p in points])
            mzs = np.array([p[1] for p in points])
            ints = np.array([p[2] for p in points])

            smoothed = ints
            if len(ints) >= 5:
                kernel = np.ones(3) / 3.0
                smoothed = np.convolve(ints, kernel, mode="same")

            prom = smoothed.max() * 0.1 if smoothed.max() > 0 else 0
            apex_idx, _ = find_peaks(smoothed, prominence=prom)
            if len(apex_idx) == 0:
                apex_idx = np.array([int(np.argmax(smoothed))])

            # split multi-apex ROIs at the midpoint between consecutive
            # apexes -- a simplification documented as such (real valley
            # detection would be more precise but adds complexity this
            # v1 doesn't need)
            bounds = [0] + [
                int((apex_idx[i] + apex_idx[i + 1]) // 2) for i in range(len(apex_idx) - 1)
            ] + [len(points)]

            for i in range(len(apex_idx)):
                lo, hi = bounds[i], bounds[i + 1]
                seg_ints, seg_rts, seg_mzs = ints[lo:hi], rts[lo:hi], mzs[lo:hi]
                if len(seg_ints) < self.min_scans:
                    continue
                apex_i = int(np.argmax(seg_ints))
                apex_intensity = float(seg_ints[apex_i])
                if apex_intensity < self.min_intensity:
                    continue
                area = float(_trapz(seg_ints, seg_rts)) if len(seg_rts) > 1 else apex_intensity
                mz_mean = float(np.average(seg_mzs, weights=seg_ints))
                peaks.append(ROIPeak(
                    mz=mz_mean, rt_apex=float(seg_rts[apex_i]),
                    rt_start=float(seg_rts[0]), rt_end=float(seg_rts[-1]),
                    apex_intensity=apex_intensity, area=area, n_scans=len(seg_ints),
                ))
        return peaks
