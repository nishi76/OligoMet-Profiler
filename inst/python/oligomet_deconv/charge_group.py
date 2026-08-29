"""Charge-state grouping / neutral-mass consensus.

Groups co-eluting ROI peaks (`roi.ROIPeak`) into charge-state envelopes by
checking whether their back-calculated neutral masses agree, using the same
formula as the R package's `envelope_consistency()` (R/ms_matching.R):

    M = z * (mz - PROTON)      # negative-ESI [M-zH]^z-, h_offset = 0

This deliberately does NOT use an isotope/averagine model -- for these
highly charged oligonucleotide species, charge-state agreement across
observed m/z (not isotope shape) is the reliable deconvolution signal.
Groups are chained by adjacent-gap ppm tolerance, the same transitive
grouping convention the R side already uses in `extract_ms1_features()`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .roi import ROIPeak

PROTON = 1.007276466  # matches R's .PROTON in R/mass_isotope.R


@dataclass
class ChargeGroup:
    mz: float
    charge: int
    neutral_mass: float
    n_charge_states: int
    mass_cv_ppm: float
    rt: float
    rt_start: float
    rt_end: float
    apex_intensity: float
    area: float
    n_scans: int


def _rt_clusters(peaks: list[ROIPeak], rt_tol: float) -> list[list[ROIPeak]]:
    ordered = sorted(peaks, key=lambda p: p.rt_apex)
    clusters, cur = [], []
    for p in ordered:
        if cur and p.rt_apex - cur[-1].rt_apex > rt_tol:
            clusters.append(cur)
            cur = []
        cur.append(p)
    if cur:
        clusters.append(cur)
    return clusters


def group_charge_states(peaks: list[ROIPeak], z_min: int = 3, z_max: int = 20,
                         rt_tol: float = 0.15, mass_tol_ppm: float = 20.0,
                         min_mass: float = 200.0, max_mass: float = 50000.0) -> list[ChargeGroup]:
    z_range = range(z_min, z_max + 1)
    groups: list[ChargeGroup] = []

    for cluster in _rt_clusters(peaks, rt_tol):
        candidates = []  # (peak_idx, z, M)
        for idx, p in enumerate(cluster):
            for z in z_range:
                M = z * (p.mz - PROTON)
                if min_mass <= M <= max_mass:
                    candidates.append((idx, z, M))
        candidates.sort(key=lambda c: c[2])

        chains, cur = [], []
        for cand in candidates:
            if cur and (cand[2] - cur[-1][2]) > cur[-1][2] * mass_tol_ppm / 1e6:
                chains.append(cur)
                cur = []
            cur.append(cand)
        if cur:
            chains.append(cur)

        for chain in chains:
            ref_M = float(np.mean([c[2] for c in chain]))
            # a peak can appear more than once in a chain under different z
            # guesses (coincidental agreement across a wide z-range) --
            # keep only its assignment closest to the chain's consensus mass
            by_peak: dict[int, tuple] = {}
            for idx, z, M in chain:
                if idx not in by_peak or abs(M - ref_M) < abs(by_peak[idx][2] - ref_M):
                    by_peak[idx] = (idx, z, M)
            entries = list(by_peak.values())
            if not entries:
                continue

            masses = np.array([e[2] for e in entries])
            intensities = np.array([cluster[e[0]].apex_intensity for e in entries])
            mean_mass = float(np.average(masses, weights=intensities))
            cv_ppm = float(np.std(masses) / mean_mass * 1e6) if len(masses) > 1 else 0.0

            rep_idx, rep_z, _ = entries[int(np.argmax(intensities))]
            rep_peak = cluster[rep_idx]

            groups.append(ChargeGroup(
                mz=rep_peak.mz, charge=rep_z, neutral_mass=mean_mass,
                n_charge_states=len(entries), mass_cv_ppm=cv_ppm,
                rt=float(np.average([cluster[e[0]].rt_apex for e in entries], weights=intensities)),
                rt_start=min(cluster[e[0]].rt_start for e in entries),
                rt_end=max(cluster[e[0]].rt_end for e in entries),
                apex_intensity=rep_peak.apex_intensity, area=rep_peak.area,
                n_scans=rep_peak.n_scans,
            ))
    return groups
