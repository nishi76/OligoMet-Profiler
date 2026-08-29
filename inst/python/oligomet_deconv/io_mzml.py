"""Streaming mzML/mzXML reader.

Wraps pyteomics' iterative parsers behind a single generator, `iter_scans`,
that yields one dict per spectrum in file order. Callers never hold more
than one scan's peak arrays in memory at a time -- the caller (deconvolve.py)
is responsible for keeping only bounded, incremental state (open ROIs,
watch-listed MS2 hits) across the iteration.
"""

from __future__ import annotations

import os
from typing import Iterator, Optional

import numpy as np


def _reader_for(path: str):
    from pyteomics import mzml, mzxml

    ext = os.path.splitext(path)[1].lower()
    if ext == ".mzxml":
        return mzxml.MzXML(path)
    return mzml.MzML(path)


def _scan_start_time(spec: dict) -> Optional[float]:
    scan_list = spec.get("scanList")
    if scan_list and scan_list.get("scan"):
        rt = scan_list["scan"][0].get("scan start time")
        if rt is not None:
            return float(rt)
    # mzXML (and some lenient mzML) expose it directly on the spectrum/scan.
    rt = spec.get("scan start time") or spec.get("retentionTime")
    return float(rt) if rt is not None else None


def _precursor_info(spec: dict) -> tuple[Optional[float], Optional[int]]:
    prec_list = spec.get("precursorList")
    if not prec_list or not prec_list.get("precursor"):
        # mzXML represents the precursor inline on the spectrum dict.
        mz = spec.get("precursorMz")
        if isinstance(mz, list) and mz:
            entry = mz[0]
            return (
                float(entry.get("precursorMz", entry)) if isinstance(entry, dict) else float(entry),
                int(entry["precursorCharge"]) if isinstance(entry, dict) and entry.get("precursorCharge") else None,
            )
        return None, None
    prec = prec_list["precursor"][0]
    sel = prec.get("selectedIonList", {}).get("selectedIon", [{}])
    sel = sel[0] if sel else {}
    mz = sel.get("selected ion m/z")
    z = sel.get("charge state")
    if mz is None:
        iso = prec.get("isolationWindow", {}) or {}
        mz = iso.get("isolation window target m/z")
    return (
        float(mz) if mz is not None else None,
        int(z) if z is not None else None,
    )


def iter_scans(path: str) -> Iterator[dict]:
    """Stream one dict per spectrum: ms_level, rt, mz, intensity, and (for
    MS2) precursor_mz/precursor_z. `mz`/`intensity` are numpy float64 arrays.
    Scans with no peaks or no ms level are skipped.
    """
    with _reader_for(path) as reader:
        for spec in reader:
            ms_level = spec.get("ms level")
            if ms_level is None:
                continue
            mz = spec.get("m/z array")
            intensity = spec.get("intensity array")
            if mz is None or intensity is None or len(mz) == 0:
                continue
            rt = _scan_start_time(spec)
            out = {
                "scan_id": spec.get("id") or spec.get("num") or str(spec.get("index")),
                "ms_level": int(ms_level),
                "rt": rt,
                "mz": np.asarray(mz, dtype=np.float64),
                "intensity": np.asarray(intensity, dtype=np.float64),
                "precursor_mz": None,
                "precursor_z": None,
            }
            if out["ms_level"] == 2:
                out["precursor_mz"], out["precursor_z"] = _precursor_info(spec)
            yield out
