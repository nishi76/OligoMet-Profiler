"""Command-line entry point: `python -m oligomet_deconv.cli ...`"""

from __future__ import annotations

import argparse
import glob
import sys

from .batch import run_batch
from .deconvolve import DeconvParams


def _expand_inputs(patterns):
    paths = []
    for pat in patterns:
        matches = sorted(glob.glob(pat))
        paths.extend(matches if matches else [pat])
    seen, out = set(), []
    for p in paths:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Parallel charge-envelope MS1 deconvolution")
    p.add_argument("--input", nargs="+", required=True, help="mzML/mzXML file(s) or glob pattern(s)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--output-file", default="combined_features.tsv")
    p.add_argument("--ms2-output-file", default=None,
                    help="defaults to combined_ms2.tsv when --precursor-watchlist is given")
    p.add_argument("--precursor-watchlist", default=None,
                    help="file of theoretical precursor m/z, one per line; enables targeted MS2 capture")
    p.add_argument("--ms2-watch-ppm", type=float, default=50.0)
    p.add_argument("--roi-ppm", type=float, default=15.0)
    p.add_argument("--rt-tol", type=float, default=0.15)
    p.add_argument("--mass-tol-ppm", type=float, default=20.0)
    p.add_argument("--z-min", type=int, default=3)
    p.add_argument("--z-max", type=int, default=20)
    p.add_argument("--min-intensity", type=float, default=1e4)
    p.add_argument("--min-scans", type=int, default=3)
    p.add_argument("--max-gap-scans", type=int, default=2)
    p.add_argument("--min-mass", type=float, default=200.0)
    p.add_argument("--max-mass", type=float, default=50000.0)
    p.add_argument("--min-charge-states", type=int, default=2,
                    help="drop charge-envelope groups confirmed by fewer than this many "
                         "charge states (default 2 -- a single-charge-state 'group' can't "
                         "actually confirm a charge, and sweeping every peak against every "
                         "z in [--z-min, --z-max] produces many such coincidental single-"
                         "observation matches; set to 1 to see everything the sweep finds)")
    p.add_argument("--n-workers", type=int, default=None)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    paths = _expand_inputs(args.input)
    if not paths:
        print("No input files matched.", file=sys.stderr)
        return 1

    params = DeconvParams(
        roi_ppm=args.roi_ppm, rt_tol=args.rt_tol, mass_tol_ppm=args.mass_tol_ppm,
        z_min=args.z_min, z_max=args.z_max, min_intensity=args.min_intensity,
        min_scans=args.min_scans, max_gap_scans=args.max_gap_scans,
        min_mass=args.min_mass, max_mass=args.max_mass, ms2_watch_ppm=args.ms2_watch_ppm,
        min_charge_states=args.min_charge_states,
    )

    ms2_output_file = args.ms2_output_file
    if args.precursor_watchlist and not ms2_output_file:
        ms2_output_file = "combined_ms2.tsv"

    feat_path, ms2_path, failures, profile_mode_files = run_batch(
        paths, params, args.output_dir, output_file=args.output_file,
        ms2_output_file=ms2_output_file, precursor_watchlist_path=args.precursor_watchlist,
        n_workers=args.n_workers,
    )

    print(f"Processed {len(paths)} file(s); {len(failures)} failed.")
    print(f"Feature table: {feat_path}")
    if ms2_path:
        print(f"MS2 table: {ms2_path}")
    for f in failures:
        print(f"  FAILED: {f['file']}", file=sys.stderr)
    if profile_mode_files:
        names = ", ".join(f["sample"] for f in profile_mode_files)
        print(f"WARNING: {len(profile_mode_files)} file(s) appear to be PROFILE mode "
              f"(not centroided): {names}. ROI/charge-envelope detection is designed for "
              f"centroided peaks -- convert with `msconvert --centroid` first for reliable, "
              f"fast results.", file=sys.stderr)

    return 1 if failures and len(failures) == len(paths) else 0


if __name__ == "__main__":
    sys.exit(main())
