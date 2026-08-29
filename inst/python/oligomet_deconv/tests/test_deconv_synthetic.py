"""Console PASS/FAIL check for the deconvolution pipeline against a small
hand-built synthetic mzML file with known charge-state ladders. Mirrors the
plain console-check style of the R package's own tests/test_*.R scripts
(no pytest dependency required).

Run with:  python -m oligomet_deconv.tests.test_deconv_synthetic
(run from inst/python/, or with inst/python/ on PYTHONPATH)
"""

from __future__ import annotations

import base64
import os
import sys
import tempfile
import zlib

PROTON = 1.007276466

_FAILED = []


def check(name, cond):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}")
    if not cond:
        _FAILED.append(name)


def _b64(vals):
    import struct
    raw = struct.pack(f"<{len(vals)}d", *vals)
    return base64.b64encode(zlib.compress(raw)).decode()


def _binary_array(kind, vals):
    accession = "MS:1000514" if kind == "mz" else "MS:1000515"
    name = "m/z array" if kind == "mz" else "intensity array"
    return (
        f'<binaryDataArray><cvParam accession="{accession}" name="{name}"/>'
        f'<cvParam accession="MS:1000523" name="64-bit float"/>'
        f'<cvParam accession="MS:1000574" name="zlib compression"/>'
        f'<binary>{_b64(vals)}</binary></binaryDataArray>'
    )


def _ms1_spectrum(idx, rt, mz_vals, int_vals):
    n = len(mz_vals)
    return (
        f'<spectrum id="scan={idx}" index="{idx - 1}" defaultArrayLength="{n}">'
        f'<cvParam accession="MS:1000511" name="ms level" value="1"/>'
        f'<scanList count="1"><scan><cvParam accession="MS:1000016" name="scan start time" '
        f'value="{rt}" unitAccession="UO:0000031" unitName="minute"/></scan></scanList>'
        f'<binaryDataArrayList count="2">{_binary_array("mz", mz_vals)}{_binary_array("int", int_vals)}'
        f'</binaryDataArrayList></spectrum>'
    )


def _ms2_spectrum(idx, rt, precursor_mz, precursor_z, mz_vals, int_vals):
    n = len(mz_vals)
    return (
        f'<spectrum id="scan={idx}" index="{idx - 1}" defaultArrayLength="{n}">'
        f'<cvParam accession="MS:1000511" name="ms level" value="2"/>'
        f'<scanList count="1"><scan><cvParam accession="MS:1000016" name="scan start time" '
        f'value="{rt}" unitAccession="UO:0000031" unitName="minute"/></scan></scanList>'
        f'<precursorList count="1"><precursor spectrumRef="scan=1">'
        f'<selectedIonList count="1"><selectedIon>'
        f'<cvParam accession="MS:1000744" name="selected ion m/z" value="{precursor_mz}"/>'
        f'<cvParam accession="MS:1000041" name="charge state" value="{precursor_z}"/>'
        f'</selectedIon></selectedIonList></precursor></precursorList>'
        f'<binaryDataArrayList count="2">{_binary_array("mz", mz_vals)}{_binary_array("int", int_vals)}'
        f'</binaryDataArrayList></spectrum>'
    )


def build_synthetic_mzml(path):
    """Two co-eluting charge envelopes (M=7000, M=9100 Da) plus one
    isolated single-charge-state case (M=5000 Da, low-confidence), and
    two MS2 scans (one on a precursor watch-list, one off it)."""
    # Deliberately non-round masses: round numbers here create spurious
    # exact collisions between unrelated (peak, wrong-z) candidate masses
    # when swept across the full z_min..z_max range (e.g. 7000/7 == 5000/5).
    mz_a = {z: 6543.21 / z + PROTON for z in (5, 6, 7, 8)}
    mz_b = {z: 8765.43 / z + PROTON for z in (6, 7, 9)}
    mz_c = 4321.98 / 10 + PROTON

    w_a = {5: 0.6, 6: 0.8, 7: 1.0, 8: 0.7}
    w_b = {6: 0.5, 7: 1.0, 9: 0.6}
    profile_ab = [4000, 10000, 18000, 25000, 18000, 10000, 4000]
    profile_c = [5000, 15000, 25000, 15000, 5000]
    rts_ab = [0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]
    rts_c = [0.90, 0.95, 1.00, 1.05, 1.10]

    spectra = []
    idx = 1
    for rt, base in zip(rts_ab, profile_ab):
        mzs = sorted(list(mz_a.values()) + list(mz_b.values()))
        ints = []
        for m in mzs:
            za = next((z for z, v in mz_a.items() if v == m), None)
            zb = next((z for z, v in mz_b.items() if v == m), None)
            if za is not None:
                ints.append(base * w_a[za])
            else:
                ints.append(base * w_b[zb])
        spectra.append(_ms1_spectrum(idx, rt, mzs, ints))
        idx += 1

    # MS2 #1: on-watchlist (matches envelope A's z=7 m/z almost exactly)
    spectra.append(_ms2_spectrum(idx, 0.56, mz_a[7] + 0.00002, 7, [100.0, 200.0], [500.0, 300.0]))
    idx += 1
    # MS2 #2: off-watchlist
    spectra.append(_ms2_spectrum(idx, 0.61, 650.0, 3, [90.0], [400.0]))
    idx += 1

    for rt, base in zip(rts_c, profile_c):
        spectra.append(_ms1_spectrum(idx, rt, [mz_c], [base]))
        idx += 1

    xml = (
        '<?xml version="1.0"?>\n<mzML xmlns="http://psi.hupo.org/ms/mzml">\n<run>\n'
        f'<spectrumList count="{idx - 1}">\n' + "\n".join(spectra) + "\n</spectrumList>\n</run>\n</mzML>\n"
    )
    with open(path, "w") as f:
        f.write(xml)
    return mz_a, mz_b, mz_c


def main():
    print("=== oligomet_deconv synthetic pipeline test ===")
    from oligomet_deconv.deconvolve import DeconvParams, process_file

    tmpdir = tempfile.mkdtemp()
    mzml_path = os.path.join(tmpdir, "synthetic.mzML")
    mz_a, mz_b, mz_c = build_synthetic_mzml(mzml_path)

    watchlist_path = os.path.join(tmpdir, "watchlist.txt")
    with open(watchlist_path, "w") as f:
        f.write(f"{mz_a[7]}\n")

    params = DeconvParams(roi_ppm=15.0, rt_tol=0.15, mass_tol_ppm=20.0,
                           z_min=3, z_max=20, min_intensity=1e4, min_scans=3,
                           max_gap_scans=2, ms2_watch_ppm=50.0)
    features, ms2 = process_file(mzml_path, params, precursor_watchlist_path=watchlist_path)

    print(f"\n--- MS1 charge-envelope groups ({len(features)}) ---")
    for feat in features:
        print(f"  neutral_mass={feat['neutral_mass']:.3f}  charge={feat['charge']}  "
              f"n_charge_states={feat['n_charge_states']}  mz={feat['mz']:.4f}")

    group_a = [f for f in features if abs(f["neutral_mass"] - 6543.21) < 0.01]
    group_b = [f for f in features if abs(f["neutral_mass"] - 8765.43) < 0.01]
    group_c = [f for f in features if abs(f["neutral_mass"] - 4321.98) < 0.01]

    check("envelope A (M=6543.21) detected exactly once", len(group_a) == 1)
    if group_a:
        check("envelope A has 4 charge states", group_a[0]["n_charge_states"] == 4)
        check("envelope A representative charge is one of 5/6/7/8", group_a[0]["charge"] in (5, 6, 7, 8))

    check("envelope B (M=8765.43) detected exactly once", len(group_b) == 1)
    if group_b:
        check("envelope B has 3 charge states (not merged with A)", group_b[0]["n_charge_states"] == 3)
        check("envelope B representative charge is one of 6/7/9", group_b[0]["charge"] in (6, 7, 9))

    check("envelope C (M=4321.98) detected exactly once", len(group_c) == 1)
    if group_c:
        check("envelope C is single-charge (low confidence)", group_c[0]["n_charge_states"] == 1)
        check("envelope C charge recovered as 10", group_c[0]["charge"] == 10)

    print(f"\n--- Targeted MS2 extraction ({len(ms2)} captured) ---")
    for row in ms2:
        print(f"  scan={row['ms2_scan_id']}  precursor_mz={row['precursor_mz']}")
    check("exactly one MS2 scan captured (on-watchlist only)", len(ms2) == 1)
    if ms2:
        check("captured MS2 scan is the on-watchlist one",
              abs(ms2[0]["precursor_mz"] - (mz_a[7] + 0.00002)) < 1e-3)

    print()
    if _FAILED:
        print(f"==== {len(_FAILED)} check(s) FAILED ====")
        return 1
    print("==== All oligomet_deconv synthetic tests passed ====")
    return 0


if __name__ == "__main__":
    sys.exit(main())
