"""Generate the bundled pre-processed data-matrix example.

Writes, next to this script:
  data_matrix_wide.csv  -- metabolites x samples (met_id, met_name, kind, <samples>)
  data_matrix_long.csv  -- one row per (metabolite, sample), with sample info columns
  sample_info.csv       -- sample, group, timepoint, sample_type, concentration

Design: 5-level calibration curve in duplicate (1-500 ng/mL) for the parent
and its 3' N-1 metabolite, three QCs, one reagent and one matrix blank, and a
two-arm incubation time course (active vs heat-inactivated, 0/4/24/48 h, n=3).
Parent decays first-order in the active arm and stays flat in the control
arm; truncation metabolites grow in step. Deterministic (fixed seed).
"""
import csv
import math
import os
import random

random.seed(7)
HERE = os.path.dirname(os.path.abspath(__file__))

mets = [
    ("M01", "inotersen", "parent"),
    ("M02", "inotersen 3' N-1", "exo_3p"),
    ("M03", "inotersen 3' N-2", "exo_3p"),
    ("M04", "inotersen 3' N-3", "exo_3p"),
    ("M12", "inotersen 5' N-1", "exo_5p"),
    ("M13", "inotersen 5' N-2", "exo_5p"),
]
# Response factors (signal per ng/mL) and intercepts for the two calibrated species.
slope = {"M01": 2000.0, "M02": 1500.0}
intercept = {"M01": 400.0, "M02": 300.0}

samples = []  # (sample, group, timepoint, sample_type, concentration)
for lvl, conc in enumerate([1, 5, 25, 100, 500], start=1):
    for rep in (1, 2):
        samples.append((f"std_L{lvl}_r{rep}", "", "", "standard", conc))
for name, conc in (("QC_low", 3), ("QC_mid", 50), ("QC_high", 400)):
    samples.append((name, "", "", "quality_control", conc))
samples.append(("RB_01", "", "", "reagent_blank", ""))
samples.append(("MB_01", "", "", "matrix_blank", ""))
for arm in ("active", "heat_inactivated"):
    for t in (0, 4, 24, 48):
        for rep in (1, 2, 3):
            samples.append((f"{arm}_t{t}_r{rep}", arm, t, "unknown", ""))


def noisy(x, cv=0.04):
    return max(0.0, x * (1 + random.gauss(0, cv)))


values = {m[0]: {} for m in mets}
for s, arm, t, stype, conc in samples:
    if stype in ("standard", "quality_control"):
        for mid in slope:
            values[mid][s] = noisy(intercept[mid] + slope[mid] * conc)
        for mid, _, _ in mets:
            if mid not in slope:
                values[mid][s] = noisy(250.0, 0.2)
    elif stype in ("reagent_blank", "matrix_blank"):
        for mid, _, _ in mets:
            values[mid][s] = noisy(200.0 if stype == "reagent_blank" else 350.0, 0.2)
    else:
        k = 0.03 if arm == "active" else 0.001
        parent_conc = 120.0 * math.exp(-k * t)
        lost = 120.0 - parent_conc
        share = {"M02": 0.40, "M03": 0.20, "M04": 0.10, "M12": 0.15, "M13": 0.05}
        values["M01"][s] = noisy(intercept["M01"] + slope["M01"] * parent_conc)
        for mid, frac in share.items():
            rf = slope.get(mid, 1200.0)
            values[mid][s] = noisy(300.0 + rf * (2.0 + lost * frac))
        if t == 0 and arm == "active":
            values["M13"].pop(s)  # not detected -> blank cell

names = [s[0] for s in samples]
with open(os.path.join(HERE, "data_matrix_wide.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["met_id", "met_name", "kind"] + names)
    for mid, mname, kind in mets:
        w.writerow([mid, mname, kind] + [
            f"{values[mid][s]:.1f}" if s in values[mid] else "" for s in names])

with open(os.path.join(HERE, "sample_info.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sample", "group", "timepoint", "sample_type", "concentration"])
    for row in samples:
        w.writerow(row)

with open(os.path.join(HERE, "data_matrix_long.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["met_id", "met_name", "kind", "sample", "intensity",
                "group", "timepoint", "sample_type", "concentration"])
    for mid, mname, kind in mets:
        for s, arm, t, stype, conc in samples:
            if s in values[mid]:
                w.writerow([mid, mname, kind, s, f"{values[mid][s]:.1f}", arm, t, stype, conc])
