#!/usr/bin/env python3
"""Generate cases.csv for the transverse-exchange study (run once; edit here, not in the CSV).

Families (W = 1 um, AR 1 and AR 5, t_end = 0.4 ms, 40 cells/W, residuals 1e-6, Courant 1):
  A   baseline: original rinse U(t), symmetry top, v = 0          (Lin, Busnaina & Suni 2002)
  B   + steady transverse inflow V = lambda * u_avg through the top (lambda = 0 checks the top-BC change)
  Bf  + V(t) = lambda * U(t): fixed flow-approach angle phi = atan(lambda)
  C   + pulsed V(t), mean lambda * u_avg, own frequency and phase relative to the rinse pulses
  V   verification: time step halved (cfl 0.5) and mesh doubled (80 cells/W) for one B case
  R   angle sweep at CONSTANT speed: whole stream enters through the top at phi = 0..90 deg,
      (u, v) = S(t) (cos phi, -sin phi); both sides open, trench centred and split upper/lower.
      Geometries: 1 um AR 1 and AR 5 (0.4 ms, steady and 200 kHz), 100 um AR 5 (0.1 s, convection-dominated).
  A0  exact copy of the original baseline mesh (no trench split) for the reproduction check.
  S   trench calibration: steady near-wall shear u = G y, G = 10 ... 2e5 1/s, AR 1 and AR 5, taller graded
      domain (H = 10 W). Gives cleaning vs local shear rate, the link to the wafer-scale flow.
  T   wafer-scale tank (one inlet, tilted wafer, steady k-omega SST, trenchflow/tank.py): jet meets the
      300 mm wafer at 15..90 deg, aimed at the centre from 30 mm; output = shear rate along the wafer.
All other cases split the trench into upper/lower halves (report zones only) for the deep-mixing metric.
Comparison grid at the angles 15/30/45/60/75 deg: steady V = tan(phi) u_avg (B, steady rinse) and
V = tan(phi) U(t) (Bf, 200 kHz) against the constant-speed rotation R at the same angles.
Walltimes: 2.5 x estimate from the measured cost of conv_f20k_n40tight (2.8 h, 7540 steps, 25 600 cells).
"""
import csv
import math
from pathlib import Path

from trenchflow.physics import plan

HERE = Path(__file__).resolve().parent
COLS = ["case", "family", "W_m", "AR", "mode", "u_avg_m_s", "freq_Hz", "t_end_s", "n_per_w", "walltime",
        "v_mode", "v_ratio", "v_freq_Hz", "v_phase_deg", "angle_deg", "split", "cfl", "res_tol", "iters", "note",
        # Case S (shear-driven trench) and Case T (wafer-scale tank); appended so field 10 stays walltime
        "H_W", "Lu_W", "Ld_W", "shear_rate", "standoff_mm", "slot_mm", "aim_mm", "wafer_mm", "refine"]
SHEAR_RATES = (("10", 10), ("100", 100), ("1k", 1e3), ("10k", 1e4), ("100k", 1e5), ("200k", 2e5))
TANK_ANGLES = (15, 30, 45, 60, 75, 90)
RINSE = {"steady": ("steady", 0), "osc200k": ("osc", 200000)}


def row(case, family, ar, rinse, note, v_mode="none", lam="", f_v="", phase="", n=40, cfl="",
        W="1e-6", t_end="4e-4", angle="", split=""):
    mode, f = RINSE[rinse]
    return dict(case=case, family=family, W_m=W, AR=ar, mode=mode, u_avg_m_s="0.15", freq_Hz=f,
                t_end_s=t_end, n_per_w=n, walltime="", v_mode=v_mode, v_ratio=lam, v_freq_Hz=f_v,
                v_phase_deg=phase, angle_deg=angle, split=split, cfl=cfl, res_tol="1e-6", iters="20", note=note)


def rows():
    out = [row("A0_ar1_osc200k", "A0", 1, "osc200k", "original baseline mesh (no split): reproduction check",
               split="0")]
    for ar in (1, 5):
        g = f"ar{ar}"
        for rinse in RINSE:
            out.append(row(f"A_{g}_{rinse}", "A", ar, rinse, f"baseline AR {ar}, {rinse} rinse"))
            for lam in (0, 0.1, 0.25, 0.5, 1, 2):
                out.append(row(f"B_{g}_{rinse}_v{lam:g}", "B", ar, rinse,
                               f"steady V = {lam:g} u_avg (phi = {math.degrees(math.atan(lam)):.0f} deg at u_avg)",
                               "steady", lam))
        for lam in (0.5, 1):
            out.append(row(f"Bf_{g}_osc200k_v{lam:g}", "Bf", ar, "osc200k",
                           f"V = {lam:g} U(t): fixed phi = {math.degrees(math.atan(lam)):.0f} deg", "follow", lam))
        for ang in (15, 30, 60, 75):             # 45 deg = lambda 1, already above
            lam = f"{math.tan(math.radians(ang)):.6g}"
            out.append(row(f"B_{g}_steady_a{ang:02d}", "B", ar, "steady",
                           f"steady V = tan({ang}) u_avg: flow at {ang} deg, U kept", "steady", lam))
            out.append(row(f"Bf_{g}_osc200k_a{ang:02d}", "Bf", ar, "osc200k",
                           f"V = tan({ang}) U(t): flow at {ang} deg, U kept", "follow", lam))
        for f_v, ph, txt in ((200000, 0, "in phase"), (200000, 180, "anti-phase"), (20000, 0, "20 kHz")):
            out.append(row(f"C_{g}_osc200k_vp{f_v // 1000}k_ph{ph}", "C", ar, "osc200k",
                           f"pulsed V, mean 0.5 u_avg, {txt}", "pulse", 0.5, f_v, ph))
    out.append(row("V_ar5_osc200k_v0.5_dt", "V", 5, "osc200k", "time-step check: B_ar5_osc200k_v0.5 with cfl 0.5",
                   "steady", 0.5, cfl="0.5"))
    out.append(row("V_ar5_osc200k_v0.5_n80", "V", 5, "osc200k", "mesh check: B_ar5_osc200k_v0.5 at 80 cells/W",
                   "steady", 0.5, n=80))
    for W, wtag, ar, rinse, t_end in (("1e-6", "w1um", 1, "steady", "4e-4"), ("1e-6", "w1um", 5, "steady", "4e-4"),
                                      ("1e-6", "w1um", 1, "osc200k", "4e-4"), ("1e-6", "w1um", 5, "osc200k", "4e-4"),
                                      ("1e-4", "w100um", 5, "steady", "0.1")):
        for ang in (0, 15, 30, 45, 60, 75, 90):
            out.append(row(f"R_{wtag}_ar{ar}_{rinse}_a{ang:02d}", "R", ar, rinse,
                           f"same speed, stream at {ang} deg from the wafer", "rotate",
                           W=W, t_end=t_end, angle=ang))
    # S: trench calibration in a near-wall shear flow u = G y (what a trench sees on a real wafer)
    for ar in (1, 5):
        for tag, G in SHEAR_RATES:
            r = row(f"S_ar{ar}_G{tag}", "S", ar, "steady", f"near-wall shear G = {G:g} 1/s", "shear")
            r.update(H_W="10", Lu_W="3", Ld_W="10", shear_rate=f"{G:g}", u_avg_m_s="0")
            out.append(r)
    # T: wafer-scale tank, one inlet, wafer tilted so the jet meets it at angle phi (steady RANS)
    for ang in TANK_ANGLES:
        out.append(tank_row(f"T_a{ang:02d}", ang, f"one inlet, jet at {ang} deg to the wafer, aimed at the centre"))
    out.append(tank_row("T_a45_r2", 45, "mesh check: T_a45 with cells refined 2x", refine="2"))
    return out


def tank_row(case, ang, note, refine="1"):
    r = {c: "" for c in COLS}
    r.update(case=case, family="T", mode="steady", u_avg_m_s="0.15", freq_Hz="0", angle_deg=ang,
             res_tol="1e-6", iters="6000", note=note, standoff_mm="30", slot_mm="10", aim_mm="0",
             wafer_mm="300", refine=refine)
    return r


def walltime(r):
    if r["family"] == "T":
        return ("12:00:00" if r["refine"] not in ("", "1") else "08:00:00"), 0.0
    p = plan(r)
    n = p["n_per_w"]
    if p.get("grade"):
        from trenchflow import mesh
        cells = mesh.build(p["W"], p["W"] * p["AR"], p["H"], p["Lu_W"] * p["W"], p["Ld_W"] * p["W"], n,
                           True, p["grade"])[5]["n_cells"]
    else:
        cells = (p["Lu_W"] + 1 + p["Ld_W"]) * 2.5 * n * n + p["AR"] * n * n
    est = 2.8 * p["n_steps"] / 7540 * cells / 25600
    return f"{min(240, max(4, math.ceil(2.5 * est))):02d}:00:00", est


if __name__ == "__main__":
    rs = rows()
    for r in rs:
        r["walltime"], est = walltime(r)
        steps = "steady" if r["family"] == "T" else f"{plan(r)['n_steps']:6d}"
        print(f"{r['case']:30s} steps={steps:>6s} est={est:5.1f} h walltime={r['walltime']}")
    with open(HERE / "cases.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, lineterminator="\n")
        w.writeheader()
        w.writerows(rs)
    print(f"{len(rs)} cases -> cases.csv")
