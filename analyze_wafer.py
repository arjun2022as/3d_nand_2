#!/usr/bin/env python3
"""Best wafer tilt with ONE inlet: wafer-scale flow (Case T) x trench calibration (Case S).

    python analyze_wafer.py        # needs matplotlib (pyfluent241 venv); no Fluent

1. Case S gives, per aspect ratio, the trench cleaning at t_end and the Sherwood number as functions of
   the local near-wall shear rate G (log-linear interpolation between the sampled G).
2. Case T gives G(x) = |tau_x(x)| / mu along the wafer for each jet angle phi (wafer tilt).
3. Every point on the wafer gets the cleaning of a trench sitting there. The best tilt maximises the
   5th percentile: the cleaning reached by 95 % of the wafer (a wafer is only as clean as its dirtiest
   dies). The single worst point is reported too, but it is the jet stagnation point, where the shear
   is zero on one face, so it depends on the mesh.
Writes results/analysis/wafer_*.png, the CSV behind each figure and wafer_verification.md.
Imposed flow only: no megasonic acoustics.
"""
import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyze import angle_metrics, cases, run  # noqa: E402

RES = HERE / "results"
OUT = RES / "analysis"
MAIN_REFINE = ("", "1")


# ---------------------------------------------------------------- calibration (Case S)
def calibration(cs):
    """{AR: sorted [(G, metrics)]} from the finished S runs."""
    cal = {}
    for name, r in cs.items():
        if r["family"] != "S":
            continue
        h, p = run(name)
        if h:
            cal.setdefault(int(float(r["AR"])), []).append((float(r["shear_rate"]), angle_metrics(h, p)))
    return {ar: sorted(v, key=lambda t: t[0]) for ar, v in cal.items()}


def interp(G, table, key):
    """Metric at shear rate G, linear in log10(G) between calibration points; clipped at the ends.
    Returns (value, clipped) so callers can report how much of the wafer lies outside the calibration."""
    pts = [(math.log10(g), m[key]) for g, m in table if g > 0 and m[key] == m[key]]
    if not pts:
        return float("nan"), True
    lg = math.log10(max(G, 1e-30))
    if lg <= pts[0][0]:
        return pts[0][1], lg < pts[0][0]
    if lg >= pts[-1][0]:
        return pts[-1][1], lg > pts[-1][0]
    for (a, va), (b, vb) in zip(pts, pts[1:]):
        if a <= lg <= b:
            return va + (vb - va) * (lg - a) / (b - a), False
    return pts[-1][1], True


# ---------------------------------------------------------------- wafer-scale flow (Case T)
def wafer_profiles(cs):
    """{phi: (x list, G list, summary)} for the main-mesh T runs; mesh-check runs separately."""
    main, extra = {}, {}
    for name, r in cs.items():
        if r["family"] != "T":
            continue
        f, s = RES / name / "wafer_shear.csv", RES / name / "summary.json"
        if not (f.exists() and s.exists()):
            continue
        with open(f, newline="") as fh:
            rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(fh)]
        prof = ([q["x_m"] for q in rows], [q["shear_rate_1_s"] for q in rows], json.loads(s.read_text()))
        if r["refine"] in MAIN_REFINE:
            main[float(r["angle_deg"])] = prof
        else:
            extra[name] = prof
    return main, extra


def weights(x):
    """Wafer length each sample represents (samples are face centres, not evenly spaced)."""
    n = len(x)
    if n == 1:
        return [1.0]
    w = [0.5 * (x[min(i + 1, n - 1)] - x[max(i - 1, 0)]) for i in range(n)]
    w[0], w[-1] = x[1] - x[0], x[-1] - x[-2]
    return w


def wafer_stats(x, vals):
    w = weights(x)
    tot = sum(w)
    srt = sorted(zip(vals, w))
    acc, p5 = 0.0, srt[0][0]
    for v, wi in srt:
        acc += wi
        if acc >= 0.05 * tot:
            p5 = v
            break
    return dict(min=srt[0][0], p5=p5, mean=sum(v * wi for v, wi in zip(vals, w)) / tot,
                x_worst=x[min(range(len(vals)), key=lambda i: vals[i])])


def best_tilt(table, key="p5", tie="mean"):
    """phi with the best 5th percentile (or `key`); ties (within 0.01 points) broken by the mean."""
    top = max(v[key] for v in table.values())
    cands = [phi for phi, v in table.items() if v[key] >= top - 0.01]
    return max(cands, key=lambda phi: table[phi][tie])


# ---------------------------------------------------------------- figures
def figures(cal, prof, plt):
    rows_cal = []
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ar, tab in sorted(cal.items()):
        G = [g for g, _ in tab]
        for ax, key, lab in zip(axes, ("eta", "eta_deep", "sh"),
                                ("Trench cleaned at t_end (%)", "Lower half cleaned (%)", "Mean Sherwood number")):
            ax.plot(G, [m[key] for _, m in tab], "o-", label=f"AR {ar}")
            ax.set(xscale="log", xlabel="Near-wall shear rate G (1/s)", ylabel=lab)
        rows_cal += [(ar, g, m["eta"], m["eta_deep"], m["sh"]) for g, m in tab]
    for ar, style in ((1, ":"), (5, "--")):
        h, p = run(f"A_ar{ar}_steady")
        if h:
            m = angle_metrics(h, p)
            for ax, key in zip(axes, ("eta", "eta_deep", "sh")):
                ax.axhline(m[key], color="0.4", ls=style, lw=1, label=f"Paper A model, AR {ar} (steady U)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Trench calibration (Case S): cleaning of a 1 um trench vs the shear rate of the flow "
                 "over it, after 0.4 ms. Imposed flow, not acoustics.", fontsize=10)
    save(fig, "wafer_calibration", ["AR", "shear_rate_1_s", "eta", "eta_deep", "sh"], rows_cal)

    fig, ax = plt.subplots(figsize=(8, 4.6))
    rows = []
    for phi, (x, G, _) in sorted(prof.items()):
        ax.plot([1e3 * v for v in x], G, label=f"jet at {phi:g} deg")
        rows += [(phi, xi, gi) for xi, gi in zip(x, G)]
    for ar, tab in cal.items():
        for g in (tab[0][0], tab[-1][0]):
            ax.axhline(g, color="0.7", lw=0.6)
    ax.set(yscale="log", xlabel="Position on wafer (mm), aim point at 0", ylabel="Near-wall shear rate G (1/s)",
           title="Wafer-scale flow (Case T): shear rate along the wafer; grey = calibrated range")
    ax.legend(fontsize=8)
    save(fig, "wafer_shear_profiles", ["phi_deg", "x_m", "shear_rate_1_s"], rows)


def wafer_maps(cal, prof, plt):
    """Cleaning map along the wafer and the best tilt, per aspect ratio."""
    summary, lines = [], []
    for ar, tab in sorted(cal.items()):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4.6))
        stats, rows = {}, []
        for phi, (x, G, _) in sorted(prof.items()):
            ev = [interp(g, tab, "eta") for g in G]
            eta = [v for v, _ in ev]
            clipped = sum(w for (v, c), w in zip(ev, weights(x)) if c) / sum(weights(x))
            stats[phi] = dict(wafer_stats(x, eta), clipped=clipped)
            ax1.plot([1e3 * v for v in x], eta, label=f"{phi:g} deg")
            rows += [(phi, xi, gi, e) for xi, gi, e in zip(x, G, eta)]
        if not stats:
            plt.close(fig)
            continue
        best = best_tilt(stats)
        phis = sorted(stats)
        wbar = 4
        for k, (key, lab) in enumerate((("min", "worst trench"), ("p5", "5th percentile"), ("mean", "mean"))):
            ax2.bar([p + (k - 1) * wbar for p in phis], [stats[p][key] for p in phis], wbar, label=lab)
        ax2.axvline(best, color="k", ls=":", lw=1)
        ax1.set(xlabel="Position on wafer (mm)", ylabel="Trench cleaned at t_end (%)",
                title=f"AR {ar}: cleaning of a trench at each point on the wafer")
        ax2.set(xlabel="Jet angle to the wafer, phi (deg)", ylabel="Trench cleaned at t_end (%)",
                xticks=phis, title=f"AR {ar}: best tilt (95 % of wafer cleanest) = {best:g} deg")
        ax1.legend(fontsize=8, title="jet angle")
        ax2.legend(fontsize=8)
        save(fig, f"wafer_cleaning_ar{ar}", ["phi_deg", "x_m", "shear_rate_1_s", "eta"], rows)
        for phi in phis:
            v = stats[phi]
            summary.append((ar, phi, 90 - phi, v["min"], v["p5"], v["mean"], v["x_worst"], v["clipped"],
                            phi == best))
        lines.append(f"AR {ar}: best jet angle {best:g} deg (wafer tilted {90 - best:g} deg from facing the "
                     f"inlet); 95 % of the wafer >= {stats[best]['p5']:.2f} %, mean {stats[best]['mean']:.2f} %, "
                     f"worst point {stats[best]['min']:.2f} %")
    if summary:
        with open(OUT / "wafer_best_tilt.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["AR", "phi_deg", "tilt_from_facing_deg", "eta_worst", "eta_p5", "eta_mean",
                        "x_worst_m", "fraction_outside_calibration", "best"])
            w.writerows(summary)
        print("\nbest tilt with one inlet (highest 5th percentile of trench cleaning):")
        for ln in lines:
            print("  " + ln)
    return lines


def verification(prof, extra, cal):
    lines = ["# Wafer-scale verification", ""]
    for phi, (x, G, s) in sorted(prof.items()):
        lines.append(f"- jet {phi:g} deg: converged {s.get('converged')}, {s.get('iterations')} iterations, "
                     f"max y+ {s.get('y_plus_max', float('nan')):.2f}, mass imbalance "
                     f"{s.get('mass_imbalance', float('nan')):.1e}, G {s.get('shear_rate_min', 0):.3g} .. "
                     f"{s.get('shear_rate_max', 0):.3g} 1/s")
    for name, (x2, G2, _) in extra.items():
        phi = float(name.split("_a")[1][:2])
        if phi in prof:
            x1, G1, _ = prof[phi]
            diffs = []
            for xi, gi in zip(x1, G1):          # nearest refined sample
                j = min(range(len(x2)), key=lambda k: abs(x2[k] - xi))
                if gi > 0:
                    diffs.append(abs(G2[j] - gi) / gi)
            diffs.sort()
            if diffs:
                lines.append(f"- mesh check {name} vs jet {phi:g} deg: median relative change in G "
                             f"{100 * diffs[len(diffs) // 2]:.1f} %, 90th percentile "
                             f"{100 * diffs[int(0.9 * (len(diffs) - 1))]:.1f} %")
    for ar, tab in sorted(cal.items()):
        lines.append(f"- calibration AR {ar}: G = " + ", ".join(f"{g:g}" for g, _ in tab) + " 1/s")
    (OUT / "wafer_verification.md").write_text("\n".join(lines) + "\n")
    print("wrote", OUT / "wafer_verification.md")


def save(fig, name, header, rows):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=150)
    with open(OUT / f"{name}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print("wrote", OUT / f"{name}.png", "+ .csv")
    fig.clf()


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)
    cs = cases()
    cal = calibration(cs)
    prof, extra = wafer_profiles(cs)
    print(f"calibration runs: {sum(len(v) for v in cal.values())}, wafer-scale runs: {len(prof)} (+{len(extra)} checks)")
    if cal or prof:
        figures(cal, prof, plt)
    if cal and prof:
        wafer_maps(cal, prof, plt)
    verification(prof, extra, cal)


if __name__ == "__main__":
    main()
