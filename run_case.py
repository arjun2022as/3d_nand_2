#!/usr/bin/env python3
"""Run one case of the transverse-exchange study (Fluent 2024 R1 via PyFluent).

Case A (v_mode = none) is the Lin, Busnaina & Suni (2002) baseline exactly as in
3d_nand/ansys_busnaina. Cases B/C add an imposed transverse inflow V(t) through the top
boundary (see trenchflow/physics.py and docs/formulation.md). Case S drives the trench with a
near-wall shear flow; Case T (wafer-scale tank, trenchflow/tank.py) is dispatched to its own runner.

    python run_case.py --case B_ar1_steady_l0.5 --procs 4
    python run_case.py --case B_ar1_steady_l0.5 --setup-only     # 5 steps into results/<case>_test
    python run_case.py --case B_ar1_steady_l0.5 --plan-only      # time-step plan, no Fluent
"""
import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from trenchflow import mesh  # noqa: E402
from trenchflow.fluent import (EXIT_CONTINUE, checkpoint, initialize, iterate, launch, log,  # noqa: E402
                               parse_slurm_time, read_case_data, reports, setup, setv,
                               settings, truncate_history, write)
from trenchflow.physics import CASES_CSV, load_case, mesh_args, plan  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, help="row name in cases.csv")
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--out", default=None, help="output folder (default results/<case>)")
    ap.add_argument("--steps-per-period", type=int, default=20)
    ap.add_argument("--min-steps", type=int, default=2000)
    ap.add_argument("--cfl", type=float, default=1.0)
    ap.add_argument("--iters", type=int, default=20)
    ap.add_argument("--res-tol", type=float, default=1e-6)
    ap.add_argument("--log-points", type=int, default=200)
    ap.add_argument("--snapshots", type=int, default=10)
    ap.add_argument("--setup-only", action="store_true")
    ap.add_argument("--plan-only", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--time-left", default="")
    ap.add_argument("--margin-min", type=float, default=30)
    ap.add_argument("--checkpoint-hours", type=float, default=6)
    a = ap.parse_args()

    row = load_case(a.case, CASES_CSV)
    if row["family"] == "T":            # wafer-scale tank flow: steady RANS, its own setup
        from trenchflow import tank
        out = Path(a.out or HERE / "results" / (f"{a.case}_test" if a.setup_only else a.case)).resolve()
        sys.exit(tank.run(a.case, row, a.procs, out, a.setup_only, a.plan_only, a.force))

    p = plan(row, a.steps_per_period, a.min_steps, a.cfl, a.res_tol)
    a.iters = p["iters"] or a.iters
    p["iters"] = a.iters
    log("case", a.case, json.dumps(p, indent=1))
    if a.plan_only:
        return

    out = Path(a.out or HERE / "results" / (f"{a.case}_test" if a.setup_only else a.case)).resolve()
    ckpt = out / "checkpoint.json"
    state = None
    if a.resume:
        if not ckpt.exists():
            raise SystemExit(f"--resume: no checkpoint in {out}")
        state = json.loads(ckpt.read_text())
    else:
        if (out / "summary.json").exists() and not a.force:
            raise SystemExit(f"{out} already has a finished run. Use --force to overwrite it.")
        if ckpt.exists() and not a.force:
            raise SystemExit(f"{out} has an unfinished run. Use --resume or --force.")
        for old in out.glob("checkpoint*"):
            old.unlink()
    time_left = parse_slurm_time(a.time_left)
    log(f"time left in this job: {a.time_left or 'unknown'}")

    out.mkdir(parents=True, exist_ok=True)
    msh = out / f"{a.case}.msh"
    if not state:
        p["mesh"] = mesh.main(mesh_args(p, msh))

    solver = launch(a.procs, out)
    t0 = time.time()
    end_time = t0 + time_left if time_left else None
    exit_code, c = 0, None
    try:
        try:
            solver.tui.file.confirm_overwrite("no")
        except Exception:  # noqa: BLE001
            pass
        hist = out / "history.csv"
        if state:
            done, names = state["steps_done"], state["report_names"]
            log(f"resuming from checkpoint at step {done}/{state['n_total']} (part {state['part'] + 1})")
            read_case_data(solver, "checkpoint.cas.h5")
            try:
                rf = settings(solver).solution.monitor.report_files["reports_rfile"]
                setv(rf.file_name, f"reports_part{state['part'] + 1}.out")
            except Exception as e:  # noqa: BLE001
                log("report file name unchanged:", e)
            truncate_history(hist, done)
            mode = "a"
        else:
            names = setup(solver, msh.name, p)
            initialize(solver, p["cavity_zones"])
            write(solver, "case-data", f"{a.case}-t0.cas.h5")
            done, mode = 0, "w"

        n_total = 5 if a.setup_only else p["n_steps"]
        chunk = max(1, n_total // a.log_points)
        snap_every = max(1, n_total // max(1, a.snapshots))
        next_snap = (done // snap_every + 1) * snap_every
        next_ckpt = time.time() + a.checkpoint_hours * 3600
        part = state["part"] + 1 if state else 1
        wall_prev = state["wall_s"] if state else 0.0
        steps_this_job = 0
        with open(hist, mode, newline="") as f:
            wr = csv.writer(f)
            if mode == "w":
                wr.writerow(["step", "time_s", "cleaning_efficiency"] + names)
                r0 = reports(solver, names)
                wr.writerow([0, 0.0, 1.0 - r0["c_cavity"]] + [r0[n] for n in names])
            while done < n_total:
                n = min(chunk, n_total - done)
                t_chunk = time.time()
                iterate(solver, n, a.iters)
                t_chunk = time.time() - t_chunk
                done += n
                steps_this_job += n
                r = reports(solver, names)
                c = r["c_cavity"]
                wr.writerow([done, done * p["dt"], 1.0 - c] + [r[k] for k in names])
                f.flush()
                if done >= next_snap and done < n_total:
                    write(solver, "data", f"{a.case}-{done:08d}.dat.h5")
                    next_snap += snap_every
                rate = (time.time() - t0) / steps_this_job
                imbalance = sum(v for k, v in r.items() if k.startswith("mdot_"))
                log(f"step {done}/{n_total}  t={done * p['dt']:.4e} s  C/C0={c:.6f}  "
                    f"Y in [{r['y_min']:.2e}, {r['y_max']:.4f}]  mass imbalance {imbalance:.2e} kg/s  "
                    f"ETA {rate * (n_total - done) / 3600:.2f} h")
                if done >= n_total or a.setup_only:
                    continue
                out_of_time = end_time and time.time() + 1.5 * t_chunk + a.margin_min * 60 > end_time
                if out_of_time or time.time() > next_ckpt:
                    checkpoint(solver, out, dict(case=a.case, steps_done=done, n_total=n_total, part=part,
                                                 wall_s=wall_prev + time.time() - t0, report_names=names))
                    next_ckpt = time.time() + a.checkpoint_hours * 3600
                if out_of_time:
                    log(f"time limit near: stopped at step {done}/{n_total}; a continuation job will resume")
                    exit_code = EXIT_CONTINUE
                    break

        if exit_code == 0:
            write(solver, "case-data", f"{a.case}-final.cas.h5")
            p.update(c_final=c, cleaning_efficiency=1.0 - c, steps_run=done, parts=part,
                     setup_only=a.setup_only, wall_s=wall_prev + time.time() - t0)
            (out / "summary.json").write_text(json.dumps(p, indent=1))
            for ck in out.glob("checkpoint*"):
                ck.unlink()
            log("done:", json.dumps({k: p[k] for k in ("c_final", "cleaning_efficiency", "steps_run")}))
    finally:
        solver.exit()
    if exit_code:
        sys.exit(exit_code)


if __name__ == "__main__":
    main()
