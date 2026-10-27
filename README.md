# 3d_nand_2: transverse trench exchange and megasonic incidence angle

This project extends the Lin, Busnaina & Suni (2002) trench-cleaning model, reproduced in the sibling project
[`../3d_nand`](https://github.com/arjun2022as/3d_nand) (`ansys_busnaina/`), to study whether transverse liquid
motion into and out of the trench, and the angle of an incoming megasonic wave, improve contaminant removal.

The original project is not modified. Case A here re-runs its baseline with the same model and settings.

## Status

| Case | What | Model level | State |
|---|---|---|---|
| A | Original baseline: rinse U(t), v = 0, symmetry top | Paper A | implemented; HPC runs pending |
| B, Bf | + imposed transverse inflow V through the top (steady; or V = λU, fixed approach angle) | Level 1: **imposed flow** | implemented; HPC runs pending |
| C | + pulsed V with its own frequency/phase | Level 1: **imposed flow** | implemented; HPC runs pending |
| D | Megasonic incidence angle θ → acoustic streaming force → flow | Level 3: perturbation acoustics | **designed, not implemented** (docs/formulation.md §7) |
| E | Comparison and verification figures | – | `analyze.py` (tested on synthetic data) |

**Cases B/C impose a transverse flow; they don't calculate megasonic streaming.** No acoustic results exist yet.

## Documentation

- [docs/formulation.md](docs/formulation.md): equations, boundary conditions, U/V, the two angle definitions (φ for the imposed flow, θ for the megasonic wave), formula-change audit, acoustic scales, Case D design
- [docs/literature.md](docs/literature.md): what each source states, abstract-only items, deductions and hypotheses
- [docs/verification.md](docs/verification.md): software tests, numerical verification, reproduction, and the absent experimental validation

## Layout

| File | Purpose |
|---|---|
| `trenchflow/mesh.py` | 2D Fluent mesh, copied from the baseline; `--top-type symmetry/velocity-inlet` |
| `trenchflow/physics.py` | case plans, U/V expressions, angle convention, acoustic scales, Sherwood number (pure Python) |
| `trenchflow/fluent.py` | Fluent setup (baseline code + velocity-component inlets + balance reports) |
| `make_cases.py` → `cases.csv` | the 40-case matrix (edit the script, not the CSV) |
| `run_case.py`, `run_case.slurm`, `submit_cases.sh` | run one case / SLURM job with auto-continuation / submit by pattern |
| `check_runs.sh`, `iter_stats.sh` | run status; iterations per time step |
| `analyze.py` | figures + CSV data + `results/analysis/verification.md` |
| `tests/` | software tests (no Fluent): `uv run --no-project --with pytest pytest -q` |

## Running on UA HPC (Puma)

This project uses the Python environment `~/venvs/pyfluent241` created for `3d_nand`. If it's missing, run
`bash ../3d_nand/ansys_busnaina/setup_env.sh` once in an `interactive` session.

```bash
cd ~ && git clone https://github.com/arjun2022as/3d_nand_2.git     # first time; later: cd ~/3d_nand_2 && git pull
cd ~/3d_nand_2

# 1. smoke test (5 time steps, writes results/<case>_test, ~5 min)
bash submit_cases.sh B_ar1_steady_v0.5 --setup-only
grep -E "ok$|failed|Error|components|residual|C/C0" $(ls -t logs/*.out | head -1)

# 2. production
bash submit_cases.sh '^A_'          # baselines (4 jobs)
bash submit_cases.sh '^B_'          # steady transverse inflow (24 jobs)
bash submit_cases.sh '^(Bf|C)_'     # fixed-angle and pulsed transverse inflow (10 jobs)
bash submit_cases.sh '^V_'          # time-step and mesh checks (2 jobs)

# 3. status, then figures + verification report (compute node, for matplotlib)
bash check_runs.sh
interactive -a krishna -t 00:30:00
cd ~/3d_nand_2 && source ~/venvs/pyfluent241/bin/activate && python analyze.py && exit
```

- **Resources:** each job uses 4 cores and 16 GB. Estimated run times are 1–7 h, plus 8 h and 33 h for the two
  verification cases, scaled from the measured baseline cost.
- **Outputs:** logs go to `logs/<case>_<jobid>.out/.err`. Results go to `results/<case>/`: history.csv, reports.out, summary.json, `.cas.h5`/`.dat.h5` (t0, 10 snapshots, final).
- **Restarts:** jobs that reach their time limit checkpoint and resubmit themselves. A stopped chain shows as
  `STOPPED` in `check_runs.sh` with the resume command. Finished cases refuse to rerun without `--force`.
- **Note:** Fig. 2 jobs of `../3d_nand` may still be running; this project doesn't touch that directory.

## Parameters and their sources

| Parameter | Value | Source |
|---|---|---|
| ρ, μ | 998.2 kg/m³, 1.003e-3 Pa·s | water |
| D | 6e-9 m²/s | back-calculated from Paper A, Table 1 |
| W, AR | 1 µm; 1 and 5 | Paper A Fig. 4 geometry; AR 5 for high aspect ratio |
| u_avg, pulses | 0.15 m/s; Us = 0, Up = π·u_avg; 200 kHz | Paper A (Up from its Fig. 4 St values) |
| V/u_avg | 0–2 | sensitivity range (no source) |
| Megasonic f, pressure | 1 MHz; 27–69 kPa, plus 172 kPa (1 W/cm²) as an upper bound | Paper B (Case D) |
| θ sweep | 0–90° in 15° steps, constant incident intensity | Case D plan |
