# Verification plan and status

The four kinds of evidence are kept separate.

## 1. Software tests (run locally: `uv run --no-project --with pytest pytest -q`)

13 tests, all passing on 2026-10-27:
- **Mesh:** zone types and cell counts for AR 1 and AR 5. The symmetry mesh is **byte-identical** to the
  baseline generator in ../3d_nand.
- **Waveform:** g(s) equals Paper A's pulse; Case A's Fluent expression evaluates to the paper's U(t).
- **Baseline plan:** Case A's time step, step count and U expression are **identical** to the baseline's `fig4_f200k`.
- **Transverse BCs:** steady, follow and pulse values; pulse mean = λ·u_avg; anti-phase pulses never overlap
  the rinse pulses; left-inlet v is 0 at the wafer and −V at the top corner.
- **Case table:** unique names, valid plans, time limits ≤ 240 h, and family A ⇔ no transverse flow.
- **Angle convention:** θ from the wafer surface (0° grazing, 90° normal).
- **Acoustic scales:** match hand values; Paper B's 27 kPa gives 0.025 W/cm².
- **Utilities:** Sherwood number of an exponential decay is constant and correct; SLURM time parsing.

These don't run Fluent. The Fluent calls (velocity components, new report types) are first exercised
by the HPC smoke test.

## 2. Numerical verification (HPC runs; `analyze.py` writes results/analysis/verification.md)

1. Baseline recovery: `A_ar1_osc200k` vs ../3d_nand `fig4_f200k`. Same model and plan, so the difference
   should be at round-off/solver-tolerance level.
2. Zero transverse amplitude: `B_*_v0` vs `A_*`. This isolates the top boundary changing from symmetry
   to a u = U velocity inlet. A difference of a few % or less is expected.
3. Fluid mass balance: |Σṁ| / ṁ_outlet at every sample.
4. Concentration bounds: Y_min ≥ 0 and Y_max ≤ 1 within ~1e-3.
5. Contaminant balance: inventory lost from the domain vs time-integrated outflow (approximate; sampled).
6. Time-step and mesh sensitivity on `B_ar5_osc200k_v0.5`. The baseline study gave ~10 % at 40 cells/W.
7. Per-step convergence: `bash iter_stats.sh` (residuals 1e-6, ≤ 20 iterations/step, as in the baseline).

## 3. Reproduction of Paper A

Inherited from the baseline project (../3d_nand):
- Fig. 3 optimum St ≈ 0.133 reproduced in all four geometries.
- Fig. 4 weak frequency effect reproduced; less total removal than the paper.

## 4. Validation against experiments

**None.** No independent measurements of dissolved-contaminant removal from 1 µm trenches with transverse
or angled megasonic flow were found. Paper B's particle removal efficiency is for particles on a **bare**
wafer, a different quantity, and can't validate this model.
