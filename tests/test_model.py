"""Software tests (no Fluent needed):  uv run --no-project --with pytest pytest -q"""
import csv
import importlib.util
import math
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from trenchflow import mesh  # noqa: E402
from trenchflow.fluent import parse_slurm_time  # noqa: E402
from trenchflow.physics import (CASES_CSV, acoustic_scales, beam_direction, half_wave,  # noqa: E402
                                inflow_rates, inlet_v_expr, plan, sherwood)

BASELINE = ROOT.parent / "3d_nand" / "ansys_busnaina"   # sibling checkout, if present


def evaluate(expr, t=0.0, y=0.0, **named):
    """Evaluate one of our Fluent expression strings in Python (units stripped)."""
    py = re.sub(r"\[[^\]]*\]", "", expr)
    env = {"sin": math.sin, "abs": abs, "t": t, "y": y}
    env.update(named)
    return eval(py, {"__builtins__": {}}, env)  # noqa: S307 - our own strings


def case_rows():
    with open(CASES_CSV, newline="") as f:
        return list(csv.DictReader(f))


def row_by(name):
    return next(r for r in case_rows() if r["case"] == name)


# ---------------------------------------------------------------- mesh
@pytest.mark.parametrize("ar", [1, 5])
def test_mesh_checks_and_top_type(tmp_path, ar):
    for top in ("symmetry", "velocity-inlet"):
        out = tmp_path / f"m_{top}.msh"
        info = mesh.main(["--W", "1e-6", "--AR", str(ar), "--n-per-w", "8", "--top-type", top, "-o", str(out)])
        text = out.read_text()
        assert f"(45 (9 {top} top)())" in text
        assert info["n_cells"] == 15 * 64 + ar * 64       # channel 6W x 2.5W + cavity W x AR W


@pytest.mark.skipif(not (BASELINE / "make_trench_mesh.py").exists(), reason="baseline checkout not found")
def test_symmetry_mesh_identical_to_baseline_generator(tmp_path):
    spec = importlib.util.spec_from_file_location("orig_mesh", BASELINE / "make_trench_mesh.py")
    orig = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(orig)
    a, b = tmp_path / "new.msh", tmp_path / "orig.msh"
    mesh.main(["--W", "1e-6", "--AR", "1", "--n-per-w", "8", "-o", str(a)])
    orig.main(["--W", "1e-6", "--AR", "1", "--n-per-w", "8", "-o", str(b)])
    assert a.read_text() == b.read_text()


# ---------------------------------------------------------------- waveforms and boundary conditions
def test_half_wave_matches_paper_definition():
    T = 1.0
    for tt in [i / 97 for i in range(97)]:
        paper = 0.0 if tt < 0.5 else math.sin(2 * math.pi * (tt / T - 0.5))
        assert half_wave(2 * math.pi * tt / T) == pytest.approx(paper, abs=1e-12)


def test_case_a_expression_reproduces_paper_pulse():
    p = plan(row_by("A_ar1_osc200k"))
    T = 1 / p["freq"]
    for k in range(40):
        t = (k + 0.5) / 40 * T
        expect = p["Us"] + p["Up"] * (0.0 if t / T < 0.5 else math.sin(2 * math.pi * (t / T - 0.5)))
        assert evaluate(p["u_expr"], t=t) == pytest.approx(expect, rel=1e-6, abs=1e-9)


@pytest.mark.skipif(not (BASELINE / "run_case.py").exists(), reason="baseline checkout not found")
def test_case_a_plan_identical_to_baseline_fig4_f200k():
    spec = importlib.util.spec_from_file_location("orig_run", BASELINE / "run_case.py")
    orig = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(BASELINE))
    spec.loader.exec_module(orig)
    base = orig.plan(orig.load_case("fig4_f200k"), 20, 2000, 1.0)
    new = plan(row_by("A0_ar1_osc200k"))
    assert new["split_cavity"] is False and new["cavity_zones"] == ["cavity"]   # original mesh exactly
    split = plan(row_by("A_ar1_osc200k"))                                        # same physics, split zones
    for k in ("dt", "n_steps", "u_expr", "top_type", "left_type", "Lu_W", "Ld_W"):
        assert split[k] == new[k]
    assert split["cavity_zones"] == ["cavity", "cavity-lower"]
    for k in ("dt", "n_steps", "Us", "Up", "St"):
        assert new[k] == pytest.approx(base[k], rel=1e-12)
    assert new["u_expr"] == base["inlet_expr"]
    assert new["top_type"] == "symmetry" and new["v_expr"] is None


def test_transverse_modes():
    p = plan(row_by("B_ar1_steady_v0.5"))
    assert p["top_type"] == "velocity-inlet"
    assert evaluate(p["v_expr"]) == pytest.approx(0.075)
    assert p["phi_deg"] == pytest.approx(math.degrees(math.atan(0.5)))

    p = plan(row_by("Bf_ar1_osc200k_v1"))
    assert p["v_expr"].strip().endswith("u_in")
    assert evaluate(p["v_expr"], u_in=0.3) == pytest.approx(0.3)

    p = plan(row_by("C_ar1_osc200k_vp200k_ph180"))
    T = 1 / p["v_freq"]
    mean = sum(evaluate(p["v_expr"], t=(k + 0.5) / 2000 * T) for k in range(2000)) / 2000
    assert mean == pytest.approx(0.5 * 0.15, rel=1e-3)          # mean V = lambda * u_avg
    # anti-phase: V pulses while U is zero and vice versa
    pu = plan(row_by("A_ar1_osc200k"))
    for k in range(20):
        t = (k + 0.5) / 20 * T
        assert min(evaluate(pu["u_expr"], t=t), evaluate(p["v_expr"], t=t)) == pytest.approx(0, abs=1e-9)


def test_rotate_keeps_speed_and_sets_open_centred_domain():
    for ang in (0, 15, 30, 45, 60, 75, 90):
        p = plan(row_by(f"R_w1um_ar5_steady_a{ang:02d}"))
        u = evaluate(p["top_u"], u_in=0.15)
        v = evaluate(p["top_v"], u_in=0.15)
        assert math.hypot(u, v) == pytest.approx(0.15)                  # same speed at every angle
        assert math.degrees(math.atan2(-v, u)) == pytest.approx(ang)    # angle from the wafer surface
        assert v <= 0                                                   # toward the wafer
        assert p["left_type"] == "pressure-outlet" and p["split_cavity"]
        assert p["Lu_W"] == p["Ld_W"] == 3.0                            # trench centred
        assert p["cavity_zones"] == ["cavity", "cavity-lower"]
    assert evaluate(plan(row_by("R_w1um_ar5_steady_a90"))["top_u"], u_in=0.15) == 0.0


def test_rotate_mesh_has_open_left_and_split_trench(tmp_path):
    from trenchflow.physics import mesh_args
    p = plan(row_by("R_w1um_ar5_steady_a45"))
    p["n_per_w"] = 8
    out = tmp_path / "r.msh"
    info = mesh.main(mesh_args(p, out))
    text = out.read_text()
    for zone in ("(45 (7 pressure-outlet left)())", "(45 (9 velocity-inlet top)())",
                 "(45 (12 fluid cavity-lower)())", "(45 (13 interior cavity-split)())"):
        assert zone in text
    assert info["n_cavity_lower"] * 2 == info["n_cavity"]


def test_optimum_finder_recovers_parabola_peak():
    sys.path.insert(0, str(ROOT))
    import analyze
    angles = [0, 15, 30, 45, 60, 75, 90]
    best, est = analyze.optimum(angles, [-(a - 37.0) ** 2 for a in angles], maximize=True)
    assert best == 30 and est == pytest.approx(37.0)
    best, est = analyze.optimum(angles, [(a - 52.0) ** 2 for a in angles], maximize=False)
    assert best == 45 and est == pytest.approx(52.0)
    best, est = analyze.optimum(angles, list(range(7)), maximize=True)   # monotonic: optimum at the end
    assert best == 90 and est == 90


def test_left_inlet_profile_matches_wall_and_top():
    H = 2.5e-6
    expr = inlet_v_expr(H)
    assert evaluate(expr, y=0.0, v_amp=0.2) == pytest.approx(0.0)          # no penetration at the wafer
    assert evaluate(expr, y=H, v_amp=0.2) == pytest.approx(-0.2)           # matches the top inflow
    q_left, q_top = inflow_rates(0.15, 0.075, H, 6e-6)
    assert q_left == pytest.approx(0.15 * H) and q_top == pytest.approx(0.075 * 6e-6)


def test_case_table_consistent():
    rows = case_rows()
    names = [r["case"] for r in rows]
    assert len(names) == len(set(names))
    for r in rows:
        h, m, s = (int(x) for x in r["walltime"].split(":"))
        assert h <= 240
        if r["family"] == "T":                 # wafer-scale steady runs: planned by trenchflow.tank
            continue
        p = plan(r)
        assert p["n_steps"] > 0 and p["dt"] > 0
        assert (p["v_mode"] == "none") == (r["family"] in ("A", "A0"))
        if r["case"].endswith(("_a15", "_a30", "_a45", "_a60", "_a75")):   # angle grid matches its name
            assert p["phi_deg"] == pytest.approx(float(r["case"][-2:]), abs=1e-4)


# ---------------------------------------------------------------- acoustic conventions and scales
def test_angle_convention_from_wafer_surface():
    assert beam_direction(0) == pytest.approx((1.0, 0.0), abs=1e-12)       # grazing: along the wafer
    assert beam_direction(90) == pytest.approx((0.0, -1.0), abs=1e-12)     # normal: toward the wafer
    kx, ky = beam_direction(30)
    assert math.degrees(math.atan2(-ky, kx)) == pytest.approx(30)


def test_acoustic_scales_against_hand_values():
    s = acoustic_scales(172e3, 1e6)                    # ~1 W/cm2 plane wave
    assert s["intensity_W_cm2"] == pytest.approx(1.0, rel=0.01)
    assert s["stokes_layer"] == pytest.approx(5.66e-7, rel=0.01)
    assert s["wavelength"] == pytest.approx(1.48e-3)
    b = acoustic_scales(27e3, 1.002e6)                 # Paper B, 1 W input (Fig. 4)
    assert b["intensity_W_cm2"] == pytest.approx(0.0247, rel=0.02)
    g = acoustic_scales(100e3, 1e6, theta_deg=0)
    assert g["u_normal"] == pytest.approx(0, abs=1e-12) and g["u_tangential"] == pytest.approx(g["U0"])


# ---------------------------------------------------------------- metrics and utilities
def test_sherwood_of_exponential_decay():
    k, depth, W = 5000.0, 1e-6, 1e-6
    t = [i * 1e-5 for i in range(50)]
    sh = sherwood(t, [math.exp(-k * x) for x in t], depth, W)
    assert all(v == pytest.approx(depth * k * W / 6e-9) for _, v in sh)


def test_parse_slurm_time():
    assert parse_slurm_time("1-02:03:04") == 93784
    assert parse_slurm_time("05:30") == 330
    assert parse_slurm_time("UNLIMITED") is None


# ---------------------------------------------------------------- Case S (shear) and Case T (tank)
def test_graded_spacing_sums_exactly_and_keeps_core_uniform():
    h = 0.025
    for L, core in ((3.0, 1.0), (10.0, 1.0), (7.5, 2.5), (0.6, 1.0)):
        s = mesh.graded(L, h, core, 1.15, 0.5)
        assert sum(s) == pytest.approx(L, rel=1e-12)
        n_core = min(len(s), round(core / h))
        assert all(x == pytest.approx(h) for x in s[:n_core - 1])        # uniform next to the trench
        assert max(s) <= 0.5 + 1e-12
        assert all(b <= a * 1.15 * (1 + 1e-9) for a, b in zip(s, s[1:]))


def test_shear_case_plan_and_graded_mesh(tmp_path):
    from trenchflow.physics import GRADE, mesh_args
    p = plan(row_by("S_ar5_G1k"))
    W, H, G = p["W"], p["H"], 1e3
    assert p["v_mode"] == "shear" and H == pytest.approx(10 * W)
    assert evaluate(p["u_expr"], y=2e-6) == pytest.approx(G * 2e-6)          # u = G y at the inlet
    assert evaluate(p["top_u"]) == pytest.approx(G * H) and evaluate(p["top_v"]) == 0   # continuous at top
    assert p["dt"] <= p["h"] / (G * GRADE["core_y"] * W) * (1 + 1e-9)
    assert p["split_cavity"] and p["cavity_zones"] == ["cavity", "cavity-lower"]
    p["n_per_w"] = 8
    info = mesh.main(mesh_args(p, tmp_path / "s.msh"))      # runs the orientation/length/area checks
    assert info["H"] == pytest.approx(H) and info["grade"]["ratio"] == GRADE["ratio"]
    with pytest.raises(SystemExit):
        plan(dict(row_by("S_ar1_G1k"), mode="osc"))


def test_tank_plan_same_flow_and_aim_at_every_angle():
    from trenchflow.tank import plan_tank
    for ang in (15, 30, 45, 60, 75, 90):
        p = plan_tank(row_by(f"T_a{ang:02d}"))
        x0, x1 = p["slot_x"]
        assert -p["v_in"] * (x1 - x0) == pytest.approx(p["U"] * p["slot"])           # same pump flow
        assert math.hypot(p["u_in"], p["v_in"]) == pytest.approx(p["U"])
        assert math.degrees(math.atan2(-p["v_in"], p["u_in"])) == pytest.approx(ang)  # angle from wafer
        xc = 0.5 * (x0 + x1)
        assert xc + p["standoff"] * p["u_in"] / -p["v_in"] == pytest.approx(p["aim"], abs=1e-12)
    with pytest.raises(SystemExit):
        plan_tank(dict(row_by("T_a15"), angle_deg="5"))             # slot would not fit above the wafer


def test_tank_mesh_checks(tmp_path):
    from trenchflow.tank import mesh_only
    for ang in (15, 90):
        row = dict(row_by(f"T_a{ang:02d}"), refine="0.25")           # coarse, for speed
        p, nx, ny = mesh_only(row, tmp_path / f"t{ang}.msh")
        text = (tmp_path / f"t{ang}.msh").read_text()
        for zone in ("velocity-inlet nozzle", "wall wafer", "pressure-outlet open-top", "fluid tank"):
            assert zone in text


def test_wafer_interp_and_best_tilt():
    import analyze_wafer as aw
    tab = [(10.0, {"eta": 1.0}), (100.0, {"eta": 3.0}), (1000.0, {"eta": 7.0})]
    assert aw.interp(10 ** 1.5, tab, "eta") == (pytest.approx(2.0), False)
    assert aw.interp(1.0, tab, "eta") == (1.0, True)                  # below the calibration: clipped
    assert aw.interp(1e4, tab, "eta") == (7.0, True)
    s = aw.wafer_stats([0.0, 1.0, 2.0, 3.0], [5.0, 1.0, 5.0, 5.0])
    assert s["min"] == 1.0 and s["x_worst"] == 1.0 and s["mean"] == pytest.approx(4.0)
    table = {15: dict(min=2.0, mean=9.0), 45: dict(min=3.0, mean=4.0), 90: dict(min=3.0, mean=5.0)}
    assert aw.best_tilt(table, key="min") == 90                       # equal worst case: higher mean wins
    table = {15: dict(p5=6.0, mean=9.0), 45: dict(p5=7.0, mean=4.0)}
    assert aw.best_tilt(table) == 45                                  # default ranking: 5th percentile
