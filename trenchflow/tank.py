r"""Case T: wafer-scale flow in a submerged tank with ONE inlet and a tilted wafer (steady, 2D).

Purpose: find how a single inlet's flow is spread over the wafer for each tilt, as the near-wall shear
rate G(x) = |tau_wall(x)| / mu along the wafer. The trench-scale Case S runs give the cleaning of a
trench as a function of G, so together they give the cleaning of every trench on the wafer.

Geometry (wafer frame: tilting the wafer = turning the inlet; the tank is taken as large, so its walls
are open far-field boundaries):

      open-top (pressure outlet)   nozzle (velocity inlet, flow at angle phi)   open-top
   y=hs +----------------------------[=====]-----------------------------------+
        |                              \  jet, speed U, angle phi from the wafer |
  open- |                               \                                        | open-
  left  |                                \                                       | right
   y=0  +------------------------ wafer (wall, no slip) ---------------------------+
      x=-L/2                         aim point x = aim                         x=+L/2

The slot of width b sits on the top boundary, centred at x_c = aim - hs / tan(phi), so the jet axis hits
the wafer at x = aim. Its opening on the top boundary is b / sin(phi) long and the inlet velocity is
U (cos phi, -sin phi): the volume flow U * b is the same at every angle (same pump).

Model: steady RANS, k-omega SST (slot Reynolds number U b / nu ~ 1500: a planar jet at this Re is
unsteady/turbulent), water, wall cells graded to y+ < 1 on the wafer. Imposed flow only: no acoustics.
"""
import csv
import json
import math
import time

from .fluent import first_ok, launch, log, option, pick, reports, set_inlet_components, settings, setv
from .mesh import graded
from .physics import MU, NU, RHO

BC = {"interior": 2, "wall": 3, "pressure-outlet": 5, "velocity-inlet": 10}
ZONES = [  # (id, name, type)
    (3, "interior-tank", "interior"),
    (4, "wafer", "wall"),
    (5, "nozzle", "velocity-inlet"),
    (6, "open-top", "pressure-outlet"),
    (7, "open-left", "pressure-outlet"),
    (8, "open-right", "pressure-outlet"),
]
FLUID_ID = 2


# ---------------------------------------------------------------- plan
def _f(row, key, default):
    v = row.get(key)
    return float(v) if v not in (None, "") else default


def plan_tank(row):
    phi = _f(row, "angle_deg", float("nan"))
    if not 0 < phi <= 90:
        raise SystemExit("Case T needs 0 < angle_deg <= 90 (the jet must reach the wafer)")
    U = _f(row, "u_avg_m_s", 0.15)
    b = _f(row, "slot_mm", 10.0) * 1e-3
    hs = _f(row, "standoff_mm", 30.0) * 1e-3
    L = _f(row, "wafer_mm", 300.0) * 1e-3
    aim = _f(row, "aim_mm", 0.0) * 1e-3
    r = _f(row, "refine", 1.0)
    s, c = math.sin(math.radians(phi)), math.cos(math.radians(phi))
    c = 0.0 if abs(c) < 1e-12 else c
    xc = aim - (hs * c / s)
    b_top = b / s
    x0, x1 = xc - b_top / 2, xc + b_top / 2
    if x0 <= -L / 2 + 1e-3 or x1 >= L / 2 - 1e-3:
        raise SystemExit(f"slot [{x0 * 1e3:.1f}, {x1 * 1e3:.1f}] mm does not fit above the wafer: "
                         "steeper angle, smaller standoff or aim point closer to the centre")
    return dict(case=row["case"], phi_deg=phi, U=U, slot=b, standoff=hs, wafer=L, aim=aim, refine=r,
                slot_x=[x0, x1], u_in=U * c, v_in=-U * s, flow_per_depth=U * b,
                Re_slot=U * b / NU, dx=2.5e-4 / r, y1=1e-5 / r, ratio=1.1 ** (1 / r), hmax=5e-4 / r,
                max_iters=int(_f(row, "iters", 0) or 6000), res_tol=_f(row, "res_tol", 1e-6))


# ---------------------------------------------------------------- mesh
def _segment(a, b, dx):
    n = max(1, math.ceil((b - a) / dx - 1e-9))
    return [a + (b - a) * i / n for i in range(n)]


def build_tank(p):
    """Structured rectangle; x breaks at the slot edges, y graded from the wafer."""
    L, (x0, x1) = p["wafer"], p["slot_x"]
    xs = _segment(-L / 2, x0, p["dx"]) + _segment(x0, x1, p["dx"]) + _segment(x1, L / 2, p["dx"]) + [L / 2]
    hy = graded(p["standoff"], p["y1"], 0.0, p["ratio"], p["hmax"])
    ys = [sum(hy[:k]) for k in range(len(hy) + 1)]
    nx, ny = len(xs) - 1, len(ys) - 1

    def node(i, j):
        return j * (nx + 1) + i + 1

    def cell(i, j):
        return j * nx + i + 1 if 0 <= i < nx and 0 <= j < ny else 0

    zones = {zid: [] for zid, _, _ in ZONES}
    tol = 1e-9
    for j in range(ny + 1):                         # horizontal faces, left-to-right: c0 = cell above
        for i in range(nx):
            above, below = cell(i, j), cell(i, j - 1)
            if above:
                f = (node(i, j), node(i + 1, j), above, below)
            else:                                   # top boundary: reverse so the cell is on the left
                f = (node(i + 1, j), node(i, j), below, 0)
            if j == 0:
                zones[4].append(f)
            elif j == ny:
                xm = 0.5 * (xs[i] + xs[i + 1])
                zones[5 if x0 - tol < xm < x1 + tol else 6].append(f)
            else:
                zones[3].append(f)
    for i in range(nx + 1):                         # vertical faces, upward: c0 = cell on the left
        for j in range(ny):
            left, right = cell(i - 1, j), cell(i, j)
            if left:
                f = (node(i, j), node(i, j + 1), left, right)
            else:
                f = (node(i, j + 1), node(i, j), right, 0)
            zones[7 if i == 0 else 8 if i == nx else 3].append(f)
    nodes = [(x, y) for y in ys for x in xs]
    return nodes, xs, ys, zones


def check_tank(nodes, xs, ys, zones, p):
    nx, ny = len(xs) - 1, len(ys) - 1

    def centroid(c):
        j, i = divmod(c - 1, nx)
        return 0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1])

    refs = [0] * (nx * ny + 1)
    for flist in zones.values():
        for a, b, c0, c1 in flist:
            (ax, ay), (bx, by) = nodes[a - 1], nodes[b - 1]
            for c, sign in ((c0, 1), (c1, -1)):
                if c:
                    cx, cy = centroid(c)
                    if sign * ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax)) <= 0:
                        raise RuntimeError(f"bad orientation: face {a}-{b}, cell {c}")
                    refs[c] += 1
    if any(r != 4 for r in refs[1:]):
        raise RuntimeError("a cell does not have exactly 4 faces")

    def length(zid):
        return sum(math.dist(nodes[a - 1], nodes[b - 1]) for a, b, _, _ in zones[zid])
    L, hs, (x0, x1) = p["wafer"], p["standoff"], p["slot_x"]
    for zid, want in ((4, L), (5, x1 - x0), (6, L - (x1 - x0)), (7, hs), (8, hs)):
        if abs(length(zid) - want) > 1e-9 * L:
            raise RuntimeError(f"zone {zid} length {length(zid)} != {want}")


def write_tank_msh(path, nodes, xs, ys, zones, title):
    nn, nc = len(nodes), (len(xs) - 1) * (len(ys) - 1)
    order = [z for z in ZONES if zones[z[0]]]
    nf = sum(len(zones[z[0]]) for z in order)
    with open(path, "w", newline="\n") as f:
        f.write(f'(0 "{title}")\n(2 2)\n')
        f.write(f"(10 (0 1 {nn:x} 0 2))\n(12 (0 1 {nc:x} 0))\n(13 (0 1 {nf:x} 0))\n")
        f.write(f"(10 (1 1 {nn:x} 1 2)(\n")
        f.writelines(f"{x:.15e} {y:.15e}\n" for x, y in nodes)
        f.write("))\n")
        f.write(f"(12 ({FLUID_ID:x} 1 {nc:x} 1 3))\n")
        first = 1
        for zid, _, ztype in order:
            flist = zones[zid]
            last = first + len(flist) - 1
            f.write(f"(13 ({zid:x} {first:x} {last:x} {BC[ztype]:x} 2)(\n")
            f.writelines(f"{a:x} {b:x} {c0:x} {c1:x}\n" for a, b, c0, c1 in flist)
            f.write("))\n")
            first = last + 1
        f.write(f"(45 ({FLUID_ID} fluid tank)())\n")
        for zid, name, ztype in order:
            f.write(f"(45 ({zid} {ztype} {name})())\n")


def make_mesh(p, path):
    nodes, xs, ys, zones = build_tank(p)
    check_tank(nodes, xs, ys, zones, p)
    title = f"Tank: wafer {p['wafer'] * 1e3:g} mm, jet {p['phi_deg']:g} deg, refine {p['refine']:g}"
    write_tank_msh(path, nodes, xs, ys, zones, title)
    info = dict(n_cells=(len(xs) - 1) * (len(ys) - 1), cells_x=len(xs) - 1, cells_y=len(ys) - 1,
                first_cell_m=ys[1], largest_dy_m=max(b - a for a, b in zip(ys, ys[1:])))
    log("tank mesh:", json.dumps(info))
    return info


# ---------------------------------------------------------------- Fluent
def setup_tank(solver, msh, p):
    s = settings(solver)
    first_ok("read mesh", lambda: s.file.read_mesh(file_name=msh),
             lambda: s.file.read(file_type="mesh", file_name=msh), lambda: solver.tui.file.read_case(msh))
    try:
        s.mesh.check()
    except Exception as e:  # noqa: BLE001
        log("mesh check call failed (continuing):", e)
    vis = s.setup.models.viscous
    first_ok("k-omega SST",
             lambda: (setv(pick(vis, "model"), "k-omega"), setv(pick(vis, "k_omega_model"), "sst")),
             lambda: setv(pick(vis, "model"), "k-omega-sst"),
             lambda: solver.tui.define.models.viscous.kw_sst("yes"))
    fl = s.setup.materials.fluid
    name = list(fl.get_object_names())[0]          # the default fluid, given water's properties
    first_ok(f"water properties on {name}",
             lambda: (setv(fl[name].density.option, "constant"), setv(fl[name].density.value, RHO),
                      setv(fl[name].viscosity.option, "constant"), setv(fl[name].viscosity.value, MU)))
    ne = s.setup.named_expressions
    ne["u_jet"] = {"definition": f"{p['u_in']:.10g} [m/s]"}
    ne["v_jet"] = {"definition": f"{p['v_in']:.10g} [m/s]"}
    set_inlet_components(s, "nozzle", "u_jet", "v_jet")
    try:
        tb = s.setup.boundary_conditions.velocity_inlet["nozzle"].turbulence
        spec = pick(tb, "turbulence_specification", "turb_specification")
        setv(spec, option(spec, "Intensity and Hydraulic Diameter", "intensity-and-hydraulic-diameter"))
        setv(pick(tb, "turbulent_intensity", "turb_intensity"), 0.05)
        setv(pick(tb, "hydraulic_diameter", "turb_hydraulic_diam"), 2 * p["slot"])
        log("nozzle turbulence: 5 %, hydraulic diameter", 2 * p["slot"])
    except Exception as e:  # noqa: BLE001
        log("nozzle turbulence left at Fluent default:", e)
    try:
        setv(s.solution.methods.p_v_coupling.flow_scheme, "Coupled")
        log("pressure-velocity coupling: Coupled")
    except Exception as e:  # noqa: BLE001
        log("Coupled scheme not set, keeping default:", e)

    rd = s.solution.report_definitions.surface
    names = []
    for rname, rtype, surf, field in (("tau_mean", ("surface-areaavg", "area-weighted-avg"), "wafer", "wall-shear"),
                                      ("mdot_nozzle", ("surface-massflowrate", "mass-flow-rate"), "nozzle", None),
                                      ("mdot_open", ("surface-massflowrate", "mass-flow-rate"),
                                       ["open-top", "open-left", "open-right"], None)):
        rd[rname] = {}
        r = rd[rname]
        setv(r.report_type, option(r.report_type, *rtype))
        if field:
            setv(r.field, field)
        setv(pick(r, "surface_names", "surfaces"), surf if isinstance(surf, list) else [surf])
        names.append(rname)
    eqs = s.solution.monitor.residual.equations
    for n_ in eqs.get_object_names():
        try:
            setv(eqs[n_].absolute_criteria, p["res_tol"])
        except Exception as e:  # noqa: BLE001
            log(f"residual criterion {n_} not set: {e}")
    return names


def export_shear(solver, out):
    """Wall shear (signed x component and magnitude) and y+ along the wafer -> wafer_shear.csv."""
    s = settings(solver)
    raw = out / "wafer_raw.csv"
    funcs = ["x-wall-shear", "wall-shear", "y-plus"]
    first_ok("export wall shear",
             lambda: s.file.export.ascii(file_name=raw.name, surface_name_list=["wafer"], delimiter="comma",
                                         cell_func_domain=funcs, location="cell-center"),
             lambda: s.file.export.ascii(file_name=raw.name, surfaces_list=["wafer"], delimiter="comma",
                                         cell_func_domain=funcs, location="cell-center"),
             lambda: solver.tui.file.export.ascii(raw.name, "wafer", "()", "yes", *funcs, "()", "no"))
    with open(raw, newline="") as f:
        rows = list(csv.reader(f))
    head = [h.strip().lower() for h in rows[0]]

    def col(*keys):
        for k in keys:
            for i, h in enumerate(head):
                if h == k or h.startswith(k):
                    return i
        raise RuntimeError(f"column {keys} not in export header {head}")
    ix, itx, it, iy = col("x-coordinate"), col("x-wall-shear"), col("wall-shear"), col("y-plus")
    data = sorted((float(r[ix]), float(r[itx]), float(r[it]), float(r[iy])) for r in rows[1:] if r)
    with open(out / "wafer_shear.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_m", "tau_x_Pa", "tau_Pa", "y_plus", "shear_rate_1_s"])
        w.writerows((x, tx, t, yp, abs(tx) / MU) for x, tx, t, yp in data)
    return data


def run(case, row, procs, out, setup_only=False, plan_only=False, force=False):
    p = plan_tank(row)
    log("tank case", case, json.dumps(p, indent=1))
    if plan_only:
        return 0
    if (out / "summary.json").exists() and not force:
        raise SystemExit(f"{out} already has a finished run. Use --force to overwrite it.")
    out.mkdir(parents=True, exist_ok=True)
    msh = out / f"{case}.msh"
    p["mesh"] = make_mesh(p, msh)
    solver = launch(procs, out)
    t0 = time.time()
    try:
        s = settings(solver)
        names = setup_tank(solver, msh.name, p)
        ini = s.solution.initialization
        first_ok("hybrid initialization", lambda: ini.hybrid_initialize(), lambda: ini.standard_initialize())
        chunk, max_it = (20, 20) if setup_only else (200, p["max_iters"])
        done, hist, prev, converged = 0, [], None, False
        while done < max_it:
            first_ok(f"iterate {chunk}", lambda: s.solution.run_calculation.iterate(iter_count=chunk),
                     lambda: solver.tui.solve.iterate(chunk))
            done += chunk
            r = reports(solver, names)
            imb = abs(r["mdot_nozzle"] + r["mdot_open"]) / max(abs(r["mdot_nozzle"]), 1e-30)
            change = abs(r["tau_mean"] - prev) / max(abs(r["tau_mean"]), 1e-30) if prev is not None else 1.0
            hist.append(dict(iter=done, **r, mass_imbalance=imb, rel_change=change))
            log(f"iter {done}: mean wafer shear {r['tau_mean']:.4e} Pa  change {change:.1e}  "
                f"mass imbalance {imb:.1e}")
            prev = r["tau_mean"]
            if len(hist) >= 3 and all(h["rel_change"] < 1e-4 for h in hist[-2:]) and imb < 1e-4:
                converged = True
                break
        with open(out / "history.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(hist[0]))
            w.writeheader()
            w.writerows(hist)
        data = export_shear(solver, out)
        first_ok("write case/data", lambda: s.file.write(file_type="case-data", file_name=f"{case}-final.cas.h5"),
                 lambda: solver.tui.file.write_case_data(f"{case}-final.cas.h5"))
        G = [abs(tx) / MU for _, tx, _, _ in data]
        p.update(converged=converged, iterations=done, setup_only=setup_only, wall_s=time.time() - t0,
                 n_wall_faces=len(data), y_plus_max=max(yp for *_, yp in data),
                 shear_rate_min=min(G), shear_rate_max=max(G), tau_mean=hist[-1]["tau_mean"],
                 mass_imbalance=hist[-1]["mass_imbalance"])
        if not converged and not setup_only:
            log("WARNING: mean wall shear still changing; the steady jet may be unsteady (see history.csv)")
        (out / "summary.json").write_text(json.dumps(p, indent=1))
        log("done:", json.dumps({k: p[k] for k in ("converged", "iterations", "y_plus_max",
                                                    "shear_rate_min", "shear_rate_max")}))
    finally:
        solver.exit()
    return 0


def mesh_only(row, path):
    """For tests / inspection: build and check the mesh without Fluent."""
    p = plan_tank(row)
    nodes, xs, ys, zones = build_tank(p)
    check_tank(nodes, xs, ys, zones, p)
    write_tank_msh(path, nodes, xs, ys, zones, "test")
    return p, len(xs) - 1, len(ys) - 1
