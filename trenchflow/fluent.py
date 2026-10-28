"""Fluent 2024 R1 (PyFluent 0.37) setup and run helpers.

Material, model, initialisation, time-stepping and checkpoint code is taken unchanged from the
converged baseline 3d_nand/ansys_busnaina/run_case.py (commit 141e41b). New here:
  * velocity-component inlets (left inlet and top) for the transverse flow of Cases B/C,
  * extra reports for verification: domain contaminant, concentration bounds, mass flow rates.
"""
import csv
import json
import os
import time

from .physics import DIFF, MU, RHO, inlet_v_expr

TRACER = "h2o"   # species of mixture-template used as the contaminant (props overridden to water)


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


# ---------------------------------------------------------------- settings helpers
def pick(obj, *names):
    """Return the first existing child among `names` (settings names differ between versions)."""
    for n in names:
        try:
            return getattr(obj, n)
        except (AttributeError, RuntimeError):
            continue
    raise AttributeError(f"none of {names} under {getattr(obj, 'path', obj)}; "
                         f"available: {getattr(obj, 'child_names', '?')}")


def setv(obj, value):
    obj.set_state(value)
    return obj


def option(obj, *candidates):
    """First candidate that this option object allows (raises with the allowed list)."""
    allowed = list(obj.allowed_values())
    for c in candidates:
        if c in allowed:
            return c
    raise ValueError(f"none of {candidates} allowed; allowed: {allowed}")


def first_ok(label, *attempts):
    """Run callables in order until one succeeds; raise with all errors otherwise."""
    errors = []
    for fn in attempts:
        try:
            out = fn()
            log(f"{label}: ok")
            return out
        except Exception as e:  # noqa: BLE001 - want every failure reported
            errors.append(f"{type(e).__name__}: {e}")
    raise RuntimeError(f"{label} failed:\n  " + "\n  ".join(errors))


def scalar(res):
    """Pull the first number out of a report_definitions.compute() result."""
    if isinstance(res, (int, float)):
        return float(res)
    items = res.values() if isinstance(res, dict) else res if isinstance(res, (list, tuple)) else []
    for v in items:
        try:
            return scalar(v)
        except ValueError:
            continue
    raise ValueError(f"no number in {res!r}")


def settings(solver):
    return getattr(solver, "settings", solver)


# ---------------------------------------------------------------- launch
def launch(procs, workdir):
    import ansys.fluent.core as pyfluent
    log(f"PyFluent {pyfluent.__version__}, launching Fluent 2D double precision on {procs} cores")
    # Inside a SLURM job PyFluent would ssh-spawn Fluent processes to the same node; hide SLURM.
    for k in [k for k in os.environ if k.startswith("SLURM_")]:
        os.environ.pop(k)
    return pyfluent.launch_fluent(
        product_version=pyfluent.FluentVersion.v241,
        dimension=pyfluent.Dimension.TWO,
        precision=pyfluent.Precision.DOUBLE,
        processor_count=procs,
        ui_mode=pyfluent.UIMode.NO_GUI_OR_GRAPHICS,
        cwd=str(workdir),
        start_timeout=900,
    )


# ---------------------------------------------------------------- boundary conditions
def set_inlet_magnitude(s, name, expr_name):
    """Baseline (Case A): velocity magnitude normal to the boundary, exactly as the original model."""
    vi = s.setup.boundary_conditions.velocity_inlet[name]
    mom = vi.momentum if "momentum" in getattr(vi, "child_names", []) else vi
    vel = pick(mom, "velocity_magnitude", "velocity", "vmag")
    first_ok(f"{name} velocity = {expr_name}",
             lambda: setv(vel.value, expr_name),
             lambda: setv(vel, expr_name),
             lambda: setv(vel, {"option": "expression", "expression": expr_name}))
    log(f"{name} velocity:", vel())


def set_inlet_components(s, name, u_name, v_name):
    """Cases B/C: Cartesian velocity components (u, v) given by named expressions."""
    vi = s.setup.boundary_conditions.velocity_inlet[name]
    mom = vi.momentum if "momentum" in getattr(vi, "child_names", []) else vi
    log(f"{name} momentum settings:", getattr(mom, "child_names", "?"))
    meth = pick(mom, "velocity_specification_method", "velocity_spec")
    setv(meth, option(meth, "Components", "components"))
    vc = pick(mom, "velocity_components", "velocity_component", "velocity")
    log(f"{name} velocity components object:", getattr(vc, "child_names", type(vc).__name__))

    def by_index(attr=None):
        for i, ex in enumerate((u_name, v_name)):
            setv(getattr(vc[i], attr) if attr else vc[i], ex)

    first_ok(f"{name} components = ({u_name}, {v_name})",
             lambda: setv(vc, [u_name, v_name]),
             lambda: by_index(),
             lambda: by_index("value"),
             lambda: setv(vc, {"x": u_name, "y": v_name}))
    log(f"{name} velocity components:", vc())


# ---------------------------------------------------------------- reports
def add_volume_report(s, name, rtype, field, zones):
    vol = s.solution.report_definitions.volume
    vol[name] = {}
    r = vol[name]
    setv(r.report_type, option(r.report_type, *rtype))
    setv(r.field, field)
    setv(r.cell_zones, zones)


def add_surface_report(s, name, rtype, surfaces, field=None):
    surf = s.solution.report_definitions.surface
    surf[name] = {}
    r = surf[name]
    setv(r.report_type, option(r.report_type, *rtype))
    if field:
        setv(r.field, field)
    setv(pick(r, "surface_names", "surfaces"), surfaces)


def setup(solver, msh, p):
    """Build the case. Returns the list of report names logged each sample."""
    s = settings(solver)
    first_ok("read mesh",
             lambda: s.file.read_mesh(file_name=msh),
             lambda: s.file.read(file_type="mesh", file_name=msh),
             lambda: solver.tui.file.read_case(msh))
    try:
        s.mesh.check()
    except Exception as e:  # noqa: BLE001
        log("mesh check call failed (continuing):", e)
    try:
        solver.tui.file.confirm_overwrite("no")
    except Exception:  # noqa: BLE001
        pass

    gen = s.setup.general
    first_ok("transient 2nd order", lambda: setv(pick(gen.solver, "time"), "unsteady-2nd-order"))
    first_ok("laminar", lambda: setv(pick(s.setup.models.viscous, "model"), "laminar"))
    sp = s.setup.models.species
    first_ok("species transport",
             lambda: (setv(sp.model.option, "species-transport"),
                      setv(sp.model.material, "mixture-template")),
             lambda: solver.tui.define.models.species.species_transport("yes", "mixture-template"))

    mix = s.setup.materials.mixture["mixture-template"]
    vs = mix.species.volumetric_species
    species = list(vs.get_object_names())
    global TRACER
    TRACER = species[0]
    log("mixture species:", species, "| contaminant tracer:", TRACER)
    setv(mix.density.option, option(mix.density.option, "constant", "volume-weighted-mixing-law"))
    setv(mix.viscosity.option, option(mix.viscosity.option, "constant", "mass-weighted-mixing-law"))
    if mix.density.option() == "constant":
        setv(mix.density.value, RHO)
    if mix.viscosity.option() == "constant":
        setv(mix.viscosity.value, MU)
    for name in species:
        for prop, val in (("density", RHO), ("viscosity", MU)):
            obj = getattr(vs[name], prop)
            if obj.is_active():
                setv(obj.option, option(obj.option, "constant"))
                setv(obj.value, val)
    try:
        setv(s.setup.models.energy.enabled, False)
    except Exception as e:  # noqa: BLE001
        log("energy equation left on:", e)
    md = mix.mass_diffusivity
    first_ok("mass diffusivity",
             lambda: (setv(md.option, "constant-dilute-appx"), setv(pick(md, "value"), DIFF)),
             lambda: (setv(md.option, "constant-dilute-appx"),
                      setv(pick(md, "constant_dilute_appx", "constant"), DIFF)))

    # ---- boundary conditions
    ne = s.setup.named_expressions
    ne["u_in"] = {"definition": p["u_expr"]}
    log("u_in =", ne["u_in"].definition())
    if p["v_mode"] == "none":
        set_inlet_magnitude(s, "inlet", "u_in")                     # Case A: unchanged baseline
    elif p["v_mode"] == "rotate":
        # Case R: the whole far-field stream enters through the top at angle phi; sides are open
        ne["u_top"] = {"definition": p["top_u"]}
        ne["v_top"] = {"definition": p["top_v"]}
        for n_ in ("u_top", "v_top"):
            log(f"{n_} =", ne[n_].definition())
        set_inlet_components(s, "top", "u_top", "v_top")
    else:
        ne["v_amp"] = {"definition": p["v_expr"]}                   # V(t) >= 0, toward the wafer
        ne["v_left"] = {"definition": inlet_v_expr(p["H"])}         # 0 at wafer -> -V at top corner
        ne["v_top"] = {"definition": "-v_amp"}
        for n_ in ("v_amp", "v_left", "v_top"):
            log(f"{n_} =", ne[n_].definition())
        set_inlet_components(s, "inlet", "u_in", "v_left")
        set_inlet_components(s, "top", "u_in", "v_top")

    try:
        s.solution.methods.p_v_coupling.flow_scheme = "PISO"
        log("pressure-velocity coupling: PISO")
    except Exception as e:  # noqa: BLE001
        log("PISO not set, keeping default:", e)

    # ---- reports: cleaning metric, in-trench mixing, verification quantities
    cav = p["cavity_zones"]
    fluid = ["channel"] + cav
    add_volume_report(s, "c_cavity", ("volume-average",), TRACER, cav)
    add_volume_report(s, "rho_avg", ("volume-average",), "density", fluid)
    add_volume_report(s, "c_domain", ("volume-average",), TRACER, fluid)
    add_volume_report(s, "y_max", ("volume-max", "volume-maximum"), TRACER, fluid)
    add_volume_report(s, "y_min", ("volume-min", "volume-minimum"), TRACER, fluid)
    add_volume_report(s, "y_max_cav", ("volume-max", "volume-maximum"), TRACER, cav)       # peak left in trench
    add_volume_report(s, "u_cav", ("volume-average",), "velocity-magnitude", cav)          # in-trench stirring
    names = ["c_cavity", "c_domain", "y_min", "y_max", "y_max_cav", "u_cav"]
    if "cavity-lower" in cav:
        add_volume_report(s, "c_deep", ("volume-average",), TRACER, ["cavity-lower"])     # lower half
        add_volume_report(s, "u_deep", ("volume-average",), "velocity-magnitude", ["cavity-lower"])
        names += ["c_deep", "u_deep"]
    if p["v_mode"] == "rotate":
        boundaries, outflows = ["top", "left", "outlet"], ["left", "outlet"]
    else:
        boundaries = ["inlet"] + ([] if p["v_mode"] == "none" else ["top"]) + ["outlet"]
        outflows = ["outlet"]
    for b in boundaries:
        add_surface_report(s, f"mdot_{b}", ("surface-massflowrate", "mass-flow-rate"), [b])
        names.append(f"mdot_{b}")
    for b in outflows:
        add_surface_report(s, f"y_{b}", ("surface-massavg", "mass-weighted-avg"), [b], TRACER)
        names.append(f"y_{b}")

    try:   # every time step, for the contaminant balance
        s.solution.monitor.report_files["reports_rfile"] = {
            "report_defs": names, "file_name": "reports.out",
            "frequency_of": "time-step", "frequency": 1}
        log("report file: reports.out (every time step)")
    except Exception as e:  # noqa: BLE001
        log("report file not created (history.csv still written):", e)

    setv(s.solution.run_calculation.transient_controls.time_step_size, p["dt"])
    eqs = s.solution.monitor.residual.equations
    n_set = 0
    for name in eqs.get_object_names():
        try:
            setv(eqs[name].absolute_criteria, p["res_tol"])
            n_set += 1
        except Exception as e:  # noqa: BLE001
            log(f"residual criterion {name} not set: {e}")
    if n_set == 0:
        raise RuntimeError("no residual criterion could be set")
    log(f"residual criterion {p['res_tol']} on {n_set} equations")
    return names


def reports(solver, names):
    """{name: value} for several report definitions."""
    s = settings(solver)
    try:
        res = s.solution.report_definitions.compute(report_defs=names)
        merged = {}
        for item in res if isinstance(res, (list, tuple)) else [res]:
            if isinstance(item, dict):
                merged.update(item)
        if all(n in merged for n in names):
            return {n: scalar(merged[n]) for n in names}
    except Exception as e:  # noqa: BLE001
        log("combined report compute failed, computing one by one:", e)
    return {n: scalar(s.solution.report_definitions.compute(report_defs=[n])) for n in names}


def initialize(solver, cavity_zones=("cavity",)):
    s = settings(solver)
    ini = s.solution.initialization
    first_ok("standard initialization", lambda: ini.standard_initialize(), lambda: ini.initialize())
    r = reports(solver, ["rho_avg", "c_cavity"])
    log(f"mixture density check: {r['rho_avg']} kg/m3")
    if abs(r["rho_avg"] - RHO) > 1e-3 * RHO:
        raise RuntimeError(f"mixture density is {r['rho_avg']}, expected {RHO}")
    if abs(r["c_cavity"]) > 1e-12:
        raise RuntimeError(f"cavity concentration after initialization is {r['c_cavity']}, expected 0")
    for var in ("species-0", "species-1", TRACER):
        try:
            ini.patch.calculate_patch(cell_zones=list(cavity_zones), variable=var, value=1.0)
        except Exception as e:  # noqa: BLE001
            log(f"patch {var}: {type(e).__name__}: {e}")
            continue
        c = reports(solver, ["c_cavity"])["c_cavity"]
        log(f"patch {var}: cavity {TRACER} = {c}")
        if abs(c - 1.0) < 1e-6:
            return
        ini.patch.calculate_patch(cell_zones=list(cavity_zones), variable=var, value=0.0)
    raise RuntimeError("could not patch contaminant into the cavity")


def iterate(solver, n, iters):
    rc = settings(solver).solution.run_calculation
    first_ok(f"advance {n} steps",
             lambda: rc.dual_time_iterate(time_step_count=n, max_iter_per_step=iters),
             lambda: (setv(rc.transient_controls.time_step_count, n),
                      setv(rc.transient_controls.max_iter_per_time_step, iters),
                      rc.calculate()),
             lambda: solver.tui.solve.dual_time_iterate(n, iters))


def write(solver, kind, name):
    s = settings(solver)
    tui = {"case": "write_case", "data": "write_data", "case-data": "write_case_data"}[kind]
    first_ok(f"write {name}",
             lambda: s.file.write(file_type=kind, file_name=name),
             lambda: getattr(solver.tui.file, tui)(name))


def read_case_data(solver, name):
    s = settings(solver)
    first_ok(f"read {name}",
             lambda: s.file.read_case_data(file_name=name),
             lambda: s.file.read(file_type="case-data", file_name=name),
             lambda: solver.tui.file.read_case_data(name))


# ---------------------------------------------------------------- checkpoint / restart
EXIT_CONTINUE = 99


def parse_slurm_time(s):
    """SLURM time ('D-HH:MM:SS', 'HH:MM:SS', 'MM:SS') -> seconds; None if not a time."""
    s = (s or "").strip()
    if not s or not s[0].isdigit():
        return None
    days, _, rest = s.partition("-") if "-" in s else ("0", "", s)
    parts = [int(x) for x in rest.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    h, m, sec = parts
    return int(days) * 86400 + h * 3600 + m * 60 + sec


def checkpoint(solver, out, state):
    write(solver, "case-data", "checkpoint_tmp.cas.h5")
    for ext in ("cas.h5", "dat.h5"):
        os.replace(out / f"checkpoint_tmp.{ext}", out / f"checkpoint.{ext}")
    (out / "checkpoint.json").write_text(json.dumps(state, indent=1))
    log(f"checkpoint written at step {state['steps_done']}")


def truncate_history(hist, steps_done):
    with open(hist, newline="") as f:
        rows = list(csv.reader(f))
    keep = [rows[0]] + [r for r in rows[1:] if r and int(r[0]) <= steps_done]
    with open(hist, "w", newline="") as f:
        csv.writer(f).writerows(keep)
