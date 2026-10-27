"""Case parameters, boundary-condition expressions and scale formulas (pure Python, no Fluent).

Coordinates (same as the baseline 3d_nand/ansys_busnaina model):
    x  along the wafer, in the rinse direction;  y  normal to the wafer, pointing up into the liquid.
    The wafer surface is y = 0; the trench occupies 0 <= x <= W, -D <= y <= 0.
    u = x-velocity (the paper's U), v = y-velocity (the paper's V; v < 0 moves toward the wafer).

Transverse flow (Cases B/C, an *imposed* far-field flow, not an acoustic calculation):
    top boundary  y = H  : (u, v) = (U(t), -V(t))
    left inlet    x = -Lu: (u, v) = (U(t), -V(t) * y / H)   (v = 0 at the wafer, -V at the top corner)
    The pressure outlet takes the extra flux; Navier-Stokes decides how the flow turns at the wafer.
    lambda = V / U and the flow-approach angle phi = atan(lambda) describe the IMPOSED FLOW.
    phi is NOT the megasonic incidence angle theta (see beam_direction()).

Acoustic angle convention (Case D, used by the scale helpers here):
    theta is measured from the wafer surface: theta = 0 deg grazing (parallel), 90 deg normal.
    The incident beam travels toward the wafer along  k_hat = (cos theta, -sin theta).
"""
import csv
import math
from pathlib import Path

RHO, MU, DIFF = 998.2, 1.003e-3, 6e-9     # water; contaminant diffusivity (from Paper A, Table 1)
C_SOUND = 1480.0                          # m/s, water
NU = MU / RHO

V_MODES = ("none", "steady", "follow", "pulse")


# ---------------------------------------------------------------- waveforms
def half_wave(s):
    """Paper A pulse shape: 0 for the first half period, sin(s - pi) in the second.
    Equal to 0.5 * (|sin s| - sin s), which is what the Fluent expressions use."""
    return 0.5 * (abs(math.sin(s)) - math.sin(s))


def half_wave_expr(omega, phase_rad=0.0):
    # no phase term when phase = 0, so Case A is character-identical to the baseline expression
    arg = f"({omega:.10g} [s^-1] * t)" if phase_rad == 0 else f"({omega:.10g} [s^-1] * t + {phase_rad:.10g})"
    return f"0.5 * (abs(sin{arg}) - sin{arg})"


# ---------------------------------------------------------------- case table
def load_case(name, path):
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["case"] == name:
                return row
    raise SystemExit(f"case '{name}' not in {path}")


def _num(row, key, default):
    v = row.get(key)
    return float(v) if v not in (None, "") else default


def plan(row, steps_per_period=20, min_steps=2000, cfl=1.0, res_tol=1e-6):
    """Time step, boundary-condition expressions and dimensionless numbers for one case."""
    cfl = _num(row, "cfl", cfl)
    res_tol = _num(row, "res_tol", res_tol)
    iters = int(_num(row, "iters", 0)) or None
    W, AR, n = float(row["W_m"]), float(row["AR"]), int(row["n_per_w"])
    H = 2.5 * W
    u_avg, f, t_end = float(row["u_avg_m_s"]), float(row["freq_Hz"]), float(row["t_end_s"])
    h = W / n
    periods = []

    # rinse velocity U(t): identical to the baseline model
    if row["mode"] == "osc":
        Us, Up = 0.0, math.pi * u_avg
        u_expr = f"{Us} [m/s] + {Up:.6g} [m/s] * {half_wave_expr(2 * math.pi * f)}"
        u_max, St = Us + Up, W * f / Up
        periods.append(1.0 / f)
    else:
        Us, Up, St = u_avg, 0.0, None
        u_expr = f"{u_avg} [m/s]"
        u_max = u_avg

    # transverse velocity magnitude V(t) >= 0 (directed toward the wafer: v = -V)
    v_mode = (row.get("v_mode") or "none").strip()
    if v_mode not in V_MODES:
        raise SystemExit(f"v_mode '{v_mode}' not one of {V_MODES}")
    lam = _num(row, "v_ratio", 0.0)
    f_v = _num(row, "v_freq_Hz", 0.0)
    phase = math.radians(_num(row, "v_phase_deg", 0.0))
    if v_mode == "none":
        v_expr, v_max, v_mean = None, 0.0, 0.0
    elif v_mode == "steady":
        v_expr, v_max, v_mean = f"{lam * u_avg:.6g} [m/s]", lam * u_avg, lam * u_avg
    elif v_mode == "follow":   # constant ratio V/U: the imposed flow keeps a fixed approach angle
        v_expr, v_max, v_mean = f"{lam:.6g} * u_in", lam * u_max, lam * u_avg
    else:                      # pulse: own frequency/phase, same half-sine shape, mean = lam * u_avg
        if f_v <= 0:
            raise SystemExit("v_mode 'pulse' needs v_freq_Hz > 0")
        vp = math.pi * lam * u_avg
        v_expr = f"{vp:.6g} [m/s] * {half_wave_expr(2 * math.pi * f_v, phase)}"
        v_max, v_mean = vp, lam * u_avg
        periods.append(1.0 / f_v)
    if v_mode != "none" and lam < 0:
        raise SystemExit("v_ratio must be >= 0 (V is the speed toward the wafer)")

    dt = min([t_end / min_steps, cfl * h / (u_max + v_max)] + [T / steps_per_period for T in periods])
    n_steps = math.ceil(t_end / dt)
    dt = t_end / n_steps
    return dict(W=W, AR=AR, H=H, n_per_w=n, h=h, u_avg=u_avg, Us=Us, Up=Up, freq=f, St=St,
                t_end=t_end, dt=dt, n_steps=n_steps, u_expr=u_expr,
                v_mode=v_mode, v_ratio=lam, v_freq=f_v, v_phase_deg=math.degrees(phase), v_expr=v_expr,
                v_max=v_max, v_mean=v_mean, phi_deg=math.degrees(math.atan(lam)) if v_mode in ("steady", "follow") else None,
                top_type="symmetry" if v_mode == "none" else "velocity-inlet",
                cfl=cfl, res_tol=res_tol, iters=iters,
                Re=RHO * u_avg * W / MU, Pe_table1=W * 0.1 * u_avg / DIFF)


def inlet_v_expr(H):
    """v on the left inlet: linear from 0 at the wafer to -V at the top corner."""
    return f"-v_amp * y / {H:.10g} [m]"


def inflow_rates(U, V, H, L_top):
    """Volume inflow per unit depth [m^2/s] through left inlet and top for the profiles above.
    The left inlet's v is tangential to it, so it carries U*H; the top carries V*L_top."""
    return U * H, V * L_top


# ---------------------------------------------------------------- acoustic conventions (Case D)
def beam_direction(theta_deg):
    """Unit propagation vector of the incident wave, theta from the wafer surface, travelling toward it."""
    t = math.radians(theta_deg)
    return math.cos(t), -math.sin(t)


def acoustic_scales(p_amp, f, theta_deg=90.0):
    """Order-of-magnitude scales of a plane wave of pressure amplitude p_amp [Pa] at frequency f [Hz]."""
    w = 2 * math.pi * f
    U0 = p_amp / (RHO * C_SOUND)
    kx, ky = beam_direction(theta_deg)
    return dict(
        intensity_W_cm2=p_amp ** 2 / (2 * RHO * C_SOUND) / 1e4,
        wavelength=C_SOUND / f,
        stokes_layer=math.sqrt(2 * NU / w),          # viscous acoustic boundary layer
        U0=U0,                                       # acoustic velocity amplitude
        displacement=U0 / w,
        u_tangential=U0 * abs(kx),                   # incident-wave velocity along the wafer
        u_normal=U0 * abs(ky),                       # incident-wave velocity toward the wafer
        streaming_scale=U0 ** 2 / C_SOUND,           # Rayleigh/Schlichting streaming velocity scale
        dispersion_scale=(U0 / w) ** 2 * w,          # oscillatory-dispersion scale, compare with DIFF
    )


# ---------------------------------------------------------------- metrics
def sherwood(times, R, depth, W, D=DIFF):
    """Instantaneous trench Sherwood number from R(t) = C/C0:
    flux out of the mouth = h * W * C_cavity  ->  h = -depth * d(ln R)/dt,  Sh = h W / D."""
    out = []
    for i in range(1, len(R)):
        if R[i] > 0 and R[i - 1] > 0 and times[i] > times[i - 1]:
            h = -depth * (math.log(R[i]) - math.log(R[i - 1])) / (times[i] - times[i - 1])
            out.append((0.5 * (times[i] + times[i - 1]), h * W / D))
    return out


CASES_CSV = Path(__file__).resolve().parent.parent / "cases.csv"
