# Mathematical formulation

## 1. Coordinates and variables

| Symbol | Meaning |
|---|---|
| x | along the wafer, in the rinse direction |
| y | normal to the wafer, pointing up into the liquid; wafer surface y = 0 |
| u | x-velocity = the paper's **U** (rinse component) |
| v | y-velocity = the paper's **V** (transverse: into/out of the trench; v < 0 is toward the wafer) |
| Y | contaminant mass fraction (dilute passive tracer); C/C0 = Y for the initial condition used |

Domain (unchanged from the baseline): channel −Lu ≤ x ≤ W + Ld, 0 ≤ y ≤ H, with Lu = 2W, Ld = 3W,
H = 2.5W; trench 0 ≤ x ≤ W, −D ≤ y ≤ 0, D = AR·W. This is 2D, so "transverse" here means **y only**;
motion along the trench axis (z) would need a 3D model.

## 2. Governing equations (all cases)

Continuity, incompressible: ∂u/∂x + ∂v/∂y = 0

Momentum, constant-property Newtonian fluid: ρ[∂**u**/∂t + (**u**·∇)**u**] = −∇p + μ∇²**u** (+ **f**_ac, Case D only)

Gravity is omitted, as in Paper A: with constant density it only adds a hydrostatic pressure.

Contaminant: ∂Y/∂t + ∇·(Y**u**) = D∇²Y, with no source (S_C = 0)

ρ = 998.2 kg/m³, μ = 1.003×10⁻³ Pa·s, D = 6×10⁻⁹ m²/s (back-calculated from Paper A, Table 1: Pe = W·0.1Up/D).
The contaminant is **dissolved** (ionic K⁺ in Paper A), not particles or surface-bound species.

## 3. Boundary conditions

| Boundary | Case A (baseline) | Cases B/C (imposed transverse inflow) |
|---|---|---|
| Left inlet x = −Lu | u = U(t), v = 0 (velocity magnitude normal to the boundary), Y = 0 | u = U(t), **v = −V(t)·y/H**, Y = 0 |
| Top y = H | symmetry: v = 0, ∂u/∂y = 0, ∂Y/∂y = 0 (stands in for the paper's free surface) | **velocity inlet: u = U(t), v = −V(t)**, Y = 0 |
| Right outlet x = W + Ld | pressure outlet, p = 0, backflow Y = 0 | same; takes the inflow U·H + V·(Lu+W+Ld) |
| Wafer, trench walls | no slip, u = v = 0; zero contaminant flux | same |
| Initial | u = v = 0; Y = 1 in the trench, 0 elsewhere | same |

U(t), the rinse waveform from Paper A, is unchanged:
- steady: U = u_avg = 0.15 m/s
- pulsed: U = Us + Up·g(2πft), where g(s) = ½(|sin s| − sin s), Us = 0 and Up = π·u_avg

g equals Paper A's sin[2π(t/T − 0.5)] in the second half period and 0 in the first. Up = π·u_avg reproduces the St values in Paper A, Fig. 4.

V(t) ≥ 0 options (`v_mode` in cases.csv):

| v_mode | V(t) | Meaning |
|---|---|---|
| none | – | Case A |
| steady | λ·u_avg | constant transverse inflow (Case B) |
| follow | λ·U(t) | fixed approach angle φ = atan λ of the imposed flow (Case B, `Bf_`) |
| pulse | π·λ·u_avg·g(2πf_v t + φ_v) | pulsed transverse inflow with its own frequency and phase; mean λ·u_avg (Case C) |

**Compatibility.** The left-inlet v is zero at the wafer and equals −V at the top corner, so it agrees with
no-penetration at y = 0 and with the top inflow. There's no discontinuity at either end. Continuity is
satisfied by the pressure-based solver, with the outlet taking the extra flux U·H + V·L. The fluid mass
balance (Σ ṁ over inlet, top and outlet) is reported at every sample and checked in `analyze.py`.

**Why not V at the existing inlet.** A uniform v ≠ 0 at x = −Lu contradicts v = 0 at both ends of that boundary (wall
and symmetry top), and decays within about H of the inlet, before reaching the trench 2W downstream.
It would show "no effect" for artificial reasons.

**Physical meaning and limits.** The top boundary is where Paper A's idealized rinse stream is defined
(2.5 µm above a 1 µm trench). Cases B/C make that stream *approach* the wafer instead of only sweeping
past it, the "normal rinse flow" Paper A's conclusion mentions without data. It's an **imposed-flow
sensitivity study**, idealized in the same way as the baseline. It doesn't predict megasonic streaming.
In a real cleaner, the near-wall flow would be a stagnation/shear boundary layer far thicker than this
domain (§6).

## 4. Angles: two different quantities

| Symbol | Definition | Used in |
|---|---|---|
| φ = atan(V/U) | approach angle of the **imposed far-field flow** | Cases B/C |
| θ | angle between the **incident megasonic wave's propagation direction** and the wafer surface; θ = 0° grazing, 90° normal; k̂ = (cos θ, −sin θ) | Case D |

φ and θ aren't interchangeable. Rotating the acoustic beam doesn't rotate the liquid velocity by the same
angle: the liquid responds to the time-averaged acoustic force, whose direction and size near walls
and inside trenches are set by viscous boundary layers and scattering. Deymier et al. (JAP 2000, 2001)
measure the angle from the wafer surface, the same as θ here, so no conversion is needed.

## 5. Formula-change audit

| Item | Paper A | Baseline code (3d_nand, 141e41b) | This project | Reason | Support | Test |
|---|---|---|---|---|---|---|
| Continuity, NS | Eqs. 1–3, 2D incompressible | Fluent pressure-based, laminar | unchanged | – | Paper A | Case A reproduces fig4_f200k |
| Species | Eq. 4 | passive tracer, D = 6e-9 | unchanged | – | Paper A, Table 1 | same |
| Inlet u | Uₛ / pulsed form | Us = 0, Up = π·u_avg | unchanged (character-identical expression) | – | Paper A, Fig. 4 St values | `test_case_a_plan_identical_to_baseline_fig4_f200k` |
| Inlet v | v = 0 | v = 0 | v = −V·y/H (Cases B/C) | transverse inflow compatible with wall and top | deduction | `test_left_inlet_profile_matches_wall_and_top` |
| Top | free surface | symmetry | velocity inlet (U, −V) (Cases B/C) | introduces V where the far field is defined | deduction; Paper A mentions normal flow | λ = 0 vs Case A (verification §2) |
| Periodic BCs | periodic cell | inlet/outlet | inlet/outlet | removed contaminant must leave | documented in baseline | – |
| Momentum source | none | none | **f**_ac in Case D only | acoustic streaming force | Bach & Bruus 2018; Muller et al. 2012 | analytic flat-wall streaming (planned) |
| Metric | average concentration in trench | R = volume average over trench / initial; η = 1 − R | same + Sherwood number Sh = h·W/D, h = −D_trench·d(ln R)/dt | exchange rate independent of the run length | deduction | `test_sherwood_of_exponential_decay` |
| St, Pe | St = Wf/Up; Pe = W·0.1Up/D | same | same; V/u_avg for the transverse amplitude | – | Paper A | – |

## 6. Acoustic scales (why Case D is a separate model)

These are plane-wave estimates in water at 1 MHz; `trenchflow.physics.acoustic_scales` computes them.

| Quantity | Paper B range (27–69 kPa, 1–10 W) | 1 W/cm² (172 kPa) |
|---|---|---|
| intensity | 0.025–0.16 W/cm² | 1 |
| acoustic velocity U₀ | 1.8–4.7 cm/s | 11.6 cm/s |
| displacement | 3–7 nm | 19 nm |
| viscous acoustic boundary layer √(2ν/ω) | 0.57 µm | 0.57 µm |
| steady streaming scale U₀²/c | 0.2–1.5 µm/s | 9 µm/s |
| oscillatory-dispersion scale ξ²ω | 0.05–0.35 ×10⁻⁹ m²/s | 2×10⁻⁹ m²/s (D = 6×10⁻⁹) |

Deductions:
1. A 1 µm trench lies **inside** the viscous acoustic boundary layer.
2. Steady streaming generated at the trench is 3–5 orders smaller than the rinse-driven velocity in the trench (mm/s–cm/s).
3. The oscillatory motion's dispersion could matter only near the 1 W/cm² end.
4. Bulk (Eckart) streaming builds up over cm-scale paths and reaches the wafer as a jet whose far-field direction follows the beam. Locally that's a near-wall flow like Cases B/C, with a magnitude that depends on the cleaner geometry, which isn't known here.

## 7. Case D design (not yet implemented)

A perturbation model in the framework of Muller et al. (2012) and Bach & Bruus (2018):
1. First-order viscous acoustics in frequency domain: an incident plane wave at angle θ plus the
   scattered field around the trench, with silicon modelled first as rigid.
2. Time-averaged second-order streaming: a Reynolds-stress body force **f**_ac, or a slip velocity on
   flat parts.
3. Couple to Fluent as a momentum source profile, keeping U(t) and the contaminant model unchanged.
4. Compare at **constant incident intensity** I₀, in the Paper B pressure range and at 1 W/cm². The power per unit wafer area then scales as sin θ, and both are reported.

Fluent can't solve frequency-domain viscous acoustics, so steps 1–2 need a separate solver
(a Python finite-element code). It would be verified first against the analytic Rayleigh/Schlichting
streaming over a flat wall.

**Missing before implementation:**
- Full texts of Deymier et al. 2000/2001 (equations, trench resonance results, silicon elastic model).
- The acoustic pressure at the wafer for an oblique beam, since Paper B is near-normal only.
- Agreement on whether trench resonances should be resolved, which needs elastic silicon.
