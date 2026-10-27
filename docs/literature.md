# Literature review

Status of each statement:
- **[S]** stated by the source
- **[A]** abstract only (full text not read)
- **[D]** our deduction
- **[H]** hypothesis this project tests

## Supplied papers

**Paper A.** H. Lin, A. A. Busnaina, I. I. Suni, "Cleaning of High Aspect Ratio Submicron Trenches,"
2002 IEEE/SEMI Advanced Semiconductor Manufacturing Conference, pp. 304–308. Full text read.
- [S] 2D, incompressible Newtonian DI water; Eqs. 1–4; u = Uₛ, v = 0 (steady) and the pulsed waveform;
  no-slip walls; free surface at the air–water interface; periodic single cavity (Fig. 1); GENFLO code.
- [S] Pulsed rinse beats steady rinse; optimum St = Wf/Up ≈ 0.133 for cavities ≳ 100 µm; for 1 µm
  cavities higher frequency gives better removal and no optimum St appears (diffusion-dominated, Table 1).
- [S] In the trench, velocity is 1–2 orders below the pulsed velocity (Fig. 5).
- [S] "Rinse flow normal to the wafer surface also shows orders of magnitude higher cleaning efficiency
  than the parallel flow". No data, geometry or BCs are given for that case.
- [S] Surface tension can hinder infiltration of submicron trenches; acoustic streaming improved Cu fill
  of 0.4 µm trenches (ref. 8).

**Paper B.** H. Kim, Y. Lee, E. Lim, "Design, Fabrication, and Performance Test of a Quartz Horn Megasonic
Waveguide for Semiconductor Cleaning," *ECS Trans.* 41(5), 115–120 (2011), doi:10.1149/1.3630834. Full text read.
- [S] 2D axisymmetric ANSYS harmonic FEM (piezoelectric); anti-resonance 1002 kHz predicted, 1003 kHz measured.
- [S] 110 mm quartz horn over a rotating wafer, DI water 1.25 L/min, 1 MHz drive.
- [S] Averaged acoustic pressure 27 → 69 kPa for 1 → 10 W input (Fig. 4); "well-distributed" pressure map at 16 W (Fig. 5).
- [S] Particle removal efficiency (≤ 1 µm particles on a bare wafer) > 90 % for ≥ 1 W; 65 % at 0.5 W (Fig. 7).
- [D] The horn's flat bottom faces the wafer (Fig. 1), so incidence is near-normal (θ ≈ 90°). The paper
  reports **no** angle study, flow field, trench geometry or transverse-flow equation.
- [D] 27–69 kPa ≈ 0.025–0.16 W/cm² for a plane wave. This sets the acoustic amplitude range used here.

## Additional literature

- P. A. Deymier et al., "Theoretical calculation of the acoustic force on a patterned silicon wafer during
  megasonic cleaning," *J. Appl. Phys.* 88, 2423 (2000), doi:10.1063/1.1287224.
  - [A] Green's-function pressure field; force vs frequency and angle *with the wafer surface*.
  - [A] Acoustic shear stress ≪ the features' shear strength.
- P. A. Deymier et al., "Streaming and removal forces due to second-order sound field during megasonic cleaning
  of silicon wafers," *J. Appl. Phys.* 88, 6821 (2000), doi:10.1063/1.1323521.
  - [A] Second-order streaming force at a planar silicon/water interface and for a slab.
  - [A] The normal force is too weak for submicron particles; the parallel force may push or roll them.
  - [A] The force is **very sensitive to the incident angle**.
- P. A. Deymier et al., "Second-order sound field during megasonic cleaning of patterned silicon wafers:
  application to ridges and trenches," *J. Appl. Phys.* 90, 4211 (2001), doi:10.1063/1.1398595.
  - [A] First-order pressure and second-order pressure gradient near ridges and trenches; energy concentrates near sub-wavelength features.
  - [A] "Grazing incident waves appear to provide a more efficient way of cleaning inside a trench."
  - [A] Trench resonances enhance the removal force. The comparison is with van der Waals adhesion of particles, not dissolved transport.
- J. S. Bach, H. Bruus, "Theory of pressure acoustics with viscous boundary layers and streaming in curved
  elastic cavities," *J. Acoust. Soc. Am.* 144, 766 (2018), doi:10.1121/1.5049579 (arXiv:1804.05802).
  - [S] Effective pressure-acoustics boundary condition and a slip-velocity condition for boundary-layer streaming.
  - Basis for Case D.
- P. B. Muller, R. Barnkob, M. J. H. Jensen, H. Bruus, "A numerical study of microparticle acoustophoresis
  driven by acoustic radiation forces and streaming-induced drag forces," *Lab Chip* 12, 4617 (2012),
  doi:10.1039/C2LC40612H.
  - [A] Perturbation method: first-order acoustics, then a time-averaged second-order streaming flow.
- Kapila, Deymier, Raghavan et al., acoustic streaming in megasonic cleaning of EUV photomasks.
  - [A] Schlichting streaming in the acoustic boundary layer; streaming velocities and drag on particles.
- Strusevich et al., *Adv. Manuf.* (2013), acoustic streaming in high-aspect-ratio PCB vias.
  - [A] In their simulations, streaming doesn't form inside micro-vias but improves transport at the via mouth.

## What the literature supports for this project

| Claim | Support |
|---|---|
| Incidence angle changes the second-order acoustic field and forces near patterned wafers | [A] Deymier 2000b, 2001 |
| Grazing incidence is more efficient *inside trenches* (for particle removal forces) | [A] Deymier 2001 |
| A transverse **flow** equation or angle-dependent **dissolved-contaminant** removal | **No source found**. Cases B/C impose the flow; Case D would compute it |
| Normal rinse flow cleans better than parallel flow | [S] Paper A conclusion (no data) |
| Megasonic pressure scale at 1 MHz in a single-wafer tool | [S] Paper B: 27–69 kPa |

## Hypotheses tested here

- [H1] Transverse inflow toward the wafer (V > 0) increases exchange between trench and channel, more for
  AR 5 than AR 1, because it pushes the dividing streamline at the mouth into the trench (Cases B, Bf).
- [H2] Pulsing V out of phase with U (Case C) breaks the closed in-trench vortex more than steady V of
  the same mean, by moving the stagnation region over the mouth each period.
- [H3] (Case D) At Paper B intensities, megasonic streaming generated at a 1 µm trench is too weak to change
  dissolved-contaminant removal; any angle effect would come from bulk streaming or from trench resonances.
