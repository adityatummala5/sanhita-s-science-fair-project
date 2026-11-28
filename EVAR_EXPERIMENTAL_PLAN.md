# DentiFlow EVAR: Experimental Plan for Finding and Testing the Optimal Texture

Companion to `DentiFlow_Proposal.pdf` (EVAR version). This document turns the proposal into a step-by-step plan:
- the physics and literature that constrain the design,
- how each design is generated,
- how it is simulated in OpenFOAM,
- how the optimum is chosen,
- how the selected designs are fabricated and tested on the bench.

The plan reuses the existing coronary pipeline in `cfd/` (Carreau blood, Poiseuille validation, Lagrangian platelets with Sorensen kinetics) wherever it can.

---

## 0. Changes to the proposal, and why

These come out of the scaling analysis in §1. Each one changes what gets measured or how. Fix them before the plan is frozen.

| # | Issue in the proposal | Why it matters | Change |
|---|---|---|---|
| 1 | Low-WSS area is defined as TAWSS < 0.4 Pa. | The proposal itself says the limb's mean WSS is about 0.2–0.3 Pa. At the worst site, *every* surface, smooth or textured, is 100% "low-WSS". The metric cannot discriminate between designs. | Use a **relative** threshold: the wetted area where TAWSS is below 50% of the smooth-surface TAWSS at the same site. Report 0.4 Pa only at the macro scale. |
| 2 | "Keep OSI within 5%" is a stated constraint. | Feature-scale flow is quasi-steady Stokes flow (§1.2). For flow along one axis, the local WSS is the driving shear times a fixed spatial pattern, so **OSI is set by the waveform, not the texture**. Passing this test is almost automatic. | Keep OSI as a reported check, but do not present it as a finding. The constraints that carry information are TAWSS ratio, relative low-WSS area and RRT. |
| 3 | The DOE factors are dimensional (h in µm, s in µm). | In shear-driven Stokes flow the hydrodynamics depend only on *shape ratios*, not on absolute size (§1.3). Absolute size matters only for platelets (~3 µm), bacteria (~1 µm) and fabrication. | Run the CFD DOE on dimensionless factors: h/s, net height/h, face angle, riblet ratio and orientation. Choose the absolute scale afterwards, from the platelet model and the printer resolution. |
| 4 | Static PRP adhesion in well plates is the main platelet test. | A ~50 µm groove-flushing mechanism only exists **under flow**. A static assay tests contact area and chemistry, not the hypothesis. | Make the **flow-chamber assay primary** for the hydrodynamic-scale designs. Keep the static assay as a secondary test and for the cell-scale (2 µm) designs. |
| 5 | Sharklet-scale (2 µm) effects are planned for the CFD. | The Lagrangian platelets are point particles. Features smaller than a platelet cannot be resolved by this model. | Treat cell-scale effects as a **bench-only** question. The CFD covers only the ≥ 10 µm hydrodynamic texture. |
| 6 | Undercut crowns are a design level. | Undercuts cannot be demolded from a rigid master and tear soft molds. They also produced the pilot's worst result. | Keep undercuts in the CFD screen to explain the pilot result. Exclude them from fabrication. |
| 7 | The CCD includes crown shape. | Heightfield vs undercut is categorical. A CCD needs continuous factors, and its axial points (α > 1) fall outside physical bounds such as h/s < 0.1. | Use a **face-centred CCD** (α = 1) on 3 continuous factors, run at each level of the categorical factor. |
| 8 | Bench stats: n ≥ 3 coupons × 3 days, one-way ANOVA. | Pooling days as independent replicates inflates n. Day-to-day donor variation is usually the largest variance component. | Use a **mixed model with day as a random block**, with Dunnett contrasts on log-transformed data. Use about 4 coupons per surface per day (§8.5). |
| 9 | The pilot's +91% platelet result is treated as a finding. | `cfd/results/platelet_comparison/` shows it is **one random seed**, dominated by aggregation (+551% aggregated; 2 clusters drive most of it). | Replicate it with ≥ 5 seeds before building a mechanism on it. Report an aggregation-excluded metric alongside it. |
| 10 | "BSL-1 surrogate per fair rules." | The common BSL-1 *S. epidermidis* strain, ATCC 12228, lacks the *ica* operon and **does not form biofilm** [Zhang 2003]. It would give a null result by design. | At the institutional lab, use a documented biofilm former (§8.4). Confirm biofilm on polystyrene before any test coupons are used. |
| 11 | PLGA coupons are incubated at 37 °C for 24 h. | PLGA's dry Tg is about 45–55 °C, and water plasticises it. Textures may relax during incubation. Check this; do not assume it. | Re-measure textures *after* a mock incubation (§7.5). Have a fallback polymer ready, such as polyurethane or PDMS. |

---

## 1. Physics and literature that constrain the design

### 1.1 Scale numbers (12 mm limb, 0.5 L/min, Carreau blood, period 0.9 s)

| Quantity | Value | How computed |
|---|---|---|
| Poiseuille mean WSS | ≈ 0.17 Pa (μ = 3.5 mPa·s) to ≈ 0.25 Pa | 4μQ/πR³ |
| Shear rate at τ = 0.25 Pa | ≈ 45 s⁻¹, where Carreau μ ≈ 5.6 mPa·s | Cho & Kensey parameters in `cfd/config/experiment.json` |
| Womersley number (limb) | ≈ 8.7 | R·√(ω/ν) |
| Mean Reynolds number (limb) | ≈ 270 | |
| Stokes-layer thickness | ≈ 1.0 mm | √(2ν/ω) |
| Feature Reynolds number Re_k = γ̇k²/ν | 0.001 (10 µm), 0.02 (50 µm), 0.09 (100 µm) | |
| Viscous diffusion time across a feature, k²/ν | 0.02 ms (10 µm) to 2 ms (100 µm) | Compare with the 900 ms cardiac period |

### 1.2 Consequence 1: feature-scale flow is quasi-steady Stokes flow

The diffusion time across a feature is at least 450× shorter than the cardiac period, and Re_k ≪ 1. So at every instant the flow inside the texture is the steady creeping-flow solution for that instant's wall shear. For a Newtonian fluid this makes it **linear**:

> **u**(x, t) = S_x(t)·**U**_x(x) + S_z(t)·**U**_z(x)

where S_x and S_z are the streamwise and spanwise driving shear from the macro model.

What follows from this:
- **Compute cost drops by about 100×.** Two steady Stokes solves per design (unit shear along the flow and across it), plus a few amplitudes to capture Carreau shear-thinning, replace a full pulsatile run. Verify this once against a full `pimpleFoam` pulsatile patch run (§4.4).
- **OSI is essentially texture-independent** under flow along one axis. When the driving flow reverses, the whole field reverses. This is why change #2 above is needed.
- **A reversing flow exposes both faces of an asymmetric crown.** Whenever OSI > 0, the "downstream" face takes a turn as the upstream face. Both face angles therefore matter (§1.4).
- **Drag reduction is off the table.** The riblet mechanism needs turbulence: optimum spacing s⁺ ≈ 15, or groove-size ℓg⁺ ≈ 11 [Bechert 1997; García-Mayoral & Jiménez 2011]. In Stokes flow, texture adds wetted area and the only lever is the *virtual origin* of the outer flow [Luchini, Manzo & Pozzi 1991].

### 1.3 Consequence 2: the hydrodynamics are scale-invariant

The patch is driven by a *shear*, not a velocity. Local shear rate is therefore γ̇_drive × F(x/s), whatever the absolute size, and the Carreau viscosity field scales the same way.

All normalised hemodynamic metrics depend only on dimensionless shape (h/s, angles, riblet ratios), provided two conditions hold:
- Re_k stays ≪ 1,
- the top boundary is far away, H ≫ s.

**Test this explicitly** as a sanity check: run one shape at s = 100 µm and s = 300 µm. The normalised TAWSS should match within the GCI.

Absolute scale enters only through:
1. platelet and bacterial size relative to the features,
2. fabrication resolution,
3. the patch-to-lumen ratio (< 1%, so negligible).

### 1.4 Literature design rules (testable predictions to pre-register)

| ID | Prediction | Source |
|---|---|---|
| **P1** | **Corner-angle rule.** In Stokes flow, a corner with interior fluid angle below **≈ 146°** contains Moffatt eddies, with WSS → 0 near the corner. For a face rising at angle θ from the floor, the base corner has interior angle 180° − θ. Eddies therefore appear once **θ > ≈ 34°**. | Moffatt 1964 |
| **P2** | **Cavity rule.** In transverse rectangular grooves, flow over a wide cavity (width/depth = 4, i.e. h/s = 0.25) has only corner eddies. At width/depth = 2 (h/s = 0.5) one central vortex fills the cavity. A stagnation step-change is expected between h/s ≈ 0.25 and 0.5. | Shen & Floryan 1985; Higdon 1985 |
| **P3** | **Orientation rule.** Grooves aligned with the flow have *no* recirculation under streamwise shear; fluid in the groove still moves, just more slowly. Transverse grooves (crimp-like) are cavities. Expected ranking for stagnant volume: aligned < 45° < transverse. This changes only if the macro near-wall shear direction rotates strongly during the cycle. | Luchini 1991 (h∥ vs h⊥) |
| **P4** | **Recession rule.** Recessing the texture moves the virtual origin below the nominal surface, so the outer flow sees the wall as slightly further away. This cannot raise local TAWSS inside the grooves. The best achievable result is TAWSS on the crests ≈ flat, with a groove deficit that shrinks as h/s → 0. | Luchini 1991 |
| **P5** | **Cell-scale rule.** Features below platelet and bacterial size reduce adhesion by cutting contact area. | Platelets on PLGA and polyurethane: Koh 2010; Milner 2006. Sharklet on TPU: about 80% less platelet coverage [May 2015]. Microtopographies on PDMS: 30–45% less bacterial adhesion [Perera-Costa 2014]. |

P1 gives a candidate mechanism for the pilot's +91% result: the steep upstream faces (> 34°) and undercuts created corner zones with near-zero WSS. One caveat: closed Moffatt streamlines exchange fluid only by diffusion. They may *exclude* platelets as much as trap them. The Lagrangian model is the way to find out.

---

## 2. Design space and design families

### 2.1 Factors (all dimensionless for the CFD)

| Factor | Symbol | Levels in the OFAT screen | Notes |
|---|---|---|---|
| Groove aspect ratio | h/s | 0.1, 0.2, 0.3, 0.5, 0.75, 1.0 | Brackets P2 |
| Net height | η = net height / h | +1 (protruding), +0.5, 0 (flush), −0.5, −1 (fully recessed) | Replaces ±25/±50 µm |
| Upstream-face angle | θ_up | 20, 30, 35, 40, 55, 70° | Dense around 34° (P1) |
| Downstream-face angle | θ_dn | = θ_up (symmetric), or 20° fixed | Matters because OSI > 0 |
| Riblets per crown | n_r | 0, 3, 5 | |
| Riblet depth ratio | d_r / h | 0, 0.25, 0.5 | |
| Orientation | φ | 0° (along flow), 45°, 90° (transverse) | Categorical in the DOE |
| Crown form | — | heightfield, undercut | Undercut is CFD-only (change #6) |
| Arrangement | — | continuous riblets, segmented diamond denticles | |

**Fixed:** a planform-coverage fraction for denticles (e.g. 50%), and the corner fillet radius (e.g. 0.05 s). A small fillet keeps meshes valid; in fabricated parts it reflects the real print or mold radius.

### 2.2 Named design families

| Code | Family | Description |
|---|---|---|
| **C0** | Smooth | Flat wall, the ePTFE-like control. Every patch metric is normalised to C0 at the same site. |
| **C1** | Polyester crimp | Transverse sinusoidal crimp. Take the crimp pitch and amplitude from a measured graft sample, or published SEM, at real scale. Expected worst case for stagnation. |
| **R** | Continuous riblets | Aligned V-, U- or trapezoid grooves (h/s, θ). The cleanest test of P1, P3 and P4. |
| **D** | Segmented denticles | Mako-type crowns with 3–5 riblets [Wen 2014; Dean & Bhushan 2010], in a diamond layout (as in `shark-denticle-diamond-layout.stl`). |
| **E** | Engraved (recessed) denticles | D with η = −1: the crown top is flush and the outline is a moat. This is the hypothesis design. |
| **S** | Sharklet-type | 2 µm ribs in a diamond array [Chung 2007]. Bench-only (change #5). |
| **H** | Hierarchical | The best E or R design with S-scale ribs on its crests. Bench-only. |

### 2.3 How the designs are generated (code)

Add a periodic-patch generator, `cfd/pipeline/patch_geometry.py`, next to the existing `geometry.py`. It can reuse `manifold3d`/`trimesh` and the `DenticleParams` idea.

1. **Heightfield designs** (R, D, E, C1). Define the surface as y = f(x, z) on a 2 × 1 mm doubly periodic tile. The pattern must be an exact integer number of periods in each direction, so the cyclic faces match.
   - Write f analytically (riblets, crimp) or by rasterising the crown outline with smooth-step edges.
   - Export: (a) an STL for snappy, (b) the raw f(x, z) array for terrain-following meshing (§4.2), and (c) a JSON of measured stats: wetted/planform area ratio, groove volume, h, s and angles.
2. **Undercut designs.** Build them with boolean unions of crown solids on a base slab, as `geometry.py` already does, then trim to the periodic tile. Remove any sliver fragment smaller than 1% of a crown's volume; this is the same clean-up the coronary sector needed.
3. **Validation geometries.** Add a 2D-extruded rectangular cavity (W/h = 4, 2 and 1) and a triangular riblet. These are the literature benchmarks in §5.
4. **Register everything in a config.** Put each design in `cfd/config/patch_experiment.json`, in the same style as `experiment.json`, so every run regenerates from the config.

---

## 3. Stage A: Macro-scale endograft CFD (months 1–2)

**Goal:** find where the limb wall shear is worst. Then extract the near-wall shear time series S(t) that will drive the patch models.

### 3.1 Geometry
- Idealised bifurcated endograft, built parametrically in Python (manifold3d sweep and union) or in CadQuery.
- Body: 24 mm diameter, about 60 mm long. Limbs: 12 mm diameter, about 90 mm long.
- Total limb divergence angle: 40°. Run a second geometry at 60° as a sensitivity case.
- Smooth crotch fillet with radius ≥ 2 mm. A sharp flow divider breaks the mesh and gives a spurious stagnation point.
- Inlet extension: 2D. Outlet extensions: 10D.
- Rigid walls (a stated limitation).
- *Optional:* add a 15–20 mm crimp-scale corrugation on the limb as a "real graft" case.

### 3.2 Boundary conditions and physics
- **Inlet:** time-varying Womersley profile built from the infrarenal flow waveform of Les et al. 2010 (rest and exercise; digitise the published curves). The mean resting infrarenal flow in AAA patients is about 17.5 mL/s ≈ 1.05 L/min [Les, CVET 2010], which gives the proposal's ≈ 0.5 L/min per limb.
- **Outlets:** zero static pressure with `pressureInletOutletVelocity`. The geometry is symmetric, so this gives a 50/50 split. Diastolic flow reversal exists in the infrarenal waveform, which is why the long outlet extensions are needed: they keep backflow away from the region of interest. An asymmetric variant would need resistance BCs (`codedFixedValue` with p = R·Q).
- **Blood:** Carreau, ρ = 1060 kg/m³, the existing `physics` block.
- **Solver:** `pimpleFoam`, adaptive Δt at Co ≤ 1, as in the existing `write_control_pulsatile`.
- **Cycles:** 4. Periodic convergence is reached when area-mean TAWSS changes < 1% between cycles 3 and 4. Analyse only the last cycle.

### 3.3 Mesh
- `snappyHexMesh` core cells of 0.4–0.6 mm.
- **Prism layers to resolve the ≈ 1 mm Stokes layer:** 10 layers, first cell 20–30 µm, growth ≤ 1.2.
- Expect about 2–5 M cells.
- Run a GCI study on 3 meshes [Celik 2008] with a refinement ratio ≥ 1.3. Metrics: limb-averaged TAWSS, OSI at the selected site, and pressure drop. Accept GCI < 5%.

### 3.4 Validation (before any endograft run)
1. **Poiseuille.** Already done in the pilot (−2.5% WSS error).
2. **Womersley.** A straight 12 mm Newtonian tube driven with the first 5–8 harmonics of the infrarenal waveform. Compare velocity profiles and WSS at 8 phases against the analytical solution (Bessel J₀ of complex argument, `scipy.special.jv`). Accept a peak-WSS error < 5% and an RMS profile error < 3%. Add the analytical solution as `cfd/pipeline/womersley.py`; it also generates the inlet profile.

### 3.5 Outputs
- Maps of TAWSS, OSI, RRT and (for reference) the 0.4 Pa low-WSS area for rest and exercise.
- **Site selection:**
  - **S1:** the limb location with the lowest TAWSS × (1 − 2·OSI), excluding the 1D near the outlets.
  - **S2:** a "typical" mid-limb location.
  - **S3:** the S1 location under exercise.
- **Extraction:** at each site, extract the wall-parallel velocity vector at a normal distance H = 0.3 mm (and 0.5 mm), with both components, over the last cycle at ≥ 200 samples per cycle. Also extract the WSS vector. These are the patch drivers. Check that H ≥ 5× the largest feature height and H < the Stokes-layer thickness.
- **Check the near-wall shear direction.** If it rotates more than ±20° during the cycle at S1, orientation becomes a real factor. If not, P3 predicts that aligned grooves win trivially.

---

## 4. Stage B: Wall-patch CFD (months 2–3)

### 4.1 Domain and BCs
- **Tile:** 2 mm streamwise (x) × 1 mm spanwise (z), extending from the deepest groove up to H above the nominal surface. H = 0.3 mm (0.5 mm in the sensitivity check).
- **Sides:** `cyclic` in x and z. Snappy breaks conformal matching, so convert with `createPatch` to `cyclicAMI` when needed.
- **Bottom:** no-slip on the textured wall.
- **Top, quasi-steady mode:** `fixedValue` U = (1, 0, 0)·u_ref (solve A) and (0, 0, 1)·u_ref (solve B), with `simpleFoam`. Because of Carreau shear-thinning, repeat each solve at 5 driving shear levels spanning the waveform: 0.02, 0.1, 0.3, 1 and 3 Pa. Build a lookup table, then reconstruct WSS(x, t) and the velocity snapshots by interpolation in |S(t)|, superposing A and B.
- **Top, full mode (verification only):** `uniformFixedValue` with a table of the extracted u(H, t) vector, solved with `pimpleFoam` for 3 cycles.
- **Why a moving top and not a pressure gradient:** inside H ≪ δ, the macro profile is close to linear. A moving lid reproduces the local shear directly. The error scales roughly with (H/δ)², about 9% at H = 0.3 mm. That is a *driving-amplitude* error, so it cancels in metrics normalised to C0. Confirm with the H = 0.5 mm case.

### 4.2 Mesh
- **Heightfield designs:** a terrain-following structured mesh. Run `blockMesh` on a flat box, then displace the bottom points to y = f(x, z) with vertical blending, using a small Python utility that writes `constant/polyMesh/points`. This gives exact periodicity and orthogonal near-wall cells.
- **Undercut designs:** `snappyHexMesh` with 2 refinement levels at the texture and 3 surface layers.
- **Resolution:** at least 15–20 cells across the narrowest groove or riblet, and at least 8 cells over riblet depth.
- **GCI** [Celik 2008] on 3 meshes for C0, one deep design (h/s = 1) and one sharp design (θ = 70°). Use the TAWSS ratio, relative low-WSS area and groove washout time; accept GCI < 5%.
  - Threshold metrics such as the low-WSS area converge noisily. Celik's least-squares GCI variant handles that.
  - If the area fraction will not converge, report the *distribution* (CDF) of the TAWSS ratio instead.

### 4.3 Metrics (define and freeze them before running)

All are normalised to C0 at the same site and driving shear.

| Metric | Definition |
|---|---|
| **Φ_TAWSS** | Planform-area-weighted TAWSS ÷ C0 TAWSS. Also report it wetted-area-weighted. |
| **A_low,rel** | Wetted area with TAWSS < 0.5 × C0 TAWSS, ÷ planform area |
| **RRT ratio** | Area-mean RRT ÷ C0 RRT (expected ≈ 1/Φ) |
| **OSI** | Area-mean OSI. Reported as a check only (§1.2). |
| **Washout time t₉₀** | Set a passive scalar c = 1 in the fluid below the crest plane at t = 0. Transport it with `scalarTransport` (D = 1 × 10⁻¹¹ m²/s, shear-enhanced platelet scale) under the reconstructed waveform. t₉₀ is the time until 90% has left the groove. This is the direct "groove flushing" metric. |
| **Stagnant volume** | Groove volume where time-averaged \|u\| < 1% of the C0 velocity at crest height |
| **h∥, h⊥** | Luchini protrusion heights from solves A and B: the virtual origin of each linear profile. Also used for validation. |
| **Platelet deposition** | Deposited platelets per **planform** mm² after 2 cycles. Mean ± SD over seeds. Also report it excluding the 2 largest aggregates. |

### 4.4 Quasi-steady verification (do this once, early)

Compare the full pulsatile run against the quasi-steady reconstruction for C0, one E design and one deep transverse design. Accept if Φ_TAWSS, OSI and t₉₀ agree within 2%.

If they don't, unsteady effects matter at feature scale. That would itself be a reportable surprise. In that case, fall back to full runs and a smaller DOE.

### 4.5 Platelet model on the patch
- **Absolute scale.** Fix s for this step, since platelets break scale invariance (§1.3). Run the top 10 shapes at s = 50, 100 and 200 µm.
- **Injection.** Inject platelets at x = 0, uniformly in 0 < y < H, and remove them at x = L. Do not recycle them through the cyclic boundary: that would artificially pre-load the near-wall layer.
- **Kinetics.** Use the existing Sorensen-kinetics code (`cfd/pipeline/platelets.py`). The flow snapshots come from the quasi-steady reconstruction, written in the `FlowField` snapshot format.
- **Seeds.** 3 per design in the screen, 5 per design for the final 6 candidates and anything near the optimum.
- **Parameter-robustness check.** Sorensen's rate constants were fitted on collagen, not on a polymer graft. For the final candidates, rerun with the surface-adhesion rate ×0.1 and ×10. The **ranking** of designs is what has to hold up, not the absolute counts. If the ranking flips, say so: it means the model cannot pick a winner.

---

## 5. Validation ladder (so the judges trust the numbers)

| Level | Case | Reference | Acceptance |
|---|---|---|---|
| V1 | Poiseuille pipe | Analytical | Done: −2.5% |
| V2 | Womersley pipe | Analytical (Womersley 1955) | Peak WSS < 5% |
| V3 | Shear flow over a rectangular cavity, W/h = 4, 2, 1 | Shen & Floryan 1985 (corner eddies at W/h = 4, single vortex at W/h = 2); Higdon 1985 | Correct eddy topology; separation points within 5% of h |
| V4 | Protrusion heights h∥, h⊥ of a triangular riblet | Luchini, Manzo & Pozzi 1991 (tabulated values) | Within 5% |
| V5 | Scale invariance: same shape at s = 100 and 300 µm | §1.3 | Normalised metrics match within the GCI |
| V6 | Quasi-steady vs full pulsatile | §4.4 | Within 2% |
| V7 *(optional, cheap)* | Physical model: a 20–50× scaled transverse groove in a glycerol channel with tracer particles filmed by phone. Keep Re_k < 0.1, so the flow is dynamically similar (§1.3). | Patch CFD of the same geometry | Eddy size and topology match |
| V8 | Flow-chamber CFD of the actual coupon in the actual chamber | Bench platelet distribution (§8.3) | Qualitative: deposits sit where the CFD predicts low shear and long residence |

V3 and V4 are the important ones. They prove the patch solver gets creeping flow over textures right, which is exactly the regime that matters.

---

## 6. Optimisation strategy (months 2–3)

The quasi-steady method makes each hemodynamic evaluation cost minutes. The platelet runs (about 1 CPU-hour each, from the pilot) are the bottleneck, so the budget below is counted in platelet runs.

**Step 1: OFAT screen (about 35 designs, hemodynamics only).**
- Sweep each factor in §2.1 around a baseline (aligned, η = 0, h/s = 0.3, θ = 30°, n_r = 3, d_r/h = 0.25, diamond D layout).
- Add C1 and the 2 pilot-like undercut shapes.
- Purpose: confirm or refute P1–P4, and drop factors whose effect is smaller than 2× the GCI.

**Step 2: Response surface (CFD + platelets).**
- Factors: the 3 continuous factors expected to survive (η, h/s, θ_up).
- Design: face-centred CCD, 8 factorial + 6 face + 1 centre = **15 points**.
- Run at the 2 most promising orientations from Step 1, giving **30 designs**. Each has hemodynamics plus 3 platelet seeds: 90 platelet runs, about 90 CPU-hours.
- Deterministic CFD needs no centre replicates. The platelet seeds give the noise estimate.
- If Step 1 leaves more than 4 important factors, use a **definitive screening design** first [Jones & Nachtsheim 2011]: 13 runs for 6 factors, and it detects curvature.
- Software: Stat-Ease Design-Expert (as planned), or Python `pyDOE3` + `statsmodels`.

**Step 3: Choose the optimum.**
- Fit quadratic models for Φ_TAWSS, A_low,rel, t₉₀ and platelet deposition.
- Plot the **Pareto front** of Φ_TAWSS against platelet deposition.
- Pick the operating point with a **Derringer–Suich desirability** [1980]:
  - **Hard constraints:** Φ_TAWSS ≥ 0.95 (the proposal's 5%) and A_low,rel ≤ C0 + 5 percentage points.
  - **Manufacturability constraint:** no undercut, minimum feature ≥ 3 printer pixels at the chosen scale (§7), and face angles ≤ the mold release limit.
  - **Objective:** minimise platelet deposition.

**Step 4: Confirmation.**
- Rerun the predicted optimum on the fine mesh with 5 seeds.
- It must land within the model's 95% prediction interval. If it doesn't, the response surface is wrong near the optimum: add points there and refit.

**Step 5: Robustness.**
- Evaluate the top 5 at S1-rest, S1-exercise and S2.
- Apply the platelet ×0.1/×10 sensitivity (§4.5).
- The design that ranks best across all conditions is "optimal". A design that wins only at one site is not.

**Step 6: Down-select for the bench (§8.1).**

**Compute estimate:** about 60 hemodynamic designs × 10 quasi-steady solves × 5–10 min, plus about 150 platelet runs × 1 h. That totals roughly 250 CPU-hours, about 3 weeks on 4 cores at night. It fits in month 2–3.

---

## 7. Fabrication (months 3–4)

### 7.1 Route by feature scale

| Feature scale | Route | Notes |
|---|---|---|
| ≥ 100 µm (and C1 crimp) | Consumer resin LCD/DLP printer, about 20–35 µm XY pixel and 10–25 µm layers | Minimum reliable feature is about 3–4 pixels. Sloped faces show layer stair-steps, so print with the slope along the build direction or tilt the part. |
| 20–100 µm | Institutional high-resolution micro-stereolithography (2–10 µm class) or grayscale lithography | Needed if the optimum lands at s < 100 µm |
| 2–5 µm (S, H) | SU-8 photolithography master in a university cleanroom, or two-photon lithography for small areas | Chung 2007 and May 2015 used lithographic masters. If none is available, a replica of a commercial **diffraction grating** (≈ 1 µm ribs) can stand in as a "cell-scale ribs" surface. It is *not* Sharklet, and must not be called that. |

**Use scale invariance (§1.3) to your advantage.** If the printer can't make the optimum at s = 50 µm, fabricate it at the smallest printable s. The hemodynamic part is unchanged. Only the platelet-to-feature size ratio changes, and it can be re-simulated at the as-built scale.

### 7.2 Replication chain (watch the polarity)

1. **Print the master as the final graft surface** (positive), so that grooves are grooves.
2. **Post-treat the resin before casting PDMS.** Most SLA resins inhibit PDMS curing: phosphine-oxide photoinitiators leach out and poison the Pt catalyst [Venzac 2021]. Follow their UV + heat + wash protocol, or coat the master (e.g. silanise or parylene). Test-cast a scrap piece first.
3. **Cast a PDMS (Sylgard 184) negative**, then silanise it for release.
4. **Cast PLGA into the PDMS negative** (solvent casting, the method Koh 2010 used for textured PLGA films).
   - Use a 10–15% w/v solution in ethyl acetate. DCM or chloroform also work, but only in a fume hood.
   - Dry under vacuum for ≥ 48 h to remove residual solvent, which is cytotoxic.
   - Alternative: hot-emboss PLGA above Tg against a rigid mold.
5. **Make every control through the same chain.** C0 is cast against flat PDMS, and C1 against a printed crimp master. Then the material, solvent history and roughness differ only in the designed texture.

### 7.3 Polymer choice
- **PLGA** is fine for a 24 h bench test (§0 #11).
- **Fallback:** medical-grade polyurethane, solution-cast. This is more clinically relevant and is what Milner 2006 and May 2015 used. Or PDMS itself.
- Decide after the mock-incubation test in §7.5.

### 7.4 Metrology and acceptance (every batch)
- **Profile:** optical profilometry, confocal microscopy or stylus profilometry on ≥ 3 coupons per batch. Measure h, s, θ and the corner radius.
  - **Accept** if h and s are within ±10% of the CAD values and θ is within ±5°.
  - **Re-run the CFD on the measured "as-built" profile.** The bench is then compared against the geometry that was actually made.
- **Imaging:** SEM on 1 coupon per design to check for defects and for transfer of the cell-scale features.
- **Wettability:** water contact angle, ≥ 3 drops per coupon. Texture changes wettability (Wenzel/Cassie), which alters protein adsorption. Record it as a covariate.

### 7.5 Stability check
- Incubate blank coupons in PBS at 37 °C for 24 h and 48 h, then re-profile them.
- Also test the sterilisation step: 70% ethanol + UV. PLGA cannot be autoclaved.
- **Accept** if h changes < 10%. Otherwise switch polymer.

---

## 8. Biological testing (months 4–6, supervised institutional lab, SRC/IRB approved)

### 8.1 Bench groups (6, ranked by priority)

| Group | Surface | Purpose |
|---|---|---|
| 1 | C0 smooth | Primary control |
| 2 | C1 polyester crimp | Clinical comparator: "the texture to beat" |
| 3 | E-opt: CFD-optimal recessed design | Hypothesis test |
| 4 | P-opt: the same crown but protruding (η = +1) | Tests recessed vs protruding directly, the pilot's main lesson |
| 5 | S: cell-scale ribs | Anti-adhesion test at the cell scale |
| 6 | H: E-opt + S | Hierarchical test |
| 7 *(only if capacity allows)* | The CFD's worst design | Tests whether the CFD ranking predicts the bench ranking (Spearman ρ across groups 1–4 and 7) |

More groups mean less power per comparison (§8.5). Drop 6 or 7 before dropping replicates.

### 8.2 Static platelet adhesion (secondary for groups 3–4; primary for groups 5–6)
- **Blood source:** citrated PRP from a source the SRC approves. Commercially sourced, de-identified blood products usually avoid human-subjects review. Volunteer draws (including the student's own blood) need IRB approval. Confirm under your fair's rules.
- **Incubation:** 60 min at 37 °C, then rinse 3× in PBS.
- **Quantification:** use both of these.
  - **LDH assay:** lyse with Triton X-100 and read against a standard curve of known platelet numbers. This gives total adhesion.
  - **Fluorescence counts:** calcein-AM or anti-CD41 staining, imaged in 10 pre-defined fields per coupon. Count in Fiji, with the counter blinded to the coupon code. Score *where* platelets sit (crest vs groove vs moat), which is the mechanistic readout.
- **SEM on a subset:** spreading or activation morphology (round, dendritic, spread).

### 8.3 Flow-chamber platelet adhesion (primary for groups 1–4)
- **Chamber:** a parallel-plate chamber with the coupon as the floor. An ibidi sticky-Slide I Luer, which attaches to a custom substrate, is one option; a 3D-printed gasket chamber is another.
- **Channel height:** choose it so the texture is < 10% of the channel height, for example the 0.8 mm channel for ≤ 50 µm features. Otherwise the texture changes the bulk channel flow.
- **Shear:** set it to the CFD-predicted near-wall shear at S1, 0.2–1 Pa, steady first. Optionally add the S1 waveform with a programmable syringe pump.
  - ibidi gives τ = η × 31.0 × Φ for the 0.8 mm Luer slide (τ in dyn/cm², η in poise, Φ in mL/min). The sticky tape and the coupon change the effective height, so **compute the actual τ with a CFD of your chamber + coupon** (this is V8).
  - Rough flow needed for 0.3 Pa: about 8 mL/min for PRP (η ≈ 1.2 mPa·s) and about 2.4 mL/min for whole blood.
- **Fluid:** reconstituted blood (washed RBC at about 40% haematocrit plus PRP with labelled platelets) or labelled whole blood. Red cells push platelets to the wall (margination). PRP alone under-delivers platelets to the wall and will under-estimate the texture effect.
- **Readout:** time-lapse fluorescence of platelet accumulation every 30 s for 10 min. Then fix and run LDH or counts. Overlay the deposit maps on the CFD shear and washout maps (V8).

### 8.4 Bacterial biofilm (groups 1, 2, 3, 5, 6)
- **Strain** (BSL-2, at the institution): *S. aureus* (e.g. ATCC 25923, or a strong biofilm former your lab already uses) or *S. epidermidis* RP62A (ATCC 35984). Staphylococci are the leading endograft pathogens.
  - **Do not use** *S. epidermidis* ATCC 12228: it is *ica*-negative and does not form biofilm [Zhang 2003].
  - If BSL-1 is unavoidable, pick a strain with *documented* biofilm formation and say clearly that it is a surrogate.
- **Biofilm growth:** 24 h, static, in TSB + 1% glucose at 37 °C, on sterile coupons (70% ethanol + UV, then sterility controls).
- **Readouts:**
  1. **Crystal violet** [O'Toole 2011]: 0.1% CV for 15 min, rinse, solubilise in 30% acetic acid, read OD590. Subtract the OD of *uninoculated coupons of the same texture*, because grooves can trap dye.
  2. **CFU per planform cm²:** sonicate and vortex, serially dilute, plate.
  3. **Live/dead imaging** (SYTO9/PI): shows where cells colonise, by analogy with Perera-Costa 2014's "upper vs lower areas".
- **Optional:** biofilm under flow in the same chamber at 0.3 Pa for 24 h.

### 8.5 Statistics and power
- **Design:** 3 independent days (separate blood donors or lots, separate bacterial cultures). On each day, 4 coupons per surface, with positions randomised across wells or the chamber sequence. Image analysis blinded by coupon codes.
- **Model:** a linear mixed model on log(response), with surface as a fixed effect and day as a random effect. Run Dunnett contrasts against C0 (α = 0.05, family-wise). Report the E-opt vs P-opt and E-opt vs C1 contrasts as pre-planned comparisons. Report effect sizes with 95% CIs, not just p-values.
- **Power:**
  - With 12 coupons per group and 5 Dunnett comparisons, there is about 85% power to detect a 30% change if the within-day CV ≤ 20%. That falls to about 65% if the CV is 25%.
  - **Run a 1-day pilot on C0 only first**, measure the CV, and adjust the number of coupons.
- **Pre-registration:** write the hypotheses, primary endpoint, metrics and decision rules (§10) into the lab notebook *before* the first bench day.

---

## 9. Timeline (6 months)

| Weeks | Work | Deliverable / gate |
|---|---|---|
| 1–2 | Literature reading (refs below). Write `womersley.py`. V2 validation. | V2 passes |
| 3–5 | Macro endograft geometry, GCI and runs at rest and exercise, extraction at S1–S3 | TAWSS/OSI/RRT maps; S(t) drivers |
| 5–6 | `patch_geometry.py` and terrain-following mesher. V3, V4 and V5 validation. | Patch solver validated |
| 7 | Quasi-steady vs full (V6). Patch GCI. | Method locked; metrics frozen |
| 7–9 | OFAT screen (P1–P4 tested). Replicate the pilot +91% case with 5 seeds. | Factor shortlist |
| 9–11 | Face-centred CCD + platelets. Desirability, Pareto, confirmation, robustness. | Optimum plus 3 contrast designs chosen |
| 11–12 | Masters printed. PDMS inhibition fix. First PLGA casts. Metrology. | Coupons pass §7.4 |
| 13 | Mock incubation and sterilisation (§7.5). As-built CFD. Polymer decision. | Final coupon batch |
| 13–14 | SRC/IRB paperwork *submitted by week 8* so approval is in hand here. Assay pilot day (CV). | Final n chosen |
| 15–20 | 3 days of static platelets, 3 days of flow chamber, 3 days of biofilm | Raw data |
| 21–22 | Mixed-model statistics, CFD-vs-bench comparison (V8), figures | Results |
| 23–26 | Write-up and board. Buffer for repeats. | Final |

**Critical path:** fabrication (week 11) → bench (week 15). Start test prints of a sample master in week 6, in parallel with the CFD. Printer and PDMS problems are the most likely cause of delay.

---

## 10. Decision rules (pre-registered)

- **H1 (hemodynamics) is supported** if E-opt has Φ_TAWSS ≥ 0.95 and A_low,rel ≤ C0 + 5 pp at S1-rest *and* S1-exercise, with each metric's GCI < 5%.
- **H2 (platelets, CFD) is supported** if E-opt has lower deposition per planform mm² than C0 *beyond seed noise* (5 seeds, Welch t-test p < 0.05) *and* keeps that ranking under the ×0.1/×10 rate sensitivity.
- **H3 (platelets, bench) is supported** if the Dunnett contrast of E-opt vs C0 in the **flow-chamber** assay shows a reduction with a 95% CI excluding 0.
- **H4 (biofilm) is supported** under the same rule as H3, using CV OD590, with CFU agreeing in direction.
- **Recessed vs protruding:** the E-opt vs P-opt contrast is reported whatever its direction.
- **The negative result is a result.** Suppose no recessed design meets H1 at h/s ≥ 0.1. The quantitative conclusion is then: "At infrarenal shear, any groove deep enough to matter for platelets creates a measurable low-shear and long-residence zone." Report the h/s threshold at which this happens. That is the useful negative finding the proposal promises.

---

## 11. Risks and mitigations

| Risk | Mitigation |
|---|---|
| The printer can't resolve the optimum scale | Use scale invariance (§7.1) and fabricate at the smallest printable s. Ask the institution for micro-SLA. |
| PDMS will not cure on the printed masters | Venzac 2021 post-treatment, or a coated master. Test on scrap in week 6. |
| PLGA textures relax at 37 °C | Mock incubation (§7.5). Switch to polyurethane or PDMS. |
| The quasi-steady assumption fails | Full pulsatile runs. Cut the CCD to one orientation (15 designs). |
| The platelet model ranking flips under rate sensitivity | Report it honestly. Rely on bench H3 and on the hemodynamic H1 for the design choice. |
| No cleanroom access for 2 µm features | Diffraction-grating replica, labelled as a substitute. Or drop groups 5–6. |
| Biosafety approval is late | Submit paperwork by week 8. Run the CFD and fabrication regardless; they are complete projects on their own. |

---

## 12. Implementation checklist for `cfd/`

New modules (extending, not replacing, the coronary pipeline):
- `pipeline/womersley.py`: analytical solution, inlet profile writer, V2 comparison.
- `pipeline/endograft_geometry.py`: parametric bifurcated graft.
- `pipeline/patch_geometry.py`: periodic heightfield and undercut tiles, plus validation cavities and riblets.
- `pipeline/patch_mesh.py`: terrain-following blockMesh displacement.
- `pipeline/patch_case.py`: cyclic patch case, lid BCs, `scalarTransport` washout, and the protrusion-height post-processing.
- `pipeline/quasi_steady.py`: the 2-direction × 5-amplitude lookup, WSS(x, t) reconstruction and snapshot writer for `platelets.py`.
- `config/endograft.json`, `config/patch_experiment.json`, `config/waveform_infrarenal_rest.csv`, `config/waveform_infrarenal_exercise.csv`.
- `doe.py`: generates the CCD from a factor table, runs it, and fits the response surfaces, desirability and Pareto front.

---

## References

Items marked * are new (not in the proposal). Numbers in square brackets in the proposal map to the proposal's own list.

**Hemodynamics and blood**
- Cho YI, Kensey KR. Effects of the non-Newtonian viscosity of blood on flows in a diseased arterial vessel. *Biorheology* 1991;28:241–262.*
- Les AS et al. Supraceliac and infrarenal aortic flow in patients with abdominal aortic aneurysms: mean flows, waveforms, and allometric scaling relationships. *Cardiovasc Eng Technol* 2010. doi:10.1007/s13239-010-0004-8*
- Les AS et al. Quantification of hemodynamics in abdominal aortic aneurysms during rest and exercise using MRI and CFD. *Ann Biomed Eng* 2010. doi:10.1007/s10439-010-9949-x
- Malek AM, Alper SL, Izumo S. Hemodynamic shear stress and its role in atherosclerosis. *JAMA* 1999;282:2035–2042.
- Womersley JR. Method for the calculation of velocity, rate of flow and viscous drag in arteries when the pressure gradient is known. *J Physiol* 1955;127:553–563.*

**Creeping flow over textures (design rules and validation)**
- Higdon JJL. Stokes flow in arbitrary two-dimensional domains: shear flow over ridges and cavities. *J Fluid Mech* 1985;159:195–226.*
- Luchini P, Manzo F, Pozzi A. Resistance of a grooved surface to parallel flow and cross-flow. *J Fluid Mech* 1991;228:87–109.*
- Moffatt HK. Viscous and resistive eddies near a sharp corner. *J Fluid Mech* 1964;18:1–18.
- Shen C, Floryan JM. Low Reynolds number flow over cavities. *Phys Fluids* 1985;28:3191–3202. doi:10.1063/1.865366*

**Riblets and shark skin**
- Bechert DW et al. Experiments on drag-reducing surfaces and their optimization with an adjustable geometry. *J Fluid Mech* 1997;338:59–87.*
- Dean B, Bhushan B. Shark-skin surfaces for fluid-drag reduction in turbulent flow: a review. *Phil Trans R Soc A* 2010;368:4775–4806.
- García-Mayoral R, Jiménez J. Drag reduction by riblets. *Phil Trans R Soc A* 2011;369:1412–1427.*
- Wen L, Weaver JC, Lauder GV. Biomimetic shark skin: design, fabrication and hydrodynamic function. *J Exp Biol* 2014;217:1656–1666.* (Mako denticles are on the order of 150 µm, with surface features down to 5–10 µm.)

**Topography, platelets and bacteria**
- Chung KK et al. Impact of engineered surface microtopography on biofilm formation of *Staphylococcus aureus*. *Biointerphases* 2007;2:89–94.
- Koh LB, Rodriguez I, Venkatraman SS. The effect of topography of polymer surfaces on platelet adhesion. *Biomaterials* 2010;31:1533–1545.* (Textured PLGA made by solvent casting.)
- May RM et al. An engineered micropattern to reduce bacterial colonization, platelet adhesion and fibrin sheath formation for improved biocompatibility of central venous catheters. *Clin Transl Med* 2015. doi:10.1186/s40169-015-0050-9* (Sharklet on TPU. Note that the work was industry-funded.)
- Milner KR, Snyder AJ, Siedlecki CA. Sub-micron texturing for reducing platelet adhesion to polyurethane biomaterials. *J Biomed Mater Res A* 2006;76:561–570.*
- Perera-Costa D et al. Studying the influence of surface topography on bacterial adhesion using spatially organized microtopographic surface patterns. *Langmuir* 2014;30:4633–4641.*
- Sorensen EN et al. Computational simulation of platelet deposition and activation (Parts I and II). *Ann Biomed Eng* 1999;27:436–448 and 449–458.
- Zhang YQ et al. Genome-based analysis of virulence genes in a non-biofilm-forming *Staphylococcus epidermidis* strain (ATCC 12228). *Mol Microbiol* 2003;49:1577–1593.*

**Methods**
- Celik IB, Ghia U, Roache PJ, Freitas CJ. Procedure for estimation and reporting of uncertainty due to discretization in CFD applications. *J Fluids Eng* 2008;130:078001.*
- Derringer G, Suich R. Simultaneous optimization of several response variables. *J Qual Technol* 1980;12:214–219.*
- ibidi Application Note 11: Shear stress and shear rates for ibidi µ-Slides; and the sticky-Slide I Luer instructions.*
- Jones B, Nachtsheim CJ. A class of three-level designs for definitive screening in the presence of second-order effects. *J Qual Technol* 2011;43:1–15.*
- O'Toole GA. Microtiter dish biofilm formation assay. *J Vis Exp* 2011;(47):2437.*
- Venzac B et al. PDMS curing inhibition on 3D-printed molds: why? Also, how to avoid it? *Anal Chem* 2021;93:7180–7187.*

> Verify page numbers and volume details against the originals before submission. I confirmed the titles, journals and key findings, but not every page range.
