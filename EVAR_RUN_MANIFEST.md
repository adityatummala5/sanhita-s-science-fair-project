# DentiFlow EVAR: Run Manifest

The specific designs to build and run, with dimensions. Companion to `EVAR_EXPERIMENTAL_PLAN.md` (method) and `DentiFlow_Proposal.pdf` (rationale).

Read the column "What it answers" as the reason the run exists. If a run has no answer, cut it.

---

## 0. The crown parametrisation (read this first)

Your existing CAD crown is a **pointed shield**: 100 µm long × 70 µm wide × 43 µm tall, 3 riblets. For a pointed crown, height and face angle are the same number:

> face angle = arctan(h / (W/2)) = arctan(43/35) ≈ **51°** across the crown, ≈ 41° along it

Two consequences:

1. **You cannot vary height and face angle independently on a pointed crown.** A 30° crown 70 µm wide can only be 20 µm tall. Any OFAT sweep over height is secretly also a sweep over angle, and the +91% result cannot be attributed to either.
2. **Your pilot crown sits above the Moffatt threshold** (eddies when the face rises steeper than ≈ 34°). That is a candidate mechanism for the +91%, and it is testable.

So every parametric design below uses a **flat-topped frustum crown**:

```
            L_top  (plateau, streamwise)
         ┌──────────────┐          ─┬─
        /                \           h   relief
       /θ_up          θ_dn\         ─┴─
  ────┘                    └────         nominal wall
      |← h·cot θ_up →|
```

| Symbol | Meaning | Baseline value |
|---|---|---|
| `L_top` × `W_top` | plateau size (streamwise × spanwise) | 100 × 70 µm (matches your CAD) |
| `h` | relief (plateau top to groove floor) | 30 µm |
| `θ_up`, `θ_dn`, `θ_side` | face angles from the floor | 30° |
| `s` | groove width = gap between crown footprints | 100 µm |
| `n_r`, `d_r/h` | riblets on the plateau, depth fraction | 3, 0.25 |
| `φ` | orientation of the crown axis vs the flow | 0° (along flow) |
| `η` | net height (see below) | −1 (fully recessed) |
| layout | staggered diamond, 50% row offset | — |

Derived footprint: `L = L_top + h(cot θ_up + cot θ_dn)`, `W = W_top + 2h·cot θ_side`.
At baseline: **L = 204 µm, W = 174 µm**, lattice pitch 304 × 274 µm.

**Net height η** is the crown base position relative to the nominal wall, in units of h:

| η | Crown top | Groove floor | Name |
|---|---|---|---|
| +1 | +h | 0 (= wall) | protruding |
| +0.5 | +0.75h | −0.25h | mostly protruding |
| 0 | +0.5h | −0.5h | flush / straddling |
| −0.5 | +0.25h | −0.75h | mostly recessed |
| −1 | 0 (= wall) | −h | **engraved / recessed** |

Your example — "0.02 mm denticles embedded" — is **η = −1, h = 20 µm**. That is design `H-R20` below.

**Tile size snapping.** The patch tile is nominally 2.0 × 1.0 mm but is snapped to a whole number of lattice periods so the cyclic faces match. At baseline: 6 × 304 = 1824 µm by 4 × 274 = 1096 µm. Record the actual tile size per design in its JSON.

---

## 1. Stage A: macro endograft runs (no texture)

These set the driving conditions. **No texture appears in any Stage A run** — features are under 1% of the limb diameter, so they cannot change vessel-scale flow. Their only job is to produce the near-wall shear waveform S(t) that drives every patch run.

| ID | Geometry | Flow | Blood | What it answers |
|---|---|---|---|---|
| **A0** | Straight tube, Ø12 mm, L = 120 mm | Infrarenal waveform, first 8 harmonics | Newtonian, 3.5 mPa·s | Womersley validation (V2). Peak WSS within 5% of analytical. |
| **A1** | Bifurcated graft: 24 mm body × 60 mm, two 12 mm limbs × 90 mm, 40° total divergence, 2 mm crotch fillet | Rest, 1.05 L/min mean inlet | Carreau | **The reference case.** Maps TAWSS/OSI/RRT; picks sites S1, S2. |
| **A2** | Same as A1 | Exercise waveform | Carreau | Does the worst site move? Gives S3. |
| **A3** | Same as A1 but 60° divergence | Rest | Carreau | Geometry sensitivity. Does the site selection survive? |
| **A4** | A1 with limb crimp: transverse corrugation, pitch 1.5 mm, amplitude 0.35 mm | Rest | Carreau | **Crimp is a macro feature, not a patch feature** (see §2.1). Its hemodynamic cost belongs here. |
| **A5** | A1 at 3 mesh densities (r ≥ 1.3) | Rest | Carreau | GCI < 5% on limb TAWSS, OSI at S1, Δp. |

**Deliverable from Stage A:** for each of S1-rest, S2-rest, S1-exercise, a CSV of the wall-parallel velocity vector at H = 0.3 mm and the WSS vector, ≥ 200 samples per cycle. Plus the angle swept by the near-wall shear direction over the cycle — if it is under ±20°, the orientation series (§2.8) becomes a formality and can be cut to 2 designs.

---

## 2. Stage B: patch designs

**Held constant across every patch run unless the row says otherwise:**
tile 2.0 × 1.0 mm (snapped), lid height H = 0.3 mm, Carreau blood (ρ = 1060, μ₀ = 0.056, μ∞ = 0.00345, λ = 3.313, n = 0.3568), driving waveform from **S1-rest**, cyclic x and z, 5 shear amplitudes (0.02–3 Pa), identical mesh settings and refinement policy, identical metric definitions, fillet radius 0.05·s.

Metrics are always normalised to **C0 at the same site and amplitude**.

### 2.0 Validation designs (run these before anything else)

| ID | Geometry | What it answers |
|---|---|---|
| **V-CAV4** | 2D transverse rectangular cavity, W = 400 µm, h = 100 µm (W/h = 4) | Shen & Floryan: corner eddies only. |
| **V-CAV2** | Same, W = 200 µm, h = 100 µm (W/h = 2) | Shen & Floryan: single central vortex. The topology must change between CAV4 and CAV2. |
| **V-CAV1** | Same, W = 100 µm, h = 100 µm | Deep-cavity limit; stacked eddies. |
| **V-RIB** | Triangular riblet, s = 100 µm, h = 50 µm, included angle 90° | Luchini h∥ and h⊥ within 5% of tabulated. |
| **V-SCALE-A** | Baseline recessed crown at s = 100 µm, h = 30 µm | Scale invariance arm 1. |
| **V-SCALE-B** | Same shape at s = 300 µm, h = 90 µm | Arm 2. Normalised metrics must match arm 1 within the GCI. If they don't, Re_k is not small enough and §1.3 of the plan is wrong. |
| **V-MESH** | C0, H-R30 and Θ-70 at 3 mesh densities | Patch GCI < 5%. |

**Gate:** do not run the design matrix until V-CAV2 reproduces the vortex transition and V-RIB matches Luchini. These two prove the solver handles creeping flow over texture, which is the entire regime of interest.

### 2.1 Controls

| ID | Design | Dimensions | What it answers |
|---|---|---|---|
| **C0** | Smooth wall | flat | **Primary control.** The ePTFE-like surface. Every number is a ratio to this. |
| **C1** | Polyester **weave** texture | transverse quasi-sinusoid, pitch 250 µm, amplitude 50 µm peak-to-trough (yarn-bundle scale) | The clinical surface at patch scale. The texture to beat. |

> **Correction to the proposal.** It lists "native polyester crimp" as a patch-scale texture. Real graft crimp is pitch ≈ 1–2 mm, amplitude ≈ 0.3–0.5 mm — *larger than the 0.3 mm lid height and comparable to the 1 mm tile width*. It cannot live in the patch model; it would be driven by the lid, not by the shear. Crimp therefore moves to **run A4** (macro), and the patch-scale polyester control becomes the **weave ripple** (C1), which is the right size. Measure the real pitch and amplitude off a graft sample or SEM before fixing C1's numbers.

### 2.2 Height / placement series — the hypothesis test

Crown shape identical in all five. **Only the vertical position of the crown changes.** h = 30 µm, θ = 30°, s = 100 µm, 3 riblets at 0.25h, φ = 0°, diamond layout.

| ID | η | Crown top | Groove floor | What it answers |
|---|---|---|---|---|
| **H-P30** | +1 | +30 µm | 0 | Protruding. The pilot geometry's placement. |
| **H-P15** | +0.5 | +22.5 µm | −7.5 µm | |
| **H-F00** | 0 | +15 µm | −15 µm | Flush / straddling. |
| **H-R15** | −0.5 | +7.5 µm | −22.5 µm | |
| **H-R30** | **−1** | 0 | −30 µm | **Fully recessed. This is the hypothesis design (family E).** |

This series alone tests the central claim: *does removing protrusion remove the hemodynamic penalty?* Expect Φ_TAWSS to rise monotonically from H-P30 to H-R30 and to approach but not exceed 1.0 (plan rule P4).

Add **H-R20** (η = −1, h = 20 µm) and **H-R50** (η = −1, h = 50 µm) — your "0.02 mm embedded" case and a deeper one — so the recessed arm has three depths.

### 2.3 Depth series (h-arm): groove aspect ratio at fixed groove width

η = −1, s = 100 µm, θ = 30°, 3 riblets at 0.25h, φ = 0°.

| ID | h (µm) | h/s | Footprint W (µm) | Coverage | What it answers |
|---|---|---|---|---|---|
| **B-010** | 10 | 0.10 | 105 | high | Shallowest; should be nearly flat-equivalent. |
| **B-020** | 20 | 0.20 | 139 | | Below the cavity transition. |
| **B-030** | 30 | 0.30 | 174 | | **Baseline** (= H-R30). |
| **B-050** | 50 | 0.50 | 243 | | **At the Shen & Floryan transition.** Expect a step change in stagnant volume and t₉₀. |
| **B-075** | 75 | 0.75 | 330 | | |
| **B-100** | 100 | 1.00 | 416 | low | Deep-groove limit. Expect trapped recirculation. |

Note the footprint grows with h at fixed θ, so coverage fraction is not constant. That is intrinsic — a deeper groove at a fixed face angle *is* a wider crown. Record coverage as a covariate and use the s-arm below to separate the two.

### 2.4 Width series (s-arm): the same h/s range reached the other way

η = −1, h = 30 µm, θ = 30°, 3 riblets at 0.25h, φ = 0°. Footprint fixed at 204 × 174 µm; only the lattice pitch changes.

| ID | s (µm) | h/s | Pitch W (µm) | What it answers |
|---|---|---|---|---|
| **S-300** | 300 | 0.10 | 474 | |
| **S-150** | 150 | 0.20 | 324 | |
| **S-100** | 100 | 0.30 | 274 | Baseline (= H-R30). |
| **S-060** | 60 | 0.50 | 234 | Transition, reached by narrowing. |
| **S-030** | 30 | 1.00 | 204 | Crowns nearly abutting (your CAD's 85 µm pitch is in here). |

**Why both arms.** If the normalised metrics from §2.3 and §2.4 collapse onto one curve against h/s, then h/s is the governing parameter and the DOE can use it as a single factor. If they don't, coverage fraction matters independently and needs its own factor. This is a cheap, decisive test and it is worth 11 hemodynamic runs.

### 2.5 Face-angle series — the Moffatt test

η = −1, h = 30 µm, s = 100 µm, symmetric faces (θ_up = θ_dn = θ_side), 3 riblets at 0.25h, φ = 0°.

| ID | θ (°) | Interior corner angle | Moffatt prediction | What it answers |
|---|---|---|---|---|
| **Θ-20** | 20 | 160° | no eddy | |
| **Θ-30** | 30 | 150° | no eddy (marginal) | Baseline. |
| **Θ-35** | 35 | 145° | **first eddy** | **The predicted threshold crossing.** |
| **Θ-40** | 40 | 140° | eddy | |
| **Θ-55** | 55 | 125° | eddy | |
| **Θ-70** | 70 | 110° | strong eddy | |
| **Θ-51** | 51 | 129° | eddy | **Your CAD crown's actual face angle.** Explains the pilot. |

Expect a near-discontinuity in stagnant volume and t₉₀ between Θ-30 and Θ-40. If you see it, P1 is confirmed and you have a publishable mechanism for the +91%. If you don't, the mechanism is something else and the asymmetry series becomes the next suspect.

### 2.6 Asymmetry series — does flow reversal punish asymmetric crowns?

η = −1, h = 30 µm, s = 100 µm, 3 riblets at 0.25h, φ = 0°.

| ID | θ_up | θ_dn | What it answers |
|---|---|---|---|
| **A-2020** | 20° | 20° | Symmetric shallow. (≈ Θ-20) |
| **A-2055** | 20° | 55° | Shallow into the flow, steep behind. |
| **A-5520** | 55° | 20° | **Steep into the flow.** The pilot's "undercut crown with steep upstream face" minus the undercut. |
| **A-5555** | 55° | 55° | Symmetric steep. (≈ Θ-55) |

Under OSI > 0 the flow reverses, so A-2055 and A-5520 are the *same shape* seen in opposite directions. If they give different platelet deposition, the difference is driven by the forward/reverse flow asymmetry in the waveform, not by shape — a useful and non-obvious result. If they give the same answer, asymmetry does not matter and θ_up/θ_dn can collapse to one factor in the DOE.

### 2.7 Riblet series

η = −1, h = 30 µm, s = 100 µm, θ = 30°, φ = 0°. Riblets run streamwise along the 70 µm plateau.

| ID | n_r | d_r/h | d_r (µm) | Riblet pitch | Groove width | What it answers |
|---|---|---|---|---|---|---|
| **R-0** | 0 | 0 | 0 | — | — | Smooth plateau. The null case. |
| **R-3L** | 3 | 0.25 | 7.5 | 23 µm | ~12 µm | Baseline (= H-R30). |
| **R-3D** | 3 | 0.50 | 15 | 23 µm | ~12 µm | Deeper riblets. |
| **R-5L** | 5 | 0.25 | 7.5 | 14 µm | ~7 µm | Finer, approaching platelet scale. |
| **R-5D** | 5 | 0.50 | 15 | 14 µm | ~7 µm | Finest and deepest. |

The pilot found riblets "slightly beat" smooth crowns. This series says whether that survives at EVAR shear and whether it is a count effect or a depth effect. **Flag:** R-5* needs ~0.5 µm cells to mesh, and 7 µm grooves need the micro-SLA route to fabricate. If the plateau riblets turn out to matter, that drives the fabrication decision.

### 2.8 Orientation series

η = −1, h = 30 µm, s = 100 µm, θ = 30°, 3 riblets at 0.25h. Only the crown axis rotates.

| ID | φ | What it answers |
|---|---|---|
| **O-00** | 0° (along flow) | Baseline. Grooves are open channels; P3 predicts no recirculation. |
| **O-45** | 45° | Intermediate. |
| **O-90** | 90° (transverse, crimp-like) | Grooves are cavities. P3 predicts the worst stagnation. |

**Cut this to O-00 and O-90 if Stage A shows the near-wall shear direction rotating less than ±20°** over the cycle. With a fixed shear direction, P3 makes aligned grooves win trivially and the 45° point adds nothing.

### 2.9 Arrangement series

η = −1, h = 30 µm, s = 100 µm, θ = 30°, φ = 0°.

| ID | Design | What it answers |
|---|---|---|
| **G-DIA** | Staggered diamond crowns (baseline) | Segmented, biomimetic layout. |
| **G-LIN** | Aligned rectangular grid of the same crowns | Does row stagger matter, or only groove geometry? |
| **G-CONT** | Continuous streamwise riblets, same h and s, no crown segmentation | **The groove-flushing test.** Continuous channels cannot trap a closed cavity; segmentation introduces transverse steps between crowns. Isolates "is segmentation the problem?" |

G-CONT vs G-DIA is the cleanest available test of the proposal's "continuous vs segmented riblets" question.

### 2.10 Undercut designs (CFD only — not fabricable)

| ID | Design | What it answers |
|---|---|---|
| **U-BAD** | Reconstruction of the pilot's **+91%** crown: pointed shield 100 × 70 × 43 µm, undercut crown, steep upstream face, η = +1 | **Replicate the pilot result at EVAR shear, with 5 seeds.** Is +91% real or one seed's aggregation? |
| **U-GOOD** | Reconstruction of the pilot's **−24%** crown (D4_dense, 41.8 µm, 65 µm spacing), η = +1 | The paired comparison. Two equal-height designs, opposite outcomes. |

These two exist to explain the preliminary result, not to be built. Per `cfd/results/platelet_comparison/`, the +91% was one random realisation and the aggregated count was +551% with two clusters dominating — rerun with 5 seeds and report both the full and aggregation-excluded metrics before drawing any mechanism from it.

### 2.11 Patch run count

| Block | Unique designs | Platelets? |
|---|---|---|
| Validation (§2.0) | 7 | no |
| Controls (§2.1) | 2 | yes, 3 seeds |
| Height (§2.2) | 5 + 2 = 7 | yes, 3 seeds |
| Depth h-arm (§2.3) | 5 new | no |
| Width s-arm (§2.4) | 4 new | no |
| Angle (§2.5) | 6 new | yes, 3 seeds on Θ-30/35/40 |
| Asymmetry (§2.6) | 2 new | yes, 3 seeds |
| Riblets (§2.7) | 4 new | yes, 3 seeds on R-0, R-3L |
| Orientation (§2.8) | 2 new | yes, 3 seeds |
| Arrangement (§2.9) | 2 new | yes, 3 seeds |
| Undercut (§2.10) | 2 | yes, **5 seeds** |
| **OFAT total** | **43 runs, ~36 unique designs** | ~60 platelet runs |

---

## 3. Stage B2: the response surface

Run **after** the OFAT screen, with factors trimmed to whatever survived.

**Factors** (assuming η, h/s, θ_up survive, which is the expected outcome):

| Factor | Low (−1) | Centre (0) | High (+1) |
|---|---|---|---|
| η | −1 (recessed) | 0 (flush) | +1 (protruding) |
| h/s | 0.10 | 0.30 | 0.50 |
| θ_up | 20° | 35° | 50° |

**Design:** face-centred CCD, α = 1 → 8 factorial + 6 face + 1 centre = **15 points**, run at φ = 0° and φ = 90° → **30 designs**, 3 platelet seeds each → **90 platelet runs** (~90 CPU-h).

Fixed throughout: s = 100 µm, L_top × W_top = 100 × 70 µm, θ_dn = θ_up, 3 riblets at 0.25h, diamond layout, S1-rest driving.

Note the h/s high level is capped at 0.5, not 1.0 — the OFAT depth series will have shown whether anything above the cavity transition is viable. If B-075 and B-100 are disasters, there is no reason to spend CCD points there.

**Confirmation:** the predicted optimum, on the fine mesh, 5 seeds. Plus the top 5 designs re-evaluated at S1-exercise and S2.

---

## 4. Bench coupons

Seven surfaces. Every one is cast through the **same** chain (print → PDMS negative → PLGA) so that material, solvent history and bulk roughness are identical and only the designed texture differs.

| ID | Surface | Dimensions | Groups it is in | What it answers |
|---|---|---|---|---|
| **BC-0** | Smooth | flat | platelets, biofilm | Primary control. |
| **BC-1** | Weave ripple | transverse, pitch 250 µm, amp 50 µm | platelets, biofilm | Clinical comparator. |
| **BC-E** | **CFD optimum, recessed** | from §3; expect η = −1, h/s ≈ 0.2–0.3, θ ≈ 20–30° | platelets, biofilm | **The hypothesis.** |
| **BC-P** | **Same crown, protruding** | identical shape, η = +1 | platelets | Recessed vs protruding, held at identical shape. The pilot's lesson, tested directly. |
| **BC-S** | Cell-scale ribs | ~2 µm ribs, diamond array | platelets, biofilm | Contact-area mechanism (Sharklet-type). |
| **BC-H** | Hierarchical | BC-E with BC-S ribs on the plateaus | platelets, biofilm | Do the two mechanisms add? |
| **BC-X** | **CFD's worst design** | whichever scored worst in §3 | platelets | **Does the CFD ranking predict the bench ranking?** Spearman ρ over BC-0/1/E/P/X. This is the run that validates the whole computational half of the project. |

**Replication:** 3 independent days × 4 coupons per surface per day = 12 per surface. Run a 1-day pilot on BC-0 alone first to measure the CV and set the final n.

**If capacity is short, cut in this order:** BC-H, then BC-S, then BC-1. Never cut BC-P or BC-X — they are the two comparisons that carry the argument. Never cut replicates to add groups.

**Per-batch acceptance** (§7.4 of the plan): h and s within ±10% of CAD, θ within ±5°, measured by profilometry on ≥ 3 coupons. Then re-run the CFD on the as-built profile, so the bench is compared against the geometry you actually made rather than the one you designed.

---

## 5. Build order

1. **A0** → Womersley passes. Nothing else is trustworthy until this does.
2. **V-CAV4 / V-CAV2 / V-RIB** → patch solver validated against literature.
3. **A1, A5** → S1 located, GCI clean, drivers extracted. Check the shear-direction rotation and decide on §2.8.
4. **V-SCALE-A/B, V-MESH, V6 (quasi-steady vs full)** → method locked, metrics frozen.
5. **C0, C1** → the baselines everything is normalised to.
6. **§2.2 height series** → the single most informative block. If recession does not recover Φ_TAWSS here, the hypothesis is in trouble and you will want to know in week 7, not week 20.
7. **§2.10 undercut pair, 5 seeds** → settle the +91%.
8. **§2.3–2.9** → the rest of the screen.
9. **A2, A3, A4** → robustness and the macro crimp cost.
10. **§3 CCD** → optimum.
11. **Masters and coupons** (start test prints in parallel at step 6, not at step 11).

Steps 1–4 are gates. Each one can fail and each failure changes the plan, so do not run the design matrix before they pass.
