# Plan: Textured stent-graft surfaces in a simulated abdominal aortic EVAR

**Status:** proposal, 2026-11-28. Nothing in this plan has been run yet.
**Replaces:** the coronary stent setup. The coronary results (`cfd/results/`) are kept as the **pilot study** that motivates the design choices below.

---

## 0. Bottom line

- **Moving to the aorta changes the physics.** The vessel is about 8× wider and the flow is still laminar at rest, but more pulsatile, and it **reverses direction** for part of each heartbeat. The device is also different: a fabric graft with a large luminal surface, not thin metal struts.
- **Shark-scale texture (~50–100 µm) is tiny compared with the near-wall flow scale in the aorta.** It will barely change the global flow. Any effect will be local, on where platelets stick. So we need a **two-scale approach** (section 3).
- **My expectation, written down now:** no micro-texture will beat a smooth graft by more than ~10% on platelet adhesion. The one candidate with a real chance of improving the hemodynamic metrics is a **helical ridge** that sets up swirling flow in the iliac limbs. That idea is only loosely shark-inspired, and the clinical evidence for it is mixed.
- **The design is built so that either outcome is a result.** Every texture is paired with a control, predictions are written before running, and a known-bad texture is included to show the method can detect harm.

---

## 1. What the coronary pilot taught us

| # | Finding (coronary, coarse mesh, 1 seed) | What it means for EVAR |
|---|---|---|
| L1 | Hemodynamic cost grew with denticle height | Keep features low, or make them recessed rather than raised |
| L2 | Separate bumps beat a uniformly thicker surface: flow passes between them | Prefer discrete or streamwise-open features over continuous transverse ones |
| L3 | Same height, opposite outcome: Sanhita's undercut denticles gave 2.7× adhesion, D4 was neutral | Shape matters. Test the undercut directly, and keep it as a negative control |
| L4 | Spacing, smoothing and orientation (D5–D8) changed very little | Don't spend runs on fine parameter sweeps of the same denticle |
| L5 | Laminar flow leaves riblets nothing to suppress | Same issue in the aorta (section 2). Test it directly rather than assume it |
| L6 | Aggregation was dominated by chance runaway clusters | The primary metric is direct adhesion. Use ≥3 seeds per design |
| L7 | Every stented mesh failed `checkMesh` | Fix mesh quality before any EVAR result counts |

---

## 2. The EVAR flow environment

Values are typical. Before the final runs, replace them with numbers from one cited source paper, as was done for the coronary LAD case.

| Quantity | Value | Basis |
|---|---|---|
| Infrarenal mean flow, rest | ~17.5 mL/s ≈ 1.05 L/min | Les et al. 2010, 43 AAA patients, PC-MRI |
| Waveform | Infrarenal resting: systolic peak, **early-diastolic reverse flow**, low diastolic flow | Digitise from Les et al. 2010 (or Moore & Ku 1994) |
| Period | 0.9 s | Same as the pilot |
| Graft main body | Ø 26 mm, ~60 mm long, tapering into the bifurcation | Idealised; adjust to a published EVAR CFD geometry |
| Iliac limbs | Ø 13 mm, ~80 mm long, ~40° total bifurcation angle | Idealised |
| Aneurysm sac | Excluded: the flow domain is the graft lumen only | Standard for EVAR CFD |
| Blood | Carreau, ρ = 1060 kg/m³: same constants as the pilot | Keeps the two studies comparable |
| Walls | Rigid, no-slip | Limitation, as in the pilot |

**Derived numbers** (my calculations from the values above):

| | Main body (Ø 26 mm) | Each limb (Ø 13 mm) |
|---|---|---|
| Mean Reynolds number | ~260 | ~260 |
| Peak Reynolds number (≈5× mean flow) | ~1,300 | ~1,300 |
| Womersley number | ~19 | ~9 |
| Near-wall oscillating layer (Stokes layer) | ~1 mm | ~1 mm |
| Viscous length ν/u_τ at a wall shear of 1 Pa | ~100 µm | ~100 µm |

What these imply:

- **Still laminar at rest.** Shark riblets reduce drag when spaced at ~15 viscous lengths in *turbulent* flow. Here that spacing would be ~1.5 mm, and there is no developed turbulence to act on. The shark-skin mechanism remains mostly irrelevant, as in the pilot.
- **Shark-scale denticles (~40 µm) are under one viscous length tall.** The bulk flow sees them as smooth. Their effect is purely local, on the near-wall flow that platelets travel through.
- **Flow reversal matters for directional shapes.** A denticle that points downstream points upstream during diastolic backflow. Symmetric or streamwise shapes avoid this.
- **The graft wall has low shear.** Much of a smooth graft will already sit below the 0.4 Pa low-shear threshold used in the pilot. So also report a **relative** threshold: area below 50% of the smooth-wall TAWSS at the same site.

---

## 3. Method: two scales

A 50 µm texture cannot be meshed across a 26 mm × 150 mm graft. Split the problem in two:

**Global model (whole idealised EVAR, smooth or with mm-scale features)**
- Purpose: realistic flow, a map of where the graft is at risk (low TAWSS, high OSI and RRT), and the near-wall shear waveform τ(t) at chosen sites.
- Also the place to test **macro-textures** (helical ridges, ~0.5–1 mm), which change the global flow.

**Patch model (a ~3 × 2 × 2 mm block of graft wall with texture on the floor)**
- Periodic front/back and left/right, so the texture repeats as it would across a whole fabric surface.
- Driven by a time-varying flow rate, tuned so that the **smooth** floor reproduces the τ(t) taken from the global model.
- Cheap (~0.5 M cells), and it resolves textures at 5–10 µm.
- **Its geometry is the same as a flat test coupon.** Any later bench test can use exactly the simulated geometry, which fixes issue 9 in `PROJECT_SUMMARY.md`.

**Patch sites**, chosen from the global run:
- **Site A (primary):** the worst graft location, expected in the iliac limb near the flow divider. Clinical CFD links low iliac TAWSS to in-graft thrombus.
- **Site B (robustness check):** main body mid-wall.

---

## 4. Textures to test, ranked by expected success

Here "success" means beating the smooth graft on the primary endpoint (section 5). The odds are my judgment from the pilot and the literature, not data.

| Rank | ID | Texture | Model | Key dimensions | Why it's in | Prediction vs smooth | Odds of success |
|---|---|---|---|---|---|---|---|
| 1 | **G1** | **Helical ridge** in the iliac limbs: a single ridge that induces swirl | Global | Height 0.5 and 1.0 mm (~4–8% of limb Ø); pitch ~3D and ~6D | Swirl raises WSS and reduces stagnation in CFD studies of spiral-flow bypass grafts; L2 (streamwise-open) | Higher limb TAWSS, smaller low-shear area; small pressure-drop penalty; a low-shear pocket in the ridge's wake | **Fair–good** for hemodynamics; **fair** for platelets |
| 2 | P5 | **Plateau-and-groove / shallow dimples**: mostly flat surface with narrow recessed streamwise grooves or shallow dimples (depth/diameter < 0.1) | Patch | Groove 30 µm wide × 20 µm deep, 150 µm pitch; or dimples Ø150 × 15 µm | L1: nothing protrudes; most of the surface keeps full shear | Within ±10%; deposits concentrate in grooves | Low–fair |
| 3 | P1 | **Shark-scale streamwise riblets**: scalloped profile, aligned with flow | Patch | s = 100 µm, h = 50 µm (h/s = 0.5) | The classic shark-skin texture; symmetric under flow reversal | Within ±10%; low shear in the valleys | Low |
| 4 | P3 | **Low parametric denticle** (D3-like, 25 µm, no undercut), run forward and reversed | Patch | Pilot `DenticleParams`, h = 25 µm | Replicates the pilot's neutral D3/D4 result; tests sensitivity to flow reversal | Within ±10%; reversed ≈ forward because flow reverses anyway | Low |
| 5 | P2 | **Flow-scaled riblets**: same shape as P1, but spaced at ~15 viscous lengths | Patch | s = 1.0 mm, h = 0.5 mm | The best case for the shark mechanism: correct spacing in wall units. P1 vs P2 isolates scale | Neutral to worse: valleys are big enough to hold slow flow | Low |
| ref | P0b | **Woven polyester weave**: approximate yarn-bundle relief | Patch | Amplitude ~20–40 µm, pitch ~200–300 µm (confirm from published SEM images) | Real grafts are already textured. Is any design better than what is used clinically? | Slightly worse than smooth | Reference |
| neg | P8 | **Transverse ridges**: same profile as P1, rotated 90° (like fabric pleats or crimps) | Patch | s = 100 µm, h = 50 µm | Orientation control for P1; shows the method can detect a bad texture | Clearly worse: recirculation in every groove | Very low (negative control) |
| neg | P6 | **Sanhita's undercut denticle**, re-meshed onto the patch floor | Patch | From her CAD, 43 µm | Does the pilot's strongest result (L3) replicate in a different flow? | Clearly worse | Very low (negative control) |
| — | P0 | **Smooth graft** (ePTFE-like) | Both | — | **Control** | — | — |

**Mechanism study (if time allows):**
- **P7, undercut parametric denticle.** Same height and footprint as P3, but the crown is wider than its stalk. If P7 ≈ P6 and P3 ≈ smooth, the undercut is the cause.

**Deliberately excluded from simulation:**
- **Sub-platelet microtexture (~2 µm, Sharklet-like).** Its proposed benefit is reduced platelet and bacteria *contact area*. The platelet model treats platelets as points, so it cannot capture that. This is the texture class with the best non-hemodynamic case, so keep it for the wet-lab phase if one happens.
- **Fine sweeps of denticle spacing, smoothing and orientation.** The pilot showed small effects (L4).

---

## 5. Endpoints and decision rule

Fix these **before** running.

**Primary endpoint (patch):** platelets directly adhered per mm² of **projected** wall area after 3 cardiac cycles, mean of 3 seeds.
- Projected area, not wetted area, because a texture that adds surface but keeps adhesion per wetted area constant still puts more platelets on the patient's graft.

**Secondary endpoints:**
- **Patch:** TAWSS, OSI, RRT; absolute (< 0.4 Pa) and relative (< 50% of smooth) low-shear area fractions; fraction of deposits in valleys or grooves.
- **Global:** TAWSS, OSI and RRT maps on the graft; low-shear area per region (body, flow divider, limbs); helicity (LNH); pressure drop; **displacement force on the graft**. Migration is a real EVAR failure mode, so a ridge must not increase drag force meaningfully.

**A texture counts as a success only if all of these hold:**
1. Primary endpoint ≥ 10% below smooth, with the difference > 2 SD of the seed-to-seed spread.
2. The result holds on the medium mesh: the endpoint changes < 10% from coarse to medium.
3. The relative low-shear area is not worse than smooth.
4. It holds at Site A, and at Site B is at least not worse.
5. **G1 only:** pressure drop +< 5% and displacement force +< 5% vs the smooth graft.

**Written predictions (record them in the project notebook now):**
- G1 is the only design to pass.
- P1, P3 and P5 land within ±10% of smooth.
- P2 and P0b come out slightly worse.
- P8 and P6 come out clearly worse.

---

## 6. Run order

| Stage | Runs | Purpose | Rough compute (4 cores)* |
|---|---|---|---|
| **0. Validation** | Smooth patch vs the exact oscillating-channel (Womersley) solution; straight smooth tube vs the exact Womersley pipe solution | Prove both models give known answers (target WSS error < 3%) | ~4 h |
| **1. Global baseline** | Smooth EVAR at coarse/medium/fine mesh, 3 cycles | Mesh independence; risk map; τ(t) at Sites A and B; compare TAWSS ranges with published EVAR CFD | ~20 h |
| **2. Patch screen** | P0, P0b, P1, P2, P3 (forward and reversed), P5, P6, P8 at Site A; flow plus platelets, 1 seed | First ranking, about 9 runs | ~30 h |
| **3. Confirmation** | Top 2 patch textures + P0 + P6: medium mesh, 3 seeds, plus Site B | Apply the decision rule | ~25 h |
| **4. Macro texture** | G1 at 2 heights (best pitch from a 2-pitch pre-screen) vs smooth global; platelets in the limb | Test the most promising design | ~25 h |
| **5. Mechanism (optional)** | P7 undercut vs P3 vs P6 | Explain the shape effect | ~8 h |

\*Scaled from the pilot's measured cost (0.35 M cells, 22 min per cycle). Treat these as ±2× estimates. The total is **~110 h**: about a week of continuous running on one 4-core machine, or a few days on a lab workstation. Stages 2 and 4 can run side by side.

---

## 7. Code work needed (existing pipeline in `cfd/`)

| Piece | Change |
|---|---|
| `pipeline/geometry.py` | New: idealised bifurcated graft generator (smooth, plus helical-ridge option). New: periodic patch generator with texture library (riblets, grooves, dimples, weave, transverse ridges). Reuse `DenticleParams` for P3, P6 and P7 |
| `pipeline/case.py` | Patch case with periodic boundaries and time-varying flow-rate forcing. Global case with Womersley inlet (or a long entrance length) and a flow split at the two outlets |
| `pipeline/waveform.py` + `config/` | Infrarenal waveform CSV, digitised from the chosen source paper |
| `pipeline/platelets.py` | Periodic wrap-around for the patch; the whole graft surface is reactive (rate constants stay collagen-based, so only *relative* comparisons are valid) |
| `pipeline/post.py` | Projected-area metrics, relative low-shear threshold, helicity, displacement force |
| Mesh quality | Find out why the stented meshes failed `checkMesh` and fix it before Stage 1 |

---

## 8. Limitations to state up front

- Rigid graft walls; idealised (not patient-specific) geometry; no aneurysm sac flow or endoleaks.
- Laminar model. It cannot capture transitional or turbulent bursts, so if a texture only helps under turbulence these runs would miss it. Resting infrarenal flow is not turbulent, so this is a minor gap.
- The patch model assumes the texture doesn't change the large-scale flow. That is justified for features under ~100 µm, but not for G1, which is why G1 runs in the global model.
- Platelet model: collagen rate constants, simplified activation, no shear-induced activation, no detachment, a short simulated window. It captures *initial* adhesion only.
- Contact-area effects of microtexture are not modelled (section 4, excluded).

---

## 9. Document changes this plan triggers

- **Target vessel:** now an **abdominal** aortic stent graft (EVAR). The title brainstorm in `2027_proj.docx` was thoracic. Pick abdominal everywhere.
- **"Stent" vs "stent graft":** resolved. EVAR devices are stent grafts, and the texture sits on the graft's luminal fabric surface.
- **Research question, suggested rewrite:** *"Can bio-inspired surface textures on an abdominal aortic stent graft reduce platelet adhesion and low-shear flow, and which texture features help or hurt?"*
- **Hypothesis:** rewrite it to match the predictions in section 5.
- **Coronary work:** present it as the pilot that produced L1–L7.
