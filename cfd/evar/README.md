# DentiFlow EVAR pipeline: textured luminal surfaces for endograft limbs

This folder holds the computational part of the EVAR project (see `proposal/`). Its question is whether a denticle-inspired texture on the inner surface of an EVAR stent graft reduces platelet deposition without worsening the slow, oscillating near-wall flow.

The vessel is 12–25 mm across, but the textures are 10–100 µm, so one mesh cannot resolve both. The work is split into two linked scales:

```
config/evar.json
   │
   ├─► run_macro.py validation   25 mm straight tube vs the exact Womersley solution  (method check)
   │
   ├─► run_macro.py graft        idealized bifurcated endograft, Carreau blood, pulsatile infrarenal flow
   │        └─ post_macro.py     TAWSS / OSI / RRT maps, per-segment stats, patch-site waveforms
   │                                   │  (shear-stress waveform at the worst limb site)
   │                                   ▼
   └─► run_patch.py              periodic textured wall patch, one texture unit cell, Stokes flow
            ├─ --stage solve     steady solves at 6 shear levels per design
            ├─ --stage post      cardiac cycle rebuilt with the site waveform -> surface metrics
            ├─ --stage verify    full time-dependent run of one design vs the quasi-steady rebuild
            └─ --stage platelets Lagrangian platelet adhesion on the reactive graft surface
```

All commands are run from `cfd/`. OpenFOAM v2406 runs in Docker, the same as for the coronary pipeline (see `cfd/README.md`).

## 1. Vessel scale (`run_macro.py`, `graft_geometry.py`, `macro_case.py`, `post_macro.py`)

**Graft geometry.** The lumen is an implicit surface turned into a mesh with marching cubes:

| Segment | Dimension |
|---|---|
| Main body | 25 mm |
| Legs (bifurcated-body outlets) | 12 mm, smooth transition over 15 mm |
| Flare | 12 → 16 mm |
| Iliac limbs | 16 mm, diverging at 8° each |

- These follow the Endurant-based idealized model of *Frontiers in Physiology* 2021 (doi:10.3389/fphys.2021.795085).
- Only a quarter of the lumen is meshed, with symmetry planes at y = 0 (between the limbs) and z = 0.

**Flow.**

| Parameter | Setting |
|---|---|
| Blood | Carreau (Cho & Kensey 1991), density 1060 kg/m³ |
| Walls | Rigid, no-slip graft fabric |
| Outlets | p = 0, backflow allowed |
| Inlet profile | Womersley, from the waveform's Fourier modes |
| Mean flow | 0.9 L/min resting infrarenal (Taylor et al., *Ann Biomed Eng* 2002;30:402) |
| Period | 1 s |
| Waveform shape | **Generic triphasic** (systolic peak, early-diastolic reversal): see the caveat below |
| Peak flow (75 mL/s) | An assumption |

> ⚠️ Replace `config/waveform_infrarenal_rest_approx.csv` with a waveform digitised from a published MRI study (e.g. Les et al. 2010) before final runs. The columns are `time_s,flow_mL_s`, and the pipeline rescales it to the mean and peak in the config.

**Solution procedure.**
1. A steady mean-flow run initialises the mean component, which develops slowly in a large vessel.
2. Three pulsatile cycles follow (pimpleFoam, maxCo 3, 2 outer correctors), and metrics are averaged over the last cycle.
3. Wall-shear probes check how much the result changes from one cycle to the next.

**Validation** (`results/evar/validation-coarse/`). A 25 mm straight tube with Newtonian blood, run on the same mesh settings, was compared with the analytical Womersley solution (Wo = 17):

| Quantity | CFD | Exact | Error |
|---|---|---|---|
| Wall shear stress over the cycle | | | **1.05% RMS** of peak |
| TAWSS | 0.191 Pa | 0.187 Pa | +2.4% |
| OSI | 0.411 | 0.409 | |
| Peak WSS | 0.970 Pa | 0.972 Pa | |

**Outputs** (`results/evar/graft-<preset>/`):

| File | Content |
|---|---|
| `graft_metrics.json` | TAWSS, low-WSS area, OSI and its > 0.25 / > 0.45 areas, RRT median / 95th percentile, per segment (body, bifurcation, leg, flare, limb, each split into inner and outer wall) |
| `graft_maps.png` | Maps of those metrics over the graft |
| `site_worst_limb_site.csv`, `site_median_limb_site.csv` | Signed shear-stress waveform along the principal shear direction at the highest-RRT and the median limb site. These drive the patch model. |

## 2. Texture scale (`textures.py`, `patch.py`, `run_patch.py`, `post_patch.py`)

**Domain.** A periodic box covering one texture unit cell:
- **Bottom:** the textured graft surface y = h(x, z).
- **Top:** a lid at 0.3 mm that moves at U(t) (Couette drive).

| Quantity | Value | Implication |
|---|---|---|
| Texture Reynolds number, γh²/ν | ~0.01–0.1 | Stokes flow |
| Patch Womersley number, H√(ω/ν) | ~0.3 | Quasi-steady |
| Lid height vs limb boundary layer | 0.3 mm vs ~1 mm | Near-wall velocity is close to linear |
| Lid height vs texture disturbance | 0.3 mm vs period/2π | The lid sees an undisturbed flow |

**Solving.** Each design is solved steadily at 6 shear levels (γ = 0.3 to 400 s⁻¹), and the cardiac cycle is rebuilt by interpolating between levels using the site waveform. `--stage verify` checks this shortcut against a full time-dependent run.

**Calibration.** For the smooth wall, the lid speed is set from the Carreau law so that it reproduces the site's wall stress exactly; CFD matches the exact Couette stress to 8 significant digits. Every textured design is driven by the same U(t), so all designs see the same outer flow.

**Mesh.**
- Structured hex cells that follow the surface, 3 µm horizontally, with 40 graded cells across the gap (the first cell is about 1 µm).
- Periodicity is exact, because the heightfield is periodic. Every design passes checkMesh.

**Designs** (`config/evar.json` → `patch.designs`). Each texture is placed in one of three ways (`placement`):
- **protrude:** stands on the surface.
- **flush:** the crown tops are level with the surface and the field between them is cut down (height-neutral).
- **engraved:** the crown shape is cut in as a pit.

| Design | What it tests |
|---|---|
| `smooth` | Control: untextured graft (ePTFE-like) |
| `pet` | Control: idealized woven-polyester yarn ridges (the clinical comparator) |
| `den_protrude`, `den_flush`, `den_engraved` | Placement of the 43 µm denticle crown |
| `den_flush_short` | Depth (20 µm, a lower depth-to-spacing ratio) |
| `den_flush_sparse` | Spacing (150 × 100 µm lattice) |
| `den_flush_reverse`, `den_flush_transverse` | Orientation relative to the oscillating flow |
| `riblets_deep`, `riblets_shallow` | Continuous streamwise riblets, depth-to-spacing ratio 0.5 and 0.25 |

Denticle crowns use the coronary pipeline's shape (100 × 70 µm, 3 riblets), with softened side walls so that the surface-following mesh stays valid. Undercut crowns, which the coronary pilot flagged, cannot be represented as a heightfield. They are a known gap.

**Metrics.** `results/evar/patch/site-<site>/patch_metrics.csv`, plus `patch_comparison.png` and `patch_tawss_maps.png`:
- TAWSS, OSI and RRT.
- Low-WSS area (< 0.4 Pa), and the area below 50% of the smooth wall's TAWSS.
- Areas with OSI > 0.25 and > 0.45.
- RRT median and 95th percentile.
- WSSG: the 95th percentile of the surface gradient of TAWSS.
- Slow-flow volume: fluid per mm² of graft moving slower than fluid 5 µm above a smooth wall.
- Stress per planform area, and wetted-area ratio.

**Platelets** (`patch_platelets.py`). The kinetics are the same as in the coronary pilot (Sorensen 1999; shear-enhanced diffusivity from Zydney & Colton). The setup differs in three ways:
- The **whole graft surface is reactive**.
- Platelets are tracked in the 70 µm layer next to the wall, with the bulk platelet count held fixed at the top of that layer.
- Each particle is one real platelet, in a tiled 0.5 × 0.5 mm domain.

Each design is run for 3 seeds and 30 s of real time. Results are reported as mean ± SD per mm² and compared with the smooth graft.

## 3. Limitations (state these in the write-up)

1. The waveform shape is generic until it is replaced with digitised data, and the peak flow is an assumption.
2. The graft is idealized, rigid and symmetric, and is simulated as a quarter model. Its fabric is smooth at the vessel scale; crimps and stent rings are not modelled.
3. The patch is driven along one principal direction. The site's cross-flow component is reported in `graft_metrics.json` but not imposed.
4. The patch is quasi-steady Stokes flow, which is valid for this regime and checked with `--stage verify`.
5. The fixed-speed lid gives protruding textures a few per cent more shear. This biases the comparison *in favour of* protruding designs, which makes it conservative against the height-neutral hypothesis.
6. The platelet model uses simplified activation and aggregation, no red cells, and no detachment. Its rate constants were measured on collagen, so compare designs rather than quoting absolute counts.
7. Textures are idealized heightfields, so undercut crowns are not represented.

## 4. Reproduce

```bash
python3 evar/run_macro.py validation
```

```bash
python3 evar/run_macro.py graft --np 4
```

```bash
python3 evar/run_patch.py --stage solve
```

```bash
python3 evar/run_patch.py --stage post --site worst
```

```bash
python3 evar/run_patch.py --stage verify --only den_flush --site worst
```

```bash
python3 evar/run_patch.py --stage platelets --site worst --seeds 3
```

| Command | Time on 4 cores |
|---|---|
| `validation` | ~30 min |
| `graft` | ~1.5 h |
| `--stage solve` | ~1.5 h for all designs |
