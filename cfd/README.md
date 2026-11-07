# DentiFlow CFD pipeline

This pipeline is an automated, reproducible workflow that takes the stent CAD files in this repo through four steps:

1. **Clean, watertight stent geometry**: the bare stent, Sanhita's denticle sector, and parametric denticle variants.
2. **High-quality volume meshes**, built with OpenFOAM `snappyHexMesh` and checked with `checkMesh`.
3. **Blood-flow simulations**, run with **OpenFOAM** (open source, the standard in academic cardiovascular CFD):
   - non-Newtonian Carreau blood
   - parabolic inlet profile
   - steady flow, or a pulsatile LAD coronary waveform
4. **Hemodynamic metrics**: TAWSS, OSI, RRT, low-WSS area, recirculation volume and pressure drop. These are written to one CSV with comparison plots.

```
config/experiment.json ──► geometry.py ──► snappyHexMesh ──► simpleFoam / pimpleFoam ──► post.py ──► results/
 (design matrix)            (STL repair,     (mesh + quality    (Carreau blood,             (TAWSS, OSI, RRT,
                             denticles)       check)             LAD inflow)                 maps, CSV)
```

---

## 1. Install (one time)

**Python 3.10+** and the Python packages:

```bash
pip install -r requirements.txt
```

**OpenFOAM.** Use Docker. This works the same on Mac, Windows and Linux:

1. Install Docker Desktop (https://www.docker.com/products/docker-desktop/) and start it.
2. That's it. The first run downloads the official `opencfd/openfoam-default:2406` image (about 500 MB).

A native OpenFOAM install (openfoam.com **v2206 or newer**) also works. Pass `--backend native` and set `FOAM_BASHRC=/path/to/openfoam/etc/bashrc` if it isn't in the default location.

> ⚠️ Do **not** use Ubuntu's `apt install openfoam` package (v1912). Every function object crashes in that build ("error in IOstream sha1"), so no wall shear stress can be computed.

---

## 2. Run

All commands are run from this `cfd/` folder.

List the designs in the experiment:

```bash
python3 run.py --list
```

Build only the cleaned stent STL files (takes seconds). They are written to `geometry_out/`:

```bash
python3 run.py --geometry-only
```

1. Check that the solver is correct against the exact Poiseuille solution:

```bash
python3 run.py --mode steady --preset coarse --only validation_poiseuille
```

2. Quick screen of all designs: steady mean flow, coarse mesh. About 6–15 min per design on 4 cores:

```bash
python3 run.py --mode steady --preset coarse
```

3. Mesh-independence check on the control and one textured design:

```bash
python3 run.py --mode steady --preset medium --only bare,D1_baseline
```

```bash
python3 run.py --mode steady --preset fine --only bare,D1_baseline
```

4. Final experiment: pulsatile flow on the chosen mesh preset. This takes hours per design:

```bash
python3 run.py --mode pulsatile --preset medium --np 8
```

Useful options:

| Option | Meaning |
|---|---|
| `--only a,b` | Run only these designs |
| `--np N` | Number of CPU cores |
| `--skip-existing` | Resume a batch without redoing finished designs |
| `--cycles N` | Number of cardiac cycles (default 2; metrics are averaged over the last one) |

**Outputs**

| Path | Content |
|---|---|
| `geometry_out/<design>.stl` + `.json` | Cleaned, watertight stent sector (mm) and its stats (denticle count, protrusion height…) |
| `runs/<preset>-<mode>/<design>/` | The full OpenFOAM case: mesh, logs (`log.*`), fields. Open `case.foam` in ParaView to visualise. |
| `runs/.../<design>/mesh.json` | Cell count and `checkMesh` quality summary |
| `results/<preset>-<mode>/metrics.csv` | **One row per design**: every metric below |
| `results/<preset>-<mode>/comparison.png` | Bar charts of the key metrics vs. the bare control |
| `results/<preset>-<mode>/wss_maps_<design>.png` | Unrolled maps of (TA)WSS on the artery wall and stent, plus OSI for pulsatile runs |

`runs/` is in `.gitignore`. The cases are large and can be fully regenerated from the config.

---

## 3. What is simulated

### Geometry
- **Stent:** `stent fr.stl` (Ø3.2 mm, 100 µm struts, 16 rings). The repo's STL files are inside-out and made of 64 overlapping bodies. The pipeline flips and unions them into one valid solid.
- **Simulated region:** a **120° sector, 2 rings long** (x = 1.253–3.760 mm). This is exactly the window of Sanhita's `stent_denticle_CFD_sector.stl`, so every design is compared on identical strut geometry.
  - I checked her sector: its struts match `stent fr.stl` to 0.00 µm.
  - Cleaning her file removed one floating denticle fragment that had been cut off at the sector edge, plus two degenerate slivers.
- **Artery:** a rigid tube, R = 1.59 mm. The struts (r = 1.50–1.60 mm) sit 10 µm into the wall and protrude 90 µm into the lumen. There is 1 mm of artery upstream and 6 mm downstream so the inlet and outlet do not affect the stent region.
- **Sector sides:** symmetry planes.

### Denticle designs (the independent variable)

Each parametric denticle is a shield-shaped crown with longitudinal riblets. The defaults are measured from Sanhita's CAD:

| Property | Default |
|---|---|
| Footprint | 100 µm × 70 µm |
| Height | 43 µm |
| Spacing along the strut centreline | 85 µm |
| Riblets | 3 |

Denticles are placed on the luminal face of every strut and fused to it with an exact boolean union. To add a design, add an entry to `config/experiment.json`. All parameters are listed in `DenticleParams` in `pipeline/geometry.py`:

| Parameter | Meaning |
|---|---|
| `length_um`, `width_um`, `height_um` | Crown size |
| `spacing_um` | Centre-to-centre spacing along the strut |
| `n_ridges`, `ridge_depth` | Riblet count, and groove depth as a fraction of height (0 = smooth) |
| `orientation` | `strut` (along strut), `flow` (along vessel axis) or `reverse` (pointing upstream) |

The default experiment is **one-factor-at-a-time** around the baseline denticle:

| Design | What changes | Role |
|---|---|---|
| `validation_poiseuille` | No stent, Newtonian blood | Proves the solver gives the exact answer |
| `bare` | Untextured stent | **Control** |
| `student_sector` | Sanhita's own denticle CAD | Her design |
| `D1_baseline` | Parametric, matched to hers | Baseline for the sweep |
| `D2_tall` / `D3_short` | Height 65 / 25 µm | Effect of height |
| `D4_dense` / `D5_sparse` | Spacing 65 / 150 µm | Effect of density |
| `D6_smooth` | No riblets | Effect of riblets |
| `D7_flow_aligned` / `D8_reversed` | Axis-aligned / upstream-pointing | Effect of orientation |

### Blood and flow (held constant)

| | Value | Source |
|---|---|---|
| Density | 1060 kg/m³ | Sanhita's notes |
| Viscosity | Carreau: μ₀ = 0.056 Pa·s, μ∞ = 0.00345 Pa·s, λ = 3.313 s, n = 0.3568 | Cho & Kensey 1991 (standard blood values) |
| Inlet | Parabolic profile | Sanhita's notes |
| Mean / peak flow | 42.2 / 79.1 mL/min (Case A) | Sanhita's notes |
| Period | 0.9 s (100 × 0.009 s) | Sanhita's notes |
| Steady mode | Mean flow (set `"steady_flow": "peak"` for peak flow) | |
| Outlet | p = 0 | |
| Walls | Rigid, no-slip | |
| Flow regime | Laminar | |

> ⚠️ **Replace the waveform before final runs.** `config/waveform_LAD_approx.csv` is a *generic* coronary shape (low systolic flow, early-diastolic peak). Its mean and peak are scaled to Case A, but its shape is not the paper's. To fix it:
> 1. Digitise the LAD waveform from the source paper with WebPlotDigitizer.
> 2. Save it as a CSV with the columns `time_s,flow_mL_min` and point `inlet.waveform_csv` at it.
>
> The pipeline rescales it to the period, mean and peak in the config.

### Metrics (dependent variables)

All metrics are computed on the artery wall under the stent (`wall_stented_*`), on the wall in the 1 mm downstream of it (`wall_wake_*`), and on the stent surface (`stent_*`):

| Metric | Definition | Better when |
|---|---|---|
| TAWSS (Pa) | Time-averaged \|wall shear stress\| | Higher |
| Low-WSS area (%) | Area with TAWSS < 0.4 Pa (Malek et al., JAMA 1999) | Lower |
| OSI | ½(1 − \|mean τ\| / TAWSS): how much the shear direction flips (0 to 0.5) | Lower |
| High-OSI area (%) | Area with OSI > 0.1 | Lower |
| RRT (1/Pa) | 1 / ((1 − 2·OSI)·TAWSS): relative residence time near the wall | Lower |
| Recirculation volume (mm³) | Fluid with reversed time-mean axial velocity, from the stent start to 1 mm past its end | Lower |
| Pressure drop (Pa) | Inlet minus outlet | Lower |

Steady mode gives only WSS-based metrics, with OSI = 0. OSI and RRT need the pulsatile mode.

---

## 4. Mesh presets and cost

The cell size at the stent sets how many cells resolve a 43 µm denticle:

| Preset | Background cell | Cell at artery wall | Cell at stent | Cells across a 43 µm denticle | Typical cells |
|---|---|---|---|---|---|
| coarse | 80 µm | 40 µm | 10 µm | ~4 | 0.35 M (bare) |
| medium | 60 µm | 30 µm | 7.5 µm | ~6 | ~1 M |
| fine | 50 µm | 12.5 µm | 3.1 µm | ~14 | several M |

Measured on a 4-core cloud machine, coarse bare stent:
- snappyHexMesh: 83 s
- steady solve: 4 min (converged in 568 iterations)

Pulsatile runs need thousands of time steps per cycle, so expect hours per design on a laptop. A university or lab workstation is a better fit for the final runs.

**Mesh-independence:** run `bare` and `D1_baseline` at coarse, medium and fine. Use the coarsest preset whose key metrics change by less than ~5% at the next finer level. Report this in the project; judges will ask.

---

## 5. Assumptions and limitations (state these in the write-up)

1. **2-ring, 120° sector with symmetry sides.** The real stent is periodic every 120°, and periodic (cyclic) sides would be exact. Symmetry planes are a close, commonly used approximation.
2. **Rigid, straight, healthy artery.** There is no wall motion, curvature, plaque, or stent expansion and deformation.
3. **Approximate inlet waveform** until the digitised one is added (see above).
4. **Carreau model:** blood is treated as a single-phase fluid. There are no red cells or platelets, so platelet adhesion itself is not simulated, only the flow conditions linked to it.
5. **Hemodynamic thresholds** (0.4 Pa, OSI 0.1) are literature conventions for "atheroprone/thrombogenic" flow, not exact biological cut-offs.
6. **Denticles are on the luminal strut face only**, and parametric denticles are idealised shapes.

---

## 6. File map

```
cfd/
├── run.py                    # orchestrator (CLI)
├── config/
│   ├── experiment.json       # design matrix + physics + flow settings
│   └── waveform_LAD_approx.csv
├── pipeline/
│   ├── geometry.py           # STL repair, sector crop, denticle generator, lumen surface
│   ├── case.py               # OpenFOAM dictionaries (mesh, BCs, Carreau, solver, function objects)
│   ├── waveform.py           # inlet waveform loading/rescaling
│   └── post.py               # TAWSS/OSI/RRT, recirculation, maps, comparison plots
├── geometry_out/             # cleaned stent STL files (generated)
├── runs/                     # OpenFOAM cases (generated, git-ignored)
└── results/                  # metrics.csv + figures (generated)
```
