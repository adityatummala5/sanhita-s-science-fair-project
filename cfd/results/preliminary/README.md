# Preliminary check: bare stent, 3 simulations

These three simulations were generated with `python3 preliminary.py` on 2026-11-07. All three use:
- the bare-stent sector (120°, 2 rings)
- the **coarse** mesh (~349k cells)
- the pulsatile LAD inflow: mean 42.2 mL/min, peak 79.1 mL/min, period 0.9 s
- 2 heartbeats, with metrics taken over the 2nd

> ⚠️ **These are preliminary.** The mesh is coarse, and the inflow waveform is the generic coronary shape, not the one digitised from the paper. See `cfd/README.md` for all limitations.

| Folder | Simulation | GIF | Data |
|---|---|---|---|
| `1_newtonian/` | Newtonian fluid with blood's density (1060 kg/m³) and high-shear viscosity (3.5 mPa·s) | `flow.gif` | `timeseries.csv`, `metrics.json`, `surface_tawss_osi_rrt.csv` |
| `2_carreau/` | Shear-thinning Carreau blood (Cho & Kensey 1991: μ₀ = 56 mPa·s, μ∞ = 3.45 mPa·s, λ = 3.313 s, n = 0.3568) | `flow.gif` | same as above |
| `3_platelets/` | Simulation 2 plus Lagrangian platelets: transport, adhesion, activation, aggregation | `platelets.gif`, `deposition_growth.gif` | `adhered_platelets.csv`, `timeseries.csv`, `summary.json` |
| *(top level)* | Newtonian vs Carreau | `comparison_newtonian_vs_carreau.png` | `comparison_newtonian_vs_carreau.csv` |

**What the flow GIFs show** (one frame per 1/20 of the heartbeat):
- the inflow waveform with the current moment marked
- wall shear stress (WSS) on the artery wall, unrolled flat
- near-wall velocity and streamlines in a slice through two struts

**Data files:**

| File | Content |
|---|---|
| `timeseries.csv` | Every snapshot: inflow, pressure drop, mean wall/stent WSS, low-WSS area, recirculation volume |
| `surface_tawss_osi_rrt.csv` | TAWSS, OSI and RRT on every wall and stent face (x, y, z in mm; face area in µm²) |
| `metrics.json` | Cycle-averaged summary, same definitions as `cfd/README.md` |

---

## 1 vs 2: Newtonian vs Carreau

Cycle-averaged values, on the artery wall under the stent unless noted:

| Metric | Newtonian | Carreau | Carreau vs Newtonian |
|---|---|---|---|
| TAWSS (Pa) | 0.49 | 0.60 | **+22%** |
| Wall area with TAWSS < 0.4 Pa | 31.1% | 23.5% | **−24%** |
| RRT (1/Pa) | 6.11 | 3.93 | −36% |
| Stent TAWSS (Pa) | 0.79 | 0.93 | +17% |
| Stent RRT (1/Pa) | 17.8 | 9.6 | −46% |
| Recirculation volume (mm³) | 0.0034 | 0.0027 | −19% |
| Pressure drop (Pa) | 9.5 | 11.4 | +20% |
| OSI | ≈ 0.0004 | ≈ 0.0003 | Both ≈ 0 |

**Interpretation:**
- **The viscosity model matters, and it matters most where it counts.** Next to the struts the shear rate is low, so the shear-thinning blood is *more* viscous there than the Newtonian fluid. That raises wall shear stress in exactly the zones that are flagged as thrombogenic.
- **A Newtonian model overstates the "bad" low-WSS area by ~30%** (31.1% vs 23.5%). The real comparison should use Carreau, which is now the default.
- **OSI ≈ 0 in both runs.** This waveform never reverses direction, and the Womersley number is low. The slow-flow pockets beside the struts therefore persist all cycle instead of oscillating. OSI will only become informative with a waveform that has flow reversal, or with stronger recirculation.
- **The danger moment is early systole** (t ≈ 0.2 s, when inflow is lowest). At that point 87% (Carreau) to 100% (Newtonian) of the wall under the stent is below 0.4 Pa.
- **Meshes:** the two runs used meshes differing by 0.06% in cell count. snappyHexMesh ran on 3 cores for the Newtonian case and 4 for the Carreau one, which is negligible.

## 3: Platelet adhesion and aggregation

The model is described in `cfd/README.md` and uses Sorensen et al. 1999 kinetics. Totals after 1.8 s (2 heartbeats):

| Quantity | Value |
|---|---|
| Platelets carried past the stent in the near-wall layer | 3.2 million |
| Adhered directly to struts | **13,500** |
| — on strut sides and corners | 8,750 (65%) |
| — on strut tops facing the flow | 4,770 (35%) |
| Aggregated onto already-deposited platelets | **2,670** (16% of all deposits, rising with time) |
| Mean deposit density | 3,950 platelets/mm², about 6% of the saturation coverage P∞ = 70,000/mm² |
| Upstream ring vs downstream ring | 13,280 vs 2,900 |

**Interpretation:**
- **Deposition is concentrated on the strut sides and corners.** These are the slow-flow, low-WSS pockets where struts meet the wall, and the same areas flagged dark blue in the WSS maps. This is the link between the hemodynamic metrics and the platelet result.
- **Deposition accelerates and turns into aggregation.** The `timeseries.csv` curve bends upward: the second heartbeat deposits ~4× more than the first. Deposited platelets activate passers-by, which then aggregate. This is the start of thrombus growth. A real thrombus takes minutes, so 1.8 s only shows the initial pattern.
- **The upstream ring gets ~80% of deposits. This is partly an artefact.** In this 2-ring section the first ring is the leading edge, and it sees fresh, undepleted near-wall blood. The densest spot (zoom panel) is the crown cut off at the inflow edge of the sector. In a full-length stent, the rings would share the load more evenly.
- **Each marker is a parcel of 21.5 platelets**, so local clustering is somewhat exaggerated. The totals are in real platelets.

## Before quoting any of these numbers
1. Replace the inflow waveform with the digitised one from the paper.
2. Run at the medium or fine mesh and check that the numbers change by less than 5%.
3. For the platelet model:
   - check the rate constants against Sorensen 1999, Table 1
   - state that the activation rule is a simplification
   - run a longer window or several seeds, so the stochastic counts have error bars

To reproduce: from `cfd/`, run `python3 preliminary.py`. This takes about 70 min of CFD per fluid model and about 60 min for the platelets on 4 cores. Re-exporting the GIFs and data alone uses `--skip-cfd --skip-platelets` and takes ~4 min.
