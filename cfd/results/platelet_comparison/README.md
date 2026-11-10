# Platelets in blood: denticle-coated stents vs bare stents with the same strut thickness

All runs were generated with `python3 compare_platelets.py` from `cfd/`, on 2026-11-10. Every design uses the same settings:
- the same 120° two-ring sector
- the same **coarse** mesh settings
- the same pulsatile LAD inflow: mean 42.2 mL/min, period 0.9 s
- the same shear-thinning Carreau blood
- the same platelet model: Sorensen 1999 kinetics, 2 heartbeats = 1.8 s, identical parameters and seed

The platelet model is described in `cfd/README.md` and is the same one used for `results/preliminary/3_platelets`. The only change is that platelets are now seeded from r = 1.40 mm, so that they also reach the 43 µm-tall denticles. Absolute counts are therefore not comparable with the preliminary run.

> ⚠️ **Preliminary.** This uses a coarse mesh, one random realisation per design, the generic inflow waveform, and a simplified aggregation rule. Read the "How much to trust this" section before quoting any number.

## Why "same strut thickness"?

The earlier steady screen found that every denticle design was worse than the bare stent, in proportion to denticle height. That result cannot separate "the denticle *shape* is bad" from "anything that sticks 40 µm further into the flow is bad". So each textured stent is compared against a bare stent whose struts are thickened inward by the denticle height:

| Textured design | Denticle height | Control | Control strut thickness |
|---|---|---|---|
| `student_sector`: Sanhita's CAD | 42.7 µm | `bare_thick143` | 143 µm |
| `D4_dense`: parametric, 65 µm spacing | 41.8 µm | `bare_thick143` | 143 µm |
| `D3_short`: parametric, 25 µm tall | 24.3 µm | `bare_thick124` | 124 µm |
| *(reference)* | n/a | `bare` | 100 µm |

The control is **stricter than the denticles**. It thickens the *whole* strut, adding 0.055 mm³ of metal at 143 µm, whereas the denticles add only 0.012–0.024 mm³ as separate bumps. So "better than its control" means *better than a uniformly thicker strut*. It does **not** mean better than the original 100 µm stent. Both comparisons are given below.

## Files

| Path | Content |
|---|---|
| `comparison.png` | Bar charts of every design: deposition, density, aggregation and WSS metrics |
| `comparison.csv` | One row per design with all metrics: platelets, flow, geometry, counting error, cluster statistics |
| `pairwise_textured_vs_control.csv` | Each textured design vs its control: values, % change, counting-noise z-scores |
| `pair_<textured>_vs_<control>.gif` | Side-by-side growth of the deposits, with cumulative curves (60 frames over 1.8 s) |
| `pair_bare_thick143_vs_bare.gif` | Effect of strut thickness alone: 143 µm vs 100 µm |
| `<design>/platelets.gif` | Free platelets (resting/activated), adhesion and aggregation, plus a zoom on the busiest spot |
| `<design>/deposition_growth.gif`, `deposition_map.png` | Where and when platelets stuck |
| `<design>/flow.gif`, `frame_peak_flow.png` | Inflow waveform, unrolled wall-WSS map, and near-wall slice with streamlines |
| `<design>/adhered_platelets.csv` | Every deposited parcel: time, x/y/z (mm), site (strut top / side / aggregate), activated flag |
| `<design>/platelet_timeseries.csv` | Free, adhered, aggregated and injected platelets vs time |
| `<design>/flow_timeseries.csv` | Per snapshot: inflow, pressure drop, wall/stent WSS, low-WSS area, recirculation |
| `<design>/surface_tawss_osi_rrt.csv` | TAWSS, OSI and RRT on every wall and stent face |
| `<design>/summary.json`, `flow_metrics.json` | Raw platelet and flow summaries |

---

## Results

### 1. Textured vs height-matched control

The numbers below are platelets after 1.8 s. "z" is the difference divided by its particle-counting noise; |z| > 2 is unlikely to be chance (see the caveats).

| | **Sanhita's denticles** vs 143 µm bare | **D4 dense** vs 143 µm bare | **D3 short** vs 124 µm bare |
|---|---|---|---|
| Directly adhered to struts | 32,900 vs 17,300: **+91%** (z = 12.3) | 13,100 vs 17,300: **−24%** (z = −4.2) | 11,300 vs 13,200: **−14%** (z = −2.1) |
| Aggregated onto deposits | 18,400 vs 2,800: **+551%** | 1,950 vs 2,800: −31% | 1,040 vs 1,920: −46% |
| — excluding the 2 largest clusters | 3,400 vs 2,500: +38% | 1,800 vs 2,500: −26% | 910 vs 1,690: −46% |
| Deposits per mm² of stent | 9,100 vs 3,750: **+143%** | 2,960 vs 3,750: −21% | 2,720 vs 3,140: −13% |
| Artery-wall TAWSS | 0.565 vs 0.519 Pa: +9% | 0.565 vs 0.519 Pa: +9% | 0.585 vs 0.552 Pa: +6% |
| Wall area with TAWSS < 0.4 Pa | 27.4 vs 32.9%: −17% | 27.5 vs 32.9%: −16% | 25.1 vs 29.2%: −14% |
| Stent TAWSS | 0.79 vs 0.92 Pa: **−15%** | 0.86 vs 0.92 Pa: −7% | 0.90 vs 0.92 Pa: −2% |
| Stent area with TAWSS < 0.4 Pa | 46.8 vs 43.2%: **+8%** | 39.8 vs 43.2%: −8% | 39.1 vs 42.0%: −7% |
| Recirculation volume | −41% | −54% | −44% |
| Pressure drop | −1.8% | −1.8% | −1.5% |

### 2. Every design vs the original 100 µm bare stent

| Design | Adhered | Aggregated | Deposits per mm² | Wall area < 0.4 Pa | Wall TAWSS |
|---|---|---|---|---|---|
| **bare, 100 µm** (reference) | 12,100 | 1,470 | 3,320 | 23.5% | 0.601 Pa |
| bare, 124 µm | 13,200 | 1,920 | 3,140 | 29.2% | 0.552 Pa |
| bare, 143 µm | 17,300 | 2,820 | 3,750 | 32.9% | 0.519 Pa |
| D3 short (24 µm) | 11,300 (−7%, z ≈ −1.0) | 1,040 | 2,720 | 25.1% | 0.585 Pa |
| D4 dense (42 µm) | 13,100 (+8%, z ≈ +1.1) | 1,950 | 2,960 | 27.5% | 0.565 Pa |
| **Sanhita's denticles (43 µm)** | **32,900 (2.7×, z ≈ +17)** | **18,400** | **9,100** | 27.4% | 0.565 Pa |

---

## What this means

**1. Height drives the flow metrics, and denticles are less bad than a uniformly thicker strut.**
- In all three pairs, the textured stent beats its height-matched control on every artery-wall metric: 14–17% less low-WSS area, 6–9% higher TAWSS, about half the recirculation.
- Separate bumps let near-wall flow pass between them, whereas a solid thicker strut blocks it.
- Against the real 100 µm stent, however, every denticle design is still worse on the wall: 25–27.5% low-WSS area vs 23.5%.
- So the earlier conclusion holds. **Adding denticles on top of a strut costs hemodynamics roughly in proportion to their height, just less than thickening the strut would.**

**2. The two parametric designs (D3, D4) do not change platelet adhesion meaningfully compared with the 100 µm stent.**
- Their adhered counts are within ±8% of bare, with |z| ≈ 1, which is within counting noise.
- Their per-area densities are 11–18% lower, because they add surface without adding proportionally more deposits.
- They beat their thick controls by 14–24%. That is real (z = −2.1 and −4.2), but it mostly reflects how bad thick struts are.
- **There is no evidence here that denticles *reduce* platelet adhesion below the plain stent.**

**3. Sanhita's own CAD is clearly worse: about 2.7× the adhesion of the 100 µm stent and 1.9× that of its height-matched control.**
- This is the one strong, statistically unambiguous result (z = 12–17).
- It is driven by the denticle shape, not its height: D4 has the same height (41.8 vs 42.7 µm) and a similar denticle surface area (1.95 vs 1.89 mm²), yet D4 is *better* than the control.
- Two things set her geometry apart:
  - **Undercuts.** 27% of her denticle surface faces back toward the artery wall, like a real shark denticle's crown on a narrower stalk. The parametric denticles are heightfields and have none.
  - **Steeper upstream faces.** 18% of the surface faces upstream, vs 9–14% for D3 and D4.
- Her design is also the only one whose *stent* surface shear stress drops (−15%), and whose low-WSS stent area grows rather than shrinks. So the blood moves more slowly over her denticles than over a flat strut top.
- 89% of her deposits sit on the denticle surfaces. Only ~1% sit in the shadow *under* a crown. The platelets stick to the exposed crowns, not in hidden pockets.
- Which shape feature causes the effect is **not pinned down by these runs**. A finer mesh, plus one or two parametric variants with undercuts, would settle it (see below).

**4. Aggregation is dominated by two runaway clusters on Sanhita's design, and those should not be over-read.**
- 81% of her aggregated platelets (≈15,000) sit in **two clusters**, both at x ≈ 1.55 mm on the upstream ring. They start at 0.38 s and 0.55 s and grow into towers ~45 µm tall.
- No other design has a cluster larger than ~200 platelets.
- Nucleating a tower is a rare, chance event, and the model has no way to stop one growing:
  - no thrombus-on-flow feedback
  - no shear erosion or embolisation
- So the "+551% aggregation" mostly measures *whether a tower happened to start*. Its z-score overstates the certainty, because clustered deposits are not independent counts.
- With those two clusters removed, aggregation on her design is +38% vs the control, which is still worse, in line with the adhesion result.
- **Use the direct-adhesion numbers as the main platelet metric.**

**5. Significance for the project's hypothesis.**
- The idea behind shark-skin riblets comes from drag reduction in **turbulent** flow, where riblets spaced at ~15 viscous wall units damp near-wall vortices.
- Coronary flow is **laminar**: Re ≈ 85 at mean flow and ≈ 160 at peak in this model, and the near-wall layer around struts is viscous-dominated.
- In laminar flow, extra roughness adds wetted area and slow-flow pockets, and there is no turbulence for it to suppress. These results match that expectation.
- So the evidence so far **does not support the hypothesis** that denticles reduce thrombogenicity on a coronary stent. That is a legitimate, explainable result for the science fair, not a failure.
- The most interesting result for the project is the shape sensitivity: two geometries of the same height behave in opposite directions against the same control.

---

## How much to trust this

| Issue | Effect | What would fix it |
|---|---|---|
| **Coarse mesh** | Stent cells are ~10 µm, so a 43 µm denticle is ~4 cells tall and D3 is ~2.5. Fine CAD features such as Sanhita's undercuts are under-resolved. This matters most for the student result. | Re-run at least `student_sector`, `D4_dense` and `bare_thick143` with `--preset medium`/`fine`, and check that adhesion changes by less than ~10%. |
| **One random realisation** | The z-scores only cover particle counting noise. They do not cover mesh or seed variation, or the clustered aggregation. | Run 3 seeds per design (`PlateletParams(seed=…)`) and report the mean ± SD. |
| **Aggregation rule** | It is simplified: activation within 10 µm of a deposit, and towers can grow without limit. | Add shear-dependent detachment or a limit on aggregate height. Until then, treat aggregation as qualitative. |
| **One-way coupling** | Deposits do not change the flow. That is fine over 1.8 s (≤ 13% of saturation coverage), but not for real thrombus growth over minutes. | Not needed for an initial-deposition comparison. |
| **Short window** | 1.8 s shows initial deposition only. | Run longer (`--platelet-cycles`), at the cost of ~45 min of CPU per cycle per design. |
| **Generic inflow waveform** | It affects absolute numbers, and probably not the ranking. | Replace it with the waveform digitised from the paper. |
| **Two-ring sector** | The upstream ring gets ~80% of deposits in every design, partly because it is the leading edge of the domain. | This affects all designs equally, so the comparisons are fair. Absolute per-ring numbers are not. |
| **Thick-strut controls** | They are stricter than the denticles (+0.055 mm³ of metal vs +0.012–0.024). | Also report vs the 100 µm stent, as done above. A volume-matched control is another option. |

**Suggested next runs, in order:**
1. A medium-mesh repeat of `student_sector` vs `D4_dense` vs `bare`, to check that the main result survives.
2. Three seeds for each of those.
3. A parametric denticle with an undercut (crown wider than its stalk), to test whether the undercut is what makes Sanhita's design worse.

## Reproduce

From `cfd/`:
- **Everything:** `python3 compare_platelets.py`. This takes about 8 h of CFD on 3 cores, then 50–90 min of platelets per design on 1 core. The two stages overlap if you run `--stage cfd` and `--stage platelets --wait` side by side.
- **Tables, chart and pair GIFs only:** `python3 compare_platelets.py --stage summary`, from existing runs, in ~5 min.
