# Literature

Reference material for the DentiFlow project. The papers cited in the proposal are listed in `proposal/proposal.html`.

| File | What it is | Relevance |
|---|---|---|
| `Coronary_Stent_CFD_Journal_Manuscript.docx` | Review manuscript, *CFD metrics for coronary stent development*. No author or journal is listed, so it is unpublished or a draft. About 2,300 words, 16 references. | Defines the hemodynamic endpoints we use (TAWSS, OSI, RRT, WSSG, recirculation, platelet shear exposure, pressure drop) and gives screening benchmarks. Also recommends a reporting and verification framework. |

## Notes on `Coronary_Stent_CFD_Journal_Manuscript.docx`

**Cite the primary papers it references, not the manuscript itself.** It has no authors or venue, so it can't be cited as a source.

Its references worth pulling:
- Beier et al. 2016, *Ann Biomed Eng* 44:315: stent design and hemodynamics
- Chesnutt & Han 2016, *Phys Biol* 13:016001: platelet simulation and malapposition
- McCarthy et al. 2025, *Sci Rep* 15:15868: sensitivity to rheology and outlet boundary conditions
- Himburg et al. 2004: definition of RRT
- He & Ku 1996: definition of OSI

**Its recommendations, compared with what our pipeline does now:**

| Item | Manuscript recommends | Our pipeline / proposal | Action |
|---|---|---|---|
| Low-TAWSS threshold | < 0.4–0.5 Pa, plus < 50% of the native value | < 0.4 Pa | Consistent. Could also report < 50% of the smooth control. |
| High-OSI threshold | Report area > 0.25, and ideally none > 0.45 | Area > 0.1 | Add > 0.25 and > 0.45 areas. This matters more for EVAR, where flow reverses. |
| WSSG (spatial gradient of WSS) | Report the 95th percentile of \|WSSG\| | Not computed | Add to `post.py`. Denticle crowns create steep gradients. |
| Statistics reported | Area fractions, median, IQR, 95th percentile, not just averages | Mostly area-averages and area fractions | Add percentiles. |
| Mesh convergence | < 5% change between the two finest meshes | < 5% (proposal) | Consistent. |
| Sensitivity analysis | Rheology, outlet BCs, inlet profile; device differences must exceed modeling uncertainty | Carreau vs Newtonian only | Add an outlet-BC check in the macro EVAR model. |
| Platelet shear exposure | Trajectory-integrated stress × time | Not computed (adhesion only) | Optional: log stress history along platelet paths. |

**Scope caveat.** The manuscript is about coronary stents: its benchmarks assume coronary flow (an unstented reference of about 2.7 Pa). The metric definitions and the reporting framework carry over to the EVAR work. The absolute benchmarks do not, because the resting shear in EVAR graft limbs is roughly an order of magnitude lower.
