# DentiFlow: Shark-Skin-Inspired Microtextured Coronary Stent

**Student:** Sanhita Tummala
**Sources:** `rp2027_.docx` (research plan), `2027_proj.docx` (project idea + abstract drafts + title brainstorm), `sf27.docx` (literature notes + CFD parameters), and the CAD files in this repo.

---

## 1. What the project is

Coronary stents keep narrowed arteries open, but they still fail through three main problems:

- **Stent thrombosis.** Platelets stick to the stent and form a clot.
- **In-stent restenosis.** The artery narrows again inside the stent.
- **Bacterial adhesion and biofilm** on the implanted device.

Shark skin is covered in microscopic **dermal denticles**. These are tooth-like scales that change how fluid flows near the surface. Shark-skin-inspired surfaces have been reported to reduce bacterial attachment, and they show some promise against platelet adhesion.

**The idea:** put denticle-inspired microtextures on a coronary stent surface and test whether they:

1. improve local blood flow (hemodynamics), and
2. reduce platelet adhesion and bacterial biofilm compared with an untextured stent.

**Research question (from `rp2027_.docx`):**
> Can shark-skin-inspired microtextured surfaces be engineered for cardiac stents to reduce platelet adhesion and bacterial biofilm formation compared with conventional, untextured surfaces?

**Hypothesis:**
> Shark-skin-inspired microtextured stent surfaces will show improved hemodynamics, reduced platelet adhesion, and reduced bacterial biofilm formation compared with untextured stent surfaces.

---

## 2. Engineering goals

1. Build multiple CAD models of denticle-textured stent surfaces.
2. Compare the geometries with CFD.
3. Fabricate physical prototypes or surface replicas of the best designs.
4. Characterize the prototype surfaces.
5. Measure platelet adhesion on each surface.
6. Measure bacterial attachment and biofilm on each surface.
7. Identify one geometry that performs well biologically **and** can be fabricated reproducibly.

---

## 3. Experimental design

### Phase 1: Background and computational foundation
- Literature review covering:
  - stent design
  - coronary blood flow and the math behind it (Poiseuille's law, non-Newtonian blood models)
  - CFD methods
  - platelet adhesion mechanisms
  - denticle biology and hydrodynamics
- Learn the CAD tools (OpenSCAD, FreeCAD, SolidWorks) and the CFD and meshing tools.

### Phase 2: CAD and CFD (computational experiment)

**Designs.** The plan is a stepwise series of designs:

| # | Design | Purpose |
|---|---|---|
| 1 | Plain stent (no texture) | **Control / baseline** |
| 2 | Stent with basic denticles | First textured test case |
| 3+ | Denticle variants | Optimization |

The notes list several ideas for the variants:
- denticle spacing and geometry
- continuous vs. segmented riblets
- struts with holes vs. without (based on the "punched strut" paper)
- a version with a graft cover

**Stent geometry (from the notes and the CAD files):**
- about 3.2 mm diameter
- about 20 mm length, which matches the reported mean stent length of 20.0 ± 7.9 mm
- 16 stacked sinusoidal rings with connector links

**Planned flow conditions (adapted from a published LAD coronary CFD study):**

| Parameter | Value |
|---|---|
| Flow | Pulsatile, using a human LAD (left anterior descending artery) waveform |
| Inlet | Parabolic velocity profile |
| Mean / peak flow, Case A | 42.2 / 79.1 mL/min |
| Mean / peak flow, Case B | 45.1 / 84.7 mL/min |
| Blood model | Carreau non-Newtonian, density 1060 kg/m³ |
| Walls | Rigid artery, no-slip on artery and stent surfaces |
| Flow regime | Laminar. Peak Re = 195 (Case A) / 260 (Case B). Womersley number = 1.9 / 1.4 |
| Time stepping | Δt = 0.009 s, 100 steps per cardiac cycle, 1 cycle simulated |

The notes also record a simpler Newtonian setup from the literature as an alternative: viscosity μ = 3.5 × 10⁻³ Pa·s.

**Variables**
- **Independent:** surface texture geometry (none vs. the denticle variants).
- **Dependent:** the hemodynamic outputs:
  - flow velocity
  - recirculation area
  - time-averaged wall shear stress (TAWSS)
  - oscillatory shear index (OSI)
  - relative residence time (RRT)
- **Constants:** vessel geometry, inlet waveform, blood properties, mesh settings, boundary conditions.

**Success criterion.** The notes say the goal is to avoid three things that promote thrombosis:
- low wall shear stress
- flow recirculation
- long residence time

So the best design has the highest TAWSS and the lowest OSI, RRT and recirculation area.

**Output:** a short list of top-performing textures to carry into fabrication.

### Phase 3: Fabrication and surface replication
- 3D print the selected designs as **negative molds**.
- Imprint the textures into **PLGA**. PLGA (poly(lactic-co-glycolic acid)) is a biodegradable, FDA-approved copolymer that is already used in drug-eluting stent coatings. Its degradation rate can be tuned through the lactic-to-glycolic ratio.
- Characterize the replicated surfaces to confirm the texture actually transferred.
- Make a flat (untextured) PLGA control.

### Phase 4: Platelet adhesion
- Expose the textured PLGA and the flat control to human platelets or blood.
- Quantify the adhered platelets, for example with a microplate reader.
- **Independent variable:** texture. **Dependent variable:** platelet adhesion. **Control:** flat PLGA.

### Phase 5: Bacterial adhesion and biofilm
- Culture bacteria on the textured and control surfaces.
- Quantify biofilm with **crystal violet staining** and absorbance on the microplate reader.

### Phase 6 (from the `2027_proj` abstract only): Endothelial cell growth
- Measure endothelial cell regrowth on the surfaces. Good regrowth is what prevents late thrombosis and restenosis.
- This phase is in the abstract drafts but not in the research plan. See issue 3 below.

### Analysis
Compare every textured design against the untextured control in each phase, using Excel, Stat-Ease or similar software.

---

## 4. Existing CAD files in this repo

All dimensions are in mm.

| File | What it appears to be |
|---|---|
| `stent fr.stl` | Full stent, 20.1 mm long, Ø3.2 mm (control geometry) |
| `graft_cover.stl` | Thin cylindrical graft sleeve, Ø3.36 mm |
| `stent+cover.stl` | Stent with the graft cover merged on |
| `plain stent with connectors.stl` | A single ring unit, 1.33 mm long |
| `shark-denticle-diamond-layout.stl` | Patch of denticles in a diamond layout |
| `cad models/.../stent_denticle_CFD_sector.zip` | Sector of the denticle-textured stent, cut down for CFD (~393k triangles) |

---

## 5. Issues found in the documents (fix before submitting)

1. **`sf27.docx` contains a different project.** Its research question, expected outcomes, procedures 4–5, risk and safety section, and Experiments 1–3 are all about *Blakeslea trispora* β-carotene production on cassava. That looks like a leftover template from an earlier project. None of it belongs in this plan, and it must not appear in anything submitted.
2. **The target vessel is inconsistent.** Every plan document and CAD file is for a **coronary** stent (Ø3.2 mm). But the title brainstorm in `2027_proj.docx` is almost entirely about **thoracic aortic stent grafts**. The thoracic aorta is roughly 25–35 mm across, with a completely different flow regime. Pick one. The work done so far supports coronary.
3. **"Stent" and "stent graft" are used interchangeably.** The hypothesis and goals in `rp2027_.docx` say "stent graft", and a graft cover model exists, but the rest of the plan describes a bare-metal-style stent. Decide whether the graft cover is part of the design or just one of the variants.
4. **Endothelial cell growth** is one of four evaluation criteria in the `2027_proj` abstract, but it is missing from the research plan's goals and procedures.
5. **The research plan has no risk and safety section.** This matters: human blood or platelets and bacterial cultures both need prior review and approval under ISEF-style fair rules. They are typically done at a supervised research institution. Check your fair's specific rules before any wet-lab work.
6. **Biological methods are unspecified:**
   - which bacterial species
   - where the blood or platelets come from
   - which platelet quantification assay (for example LDH, microscopy or SEM)
   - sample sizes and replicates
   - which statistical test
7. **Fabrication resolution is unaddressed.** Real shark denticles are on the order of 100–200 µm, and stent struts are about 80–100 µm thick. A typical FDM 3D printer cannot make molds at that scale. Either:
   - test the textures on **scaled-up flat coupons**, and say so explicitly, or
   - use high-resolution resin (SLA/DLP) printing or soft lithography.
8. **Stent expansion is not modeled.** The notes flag this ("Did not account for the expansion"). It is acceptable as a stated limitation. Also list the rigid-wall assumption as a limitation.
9. **The physical tests use flat PLGA surfaces, but the CFD uses a full stent in an artery.** That is fine, but state clearly that the bench tests isolate surface texture effects, not whole-stent performance.
