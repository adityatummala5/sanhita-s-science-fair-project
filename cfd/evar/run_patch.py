#!/usr/bin/env python3
"""Textured wall-patch CFD for the EVAR graft surface.

    python3 evar/run_patch.py --list
    python3 evar/run_patch.py --stage solve                 # mesh + steady solves at every shear level, all designs
    python3 evar/run_patch.py --stage post --site worst     # rebuild the cardiac cycle with a site waveform, metrics
    python3 evar/run_patch.py --only smooth,den_flush       # subset

Designs and shear levels are in config/evar.json ("patch"). Cases go to
cfd/runs/evar/patch/<design>/, results to cfd/results/evar/patch/.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import glob
import json
import math
import os
import re
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = os.path.dirname(HERE)
sys.path.insert(0, CFD)
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402

import run as R  # noqa: E402
import patch as P  # noqa: E402
import textures as TX  # noqa: E402
sys.path.insert(0, os.path.join(CFD, "pipeline"))
import case as fc  # noqa: E402

RUNS = os.path.join(CFD, "runs", "evar", "patch")
RES = os.path.join(CFD, "results", "evar", "patch")


def cfg_all():
    return json.load(open(os.path.join(CFD, "config", "evar.json")))


def build_mesh(design: TX.Texture, pm: P.PatchMesh) -> dict:
    base = os.path.join(RUNS, design.name, "mesh")
    if os.path.exists(base):
        shutil.rmtree(base)
    os.makedirs(base)
    Lx, Lz = design.domain()
    if design.mesh_dx_um or design.mesh_ny:
        pm = P.PatchMesh(**{**pm.__dict__, **({"dx_um": design.mesh_dx_um} if design.mesh_dx_um else {}),
                            **({"ny": design.mesh_ny} if design.mesh_ny else {})})
    nx, nz = P.write_block(base, Lx, Lz, pm)
    P.write_control(base, 1, ascii_points=True)
    P.write_common(base, fc.Physics())
    R.sh(["blockMesh"], base, "blockMesh")
    P.map_points(base, design, pm)
    json.dump(pm.__dict__, open(os.path.join(base, "patchmesh.json"), "w"))
    out = R.sh(["checkMesh", "-constant"], base, "checkMesh")
    info = {"cells": int(re.search(r"cells:\s+(\d+)", out).group(1)), "nx": nx, "nz": nz, "ny": pm.ny,
            "mesh_ok": "Mesh OK" in out,
            "failed_checks": int(m.group(1)) if (m := re.search(r"Failed (\d+) mesh checks", out)) else 0,
            "max_non_ortho": float(m.group(1)) if (m := re.search(r"Mesh non-orthogonality Max: ([\d.eE+-]+)", out)) else None,
            "max_skewness": float(m.group(1)) if (m := re.search(r"Max skewness = ([\d.eE+-]+)", out)) else None,
            **design.summary()}
    json.dump(info, open(os.path.join(RUNS, design.name, "mesh.json"), "w"), indent=2)
    return info


def solve_level(design: TX.Texture, k: int, gamma: float, pm: P.PatchMesh, ph, n_iter: int,
                init_from: str | None = None) -> str:
    cdir = os.path.join(RUNS, design.name, f"L{k}")
    if os.path.exists(cdir):
        shutil.rmtree(cdir)
    os.makedirs(os.path.join(cdir, "constant"))
    shutil.copytree(os.path.join(RUNS, design.name, "mesh", "constant", "polyMesh"),
                    os.path.join(cdir, "constant", "polyMesh"))
    P.write_common(cdir, ph)
    u_lid = gamma * pm.h_lid_mm * 1e-3
    eta = np.load(os.path.join(RUNS, design.name, "mesh", "eta_cells.npy"))
    if init_from:
        last = sorted(glob.glob(os.path.join(init_from, "[1-9]*")), key=lambda q: float(os.path.basename(q)))[-1]
        u_prev = json.load(open(os.path.join(init_from, "level.json")))["u_lid_m_s"]
        P.write_fields(cdir, u_lid, internal=P.read_internal_U(os.path.join(last, "U")) * (u_lid / u_prev))
    else:
        P.write_fields(cdir, u_lid, eta)
    P.write_control(cdir, n_iter)
    R.sh(["simpleFoam"], cdir, "simpleFoam")
    R.sh(["simpleFoam", "-postProcess", "-func", "wallShearStress", "-latestTime"], cdir, "wallShearStress")
    R.sh(["simpleFoam", "-postProcess", "-func", "grad(U)", "-latestTime"], cdir, "gradU")
    conv = None
    f = glob.glob(os.path.join(cdir, "postProcessing/wallStress/*/surfaceFieldValue.dat"))
    if f:
        v = [float(ln.split()[1].strip("(")) for ln in open(f[0]) if ln.strip() and not ln.startswith("#")]
        if len(v) >= 5:
            conv = abs(v[-1] - v[-5]) / max(abs(v[-1]), 1e-30)       # relative change over the last 100 iterations
    json.dump({"gamma_nominal_1_s": gamma, "u_lid_m_s": u_lid, "init_from": init_from,
               "wall_stress_rel_change_last100": conv}, open(os.path.join(cdir, "level.json"), "w"))
    return cdir


def extract_level(cdir: str, rho: float) -> dict:
    import pyvista as pv
    import logging
    foam = os.path.join(cdir, "case.foam")
    open(foam, "a").close()
    r = pv.OpenFOAMReader(foam)
    r.enable_all_patch_arrays()
    r.cell_to_point_creation = False
    r.set_active_time_value(r.time_values[-1])
    d = r.read()
    vol = d["internalMesh"]
    b = d["boundary"]
    gr = b[[k for k in b.keys() if k.split("/")[-1] == "graft"][0]]
    logging.disable(logging.CRITICAL)
    sized = gr.compute_cell_sizes(length=False, volume=False)
    logging.disable(logging.NOTSET)
    # local shear rate gamma = sqrt(2 S:S) from OpenFOAM's cell gradient grad(U)
    Gc = np.asarray(vol.cell_data["grad(U)"]).reshape(-1, 3, 3)
    Sc = 0.5 * (Gc + Gc.transpose(0, 2, 1))
    gam = np.sqrt(2 * np.einsum("nij,nij->n", Sc, Sc))
    out = {"face_c": np.asarray(gr.cell_centers().points), "face_A": np.asarray(sized.cell_data["Area"]),
           "gamma": gam,
           "tau": rho * np.asarray(gr.cell_data["wallShearStress"]),
           "cell_c": np.asarray(vol.cell_data["C"]) if "C" in vol.cell_data else np.asarray(vol.cell_centers().points),
           "cell_V": np.asarray(vol.cell_data["V"]), "U": np.asarray(vol.cell_data["U"]),
           "u_lid": json.load(open(os.path.join(cdir, "level.json")))["u_lid_m_s"],
           "iters": int(float(r.time_values[-1]))}
    np.savez_compressed(os.path.join(cdir, "fields.npz"), **out)
    return out


def regrad(d, rho):
    """Add grad(U) to already-solved levels and re-extract (for runs made before
    the gradient step was added)."""
    for c in sorted(glob.glob(os.path.join(RUNS, d.name, "L*"))):
        R.sh(["simpleFoam", "-postProcess", "-func", "grad(U)", "-latestTime"], c, "gradU")
        extract_level(c, rho)


def solve_design(d, pm, ph, levels, pc):
    """Mid level from a linear profile, then march outward, each level starting
    from the scaled neighbouring solution."""
    t0 = time.time()
    info = build_mesh(d, pm)
    mid = len(levels) // 2
    order = [mid] + list(range(mid + 1, len(levels))) + list(range(mid - 1, -1, -1))
    done = {}
    for k in order:
        nb = k - 1 if k > mid else k + 1
        init = done.get(nb) if k != mid else None
        done[k] = solve_level(d, k, levels[k], pm, ph, pc["n_iter"] if k == mid else pc["n_iter_from_neighbour"],
                              init_from=init)
        extract_level(done[k], ph.rho)
    print(f"=== {d.name}: {info['cells']:,} cells, checkMesh {'OK' if info['mesh_ok'] else info['failed_checks']} "
          f"(non-ortho {info['max_non_ortho']:.0f}, skew {info['max_skewness']:.2f}); "
          f"{len(levels)} levels in {(time.time() - t0) / 60:.1f} min", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", default="solve", choices=["solve", "post", "platelets", "verify", "all"])
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--duration", type=float, default=None, help="platelet simulated time (s)")
    ap.add_argument("--only", default="")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--site", default="worst", help="worst | median | analytical_limb | path/to/site.csv")
    ap.add_argument("--layer-ref", default="nominal", choices=["nominal", "crest"],
                    help="platelet reservoir height measured from the nominal surface or from the texture crest")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    args = ap.parse_args()

    cfg = cfg_all()
    pc = cfg["patch"]
    designs = [TX.from_dict(d) for d in pc["designs"]]
    extra = [TX.from_dict(d) for d in pc.get("mesh_check", [])]
    if args.list:
        for d in designs:
            print(f"{d.name:<22s} {d.kind:<10s} {d.placement:<9s} {d.description}")
        return
    if args.only:
        want = set(args.only.split(","))
        designs = [d for d in designs + extra if d.name in want]
    pm = P.PatchMesh(**pc.get("mesh", {}))
    ph = fc.Physics(**cfg["physics"])
    levels = pc["gamma_levels_1_s"]

    if args.stage in ("solve", "all"):
        R.setup_backend("auto")
        with cf.ThreadPoolExecutor(args.jobs) as ex:
            list(ex.map(lambda d: solve_design(d, pm, ph, levels, pc), designs))


    if args.stage in ("post", "all"):
        import post_patch
        post_patch.run(designs, cfg, pm, ph, args.site, RUNS, RES)

    if args.stage in ("platelets", "all"):
        run_platelets(designs, cfg, pm, ph, args)

    if args.stage == "verify":
        R.setup_backend("auto")
        for d in designs:
            verify_quasi_steady(d, cfg, pm, ph, args.site)


def verify_quasi_steady(d, cfg, pm, ph, site, n_cycles=2):
    """Full time-dependent patch run vs the quasi-steady reconstruction."""
    import post_patch
    t, tau_site, src = post_patch.site_waveform(site, cfg, ph, os.path.dirname(RES))
    U_t = np.sign(tau_site) * pm.h_lid_mm * 1e-3 * P.gamma_from_tau(tau_site, ph)
    T = cfg["inflow"]["period_s"]
    cdir = os.path.join(RUNS, d.name, "pulsatile")
    if os.path.exists(cdir):
        shutil.rmtree(cdir)
    os.makedirs(os.path.join(cdir, "constant"))
    shutil.copytree(os.path.join(RUNS, d.name, "mesh", "constant", "polyMesh"), os.path.join(cdir, "constant", "polyMesh"))
    P.write_common(cdir, ph)
    eta = np.load(os.path.join(RUNS, d.name, "mesh", "eta_cells.npy"))
    P.write_pulsatile(cdir, t, U_t, T, n_cycles, eta)
    t0 = time.time()
    R.sh(["pimpleFoam"], cdir, "pimpleFoam")
    import pyvista as pv
    foam = os.path.join(cdir, "case.foam")
    open(foam, "a").close()
    r = pv.OpenFOAMReader(foam)
    r.enable_all_patch_arrays()
    r.cell_to_point_creation = False
    r.set_active_time_value(r.time_values[-1])
    b = r.read()["boundary"]
    gr = b[[k for k in b.keys() if k.split("/")[-1] == "graft"][0]]
    taw_p = ph.rho * np.asarray(gr.cell_data["magWSSMean"])
    mv_p = ph.rho * np.linalg.norm(np.asarray(gr.cell_data["wallShearStressMean"]), axis=1)
    osi_p = 0.5 * (1 - mv_p / np.maximum(taw_p, 1e-12))
    lv = post_patch.load_levels(os.path.join(RUNS, d.name))
    sm = post_patch.load_levels(os.path.join(RUNS, "smooth"))
    _, f, _ = post_patch.design_metrics(d.name, d, lv, sm, t, tau_site, ph, pm)
    A = lv[0]["face_A"]
    out = {"design": d.name, "site": src, "runtime_min": round((time.time() - t0) / 60, 1),
           "TAWSS_pulsatile_Pa": float(np.sum(taw_p * A) / A.sum()),
           "TAWSS_quasisteady_Pa": float(np.sum(f["tawss"] * A) / A.sum()),
           "TAWSS_facewise_rms_diff_pct": float(100 * np.sqrt(np.sum((taw_p - f["tawss"]) ** 2 * A) / A.sum())
                                                / (np.sum(taw_p * A) / A.sum())),
           "OSI_pulsatile": float(np.sum(osi_p * A) / A.sum()), "OSI_quasisteady": float(np.sum(f["osi"] * A) / A.sum())}
    out["TAWSS_diff_pct"] = 100 * (out["TAWSS_quasisteady_Pa"] / out["TAWSS_pulsatile_Pa"] - 1)
    os.makedirs(RES, exist_ok=True)
    json.dump(out, open(os.path.join(RES, f"verify_quasisteady_{d.name}.json"), "w"), indent=2)
    print(json.dumps(out, indent=2))


def _platelet_job(a):
    import post_patch
    import patch_platelets as PP
    dd, seed, U_t, t, src, base, prm = a
    d = TX.from_dict(dd)
    done = os.path.join(base, d.name, f"seed{seed}", "summary.json")
    if os.path.exists(done) and json.load(open(done))["duration_s"] == prm.get("duration_s", 30.0):
        return json.load(open(done))                 # resume: this run already finished
    lv = post_patch.load_levels(os.path.join(RUNS, d.name))
    p = PP.PatchPlateletParams(**{**prm, "seed": seed})
    return PP.simulate(d, lv, U_t, t, os.path.join(base, d.name, f"seed{seed}"), p, label=src)


def run_platelets(designs, cfg, pm, ph, args):
    import post_patch
    import patch_platelets as PP
    t, tau_site, src = post_patch.site_waveform(args.site, cfg, ph, os.path.dirname(RES))
    U_t = np.sign(tau_site) * pm.h_lid_mm * 1e-3 * P.gamma_from_tau(tau_site, ph)
    site_tag = f"site-{os.path.basename(args.site).replace('.csv', '')}"
    base = os.path.join(RES, site_tag, "platelets" if args.layer_ref == "nominal" else f"platelets_{args.layer_ref}ref")
    pp_cfg = {k: v for k, v in cfg["patch"].get("platelets", {}).items() if not k.startswith("_")}
    pp_cfg["layer_ref"] = args.layer_ref
    rows = []

    jobs = [(d, s) for d in designs for s in range(1, args.seeds + 1)]
    extra = {"duration_s": args.duration} if args.duration else {}
    with cf.ProcessPoolExecutor(args.jobs) as ex:
        res = list(ex.map(_platelet_job, [(d.__dict__, s, U_t, t, src, base, {**pp_cfg, **extra}) for d, s in jobs]))
    by = {}
    for (d, s), r in zip(jobs, res):
        by.setdefault(d.name, []).append(r)
    ref = None
    for d in designs:
        rs = by[d.name]
        v = np.array([r["deposited_per_mm2_planform"] for r in rs])
        a = np.array([r["adhered"] / r["planform_mm2"] for r in rs])
        g = np.array([r["aggregated"] / r["planform_mm2"] for r in rs])
        rl = np.array([r["rate_late_per_mm2_s"] for r in rs])
        row = {"design": d.name, "seeds": len(rs),
               "rate_late_per_mm2_s_mean": round(float(rl.mean()), 2),
               "rate_late_per_mm2_s_sd": round(float(rl.std(ddof=1)) if len(rl) > 1 else 0.0, 2),
               "deposited_per_mm2_mean": round(float(v.mean()), 1),
               "deposited_per_mm2_sd": round(float(v.std(ddof=1)) if len(v) > 1 else 0.0, 1),
               "adhered_per_mm2_mean": round(float(a.mean()), 1), "aggregated_per_mm2_mean": round(float(g.mean()), 1),
               "adhered_per_mm2_wetted_mean": round(float(np.mean([r["adhered_per_mm2_wetted"] for r in rs])), 1),
               "pct_in_lower_third": rs[0]["pct_in_lower_third"], "duration_s": rs[0]["duration_s"]}
        if d.name == "smooth":
            ref = row
        rows.append(row)
    if ref:
        for row in rows:
            row["deposited_vs_smooth_pct"] = round(100 * (row["deposited_per_mm2_mean"] / ref["deposited_per_mm2_mean"] - 1), 1)
            row["rate_late_vs_smooth_pct"] = round(100 * (row["rate_late_per_mm2_s_mean"] / ref["rate_late_per_mm2_s_mean"] - 1), 1)
            se = math.hypot(row["rate_late_per_mm2_s_sd"], ref["rate_late_per_mm2_s_sd"]) / math.sqrt(row["seeds"])
            row["z_rate_vs_smooth"] = round((row["rate_late_per_mm2_s_mean"] - ref["rate_late_per_mm2_s_mean"]) / se, 2) if se else None
    import csv as _csv
    with open(os.path.join(base, "platelet_summary.csv"), "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=list(rows[-1].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
