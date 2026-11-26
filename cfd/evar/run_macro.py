#!/usr/bin/env python3
"""Vessel-scale EVAR CFD: Womersley validation and the bifurcated endograft.

    python3 evar/run_macro.py validation            # straight 25 mm tube, Newtonian, vs analytical Womersley
    python3 evar/run_macro.py graft                 # idealized bifurcated endograft, Carreau blood
    python3 evar/run_macro.py graft --preset medium # mesh-independence repeat

Each run: lumen surface -> snappyHexMesh (quarter model) -> steady mean-flow
initialisation (simpleFoam) -> pulsatile pimpleFoam (Womersley inlet) ->
evar/post_macro.py metrics. Cases go to cfd/runs/evar/<name>-<preset>/,
results to cfd/results/evar/<name>-<preset>/.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = os.path.dirname(HERE)
sys.path.insert(0, CFD)
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402

import run as R  # noqa: E402  (OpenFOAM backend + shell helper from the coronary pipeline)
import flow  # noqa: E402
import graft_geometry as gg  # noqa: E402
import macro_case as mc  # noqa: E402
sys.path.insert(0, os.path.join(CFD, "pipeline"))
import case as fc  # noqa: E402


def load_cfg():
    return json.load(open(os.path.join(CFD, "config", "evar.json")))


def waveform(cfg):
    inl = cfg["inflow"]
    path = os.path.join(CFD, inl["waveform_csv"])
    if not os.path.exists(path):
        t, f = flow.infrarenal_rest_shape(period=inl["period_s"])
        flow.write_csv(path, t, flow.rescale(f, inl["mean_mL_s"], inl["peak_mL_s"]))
    return flow.load_csv(path, inl["period_s"], inl["mean_mL_s"], inl["peak_mL_s"], samples=200)


def mesh(cdir, np_):
    R.sh(["blockMesh"], cdir, "blockMesh")
    R.sh(["decomposePar", "-force"], cdir, "decomposePar.mesh")
    R.sh(["snappyHexMesh", "-overwrite"], cdir, "snappyHexMesh", np_)
    R.sh(["reconstructParMesh", "-constant"], cdir, "reconstructParMesh")
    for d in glob.glob(os.path.join(cdir, "processor*")):
        shutil.rmtree(d)
    R.sh(["transformPoints", "-scale", "(0.001 0.001 0.001)"], cdir, "transformPoints")
    out = R.sh(["checkMesh", "-constant"], cdir, "checkMesh")
    import re
    info = {"cells": int(re.search(r"cells:\s+(\d+)", out).group(1)),
            "mesh_ok": "Mesh OK" in out,
            "failed_checks": int(m.group(1)) if (m := re.search(r"Failed (\d+) mesh checks", out)) else 0,
            "max_non_ortho": float(m.group(1)) if (m := re.search(r"Mesh non-orthogonality Max: ([\d.eE+-]+)", out)) else None,
            "max_skewness": float(m.group(1)) if (m := re.search(r"Max skewness = ([\d.eE+-]+)", out)) else None}
    lay = open(os.path.join(cdir, "log.snappyHexMesh")).read()
    import re as _re
    m = _re.findall(r"wall\s+\d+\s+[\d.]+\s+([\d.e+-]+)\s+([\d.]+)", lay)
    if m:
        info["layer_thickness_m"], info["layer_pct"] = float(m[-1][0]), float(m[-1][1])
    bnd = open(os.path.join(cdir, "constant/polyMesh/boundary")).read()
    info["patch_faces"] = {n: int(f) for n, f in _re.findall(r"(\w+)\s*\{[^}]*?nFaces\s+(\d+);", bnd)}
    return info


def solve(cdir, app, np_, all_times=False):
    R.sh(["decomposePar", "-force"], cdir, f"decomposePar.{app}")
    for pd in glob.glob(os.path.join(cdir, "processor*")):
        shutil.copytree(os.path.join(cdir, "constant/boundaryData"), os.path.join(pd, "constant/boundaryData"),
                        dirs_exist_ok=True)
    R.sh([app], cdir, app, np_)
    R.sh(["reconstructPar"] + ([] if all_times else ["-latestTime"]), cdir, f"reconstructPar.{app}")
    for pd in glob.glob(os.path.join(cdir, "processor*")):
        shutil.rmtree(pd)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=["validation", "graft"])
    ap.add_argument("--preset", default="coarse", choices=list(mc.PRESETS))
    ap.add_argument("--np", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--cycles", type=int, default=None)
    ap.add_argument("--stage", default="all", choices=["all", "mesh", "solve", "post"])
    args = ap.parse_args()

    cfg = load_cfg()
    R.setup_backend("auto")
    inl = cfg["inflow"]
    T = inl["period_s"]
    ts, q_mL_s = waveform(cfg)
    name = f"{args.case}-{args.preset}"
    cdir = os.path.join(CFD, "runs", "evar", name)
    rdir = os.path.join(CFD, "results", "evar", name)
    os.makedirs(rdir, exist_ok=True)

    if args.case == "validation":
        Rmm, L = cfg["validation"]["radius_mm"], cfg["validation"]["length_mm"]
        phys = fc.Physics(newtonian=True, mu_newtonian=cfg["validation"]["mu_Pa_s"])
        surf = lambda: gg.tube_surface(Rmm, L)                                  # noqa: E731
        box = ((0.0, L), (0.0, Rmm + 0.8), (0.0, Rmm + 0.8))
        loc = (0.37 * L, 0.41 * Rmm, 0.43 * Rmm)
        x_in = 0.0
        probes = [(x * 1e-3, Rmm * 1e-3 * math.cos(a), Rmm * 1e-3 * math.sin(a))
                  for x in (0.5 * L, 0.75 * L) for a in (math.radians(20), math.radians(45), math.radians(70))]
        n_cycles = args.cycles or cfg["validation"]["n_cycles"]
    else:
        g = gg.GraftParams(**cfg.get("graft", {}))
        Rmm = g.D_body / 2
        phys = fc.Physics(**cfg.get("physics", {}))
        surf = lambda: gg.graft_surface(g)                                       # noqa: E731
        c_end, r_end = gg.centreline(np.array([g.x_end]), g)
        box = ((-g.L_body, g.x_end), (0.0, float(c_end[0] + r_end[0]) + 0.8), (0.0, Rmm + 0.8))
        loc = (-0.5 * g.L_body + 0.013, 0.37 * Rmm, 0.41 * Rmm)
        x_in = -g.L_body
        probes = []
        for x in (-30.0, 30.0, 45.0, 62.0, 100.0):
            c, r = gg.centreline(np.array([x]), g)
            c, r = float(c[0]), float(r[0])
            probes += [(x, c + r, 0.0), (x, c, r)] + ([(x, max(c - r, 0.0), 0.0)] if x > g.L_trans else [])
        probes = [tuple(v * 1e-3 for v in p) for p in probes]
        n_cycles = args.cycles or cfg["pulsatile"]["n_cycles"]

    W = flow.Womersley(ts, q_mL_s * 1e-6, T, Rmm * 1e-3, mu=phys.mu_newtonian if phys.newtonian else 0.0035,
                       rho=phys.rho)
    meta = {"case": args.case, "preset": args.preset, "womersley_number_inlet": W.womersley_number,
            "Re_mean_inlet": flow.reynolds(np.mean(q_mL_s) * 1e-6, 2 * Rmm * 1e-3),
            "Re_peak_inlet": flow.reynolds(np.max(q_mL_s) * 1e-6, 2 * Rmm * 1e-3),
            "mean_mL_s": float(np.mean(q_mL_s)), "peak_mL_s": float(np.max(q_mL_s)),
            "min_mL_s": float(np.min(q_mL_s)), "period_s": T, "n_cycles": n_cycles}
    t0 = time.time()

    if args.stage in ("all", "mesh"):
        if os.path.exists(cdir):
            shutil.rmtree(cdir)
        os.makedirs(os.path.join(cdir, "constant/triSurface"))
        tm = surf()
        tm.export(os.path.join(cdir, "constant/triSurface/lumen.stl"))
        minfo = mc.write_mesh_dicts(cdir, box, args.preset, args.np, loc)
        mc.write_control_steady(cdir, 10)          # blockMesh/snappy need controlDict, fvSchemes, fvSolution
        mc.write_schemes(cdir, steady=True)
        minfo.update(mesh(cdir, args.np))
        print(f"    mesh: {minfo['cells']:,} cells, checkMesh "
              f"{'OK' if minfo['mesh_ok'] else str(minfo['failed_checks']) + ' failed checks'}")
        json.dump(minfo, open(os.path.join(cdir, "mesh.json"), "w"), indent=2)

    if args.stage in ("all", "solve"):
        for d in glob.glob(os.path.join(cdir, "[0-9]*")) + glob.glob(os.path.join(cdir, "postProcessing")):
            shutil.rmtree(d)
        shutil.rmtree(os.path.join(cdir, "constant/boundaryData"), ignore_errors=True)
        # 1) steady mean-flow initialisation (the slow-to-develop mean component)
        mc.write_fields(cdir)
        mc.write_physics(cdir, phys)
        mc.write_schemes(cdir, steady=True)
        mc.write_inlet_poiseuille(cdir, float(np.mean(q_mL_s)) * 1e-6, x_in, Rmm)
        n_it = cfg["steady_init_iter"]
        mc.write_control_steady(cdir, n_it)
        solve(cdir, "simpleFoam", args.np)
        last = sorted(glob.glob(os.path.join(cdir, "[1-9]*")), key=lambda p: float(os.path.basename(p)))[-1]
        for f in ("U", "p"):
            shutil.copy(os.path.join(last, f), os.path.join(cdir, "0", f))
        shutil.rmtree(last)
        shutil.rmtree(os.path.join(cdir, "postProcessing"), ignore_errors=True)
        # 2) pulsatile
        shutil.rmtree(os.path.join(cdir, "constant/boundaryData"))
        Tall = np.concatenate([ts + k * T for k in range(n_cycles)] + [[n_cycles * T]])
        mc.write_inlet_womersley(cdir, W, Tall, x_in, Rmm)
        mc.write_schemes(cdir, steady=False)
        pc = cfg["pulsatile"]
        mc.write_control_pulsatile(cdir, T, n_cycles, pc["max_co"], pc["snapshots_per_cycle"], probes_m=probes)
        solve(cdir, "pimpleFoam", args.np, all_times=True)
        meta["solve_min"] = round((time.time() - t0) / 60, 1)

    json.dump(meta, open(os.path.join(rdir, "run_meta.json"), "w"), indent=2)
    if args.stage in ("all", "post", "solve"):
        import post_macro
        post_macro.analyse(cdir, rdir, args.case, cfg, W, phys)
    print(f"done in {(time.time() - t0) / 60:.1f} min -> {rdir}")


if __name__ == "__main__":
    main()
