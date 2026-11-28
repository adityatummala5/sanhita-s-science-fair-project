#!/usr/bin/env python3
"""EVAR helical-ridge experiment (design G1 in EVAR_TEXTURE_PLAN.md).

Pipeline per design:
    geometry  lumen.stl (+ ridge.stl) from evar/geometry_evar.py
    mesh      blockMesh + snappyHexMesh (inlet, outlet and symmetry plane are
              background-block faces), checkMesh
    steady    simpleFoam at the mean flow, used only as the initial condition
    pulsatile pimpleFoam, Womersley inlet, n_cycles cardiac cycles; time
              averages over the last cycle
    post      evar/post_evar.py

Usage (from cfd/):
    python3 run_evar.py                         # every design in config/evar_helical.json
    python3 run_evar.py --only smooth,ridge_h1_p78
    python3 run_evar.py --stage summary         # tables and figures from finished runs
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

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "pipeline"))
sys.path.insert(0, os.path.join(HERE, "evar"))

import case as fc              # noqa: E402
import geometry_evar as ge     # noqa: E402
import run as base             # noqa: E402  (OpenFOAM backend helpers)
import waveform_evar as wv     # noqa: E402

H = fc.HEADER
RUN_ROOT = os.path.join(HERE, "runs", "evar")
RES_ROOT = os.path.join(HERE, "results", "evar_helical")


def ridge_for(cfg_ridge: dict | None) -> ge.Ridge | None:
    return None if cfg_ridge is None else ge.Ridge(**cfg_ridge)


# --------------------------------------------------------------------------
# meshing
# --------------------------------------------------------------------------
def write_mesh_dicts(cdir: str, g: ge.EvarGeom, mesh_cfg: dict, has_ridge: bool, has_refine: bool, nproc: int):
    bg = mesh_cfg["bg_mm"]
    lw, lr = mesh_cfg["wall_level"], mesh_cfg["ridge_level"]
    x0, x1 = 0.0, g.x_out
    y0 = 0.0
    y1 = g.A[1] + (g.x_out - g.x_div) * math.tan(math.radians(g.theta_deg)) + g.R_limb / math.cos(
        math.radians(g.theta_deg)) + 4 * bg
    z1 = g.R_body + 4 * bg
    nx, ny, nz = (math.ceil((x1 - x0) / bg), math.ceil((y1 - y0) / bg), math.ceil(2 * z1 / bg))
    y1, z1 = y0 + ny * bg, nz * bg / 2
    fc._w(f"{cdir}/system/blockMeshDict", H.format(cls="dictionary", obj="blockMeshDict") + f"""
convertToMeters 1;   // mm; scaled to m after meshing
vertices
(
    ({x0} {y0} {-z1}) ({x1} {y0} {-z1}) ({x1} {y1} {-z1}) ({x0} {y1} {-z1})
    ({x0} {y0} {z1})  ({x1} {y0} {z1})  ({x1} {y1} {z1})  ({x0} {y1} {z1})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );
boundary
(
    inlet  {{ type patch; faces ((0 4 7 3)); }}
    outlet {{ type patch; faces ((1 2 6 5)); }}
    sym    {{ type symmetryPlane; faces ((0 1 5 4)); }}
    box    {{ type wall; faces ((3 7 6 2) (0 3 2 1) (4 5 6 7)); }}
);
""")
    geo = ["lumen.stl { type triSurfaceMesh; name lumen; }"]
    surf = [f"lumen {{ level ({lw} {lw}); patchInfo {{ type wall; }} }}"]
    regs = []
    if has_ridge:
        geo.append("ridge.stl { type triSurfaceMesh; name ridge; }")
        surf.append(f"ridge {{ level ({lr} {lr}); patchInfo {{ type wall; }} }}")
        regs.append(f"ridge {{ mode distance; levels ((0.15 {lr}) (0.6 {lr - 1})); }}")
    if has_refine:
        geo.append("ridgeRefine.stl { type triSurfaceMesh; name ridgeRefine; }")
        regs.append(f"ridgeRefine {{ mode distance; levels ((0.15 {lr}) (0.6 {lr - 1})); }}")
    nl = "\n        "
    fc._w(f"{cdir}/system/snappyHexMeshDict", H.format(cls="dictionary", obj="snappyHexMeshDict") + f"""
castellatedMesh true;
snap            true;
addLayers       true;

geometry
{{
        {nl.join(geo)}
}}

castellatedMeshControls
{{
    maxLocalCells 6000000;
    maxGlobalCells 30000000;
    minRefinementCells 0;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 3;
    features ( );
    refinementSurfaces
    {{
        {nl.join(surf)}
    }}
    resolveFeatureAngle 30;
    refinementRegions
    {{
        {nl.join(regs)}
    }}
    locationInMesh (10.0013 5.0017 0.3011);
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 100;
    nRelaxIter 5;
    nFeatureSnapIter 10;
    implicitFeatureSnap true;
    explicitFeatureSnap false;
    multiRegionFeatureSnap false;
}}

addLayersControls
{{
    relativeSizes true;
    layers {{ lumen {{ nSurfaceLayers {mesh_cfg.get('wall_layers', 3)}; }} }}
    expansionRatio 1.3;
    finalLayerThickness 0.5;     // last layer = half the local wall cell (0.25 mm -> 0.125 mm)
    minThickness 0.05;
    nGrow 0;
    featureAngle 100;
    slipFeatureAngle 30;
    nRelaxIter 5;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}}

meshQualityControls
{{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-30;
    minTetQuality 1e-15;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.05;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
    relaxed {{ maxNonOrtho 75; }}
}}

mergeTolerance 1e-6;
""")
    fc._w(f"{cdir}/system/decomposeParDict", H.format(cls="dictionary", obj="decomposeParDict") + f"""
numberOfSubdomains {nproc};
method scotch;
""")
    return {"background_cells": nx * ny * nz, "bg_mm": bg, "wall_cell_um": 1e3 * bg / 2 ** lw,
            "ridge_cell_um": 1e3 * bg / 2 ** lr}


def mesh(cdir: str, np_: int) -> dict:
    import re
    base.sh(["blockMesh"], cdir, "blockMesh")
    if np_ > 1:
        base.sh(["decomposePar", "-force"], cdir, "decomposePar.mesh")
        base.sh(["snappyHexMesh", "-overwrite"], cdir, "snappyHexMesh", np_)
        base.sh(["reconstructParMesh", "-constant"], cdir, "reconstructParMesh")
        for d in glob.glob(os.path.join(cdir, "processor*")):
            shutil.rmtree(d)
    else:
        base.sh(["snappyHexMesh", "-overwrite"], cdir, "snappyHexMesh")
    base.sh(["transformPoints", "-scale", "(0.001 0.001 0.001)"], cdir, "transformPoints")
    out = base.sh(["checkMesh", "-constant"], cdir, "checkMesh")
    info = {"cells": int(re.search(r"cells:\s+(\d+)", out).group(1)),
            "mesh_ok": "Mesh OK" in out,
            "failed_checks": int(m.group(1)) if (m := re.search(r"Failed (\d+) mesh checks", out)) else 0,
            "max_non_ortho": float(m.group(1)) if (m := re.search(r"Mesh non-orthogonality Max: ([\d.eE+-]+)", out)) else None,
            "max_skewness": float(m.group(1)) if (m := re.search(r"Max skewness = ([\d.eE+-]+)", out)) else None}
    bnd = open(os.path.join(cdir, "constant/polyMesh/boundary")).read()
    info["patch_faces"] = {n: int(f) for n, f in re.findall(r"(\w+)\s*\{[^}]*?nFaces\s+(\d+);", bnd)}
    return info


# --------------------------------------------------------------------------
# boundary conditions and inlet data
# --------------------------------------------------------------------------
def inlet_points(g: ge.EvarGeom):
    R = g.R_body * 1e-3
    pts = [(0.0, 0.0, 0.0)]
    for r in np.linspace(R / 40, R * 1.04, 42):
        for t in np.radians(np.linspace(-96, 96, 97)):
            pts.append((0.0, r * math.cos(t), r * math.sin(t)))
    pts = np.array(pts)
    return pts, np.hypot(pts[:, 1], pts[:, 2])


def write_inlet(cdir: str, g: ge.EvarGeom, times: np.ndarray, u_of_r):
    """u_of_r: array (n_times, n_points) of axial velocity (m/s)."""
    pts, _ = inlet_points(g)
    base_dir = f"{cdir}/constant/boundaryData/inlet"
    if os.path.isdir(base_dir):
        shutil.rmtree(base_dir)
    os.makedirs(base_dir)
    with open(f"{base_dir}/points", "w") as fh:
        fh.write(f"{len(pts)}\n(\n" + "\n".join(f"({p[0]:.9e} {p[1]:.9e} {p[2]:.9e})" for p in pts) + "\n)\n")
    for t, u in zip(times, u_of_r):
        d = f"{base_dir}/{fc._tname(t)}"
        os.makedirs(d)
        with open(f"{d}/U", "w") as fh:
            fh.write(f"{len(u)}\n(\n" + "\n".join(f"({v:.9e} 0 0)" for v in u) + "\n)\n")


def write_fields(cdir: str, has_ridge: bool):
    walls = ["graftWall", "box"] + (["ridge"] if has_ridge else [])
    U_w = "\n".join(f"    {w} {{ type noSlip; }}" for w in walls)
    p_w = "\n".join(f"    {w} {{ type zeroGradient; }}" for w in walls)
    fc._w(f"{cdir}/0/U", H.format(cls="volVectorField", obj="U") + f"""
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
    inlet
    {{
        type            timeVaryingMappedFixedValue;
        offset          (0 0 0);
        setAverage      off;
        mapMethod       planarInterpolation;
        value           uniform (0 0 0);
    }}
    outlet {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
    sym    {{ type symmetryPlane; }}
{U_w}
}}
""")
    fc._w(f"{cdir}/0/p", H.format(cls="volScalarField", obj="p") + f"""
dimensions [0 2 -2 0 0 0 0];   // kinematic pressure p/rho
internalField uniform 0;
boundaryField
{{
    inlet  {{ type zeroGradient; }}
    outlet {{ type fixedValue; value uniform 0; }}
    sym    {{ type symmetryPlane; }}
{p_w}
}}
""")


def write_fvsolution(cdir: str):
    fc._w(f"{cdir}/system/fvSolution", H.format(cls="dictionary", obj="fvSolution") + """
solvers
{
    p      { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.01; }
    pFinal { $p; relTol 0; }
    "(U|UFinal)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
}
PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 1;
    nCorrectors 2;
    nNonOrthogonalCorrectors 1;
}
SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent no;
    residualControl { p 1e-4; U 1e-5; }
}
relaxationFactors
{
    fields { p 0.3; }
    equations { U 0.7; }
}
""")


def function_objects(walls: str, avg_start: float | None, rho: float) -> str:
    fos = f"""
    wallShearStress
    {{
        type wallShearStress; libs ("libfieldFunctionObjects.so");
        patches {walls};
        executeControl timeStep; writeControl writeTime;
    }}
    magWSS
    {{
        type mag; libs ("libfieldFunctionObjects.so");
        field wallShearStress; result magWSS;
        executeControl timeStep; writeControl writeTime;
    }}
    inletFlow
    {{
        type surfaceFieldValue; libs ("libfieldFunctionObjects.so");
        regionType patch; name inlet; operation sum; fields (phi);
        writeFields false; writeControl timeStep; writeInterval 2;
    }}
    outletFlow {{ $inletFlow; name outlet; }}
    inletPressure
    {{
        type surfaceFieldValue; libs ("libfieldFunctionObjects.so");
        regionType patch; name inlet; operation areaAverage; fields (p);
        writeFields false; writeControl timeStep; writeInterval 2;
    }}
    wallWSS
    {{
        type surfaceFieldValue; libs ("libfieldFunctionObjects.so");
        regionType patch; name graftWall; operation areaAverage; fields (magWSS);
        writeFields false; writeControl timeStep; writeInterval 2;
    }}
    graftForces
    {{
        type forces; libs ("libforces.so");
        patches {walls};
        rho rhoInf; rhoInf {rho};
        CofR (0 0 0);
        writeControl timeStep; writeInterval 2;
    }}
"""
    if avg_start is not None:
        fos += f"""
    vorticity
    {{
        type vorticity; libs ("libfieldFunctionObjects.so");
        executeControl writeTime; writeControl writeTime;
    }}
    cycleAverage
    {{
        type fieldAverage; libs ("libfieldFunctionObjects.so");
        timeStart {avg_start};
        restartOnRestart false; restartOnOutput false; periodicRestart false;
        executeControl timeStep; writeControl writeTime;
        fields
        (
            wallShearStress {{ mean on; prime2Mean off; base time; }}
            magWSS          {{ mean on; prime2Mean off; base time; }}
            U               {{ mean on; prime2Mean off; base time; }}
            p               {{ mean on; prime2Mean off; base time; }}
        );
    }}
"""
    return "functions\n{" + fos + "}\n"


def write_control(cdir: str, app: str, walls: str, rho: float, **kw):
    if app == "simpleFoam":
        body = f"""
application     simpleFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {kw['n_iter']};
deltaT          1;
writeControl    timeStep;
writeInterval   {kw['n_iter']};
purgeWrite      1;
writeFormat     binary;
writePrecision  8;
timeFormat      general;
runTimeModifiable true;
"""
        avg = None
    else:
        T, n, snaps = kw["period"], kw["n_cycles"], kw["snapshots"]
        body = f"""
application     pimpleFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {T * n:.6f};
deltaT          1e-4;
adjustTimeStep  yes;
maxCo           {kw['max_co']};
maxDeltaT       {T / 200:.6f};
writeControl    adjustableRunTime;
writeInterval   {T / snaps:.6f};
purgeWrite      {snaps};
writeFormat     binary;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable true;
"""
        avg = T * (n - 1)
    fc._w(f"{cdir}/system/controlDict", H.format(cls="dictionary", obj="controlDict") + body
          + function_objects(walls, avg, rho))


def solve(cdir: str, app: str, np_: int):
    if np_ > 1:
        base.sh(["decomposePar", "-force"], cdir, f"decomposePar.{app}")
        for pd in glob.glob(os.path.join(cdir, "processor*")):
            shutil.copytree(os.path.join(cdir, "constant/boundaryData"),
                            os.path.join(pd, "constant/boundaryData"), dirs_exist_ok=True)
        base.sh([app], cdir, app, np_)
        base.sh(["reconstructPar", "-newTimes"], cdir, f"reconstructPar.{app}")
        for pd in glob.glob(os.path.join(cdir, "processor*")):
            shutil.rmtree(pd)
    else:
        base.sh([app], cdir, app)


def time_dirs(cdir: str):
    out = []
    for d in os.listdir(cdir):
        try:
            out.append((float(d), d))
        except ValueError:
            pass
    return sorted(out)


# --------------------------------------------------------------------------
def run_design(d: dict, cfg: dict, designs: dict, np_: int, stages: set):
    g = ge.EvarGeom()
    cdir = os.path.join(RUN_ROOT, d["name"])
    os.makedirs(cdir, exist_ok=True)
    ridge = ridge_for(d.get("ridge"))
    refine = ridge_for(designs[d["refine_like"]]["ridge"]) if d.get("refine_like") else None
    walls = "(graftWall ridge)" if ridge else "(graftWall)"
    ph = fc.Physics(**cfg["physics"])
    T, ncyc = cfg["period_s"], cfg["n_cycles"]
    meta_path = os.path.join(cdir, "run.json")
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {"design": d}
    print(f"\n=== {d['name']}: {d['description']}")

    if "mesh" in stages:
        for sub in ("constant/polyMesh", "0", "processor0"):
            shutil.rmtree(os.path.join(cdir, sub), ignore_errors=True)
        tri = os.path.join(cdir, "constant/triSurface")
        meta["geometry"] = ge.write_all(tri, g, ridge, refine)
        # snappy names the lumen patch after the surface; rename it graftWall below
        meta["mesh_dicts"] = write_mesh_dicts(cdir, g, cfg["mesh"], ridge is not None, refine is not None, np_)
        fc.write_physics(cdir, ph)
        fc.write_schemes(cdir, steady=True)
        write_fvsolution(cdir)
        write_control(cdir, "simpleFoam", walls.replace("graftWall", "lumen"), ph.rho, n_iter=10)
        t0 = time.time()
        meta["mesh"] = mesh(cdir, np_)
        meta["mesh"]["minutes"] = (time.time() - t0) / 60
        bpath = os.path.join(cdir, "constant/polyMesh/boundary")
        txt = open(bpath).read().replace("    lumen\n", "    graftWall\n")
        open(bpath, "w").write(txt)
        json.dump(meta, open(meta_path, "w"), indent=2)
        print(f"    mesh: {meta['mesh']['cells']:,} cells, patches {meta['mesh']['patch_faces']}, "
              f"checkMesh {'OK' if meta['mesh']['mesh_ok'] else str(meta['mesh']['failed_checks']) + ' failed'}")

    if "steady" in stages:
        for _, td in time_dirs(cdir):
            shutil.rmtree(os.path.join(cdir, td))
        shutil.rmtree(os.path.join(cdir, "postProcessing"), ignore_errors=True)
        write_fields(cdir, ridge is not None)
        ts, q = wv.load_cycle(os.path.join(HERE, cfg["waveform_csv"]), T)
        _, r = inlet_points(g)
        nu = ph.mu_inf / ph.rho
        u_mean = wv.womersley_profile(np.full(4, q.mean() * 1e-6), T, g.R_body * 1e-3, nu, r)[0]
        write_inlet(cdir, g, np.array([0.0, 1e6]), [u_mean, u_mean])
        fc.write_schemes(cdir, steady=True)
        write_control(cdir, "simpleFoam", walls, ph.rho, n_iter=600)
        t0 = time.time()
        solve(cdir, "simpleFoam", np_)
        meta["steady_minutes"] = (time.time() - t0) / 60
        json.dump(meta, open(meta_path, "w"), indent=2)

    if "pulsatile" in stages:
        tds = time_dirs(cdir)
        last = tds[-1][1]
        assert tds[-1][0] > 0, "run the steady stage first"
        # the steady solution becomes the initial condition
        shutil.rmtree(os.path.join(cdir, "0"), ignore_errors=True)
        os.rename(os.path.join(cdir, last), os.path.join(cdir, "0"))
        for junk in ("uniform", "phi", "wallShearStress", "magWSS"):
            pth = os.path.join(cdir, "0", junk)
            shutil.rmtree(pth, ignore_errors=True) if os.path.isdir(pth) else (os.path.exists(pth) and os.remove(pth))
        shutil.move(os.path.join(cdir, "postProcessing"), os.path.join(cdir, "postProcessing_steady")) \
            if os.path.isdir(os.path.join(cdir, "postProcessing")) else None
        ts, q = wv.load_cycle(os.path.join(HERE, cfg["waveform_csv"]), T, n=200)
        _, r = inlet_points(g)
        u = wv.womersley_profile(q * 1e-6, T, g.R_body * 1e-3, ph.mu_inf / ph.rho, r)
        times = np.concatenate([ts + k * T for k in range(ncyc)] + [[ncyc * T]])
        U = np.concatenate([u] * ncyc + [u[:1]])
        write_inlet(cdir, g, times, U)
        fc.write_schemes(cdir, steady=False)
        write_control(cdir, "pimpleFoam", walls, ph.rho, period=T, n_cycles=ncyc,
                      snapshots=cfg["snapshots_per_cycle"], max_co=cfg["max_co"])
        t0 = time.time()
        solve(cdir, "pimpleFoam", np_)
        meta["pulsatile_minutes"] = (time.time() - t0) / 60
        json.dump(meta, open(meta_path, "w"), indent=2)

    if "post" in stages:
        import post_evar
        res = post_evar.analyse(cdir, g, ph, T, d["name"], os.path.join(RES_ROOT, d["name"]))
        res.update({"cells": meta.get("mesh", {}).get("cells"), "mesh_ok": meta.get("mesh", {}).get("mesh_ok"),
                    "pulsatile_minutes": meta.get("pulsatile_minutes")})
        json.dump(res, open(os.path.join(RES_ROOT, d["name"], "metrics.json"), "w"), indent=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config/evar_helical.json")
    ap.add_argument("--only", default="")
    ap.add_argument("--np", type=int, default=4)
    ap.add_argument("--stage", default="mesh,steady,pulsatile,post",
                    help="comma list of mesh, steady, pulsatile, post, summary")
    ap.add_argument("--backend", default="auto", choices=["auto", "docker", "native"])
    a = ap.parse_args()
    cfg = json.load(open(os.path.join(HERE, a.config)))
    designs = {d["name"]: d for d in cfg["designs"]}
    sel = [designs[n] for n in a.only.split(",")] if a.only else cfg["designs"]
    stages = set(a.stage.split(","))
    os.makedirs(RES_ROOT, exist_ok=True)
    if stages - {"summary"}:
        base.setup_backend(a.backend)
        for d in sel:
            run_design(d, cfg, designs, a.np, stages)
    if "summary" in stages:
        import post_evar
        post_evar.summarise(RES_ROOT, [d["name"] for d in cfg["designs"]])


if __name__ == "__main__":
    main()
