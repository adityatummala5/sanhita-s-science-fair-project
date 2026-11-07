#!/usr/bin/env python3
"""DentiFlow CFD pipeline: geometry -> mesh -> OpenFOAM solve -> hemodynamic metrics.

Examples
--------
  # list designs in the config
  python3 run.py --list

  # quick smoke test: steady flow, coarse mesh, two designs
  python3 run.py --preset coarse --mode steady --only bare,student_sector

  # full experiment: pulsatile LAD flow, every design, 4 cores
  python3 run.py --preset medium --mode pulsatile --np 4

  # only build the cleaned stent STL files (no CFD)
  python3 run.py --geometry-only

Outputs
-------
  geometry_out/<design>.stl           cleaned, watertight stent sector (mm)
  runs/<preset>-<mode>/<design>/      OpenFOAM case (mesh, logs, results)
  results/<preset>-<mode>/metrics.csv one row of hemodynamic metrics per design
  results/<preset>-<mode>/*.png       wall-shear-stress maps and comparison charts
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "pipeline"))

import case as fc        # noqa: E402
import geometry as geo   # noqa: E402
import waveform as wf    # noqa: E402

DOCKER_IMAGE = os.environ.get("DENTIFLOW_IMAGE", "opencfd/openfoam-default:2406")
MPI_ENV = {"OMPI_ALLOW_RUN_AS_ROOT": "1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
           "OMPI_MCA_btl_vader_single_copy_mechanism": "none"}
BACKEND = {"kind": None, "env": None}


def setup_backend(kind: str):
    """'docker' runs every OpenFOAM command inside the official OpenFOAM image
    (works on Linux, macOS and Windows). 'native' uses an OpenFOAM install on
    this machine (must be openfoam.com v2206 or newer; set FOAM_BASHRC if it
    is not on the default path)."""
    if kind == "auto":
        ok = shutil.which("docker") and subprocess.run(["docker", "info"], capture_output=True).returncode == 0
        kind = "docker" if ok else "native"
    if kind == "docker":
        if subprocess.run(["docker", "image", "inspect", DOCKER_IMAGE], capture_output=True).returncode != 0:
            print(f"pulling {DOCKER_IMAGE} (one-time, ~500 MB download)...")
            subprocess.run(["docker", "pull", DOCKER_IMAGE], check=True)
    else:
        rcs = [os.environ.get("FOAM_BASHRC", ""),
               *sorted(glob.glob("/usr/lib/openfoam/openfoam*/etc/bashrc"), reverse=True),
               *sorted(glob.glob("/opt/openfoam*/etc/bashrc"), reverse=True)]
        rc = next((p for p in rcs if p and os.path.exists(p)), None)
        if rc is None:
            sys.exit("No usable OpenFOAM found. Install Docker (recommended) or OpenFOAM v2206+, "
                     "or set FOAM_BASHRC=/path/to/openfoam/etc/bashrc")
        dump = subprocess.run(["bash", "-c", f"source '{rc}' >/dev/null 2>&1; env -0"],
                              capture_output=True, check=True).stdout.decode()
        env = dict(kv.split("=", 1) for kv in dump.split("\0") if "=" in kv)
        env.update(MPI_ENV)
        BACKEND["env"] = env
    BACKEND["kind"] = kind
    print(f"OpenFOAM backend: {kind}" + (f" ({DOCKER_IMAGE})" if kind == "docker" else ""))


def sh(cmd: list[str], cwd: str, log: str, np_: int = 1) -> str:
    """Run an OpenFOAM utility in case directory cwd, logging to cwd/log.<log>."""
    if np_ > 1:
        cmd = ["mpirun", "-np", str(np_), "--oversubscribe"] + cmd + ["-parallel"]
    if BACKEND["kind"] == "docker":
        inner = ("source /usr/lib/openfoam/openfoam*/etc/bashrc >/dev/null 2>&1; cd /case; "
                 + " ".join(shlex.quote(c) for c in cmd))
        full = ["docker", "run", "--rm", "--shm-size=2g", "-v", f"{os.path.abspath(cwd)}:/case"]
        if hasattr(os, "getuid") and os.getuid() != 0:
            full += ["--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp"]
        for k, v in MPI_ENV.items():
            full += ["-e", f"{k}={v}"]
        full += [DOCKER_IMAGE, "bash", "-c", inner]
        env = None
    else:
        full, env = cmd, BACKEND["env"]
    path = os.path.join(cwd, f"log.{log}")
    t0 = time.time()
    with open(path, "w") as fh:
        r = subprocess.run(full, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    out = open(path).read()
    if r.returncode != 0 or "FOAM FATAL" in out:
        tail = "\n".join(out.splitlines()[-30:])
        raise RuntimeError(f"{' '.join(cmd)} failed (see {path}):\n{tail}")
    print(f"    {log:<24s} ok  ({time.time() - t0:6.1f} s)")
    return out


def rel(p: str) -> str:
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(HERE, p))


# --------------------------------------------------------------------------
def build_geometry(design: dict, cfg: dict, win: geo.Window, cache: dict) -> geo.StentResult:
    kind = design["type"]
    if kind == "empty":
        return None
    if "bare" not in cache:
        cache["bare"] = geo.bare_sector(rel(cfg["stent_stl"]), win)
    if kind == "bare":
        return cache["bare"]
    if kind == "parametric":
        p = geo.DenticleParams(**design.get("params", {}))
        return geo.denticle_sector(cache["bare"], p, win)
    if kind == "import":
        src = rel(design["stl"])
        if src.endswith(".zip"):
            dest = os.path.join(HERE, "geometry_out", "_unzipped")
            with zipfile.ZipFile(src) as z:
                name = [n for n in z.namelist() if n.lower().endswith(".stl")][0]
                z.extract(name, dest)
            src = os.path.join(dest, name)
        return geo.imported_sector(src)
    raise ValueError(f"unknown design type {kind}")


def mesh_case(cdir: str, preset: str, np_: int, parallel_snappy: bool) -> dict:
    sh(["blockMesh"], cdir, "blockMesh")
    sh(["surfaceFeatureExtract"], cdir, "surfaceFeatureExtract")
    if parallel_snappy and np_ > 1:
        sh(["decomposePar", "-force"], cdir, "decomposePar.mesh")
        sh(["snappyHexMesh", "-overwrite"], cdir, "snappyHexMesh", np_)
        sh(["reconstructParMesh", "-constant"], cdir, "reconstructParMesh")
        for d in glob.glob(os.path.join(cdir, "processor*")):
            shutil.rmtree(d)
    else:
        sh(["snappyHexMesh", "-overwrite"], cdir, "snappyHexMesh")
    sh(["transformPoints", "-scale", "(0.001 0.001 0.001)"], cdir, "transformPoints")
    out = sh(["checkMesh", "-constant"], cdir, "checkMesh")
    info = {"cells": int(re.search(r"cells:\s+(\d+)", out).group(1)),
            "mesh_ok": "Mesh OK" in out,
            "failed_checks": int(m.group(1)) if (m := re.search(r"Failed (\d+) mesh checks", out)) else 0,
            "max_non_ortho": float(m.group(1)) if (m := re.search(r"Mesh non-orthogonality Max: ([\d.eE+-]+)", out)) else None,
            "max_skewness": float(m.group(1)) if (m := re.search(r"Max skewness = ([\d.eE+-]+)", out)) else None}
    # boundary summary
    bnd = open(os.path.join(cdir, "constant/polyMesh/boundary")).read()
    info["patch_faces"] = {n: int(f) for n, f in re.findall(r"(\w+)\s*\{[^}]*?nFaces\s+(\d+);", bnd)}
    return info


def solve_case(cdir: str, mode: str, np_: int):
    app = "pimpleFoam" if mode == "pulsatile" else "simpleFoam"
    if np_ > 1:
        sh(["decomposePar", "-force"], cdir, "decomposePar")
        for pd in glob.glob(os.path.join(cdir, "processor*")):
            shutil.copytree(os.path.join(cdir, "constant/boundaryData"),
                            os.path.join(pd, "constant/boundaryData"), dirs_exist_ok=True)
        sh([app], cdir, app, np_)
        sh(["reconstructPar", "-latestTime"], cdir, "reconstructPar")
    else:
        sh([app], cdir, app)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config/experiment.json")
    ap.add_argument("--preset", default="coarse", choices=list(fc.PRESETS))
    ap.add_argument("--mode", default="steady", choices=["steady", "pulsatile"])
    ap.add_argument("--np", type=int, default=max(1, (os.cpu_count() or 2)), help="MPI processes")
    ap.add_argument("--only", default="", help="comma-separated design names")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--geometry-only", action="store_true")
    ap.add_argument("--skip-existing", action="store_true", help="skip designs that already have results")
    ap.add_argument("--serial-snappy", action="store_true", help="run snappyHexMesh on one core")
    ap.add_argument("--backend", default="auto", choices=["auto", "docker", "native"])
    ap.add_argument("--cycles", type=int, default=None, help="override number of cardiac cycles")
    ap.add_argument("--end-time", type=float, default=None,
                    help="SMOKE TEST ONLY: stop the pulsatile run early at this time (s)")
    args = ap.parse_args()

    cfg = json.load(open(rel(args.config)))
    designs = cfg["designs"]
    if args.list:
        for d in designs:
            print(f"{d['name']:<22s} {d['type']:<11s} {d.get('description', '')}")
        return
    if args.only:
        want = set(args.only.split(","))
        designs = [d for d in designs if d["name"] in want]
        missing = want - {d["name"] for d in designs}
        if missing:
            sys.exit(f"unknown design(s): {', '.join(sorted(missing))}")

    if not args.geometry_only:
        setup_backend(args.backend)
    win = geo.Window(**cfg["window"])
    dom = fc.Domain(x0=win.x0, x1=win.x1, theta0_deg=win.theta0_deg, theta1_deg=win.theta1_deg,
                    **cfg.get("domain", {}))
    phys = fc.Physics(**cfg.get("physics", {}))
    inl = cfg["inlet"]
    tag = f"{args.preset}-{args.mode}"
    run_root = os.path.join(HERE, "runs", tag)
    res_root = os.path.join(HERE, "results", tag)
    os.makedirs(os.path.join(HERE, "geometry_out"), exist_ok=True)
    os.makedirs(res_root, exist_ok=True)

    cache: dict = {}
    for d in designs:
        name = d["name"]
        print(f"\n=== {name} ({d['type']}) ===")
        t0 = time.time()
        g = build_geometry(d, cfg, win, cache)
        if g is not None:
            g.mesh.export(os.path.join(HERE, "geometry_out", f"{name}.stl"))
            ginfo = {"design": d, **geo.describe(g.mesh),
                     **{k: v for k, v in g.info.items() if k != "params"}}
            json.dump(ginfo, open(os.path.join(HERE, "geometry_out", f"{name}.json"), "w"), indent=2)
            print(f"    geometry                 ok  ({time.time() - t0:6.1f} s)  "
                  f"{ginfo['triangles']} tris, watertight={ginfo['watertight']}, "
                  f"protrusion={ginfo['max_protrusion_um']:.1f} um"
                  + (f", {ginfo['n_denticles']} denticles" if "n_denticles" in ginfo else ""))
        else:
            print("    geometry                 none (empty-vessel validation case)")
        if args.geometry_only:
            continue
        phys_d = fc.Physics(**{**cfg.get("physics", {}), **d.get("physics", {})})

        cdir = os.path.join(run_root, name)
        if args.skip_existing and os.path.exists(os.path.join(cdir, "metrics.json")):
            print("    results exist, skipping")
            continue
        if os.path.exists(cdir):
            shutil.rmtree(cdir)
        os.makedirs(os.path.join(cdir, "constant/triSurface"))
        if g is not None:
            g.mesh.export(os.path.join(cdir, "constant/triSurface/stent.stl"))
        geo.write_multisolid_stl(os.path.join(cdir, "constant/triSurface/lumen.stl"),
                                 geo.lumen_surface(win, dom.xa - 0.3, dom.xb + 0.3, dom.radius))
        minfo = fc.write_mesh_dicts(cdir, dom, args.preset, args.np, has_stent=g is not None)
        fc.write_physics(cdir, phys_d)
        fc.write_schemes(cdir, steady=(args.mode == "steady"))

        if args.mode == "pulsatile":
            pc = dict(cfg["pulsatile"])
            if args.cycles:
                pc["n_cycles"] = args.cycles
            T, Q = wf.load(rel(inl["waveform_csv"]), inl["period_s"], inl.get("mean_mL_min"),
                           inl.get("peak_mL_min"), pc["n_cycles"])
            fc.write_inlet_data(cdir, dom, T, Q)
            fc.write_control_pulsatile(cdir, inl["period_s"], pc["n_cycles"], pc["max_co"],
                                       pc["snapshots_per_cycle"])
            if args.end_time:
                cd = os.path.join(cdir, "system/controlDict")
                txt = open(cd).read()
                txt = re.sub(r"\nendTime\s+[^;]+;", f"\nendTime         {args.end_time};", txt)
                txt = re.sub(r"\nwriteInterval\s+[^;]+;", f"\nwriteInterval   {args.end_time / 2};", txt)
                open(cd, "w").write(txt)
        else:
            q = inl["mean_mL_min"] if inl.get("steady_flow", "mean") == "mean" else inl["peak_mL_min"]
            fc.write_inlet_data(cdir, dom, [0.0, 1e6], [q, q])
            fc.write_control_steady(cdir, cfg["steady"]["n_iter"])

        mesh_info = mesh_case(cdir, args.preset, args.np, not args.serial_snappy)
        fc.write_fields(cdir, has_stent=g is not None)   # after meshing, so snappy cannot overwrite BCs
        mesh_info.update(minfo)
        print(f"    mesh: {mesh_info['cells']:,} cells, checkMesh "
              f"{'OK' if mesh_info['mesh_ok'] else str(mesh_info['failed_checks']) + ' failed checks'}")
        json.dump(mesh_info, open(os.path.join(cdir, "mesh.json"), "w"), indent=2)

        solve_case(cdir, args.mode, args.np)

        import post  # noqa: E402  (imports pyvista; only needed here)
        m = post.analyse(cdir, dom, phys_d, args.mode, name, inl)
        m.update({"cells": mesh_info["cells"], "mesh_ok": mesh_info["mesh_ok"],
                  "runtime_min": round((time.time() - t0) / 60, 2)})
        json.dump(m, open(os.path.join(cdir, "metrics.json"), "w"), indent=2)
        print(f"    done in {(time.time() - t0) / 60:.1f} min")

    import preview  # noqa: E402
    preview.preview(os.path.join(HERE, "geometry_out"), os.path.join(HERE, "geometry_out", "preview.png"))
    if not args.geometry_only:
        import post  # noqa: E402
        post.summarise(run_root, res_root, dom)
        print(f"\nResults: {res_root}")


if __name__ == "__main__":
    main()
