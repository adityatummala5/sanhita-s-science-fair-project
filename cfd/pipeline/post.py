"""Post-processing: hemodynamic metrics and figures from a finished case.

Definitions (Ku 1985; He & Ku 1996; Himburg 2004):
    TAWSS = (1/T) * integral |tau_w| dt                          [Pa]
    OSI   = 0.5 * (1 - |(1/T) integral tau_w dt| / TAWSS)        [0 .. 0.5]
    RRT   = 1 / ((1 - 2*OSI) * TAWSS) = 1 / |time-mean tau_w|   [1/Pa]
OpenFOAM's incompressible solvers report kinematic stress (m2/s2); it is
multiplied by rho to get Pa. In steady mode TAWSS is just |WSS| and OSI = 0.

Regions used for the statistics:
    wall_stented   artery wall between the first and last strut (x0..x1)
    wall_wake      artery wall from x1 to x1 + 1 mm (downstream of the stent)
    stent          all strut (+ denticle) surfaces
Low-WSS threshold 0.4 Pa (Malek, Alper & Izumo, JAMA 1999).
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402

LOW_WSS_PA = 0.4
HIGH_OSI = 0.1


def _read(cdir: str):
    foam = os.path.join(cdir, "case.foam")
    open(foam, "a").close()
    r = pv.OpenFOAMReader(foam)
    r.set_active_time_value(r.time_values[-1])
    r.enable_all_patch_arrays()
    r.cell_to_point_creation = False
    data = r.read()
    patches = {}
    bnd = data["boundary"]
    for k in bnd.keys():
        patches[k.split("/")[-1]] = bnd[k]
    return data["internalMesh"], patches, r.time_values[-1]


def _face_table(mesh: pv.DataSet, rho: float, mode: str):
    sized = mesh.compute_cell_sizes(length=False, volume=False)
    area = np.asarray(sized.cell_data["Area"])
    c = np.asarray(mesh.cell_centers().points)
    cd = mesh.cell_data
    if mode == "pulsatile" and "magWSSMean" in cd:
        tawss = rho * np.asarray(cd["magWSSMean"])
        mean_vec = rho * np.linalg.norm(np.asarray(cd["wallShearStressMean"]), axis=1)
    else:
        tawss = rho * np.linalg.norm(np.asarray(cd["wallShearStress"]), axis=1)
        mean_vec = tawss.copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        osi = np.where(tawss > 0, 0.5 * (1.0 - mean_vec / tawss), 0.0)
        osi = np.clip(osi, 0.0, 0.5)
        rrt = np.where(mean_vec > 1e-9, 1.0 / mean_vec, np.nan)
    return dict(c=c, area=area, tawss=tawss, osi=osi, rrt=rrt)


def _stats(t: dict, sel: np.ndarray, prefix: str) -> dict:
    a = t["area"][sel]
    if a.sum() == 0:
        return {}
    w = a / a.sum()
    rrt = t["rrt"][sel]
    ok = np.isfinite(rrt)
    order = np.argsort(rrt[ok])
    cw = np.cumsum(w[ok][order])
    p95 = float(rrt[ok][order][np.searchsorted(cw, 0.95 * cw[-1])]) if ok.any() else np.nan
    return {
        f"{prefix}_area_mm2": float(a.sum() * 1e6),
        f"{prefix}_TAWSS_mean_Pa": float((w * t["tawss"][sel]).sum()),
        f"{prefix}_lowWSS_area_pct": float(100 * w[t["tawss"][sel] < LOW_WSS_PA].sum()),
        f"{prefix}_OSI_mean": float((w * t["osi"][sel]).sum()),
        f"{prefix}_highOSI_area_pct": float(100 * w[t["osi"][sel] > HIGH_OSI].sum()),
        f"{prefix}_RRT_mean_perPa": float((w[ok] * rrt[ok]).sum() / w[ok].sum()) if ok.any() else np.nan,
        f"{prefix}_RRT_p95_perPa": p95,
    }


def _last_value(path_glob: str) -> float | None:
    files = sorted(glob.glob(path_glob))
    if not files:
        return None
    rows = [ln.split() for ln in open(files[-1]) if ln.strip() and not ln.startswith("#")]
    return float(rows[-1][-1]) if rows else None


def analyse(cdir: str, dom, phys, mode: str, name: str, inlet: dict | None = None) -> dict:
    vol, patches, t_end = _read(cdir)
    rho = phys.rho
    x0, x1 = dom.x0 * 1e-3, dom.x1 * 1e-3
    out = {"design": name, "mode": mode, "time_end": float(t_end)}

    wall = _face_table(patches["arteryWall"], rho, mode)
    stent = _face_table(patches["stent"], rho, mode) if "stent" in patches else None
    xw = wall["c"][:, 0]
    out.update(_stats(wall, (xw >= x0) & (xw <= x1), "wall_stented"))
    out.update(_stats(wall, (xw > x1) & (xw <= x1 + 1e-3), "wall_wake"))
    if stent is not None:
        out.update(_stats(stent, np.ones(len(stent["area"]), bool), "stent"))
    elif mode == "steady" and phys.newtonian and inlet:
        # validation against Poiseuille flow: tau_w = 4 mu Q / (pi R^3)
        q = (inlet["mean_mL_min"] if inlet.get("steady_flow", "mean") == "mean" else inlet["peak_mL_min"]) * 1e-6 / 60
        R = dom.radius * 1e-3
        tau = 4 * phys.mu_newtonian * q / (math.pi * R ** 3)
        dp = 8 * phys.mu_newtonian * (dom.xb - dom.xa) * 1e-3 * q / (math.pi * R ** 4)
        sel = (xw >= x0) & (xw <= x1)
        w = wall["area"][sel] / wall["area"][sel].sum()
        sim = float((w * wall["tawss"][sel]).sum())
        out.update({"analytic_WSS_Pa": tau, "simulated_WSS_Pa": sim,
                    "WSS_error_pct": 100 * (sim - tau) / tau, "analytic_pressure_drop_Pa": dp})

    # recirculation: cells with reversed (time-mean) axial velocity near the stent
    cd = vol.cell_data
    U = np.asarray(cd["UMean"] if (mode == "pulsatile" and "UMean" in cd) else cd["U"])
    cv = np.asarray(vol.compute_cell_sizes(length=False, area=False).cell_data["Volume"])
    xc = np.asarray(vol.cell_centers().points)[:, 0]
    reg = (xc >= x0) & (xc <= x1 + 1e-3)
    rev = reg & (U[:, 0] < 0)
    out["recirc_volume_mm3"] = float(cv[rev].sum() * 1e9)
    out["recirc_volume_pct"] = float(100 * cv[rev].sum() / cv[reg].sum())

    # pressure drop across the domain (outlet p = 0)
    pin = _last_value(os.path.join(cdir, "postProcessing/inletPressure/*/surfaceFieldValue.dat"))
    if mode == "pulsatile" and "pMean" in patches["inlet"].cell_data:
        pa = patches["inlet"]
        ar = np.asarray(pa.compute_cell_sizes(length=False, volume=False).cell_data["Area"])
        pin = float((np.asarray(pa.cell_data["pMean"]) * ar).sum() / ar.sum())
    out["pressure_drop_Pa"] = float(rho * pin) if pin is not None else None
    qin = _last_value(os.path.join(cdir, "postProcessing/inletFlow/*/surfaceFieldValue.dat"))
    qout = _last_value(os.path.join(cdir, "postProcessing/outletFlow/*/surfaceFieldValue.dat"))
    if qin is not None and qout is not None:
        out["inlet_flow_mL_min_sector"] = float(-qin * 6e7)
        out["mass_imbalance_pct"] = float(100 * abs(qin + qout) / abs(qin))

    if stent is not None:
        _plot_maps(cdir, name, wall, stent, dom, mode)
    if "analytic_pressure_drop_Pa" in out and out.get("pressure_drop_Pa"):
        out["pressure_drop_error_pct"] = 100 * (out["pressure_drop_Pa"] - out["analytic_pressure_drop_Pa"]) \
            / out["analytic_pressure_drop_Pa"]
    return out


def _plot_maps(cdir, name, wall, stent, dom, mode):
    lab = "TAWSS (Pa)" if mode == "pulsatile" else "WSS (Pa)"
    fig, axs = plt.subplots(1, 3 if mode == "pulsatile" else 2, figsize=(18 if mode == "pulsatile" else 12, 6))
    lim = (dom.x0 - 0.3, dom.x1 + 1.5)

    def unroll(c, r_ref):
        return c[:, 0] * 1e3, r_ref * np.arctan2(c[:, 2], c[:, 1])

    for ax, (tab, title, r_ref) in zip(axs[:2], [(wall, "artery wall", dom.radius), (stent, "stent surface", 1.5)]):
        x, s = unroll(tab["c"], r_ref)
        sel = (x > lim[0]) & (x < lim[1])
        sc = ax.scatter(x[sel], s[sel], c=tab["tawss"][sel], s=0.4, cmap="turbo", vmin=0, vmax=2.0, rasterized=True)
        ax.set_aspect("equal")
        ax.set_title(f"{name}: {lab} on {title}")
        ax.set_xlabel("axial position (mm)")
        ax.set_ylabel("circumferential arc (mm)")
        fig.colorbar(sc, ax=ax, shrink=0.8, label=lab)
    if mode == "pulsatile":
        x, s = unroll(wall["c"], dom.radius)
        sel = (x > lim[0]) & (x < lim[1])
        sc = axs[2].scatter(x[sel], s[sel], c=wall["osi"][sel], s=0.4, cmap="magma", vmin=0, vmax=0.5, rasterized=True)
        axs[2].set_aspect("equal")
        axs[2].set_title(f"{name}: OSI on artery wall")
        axs[2].set_xlabel("axial position (mm)")
        fig.colorbar(sc, ax=axs[2], shrink=0.8, label="OSI")
    plt.tight_layout()
    plt.savefig(os.path.join(cdir, "wss_maps.png"), dpi=110)
    plt.close(fig)


KEY_METRICS = [
    ("wall_stented_TAWSS_mean_Pa", "Artery wall TAWSS (Pa)\nhigher = better"),
    ("wall_stented_lowWSS_area_pct", "Wall area with TAWSS < 0.4 Pa (%)\nlower = better"),
    ("wall_stented_OSI_mean", "Artery wall OSI\nlower = better"),
    ("wall_stented_RRT_mean_perPa", "Artery wall RRT (1/Pa)\nlower = better"),
    ("stent_lowWSS_area_pct", "Stent area with TAWSS < 0.4 Pa (%)\nlower = better"),
    ("recirc_volume_mm3", "Recirculation volume (mm^3)\nlower = better"),
]


def summarise(run_root: str, res_root: str, dom):
    rows = []
    for mpath in sorted(glob.glob(os.path.join(run_root, "*", "metrics.json"))):
        rows.append(json.load(open(mpath)))
        png = os.path.join(os.path.dirname(mpath), "wss_maps.png")
        if os.path.exists(png):
            shutil.copy(png, os.path.join(res_root, f"wss_maps_{rows[-1]['design']}.png"))
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(res_root, "metrics.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    names = [r["design"] for r in rows]
    metrics = [(k, t) for k, t in KEY_METRICS if any(r.get(k) is not None for r in rows)]
    n = len(metrics)
    fig, axs = plt.subplots(2, (n + 1) // 2, figsize=(5.5 * ((n + 1) // 2), 9))
    for ax, (k, title) in zip(axs.ravel(), metrics):
        vals = [r.get(k) if r.get(k) is not None else np.nan for r in rows]
        cols = ["0.55" if nm == "bare" else ("tab:orange" if nm == "student_sector" else "tab:blue") for nm in names]
        ax.bar(range(len(vals)), vals, color=cols)
        ax.set_xticks(range(len(vals)), names, rotation=60, ha="right", fontsize=8)
        ax.set_title(title, fontsize=10)
        if "bare" in names:
            ax.axhline(vals[names.index("bare")], color="0.3", ls="--", lw=0.8)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    plt.suptitle("DentiFlow CFD comparison (grey = bare control, orange = Sanhita's sector)")
    plt.tight_layout()
    plt.savefig(os.path.join(res_root, "comparison.png"), dpi=110)
    plt.close(fig)
