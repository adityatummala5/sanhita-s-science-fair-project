"""Post-processing for the EVAR helical-ridge runs.

Wall metrics (Ku 1985; He & Ku 1996; Himburg 2004), time-averaged over the
LAST cardiac cycle:
    TAWSS = mean |tau_w|                [Pa]
    OSI   = 0.5 (1 - |mean tau_w| / TAWSS)
    RRT   = 1 / |mean tau_w|            [1/Pa]
    low-WSS area: TAWSS < 0.4 Pa (Malek et al. 1999)
    high-OSI area: OSI > 0.1

Graft-wall regions (graftWall + ridge faces):
    body          x < x_taper (straight Ø26 mm main body)
    bifurcation   taper, flow divider and the first 12 mm of the limb
    limb_ridge    limb, 12-87 mm from the divider (where the ridge is)
    limb_distal   limb, 87 mm to 5 mm before the outlet
    limb          limb_ridge + limb_distal

Flow metrics in the limb lumen (12 mm from the divider to 5 mm before the
outlet), averaged over the last-cycle snapshots:
    h1   mean helicity density u.w            [m/s^2]  (signed: sense of rotation)
    h2   mean |u.w|                           [m/s^2]  (helicity intensity, Gallo et al. 2012)
    LNH  mean |u.w| / (|u||w|)                (0 = no helical structure, 1 = fully helical)
    secondary-flow ratio  mean |u_inplane| / mean |u_axial|

Global: pressure drop (inlet minus outlet, time mean and peak) and the
displacement force on the graft (pressure + shear, x component, both limbs
= 2 x half model). Pressures are gauge values relative to the outlet, so the
force is the flow-induced part only; the static ~100 mmHg part is the same
for every design and cancels in the comparison.
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402

LOW_WSS_PA = 0.4
HIGH_OSI = 0.1
RIDGE_S = (12.0, 87.0)


# --------------------------------------------------------------------------
def _reader(cdir):
    foam = os.path.join(cdir, "case.foam")
    open(foam, "a").close()
    r = pv.OpenFOAMReader(foam)
    r.cell_to_point_creation = False
    return r


def _patch_table(block, rho):
    sized = block.compute_cell_sizes(length=False, volume=False)
    area = np.asarray(sized.cell_data["Area"]) * 1e6                 # mm2
    c = np.asarray(block.cell_centers().points) * 1e3                # mm
    cd = block.cell_data
    tawss = rho * np.asarray(cd["magWSSMean"])
    mean_vec = rho * np.linalg.norm(np.asarray(cd["wallShearStressMean"]), axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        osi = np.clip(np.where(tawss > 0, 0.5 * (1 - mean_vec / tawss), 0.0), 0, 0.5)
        rrt = np.where(mean_vec > 1e-9, 1.0 / mean_vec, np.nan)
    return dict(c=c, area=area, tawss=tawss, osi=osi, rrt=rrt)


def _cat(a, b):
    return {k: np.concatenate([a[k], b[k]]) for k in a}


def limb_coords(g, c):
    d = g.d
    e1 = np.array([-d[1], d[0], 0.0])
    e2 = np.array([0.0, 0.0, 1.0])
    v = c - g.A
    s = v @ d
    rv = v - np.outer(s, d)
    rl = np.linalg.norm(rv, axis=1)
    phi = np.degrees(np.arctan2(rv @ e2, rv @ e1))
    return s, rl, phi


def regions(g, c):
    s, rl, _ = limb_coords(g, c)
    x = c[:, 0]
    in_limb = (rl < g.R_limb + 1.5) & (s >= RIDGE_S[0])
    reg = {
        "body": x < g.x_taper,
        "bifurcation": (x >= g.x_taper) & ~in_limb,
        "limb_ridge": in_limb & (s <= RIDGE_S[1]),
        "limb_distal": in_limb & (s > RIDGE_S[1]) & (x <= g.x_out - 5.0),
    }
    reg["limb"] = reg["limb_ridge"] | reg["limb_distal"]
    reg["all"] = x <= g.x_out - 5.0
    return reg


def stats(t, sel, prefix):
    a = t["area"][sel]
    if a.sum() == 0:
        return {}
    w = a / a.sum()
    rrt = t["rrt"][sel]
    ok = np.isfinite(rrt)
    return {
        f"{prefix}_area_mm2": float(a.sum()),
        f"{prefix}_TAWSS_Pa": float((w * t["tawss"][sel]).sum()),
        f"{prefix}_TAWSS_p05_Pa": float(_wpct(t["tawss"][sel], w, 5)),
        f"{prefix}_lowWSS_area_pct": float(100 * w[t["tawss"][sel] < LOW_WSS_PA].sum()),
        f"{prefix}_OSI": float((w * t["osi"][sel]).sum()),
        f"{prefix}_highOSI_area_pct": float(100 * w[t["osi"][sel] > HIGH_OSI].sum()),
        f"{prefix}_RRT_perPa": float((w[ok] * rrt[ok]).sum() / w[ok].sum()),
        f"{prefix}_RRT_p95_perPa": float(_wpct(rrt[ok], w[ok], 95)),
    }


def _wpct(v, w, q):
    o = np.argsort(v)
    cw = np.cumsum(w[o])
    return v[o][min(np.searchsorted(cw, q / 100 * cw[-1]), len(v) - 1)]


# --------------------------------------------------------------------------
def _dat(path_glob):
    """Concatenate OpenFOAM function-object .dat files (restarts give several)."""
    rows = []
    for f in sorted(glob.glob(path_glob)):
        for ln in open(f):
            if ln.startswith("#") or not ln.strip():
                continue
            vals = [float(v) for v in re.sub(r"[()]", " ", ln).split()]
            rows.append(vals)
    if not rows:
        return None
    a = np.array(rows)
    _, idx = np.unique(a[:, 0], return_index=True)
    return a[idx]


def _cycle_mean(a, col, t0, t1):
    m = (a[:, 0] >= t0) & (a[:, 0] <= t1)
    t, v = a[m, 0], a[m, col]
    return float(np.trapezoid(v, t) / (t[-1] - t[0])), v, t


# --------------------------------------------------------------------------
def analyse(cdir: str, g, ph, period: float, name: str, out_dir: str) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    rho = ph.rho
    r = _reader(cdir)
    times = np.array(r.time_values)
    t_end = times[-1]
    t0 = t_end - period
    out = {"design": name, "time_end_s": float(t_end), "cycles_simulated": float(t_end / period)}

    # ---- wall metrics (last time holds the last-cycle averages) ----
    r.set_active_time_value(t_end)
    r.enable_all_patch_arrays()
    data = r.read()
    bnd = {k.split("/")[-1]: data["boundary"][k] for k in data["boundary"].keys()}
    wall = _patch_table(bnd["graftWall"], rho)
    has_ridge = "ridge" in bnd and bnd["ridge"].n_cells > 0
    surf = _cat(wall, _patch_table(bnd["ridge"], rho)) if has_ridge else wall
    reg = regions(g, surf["c"])
    for k, sel in reg.items():
        out.update(stats(surf, sel, k))
    if has_ridge:
        rt = _patch_table(bnd["ridge"], rho)
        out.update(stats(rt, np.ones(len(rt["area"]), bool), "ridgeSurface"))
    s, rl, phi = limb_coords(g, surf["c"])
    np.savez_compressed(os.path.join(out_dir, "wall_faces.npz"), c=surf["c"], area=surf["area"],
                        tawss=surf["tawss"], osi=surf["osi"], rrt=surf["rrt"], s=s, phi=phi,
                        is_ridge=np.r_[np.zeros(len(wall["area"]), bool),
                                       np.ones(len(surf["area"]) - len(wall["area"]), bool)])
    with open(os.path.join(out_dir, "surface_tawss_osi_rrt.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["x_mm", "y_mm", "z_mm", "area_mm2", "TAWSS_Pa", "OSI", "RRT_perPa", "s_limb_mm", "phi_deg", "region"])
        names = np.full(len(surf["area"]), "", dtype=object)
        for k in ("body", "bifurcation", "limb_ridge", "limb_distal"):
            names[reg[k]] = k
        for i in range(len(surf["area"])):
            w.writerow([f"{surf['c'][i, 0]:.3f}", f"{surf['c'][i, 1]:.3f}", f"{surf['c'][i, 2]:.3f}",
                        f"{surf['area'][i]:.5f}", f"{surf['tawss'][i]:.4f}", f"{surf['osi'][i]:.4f}",
                        f"{surf['rrt'][i]:.3f}", f"{s[i]:.2f}", f"{phi[i]:.1f}", names[i]])

    # ---- helicity and secondary flow over the last-cycle snapshots ----
    snaps = [t for t in times if t > t0 + 1e-6]
    r.disable_all_patch_arrays()
    hel = []
    sec_profiles = []
    for t in snaps:
        r.set_active_time_value(t)
        r.disable_all_cell_arrays()
        for f in ("U", "vorticity"):
            r.enable_cell_array(f)
        vol = r.read()["internalMesh"]
        c = np.asarray(vol.cell_centers().points) * 1e3
        sl, rll, _ = limb_coords(g, c)
        sel = (sl >= RIDGE_S[0]) & (rll < g.R_limb) & (c[:, 0] <= g.x_out - 5.0)
        V = np.asarray(vol.compute_cell_sizes(length=False, area=False).cell_data["Volume"])[sel]
        U = np.asarray(vol.cell_data["U"])[sel]
        Wv = np.asarray(vol.cell_data["vorticity"])[sel]
        h = np.einsum("ij,ij->i", U, Wv)
        nu, nw = np.linalg.norm(U, axis=1), np.linalg.norm(Wv, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            lnh = np.where(nu * nw > 1e-12, np.abs(h) / (nu * nw), 0.0)
        uax = U @ g.d
        usec = np.linalg.norm(U - np.outer(uax, g.d), axis=1)
        Vs = V.sum()
        hel.append([t, (h * V).sum() / Vs, (np.abs(h) * V).sum() / Vs, (lnh * V).sum() / Vs,
                    (usec * V).sum() / max((np.abs(uax) * V).sum(), 1e-30)])
        bins = np.arange(RIDGE_S[0], 130.0 + 1e-9, 5.0)
        idx = np.digitize(sl[sel], bins) - 1
        prof = [(np.abs(h)[idx == k] * V[idx == k]).sum() / max(V[idx == k].sum(), 1e-30)
                for k in range(len(bins) - 1)]
        sec_profiles.append(prof)
    hel = np.array(hel)
    out.update({"limb_h1_signed_m_s2": float(hel[:, 1].mean()), "limb_h2_m_s2": float(hel[:, 2].mean()),
                "limb_LNH": float(hel[:, 3].mean()), "limb_secondary_flow_ratio": float(hel[:, 4].mean())})
    np.savetxt(os.path.join(out_dir, "helicity_timeseries.csv"), hel, delimiter=",",
               header="time_s,h1_signed_m_s2,h2_m_s2,LNH,secondary_flow_ratio", comments="")
    bins = np.arange(RIDGE_S[0], 130.0 + 1e-9, 5.0)
    np.savetxt(os.path.join(out_dir, "h2_along_limb.csv"),
               np.c_[0.5 * (bins[1:] + bins[:-1]), np.mean(sec_profiles, axis=0)], delimiter=",",
               header="s_limb_mm,h2_m_s2", comments="")

    # ---- integral quantities from function objects ----
    pp = os.path.join(cdir, "postProcessing")
    a = _dat(f"{pp}/inletPressure/*/surfaceFieldValue.dat")
    if a is not None:
        m, v, _ = _cycle_mean(a, 1, t0, t_end)
        out["pressure_drop_mean_Pa"] = rho * m
        out["pressure_drop_peak_Pa"] = rho * float(v.max())
    a = _dat(f"{pp}/inletFlow/*/surfaceFieldValue.dat")
    if a is not None:
        m, _, _ = _cycle_mean(a, 1, t0, t_end)
        out["inflow_mean_mL_s_full_model"] = -2 * m * 1e6
    a = _dat(f"{pp}/graftForces/*/force.dat")
    if a is not None:
        out["force_x_total_mN_full"] = 2e3 * _cycle_mean(a, 1, t0, t_end)[0]
        out["force_x_pressure_mN_full"] = 2e3 * _cycle_mean(a, 4, t0, t_end)[0]
        out["force_x_viscous_mN_full"] = 2e3 * _cycle_mean(a, 7, t0, t_end)[0]
        out["force_x_peak_mN_full"] = 2e3 * float(np.max(a[(a[:, 0] >= t0), 1]))
    a = _dat(f"{pp}/wallWSS/*/surfaceFieldValue.dat")
    if a is not None and t_end >= 2 * period - 1e-6:
        ph_grid = np.linspace(0, period, 200, endpoint=False)
        last = np.interp(ph_grid + t0, a[:, 0], a[:, 1])
        prev = np.interp(ph_grid + t0 - period, a[:, 0], a[:, 1])
        out["cycle_to_cycle_WSS_change_pct"] = float(100 * np.abs(last - prev).mean() / np.abs(last).mean())

    _plot_maps(out_dir, name, s, phi, surf, reg, has_ridge)
    _plot_sections(r, g, snaps, t0, period, out_dir, name)
    return out


# --------------------------------------------------------------------------
def _plot_maps(out_dir, name, s, phi, surf, reg, has_ridge):
    sel = reg["limb"]
    fig, ax = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True)
    for a, key, lab, vmax, cmap in ((ax[0], "tawss", "TAWSS (Pa)", 1.5, "viridis"),
                                    (ax[1], "osi", "OSI", 0.3, "magma")):
        sc = a.scatter(s[sel], phi[sel], c=surf[key][sel], s=0.4, vmin=0, vmax=vmax, cmap=cmap, rasterized=True)
        plt.colorbar(sc, ax=a, label=lab)
        a.set_ylabel("angle around limb (°)")
        a.set_ylim(-180, 180)
        a.axvspan(*RIDGE_S, color="w", alpha=0.0, ec="w", ls="--")
    ax[0].set_title(f"{name}: iliac limb wall, unrolled (flow left to right; ridge zone {RIDGE_S[0]:.0f}–{RIDGE_S[1]:.0f} mm)")
    ax[1].set_xlabel("distance from flow divider along the limb (mm)")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "limb_tawss_osi_map.png"), dpi=110)
    plt.close(fig)


def _plot_sections(r, g, snaps, t0, period, out_dir, name):
    """Axial velocity + in-plane (secondary) velocity on two limb cross-sections,
    at peak systole and in late systolic deceleration of the last cycle."""
    d = g.d
    e1 = np.array([-d[1], d[0], 0.0])
    e2 = np.array([0.0, 0.0, 1.0])
    phases = {"peak systole": 0.126, "deceleration": 0.22, "reverse flow": 0.30}
    secs = (50.0, 95.0)
    fig, ax = plt.subplots(len(secs), len(phases), figsize=(14, 9))
    for j, (lab, phs) in enumerate(phases.items()):
        t = min(snaps, key=lambda tt: abs((tt - t0) - phs))
        r.set_active_time_value(t)
        r.disable_all_cell_arrays()
        r.enable_cell_array("U")
        vol = r.read()["internalMesh"]
        for i, sv in enumerate(secs):
            o = (g.A + sv * d) * 1e-3
            sl = vol.slice(normal=d, origin=o)
            if sl.n_cells == 0:
                continue
            c = (np.asarray(sl.cell_centers().points) - o) * 1e3
            U = np.asarray(sl.cell_data["U"])
            a1, a2 = c @ e1, c @ e2
            ua = U @ d
            u1, u2 = U @ e1, U @ e2
            a = ax[i, j]
            sc = a.tricontourf(a1, a2, ua, levels=21, cmap="RdBu_r",
                               vmin=-max(abs(ua).max(), 1e-6), vmax=max(abs(ua).max(), 1e-6))
            grid = np.linspace(-g.R_limb, g.R_limb, 15)
            GX, GY = np.meshgrid(grid, grid)
            from scipy.interpolate import griddata
            G1 = griddata((a1, a2), u1, (GX, GY))
            G2 = griddata((a1, a2), u2, (GX, GY))
            a.quiver(GX, GY, G1, G2, color="k", scale=None)
            plt.colorbar(sc, ax=a, label="axial velocity (m/s)")
            a.set_aspect("equal")
            a.set_title(f"s = {sv:.0f} mm, {lab} (t = {t - t0:.3f} s)", fontsize=9)
            a.set_xlabel("lateral (mm)"); a.set_ylabel("anterior (mm)")
    fig.suptitle(f"{name}: limb cross-sections (arrows = secondary, in-plane flow)")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "limb_sections.png"), dpi=100)
    plt.close(fig)


# --------------------------------------------------------------------------
KEYS = [
    ("limb_TAWSS_Pa", "Limb TAWSS (Pa)", "higher"),
    ("limb_lowWSS_area_pct", "Limb area TAWSS < 0.4 Pa (%)", "lower"),
    ("limb_relLowWSS_area_pct", "Limb area TAWSS < 50% of smooth mean (%)", "lower"),
    ("limb_OSI", "Limb OSI", "lower"),
    ("limb_RRT_perPa", "Limb RRT (1/Pa)", "lower"),
    ("bifurcation_TAWSS_Pa", "Bifurcation TAWSS (Pa)", "higher"),
    ("limb_h2_m_s2", "Limb helicity intensity h2 (m/s²)", "higher"),
    ("pressure_drop_mean_Pa", "Pressure drop, cycle mean (Pa)", "lower"),
    ("force_x_total_mN_full", "Axial displacement force, flow part (mN)", "lower"),
]


def summarise(res_root: str, order: list[str]):
    rows = []
    for n in order:
        p = os.path.join(res_root, n, "metrics.json")
        if os.path.exists(p):
            rows.append(json.load(open(p)))
    if not rows:
        print("no finished designs")
        return
    names = [r["design"] for r in rows]
    ref = next((r for r in rows if r["design"] == "smooth"), None)
    if ref is not None:          # relative low-shear threshold from the smooth graft
        thr = 0.5 * ref["limb_TAWSS_Pa"]
        for r in rows:
            f = np.load(os.path.join(res_root, r["design"], "wall_faces.npz"))
            import geometry_evar as ge
            reg = regions(ge.EvarGeom(), f["c"])["limb"]
            w = f["area"][reg] / f["area"][reg].sum()
            r["limb_relLowWSS_area_pct"] = float(100 * w[f["tawss"][reg] < thr].sum())
            r["relLowWSS_threshold_Pa"] = thr
    keys = sorted({k for r in rows for k in r})
    with open(os.path.join(res_root, "metrics.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["design"] + [k for k in keys if k != "design"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    if ref is not None:
        with open(os.path.join(res_root, "vs_smooth.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["metric", "better_when"] + [f"{n}" for n in names] + [f"{n}_change_pct" for n in names])
            for k, lab, better in KEYS:
                vals = [r.get(k) for r in rows]
                ch = [None if (v is None or ref.get(k) in (None, 0)) else 100 * (v - ref[k]) / abs(ref[k]) for v in vals]
                w.writerow([lab, better] + [f"{v:.4g}" if v is not None else "" for v in vals]
                           + [f"{c:+.1f}" if c is not None else "" for c in ch])
    fig, axs = plt.subplots(3, 3, figsize=(15, 11))
    for a, (k, lab, better) in zip(axs.flat, KEYS):
        vals = [r.get(k, np.nan) or np.nan for r in rows]
        cols = ["0.5" if n == "smooth" else "C0" for n in names]
        a.bar(range(len(names)), vals, color=cols)
        a.set_xticks(range(len(names)))
        a.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
        a.set_title(f"{lab}\n({better} is better)", fontsize=9)
        if ref is not None and ref.get(k) is not None:
            a.axhline(ref[k], color="k", lw=0.8, ls="--")
    plt.tight_layout()
    plt.savefig(os.path.join(res_root, "comparison.png"), dpi=100)
    plt.close(fig)
    print(f"summary written to {res_root}")
