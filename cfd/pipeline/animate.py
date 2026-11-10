"""Animated GIFs + data export for a finished pulsatile case.

For every snapshot of the last cardiac cycle it records:
  * inlet flow (whole vessel, mL/min) and pressure drop
  * wall shear stress statistics on the artery wall under the stent and on the stent
  * reversed-flow (recirculation) volume
and draws one frame:
  (a) inlet waveform with the current phase
  (b) unrolled wall-shear-stress map of the artery wall
  (c) longitudinal slice through the stent with velocity magnitude + streamlines
"""
from __future__ import annotations

import csv
import io
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402
from PIL import Image  # noqa: E402

LOW_WSS_PA = 0.4


def _reader(cdir):
    foam = os.path.join(cdir, "case.foam")
    open(foam, "a").close()
    r = pv.OpenFOAMReader(foam)
    r.enable_all_patch_arrays()
    return r


def snapshot_times(cdir, period):
    r = _reader(cdir)
    ts = [t for t in r.time_values if t > 0]
    if not ts:
        return []
    last_cycle_start = (math.ceil(ts[-1] / period - 1e-6) - 1) * period
    return [t for t in ts if t > last_cycle_start + 1e-9]


def _areas(mesh):
    return np.asarray(mesh.compute_cell_sizes(length=False, volume=False).cell_data["Area"])


def save_gif(frames: list[Image.Image], path: str, fps: float = 6):
    pal = [f.convert("RGB").quantize(colors=200, method=Image.Quantize.MEDIANCUT) for f in frames]
    pal[0].save(path, save_all=True, append_images=pal[1:], duration=int(1000 / fps), loop=0, optimize=True)


def _fig_to_image(fig) -> Image.Image:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=80)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).copy()


def export_flow(cdir: str, out_dir: str, dom, rho: float, period: float, title: str,
                slice_theta_deg: float = 30.0, fps: float = 6, timeseries_name: str = "timeseries.csv",
                make_gif: bool = True) -> dict:
    """Write <timeseries_name>, surface_tawss_osi_rrt.csv and (if make_gif) flow.gif into out_dir."""
    os.makedirs(out_dir, exist_ok=True)
    r = _reader(cdir)
    times = snapshot_times(cdir, period)
    if not times:
        raise RuntimeError(f"no snapshots in {cdir}; run with --keep-snapshots")
    x0, x1, R = dom.x0 * 1e-3, dom.x1 * 1e-3, dom.radius * 1e-3
    th = math.radians(slice_theta_deg)
    normal = (0.0, -math.sin(th), math.cos(th))

    rows, cache = [], []
    for t in times:
        r.set_active_time_value(t)
        d = r.read()
        vol, bnd = d["internalMesh"], d["boundary"]
        wall, stent, inlet = bnd["arteryWall"], bnd["stent"], bnd["inlet"]
        wc = np.asarray(wall.cell_centers().points)
        wsel = (wc[:, 0] >= x0) & (wc[:, 0] <= x1)
        wa = _areas(wall)
        wwss = rho * np.linalg.norm(np.asarray(wall.cell_data["wallShearStress"]), axis=1)
        sa = _areas(stent)
        swss = rho * np.linalg.norm(np.asarray(stent.cell_data["wallShearStress"]), axis=1)
        ia = _areas(inlet)
        q_sector = float((np.asarray(inlet.cell_data["U"])[:, 0] * ia).sum())
        p_in = float((np.asarray(inlet.cell_data["p"]) * ia).sum() / ia.sum()) * rho
        U = np.asarray(vol.cell_data["U"])
        cv = np.asarray(vol.compute_cell_sizes(length=False, area=False).cell_data["Volume"])
        xc = np.asarray(vol.cell_centers().points)[:, 0]
        reg = (xc >= x0) & (xc <= x1 + 1e-3)
        rev = reg & (U[:, 0] < 0)
        rows.append({
            "time_s": t, "phase": (t % period) / period,
            "inlet_flow_mL_min_full_vessel": q_sector * 360.0 / (dom.theta1_deg - dom.theta0_deg) * 6e7,
            "pressure_drop_Pa": p_in,
            "wall_WSS_mean_Pa": float((wwss[wsel] * wa[wsel]).sum() / wa[wsel].sum()),
            "wall_lowWSS_area_pct": float(100 * wa[wsel][wwss[wsel] < LOW_WSS_PA].sum() / wa[wsel].sum()),
            "stent_WSS_mean_Pa": float((swss * sa).sum() / sa.sum()),
            "stent_lowWSS_area_pct": float(100 * sa[swss < LOW_WSS_PA].sum() / sa.sum()),
            "recirc_volume_mm3": float(cv[rev].sum() * 1e9),
            "max_velocity_m_s": float(np.linalg.norm(U, axis=1).max()),
        })
        if not make_gif:
            continue
        sl = vol.slice(normal=normal, origin=(0, 0, 0)).triangulate()
        sp = np.asarray(sl.points)
        rr = sp[:, 1] * math.cos(th) + sp[:, 2] * math.sin(th)
        Up = np.asarray(sl.point_data["U"])
        ur = Up[:, 1] * math.cos(th) + Up[:, 2] * math.sin(th)
        tris = sl.faces.reshape(-1, 4)[:, 1:]
        keep = rr[tris].min(1) > 0
        cache.append(dict(wc=wc, wwss=wwss, sx=sp[:, 0], sr=rr, ux=Up[:, 0], ur=ur, tris=tris[keep]))

    # ---- data -----------------------------------------------------------
    with open(os.path.join(out_dir, timeseries_name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    _export_surface(r, times[-1], out_dir, rho)
    if not make_gif:
        return {"n_frames": 0, "times": [float(t) for t in times]}

    # ---- frames ---------------------------------------------------------
    tt = np.array([q["time_s"] for q in rows])
    qq = np.array([q["inlet_flow_mL_min_full_vessel"] for q in rows])
    order = np.argsort(tt % period)
    xg = np.linspace(dom.x0 - 0.2, dom.x1 + 0.7, 320) * 1e-3
    rg = np.linspace(1.38, dom.radius, 70) * 1e-3
    XG, RG = np.meshgrid(xg, rg)
    frames = []
    for k, (t, c) in enumerate(zip(times, cache)):
        fig = plt.figure(figsize=(13, 9))
        gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[1.35, 1], hspace=0.32, wspace=0.18)
        ax0 = fig.add_subplot(gs[0, 0])
        ax0.plot((tt % period)[order], qq[order], "k-", lw=1.5)
        ax0.plot([t % period], [rows[k]["inlet_flow_mL_min_full_vessel"]], "o", color="crimson", ms=10)
        ax0.set_xlim(0, period)
        ax0.set_xlabel("time in cardiac cycle (s)")
        ax0.set_ylabel("inflow (mL/min)")
        ax0.set_title(f"{title}\nt = {t % period:.3f} s", fontsize=12, loc="left")
        rk = rows[k]
        ax0.text(0.02, 0.97, f"wall WSS under stent: {rk['wall_WSS_mean_Pa']:.2f} Pa\n"
                 f"wall area < 0.4 Pa: {rk['wall_lowWSS_area_pct']:.0f} %\n"
                 f"pressure drop: {rk['pressure_drop_Pa']:.1f} Pa",
                 transform=ax0.transAxes, va="top", fontsize=9, family="monospace",
                 bbox=dict(fc="white", ec="0.7"))

        ax1 = fig.add_subplot(gs[0, 1])
        xs = c["wc"][:, 0] * 1e3
        ss = dom.radius * np.arctan2(c["wc"][:, 2], c["wc"][:, 1])
        m = (xs > dom.x0 - 0.3) & (xs < dom.x1 + 1.2)
        sc = ax1.scatter(xs[m], ss[m], c=c["wwss"][m], s=1.5, cmap="turbo", vmin=0, vmax=3.0, rasterized=True)
        ax1.axhline(dom.radius * th, color="w", ls="--", lw=1)
        ax1.set_aspect("equal")
        ax1.set_xlim(dom.x0 - 0.3, dom.x1 + 1.2)
        ax1.set_xlabel("axial (mm)  → flow")
        ax1.set_ylabel("circumferential (mm)")
        ax1.set_title("WSS on artery wall (Pa). White = under struts;\ndashed = slice below", fontsize=10, loc="left")
        fig.colorbar(sc, ax=ax1, fraction=0.04, pad=0.02)

        ax2 = fig.add_subplot(gs[1, :])
        speed = np.hypot(c["ux"], c["ur"])
        tp = ax2.tripcolor(c["sx"] * 1e3, c["sr"] * 1e3, c["tris"], speed, cmap="magma",
                           vmin=0, vmax=0.08, shading="gouraud", rasterized=True)
        UX, UR = _slice_velocity_on_grid(c, XG, RG)
        try:
            ax2.streamplot(XG * 1e3, RG * 1e3, UX, UR, density=(4.0, 1.6), color="cyan",
                           linewidth=0.45, arrowsize=0.5)
        except ValueError:
            pass
        ax2.set_xlim(xg[0] * 1e3, xg[-1] * 1e3)
        ax2.set_ylim(rg[0] * 1e3, rg[-1] * 1e3)
        ax2.set_xlabel("axial position (mm)  → flow")
        ax2.set_ylabel("radius (mm)")
        ax2.set_title(f"near-wall velocity in a longitudinal slice (θ = {slice_theta_deg:.0f}°); artery wall at top, "
                      "white = strut cross-sections, radial scale stretched ×4", fontsize=10, loc="left")
        fig.colorbar(tp, ax=ax2, fraction=0.025, pad=0.01, label="|u| (m/s)")
        frames.append(_fig_to_image(fig))
    save_gif(frames, os.path.join(out_dir, "flow.gif"), fps=fps)
    frames[int(np.argmax(qq))].save(os.path.join(out_dir, "frame_peak_flow.png"))
    return {"n_frames": len(frames), "times": [float(t) for t in times]}


def _slice_velocity_on_grid(c, XG, RG, bridge_factor: float = 2.5):
    """Interpolate the slice velocity onto the streamline grid. Uses the slice's own triangles
    (struts stay holes). Some meshes (e.g. dense denticles) give overlapping triangles that
    matplotlib rejects; then re-triangulate with Delaunay and mask the triangles that bridge a
    hole (longest edge > bridge_factor x the local cell size of their corners)."""
    try:
        tri = mtri.Triangulation(c["sx"], c["sr"], c["tris"])
        UX = mtri.LinearTriInterpolator(tri, c["ux"])(XG, RG).filled(np.nan)
        UR = mtri.LinearTriInterpolator(tri, c["ur"])(XG, RG).filled(np.nan)
        return UX, UR
    except RuntimeError:
        pass
    x, r, t0 = c["sx"], c["sr"], c["tris"]

    def edges(t):
        e = [np.hypot(x[t[:, i]] - x[t[:, j]], r[t[:, i]] - r[t[:, j]]) for i, j in ((0, 1), (1, 2), (2, 0))]
        return np.stack(e, 1)
    e0 = edges(t0)
    local = np.zeros(len(x))
    np.maximum.at(local, t0.ravel(), np.repeat(e0.max(1), 3))
    d = mtri.Triangulation(x, r)
    td = d.triangles
    d.set_mask(edges(td).max(1) > bridge_factor * local[td].max(1))
    UX = mtri.LinearTriInterpolator(d, c["ux"])(XG, RG).filled(np.nan)
    UR = mtri.LinearTriInterpolator(d, c["ur"])(XG, RG).filled(np.nan)
    return UX, UR


def _export_surface(r, t, out_dir, rho):
    """Per-face TAWSS / OSI / RRT over the last cycle on artery wall and stent."""
    r.set_active_time_value(t)
    d = r.read()
    with open(os.path.join(out_dir, "surface_tawss_osi_rrt.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["surface", "x_mm", "y_mm", "z_mm", "area_um2", "TAWSS_Pa", "OSI", "RRT_perPa"])
        for name in ("arteryWall", "stent"):
            p = d["boundary"][name]
            c = np.asarray(p.cell_centers().points) * 1e3
            a = _areas(p) * 1e12
            taw = rho * np.asarray(p.cell_data["magWSSMean"])
            mv = rho * np.linalg.norm(np.asarray(p.cell_data["wallShearStressMean"]), axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                osi = np.clip(np.where(taw > 0, 0.5 * (1 - mv / taw), 0), 0, 0.5)
                rrt = np.where(mv > 1e-9, 1 / mv, np.nan)
            for i in range(len(a)):
                w.writerow([name, f"{c[i, 0]:.5f}", f"{c[i, 1]:.5f}", f"{c[i, 2]:.5f}", f"{a[i]:.2f}",
                            f"{taw[i]:.5f}", f"{osi[i]:.5f}", f"{rrt[i]:.4f}"])
