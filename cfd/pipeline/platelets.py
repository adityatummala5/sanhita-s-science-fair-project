"""Lagrangian platelet transport, adhesion and aggregation on a computed blood-flow field.

Model (one-way coupled: platelets are moved by the flow, they do not change it)
------------------------------------------------------------------------------
Transport  dx = u(x,t) dt + sqrt(2 D dt) * N(0,1)
           u: velocity from the OpenFOAM pulsatile Carreau solution (last cycle, periodic)
           D = D_B + alpha * shear_rate      shear-enhanced diffusivity from red-cell
               collisions (Zydney & Colton 1988), alpha = 0.15 a_rbc^2 Hct (1-Hct)^0.8,
               D_B = 1.58e-13 m2/s Brownian (Sorensen et al. 1999)
States     resting -> activated -> adhered (to the stent) / aggregated (to an adhered platelet)
Adhesion   first-order surface reaction, flux J = k C (Sorensen et al. 1999, Ann Biomed Eng 27:436).
           For a particle inside a reaction layer of thickness delta next to a reactive
           surface this gives a probability per step  P = 1 - exp(-k dt / delta).
           k_rs (resting -> surface), k_as (activated -> surface),
           k_aa (activated -> already-deposited platelet = aggregation).
Activation deposited platelets release agonists (ADP, TxA2). Sorensen solves agonist
           transport; here this is SIMPLIFIED to: a free platelet passing within
           activation_radius of a deposited platelet becomes activated.
Surfaces   stent struts are reactive (foreign surface). The artery wall is treated
           as intact endothelium (non-adhesive) unless artery_wall_reactive = True.
Parcels    each computational particle represents `weight` real platelets so the true
           platelet flux (2.5e14 /m3) can be simulated; reported counts are in real platelets.

Main simplifications (state them): no red cells (margination only via near-wall seeding),
no agonist transport field, no detachment/embolisation, no thrombus feedback on the flow,
no shear-induced (mechanical) activation, real time scale of ~2 s (thrombus growth takes minutes).
"""
from __future__ import annotations

import csv
import io
import json
import math
import os
import time
from dataclasses import asdict, dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402
import trimesh  # noqa: E402
from PIL import Image  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

try:
    import vtk  # noqa: E402
    vtk.vtkLogger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_OFF)
except Exception:  # pragma: no cover
    pass

RESTING, ACTIVATED = 0, 1


@dataclass
class PlateletParams:
    concentration_per_uL: float = 250_000.0    # normal count 150-450 x10^3 /uL
    diameter_um: float = 3.0                   # platelet diameter 2-3 um
    D_brownian_m2_s: float = 1.58e-13          # Sorensen 1999
    hematocrit: float = 0.40
    rbc_radius_um: float = 4.2                 # Zydney & Colton 1988
    k_rs_m_s: float = 3.7e-6                   # resting platelet -> surface (Sorensen 1999)
    k_as_m_s: float = 4.6e-5                   # activated platelet -> surface (Sorensen 1999)
    k_aa_m_s: float = 4.6e-5                   # activated platelet -> deposited platelet (Sorensen 1999)
    reaction_layer_um: float = 3.0             # delta: ~1 platelet diameter
    activation_radius_um: float = 10.0         # simplification of agonist release (see docstring)
    p_inf_per_m2: float = 7.0e10               # max surface coverage (Sorensen 1999)
    artery_wall_reactive: bool = False
    seed_r_min_mm: float = 1.40                # >= 50 um beyond the tallest strut/denticle tip (r = 1.457)
    track_r_min_mm: float = 1.30
    seed_offset_mm: float = 0.5                # seed plane this far upstream of the stent
    track_end_mm: float = 1.5                  # stop tracking this far downstream of the stent
    target_parcels: int = 150_000              # sets the parcel weight
    n_cycles: float = 2.0
    dt_s: float = 1e-4
    frame_interval_s: float = 0.02
    seed: int = 1


# ---------------------------------------------------------------------------
class FlowField:
    """Periodic, time-interpolated velocity and shear rate, interpolated from the
    OpenFOAM mesh points near the wall by inverse-distance weighting of the 8
    nearest points (wall points carry u = 0). Units: mm, mm/s, 1/s."""

    def __init__(self, cdir, times, period, dom, p: PlateletParams):
        foam = os.path.join(cdir, "case.foam")
        open(foam, "a").close()
        r = pv.OpenFOAMReader(foam)
        self.period = period
        self.phase = np.array([(t % period) if (t % period) < period - 1e-6 else 0.0 for t in times])
        order = np.argsort(self.phase)
        self.phase = self.phase[order]
        times = [times[i] for i in order]
        r.set_active_time_value(times[0])
        vol = r.read()["internalMesh"]
        P = np.asarray(vol.points) * 1e3
        rad = np.hypot(P[:, 1], P[:, 2])
        x_lo = dom.x0 - p.seed_offset_mm - 0.15
        x_hi = dom.x1 + p.track_end_mm + 0.15
        self.mask = (P[:, 0] > x_lo) & (P[:, 0] < x_hi) & (rad > p.track_r_min_mm - 0.08)
        self.pts = P[self.mask]
        self.tree = cKDTree(self.pts)
        print(f"    platelet flow mesh: {len(self.pts):,} points")
        self.U, self.G = [], []
        ids = np.where(self.mask)[0]
        for t in times:
            r.set_active_time_value(t)
            v = r.read()["internalMesh"]
            sub = v.extract_points(ids, adjacent_cells=False, include_cells=True)
            import logging
            logging.disable(logging.CRITICAL)        # VTK warns about a few degenerate cells
            sub = sub.compute_derivative(scalars="U", gradient="gradU")
            logging.disable(logging.NOTSET)
            # map extracted points back to our ordering
            orig = np.asarray(sub.point_data["vtkOriginalPointIds"])
            Gs = np.asarray(sub.point_data["gradU"]).reshape(-1, 3, 3)
            S = 0.5 * (Gs + Gs.transpose(0, 2, 1))
            gam = np.nan_to_num(np.sqrt(2.0 * np.einsum("nij,nij->n", S, S)))
            g_full = np.zeros(len(P))
            g_full[orig] = gam
            self.U.append(np.asarray(v.point_data["U"])[self.mask] * 1e3)
            self.G.append(g_full[self.mask])
        self.U = np.array(self.U)
        self.G = np.array(self.G)

    def _weights(self, t):
        ph = t % self.period
        k = np.searchsorted(self.phase, ph) - 1
        k0 = k % len(self.phase)
        k1 = (k + 1) % len(self.phase)
        p0, p1 = self.phase[k0], self.phase[k1]
        span = (p1 - p0) % self.period or self.period
        a = ((ph - p0) % self.period) / span
        return k0, k1, a

    def sample(self, x, t, k=8):
        d, vi = self.tree.query(x, k=k)
        ok = d[:, 0] < 0.05                       # inside the tracked region
        w = 1.0 / np.maximum(d, 1e-5) ** 2
        w /= w.sum(1, keepdims=True)
        k0, k1, a = self._weights(t)
        Uk = (1 - a) * self.U[k0][vi] + a * self.U[k1][vi]
        Gk = (1 - a) * self.G[k0][vi] + a * self.G[k1][vi]
        u = np.einsum("ni,nij->nj", w, Uk)
        g = np.einsum("ni,ni->n", w, Gk)
        u[~ok] = 0.0
        return u, g, ok


class Surface:
    """Dense point sampling of the stent surface with outward normals:
    signed distance = (x - q) . n  (negative = inside a strut)."""

    def __init__(self, stl, lumen_radius_mm, spacing_um=2.5):
        m = trimesh.load(stl)
        n = int(m.area / (spacing_um * 1e-3) ** 2)
        pts, fid = trimesh.sample.sample_surface_even(m, n, seed=0)
        pts = np.asarray(pts)
        wet = np.hypot(pts[:, 1], pts[:, 2]) < lumen_radius_mm   # only faces in contact with blood
        self.wetted_area = float(m.area * wet.mean())
        self.q = pts[wet]
        self.n = m.face_normals[fid][wet]
        self.tree = cKDTree(self.q)
        radial = (self.n[:, 1] * self.q[:, 1] + self.n[:, 2] * self.q[:, 2]) / np.hypot(self.q[:, 1], self.q[:, 2])
        self.kind = np.where(radial < -0.7, 1, 2)    # 1 = luminal (faces the blood stream), 2 = strut side

    def signed(self, x):
        d, i = self.tree.query(x, distance_upper_bound=0.05)
        far = ~np.isfinite(d)
        i = np.where(far, 0, i)
        sd = np.einsum("ni,ni->n", x - self.q[i], self.n[i])
        sd[far] = 1.0
        return sd, i


# ---------------------------------------------------------------------------
def simulate(cdir, times, period, dom, stent_stl, out_dir, p: PlateletParams, title="bare stent",
             timeseries_name: str = "timeseries.csv"):
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(p.seed)
    ff = FlowField(cdir, times, period, dom, p)
    surf = Surface(stent_stl, dom.radius)
    R = dom.radius
    th0, th1 = math.radians(dom.theta0_deg), math.radians(dom.theta1_deg)
    alpha = 0.15 * (p.rbc_radius_um * 1e-3) ** 2 * p.hematocrit * (1 - p.hematocrit) ** 0.8   # mm2
    D_B = p.D_brownian_m2_s * 1e6                                                              # mm2/s
    conc = p.concentration_per_uL                                                              # per mm3
    delta = p.reaction_layer_um * 1e-3
    d_p = p.diameter_um * 1e-3
    r_act = p.activation_radius_um * 1e-3
    pr = {k: 1 - math.exp(-getattr(p, k) * 1e3 * p.dt_s / delta) for k in ("k_rs_m_s", "k_as_m_s", "k_aa_m_s")}
    x_seed = dom.x0 - p.seed_offset_mm
    x_end = dom.x1 + p.track_end_mm

    # seed table per snapshot: candidate points on the seed plane annulus, weighted by axial velocity
    nc = 40000
    rs = np.sqrt(rng.uniform(p.seed_r_min_mm ** 2, R ** 2, nc))
    ts = rng.uniform(th0, th1, nc)
    cand = np.c_[np.full(nc, x_seed), rs * np.cos(ts), rs * np.sin(ts)]
    area = 0.5 * (th1 - th0) * (R ** 2 - p.seed_r_min_mm ** 2)
    seed_u = []
    for k in range(len(ff.phase)):
        u, _, ok = ff.sample(cand, ff.phase[k])
        seed_u.append(np.where(ok, np.clip(u[:, 0], 0, None), 0.0))
    seed_u = np.array(seed_u)
    q_seed = seed_u.mean(1) * area                            # mm3/s per snapshot
    total_platelets = conc * q_seed.mean() * p.n_cycles * period
    weight = max(1.0, total_platelets / p.target_parcels)
    print(f"    platelets reaching the seed annulus: {conc * q_seed.mean():.3g}/s; "
          f"parcel weight = {weight:.1f} platelets")

    pos = np.zeros((0, 3))
    state = np.zeros(0, int)
    pid = np.zeros(0, int)
    ads_pos, ads_kind, ads_rec = [], [], []
    next_id = 0
    t, t_end = 0.0, p.n_cycles * period
    frames, series = [], []
    next_frame = 0.0
    carry = 0.0
    counts = dict(injected=0, removed=0, activated=0)
    ads_tree, ads_dirty = None, False
    t_start = time.time()
    while t < t_end:
        # --- inject
        k0, k1, a = ff._weights(t)
        qs = (1 - a) * q_seed[k0] + a * q_seed[k1]
        carry += conc * qs * p.dt_s / weight
        n_new = int(carry)
        carry -= n_new
        if n_new:
            wts = (1 - a) * seed_u[k0] + a * seed_u[k1]
            pick = rng.choice(nc, n_new, p=wts / wts.sum())
            pos = np.vstack([pos, cand[pick]])
            state = np.r_[state, np.zeros(n_new, int)]
            pid = np.r_[pid, np.arange(next_id, next_id + n_new)]
            next_id += n_new
            counts["injected"] += n_new
        if len(pos):
            u, g, ok = ff.sample(pos, t)
            D = D_B + alpha * g
            new = pos + u * p.dt_s + np.sqrt(2 * D * p.dt_s)[:, None] * rng.standard_normal(pos.shape)
            # artery wall (reflect) and symmetry planes (reflect)
            rr = np.hypot(new[:, 1], new[:, 2])
            th = np.arctan2(new[:, 2], new[:, 1])
            rr = np.where(rr > R, 2 * R - rr, rr)
            th = np.where(th < th0, 2 * th0 - th, th)
            th = np.where(th > th1, 2 * th1 - th, th)
            new = np.c_[new[:, 0], rr * np.cos(th), rr * np.sin(th)]
            # struts: reject moves into the solid
            sd, si = surf.signed(new)
            inside = sd < 0
            new[inside] = pos[inside]
            sd[inside], si[inside] = surf.signed(new[inside]) if inside.any() else (sd[inside], si[inside])
            pos = new
            # --- adhesion to the stent
            near = sd < delta
            prob = np.where(state == ACTIVATED, pr["k_as_m_s"], pr["k_rs_m_s"])
            stick = near & (rng.random(len(pos)) < prob)
            if p.artery_wall_reactive:
                nearw = (R - rr) < delta
                stick |= nearw & (rng.random(len(pos)) < prob)
            # --- aggregation onto deposited platelets + activation by them
            agg = np.zeros(len(pos), bool)
            if ads_pos:
                if ads_dirty or ads_tree is None:
                    ads_tree = cKDTree(np.array(ads_pos))
                    ads_dirty = False
                dist, j = ads_tree.query(pos, distance_upper_bound=r_act)
                close = np.isfinite(dist)
                newly = close & (state == RESTING)
                counts["activated"] += int(newly.sum())
                state[newly] = ACTIVATED
                contact = (dist < d_p + delta) & (state == ACTIVATED) & ~stick
                agg = contact & (rng.random(len(pos)) < pr["k_aa_m_s"])
            for idx in np.where(stick | agg)[0]:
                if stick[idx]:
                    qi = si[idx]
                    xp = surf.q[qi] + surf.n[qi] * (d_p / 2)
                    kind = "stent_luminal" if surf.kind[qi] == 1 else "stent_side"
                else:
                    _, jj = ads_tree.query(pos[idx])
                    base = np.array(ads_pos[jj])
                    v = pos[idx] - base
                    xp = base + v / (np.linalg.norm(v) + 1e-12) * d_p
                    kind = "aggregate"
                ads_pos.append(xp)
                ads_kind.append(kind)
                ads_rec.append((int(pid[idx]), t, *xp, kind, int(state[idx] == ACTIVATED)))
                ads_dirty = True
            gone = stick | agg | ~ok | (pos[:, 0] > x_end) | (np.hypot(pos[:, 1], pos[:, 2]) < p.track_r_min_mm)
            counts["removed"] += int((gone & ~(stick | agg)).sum())
            pos, state, pid = pos[~gone], state[~gone], pid[~gone]
        t += p.dt_s
        if t >= next_frame:
            n_surf = sum(k != "aggregate" for k in ads_kind)
            n_agg = len(ads_kind) - n_surf
            series.append({"time_s": round(t, 5), "phase": (t % period) / period,
                           "free_parcels": len(pos), "free_activated_parcels": int((state == ACTIVATED).sum()),
                           "adhered_platelets": n_surf * weight, "aggregated_platelets": n_agg * weight,
                           "injected_platelets": counts["injected"] * weight})
            sub = rng.choice(len(pos), min(len(pos), 6000), replace=False) if len(pos) else np.zeros(0, int)
            frames.append((t, pos[sub].copy(), state[sub].copy(), np.array(ads_pos), list(ads_kind)))
            next_frame += p.frame_interval_s
            if len(frames) % 10 == 0:
                print(f"      t = {t:.2f} s  free parcels {len(pos):6d}  adhered {n_surf:5d}  aggregated {n_agg:5d}"
                      f"  ({time.time() - t_start:.0f} s)")

    # ---- data -----------------------------------------------------------
    with open(os.path.join(out_dir, "adhered_platelets.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["parcel_id", "platelets_represented", "time_s", "x_mm", "y_mm", "z_mm", "site", "was_activated"])
        for rec in ads_rec:
            w.writerow([rec[0], f"{weight:.2f}", f"{rec[1]:.5f}", f"{rec[2]:.5f}", f"{rec[3]:.5f}",
                        f"{rec[4]:.5f}", rec[5], rec[6]])
    with open(os.path.join(out_dir, timeseries_name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(series[0]))
        w.writeheader()
        w.writerows(series)
    lum_area = float((surf.kind == 1).mean() * surf.wetted_area)
    kinds = np.array(ads_kind)
    xa = np.array([r[2] for r in ads_rec]) if ads_rec else np.zeros(0)
    ring_split = 0.5 * (dom.x0 + dom.x1)
    summary = {
        "title": title, "parameters": asdict(p), "parcel_weight_platelets": weight,
        "simulated_time_s": t_end, "flow_case": cdir,
        "platelets_injected": counts["injected"] * weight,
        "platelets_adhered_to_stent": float((kinds != "aggregate").sum() * weight),
        "platelets_aggregated": float((kinds == "aggregate").sum() * weight),
        "platelets_adhered_luminal_face": float((kinds == "stent_luminal").sum() * weight),
        "platelets_adhered_strut_sides": float((kinds == "stent_side").sum() * weight),
        "stent_wetted_area_mm2": surf.wetted_area,
        "stent_luminal_area_mm2": lum_area,
        "mean_deposit_density_per_mm2": float(len(kinds) * weight / surf.wetted_area),
        "platelets_on_upstream_ring": float((xa < ring_split).sum() * weight),
        "platelets_on_downstream_ring": float((xa >= ring_split).sum() * weight),
        "p_inf_per_mm2": p.p_inf_per_m2 * 1e-6,
        "runtime_min": (time.time() - t_start) / 60,
    }
    json.dump(summary, open(os.path.join(out_dir, "summary.json"), "w"), indent=2)
    np.savez_compressed(os.path.join(out_dir, "frames.npz"),
                        t=np.array([f[0] for f in frames]),
                        free_pos=np.array([f[1] for f in frames], dtype=object),
                        free_state=np.array([f[2] for f in frames], dtype=object),
                        series=json.dumps(series))
    _gif(frames, series, weight, dom, stent_stl, period, os.path.join(out_dir, "platelets.gif"), title)
    _deposit_map(ads_rec, weight, dom, stent_stl, os.path.join(out_dir, "deposition_map.png"), t_end)
    return summary


# ---------------------------------------------------------------------------
def _stent_outline(stl, dom):
    """Unrolled strut luminal faces (mm) for drawing."""
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from geometry import _luminal_polygon
    m = trimesh.load(stl)
    poly = _luminal_polygon(m)
    if poly.is_empty:    # thickened struts: luminal face sits inside R_LUMINAL
        fc = m.triangles_center
        rc = np.hypot(fc[:, 1], fc[:, 2])
        radial = (m.face_normals[:, 1] * fc[:, 1] + m.face_normals[:, 2] * fc[:, 2]) / rc
        poly = _luminal_polygon(m, float(rc[radial < -0.98].min()))
    geoms = [poly] if poly.geom_type == "Polygon" else list(poly.geoms)
    return [np.asarray(g.exterior.coords) for g in geoms]


def _unroll(x):
    return x[:, 0], 1.5 * np.arctan2(x[:, 2], x[:, 1])


def _gif(frames, series, weight, dom, stl, period, path, title):
    outl = _stent_outline(stl, dom)
    ts = np.array([s["time_s"] for s in series])
    adh = np.array([s["adhered_platelets"] for s in series])
    agg = np.array([s["aggregated_platelets"] for s in series])
    zoom = _busiest_window(frames[-1][3], dom)
    imgs = []
    step = max(1, len(frames) // 90)
    for fi in list(range(0, len(frames), step)) + [len(frames) - 1]:
        t, P, S, A, K = frames[fi]
        fig = plt.figure(figsize=(13, 8.5))
        gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1], height_ratios=[3, 1.1], hspace=0.35, wspace=0.15)
        for ax, lim, ms in ((fig.add_subplot(gs[0, 0]), None, 1.0), (fig.add_subplot(gs[0, 1]), zoom, 6.0)):
            for o in outl:
                ax.fill(o[:, 0], o[:, 1], color="0.82", lw=0)
            near = np.hypot(P[:, 1], P[:, 2]) > 1.45 if len(P) else np.zeros(0, bool)
            if near.any():
                x, s = _unroll(P[near])
                st = S[near]
                ax.scatter(x[st == RESTING], s[st == RESTING], s=ms * 0.6, c="tab:blue", alpha=0.5, lw=0,
                           label="free, resting")
                ax.scatter(x[st == ACTIVATED], s[st == ACTIVATED], s=ms * 1.2, c="orange", lw=0,
                           label="free, activated")
            if len(A):
                x, s = _unroll(A)
                k = np.array(K)
                ax.scatter(x[k != "aggregate"], s[k != "aggregate"], s=ms * 4, c="red", lw=0,
                           label="adhered to stent")
                ax.scatter(x[k == "aggregate"], s[k == "aggregate"], s=ms * 4, c="darkred", marker="s", lw=0,
                           label="aggregated")
            ax.set_aspect("equal")
            if lim:
                ax.set_xlim(lim[0], lim[1])
                ax.set_ylim(lim[2], lim[3])
                ax.set_title("zoom: region with the most deposition", fontsize=10)
            else:
                ax.set_xlim(dom.x0 - 0.3, dom.x1 + 0.5)
                ax.set_ylim(1.5 * math.radians(dom.theta0_deg), 1.5 * math.radians(dom.theta1_deg))
                ax.set_title(f"{title}: platelets near the stent (view from the vessel axis, unrolled)   "
                             f"t = {t:.2f} s", fontsize=10, loc="left")
                ax.plot([zoom[0], zoom[1], zoom[1], zoom[0], zoom[0]],
                        [zoom[2], zoom[2], zoom[3], zoom[3], zoom[2]], "k--", lw=0.8)
                ax.legend(loc="upper right", fontsize=7, markerscale=3, framealpha=0.9)
            ax.set_xlabel("axial (mm)  → flow", fontsize=8)
        ax = fig.add_subplot(gs[1, :])
        m = ts <= t + 1e-9
        ax.plot(ts[m], adh[m], "r-", label="adhered to stent")
        ax.plot(ts[m], agg[m], color="darkred", ls="--", label="aggregated onto platelets")
        ax.set_xlim(0, ts[-1])
        ax.set_ylim(0, max(1.0, (adh + agg).max()) * 1.05)
        ax.set_xlabel("time (s)")
        ax.set_ylabel("platelets (real count)")
        ax.legend(loc="upper left", fontsize=8)
        for c in range(1, int(ts[-1] / period) + 1):
            ax.axvline(c * period, color="0.7", lw=0.8)
        ax.set_title(f"cumulative deposition (each marker = {weight:.0f} platelets); grey lines = end of a heartbeat",
                     fontsize=9, loc="left")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=72)
        plt.close(fig)
        buf.seek(0)
        imgs.append(Image.open(buf).convert("RGB"))
    pal = [im.quantize(colors=128, method=Image.Quantize.MEDIANCUT) for im in imgs]
    pal[0].save(path, save_all=True, append_images=pal[1:], duration=110, loop=0, optimize=True)
    imgs[-1].save(path.replace(".gif", "_final.png"))


def _busiest_window(A, dom, w=0.45):
    """Zoom box (x0, x1, s0, s1) of size w mm centred on the densest deposits."""
    if len(A) == 0:
        return (2.30, 2.30 + w, 0.85, 0.85 + w)
    x, s = _unroll(np.asarray(A))
    H, xe, se = np.histogram2d(x, s, bins=[np.arange(dom.x0, dom.x1, 0.05), np.arange(-0.6, 2.8, 0.05)])
    from scipy.ndimage import uniform_filter
    i, j = np.unravel_index(np.argmax(uniform_filter(H, size=int(w / 0.05))), H.shape)
    cx, cs = 0.5 * (xe[i] + xe[i + 1]), 0.5 * (se[j] + se[j + 1])
    return (cx - w / 2, cx + w / 2, cs - w / 2, cs + w / 2)


def growth_gif(adhered_csv, dom, stl, path, t_end, n_frames=60):
    """Animate deposit accumulation from adhered_platelets.csv (full view + zoom)."""
    rows = list(csv.DictReader(open(adhered_csv)))
    X = np.array([[float(r["x_mm"]), float(r["y_mm"]), float(r["z_mm"])] for r in rows])
    T = np.array([float(r["time_s"]) for r in rows])
    K = np.array([r["site"] for r in rows])
    wgt = float(rows[0]["platelets_represented"]) if rows else 1.0
    outl = _stent_outline(stl, dom)
    zoom = _busiest_window(X, dom)
    x, s = _unroll(X) if len(X) else (np.zeros(0), np.zeros(0))
    imgs = []
    for t in np.linspace(t_end / n_frames, t_end, n_frames):
        fig, axs = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [1.5, 1]})
        for ax, lim, ms in ((axs[0], None, 6), (axs[1], zoom, 40)):
            for o in outl:
                ax.fill(o[:, 0], o[:, 1], color="0.82", lw=0)
            m = T <= t
            for site, col, mk, lab in (("stent_luminal", "red", "o", "on strut top (faces the flow)"),
                                       ("stent_side", "darkorange", "o", "on strut side / corner"),
                                       ("aggregate", "darkred", "s", "aggregated onto platelets")):
                mm = m & (K == site)
                ax.scatter(x[mm], s[mm], s=ms, c=col, marker=mk, lw=0, label=lab)
            ax.set_aspect("equal")
            if lim:
                ax.set_xlim(lim[0], lim[1])
                ax.set_ylim(lim[2], lim[3])
                ax.set_title("zoom: densest deposits", fontsize=10)
            else:
                ax.set_xlim(dom.x0 - 0.2, dom.x1 + 0.3)
                ax.plot([zoom[0], zoom[1], zoom[1], zoom[0], zoom[0]],
                        [zoom[2], zoom[2], zoom[3], zoom[3], zoom[2]], "k--", lw=0.8)
                ax.set_title(f"platelet deposits on the bare stent   t = {t:.2f} s   "
                             f"({int(m.sum() * wgt):,} platelets)", fontsize=10, loc="left")
                ax.legend(loc="upper right", fontsize=8)
            ax.set_xlabel("axial (mm)  → flow")
        axs[0].set_ylabel("circumferential (mm)")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=72)
        plt.close(fig)
        buf.seek(0)
        imgs.append(Image.open(buf).convert("RGB"))
    pal = [im.quantize(colors=96, method=Image.Quantize.MEDIANCUT) for im in imgs]
    pal[0].save(path, save_all=True, append_images=pal[1:], duration=120, loop=0, optimize=True)


def _deposit_map(recs, weight, dom, stl, path, t_end):
    outl = _stent_outline(stl, dom)
    fig, ax = plt.subplots(figsize=(10, 6))
    for o in outl:
        ax.fill(o[:, 0], o[:, 1], color="0.85", lw=0)
    if recs:
        X = np.array([[r[2], r[3], r[4]] for r in recs])
        x, s = _unroll(X)
        k = np.array([r[5] for r in recs])
        tt = np.array([r[1] for r in recs])
        sc = ax.scatter(x, s, c=tt, s=np.where(k == "aggregate", 9, 5), cmap="plasma", vmin=0, vmax=t_end,
                        marker="o", lw=0)
        fig.colorbar(sc, ax=ax, shrink=0.7, label="time of adhesion (s)")
    ax.set_aspect("equal")
    ax.set_xlim(dom.x0 - 0.2, dom.x1 + 0.3)
    ax.set_xlabel("axial (mm)  → flow")
    ax.set_ylabel("circumferential (mm)")
    ax.set_title(f"where platelets deposited after {t_end:.1f} s (each dot = {weight:.0f} platelets)")
    plt.tight_layout()
    plt.savefig(path, dpi=120)
    plt.close(fig)


def rerender(out_dir, dom, stent_stl, period, title):
    """Redraw platelets.gif, deposition_map.png and deposition_growth.gif from a finished run
    (frames.npz + adhered_platelets.csv) without re-running the simulation."""
    z = np.load(os.path.join(out_dir, "frames.npz"), allow_pickle=True)
    series = json.loads(str(z["series"]))
    rows = list(csv.DictReader(open(os.path.join(out_dir, "adhered_platelets.csv"))))
    weight = float(rows[0]["platelets_represented"]) if rows else 1.0
    recs = [(int(r["parcel_id"]), float(r["time_s"]), float(r["x_mm"]), float(r["y_mm"]), float(r["z_mm"]),
             r["site"], int(r["was_activated"])) for r in rows]
    rt = np.array([r[1] for r in recs])
    rx = np.array([r[2:5] for r in recs]).reshape(-1, 3)
    rk = [r[5] for r in recs]
    frames = []
    for t, P, S in zip(z["t"], z["free_pos"], z["free_state"]):
        m = rt <= t + 1e-12
        frames.append((float(t), np.asarray(P), np.asarray(S), rx[m], [k for k, mm in zip(rk, m) if mm]))
    t_end = float(series[-1]["time_s"])
    _gif(frames, series, weight, dom, stent_stl, period, os.path.join(out_dir, "platelets.gif"), title)
    _deposit_map(recs, weight, dom, stent_stl, os.path.join(out_dir, "deposition_map.png"), t_end)
    growth_gif(os.path.join(out_dir, "adhered_platelets.csv"), dom, stent_stl,
               os.path.join(out_dir, "deposition_growth.gif"), t_end)
