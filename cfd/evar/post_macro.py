"""Post-processing of the vessel-scale EVAR runs.

validation: compares the simulated wall shear stress in the straight tube with
            the analytical Womersley solution over the last cycle.
graft:      TAWSS / OSI / RRT on the graft wall, statistics per segment, maps,
            cycle-to-cycle periodicity check, and selection of the wall-patch
            site(s) whose shear waveform drives the texture model.

Definitions (He & Ku 1996; Himburg et al. 2004), rho * kinematic stress -> Pa:
    TAWSS = (1/T) int |tau| dt,  OSI = 0.5 (1 - |(1/T) int tau dt| / TAWSS),
    RRT = 1 / ((1 - 2 OSI) TAWSS)
"""
from __future__ import annotations

import csv
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

LOW_WSS = 0.4


def _reader(cdir):
    foam = os.path.join(cdir, "case.foam")
    open(foam, "a").close()
    r = pv.OpenFOAMReader(foam)
    r.enable_all_patch_arrays()
    r.cell_to_point_creation = False
    return r


def _wall(r, t):
    r.set_active_time_value(t)
    d = r.read()
    b = d["boundary"]
    key = [k for k in b.keys() if k.split("/")[-1] == "wall"][0]
    return b[key]


def read_wall_series(cdir, rho):
    """Face centres (mm), areas (mm2), time-averaged fields and the snapshot series."""
    r = _reader(cdir)
    times = [t for t in r.time_values if t > 0]
    w_last = _wall(r, times[-1])
    sized = w_last.compute_cell_sizes(length=False, volume=False)
    area = np.asarray(sized.cell_data["Area"]) * 1e6
    cen = np.asarray(w_last.cell_centers().points) * 1e3
    cd = w_last.cell_data
    # OpenFOAM's wallShearStress points opposite to the near-wall flow; flip it so that
    # positive = stress in the flow direction (the convention of the Womersley solution)
    tawss = rho * np.asarray(cd["magWSSMean"])
    mean_vec = -rho * np.asarray(cd["wallShearStressMean"])
    snaps = np.stack([-rho * np.asarray(_wall(r, t).cell_data["wallShearStress"]) for t in times])
    return dict(times=np.array(times), area=area, cen=cen, tawss=tawss, mean_vec=mean_vec, snaps=snaps)


def metrics(tawss, mean_vec):
    mean_mag = np.linalg.norm(mean_vec, axis=1)
    osi = np.clip(0.5 * (1.0 - mean_mag / np.maximum(tawss, 1e-12)), 0.0, 0.5)
    rrt = 1.0 / np.maximum(mean_mag, 1e-9)
    return osi, rrt


def _stats(mask, area, tawss, osi, rrt):
    a = area[mask]
    A = a.sum()
    wq = lambda v, q: float(np.interp(q, np.cumsum(a[np.argsort(v)]) / A, np.sort(v)))  # noqa: E731
    return {"area_mm2": round(float(A), 1),
            "TAWSS_Pa": round(float((tawss[mask] * a).sum() / A), 4),
            "TAWSS_p05_Pa": round(wq(tawss[mask], 0.05), 4),
            "lowWSS_area_pct": round(100 * float(a[tawss[mask] < LOW_WSS].sum() / A), 2),
            "OSI": round(float((osi[mask] * a).sum() / A), 4),
            "OSI_gt_0.25_area_pct": round(100 * float(a[osi[mask] > 0.25].sum() / A), 2),
            "OSI_gt_0.45_area_pct": round(100 * float(a[osi[mask] > 0.45].sum() / A), 2),
            "RRT_median_perPa": round(wq(rrt[mask], 0.5), 3),
            "RRT_p95_perPa": round(wq(rrt[mask], 0.95), 3)}


def _probe_series(cdir):
    f = sorted(glob.glob(os.path.join(cdir, "postProcessing/wssProbes/*/wallShearStress")))
    if not f:
        return None, None
    import re
    rows = [ln for ln in open(f[0]) if not ln.startswith("#") and ln.strip()]
    t = np.array([float(r.split()[0]) for r in rows])
    num = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
    vals = np.array([[float(v) for v in num.findall(r.split(None, 1)[1])] for r in rows]).reshape(len(t), -1, 3)
    return t, vals


# ---------------------------------------------------------------------------
def analyse(cdir, rdir, kind, cfg, W, phys):
    os.makedirs(rdir, exist_ok=True)
    rho = phys.rho
    S = read_wall_series(cdir, rho)
    osi, rrt = metrics(S["tawss"], S["mean_vec"])
    T = cfg["inflow"]["period_s"]
    if kind == "validation":
        return _validation(cdir, rdir, cfg, W, S, osi, rho, T)
    return _graft(cdir, rdir, cfg, S, osi, rrt, rho, T)


def _validation(cdir, rdir, cfg, W, S, osi, rho, T):
    L = cfg["validation"]["length_mm"]
    sel = (S["cen"][:, 0] > 0.4 * L) & (S["cen"][:, 0] < 0.85 * L)
    t_an = np.linspace(0, T, 400, endpoint=False)
    tau_an = W.wall_shear(t_an)
    taw_an = float(np.mean(np.abs(tau_an)))
    osi_an = float(0.5 * (1 - abs(np.mean(tau_an)) / taw_an))
    # snapshot comparison (axial component, area-averaged over the developed region)
    a = S["area"][sel]
    tau_cfd = (S["snaps"][:, sel, 0] * a).sum(1) / a.sum()
    ph = S["times"] % T
    tau_an_s = W.wall_shear(ph)
    peak = np.max(np.abs(tau_an))
    out = {"womersley_number": W.womersley_number,
           "TAWSS_analytical_Pa": round(taw_an, 4),
           "TAWSS_cfd_Pa": round(float((S["tawss"][sel] * a).sum() / a.sum()), 4),
           "OSI_analytical": round(osi_an, 4),
           "OSI_cfd": round(float((osi[sel] * a).sum() / a.sum()), 4),
           "peak_tau_analytical_Pa": round(float(peak), 4),
           "peak_tau_cfd_Pa": round(float(np.max(np.abs(tau_cfd))), 4),
           "rms_error_pct_of_peak": round(100 * float(np.sqrt(np.mean((tau_cfd - tau_an_s) ** 2)) / peak), 2)}
    out["TAWSS_error_pct"] = round(100 * (out["TAWSS_cfd_Pa"] / out["TAWSS_analytical_Pa"] - 1), 2)
    tp, vp = _probe_series(cdir)
    fig, ax = plt.subplots(figsize=(7, 3.6))
    ax.plot(t_an, tau_an, "k-", lw=2, label="analytical (Womersley)")
    ax.plot(ph, tau_cfd, "o", color="C3", ms=5, label="OpenFOAM, area-average x = 24-51 mm")
    if tp is not None:
        last = tp > tp.max() - T
        vp = -vp
        o = np.argsort(tp[last] % T)
        ax.plot((tp[last] % T)[o], (rho * vp[last, 0, 0])[o], "-", color="C0", lw=1, alpha=0.7,
                label="OpenFOAM probe (x = 30 mm)")
        prev = (tp > tp.max() - 2 * T) & ~last
        if prev.any():
            a2 = np.interp(tp[last] % T, tp[prev] % T, rho * vp[prev, 0, 0], period=T)
            out["probe_cycle_to_cycle_rms_pct_of_peak"] = round(
                100 * float(np.sqrt(np.mean((rho * vp[last, 0, 0] - a2) ** 2)) / peak), 2)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel("time in cycle (s)")
    ax.set_ylabel("axial wall shear stress (Pa)")
    ax.set_title(f"Validation: 25 mm tube, Wo = {W.womersley_number:.1f}, RMS error "
                 f"{out['rms_error_pct_of_peak']}% of peak")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(rdir, "validation_womersley.png"), dpi=150)
    json.dump(out, open(os.path.join(rdir, "validation.json"), "w"), indent=2)
    print(json.dumps(out, indent=2))
    return out


def _segments(cen, g):
    import graft_geometry as gg
    x = cen[:, 0]
    seg = np.full(len(x), "", object)
    seg[x < -5] = "body"
    seg[(x >= -5) & (x < g.L_trans + 5)] = "bifurcation"
    seg[(x >= g.L_trans + 5) & (x < g.x_flare0)] = "leg"
    seg[(x >= g.x_flare0) & (x < g.x_flare1)] = "flare"
    seg[(x >= g.x_flare1) & (x < g.x_end - 15)] = "limb"
    seg[x >= g.x_end - 15] = "outlet_zone"
    c, r = gg.centreline(x, g)
    side = np.where(cen[:, 1] < c, "inner", "outer")      # inner = facing the other limb
    side[x < 0] = ""
    return seg, side


def _graft(cdir, rdir, cfg, S, osi, rrt, rho, T):
    import graft_geometry as gg
    g = gg.GraftParams(**cfg.get("graft", {}))
    cen, area, tawss = S["cen"], S["area"], S["tawss"]
    seg, side = _segments(cen, g)
    stats = {"whole_graft (excl. outlet zone)": _stats(seg != "outlet_zone", area, tawss, osi, rrt)}
    for s in ("body", "bifurcation", "leg", "flare", "limb"):
        stats[s] = _stats(seg == s, area, tawss, osi, rrt)
        for sd in ("inner", "outer"):
            m = (seg == s) & (side == sd)
            if m.sum() > 20:
                stats[f"{s}_{sd}"] = _stats(m, area, tawss, osi, rrt)
    # periodicity from probes
    tp, vp = _probe_series(cdir)
    period_check = {}
    if tp is not None and tp.max() > 2 * T:
        last = tp > tp.max() - T
        prev = (tp > tp.max() - 2 * T) & ~last
        for k in range(vp.shape[1]):
            a = rho * np.linalg.norm(vp[last, k], axis=1)
            b = np.interp(tp[last] % T, tp[prev] % T, rho * np.linalg.norm(vp[prev, k], axis=1), period=T)
            period_check[f"probe{k}"] = round(100 * float(abs(a.mean() - b.mean()) / max(a.mean(), 1e-9)), 2)
    # ---- patch sites: neighbourhood-averaged RRT in the legs/flare/limbs
    tree = cKDTree(cen)
    cand = np.where(np.isin(seg, ["leg", "flare", "limb"]))[0]
    nb = tree.query_ball_point(cen[cand], 2.0)
    rrt_s = np.array([np.average(rrt[n], weights=area[n]) for n in nb])
    taw_s = np.array([np.average(tawss[n], weights=area[n]) for n in nb])
    sites = {"worst_limb_site": cand[np.argmax(rrt_s)],
             "median_limb_site": cand[np.argsort(taw_s)[len(taw_s) // 2]]}
    site_info = {}
    for nm, i in sites.items():
        n = tree.query_ball_point(cen[i], 2.0)
        ser = (S["snaps"][:, n, :] * area[n][None, :, None]).sum(1) / area[n].sum()     # (nt, 3)
        # principal shear direction over the cycle (sign kept)
        u, sv, vt = np.linalg.svd(ser - 0.0, full_matrices=False)
        e1 = vt[0] if np.dot(vt[0], ser.mean(0)) >= 0 else -vt[0]
        tau_s = ser @ e1
        tau_perp = np.linalg.norm(ser - np.outer(tau_s, e1), axis=1)
        ph = S["times"] % T
        o = np.argsort(ph)
        site_info[nm] = {"centre_mm": cen[i].round(2).tolist(), "segment": str(seg[i]), "side": str(side[i]),
                         "TAWSS_Pa": round(float(np.average(tawss[n], weights=area[n])), 4),
                         "OSI": round(float(np.average(osi[n], weights=area[n])), 4),
                         "RRT_perPa": round(float(np.average(rrt[n], weights=area[n])), 3),
                         "principal_direction": e1.round(4).tolist(),
                         "perp_to_principal_rms_pct": round(100 * float(np.sqrt(np.mean(tau_perp ** 2))
                                                                        / np.sqrt(np.mean(tau_s ** 2))), 1)}
        with open(os.path.join(rdir, f"site_{nm}.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["phase_s", "tau_principal_Pa", "tau_perp_Pa"])
            for k in o:
                w.writerow([f"{ph[k]:.5f}", f"{tau_s[k]:.6f}", f"{tau_perp[k]:.6f}"])
    out = {"segments": stats, "patch_sites": site_info, "periodicity_TAWSS_change_last2cycles_pct": period_check}
    json.dump(out, open(os.path.join(rdir, "graft_metrics.json"), "w"), indent=2)
    _maps(rdir, cen, tawss, osi, rrt, site_info, g)
    _site_plot(rdir, site_info)
    print(json.dumps({k: stats[k] for k in ("whole_graft (excl. outlet zone)", "leg", "limb")}, indent=2))
    print(json.dumps(site_info, indent=2))
    return out


def _maps(rdir, cen, tawss, osi, rrt, sites, g):
    fig, axs = plt.subplots(3, 2, figsize=(12, 7.5), gridspec_kw={"width_ratios": [1, 1]})
    views = [("top view (looking down z)", cen[:, 2] > 0.5, 0, 1), ("side view (looking along y)", cen[:, 1] > 0, 0, 2)]
    fields = [("TAWSS (Pa)", tawss, "viridis", (0, 0.6)), ("OSI", osi, "magma", (0, 0.5)),
              ("RRT (1/Pa)", np.log10(rrt), "inferno", (0, 3))]
    for row, (lab, f, cm, lim) in enumerate(fields):
        for col, (vt, m, i, j) in enumerate(views):
            ax = axs[row, col]
            o = np.argsort(cen[m, 3 - i - j] if col == 0 else cen[m, 1])
            sc = ax.scatter(cen[m, i][o], cen[m, j][o], c=f[m][o], s=1.2, cmap=cm, vmin=lim[0], vmax=lim[1])
            ax.set_aspect("equal")
            ax.set_title(f"{lab}: {vt}" if lab[:3] != "RRT" else f"log10 RRT: {vt}", fontsize=9)
            if col == 0:
                for nm, s in sites.items():
                    ax.plot(s["centre_mm"][0], s["centre_mm"][1], "c*" if "worst" in nm else "w*", ms=10, mec="k")
            plt.colorbar(sc, ax=ax, shrink=0.8)
    fig.suptitle("Idealized EVAR endograft (quarter model, mirrored halves not shown). "
                 "Cyan star: worst limb site; white star: median limb site", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(rdir, "graft_maps.png"), dpi=140)
    plt.close(fig)


def _site_plot(rdir, sites):
    fig, ax = plt.subplots(figsize=(7, 3.4))
    for nm in sites:
        rows = list(csv.DictReader(open(os.path.join(rdir, f"site_{nm}.csv"))))
        t = [float(r["phase_s"]) for r in rows]
        ax.plot(t, [float(r["tau_principal_Pa"]) for r in rows], "-o", ms=3,
                label=f"{nm.replace('_', ' ')} (TAWSS {sites[nm]['TAWSS_Pa']:.2f} Pa, OSI {sites[nm]['OSI']:.2f})")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.axhline(LOW_WSS, color="C3", lw=0.6, ls="--")
    ax.axhline(-LOW_WSS, color="C3", lw=0.6, ls="--")
    ax.set_xlabel("time in cycle (s)")
    ax.set_ylabel("wall shear stress along principal direction (Pa)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(rdir, "patch_site_waveforms.png"), dpi=150)
    plt.close(fig)
