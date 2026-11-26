"""Rebuild the cardiac cycle on each textured patch and compute surface metrics.

Quasi-steady reconstruction (see patch.py for why this is valid here):
  1. Site waveform tau_site(t) (signed, along the principal direction) from the
     vessel-scale run. For the smooth wall the lid speed that produces it is
     U(t) = sign(tau) * Hlid * gamma(|tau|), gamma from the Carreau law.
  2. Each design was solved steadily at lid speeds U_k. Wall-stress and
     velocity fields are normalised by the smooth-wall stress / lid speed at the
     same U_k (an 'amplification' field that varies only weakly with U through
     shear thinning), interpolated in log|U|, and multiplied back by
     tau_site(t) / U(t). Stokes flow is reversible, so negative U flips the
     field's sign.

Metrics (per design, per site; area statistics use wetted area):
  TAWSS, OSI, RRT (He & Ku 1996; Himburg 2004)
  low-WSS area: TAWSS < 0.4 Pa (Malek 1999) and TAWSS < 50% of the smooth wall
  OSI > 0.25 / > 0.45 area, RRT median / 95th percentile
  WSSG: 95th percentile of the surface gradient of TAWSS (Pa/mm)
  slow-flow volume: fluid per mm2 of graft moving slower (cycle-mean speed) than
     fluid 5 um above a smooth wall - i.e. near-wall fluid with platelet-scale
     residence. Smooth wall = 5.0e-3 mm3/mm2 by construction.
  stress per planform area: cycle-mean |sum tau_x dA| / planform area, relative
     to smooth (the drag the texture exerts on the blood).
"""
from __future__ import annotations

import csv
import json
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402

import patch as P  # noqa: E402

LOW_WSS = 0.4
NT = 200


def site_waveform(site, cfg, ph, res_root):
    """Return (t, tau_signed_Pa) over one period."""
    T = cfg["inflow"]["period_s"]
    if site == "analytical_limb":
        import flow
        import run_macro
        ts, q = run_macro.waveform(cfg)
        W = flow.Womersley(ts, 0.5 * q * 1e-6, T, cfg["graft"]["D_limb"] / 2 * 1e-3)
        t = np.linspace(0, T, NT, endpoint=False)
        return t, W.wall_shear(t), "Womersley, 16 mm limb, Newtonian (analytical)"
    path = site if site.endswith(".csv") else os.path.join(res_root, "graft-coarse", f"site_{site}_limb_site.csv")
    rows = list(csv.DictReader(open(path)))
    tp = np.array([float(r["phase_s"]) for r in rows])
    tau = np.array([float(r["tau_principal_Pa"]) for r in rows])
    t = np.linspace(0, T, NT, endpoint=False)
    # periodic cubic interpolation via FFT-free approach: linear is adequate at 40 samples/cycle
    return t, np.interp(t, tp, tau, period=T), os.path.relpath(path, res_root)


def load_levels(run_dir):
    import glob
    out = []
    for d in sorted(glob.glob(os.path.join(run_dir, "L*")), key=lambda p: int(os.path.basename(p)[1:])):
        f = os.path.join(d, "fields.npz")
        if os.path.exists(f):
            out.append(dict(np.load(f)))
    return out


def reconstruct(levels, smooth_levels, t, tau_site, ph, pm):
    """tau(face, t, 3), speed(cell, t)-mean, and the lid speed series."""
    H = pm.h_lid_mm * 1e-3
    U_t = np.sign(tau_site) * H * P.gamma_from_tau(tau_site, ph)
    Uk = np.array([L["u_lid"] for L in levels])
    # smooth-wall stress at each level (area-weighted mean x-stress)
    tau_s_k = np.array([np.sum(S["tau"][:, 0] * S["face_A"]) / S["face_A"].sum() for S in smooth_levels])
    amp = np.stack([L["tau"] / ts for L, ts in zip(levels, tau_s_k)])          # (K, F, 3)
    vel = np.stack([L["U"] / u for L, u in zip(levels, Uk)])                    # (K, C, 3)
    lu = np.log(Uk)
    a = np.abs(U_t)
    la = np.log(np.clip(a, Uk[0], Uk[-1]))
    idx = np.clip(np.searchsorted(lu, la) - 1, 0, len(Uk) - 2)
    w = (la - lu[idx]) / (lu[idx + 1] - lu[idx])
    F, C = amp.shape[1], vel.shape[1]
    tau = np.empty((len(t), F, 3))
    speed_mean = np.zeros(C)
    for i in range(len(t)):
        A = (1 - w[i]) * amp[idx[i]] + w[i] * amp[idx[i] + 1]
        tau[i] = A * tau_site[i]
        V = (1 - w[i]) * vel[idx[i]] + w[i] * vel[idx[i] + 1]
        speed_mean += np.linalg.norm(V * U_t[i], axis=1)
    speed_mean /= len(t)
    return tau, speed_mean, U_t


def _wq(v, a, q):
    o = np.argsort(v)
    return float(np.interp(q, np.cumsum(a[o]) / a.sum(), v[o]))


def surface_gradient(face_c, f, tex, n=160):
    """95th percentile |grad_s f| (Pa/mm) on a regular grid over the heightfield."""
    Lx, Lz = tex.domain()
    xs = (np.arange(n) + 0.5) * Lx / n
    zs = (np.arange(n) + 0.5) * Lz / n
    X, Z = np.meshgrid(xs, zs, indexing="ij")
    pts = face_c[:, [0, 2]] * 1e3
    G = griddata(pts, f, (X, Z), method="linear")
    G = np.where(np.isnan(G), griddata(pts, f, (X, Z), method="nearest"), G)
    h = tex.height(X, Z)
    hx, hz = np.gradient(h, Lx / n, Lz / n)
    fx, fz = np.gradient(G, Lx / n, Lz / n)
    # |grad_s f|^2 = grad f . (I + grad h grad h^T)^-1 . grad f
    det = 1 + hx ** 2 + hz ** 2
    g2 = (fx ** 2 * (1 + hz ** 2) + fz ** 2 * (1 + hx ** 2) - 2 * fx * fz * hx * hz) / det
    return float(np.percentile(np.sqrt(np.maximum(g2, 0)), 95))


def design_metrics(name, tex, levels, smooth_levels, t, tau_site, ph, pm, smooth_ref=None):
    tau, speed_mean, U_t = reconstruct(levels, smooth_levels, t, tau_site, ph, pm)
    L0 = levels[0]
    A = L0["face_A"]
    mag = np.linalg.norm(tau, axis=2)
    tawss = mag.mean(0)
    mvec = tau.mean(0)
    mmag = np.linalg.norm(mvec, axis=1)
    osi = np.clip(0.5 * (1 - mmag / np.maximum(tawss, 1e-12)), 0, 0.5)
    rrt = 1 / np.maximum(mmag, 1e-9)
    Lx, Lz = tex.domain()
    planform = Lx * Lz                                     # mm2
    Aw = A.sum() * 1e6                                     # mm2 wetted
    # slow-flow volume
    gam_abs_mean = np.mean(P.gamma_from_tau(tau_site, ph))
    u_ref = gam_abs_mean * 5e-6                            # speed 5 um above a smooth wall (m/s)
    V = L0["cell_V"]
    slow = V[speed_mean < u_ref].sum() * 1e9 / planform    # mm3 per mm2
    drag = abs(np.mean(np.sum(tau[:, :, 0] * A[None, :], axis=1))) / (planform * 1e-6)
    m = {"design": name,
         "wetted_area_ratio": round(Aw / planform, 4),
         "TAWSS_Pa": round(float(np.sum(tawss * A) / A.sum()), 4),
         "TAWSS_p05_Pa": round(_wq(tawss, A, 0.05), 4),
         "lowWSS_area_pct": round(100 * float(A[tawss < LOW_WSS].sum() / A.sum()), 2),
         "OSI": round(float(np.sum(osi * A) / A.sum()), 4),
         "OSI_gt_0.25_area_pct": round(100 * float(A[osi > 0.25].sum() / A.sum()), 2),
         "OSI_gt_0.45_area_pct": round(100 * float(A[osi > 0.45].sum() / A.sum()), 2),
         "RRT_median_perPa": round(_wq(rrt, A, 0.5), 3),
         "RRT_p95_perPa": round(_wq(rrt, A, 0.95), 3),
         "WSSG_p95_Pa_mm": round(surface_gradient(L0["face_c"], tawss, tex), 2),
         "slow_flow_volume_mm3_per_mm2": round(float(slow), 5),
         "planform_stress_Pa": round(float(drag), 4)}
    if smooth_ref is not None:
        m["area_below_half_smooth_TAWSS_pct"] = round(
            100 * float(A[tawss < 0.5 * smooth_ref["TAWSS_Pa"]].sum() / A.sum()), 2)
        for k in ("TAWSS_Pa", "RRT_median_perPa", "slow_flow_volume_mm3_per_mm2", "planform_stress_Pa", "OSI"):
            m[f"{k}_vs_smooth_pct"] = round(100 * (m[k] / smooth_ref[k] - 1), 1) if smooth_ref[k] else None
    fields = {"face_c": L0["face_c"], "face_A": A, "tawss": tawss, "osi": osi, "rrt": rrt}
    return m, fields, U_t


def run(designs, cfg, pm, ph, site, runs_root, res_root):
    evar_res = os.path.dirname(res_root)
    t, tau_site, src = site_waveform(site, cfg, ph, evar_res)
    out_dir = os.path.join(res_root, f"site-{os.path.basename(site).replace('.csv', '')}")
    os.makedirs(out_dir, exist_ok=True)
    smooth_levels = load_levels(os.path.join(runs_root, "smooth"))
    if not smooth_levels:
        raise SystemExit("run the 'smooth' design first: it calibrates every other design")
    # check smooth-wall calibration: Carreau Couette should give tau = mu(gamma) gamma exactly
    calib = []
    for S in smooth_levels:
        g = S["u_lid"] / (pm.h_lid_mm * 1e-3)
        tau_cfd = -float(np.sum(S["tau"][:, 0] * S["face_A"]) / S["face_A"].sum())   # OpenFOAM sign -> flow direction
        calib.append({"gamma_1_s": g, "tau_exact_Pa": float(P.carreau_mu(g, ph) * g), "tau_cfd_Pa": tau_cfd,
                      "error_pct": round(100 * (tau_cfd / (P.carreau_mu(g, ph) * g) - 1), 3)})
    rows, allf = [], {}
    sm, _, U_t = design_metrics("smooth", designs[0], smooth_levels, smooth_levels, t, tau_site, ph, pm)
    for tex in designs:
        lv = load_levels(os.path.join(runs_root, tex.name))
        if len(lv) != len(smooth_levels):
            print(f"  {tex.name}: missing levels, skipped")
            continue
        m, f, _ = design_metrics(tex.name, tex, lv, smooth_levels, t, tau_site, ph, pm, smooth_ref=sm)
        m["description"] = tex.description
        rows.append(m)
        allf[tex.name] = (tex, f)
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(out_dir, "patch_metrics.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    json.dump({"site": src, "TAWSS_site_Pa": float(np.mean(np.abs(tau_site))),
               "OSI_site": float(0.5 * (1 - abs(np.mean(tau_site)) / np.mean(np.abs(tau_site)))),
               "smooth_calibration": calib, "lid_speed_max_mm_s": float(np.max(np.abs(U_t)) * 1e3)},
              open(os.path.join(out_dir, "site_and_calibration.json"), "w"), indent=2)
    _plots(out_dir, rows, allf, t, tau_site)
    print(f"site: {src}")
    cols = ["design", "TAWSS_Pa", "TAWSS_Pa_vs_smooth_pct", "lowWSS_area_pct", "area_below_half_smooth_TAWSS_pct",
            "OSI", "RRT_median_perPa", "WSSG_p95_Pa_mm", "slow_flow_volume_mm3_per_mm2", "wetted_area_ratio"]
    print("  ".join(f"{c[:14]:>14s}" for c in cols))
    for r in rows:
        print("  ".join(f"{str(r.get(c, '')):>14s}"[:14 if c != 'design' else 22] for c in cols))
    return rows


def _plots(out_dir, rows, allf, t, tau_site):
    n = len(allf)
    nc = 4
    nr = math.ceil(n / nc)
    vmax = max(np.percentile(f["tawss"], 99) for _, f in allf.values())
    fig, axs = plt.subplots(nr, nc, figsize=(4 * nc, 3.4 * nr))
    for ax, (name, (tex, f)) in zip(axs.ravel(), allf.items()):
        c = f["face_c"] * 1e3
        sc = ax.scatter(c[:, 0], c[:, 2], c=f["tawss"], s=2.0, cmap="viridis", vmin=0, vmax=vmax, marker="s")
        ax.set_aspect("equal")
        r = next(r for r in rows if r["design"] == name)
        ax.set_title(f"{name}\nTAWSS {r['TAWSS_Pa']:.3f} Pa, slow vol {r['slow_flow_volume_mm3_per_mm2'] * 1e3:.1f} um",
                     fontsize=8)
        ax.set_xlabel("x (mm), flow ->", fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.colorbar(sc, ax=axs, shrink=0.6, label="TAWSS (Pa)")
    fig.savefig(os.path.join(out_dir, "patch_tawss_maps.png"), dpi=130)
    plt.close(fig)
    # bar chart vs smooth
    keys = [("TAWSS_Pa_vs_smooth_pct", "TAWSS vs smooth (%)"),
            ("area_below_half_smooth_TAWSS_pct", "area with TAWSS < 50% of smooth (%)"),
            ("slow_flow_volume_mm3_per_mm2_vs_smooth_pct", "slow-flow volume vs smooth (%)"),
            ("RRT_median_perPa_vs_smooth_pct", "median RRT vs smooth (%)")]
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.2))
    names = [r["design"] for r in rows]
    for ax, (k, lab) in zip(axs, keys):
        v = [r.get(k) or 0 for r in rows]
        ax.barh(names, v, color=["0.5" if nm in ("smooth", "pet") else "C0" for nm in names])
        ax.axvline(0, color="k", lw=0.6)
        ax.set_title(lab, fontsize=9)
        ax.invert_yaxis()
        ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "patch_comparison.png"), dpi=140)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(t, tau_site, "k-")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel("time in cycle (s)")
    ax.set_ylabel("site wall shear stress (Pa)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "site_waveform.png"), dpi=130)
    plt.close(fig)
