#!/usr/bin/env python3
"""Check of the quasi-steady patch assumption with an exact 1D unsteady solve.

The patch model treats the 0.3 mm gap as quasi-steady: at every instant the
wall stress is the steady Carreau-Couette value for the current lid speed.
The neglected effect is fluid inertia across the gap (patch Womersley number
~0.4 for the fundamental, higher for harmonics). Inertia across the much
smaller texture features (10-100 um) is smaller still, so the gap-scale 1D
problem bounds the error.

Solves  rho du/dt = d/dy( mu(|du/dy|) du/dy ),  u(0) = 0,  u(H) = U(t)
(implicit in time, Picard iteration on the Carreau viscosity), three cycles,
and compares the wall stress of the last cycle with the quasi-steady target.

    python3 evar/verify_quasisteady.py --site worst
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(CFD, "pipeline"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.linalg import solve_banded  # noqa: E402

import case as fc  # noqa: E402
import patch as P  # noqa: E402
import post_patch  # noqa: E402


def run(site, cycles=3, ny=200, steps=4000):
    cfg = json.load(open(os.path.join(CFD, "config", "evar.json")))
    ph = fc.Physics(**cfg["physics"])
    pm = P.PatchMesh(**cfg["patch"]["mesh"])
    T = cfg["inflow"]["period_s"]
    t, tau_site, src = post_patch.site_waveform(site, cfg, ph, os.path.join(CFD, "results", "evar"))
    H = pm.h_lid_mm * 1e-3
    U_t = np.sign(tau_site) * H * P.gamma_from_tau(tau_site, ph)
    y = np.linspace(0, H, ny + 1)
    dy = y[1] - y[0]
    dt = T / steps
    u = U_t[0] * y / H
    tw, tt = [], []
    for n in range(cycles * steps):
        tn = (n + 1) * dt
        Ub = np.interp(tn % T, t, U_t, period=T)
        un = u.copy()
        for _ in range(6):                               # Picard on viscosity
            g = np.abs(np.diff(un)) / dy
            mu = P.carreau_mu(g, ph)                     # face viscosities (ny)
            a = mu / dy ** 2 / ph.rho * dt
            # interior nodes 1..ny-1:  u_i - dt/rho * d/dy(mu du/dy) = u_old
            m = ny - 1
            ab = np.zeros((3, m))
            ab[0, 1:] = -a[1:m]
            ab[1, :] = 1 + a[:m] + a[1:m + 1]
            ab[2, :-1] = -a[1:m]
            rhs = u[1:ny].copy()
            rhs[-1] += a[m] * Ub
            un = np.concatenate([[0.0], solve_banded((1, 1), ab, rhs), [Ub]])
        u = un
        if n >= (cycles - 1) * steps:
            g0 = (u[1] - u[0]) / dy
            tw.append(P.carreau_mu(abs(g0), ph) * g0)
            tt.append(tn % T)
    tw, tt = np.array(tw), np.array(tt)
    o = np.argsort(tt)
    tw, tt = tw[o], tt[o]
    target = np.interp(tt, t, tau_site, period=T)
    peak = np.max(np.abs(target))
    taw_u, taw_q = np.mean(np.abs(tw)), np.mean(np.abs(target))
    osi = lambda v: 0.5 * (1 - abs(np.mean(v)) / np.mean(np.abs(v)))  # noqa: E731
    out = {"site": src, "gap_mm": pm.h_lid_mm,
           "rms_error_pct_of_peak": round(100 * float(np.sqrt(np.mean((tw - target) ** 2)) / peak), 2),
           "TAWSS_unsteady_Pa": round(float(taw_u), 4), "TAWSS_quasisteady_Pa": round(float(taw_q), 4),
           "TAWSS_diff_pct": round(100 * float(taw_u / taw_q - 1), 2),
           "OSI_unsteady": round(float(osi(tw)), 4), "OSI_quasisteady": round(float(osi(target)), 4)}
    res = os.path.join(CFD, "results", "evar", "patch")
    os.makedirs(res, exist_ok=True)
    json.dump(out, open(os.path.join(res, f"verify_quasisteady_1d_{site}.json"), "w"), indent=2)
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    ax.plot(tt, target, "k-", lw=2, label="quasi-steady (used by the patch model)")
    ax.plot(tt, tw, "C3--", lw=1.4, label="unsteady 1D Carreau Couette (exact inertia)")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel("time in cycle (s)")
    ax.set_ylabel("smooth-wall shear stress (Pa)")
    ax.set_title(f"Quasi-steady check, {site} site: RMS error {out['rms_error_pct_of_peak']}% of peak", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(res, f"verify_quasisteady_1d_{site}.png"), dpi=150)
    print(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", default="worst")
    run(ap.parse_args().site)
