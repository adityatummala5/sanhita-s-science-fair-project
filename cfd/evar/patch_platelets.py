"""Lagrangian platelet adhesion on a textured graft patch (reactive graft surface).

Same kinetics as the coronary pilot (pipeline/platelets.py; Sorensen et al.
1999, Ann Biomed Eng 27:436), adapted to the periodic wall patch:

Transport   dx = u(x, t) dt + sqrt(2 D dt) N(0, 1), with u rebuilt quasi-steadily
            from the steady shear-level solutions and the site waveform, and
            D = D_B + 0.15 a^2 Hct (1 - Hct)^0.8 * gamma_local (shear-enhanced
            diffusivity, Zydney & Colton 1988).
Domain      periodic in x and z; platelets are tracked in the near-wall layer
            from the wall to `layer_um` above the nominal surface. The top
            `reservoir_um` of that layer is held at the bulk platelet count
            (Dirichlet condition: blood above the layer is well mixed).
            Each particle is one real platelet.
Adhesion    the whole graft surface is reactive (foreign polymer). A platelet
            within delta of the wall adheres with probability
            1 - exp(-k dt / delta) * (1 - local coverage / p_inf),
            k = k_rs (resting) or k_as (activated).
Activation  SIMPLIFIED agonist release: a free platelet passing within 10 um of a
            deposited platelet becomes activated.
Aggregation an activated platelet within one diameter of a deposit binds with k_aa.
Not modelled: red cells / margination, agonist fields, detachment, thrombus
feedback on the flow, shear-induced activation. Rate constants were measured
on collagen, so absolute counts are an upper bound; use the comparison
between designs.
"""
from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

import patch as P  # noqa: E402

RESTING, ACTIVATED = 0, 1


@dataclass
class PatchPlateletParams:
    concentration_per_uL: float = 250_000.0
    diameter_um: float = 3.0
    D_brownian_m2_s: float = 1.58e-13
    hematocrit: float = 0.40
    rbc_radius_um: float = 4.2
    k_rs_m_s: float = 3.7e-6
    k_as_m_s: float = 4.6e-5
    k_aa_m_s: float = 4.6e-5
    reaction_layer_um: float = 3.0
    activation_radius_um: float = 10.0
    p_inf_per_m2: float = 7.0e10
    layer_um: float = 70.0            # tracked layer above the nominal surface
    reservoir_um: float = 10.0
    layer_ref: str = "nominal"        # nominal: layer measured from y = 0; crest: from the texture's highest point
    duration_s: float = 30.0
    dt_s: float = 0.005
    coverage_bin_um: float = 10.0
    ds_max_um: float = 4.0            # max advective travel per sub-step
    seed: int = 1


def _alpha(p):
    a = p.rbc_radius_um * 1e-6
    return 0.15 * a * a * p.hematocrit * (1 - p.hematocrit) ** 0.8


class QuasiSteadyFlow:
    """Velocity (mm/s) and shear rate (1/s) at arbitrary points and times."""

    def __init__(self, levels, U_t, t, tex):
        self.Lx, self.Lz = tex.domain()
        self.t, self.T = t, t[-1] + (t[1] - t[0])
        self.U_t = U_t
        Uk = np.array([L["u_lid"] for L in levels])
        self.lu = np.log(Uk)
        c = levels[0]["cell_c"] * 1e3
        self.off = 1.0
        self.tree = cKDTree(np.c_[c[:, 0] % self.Lx, c[:, 1] + self.off, c[:, 2] % self.Lz],
                            boxsize=[self.Lx, 1e3, self.Lz])
        self.V = np.stack([L["U"] / u for L, u in zip(levels, Uk)]) * 1e3          # mm/s per (m/s of lid)
        self.G = np.stack([L["gamma"] / u for L, u in zip(levels, Uk)])            # 1/s per (m/s of lid)

    def _w(self, ti):
        a = abs(self.U_t[ti])
        la = math.log(min(max(a, math.exp(self.lu[0])), math.exp(self.lu[-1])))
        k = int(np.clip(np.searchsorted(self.lu, la) - 1, 0, len(self.lu) - 2))
        return k, (la - self.lu[k]) / (self.lu[k + 1] - self.lu[k])

    def sample(self, x, ti):
        q = np.c_[x[:, 0] % self.Lx, x[:, 1] + self.off, x[:, 2] % self.Lz]
        d, idx = self.tree.query(q, k=6)
        w = 1.0 / np.maximum(d, 1e-5) ** 2
        w /= w.sum(1, keepdims=True)
        k, f = self._w(ti)
        V = (1 - f) * self.V[k] + f * self.V[k + 1]
        G = (1 - f) * self.G[k] + f * self.G[k + 1]
        U = self.U_t[ti]
        u = np.einsum("nk,nkj->nj", w, V[idx]) * U
        g = np.einsum("nk,nk->n", w, G[idx]) * abs(U)
        return u, g


class Wall:
    """Heightfield wall: inside test y > h(x, z), distance to a dense surface sample."""

    def __init__(self, tex, spacing_um=1.5):
        self.tex = tex
        self.Lx, self.Lz = tex.domain()
        nx = max(8, int(self.Lx * 1e3 / spacing_um))
        nz = max(8, int(self.Lz * 1e3 / spacing_um))
        X, Z = np.meshgrid((np.arange(nx) + 0.5) * self.Lx / nx, (np.arange(nz) + 0.5) * self.Lz / nz, indexing="ij")
        Y = tex.height(X, Z)
        self.h_min, self.h_max = float(Y.min()), float(Y.max())
        self.pts = np.c_[X.ravel(), Y.ravel() + 1.0, Z.ravel()]
        self.tree = cKDTree(self.pts, boxsize=[self.Lx, 1e3, self.Lz])
        # cached heightfield on a periodic grid (0.75 um) for fast bilinear lookup
        self.gx = max(16, int(self.Lx * 1e3 / 0.75))
        self.gz = max(16, int(self.Lz * 1e3 / 0.75))
        Xg, Zg = np.meshgrid(np.arange(self.gx) * self.Lx / self.gx, np.arange(self.gz) * self.Lz / self.gz,
                             indexing="ij")
        self.H = tex.height(Xg, Zg)

    def height(self, x):
        fx = (x[:, 0] % self.Lx) / self.Lx * self.gx
        fz = (x[:, 2] % self.Lz) / self.Lz * self.gz
        i0, k0 = np.floor(fx).astype(int) % self.gx, np.floor(fz).astype(int) % self.gz
        i1, k1 = (i0 + 1) % self.gx, (k0 + 1) % self.gz
        wx, wz = fx - np.floor(fx), fz - np.floor(fz)
        H = self.H
        return ((1 - wx) * (1 - wz) * H[i0, k0] + wx * (1 - wz) * H[i1, k0]
                + (1 - wx) * wz * H[i0, k1] + wx * wz * H[i1, k1])

    def distance(self, x):
        d, i = self.tree.query(np.c_[x[:, 0] % self.Lx, x[:, 1] + 1.0, x[:, 2] % self.Lz])
        return d, i


def simulate(tex, levels, U_t, t, out_dir, p: PatchPlateletParams, label=""):
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(p.seed)
    flow = QuasiSteadyFlow(levels, U_t, t, tex)
    wall = Wall(tex)
    Lx, Lz = tex.domain(cfd=False)          # tiled domain; flow and wall lookups wrap the unit cell
    planform = Lx * Lz
    c_mm3 = p.concentration_per_uL                 # 1 uL = 1 mm3
    top = p.layer_um / 1e3 + (max(wall.h_max, 0.0) if p.layer_ref == "crest" else 0.0)
    res_lo = top - p.reservoir_um / 1e3
    delta = p.reaction_layer_um / 1e3
    alpha = _alpha(p)
    dt = p.dt_s
    nsteps = int(round(p.duration_s / dt))
    nt = len(t)

    def seed_uniform(n, y_lo, y_hi):
        out = np.empty((0, 3))
        while len(out) < n:
            m = 2 * (n - len(out)) + 16
            c = np.c_[rng.uniform(0, Lx, m), rng.uniform(y_lo, y_hi, m), rng.uniform(0, Lz, m)]
            c = c[c[:, 1] > wall.height(c) + 0.5 * p.diameter_um / 1e3]
            out = np.vstack([out, c])
        return out[:n]

    # initial: uniform bulk concentration in the fluid of the layer
    fluid_vol = planform * (top - np.mean(wall.height(np.c_[rng.uniform(0, Lx, 20000), np.zeros(20000),
                                                              rng.uniform(0, Lz, 20000)])))
    n0 = rng.poisson(c_mm3 * fluid_vol)
    x = seed_uniform(n0, wall.h_min, top)
    state = np.zeros(len(x), np.int8)
    n_res_target = c_mm3 * planform * (top - res_lo)
    # coverage bins
    nbx, nbz = max(1, int(Lx * 1e3 / p.coverage_bin_um)), max(1, int(Lz * 1e3 / p.coverage_bin_um))
    cov = np.zeros((nbx, nbz))
    bin_area_m2 = None
    sx, sz = np.meshgrid(np.arange(nbx), np.arange(nbz), indexing="ij")
    # wetted area per bin (from the heightfield)
    xs = (np.arange(nbx * 8) + 0.5) * Lx / (nbx * 8)
    zs = (np.arange(nbz * 8) + 0.5) * Lz / (nbz * 8)
    XX, ZZ = np.meshgrid(xs, zs, indexing="ij")
    hh = tex.height(XX % wall.Lx, ZZ % wall.Lz)
    gx, gz = np.gradient(hh, xs[1] - xs[0], zs[1] - zs[0])
    wa = np.sqrt(1 + gx ** 2 + gz ** 2) * (xs[1] - xs[0]) * (zs[1] - zs[0])
    bin_area_m2 = wa.reshape(nbx, 8, nbz, 8).sum((1, 3)) * 1e-6
    dep = []                                    # (t, x, y, z, kind, activated)
    dep_tree, dep_pts = None, np.empty((0, 3))
    series = []
    p_rs = 1 - math.exp(-p.k_rs_m_s * dt / (delta * 1e-3))
    p_as = 1 - math.exp(-p.k_as_m_s * dt / (delta * 1e-3))
    p_aa = 1 - math.exp(-p.k_aa_m_s * dt / (delta * 1e-3))
    n_sub_total = 0
    for step in range(nsteps):
        ti = step % nt
        tt = step * dt
        # advection-diffusion with adaptive midpoint (RK2) sub-steps: each sub-step moves a
        # platelet at most ds_max, so particles follow curved streamlines over the texture
        # (a single explicit 5 ms step would ratchet them away from the wall)
        u, g = flow.sample(x, ti)
        nsub = max(1, int(math.ceil(np.abs(u).max() * dt / (p.ds_max_um / 1e3)))) if len(x) else 1
        h = dt / nsub
        for _ in range(nsub):
            u, g = flow.sample(x, ti)
            um, _g = flow.sample(x + 0.5 * h * u, ti)
            D = (p.D_brownian_m2_s + alpha * g) * 1e6                  # mm2/s
            xn = x + h * um + np.sqrt(2 * D * h)[:, None] * rng.standard_normal(x.shape)
            bad = xn[:, 1] < wall.height(xn) + 0.5 * p.diameter_um / 1e3       # no wall penetration
            xn[bad] = x[bad]
            x = xn
            x[:, 0] %= Lx
            x[:, 2] %= Lz
        n_sub_total += nsub
        # activation near deposits, aggregation onto deposits
        if len(dep_pts):
            dd, di = dep_tree.query(np.c_[x[:, 0], x[:, 1] + 1.0, x[:, 2]],
                                    distance_upper_bound=p.activation_radius_um / 1e3)
            state[np.isfinite(dd)] = ACTIVATED
            close = (dd < p.diameter_um / 1e3) & (state == ACTIVATED)
            agg = close & (rng.random(len(x)) < p_aa)
        else:
            agg = np.zeros(len(x), bool)
        # adhesion to the graft surface
        dw, _ = wall.distance(x)
        near = (dw < delta + 0.5 * p.diameter_um / 1e3) & ~agg
        bi = np.clip((x[:, 0] / Lx * nbx).astype(int), 0, nbx - 1)
        bk = np.clip((x[:, 2] / Lz * nbz).astype(int), 0, nbz - 1)
        free_frac = np.clip(1 - cov[bi, bk] / (p.p_inf_per_m2 * bin_area_m2[bi, bk]), 0, 1)
        prob = np.where(state == ACTIVATED, p_as, p_rs) * free_frac
        adh = near & (rng.random(len(x)) < prob)
        new = adh | agg
        if new.any():
            for k in np.where(new)[0]:
                dep.append((tt, x[k, 0], x[k, 1], x[k, 2], "aggregate" if agg[k] else "surface", int(state[k])))
            np.add.at(cov, (bi[adh], bk[adh]), 1)
            dep_pts = np.vstack([dep_pts, x[new]])
            dep_tree = cKDTree(np.c_[dep_pts[:, 0] % Lx, dep_pts[:, 1] + 1.0, dep_pts[:, 2] % Lz],
                               boxsize=[Lx, 1e3, Lz])
            keep = ~new
            x, state = x[keep], state[keep]
        # top boundary + reservoir (Dirichlet: bulk concentration)
        keep = x[:, 1] < top
        x, state = x[keep], state[keep]
        inres = x[:, 1] > res_lo
        n_res = int(inres.sum())
        want = rng.poisson(n_res_target)
        if n_res < want:
            x = np.vstack([x, seed_uniform(want - n_res, res_lo, top)])
            state = np.concatenate([state, np.zeros(want - n_res, np.int8)])
        elif n_res > want:
            drop = rng.choice(np.where(inres)[0], n_res - want, replace=False)
            m = np.ones(len(x), bool)
            m[drop] = False
            x, state = x[m], state[m]
        if step % 20 == 0:
            n_surf = sum(1 for d in dep if d[4] == "surface")
            series.append((tt, len(x), n_surf, len(dep) - n_surf, int((state == ACTIVATED).sum())))
    # ---- outputs
    with open(os.path.join(out_dir, "deposits.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t_s", "x_mm", "y_mm", "z_mm", "kind", "was_activated"])
        w.writerows(dep)
    with open(os.path.join(out_dir, "timeseries.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t_s", "free", "adhered", "aggregated", "activated_free"])
        w.writerows(series)
    D = np.array([(d[1], d[2], d[3]) for d in dep]) if dep else np.empty((0, 3))
    n_surf = sum(1 for d in dep if d[4] == "surface")
    hrel = None
    if len(D):
        h = wall.height(D)
        rng_h = max(wall.h_max - wall.h_min, 1e-9)
        hrel = (h - wall.h_min) / rng_h
    wet = float(bin_area_m2.sum() * 1e6)
    t_half = 0.5 * p.duration_s
    late = [d for d in dep if d[0] >= t_half]
    summ = {"design": tex.name, "label": label, "duration_s": p.duration_s, "planform_mm2": planform,
            "wetted_mm2": wet, "deposited_total": len(dep), "adhered": n_surf, "aggregated": len(dep) - n_surf,
            "deposited_per_mm2_planform": round(len(dep) / planform, 1),
            "adhered_per_mm2_wetted": round(n_surf / wet, 1),
            "counting_error_per_mm2_planform": round(math.sqrt(max(len(dep), 1)) / planform, 1),
            # primary metric: deposition rate once the start-up (initial near-wall platelets) is over
            "rate_late_per_mm2_s": round(len(late) / planform / (p.duration_s - t_half), 2),
            "rate_late_counting_error": round(math.sqrt(max(len(late), 1)) / planform / (p.duration_s - t_half), 2),
            "pct_in_lower_third": round(100 * float((hrel < 1 / 3).mean()), 1) if hrel is not None and wall.h_max > wall.h_min else None,
            "pct_in_upper_third": round(100 * float((hrel > 2 / 3).mean()), 1) if hrel is not None and wall.h_max > wall.h_min else None,
            "mean_substeps_per_step": round(n_sub_total / max(nsteps, 1), 2),
            "params": asdict(p)}
    json.dump(summ, open(os.path.join(out_dir, "summary.json"), "w"), indent=2)
    # map
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    n = 300
    X, Z = np.meshgrid((np.arange(n) + .5) * Lx / n, (np.arange(n) + .5) * Lz / n, indexing="ij")
    Lxc, Lzc = tex.domain()
    ax.pcolormesh(X, Z, tex.height(X % Lxc, Z % Lzc) * 1e3, cmap="Greys_r", shading="auto")
    if len(D):
        kinds = np.array([d[4] for d in dep])
        ax.scatter(D[kinds == "surface", 0], D[kinds == "surface", 2], s=1.5, c="C3", label="adhered")
        ax.scatter(D[kinds == "aggregate", 0], D[kinds == "aggregate", 2], s=1.5, c="C0", label="aggregated")
        ax.legend(fontsize=7, loc="upper right", markerscale=4)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm), flow ->")
    ax.set_title(f"{tex.name}: {summ['deposited_per_mm2_planform']:.0f} platelets/mm2 in {p.duration_s:.0f} s",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "deposition_map.png"), dpi=140)
    plt.close(fig)
    return summ
