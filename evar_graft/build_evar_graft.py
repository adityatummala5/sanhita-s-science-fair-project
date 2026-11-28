"""Parametric CAD generator: bifurcated abdominal-aortic EVAR stent graft whose
luminal (blood-contacting) surface is covered in overlapping shark denticles.

Units: millimetres. Run from anywhere:
    python3 evar_graft/build_evar_graft.py            # writes all STLs + preview.png
    python3 evar_graft/build_evar_graft.py --no-render

Outputs (next to this script):
    denticle_unit_true_scale.stl        one denticle, real size (200 x 150 x 40 um)
    denticle_patch_true_scale.stl       2 x 2 mm flat patch of graft wall, real-size overlapping denticles
    evar_graft_smooth.stl               full bifurcated graft, true dimensions, no texture (CFD control)
    evar_graft_cutaway_illustrative.stl full graft, front half removed, denticles ENLARGED (x10) so they
                                        are visible at whole-device scale. NOT to scale - see README.
    preview.png                         renders of all of the above

Coordinate system for the graft: z = device axis, proximal (renal) end at z = 0,
blood flows toward -z (into the iliac limbs). Patch: flow along +x, denticles on +z face.
"""
from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, replace

import manifold3d as mf
import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------
@dataclass
class Denticle:
    """Overlapping (imbricated) denticle. Lengths in mm."""
    length: float = 0.200      # streamwise crown length
    width: float = 0.150       # spanwise crown width
    height: float = 0.040      # crown height at the trailing edge (above the graft surface)
    lead_height: float = 0.12  # leading-edge height as a fraction of `height` (crown is inclined)
    n_ridges: int = 3          # longitudinal riblets on the crown
    ridge_depth: float = 0.35  # groove depth, fraction of local height
    tip_notch: float = 0.12    # trailing-edge scallop depth, fraction of length (gives 3 points)
    pitch_stream: float = 0.120  # row-to-row spacing along flow  (< length  -> rows overlap)
    pitch_span: float = 0.130    # spacing within a row            (< width   -> neighbours overlap)
    embed: float = 0.005       # base sunk into the wall for a clean boolean union
    grid: tuple = (16, 19)     # crown surface resolution (streamwise, spanwise)

    def scaled(self, k: float, grid=None) -> "Denticle":
        return replace(self, length=self.length * k, width=self.width * k, height=self.height * k,
                       pitch_stream=self.pitch_stream * k, pitch_span=self.pitch_span * k,
                       embed=self.embed * k, grid=grid or self.grid)


@dataclass
class Graft:
    """Generic bifurcated EVAR main body (not modelled on any one commercial device)."""
    d_main: float = 26.0       # proximal main-body OD (fits a ~22-23 mm neck at ~15% oversize)
    l_main: float = 50.0       # straight main-body length
    l_trans: float = 10.0      # taper from round body to the two-limb "stadium"
    d_limb: float = 14.0       # iliac limb OD
    l_limb: float = 70.0       # limb length (along its own axis)
    limb_angle: float = 6.0    # each limb diverges this many degrees from the device axis
    wall: float = 0.30         # graft fabric wall thickness
    wire_d: float = 0.50       # external Z-stent wire diameter
    segs: int = 256            # circular resolution


# --------------------------------------------------------------------------
# Mesh helpers
# --------------------------------------------------------------------------
def to_manifold(v, f) -> mf.Manifold:
    m = mf.Manifold(mf.Mesh(vert_properties=np.asarray(v, np.float32),
                            tri_verts=np.asarray(f, np.uint32)))
    if m.status() != mf.Error.NoError:
        raise ValueError(f"invalid manifold: {m.status()}")
    return m


def to_trimesh(m: mf.Manifold) -> trimesh.Trimesh:
    mesh = m.to_mesh()
    return trimesh.Trimesh(np.asarray(mesh.vert_properties)[:, :3], np.asarray(mesh.tri_verts),
                           process=False)


def affine(R: np.ndarray, t) -> list:
    """3x4 matrix for manifold.transform."""
    return np.c_[R, np.asarray(t, float)].tolist()


def rot_y(deg: float) -> np.ndarray:
    a = math.radians(deg)
    return np.array([[math.cos(a), 0, math.sin(a)], [0, 1, 0], [-math.sin(a), 0, math.cos(a)]])


# --------------------------------------------------------------------------
# Denticle template (local frame: a = downstream, b = spanwise, c = height)
# --------------------------------------------------------------------------
def _crown_fields(P, Q, p: Denticle):
    phase = np.pi * (p.n_ridges * (Q + 1) / 2 - 0.5)
    groove = np.sin(phase) ** 2 if p.n_ridges > 0 else np.zeros_like(Q)
    p_end = 1.0 - p.tip_notch * groove
    # narrow, rounded leading edge widening to a broad cusped tail (mako-like crown);
    # flat near the tail so the (P, Q) -> (a, b) map never folds
    outline = 0.55 + 0.45 * np.sin(0.5 * np.pi * np.clip(P / 0.7, 0, 1))
    ramp = p.lead_height + (1 - p.lead_height) * np.sin(0.5 * np.pi * np.clip(P / 0.85, 0, 1)) ** 1.2
    ramp *= 1 - 0.25 * np.clip((P - 0.9) / 0.1, 0, 1) ** 2
    prof = np.clip(1 - Q ** 2, 0, 1) ** 0.3
    c = p.height * ramp * prof * (1 - p.ridge_depth * groove)
    return p_end, outline, c


def denticle_template(p: Denticle) -> tuple[np.ndarray, np.ndarray]:
    L, W, H, e = p.length, p.width, p.height, p.embed
    nu, nv = p.grid
    P, Q = np.meshgrid(np.linspace(0, 1, nu), np.linspace(-1, 1, nv), indexing="ij")

    unit = replace(p, length=1.0, width=1.0, height=1.0)
    p_end, outline, c = _crown_fields(P, Q, unit)
    a = L * (P * p_end - 0.5)
    b = 0.5 * W * Q * outline
    c = H * c
    top = np.stack([a, b, c], -1).reshape(-1, 3)
    idx = lambda i, j: i * nv + j
    faces = []
    for i in range(nu - 1):
        for j in range(nv - 1):
            q = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            faces += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
    ring = ([idx(i, 0) for i in range(nu)] + [idx(nu - 1, j) for j in range(1, nv)]
            + [idx(i, nv - 1) for i in range(nu - 2, -1, -1)] + [idx(0, j) for j in range(nv - 2, 0, -1)])
    nt = len(top)
    bot = top[ring].copy()
    bot[:, 2] = -e
    centre = nt + len(ring)
    verts = np.vstack([top, bot, [[0, 0, -e]]])
    n = len(ring)
    for k in range(n):
        t0, t1 = ring[k], ring[(k + 1) % n]
        b0, b1 = nt + k, nt + (k + 1) % n
        faces += [(t0, b1, t1), (t0, b0, b1), (b0, centre, b1)]   # side wall + flat base fan
    tm = trimesh.Trimesh(verts, np.array(faces), process=True)
    if tm.volume < 0:
        tm.invert()
    assert tm.is_watertight and tm.is_winding_consistent, "denticle template not watertight"
    return np.asarray(tm.vertices), np.asarray(tm.faces)


_CROWN_CACHE: dict = {}


def crown_height(a, b, p: Denticle):
    """Height of one crown at local (a, b); 0 outside its footprint. Evaluated by linear
    interpolation on a fine (241 x 241) forward-mapped copy of the template surface, in
    coordinates normalised by length/width so one interpolator serves every scale."""
    import matplotlib.tri as mtri
    key = (p.lead_height, p.n_ridges, p.ridge_depth, p.tip_notch)
    if key not in _CROWN_CACHE:
        n = 241
        P, Q = np.meshgrid(np.linspace(0, 1, n), np.linspace(-1, 1, n), indexing="ij")
        unit = replace(p, length=1.0, width=1.0, height=1.0)
        p_end, outline, c = _crown_fields(P, Q, unit)
        A, B = P * p_end - 0.5, 0.5 * Q * outline
        g = np.arange(n * n).reshape(n, n)
        q = np.stack([g[:-1, :-1], g[1:, :-1], g[1:, 1:], g[:-1, 1:]], -1).reshape(-1, 4)
        tri = mtri.Triangulation(A.ravel(), B.ravel(), np.r_[q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
        _CROWN_CACHE[key] = mtri.LinearTriInterpolator(tri, c.ravel())
    an, bn = np.asarray(a) / p.length, np.asarray(b) / p.width
    out = np.zeros(an.shape)
    near = (np.abs(an) <= 0.51) & (np.abs(bn) <= 0.51)
    if near.any():
        v = _CROWN_CACHE[key](an[near], bn[near])
        out[near] = np.ma.filled(v, 0.0) * p.height
    return np.maximum(out, 0.0)


def texture_height(u, s, p: Denticle, n_rows: int, n_span: int, arc: float):
    """Union (max) of all crowns of the staggered layout at surface coords
    u (downstream distance) and s (arc length); periodic in s with period n_span*arc."""
    circ = n_span * arc
    h = np.zeros_like(u)
    i0 = np.round(u / p.pitch_stream).astype(int)
    for di in range(-3, 4):
        i = i0 + di
        ok = (i >= 0) & (i < n_rows)
        off = 0.5 * (i % 2)
        a = u - i * p.pitch_stream
        j0 = np.round(s / arc - off)
        for dj in (-1, 0, 1):
            b = s - (j0 + dj + off) * arc
            b = (b + 0.5 * circ) % circ - 0.5 * circ
            h = np.maximum(h, np.where(ok, crown_height(a, b, p), 0.0))
    return h


def textured_liner(r_in, z_top, z_bot, p: Denticle, th0, th1, cell, embed, lift,
                   R=np.eye(3), t=(0, 0, 0)) -> mf.Manifold:
    """Closed solid lining the inside of a tube (axis z, flow -z) over angles [th0, th1]:
    inner face = wall displaced inward by the denticle height field, outer face buried
    `embed` into the wall. `lift` keeps the gaps between crowns just off the wall surface
    so the boolean union has no coplanar faces."""
    n_span = max(3, int(round(2 * math.pi * r_in / p.pitch_span)))
    arc = 2 * math.pi * r_in / n_span
    length = z_top - z_bot
    n_rows = int(math.floor(length / p.pitch_stream)) + 1
    nt = int(math.ceil((th1 - th0) * r_in / cell)) + 1
    nz = int(math.ceil(length / cell)) + 1
    TH, U = np.meshgrid(np.linspace(th0, th1, nt), np.linspace(0, length, nz), indexing="ij")
    H = texture_height(U, (TH % (2 * math.pi)) * r_in, p, n_rows, n_span, arc)
    def ring(r):
        return np.stack([r * np.cos(TH), r * np.sin(TH), z_top - U], -1).reshape(-1, 3)
    inner = ring(r_in - lift - H)
    outer = ring(np.full_like(H, r_in + embed))
    V = np.vstack([inner, outer])
    N = nt * nz
    g = np.arange(N).reshape(nt, nz)
    q = np.stack([g[:-1, :-1], g[1:, :-1], g[1:, 1:], g[:-1, 1:]], -1).reshape(-1, 4)
    F = [q[:, [0, 1, 2]], q[:, [0, 2, 3]], q[:, [0, 2, 1]] + N, q[:, [0, 3, 2]] + N]
    border = np.r_[g[:, 0], g[-1, 1:], g[-2::-1, -1], g[0, -2:0:-1]]
    b1 = np.roll(border, -1)
    F += [np.c_[border, b1 + N, b1], np.c_[border, border + N, b1 + N]]
    tm = trimesh.Trimesh(V @ R.T + t, np.vstack(F), process=False)
    if tm.volume < 0:
        tm.invert()
    return to_manifold(tm.vertices, tm.faces)


def staggered_layout(u0, u1, n_span: int, p: Denticle) -> list[tuple[float, float]]:
    """Imbricated layout: rows every pitch_stream, alternate rows shifted half a pitch
    spanwise, so each crown's raised tail sits between two crowns of the next row."""
    out = []
    rows = int(math.floor((u1 - u0) / p.pitch_stream)) + 1
    for i in range(rows):
        u = u0 + i * p.pitch_stream
        off = 0.5 * (i % 2)
        out += [(u, (j + off)) for j in range(n_span)]       # spanwise index (scaled by caller)
    return out


# --------------------------------------------------------------------------
# Flat patch (true scale)
# --------------------------------------------------------------------------
def build_patch(p: Denticle, size=(2.0, 2.0), base=0.30) -> mf.Manifold:
    lx, ly = size
    lv, lf = denticle_template(p)
    n_span = int(math.ceil(ly / p.pitch_span)) + 2
    parts = []
    for u, j in staggered_layout(-p.length, lx + p.length, n_span, p):
        v = (j - 1) * p.pitch_span
        parts.append(to_manifold(lv + [u, v, 0.0], lf))
    dent = mf.Manifold.batch_boolean(parts, mf.OpType.Add)
    clip = mf.Manifold.cube([lx, ly, base + 1.0]).translate([0, 0, -base])
    slab = mf.Manifold.cube([lx, ly, base]).translate([0, 0, -base])
    return slab + (dent ^ clip)


# --------------------------------------------------------------------------
# Graft body
# --------------------------------------------------------------------------
def disk(r, z, x=0.0, segs=256):
    return mf.Manifold.cylinder(1e-3, r, r, segs).translate([x, 0, z])


def limb_frames(g: Graft):
    """(rotation, start point) for left (-x) and right (+x) limbs. Canonical limb: axis -z from 0."""
    xc = 0.5 * g.d_main - 0.5 * g.d_limb                      # limb centres so the stadium is d_main wide
    z0 = -(g.l_main + g.l_trans) + 2.0                        # limbs start inside the taper
    return [(rot_y(+g.limb_angle), (-xc, 0, z0)), (rot_y(-g.limb_angle), (xc, 0, z0))]


def graft_solid(g: Graft) -> mf.Manifold:
    ro, ri = g.d_main / 2, g.d_main / 2 - g.wall
    lo, li = g.d_limb / 2, g.d_limb / 2 - g.wall
    xc = ro - lo
    zt0, zt1 = -g.l_main, -(g.l_main + g.l_trans)
    s = g.segs
    outer = mf.Manifold.cylinder(g.l_main, ro, ro, s).translate([0, 0, zt0])
    outer += mf.Manifold.batch_hull([disk(ro, zt0, segs=s), disk(lo, zt1 - g.wall, -xc, s // 2),
                                     disk(lo, zt1 - g.wall, xc, s // 2)])
    inner = mf.Manifold.cylinder(g.l_main + 1.0, ri, ri, s).translate([0, 0, zt0])
    inner += mf.Manifold.batch_hull([disk(ri, zt0 + 1e-3, segs=s), disk(li, zt1, -xc, s // 2),
                                     disk(li, zt1, xc, s // 2)])
    for R, t in limb_frames(g):
        outer += mf.Manifold.cylinder(g.l_limb, lo, lo, s // 2).translate([0, 0, -g.l_limb]).transform(affine(R, t))
        inner += mf.Manifold.cylinder(g.l_limb + 1.0, li, li, s // 2).translate([0, 0, -g.l_limb - 0.5]).transform(affine(R, t))
    return outer - inner


def zstent_ring(R_c, z0, amp, n_apex, wire_r, R=np.eye(3), t=(0, 0, 0), n=None, m=10):
    """Closed sinusoidal (Z-stent) wire ring, swept circular section, as a watertight tube."""
    n = n or n_apex * 24
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    C = np.c_[R_c * np.cos(th), R_c * np.sin(th), z0 + amp * np.sin(n_apex * th)]
    T = np.roll(C, -1, 0) - np.roll(C, 1, 0)
    T /= np.linalg.norm(T, axis=1)[:, None]
    N1 = np.c_[np.cos(th), np.sin(th), np.zeros(n)]
    N1 -= (N1 * T).sum(1)[:, None] * T
    N1 /= np.linalg.norm(N1, axis=1)[:, None]
    N2 = np.cross(T, N1)
    ph = np.linspace(0, 2 * np.pi, m, endpoint=False)
    V = (C[:, None, :] + wire_r * (np.cos(ph)[None, :, None] * N1[:, None, :]
                                   + np.sin(ph)[None, :, None] * N2[:, None, :])).reshape(-1, 3)
    F = []
    for i in range(n):
        for k in range(m):
            a, b = i * m + k, i * m + (k + 1) % m
            c, d = ((i + 1) % n) * m + k, ((i + 1) % n) * m + (k + 1) % m
            F += [(a, c, d), (a, d, b)]
    tm = trimesh.Trimesh(V @ R.T + t, np.array(F), process=True)
    if tm.volume < 0:
        tm.invert()
    return to_manifold(tm.vertices, tm.faces)


def stents(g: Graft) -> list[mf.Manifold]:
    w = g.wire_d / 2
    rc_main = g.d_main / 2 + 0.8 * w           # slightly embedded so the union is one body
    rc_limb = g.d_limb / 2 + 0.8 * w
    rings = [zstent_ring(rc_main, z, 4.0, 8, w) for z in (-5.5, -17.5, -29.5, -41.5)]
    for R, t in limb_frames(g):
        rings += [zstent_ring(rc_limb, -z, 3.0, 6, w, R, t) for z in (14.0, 25.5, 37.0, 48.5, 60.0)]
    return rings


def build_graft(g: Graft, p: Denticle | None, cutaway: bool, with_stents: bool,
                cell: float = 0.1) -> mf.Manifold:
    """p = None -> smooth graft. With p, the luminal surface of the main body and both
    limbs is lined with the denticle height field (the transition/flow divider is left smooth)."""
    parts = [graft_solid(g)]
    if p is not None:
        assert cutaway, "the height-field lining is only built for the half-graft cutaway"
        ri, li = g.d_main / 2 - g.wall, g.d_limb / 2 - g.wall
        th0, th1 = (math.pi - 0.05, 2 * math.pi + 0.05) if cutaway else (0.0, 2 * math.pi)
        kw = dict(cell=cell, embed=0.5 * g.wall, lift=0.02 * p.height / 0.4)
        m = 1.0
        parts.append(textured_liner(ri, -m, -g.l_main + m, p, th0, th1, **kw))
        for R, t in limb_frames(g):
            parts.append(textured_liner(li, -10.0, -g.l_limb + m, p, th0, th1, R=R,
                                        t=np.asarray(t, float), **kw))
    if with_stents:
        parts += stents(g)
    body = mf.Manifold.batch_boolean(parts, mf.OpType.Add)
    if cutaway:   # remove the front (+y) half so the luminal texture is visible
        body = body ^ mf.Manifold.cube([400, 200, 400]).translate([-200, -200, -300])
    return body


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------
def render(paths: dict, out_png: str):
    import pyvista as pv
    pv.global_theme.background = "white"
    p = pv.Plotter(off_screen=True, shape=(2, 3), window_size=(2700, 1900))
    graft = "#c9d6df"
    specs = [
        ("smooth", "1. EVAR graft, true dimensions\n26 mm body, 14 mm iliac limbs, Z-stents outside",
         [(170, 210, 10), (0, 0, -64), (0, 0, 1)], 1.0),
        ("cutaway", "2. Front half removed: luminal denticle lining\n(denticles drawn x10 size - illustrative)",
         [(0, 250, -64), (0, 0, -64), (0, 0, 1)], 1.0),
        ("cutaway", "3. Close-up of 2 (flow is downward)\noverlapping crowns, tails pointing downstream",
         [(5, -2, -16), (0, -12.7, -25), (0, 0, 1)], 1.0),
        ("patch", "4. True-scale luminal wall, 2 x 2 mm\n(flow -> +x; 0.3 mm graft wall)",
         [(-0.9, -2.4, 1.6), (1.0, 1.0, 0.0), (0, 0, 1)], 1.0),
        ("patch", "5. True-scale patch, top view (flow -> right)\nstaggered rows: 120 um pitch, 130 um span",
         [(1.0, 1.0, 5), (1.0, 1.0, 0.0), (0, 1, 0)], 2.2),
        ("unit", "6. One denticle, true scale\n200 x 150 x 40 um, 3 riblets, 3-cusp tail",
         [(-0.30, -0.33, 0.22), (0.0, 0.0, 0.01), (0, 0, 1)], 1.0),
    ]
    cache = {}
    for k, (key, title, cpos, zoom) in enumerate(specs):
        p.subplot(k // 3, k % 3)
        mesh = cache.setdefault(key, pv.read(paths[key]))
        p.add_mesh(mesh, color=graft, specular=0.4, specular_power=20)
        p.add_text(title, font_size=11, color="black")
        p.camera_position = cpos
        p.camera.zoom(zoom)
    p.screenshot(out_png)


# --------------------------------------------------------------------------
def save(m: mf.Manifold, name: str) -> str:
    tm = to_trimesh(m)
    assert tm.is_watertight, f"{name} is not watertight"
    path = os.path.join(HERE, name)
    tm.export(path)
    print(f"{name:40s} {len(tm.faces):>9,d} tris  {os.path.getsize(path)/1e6:6.1f} MB  "
          f"extent {np.round(tm.extents, 3)} mm")
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--illustrative-scale", type=float, default=10.0,
                    help="denticle enlargement for the whole-graft cutaway (default 10)")
    ap.add_argument("--patch-size", type=float, default=2.0, help="true-scale patch edge (mm)")
    ap.add_argument("--cell", type=float, default=0.1,
                    help="height-field resolution of the cutaway's denticle lining (mm)")
    ap.add_argument("--no-render", action="store_true")
    a = ap.parse_args()

    d, g = Denticle(), Graft()
    paths = {
        "unit": save(to_manifold(*denticle_template(d)), "denticle_unit_true_scale.stl"),
        "patch": save(build_patch(d, (a.patch_size, a.patch_size), g.wall), "denticle_patch_true_scale.stl"),
        "smooth": save(build_graft(g, None, cutaway=False, with_stents=True), "evar_graft_smooth.stl"),
        "cutaway": save(build_graft(g, d.scaled(a.illustrative_scale), cutaway=True, with_stents=True,
                                    cell=a.cell), "evar_graft_cutaway_illustrative.stl"),
    }
    if not a.no_render:
        render(paths, os.path.join(HERE, "preview.png"))
        print("preview.png written")


if __name__ == "__main__":
    main()
