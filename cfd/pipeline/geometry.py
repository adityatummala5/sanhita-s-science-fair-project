"""Geometry stage: build clean, watertight stent surfaces and the artery lumen.

All geometry is in millimetres (the units of the CAD files). The mesh is
converted to metres after meshing (see case.py).

Coordinate system (same as the repo STL files):
    x      = stent axis = flow direction
    (y, z) = cross-section; theta = atan2(z, y)
    Strut luminal (inner) surface at r = 1.50 mm, abluminal (outer) at r = 1.60 mm.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import manifold3d as mf
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from shapely import contains_xy
from shapely.geometry import Polygon
from shapely.ops import unary_union

R_LUMINAL = 1.50   # strut inner radius (mm)
R_ABLUMINAL = 1.60  # strut outer radius (mm)

# Map extrusion frame (u, v, w) -> STL frame (x=w, y=u, z=v)
_W_TO_X = np.array([[0.0, 0.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0]])


# --------------------------------------------------------------------------
# trimesh <-> manifold helpers
# --------------------------------------------------------------------------
def to_manifold(tm: trimesh.Trimesh) -> mf.Manifold:
    mesh = mf.Mesh(vert_properties=np.asarray(tm.vertices, dtype=np.float32),
                   tri_verts=np.asarray(tm.faces, dtype=np.uint32))
    m = mf.Manifold(mesh)
    if m.status() != mf.Error.NoError:
        raise ValueError(f"mesh is not a valid closed manifold: {m.status()}")
    return m


def to_trimesh(m: mf.Manifold) -> trimesh.Trimesh:
    mesh = m.to_mesh()
    return trimesh.Trimesh(vertices=np.asarray(mesh.vert_properties)[:, :3],
                           faces=np.asarray(mesh.tri_verts), process=False)


def load_solid(path: str) -> mf.Manifold:
    """Load an STL made of (possibly overlapping, possibly inside-out) closed
    bodies and return their boolean union as one valid solid."""
    tm = trimesh.load(path, process=True)
    bodies = tm.split(only_watertight=False)
    parts = []
    for b in bodies:
        if not b.is_watertight:
            raise ValueError(f"{path}: a body is not watertight; cannot repair automatically")
        if b.volume < 0:          # repo stent files are inside-out
            b.invert()
        parts.append(to_manifold(b))
    return mf.Manifold.batch_boolean(parts, mf.OpType.Add)


def drop_fragments(m: mf.Manifold, min_volume: float) -> tuple[mf.Manifold, int]:
    """Remove disconnected pieces smaller than min_volume (mm^3)."""
    pieces = m.decompose()
    keep = [p for p in pieces if p.volume() >= min_volume]
    dropped = len(pieces) - len(keep)
    return mf.Manifold.batch_boolean(keep, mf.OpType.Add), dropped


# --------------------------------------------------------------------------
# Sector window (the region of the stent that is simulated)
# --------------------------------------------------------------------------
@dataclass
class Window:
    x0: float = 1.253          # matches stent_denticle_CFD_sector.stl
    x1: float = 3.760
    theta0_deg: float = -20.0
    theta1_deg: float = 100.0

    def solid(self, r_out: float = 2.5, n_arc: int = 480) -> mf.Manifold:
        t = np.radians(np.linspace(self.theta0_deg, self.theta1_deg, n_arc))
        poly = [(0.0, 0.0)] + [(r_out * math.cos(a), r_out * math.sin(a)) for a in t]
        cs = mf.CrossSection([poly])
        ext = mf.Manifold.extrude(cs, self.x1 - self.x0)
        return ext.transform(_W_TO_X).translate((self.x0, 0.0, 0.0))


# --------------------------------------------------------------------------
# Parametric shark-denticle generator
# --------------------------------------------------------------------------
@dataclass
class DenticleParams:
    """Dimensions in micrometres. Defaults are measured from Sanhita's
    stent_denticle_CFD_sector.stl (footprint ~100 x 73 um, height ~43 um,
    spacing ~81 um along the strut)."""
    length_um: float = 100.0     # streamwise footprint length
    width_um: float = 70.0       # spanwise footprint width
    height_um: float = 43.0      # crown height above the strut surface
    n_ridges: int = 3            # longitudinal riblets on the crown
    ridge_depth: float = 0.35    # groove depth as a fraction of local height (0 = smooth)
    spacing_um: float = 85.0     # min centre-to-centre spacing along struts
    orientation: str = "strut"   # 'strut' (along strut), 'flow' (axial) or 'reverse'
    embed_um: float = 5.0        # how far the base is sunk into the strut (for a clean union)
    grid: tuple = (28, 18)       # surface resolution of each denticle (streamwise, spanwise)


def denticle_local(p: DenticleParams) -> tuple[np.ndarray, np.ndarray]:
    """Closed triangle mesh of one denticle in local coords (a, b, c) in mm:
    a = streamwise, b = spanwise, c = height (outward from strut surface).
    Shield-shaped crown, blunt upstream edge, pointed downstream tip, ridges."""
    L, W, H = p.length_um / 1e3, p.width_um / 1e3, p.height_um / 1e3
    e = p.embed_um / 1e3
    nu, nv = p.grid
    P, Q = np.meshgrid(np.linspace(0.0, 1.0, nu), np.linspace(-1.0, 1.0, nv), indexing="ij")

    width = np.maximum(0.06, np.sin(np.pi * np.clip(P, 0, 1) ** 0.8) ** 0.55)  # shield outline
    env_s = np.sin(np.pi * P ** 0.85) ** 0.5 * (0.55 + 0.45 * P)                 # rises downstream
    env_s /= env_s.max()
    env_q = np.sqrt(np.clip(1.0 - Q ** 2, 0.0, 1.0))
    if p.n_ridges > 0 and p.ridge_depth > 0:
        ridge = 1.0 - p.ridge_depth * np.sin(np.pi * (p.n_ridges * (Q + 1.0) / 2.0 - 0.5)) ** 2
    else:
        ridge = np.ones_like(Q)

    a = L * (P - 0.5)
    b = 0.5 * W * Q * width
    c_top = H * env_s * env_q * ridge
    top = np.stack([a, b, c_top], -1).reshape(-1, 3)
    bot = np.stack([a, b, np.full_like(a, -e)], -1).reshape(-1, 3)
    verts = np.vstack([top, bot])
    nt = nu * nv

    def idx(i, j, base=0):
        return base + i * nv + j

    faces = []
    for i in range(nu - 1):
        for j in range(nv - 1):
            a0, a1, a2, a3 = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            faces += [(a0, a1, a2), (a0, a2, a3)]                      # top (normal +c)
            b0, b1, b2, b3 = (k + nt for k in (a0, a1, a2, a3))
            faces += [(b0, b2, b1), (b0, b3, b2)]                      # bottom (normal -c)
    ring = ([idx(i, 0) for i in range(nu)] + [idx(nu - 1, j) for j in range(1, nv)]
            + [idx(i, nv - 1) for i in range(nu - 2, -1, -1)] + [idx(0, j) for j in range(nv - 2, 0, -1)])
    for k in range(len(ring)):
        t0, t1 = ring[k], ring[(k + 1) % len(ring)]
        faces += [(t0, t1 + nt, t1), (t0, t0 + nt, t1 + nt)]
    tm = trimesh.Trimesh(verts, np.array(faces), process=True)
    if tm.volume < 0:
        tm.invert()
    assert tm.is_watertight, "denticle template is not watertight"
    return np.asarray(tm.vertices), np.asarray(tm.faces)


def _luminal_polygon(stent: trimesh.Trimesh, r_luminal: float | None = None) -> Polygon:
    """Unrolled (x, s = R*theta) outline of the strut luminal surface (at r_luminal, default
    R_LUMINAL; pass the inner radius for thickened struts)."""
    r_lum = R_LUMINAL if r_luminal is None else r_luminal
    fc = stent.triangles_center
    rc = np.hypot(fc[:, 1], fc[:, 2])
    n = stent.face_normals
    radial = (n[:, 1] * fc[:, 1] + n[:, 2] * fc[:, 2]) / rc
    sel = (np.abs(rc - r_lum) < 2e-3) & (radial < -0.98)
    tri = stent.vertices[stent.faces[sel]]
    s = R_LUMINAL * np.arctan2(tri[..., 2], tri[..., 1])
    polys = [Polygon(np.c_[tri[k, :, 0], s[k]]) for k in range(len(tri))]
    return unary_union([p.buffer(1e-5) for p in polys if p.area > 0])


def place_denticles(stent: trimesh.Trimesh, p: DenticleParams, win: Window,
                    grid_step: float = 0.002) -> np.ndarray:
    """Return (N, 3) array of (x, s, phi): denticle centres on the strut
    centreline (unrolled coords, mm) and their orientation angle phi (rad)."""
    poly = _luminal_polygon(stent)
    minx, mins, maxx, maxs = poly.bounds
    X, S = np.meshgrid(np.arange(minx, maxx, grid_step), np.arange(mins, maxs, grid_step), indexing="ij")
    X, S = X.ravel(), S.ravel()
    inside = contains_xy(poly, X, S)
    X, S = X[inside], S[inside]
    # distance to the strut edge; keep the centreline band
    bnd = np.asarray(poly.boundary.segmentize(grid_step / 2).coords if poly.geom_type == "Polygon"
                     else np.vstack([np.asarray(g.boundary.segmentize(grid_step / 2).coords) for g in poly.geoms]))
    d, _ = cKDTree(bnd).query(np.c_[X, S])
    band = d >= 0.8 * np.percentile(d, 99)
    C = np.c_[X[band], S[band]]
    # keep only centres inside the window (denticles straddling the edge get clipped later)
    s0, s1 = R_LUMINAL * math.radians(win.theta0_deg), R_LUMINAL * math.radians(win.theta1_deg)
    C = C[(C[:, 0] > win.x0) & (C[:, 0] < win.x1) & (C[:, 1] > s0) & (C[:, 1] < s1)]
    # greedy Poisson-disk selection along the band
    sp = p.spacing_um / 1e3
    tree = cKDTree(C)
    order = np.lexsort((C[:, 1], C[:, 0]))
    taken = np.zeros(len(C), bool)
    blocked = np.zeros(len(C), bool)
    for k in order:
        if blocked[k]:
            continue
        taken[k] = True
        blocked[tree.query_ball_point(C[k], sp)] = True
    centres = C[taken]
    # local strut direction from centreline points nearby (PCA)
    out = []
    for cx, cs in centres:
        nb = C[tree.query_ball_point((cx, cs), 0.06)]
        if p.orientation == "flow" or len(nb) < 3:
            phi = 0.0
        else:
            w, V = np.linalg.eigh(np.cov((nb - nb.mean(0)).T))
            dvec = V[:, np.argmax(w)]
            if dvec[0] < 0:
                dvec = -dvec               # always point downstream
            phi = math.atan2(dvec[1], dvec[0])
        if p.orientation == "reverse":
            phi += math.pi
        out.append((cx, cs, phi))
    return np.array(out)


def build_denticles(placements: np.ndarray, p: DenticleParams) -> list[mf.Manifold]:
    lv, lf = denticle_local(p)
    a, b, c = lv[:, 0], lv[:, 1], lv[:, 2]
    parts = []
    for cx, cs, phi in placements:
        x = cx + a * math.cos(phi) - b * math.sin(phi)
        s = cs + a * math.sin(phi) + b * math.cos(phi)
        r = R_LUMINAL - c                    # crown points inward, into the lumen
        th = s / R_LUMINAL
        v = np.c_[x, r * np.cos(th), r * np.sin(th)]
        parts.append(to_manifold(trimesh.Trimesh(v, lf, process=False)))
    return parts


# --------------------------------------------------------------------------
# High-level builders
# --------------------------------------------------------------------------
@dataclass
class StentResult:
    mesh: trimesh.Trimesh
    info: dict = field(default_factory=dict)


def bare_sector(stent_stl: str, win: Window) -> StentResult:
    full = load_solid(stent_stl)
    cut = full ^ win.solid()            # boolean intersection
    cut, dropped = drop_fragments(cut, 1e-3)
    tm = to_trimesh(cut)
    return StentResult(tm, {"source": stent_stl, "fragments_removed": dropped})


def thickened_sector(bare: StentResult, added_um: float) -> StentResult:
    """Height-matched control: the same bare stent with every strut made
    `added_um` thicker, by moving the luminal face radially inward. The
    abluminal face (r = 1.60 mm, in the artery wall) and strut paths are unchanged."""
    v = bare.mesh.vertices.copy()
    r = np.hypot(v[:, 1], v[:, 2])
    t = R_ABLUMINAL - R_LUMINAL
    h = added_um / 1e3
    r_new = R_ABLUMINAL - (R_ABLUMINAL - r) * (t + h) / t
    v[:, 1:] *= (r_new / r)[:, None]
    tm = trimesh.Trimesh(v, bare.mesh.faces.copy(), process=False)
    assert tm.is_watertight
    return StentResult(tm, {"strut_thickness_um": 1e3 * t + added_um, "added_um": added_um})


def denticle_sector(bare: StentResult, p: DenticleParams, win: Window) -> StentResult:
    pl = place_denticles(bare.mesh, p, win)
    dent = build_denticles(pl, p)
    solid = mf.Manifold.batch_boolean([to_manifold(bare.mesh)] + dent, mf.OpType.Add)
    solid = solid ^ win.solid()
    solid, dropped = drop_fragments(solid, 1e-3)
    tm = to_trimesh(solid)
    return StentResult(tm, {"n_denticles": int(len(pl)), "fragments_removed": dropped,
                            "params": p.__dict__})


def imported_sector(stl: str) -> StentResult:
    """Use an existing textured-sector STL as-is (after cleaning)."""
    tm = trimesh.load(stl, process=True)
    bodies = tm.split(only_watertight=False)
    keep, dropped = [], 0
    for b in bodies:
        if b.is_watertight and abs(b.volume) >= 1e-3:
            if b.volume < 0:
                b.invert()
            keep.append(to_manifold(b))
        else:
            dropped += 1
    solid = mf.Manifold.batch_boolean(keep, mf.OpType.Add)
    return StentResult(to_trimesh(solid), {"source": stl, "fragments_removed": dropped})


def lumen_surface(win: Window, x_start: float, x_end: float, radius: float,
                  n_arc: int = 241, n_ax: int = 60, n_rad: int = 30) -> dict[str, trimesh.Trimesh]:
    """Wedge-shaped artery lumen, returned as separate named patches so
    snappyHexMesh can assign boundary types: wall, symA, symB, caps."""
    th = np.radians(np.linspace(win.theta0_deg, win.theta1_deg, n_arc))
    xs = np.linspace(x_start, x_end, n_ax)
    # curved artery wall
    TH, XX = np.meshgrid(th, xs, indexing="ij")
    V = np.c_[XX.ravel(), radius * np.cos(TH).ravel(), radius * np.sin(TH).ravel()]
    F = []
    for i in range(n_arc - 1):
        for j in range(n_ax - 1):
            a0, a1, a2, a3 = i * n_ax + j, (i + 1) * n_ax + j, (i + 1) * n_ax + j + 1, i * n_ax + j + 1
            F += [(a0, a2, a1), (a0, a3, a2)]
    wall = trimesh.Trimesh(V, F, process=True)

    def plane(theta):
        rr = np.linspace(0.0, radius, n_rad)
        R_, X_ = np.meshgrid(rr, xs, indexing="ij")
        v = np.c_[X_.ravel(), (R_ * math.cos(theta)).ravel(), (R_ * math.sin(theta)).ravel()]
        f = []
        for i in range(n_rad - 1):
            for j in range(n_ax - 1):
                a0, a1, a2, a3 = i * n_ax + j, (i + 1) * n_ax + j, (i + 1) * n_ax + j + 1, i * n_ax + j + 1
                f += [(a0, a1, a2), (a0, a2, a3)]
        return trimesh.Trimesh(v, f, process=True)

    def cap(x):
        v = [(x, 0.0, 0.0)] + [(x, radius * math.cos(t), radius * math.sin(t)) for t in th]
        f = [(0, k, k + 1) for k in range(1, n_arc)]
        return trimesh.Trimesh(np.array(v), f, process=True)

    caps = trimesh.util.concatenate([cap(x_start), cap(x_end)])
    return {"wall": wall, "symA": plane(th[0]), "symB": plane(th[-1]), "caps": caps}


def write_multisolid_stl(path: str, patches: dict[str, trimesh.Trimesh]) -> None:
    """ASCII STL with one named 'solid' per region (snappyHexMesh regions)."""
    with open(path, "w") as fh:
        for name, m in patches.items():
            fh.write(f"solid {name}\n")
            for tri, n in zip(m.triangles, m.face_normals):
                fh.write(f" facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n  outer loop\n")
                for v in tri:
                    fh.write(f"   vertex {v[0]:.8e} {v[1]:.8e} {v[2]:.8e}\n")
                fh.write("  endloop\n endfacet\n")
            fh.write(f"endsolid {name}\n")


def describe(tm: trimesh.Trimesh) -> dict:
    v = tm.vertices
    r = np.hypot(v[:, 1], v[:, 2])
    return {"triangles": int(len(tm.faces)), "watertight": bool(tm.is_watertight),
            "volume_mm3": float(tm.volume), "area_mm2": float(tm.area),
            "x_range": [float(v[:, 0].min()), float(v[:, 0].max())],
            "r_min_mm": float(r.min()), "r_max_mm": float(r.max()),
            "max_protrusion_um": float((R_LUMINAL - r.min()) * 1e3)}
