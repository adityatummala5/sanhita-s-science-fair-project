"""Idealised bifurcated EVAR stent-graft lumen, plus an optional helical ridge.

Coordinates in mm. x = flow direction (inlet at x = 0), y = lateral (the two
iliac limbs diverge in the x-y plane), z = antero-posterior.

The full lumen is symmetric about y = 0, so only the half y >= 0 (one limb)
is simulated; the y = 0 plane is a symmetry plane. A helical ridge in the
simulated limb therefore stands for a mirror-image ridge in the other limb.

Lumen = smooth union of
    main body   circular (radius R_body) up to x_taper, then flattening into
                a superellipse that just wraps the two limbs at the divider
                (x_div): the "trouser" shape of a bifurcated graft
    two limbs   capsules of radius R_limb whose axes pass through
                (x_div, +-y_div, 0) at +-theta to the x axis
built as a signed-distance field and triangulated with marching cubes
(pyvista / VTK flying edges). The inlet (x = 0), outlet (x = x_out) and
symmetry (y = 0) planes are faces of the OpenFOAM background block, so the
lumen surface simply runs through them.

The ridge is a closed swept solid: a cosine-bump cross-section of height h
and base width w (= 3 h by default) swept along a helix lying on the limb
wall, with 10 mm height ramps at both ends. Its base is buried 0.6 mm into the
wall so that snappyHexMesh sees a clean intersection.
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import asdict, dataclass

import numpy as np
import pyvista as pv
import trimesh


@dataclass
class EvarGeom:
    R_body: float = 13.0        # main body radius (Ø26 mm)
    R_limb: float = 6.5         # iliac limb radius (Ø13 mm)
    x_taper: float = 25.0       # end of straight main body
    x_div: float = 50.0         # axial position of the flow divider
    y_div: float = 7.0          # limb axis offset from the mid-plane at x_div
    theta_deg: float = 15.0     # limb divergence angle (each side; 30° total)
    blend: float = 4.0          # smooth-union blending length
    limb_s0: float = -15.0      # limb capsule start along its axis (inside the body)
    limb_s1: float = 115.0      # limb capsule end along its axis (beyond the outlet)
    x_out: float = 145.0        # outlet plane

    @property
    def d(self):                # limb axis direction (simulated limb, y > 0)
        t = math.radians(self.theta_deg)
        return np.array([math.cos(t), math.sin(t), 0.0])

    @property
    def A(self):                # limb axis origin (at the divider)
        return np.array([self.x_div, self.y_div, 0.0])

    def limb_point(self, s):
        return self.A + s * self.d


@dataclass
class Ridge:
    height: float = 1.0         # mm
    pitch: float = 78.0         # mm of limb axis per full turn
    width_ratio: float = 3.0    # base width / height
    s_start: float = 12.0       # ridge start along the limb axis (mm from the divider)
    s_end: float = 87.0         # ridge end
    ramp: float = 10.0          # height ramp length at each end
    phase_deg: float = 90.0     # angular start position (90° = +z, anterior)
    handed: int = 1             # +1 / -1: twist direction
    bury: float = 0.6           # how far the base is buried in the wall
    dip: float = 0.1            # foot of the bump dips this far below the wall

    @property
    def width(self):
        return self.width_ratio * self.height


# --------------------------------------------------------------------------
# signed-distance primitives (negative inside)
# --------------------------------------------------------------------------
def sd_capsule(p, a, b, r):
    pa, ba = p - a, b - a
    h = np.clip((pa @ ba) / (ba @ ba), 0.0, 1.0)
    return np.linalg.norm(pa - h[:, None] * ba, axis=1) - r


def smin(a, b, k):
    h = np.clip(k - np.abs(a - b), 0.0, None) / k
    return np.minimum(a, b) - h * h * k * 0.25


def smax(a, b, k):
    return -smin(-a, -b, k)


def _smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3 - 2 * t)


def sd_body(g: EvarGeom, p: np.ndarray) -> np.ndarray:
    """Main body: circle (radius R_body) for x < x_taper, morphing into a
    flattened superellipse (half-width y_div + R_limb, half-height R_limb,
    exponent 2 -> 3) at the divider, then capped just past x_div."""
    x, y, z = p[:, 0], np.abs(p[:, 1]), np.abs(p[:, 2])
    t = _smoothstep((x - g.x_taper) / (g.x_div - g.x_taper))
    ay = g.R_body + t * (g.y_div + g.R_limb - g.R_body)
    az = g.R_body + t * (g.R_limb - g.R_body)
    m = 2.0 + t
    rho = ((y / ay) ** m + (z / az) ** m) ** (1.0 / m)
    f = (rho - 1.0) * np.minimum(ay, az)          # approximate distance
    return smax(f, x - (g.x_div + 2.0), g.blend)


def lumen_sdf(g: EvarGeom, p: np.ndarray) -> np.ndarray:
    d = g.d
    dm = np.array([d[0], -d[1], 0.0])
    Am = np.array([g.x_div, -g.y_div, 0.0])
    limb = sd_capsule(p, g.A + g.limb_s0 * d, g.A + g.limb_s1 * d, g.R_limb)
    limb_m = sd_capsule(p, Am + g.limb_s0 * dm, Am + g.limb_s1 * dm, g.R_limb)
    body = smax(sd_body(g, p), -(p[:, 0] + 20.0), 1.0)
    return smin(body, smin(limb, limb_m, g.blend), g.blend)


def build_lumen(g: EvarGeom, h: float = 0.25) -> pv.PolyData:
    lo = np.array([-2.0, -2.0, -g.R_body - 2.0])
    top_y = g.limb_point(g.limb_s1)[1] + g.R_limb + 2.0
    hi = np.array([g.x_out + 3.0, top_y, g.R_body + 2.0])
    n = np.ceil((hi - lo) / h).astype(int) + 1
    grid = pv.ImageData(dimensions=n, spacing=(h, h, h), origin=lo)
    pts = grid.points
    vals = np.empty(len(pts))
    step = 2_000_000
    for i in range(0, len(pts), step):
        vals[i:i + step] = lumen_sdf(g, pts[i:i + step])
    grid.point_data["sdf"] = vals
    surf = grid.contour([0.0], scalars="sdf", method="flying_edges").triangulate().clean()
    return surf


# --------------------------------------------------------------------------
# helical ridge
# --------------------------------------------------------------------------
def _profile(r: Ridge, n_pts: int = 25):
    """Closed cross-section polygon (b, n): b across the ridge, n = height
    above the wall (negative = buried in the wall). Counter-clockwise."""
    w, h, e = r.width, r.height, r.dip
    b = np.linspace(-w / 2, w / 2, n_pts)
    n = (h + e) * 0.5 * (1 + np.cos(2 * np.pi * b / w)) - e
    bump = np.stack([b, n], 1)[::-1]               # from +w/2 to -w/2 over the crest
    base = np.array([[-w / 2 - 0.2, -r.bury], [w / 2 + 0.2, -r.bury]])
    return np.vstack([bump, base])


def build_ridge(g: EvarGeom, r: Ridge, ds: float = 0.12) -> trimesh.Trimesh:
    d = g.d
    e1 = np.array([-d[1], d[0], 0.0])               # in-plane lateral
    e2 = np.array([0.0, 0.0, 1.0])
    R = g.R_limb
    prof = _profile(r)
    k = 2 * np.pi / r.pitch * r.handed
    phi0 = math.radians(r.phase_deg)
    # arc length of the helix per mm of axis
    stretch = math.sqrt(1 + (k * R) ** 2)
    ns = int(math.ceil((r.s_end - r.s_start) * stretch / ds)) + 1
    s_all = np.linspace(r.s_start, r.s_end, ns)
    rings = []
    for s in s_all:
        phi = phi0 + k * (s - r.s_start)
        rad = math.cos(phi) * e1 + math.sin(phi) * e2
        tang = d + R * k * (-math.sin(phi) * e1 + math.cos(phi) * e2)
        tang /= np.linalg.norm(tang)
        bvec = np.cross(tang, rad)
        ramp = max(min(1.0, (s - r.s_start) / r.ramp, (r.s_end - s) / r.ramp), 0.03)
        c = g.limb_point(s)
        ring = []
        for b, nn in prof:
            nn_eff = nn * ramp if nn > 0 else nn   # only the protruding part is ramped
            q = c + R * rad + b * bvec
            a = (q - g.A) @ d
            rq = q - g.A - a * d
            q = g.A + a * d + (R - nn_eff) * rq / np.linalg.norm(rq)
            ring.append(q)
        rings.append(ring)
    V = np.array(rings).reshape(-1, 3)
    m = len(prof)
    F = []
    for i in range(ns - 1):
        for j in range(m):
            a0, a1 = i * m + j, i * m + (j + 1) % m
            b0, b1 = a0 + m, a1 + m
            F += [[a0, b0, a1], [a1, b0, b1]]
    # end caps (fan from the polygon centroid)
    V = np.vstack([V, V[:m].mean(0), V[-m:].mean(0)])
    c0, c1 = len(V) - 2, len(V) - 1
    last = (ns - 1) * m
    for j in range(m):
        F.append([c0, (j + 1) % m, j])
        F.append([c1, last + j, last + (j + 1) % m])
    mesh = trimesh.Trimesh(V, np.array(F), process=True)
    if mesh.volume < 0:
        mesh.invert()
    trimesh.repair.fix_normals(mesh)
    return mesh


def ridge_stats(g: EvarGeom, r: Ridge, mesh: trimesh.Trimesh) -> dict:
    k = 2 * np.pi / r.pitch
    helix_angle = math.degrees(math.atan(k * g.R_limb))
    return {"height_mm": r.height, "pitch_mm": r.pitch, "pitch_over_D": r.pitch / (2 * g.R_limb),
            "base_width_mm": r.width, "height_over_D_pct": 100 * r.height / (2 * g.R_limb),
            "turns": (r.s_end - r.s_start) / r.pitch, "helix_angle_deg": helix_angle,
            "watertight": bool(mesh.is_watertight), "solid_volume_mm3": float(mesh.volume)}


def write_all(out_dir: str, g: EvarGeom, ridge: Ridge | None, refine_ridge: Ridge | None = None):
    """Write lumen.stl and (optionally) ridge.stl / ridgeRefine.stl in mm."""
    os.makedirs(out_dir, exist_ok=True)
    lumen = build_lumen(g)
    lumen.save(os.path.join(out_dir, "lumen.stl"))
    info = {"geometry": asdict(g), "lumen_triangles": int(lumen.n_cells)}
    if ridge is not None:
        rm = build_ridge(g, ridge)
        rm.export(os.path.join(out_dir, "ridge.stl"))
        info["ridge"] = {**asdict(ridge), **ridge_stats(g, ridge, rm)}
    if refine_ridge is not None:      # smooth control: same mesh refinement, no ridge wall
        rm = build_ridge(g, refine_ridge)
        rm.export(os.path.join(out_dir, "ridgeRefine.stl"))
        info["refinement_path"] = asdict(refine_ridge)
    with open(os.path.join(out_dir, "geometry.json"), "w") as fh:
        json.dump(info, fh, indent=2)
    return info
