"""Idealized EVAR endograft lumen (vessel scale), built as an implicit surface.

Geometry (mm), flow along +x, bifurcation in the x-y plane:
    main body     D_body = 25 mm, from x = -L_body to 0
    transition    0 .. L_trans: the circular body splits into two legs
                  (each cross-section is the smooth union of two circles)
    legs          D_leg = 12 mm  (bifurcated main body outlet)
    flare         D_leg -> D_limb over [x_flare0, x_flare1]
    iliac limbs   D_limb = 16 mm, to x = x_end
    The limbs diverge from the mid-plane with half-angle `half_angle_deg`.
Dimensions follow the Endurant-based idealized model of Frontiers in
Physiology 2021 (doi:10.3389/fphys.2021.795085): 25 mm body inlet, 12 mm
bifurcated-body outlets, 16 mm iliac limbs.

The lumen is symmetric about y = 0 and z = 0, so only one quarter is meshed;
the symmetry planes, inlet and outlet come from the background blockMesh box
(see case.py). This file only has to supply the closed lumen wall surface.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import trimesh
from skimage.measure import marching_cubes


def smoothstep(u):
    u = np.clip(u, 0.0, 1.0)
    return u * u * (3.0 - 2.0 * u)


def smin(a, b, k):
    """Polynomial smooth minimum (blend width k)."""
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b * (1 - h) + a * h - k * h * (1 - h)


@dataclass
class GraftParams:
    D_body: float = 25.0
    D_leg: float = 12.0
    D_limb: float = 16.0
    L_body: float = 60.0
    L_trans: float = 15.0
    x_flare0: float = 55.0
    x_flare1: float = 70.0
    x_end: float = 130.0
    half_angle_deg: float = 8.0
    crotch_blend: float = 1.5        # fillet width at the flow divider (mm)

    @property
    def leg_offset(self) -> float:   # leg centre at the end of the transition
        return 0.5 * self.D_body - 0.5 * self.D_leg - 0.25


def centreline(x: np.ndarray, g: GraftParams) -> tuple[np.ndarray, np.ndarray]:
    """Leg/limb centre offset c(x) >= 0 and radius r(x); C1-smooth in x."""
    xs = np.linspace(-g.L_body - 5, g.x_end + 5, 8001)
    dx = xs[1] - xs[0]
    u = (xs - 0.0) / g.L_trans
    bump = 6.0 * np.clip(u, 0, 1) * (1 - np.clip(u, 0, 1)) / g.L_trans        # d/dx smoothstep
    slope = g.leg_offset * bump + math.tan(math.radians(g.half_angle_deg)) * smoothstep((xs - 0.3 * g.L_trans) / g.L_trans)
    c = np.cumsum(slope) * dx
    c -= c[np.searchsorted(xs, 0.0)]
    c = np.maximum(c, 0.0)
    R_body, R_leg, R_limb = g.D_body / 2, g.D_leg / 2, g.D_limb / 2
    r = R_body - (R_body - R_leg) * smoothstep(xs / g.L_trans)
    r = r + (R_limb - R_leg) * smoothstep((xs - g.x_flare0) / (g.x_flare1 - g.x_flare0))
    return np.interp(x, xs, c), np.interp(x, xs, r)


def lumen_sdf(X, Y, Z, g: GraftParams):
    """Approximate signed distance (negative inside the lumen)."""
    c, r = centreline(X, g)
    d1 = np.hypot(Y - c, Z) - r
    d2 = np.hypot(Y + c, Z) - r
    return smin(d1, d2, g.crotch_blend)


def straight_tube_sdf(X, Y, Z, R):
    return np.hypot(Y, Z) - R


def surface_from_sdf(sdf_fn, x_range, y_range, z_range, h: float) -> trimesh.Trimesh:
    """Marching-cubes surface of sdf = 0 on a grid with spacing h. The grid is
    padded with 'outside' values so the surface is closed."""
    xs = np.arange(x_range[0], x_range[1] + h, h)
    ys = np.arange(y_range[0], y_range[1] + h, h)
    zs = np.arange(z_range[0], z_range[1] + h, h)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    F = sdf_fn(X, Y, Z).astype(np.float32)
    F = np.pad(F, 1, constant_values=10.0)
    verts, faces, _, _ = marching_cubes(F, 0.0, spacing=(h, h, h))
    verts += np.array([xs[0] - h, ys[0] - h, zs[0] - h])
    tm = trimesh.Trimesh(verts, faces, process=True)
    # marching_cubes orients normals toward increasing F (outside) -> outward normals already
    if tm.volume < 0:
        tm.invert()
    return tm


def graft_surface(g: GraftParams, h: float = 0.25, margin: float = 2.0) -> trimesh.Trimesh:
    """Lumen wall for the +y, +z quarter (plus margin across the symmetry planes)."""
    c_end, r_end = centreline(np.array([g.x_end]), g)
    y_max = float(c_end[0] + r_end[0]) + 1.0
    z_max = g.D_body / 2 + 1.0
    return surface_from_sdf(lambda X, Y, Z: lumen_sdf(X, Y, Z, g),
                            (-g.L_body - margin, g.x_end + margin), (-margin, y_max), (-margin, z_max), h)


def tube_surface(R: float, L: float, h: float = 0.25, margin: float = 2.0) -> trimesh.Trimesh:
    return surface_from_sdf(lambda X, Y, Z: straight_tube_sdf(X, Y, Z, R),
                            (-margin, L + margin), (-margin, R + 1.0), (-margin, R + 1.0), h)


def describe(g: GraftParams) -> dict:
    xs = np.array([-g.L_body, 0.0, g.L_trans, g.x_flare0, g.x_flare1, g.x_end])
    c, r = centreline(xs, g)
    return {"params": g.__dict__, "stations_x_mm": xs.tolist(), "centre_offset_mm": c.round(3).tolist(),
            "radius_mm": r.round(3).tolist()}
