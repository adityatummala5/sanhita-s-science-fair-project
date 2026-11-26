"""Periodic surface textures for the wall-patch model (heightfields, in mm).

Every texture is a periodic heightfield y = h(x, z) over the graft surface:
x = local principal flow direction, z = spanwise, y = into the lumen.
y = 0 is the nominal graft surface.

placement
    'protrude'  features stand on the nominal surface (base at y = 0, tops at +H)
    'flush'     feature tops are level with the nominal surface (tops at 0, the
                field between features is cut down to -H): a height-neutral
                texture embossed INTO the graft coating
    'engraved'  the feature shape is cut into the surface as a pit (y = -h)

kinds
    smooth      untextured graft (ePTFE-like control)
    pet         idealized woven-polyester texture: transverse yarn ridges with a
                weaker spanwise weave modulation (clinical comparator)
    denticles   shield-shaped crowns with streamwise riblets on a staggered
                (diamond) lattice; crown shape follows the coronary pipeline
                (pipeline/geometry.py) but with softened side walls so a
                surface-following mesh stays valid
    riblets     continuous streamwise riblets (triangular or scalloped) - the
                classic shark-skin drag-reduction surface
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np


@dataclass
class Texture:
    name: str
    kind: str = "smooth"                 # smooth | pet | denticles | riblets
    placement: str = "protrude"          # protrude | flush | engraved
    # denticles (um)
    length_um: float = 100.0
    width_um: float = 70.0
    height_um: float = 43.0
    n_ridges: int = 3
    ridge_depth: float = 0.35
    spacing_um: float = 85.0             # streamwise centre-to-centre
    row_spacing_um: float = 70.0         # spanwise distance between staggered rows
    orientation: str = "flow"            # flow | reverse | transverse
    edge_exp: float = 0.8                # side-wall softness (0.5 = vertical walls, as in the CAD)
    # riblets (um)
    riblet_spacing_um: float = 60.0
    riblet_shape: str = "triangle"       # triangle | scallop
    # pet (um)
    pet_wavelength_um: float = 200.0
    pet_span_wavelength_um: float = 400.0
    description: str = ""
    target_size_mm: float = 0.5          # platelet-simulation domain (tiled unit cells)
    mesh_dx_um: float = 0.0              # per-design override of the patch cell size (0 = default)

    # ------------------------------------------------------------------
    @property
    def H(self) -> float:
        return self.height_um / 1e3

    def periods(self) -> tuple[float, float]:
        """Fundamental periods (mm) in x and z."""
        if self.kind == "denticles":
            sx, sz = self.spacing_um / 1e3, 2 * self.row_spacing_um / 1e3
            if self.orientation == "transverse":
                return sx, sz      # lattice is not rotated; only the crowns are
            return sx, sz
        if self.kind == "riblets":
            return 0.03, self.riblet_spacing_um / 1e3       # invariant along x: any x-length works
        if self.kind == "pet":
            return self.pet_wavelength_um / 1e3, self.pet_span_wavelength_um / 1e3
        return 0.03, 0.03

    def domain(self, cfd: bool = True) -> tuple[float, float]:
        """cfd=True: one texture period (the flow is periodic, so one unit cell gives
        exact per-area statistics). cfd=False: the tiled platelet-simulation domain."""
        px, pz = self.periods()
        if cfd:
            return px, pz
        nx = max(1, round(self.target_size_mm / px))
        nz = max(1, round(self.target_size_mm / pz))
        return nx * px, nz * pz

    # ------------------------------------------------------------------
    def _crown(self, a, b):
        """Height (mm) of one crown at local coords a (streamwise), b (spanwise), mm."""
        L, W, H = self.length_um / 1e3, self.width_um / 1e3, self.H
        P = a / L + 0.5
        inside_p = (P > 0) & (P < 1)
        Pc = np.clip(P, 1e-6, 1 - 1e-6)
        width = np.maximum(0.06, np.sin(np.pi * Pc ** 0.8) ** 0.55)
        Q = b / (0.5 * W * width)
        inside = inside_p & (np.abs(Q) < 1)
        env_s = np.sin(np.pi * Pc ** 0.85) ** self.edge_exp * (0.55 + 0.45 * Pc)
        env_s /= (np.sin(np.pi * np.linspace(0.01, 0.99, 400) ** 0.85) ** self.edge_exp
                  * (0.55 + 0.45 * np.linspace(0.01, 0.99, 400))).max()
        env_q = np.clip(1 - Q ** 2, 0, 1) ** self.edge_exp
        if self.n_ridges > 0 and self.ridge_depth > 0:
            ridge = 1.0 - self.ridge_depth * np.sin(np.pi * (self.n_ridges * (Q + 1.0) / 2.0 - 0.5)) ** 2
        else:
            ridge = 1.0
        return np.where(inside, H * env_s * env_q * ridge, 0.0)

    def raw(self, X, Z):
        """Feature heightfield >= 0 (mm), before placement."""
        if self.kind == "smooth":
            return np.zeros_like(X)
        Lx, Lz = self.domain()
        X = np.where(np.isclose(X, Lx), X, np.mod(X, Lx)) if np.ndim(X) else X
        Z = np.where(np.isclose(Z, Lz), Z, np.mod(Z, Lz)) if np.ndim(Z) else Z
        if self.kind == "pet":
            lx, lz = self.pet_wavelength_um / 1e3, self.pet_span_wavelength_um / 1e3
            yarn = np.cos(np.pi * X / lx) ** 2                       # transverse ridges
            weave = 0.8 + 0.2 * np.cos(2 * np.pi * Z / lz)
            return self.H * yarn * weave
        if self.kind == "riblets":
            s = self.riblet_spacing_um / 1e3
            u = (Z / s) % 1.0
            if self.riblet_shape == "scallop":
                prof = 1.0 - np.sqrt(np.clip(1 - (2 * u - 1) ** 2, 0, 1))
            else:
                prof = np.abs(2 * u - 1)                              # peaks at u = 0, 1
            return self.H * prof
        if self.kind == "denticles":
            sx, sz = self.spacing_um / 1e3, self.row_spacing_um / 1e3
            out = np.zeros_like(X)
            nrow = int(round(Lz / sz))
            ncol = int(round(Lx / sx))
            reach = 0.6 * max(self.length_um, self.width_um) / 1e3
            # explicit lattice images (crowns can be larger than the unit cell)
            for j in range(-3, nrow + 3):
                zc = (j + 0.5) * sz
                xoff = (j % 2) * 0.5 * sx
                for i in range(-3, ncol + 3):
                    xc = (i + 0.5) * sx + xoff
                    dx = X - xc
                    dz = Z - zc
                    near = (np.abs(dx) < reach) & (np.abs(dz) < reach)
                    if not near.any():
                        continue
                    if self.orientation == "flow":
                        a, b = dx, dz
                    elif self.orientation == "reverse":
                        a, b = -dx, dz
                    elif self.orientation == "transverse":
                        a, b = dz, -dx
                    else:
                        raise ValueError(self.orientation)
                    c = np.zeros_like(X)
                    c[near] = self._crown(a[near], b[near])
                    out = np.maximum(out, c)
            return out
        raise ValueError(self.kind)

    def height(self, X, Z):
        r = self.raw(X, Z)
        if self.placement == "protrude":
            return r
        if self.placement == "flush":
            return r - self.H
        if self.placement == "engraved":
            return -r
        raise ValueError(self.placement)

    def summary(self, n: int = 400) -> dict:
        Lx, Lz = self.domain()
        x = (np.arange(n) + 0.5) * Lx / n
        z = (np.arange(n) + 0.5) * Lz / n
        X, Z = np.meshgrid(x, z, indexing="ij")
        h = self.height(X, Z)
        gx, gz = np.gradient(h, Lx / n, Lz / n)
        wet = np.mean(np.sqrt(1 + gx ** 2 + gz ** 2))
        return {"domain_mm": [round(Lx, 4), round(Lz, 4)], "h_min_um": round(1e3 * h.min(), 2),
                "h_max_um": round(1e3 * h.max(), 2), "h_mean_um": round(1e3 * h.mean(), 2),
                "wetted_area_ratio": round(float(wet), 4),
                "coverage_pct": round(100 * float((self.raw(X, Z) > 1e-4).mean()), 1),
                "max_slope_deg": round(float(np.degrees(np.arctan(np.sqrt(gx ** 2 + gz ** 2).max()))), 1)}


def from_dict(d: dict) -> Texture:
    return Texture(**d)
