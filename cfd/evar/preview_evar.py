"""Geometry preview: mid-plane outline, limb cross-sections, ridge path."""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pyvista as pv  # noqa: E402


def preview(geo_dir: str, out_png: str, title: str = ""):
    L = pv.read(f"{geo_dir}/lumen.stl").clip_box([0, 145, 0, 80, -20, 20], invert=False)
    try:
        R = pv.read(f"{geo_dir}/ridge.stl")
    except Exception:
        R = None
    fig, ax = plt.subplots(2, 3, figsize=(16, 7.5))
    for a, z in zip(ax[0], [0.0, 6.0, 10.0]):
        s = L.slice(normal=(0, 0, 1), origin=(0, 0, z))
        a.plot(s.points[:, 0], s.points[:, 1], "k.", ms=0.6)
        if R is not None:
            r = R.slice(normal=(0, 0, 1), origin=(0, 0, z))
            if r.n_points:
                a.plot(r.points[:, 0], r.points[:, 1], "r.", ms=0.6)
        a.set_aspect("equal"); a.set_title(f"section z = {z} mm (y = 0 is the symmetry plane)")
        a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
    for a, x in zip(ax[1], [20.0, 45.0, 100.0]):
        s = L.slice(normal=(1, 0, 0), origin=(x, 0, 0))
        a.plot(s.points[:, 1], s.points[:, 2], "k.", ms=0.6)
        if R is not None:
            r = R.slice(normal=(1, 0, 0), origin=(x, 0, 0))
            if r.n_points:
                a.plot(r.points[:, 1], r.points[:, 2], "r.", ms=0.6)
        a.set_aspect("equal"); a.set_title(f"cross-section x = {x} mm")
        a.set_xlabel("y (mm)"); a.set_ylabel("z (mm)")
    fig.suptitle(title or geo_dir)
    plt.tight_layout()
    plt.savefig(out_png, dpi=90)
    plt.close(fig)


if __name__ == "__main__":
    preview(sys.argv[1], sys.argv[2])
