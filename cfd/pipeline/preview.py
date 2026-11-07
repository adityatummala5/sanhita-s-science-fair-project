"""Visual check of the stent designs: zoomed, unrolled height maps of the
strut luminal surface (colour = how far the surface sticks into the lumen)."""
from __future__ import annotations

import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import trimesh  # noqa: E402

R_LUMINAL = 1.50


def preview(geometry_dir: str, out_png: str, xlim=(1.95, 2.45), slim=(0.70, 1.05)):
    stls = sorted(glob.glob(os.path.join(geometry_dir, "*.stl")))
    stls.sort(key=lambda p: (not os.path.basename(p).startswith("bare"),
                             not os.path.basename(p).startswith("student"), p))
    if not stls:
        return
    n = len(stls)
    cols = min(n, 4)
    rows = (n + cols - 1) // cols
    fig, axs = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.4 * rows), squeeze=False)
    for ax, path in zip(axs.ravel(), stls):
        m = trimesh.load(path)
        v = m.vertices
        s = R_LUMINAL * np.arctan2(v[:, 2], v[:, 1])
        h = (R_LUMINAL - np.hypot(v[:, 1], v[:, 2])) * 1e3   # um above strut luminal surface
        sel = (v[:, 0] > xlim[0]) & (v[:, 0] < xlim[1]) & (s > slim[0]) & (s < slim[1]) & (h > -1)
        order = np.argsort(h[sel])
        sc = ax.scatter(v[sel, 0][order], s[sel][order], c=h[sel][order], s=0.6, cmap="viridis",
                        vmin=0, vmax=65, rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlim(*xlim)
        ax.set_ylim(*slim)
        ax.set_title(os.path.splitext(os.path.basename(path))[0], fontsize=10)
        ax.set_xlabel("axial (mm)  → flow", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.colorbar(sc, ax=axs, shrink=0.6, label="height into lumen (µm)")
    fig.suptitle("Strut luminal surface, unrolled (same strut segment in every design)")
    plt.savefig(out_png, dpi=120)
    plt.close(fig)
