#!/usr/bin/env python3
"""Platelet adhesion/aggregation: denticle-coated stents vs height-matched bare stents.

Pairs (textured design  vs  bare stent whose struts are as thick as strut + denticle):
    student_sector (42.7 um denticles)  vs  bare_thick143 (143 um struts)
    D4_dense       (41.8 um, dense)     vs  bare_thick143
    D3_short       (24.3 um)            vs  bare_thick124 (124 um struts)
plus the standard 100 um bare stent as a reference.

Every design gets: pulsatile Carreau CFD (coarse mesh, 2 heartbeats) and the same
platelet simulation (pipeline/platelets.py, identical parameters and random seed),
exported as GIFs + data to results/platelet_comparison/.

Stages (run from cfd/):
  python3 compare_platelets.py --stage cfd        # OpenFOAM runs that are missing
  python3 compare_platelets.py --stage platelets  # platelet runs (waits for each CFD case)
  python3 compare_platelets.py --stage summary    # tables, comparison figures and GIFs
  python3 compare_platelets.py                    # all three, one after another
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "pipeline"))

import case as fc        # noqa: E402

PAIRS = [("student_sector", "bare_thick143"), ("D4_dense", "bare_thick143"), ("D3_short", "bare_thick124")]
DESIGNS = ["bare", "bare_thick124", "bare_thick143", "student_sector", "D3_short", "D4_dense"]   # CFD finishing order
LABEL = {"bare": "bare, 100 µm struts", "bare_thick143": "bare, 143 µm struts", "bare_thick124": "bare, 124 µm struts",
         "student_sector": "Sanhita's denticles (43 µm)", "D4_dense": "D4 dense denticles (42 µm)",
         "D3_short": "D3 short denticles (24 µm)"}
RUN_ROOT = os.path.join(HERE, "runs", "coarse-pulsatile")
OUT = os.path.join(HERE, "results", "platelet_comparison")


def setup():
    cfg = json.load(open(os.path.join(HERE, "config/experiment.json")))
    w = cfg["window"]
    dom = fc.Domain(x0=w["x0"], x1=w["x1"], theta0_deg=w["theta0_deg"], theta1_deg=w["theta1_deg"],
                    **cfg.get("domain", {}))
    return cfg, dom


def stage_cfd(np_):
    need = [d for d in DESIGNS if not os.path.exists(os.path.join(RUN_ROOT, d, "metrics.json"))]
    if need:
        subprocess.run([sys.executable, os.path.join(HERE, "run.py"), "--preset", "coarse", "--mode", "pulsatile",
                        "--only", ",".join(need), "--np", str(np_), "--cycles", "2", "--snapshots", "20",
                        "--keep-snapshots", "--skip-existing"], check=True)


def stage_platelets(cfg, dom, wait):
    import animate
    import platelets as pl
    period = cfg["inlet"]["period_s"]
    for d in DESIGNS:
        out = os.path.join(OUT, d)
        if os.path.exists(os.path.join(out, "summary.json")):
            print(f"{d}: platelet results exist, skipping")
            continue
        cdir = os.path.join(RUN_ROOT, d)
        while wait and not os.path.exists(os.path.join(cdir, "metrics.json")):
            time.sleep(30)
        print(f"=== {d}: flow export + platelet simulation")
        animate.export_flow(cdir, out, dom, cfg["physics"]["rho"], period, LABEL[d])
        json.dump(json.load(open(os.path.join(cdir, "metrics.json"))), open(os.path.join(out, "flow_metrics.json"), "w"),
                  indent=2)
        times = animate.snapshot_times(cdir, period)
        pl.simulate(cdir, times, period, dom, os.path.join(HERE, "geometry_out", f"{d}.stl"), out,
                    pl.PlateletParams(n_cycles=2.0), title=LABEL[d])
        pl.growth_gif(os.path.join(out, "adhered_platelets.csv"), dom, os.path.join(HERE, "geometry_out", f"{d}.stl"),
                      os.path.join(out, "deposition_growth.gif"), 2 * period)


# --------------------------------------------------------------------------
def _load(d):
    s = json.load(open(os.path.join(OUT, d, "summary.json")))
    f = json.load(open(os.path.join(OUT, d, "flow_metrics.json")))
    g = json.load(open(os.path.join(HERE, "geometry_out", f"{d}.json")))
    rows = list(csv.DictReader(open(os.path.join(OUT, d, "adhered_platelets.csv"))))
    w = s["parcel_weight_platelets"]
    n_par = len(rows)
    total = n_par * w
    return {
        "design": d, "label": LABEL[d],
        "max_height_above_100um_strut_um": round(g["max_protrusion_um"], 1),
        "strut_volume_mm3": round(g["volume_mm3"], 4),
        "wetted_stent_area_mm2": round(s["stent_wetted_area_mm2"], 3),
        "platelets_total_deposited": round(total),
        "platelets_total_counting_error": round(math.sqrt(max(n_par, 1)) * w),
        "platelets_adhered_to_stent": round(s["platelets_adhered_to_stent"]),
        "platelets_aggregated": round(s["platelets_aggregated"]),
        "aggregated_fraction_pct": round(100 * s["platelets_aggregated"] / max(total, 1), 1),
        "deposit_density_per_mm2": round(total / s["stent_wetted_area_mm2"]),
        "deposits_on_strut_sides_pct": round(100 * s["platelets_adhered_strut_sides"] / max(s["platelets_adhered_to_stent"], 1), 1),
        "platelets_passing_per_s": round(s["platelets_injected"] / s["simulated_time_s"]),
        "deposited_per_million_passing": round(1e6 * total / s["platelets_injected"], 1),
        "parcel_weight": round(w, 2),
        "wall_TAWSS_Pa": round(f["wall_stented_TAWSS_mean_Pa"], 3),
        "wall_lowWSS_area_pct": round(f["wall_stented_lowWSS_area_pct"], 1),
        "wall_RRT_perPa": round(f["wall_stented_RRT_mean_perPa"], 2),
        "stent_TAWSS_Pa": round(f["stent_TAWSS_mean_Pa"], 3),
        "stent_lowWSS_area_pct": round(f["stent_lowWSS_area_pct"], 1),
        "recirc_volume_mm3": round(f["recirc_volume_mm3"], 5),
        "pressure_drop_Pa": round(f["pressure_drop_Pa"], 2),
    }


def stage_summary(cfg, dom):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    rows = {d: _load(d) for d in DESIGNS if os.path.exists(os.path.join(OUT, d, "summary.json"))}
    with open(os.path.join(OUT, "comparison.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(next(iter(rows.values()))))
        w.writeheader()
        for d in DESIGNS:
            if d in rows:
                w.writerow(rows[d])
    pair_rows = []
    for tex, ctl in PAIRS:
        if tex in rows and ctl in rows:
            a, b = rows[tex], rows[ctl]
            pr = {"textured": tex, "control": ctl}
            for k in ("platelets_total_deposited", "platelets_adhered_to_stent", "platelets_aggregated",
                      "deposit_density_per_mm2", "wall_TAWSS_Pa", "wall_lowWSS_area_pct", "stent_TAWSS_Pa",
                      "stent_lowWSS_area_pct", "recirc_volume_mm3", "pressure_drop_Pa"):
                pr[f"{k}_textured"] = a[k]
                pr[f"{k}_control"] = b[k]
                pr[f"{k}_change_pct"] = round(100 * (a[k] - b[k]) / b[k], 1) if b[k] else ""
            pair_rows.append(pr)
    with open(os.path.join(OUT, "pairwise_textured_vs_control.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pair_rows[0]))
        w.writeheader()
        w.writerows(pair_rows)

    # ---- bar chart ---------------------------------------------------------
    order = [d for d in ["bare", "bare_thick124", "D3_short", "bare_thick143", "student_sector", "D4_dense"] if d in rows]
    col = {"bare": "0.6", "bare_thick124": "0.35", "bare_thick143": "0.2",
           "D3_short": "tab:green", "student_sector": "tab:orange", "D4_dense": "tab:blue"}
    panels = [("platelets_total_deposited", "platelets deposited in 1.8 s\n(lower = better)", True),
              ("deposit_density_per_mm2", "deposits per mm² of stent\n(lower = better)", False),
              ("platelets_aggregated", "platelets aggregated\n(lower = better)", False),
              ("wall_lowWSS_area_pct", "artery wall area, TAWSS < 0.4 Pa (%)\n(lower = better)", False),
              ("stent_lowWSS_area_pct", "stent area, TAWSS < 0.4 Pa (%)\n(lower = better)", False),
              ("wall_TAWSS_Pa", "artery wall TAWSS (Pa)\n(higher = better)", False)]
    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    for ax, (k, title, err) in zip(axs.ravel(), panels):
        vals = [rows[d][k] for d in order]
        e = [rows[d]["platelets_total_counting_error"] for d in order] if err else None
        ax.bar(range(len(order)), vals, color=[col[d] for d in order], yerr=e, capsize=3)
        ax.set_xticks(range(len(order)), [LABEL[d] for d in order], rotation=35, ha="right", fontsize=8)
        ax.set_title(title, fontsize=10)
    fig.suptitle("Denticle-coated stents vs bare stents with the same strut thickness "
                 "(grey = bare controls; error bars = particle counting noise)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "comparison.png"), dpi=110)
    plt.close(fig)

    # ---- per-pair side-by-side deposition-growth GIFs ------------------------
    for tex, ctl in PAIRS + [("bare_thick143", "bare")]:
        if tex in rows and ctl in rows:
            _pair_gif(tex, ctl, dom, os.path.join(OUT, f"pair_{tex}_vs_{ctl}.gif"), 2 * cfg["inlet"]["period_s"])
    print(f"summary written to {OUT}")


def _pair_gif(a, b, dom, path, t_end, n_frames=60):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import trimesh
    from PIL import Image
    import platelets as pl

    data = {}
    for d in (a, b):
        rows = list(csv.DictReader(open(os.path.join(OUT, d, "adhered_platelets.csv"))))
        X = np.array([[float(r["x_mm"]), float(r["y_mm"]), float(r["z_mm"])] for r in rows]).reshape(-1, 3)
        data[d] = dict(x=X[:, 0], s=1.5 * np.arctan2(X[:, 2], X[:, 1]), t=np.array([float(r["time_s"]) for r in rows]),
                       k=np.array([r["site"] for r in rows]), w=float(rows[0]["platelets_represented"]) if rows else 1,
                       outline=pl._stent_outline(os.path.join(HERE, "geometry_out", f"{d}.stl"), dom))
    imgs = []
    tt = np.linspace(t_end / n_frames, t_end, n_frames)
    for t in tt:
        fig = plt.figure(figsize=(14, 8))
        gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.1], hspace=0.3)
        for j, d in enumerate((a, b)):
            ax = fig.add_subplot(gs[0, j])
            q = data[d]
            for o in q["outline"]:
                ax.fill(o[:, 0], o[:, 1], color="0.82", lw=0)
            m = q["t"] <= t
            for site, c, mk in (("stent_luminal", "red", "o"), ("stent_side", "darkorange", "o"),
                                ("aggregate", "darkred", "s")):
                mm = m & (q["k"] == site)
                ax.scatter(q["x"][mm], q["s"][mm], s=7, c=c, marker=mk, lw=0)
            ax.set_aspect("equal")
            ax.set_xlim(dom.x0 - 0.1, dom.x1 + 0.1)
            ax.set_ylim(1.5 * math.radians(dom.theta0_deg), 1.5 * math.radians(dom.theta1_deg))
            ax.set_title(f"{LABEL[d]}\n{int(m.sum() * q['w']):,} platelets deposited", fontsize=11)
            ax.set_xlabel("axial (mm)  → flow")
        ax = fig.add_subplot(gs[1, :])
        for d, c in ((a, "tab:orange"), (b, "0.3")):
            q = data[d]
            ts = np.sort(q["t"])
            ax.plot(ts, np.arange(1, len(ts) + 1) * q["w"], color=c, label=LABEL[d])
        ax.axvline(t, color="crimson", lw=1)
        ax.set_xlim(0, t_end)
        ax.set_xlabel("time (s)")
        ax.set_ylabel("platelets deposited")
        ax.legend(fontsize=9, loc="upper left")
        ax.set_title("red = strut top, orange = strut side/corner, dark squares = aggregated", fontsize=9, loc="left")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=70)
        plt.close(fig)
        buf.seek(0)
        imgs.append(Image.open(buf).convert("RGB"))
    pal = [im.quantize(colors=96, method=Image.Quantize.MEDIANCUT) for im in imgs]
    pal[0].save(path, save_all=True, append_images=pal[1:], duration=120, loop=0, optimize=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=["all", "cfd", "platelets", "summary"], default="all")
    ap.add_argument("--np", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--wait", action="store_true", help="platelet stage: wait for CFD cases still running")
    args = ap.parse_args()
    cfg, dom = setup()
    os.makedirs(OUT, exist_ok=True)
    if args.stage in ("all", "cfd"):
        stage_cfd(args.np)
    if args.stage in ("all", "platelets"):
        stage_platelets(cfg, dom, args.wait)
    if args.stage in ("all", "summary"):
        stage_summary(cfg, dom)


if __name__ == "__main__":
    main()
