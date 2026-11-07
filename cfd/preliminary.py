#!/usr/bin/env python3
"""Preliminary check on the bare stent: three simulations, each exported as GIF + data.

  1_newtonian   pulsatile flow, Newtonian fluid with blood's density and high-shear
                viscosity (rho 1060 kg/m3, mu 3.5 mPa s)
  2_carreau     pulsatile flow, shear-thinning Carreau blood (Cho & Kensey 1991)
  3_platelets   platelet transport, adhesion and aggregation (Lagrangian, Sorensen 1999
                kinetics) driven by the Carreau flow of simulation 2

Usage (from the cfd/ folder):
  python3 preliminary.py                 # runs the CFD if needed, then exports everything
  python3 preliminary.py --skip-cfd      # only re-export / re-run platelets from existing cases
  python3 preliminary.py --platelet-cycles 3
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "pipeline"))

import case as fc        # noqa: E402

PRESET, CYCLES, SNAPS = "coarse", 2, 20
CASES = {"1_newtonian": "bare_newtonian", "2_carreau": "bare"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-cfd", action="store_true")
    ap.add_argument("--skip-platelets", action="store_true")
    ap.add_argument("--np", type=int, default=os.cpu_count() or 2)
    ap.add_argument("--platelet-cycles", type=float, default=2.0)
    args = ap.parse_args()

    cfg = json.load(open(os.path.join(HERE, "config/experiment.json")))
    w = cfg["window"]
    dom = fc.Domain(x0=w["x0"], x1=w["x1"], theta0_deg=w["theta0_deg"], theta1_deg=w["theta1_deg"],
                    **cfg.get("domain", {}))
    period = cfg["inlet"]["period_s"]
    rho = cfg["physics"]["rho"]
    run_root = os.path.join(HERE, "runs", f"{PRESET}-pulsatile")
    out_root = os.path.join(HERE, "results", "preliminary")
    os.makedirs(out_root, exist_ok=True)

    if not args.skip_cfd:
        need = [d for d in CASES.values() if not os.path.exists(os.path.join(run_root, d, "metrics.json"))]
        if need:
            subprocess.run([sys.executable, os.path.join(HERE, "run.py"), "--preset", PRESET, "--mode", "pulsatile",
                            "--only", ",".join(need), "--np", str(args.np), "--cycles", str(CYCLES),
                            "--snapshots", str(SNAPS), "--keep-snapshots"], check=True)

    import animate  # noqa: E402
    titles = {"1_newtonian": "1) Newtonian (μ = 3.5 mPa·s)", "2_carreau": "2) Carreau blood (shear-thinning)"}
    metrics = {}
    for key, design in CASES.items():
        cdir = os.path.join(run_root, design)
        print(f"exporting {key} from {cdir}")
        animate.export_flow(cdir, os.path.join(out_root, key), dom, rho, period, titles[key])
        metrics[key] = json.load(open(os.path.join(cdir, "metrics.json")))
        json.dump(metrics[key], open(os.path.join(out_root, key, "metrics.json"), "w"), indent=2)
    _compare(out_root, metrics)

    if not args.skip_platelets:
        import platelets as pl  # noqa: E402
        cdir = os.path.join(run_root, CASES["2_carreau"])
        times = animate.snapshot_times(cdir, period)
        p = pl.PlateletParams(n_cycles=args.platelet_cycles)
        print("platelet simulation on the Carreau flow")
        s = pl.simulate(cdir, times, period, dom, os.path.join(HERE, "geometry_out", "bare.stl"),
                        os.path.join(out_root, "3_platelets"), p, title="3) bare stent, Carreau blood + platelets")
        print(json.dumps({k: v for k, v in s.items() if k != "parameters"}, indent=2, default=float))
    print(f"\nResults in {out_root}")


def _compare(out_root, metrics):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    keys = [k for k in metrics["1_newtonian"] if isinstance(metrics["1_newtonian"][k], (int, float))
            and not isinstance(metrics["1_newtonian"][k], bool)]
    with open(os.path.join(out_root, "comparison_newtonian_vs_carreau.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["metric", "newtonian", "carreau", "carreau_vs_newtonian_pct"])
        for k in keys:
            a, b = metrics["1_newtonian"][k], metrics["2_carreau"].get(k)
            pct = 100 * (b - a) / a if (b is not None and a) else ""
            wr.writerow([k, a, b, f"{pct:.2f}" if pct != "" else ""])

    series = {}
    for key in ("1_newtonian", "2_carreau"):
        rows = list(csv.DictReader(open(os.path.join(out_root, key, "timeseries.csv"))))
        series[key] = {c: np.array([float(r[c]) for r in rows]) for c in rows[0]}
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.5))
    plots = [("inlet_flow_mL_min_full_vessel", "inflow (mL/min)"), ("wall_WSS_mean_Pa", "artery-wall WSS under stent (Pa)"),
             ("wall_lowWSS_area_pct", "wall area with WSS < 0.4 Pa (%)"), ("pressure_drop_Pa", "pressure drop (Pa)")]
    for ax, (col, lab) in zip(axs.ravel(), plots):
        for key, st in (("1_newtonian", "b-o"), ("2_carreau", "r-s")):
            s = series[key]
            o = np.argsort(s["phase"])
            ax.plot(s["phase"][o] * 0.9, s[col][o], st, ms=3, label=key)
        ax.set_xlabel("time in cardiac cycle (s)")
        ax.set_ylabel(lab)
        ax.legend(fontsize=8)
    fig.suptitle("Bare stent: Newtonian vs Carreau blood over one cardiac cycle")
    plt.tight_layout()
    plt.savefig(os.path.join(out_root, "comparison_newtonian_vs_carreau.png"), dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
