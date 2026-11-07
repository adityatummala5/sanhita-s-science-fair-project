"""Inlet flow waveform for the LAD (left anterior descending) coronary artery.

The pipeline reads a CSV with columns  time_s, flow_mL_min  covering ONE
cardiac cycle. It is resampled, repeated for each cycle and rescaled so its
mean and peak match the targets in the config (Case A in Sanhita's notes:
mean 42.2 mL/min, peak 79.1 mL/min, period 0.9 s = 100 steps x 0.009 s).

IMPORTANT: config/waveform_LAD_approx.csv is a generic, hand-built coronary
shape (low flow in systole, peak in early diastole). Replace it with the
waveform digitised from the source paper (e.g. using WebPlotDigitizer) before
reporting final results.
"""
from __future__ import annotations

import csv
import numpy as np


def generic_lad_shape(n: int = 200, period: float = 0.9) -> tuple[np.ndarray, np.ndarray]:
    """Dimensionless coronary-like flow shape: small systolic bump, systolic
    dip, dominant early-diastolic peak, then slow diastolic decay."""
    t = np.linspace(0.0, period, n, endpoint=False)
    p = t / period

    def g(c, w, a):  # periodic Gaussian bump
        d = (p - c + 0.5) % 1.0 - 0.5
        return a * np.exp(-0.5 * (d / w) ** 2)

    f = 0.30 + g(0.10, 0.05, 0.25) - g(0.25, 0.05, 0.12) + g(0.47, 0.09, 0.95) + g(0.72, 0.16, 0.35)
    return t, f


def write_default_csv(path: str, mean: float = 42.2, peak: float = 79.1, period: float = 0.9):
    t, f = generic_lad_shape(period=period)
    q = rescale(f, mean, peak)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_s", "flow_mL_min"])
        for a, b in zip(t, q):
            w.writerow([f"{a:.5f}", f"{b:.4f}"])


def rescale(q: np.ndarray, mean: float, peak: float) -> np.ndarray:
    """Affine rescale so that mean(q) = mean and max(q) = peak."""
    a = (peak - mean) / (q.max() - q.mean())
    out = mean + a * (q - q.mean())
    if out.min() <= 0:
        raise ValueError("rescaled waveform has non-positive flow; check targets")
    return out


def load(path: str, period: float, mean: float | None, peak: float | None,
         n_cycles: int, samples_per_cycle: int = 100) -> tuple[np.ndarray, np.ndarray]:
    """Return (times, flow mL/min) sampled over n_cycles (+1 sample past the end)."""
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh)]
    t = np.array([float(r["time_s"]) for r in rows])
    q = np.array([float(r["flow_mL_min"]) for r in rows])
    t = (t - t[0]) * period / (t[-1] - t[0] + (t[1] - t[0]))   # stretch to the period
    if mean is not None and peak is not None:
        q = rescale(q, mean, peak)
    ts = np.linspace(0.0, period, samples_per_cycle, endpoint=False)
    qs = np.interp(ts, t, q, period=period)
    T = np.concatenate([ts + k * period for k in range(n_cycles)] + [[n_cycles * period]])
    Q = np.concatenate([qs] * n_cycles + [[qs[0]]])
    return T, Q
