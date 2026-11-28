"""Infrarenal aortic inflow for the EVAR model.

waveform: a GENERIC triphasic resting infrarenal shape (systolic peak, early
diastolic reverse flow, small late-diastolic forward wave), rescaled to
    mean 17.5 mL/s   (Les et al. 2010, Cardiovasc Eng Technol 1:39, 43 AAA
                      patients, phase-contrast MRI, cycle-averaged IR flow)
    peak 90 mL/s     (assumed: ~5x mean, typical of resting infrarenal flow)
The shape is NOT digitised from a paper. Replace config/waveform_infrarenal_generic.csv
with a digitised waveform (e.g. Les et al. 2010 Fig. 3) before final results.

inlet profile: Womersley solution for fully developed pulsatile flow in a
straight rigid pipe of radius R (Newtonian, nu = mu_inf / rho), harmonic by
harmonic. With Womersley numbers of ~19 a parabolic profile would be badly
wrong (it cannot carry the near-wall reverse flow of a real aortic pulse).
"""
from __future__ import annotations

import csv
import math

import numpy as np
from scipy.special import jv


def generic_infrarenal_shape(n: int = 400, period: float = 0.9):
    t = np.linspace(0.0, period, n, endpoint=False)
    p = t / period

    def g(c, w, a):
        d = (p - c + 0.5) % 1.0 - 0.5
        return a * np.exp(-0.5 * (d / w) ** 2)

    f = g(0.14, 0.045, 1.0) - g(0.33, 0.05, 0.45) + g(0.50, 0.07, 0.10)
    return t, f


def rescale_signed(f: np.ndarray, mean: float, peak: float) -> np.ndarray:
    """Affine rescale to the target mean and peak; negative values allowed."""
    a = (peak - mean) / (f.max() - f.mean())
    return mean + a * (f - f.mean())


def write_default_csv(path: str, mean_mL_s: float = 17.5, peak_mL_s: float = 90.0, period: float = 0.9):
    t, f = generic_infrarenal_shape(period=period)
    q = rescale_signed(f, mean_mL_s, peak_mL_s)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_s", "flow_mL_s"])
        for a, b in zip(t, q):
            w.writerow([f"{a:.5f}", f"{b:.4f}"])


def load_cycle(path: str, period: float, n: int = 400):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(r["time_s"]) for r in rows])
    q = np.array([float(r["flow_mL_s"]) for r in rows])
    t = (t - t[0]) * period / (t[-1] - t[0] + (t[1] - t[0]))
    ts = np.linspace(0.0, period, n, endpoint=False)
    return ts, np.interp(ts, t, q, period=period)


def womersley_profile(q_cycle: np.ndarray, period: float, R: float, nu: float,
                      r: np.ndarray, n_harm: int = 25) -> np.ndarray:
    """Axial velocity u(t_k, r_j) (m/s) for one cycle sampled at len(q_cycle)
    equally spaced times. q_cycle in m^3/s, R and r in m."""
    N = len(q_cycle)
    Q = np.fft.rfft(q_cycle) / N                     # Q(t) = sum Re(c_n e^{i n w t}) with c_0 = Q0, c_n = 2 Q_n
    w = 2 * np.pi / period
    t = np.arange(N) * period / N
    A = math.pi * R ** 2
    xi = np.clip(r / R, 0.0, 1.0)
    u = np.outer(np.full(N, 2 * Q[0].real / A), 1 - xi ** 2)   # mean part: Poiseuille
    for n in range(1, min(n_harm, len(Q) - 1) + 1):
        alpha = R * math.sqrt(n * w / nu)
        z = alpha * 1j ** 1.5
        J0z, J1z = jv(0, z), jv(1, z)
        shape = (1 - jv(0, z * xi) / J0z) / (1 - 2 * J1z / (z * J0z))   # mean over the area = 1
        cn = 2 * Q[n]
        u += np.real(np.outer(cn * np.exp(1j * n * w * t), shape)) / A
    return u


def womersley_number(R: float, period: float, nu: float) -> float:
    return R * math.sqrt(2 * math.pi / period / nu)
