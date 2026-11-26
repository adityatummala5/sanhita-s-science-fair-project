"""Infrarenal inflow waveform and the analytical Womersley solution.

Waveform
--------
The resting infrarenal aorta has a triphasic flow waveform: a short systolic
forward peak, a period of reversed flow in early diastole, and little flow in
late diastole. In vitro, near-wall flow there reverses for most of the cycle
(Moore & Ku, J Biomech Eng 1994;116:337). The mean infrarenal flow at rest is
0.9 +/- 0.4 L/min (Taylor et al., Ann Biomed Eng 2002;30:402).

``infrarenal_rest_shape`` is a GENERIC, hand-built triphasic shape. Its mean is
scaled to 0.9 L/min and its peak to the value in the config. Replace it with a
waveform digitised from a published MRI study (e.g. Les et al. 2010) before
final runs: save ``time_s,flow_mL_s`` as CSV and point the config at it.

Womersley solution
------------------
Fully developed oscillatory flow of a Newtonian fluid in a rigid straight tube
(Womersley 1955). For harmonic n with angular frequency w_n and complex flow
amplitude Q_n:

    k  = i^{3/2} sqrt(w_n / nu)
    u_n(r) = Q_n/(pi R^2) * [1 - J0(k r)/J0(k R)] / [1 - 2 J1(kR)/(kR J0(kR))]
    tau_n  = -mu * Q_n/(pi R^2) * [k J1(kR)/J0(kR)] / [1 - 2 J1(kR)/(kR J0(kR))]

with the n = 0 term being Poiseuille flow (u = 2 U (1 - r^2/R^2), tau = 4 mu U / R).
tau is the wall shear stress in the flow direction (positive for forward flow).
"""
from __future__ import annotations

import csv
import math

import numpy as np
from scipy.special import jv


def infrarenal_rest_shape(n: int = 400, period: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Dimensionless triphasic infrarenal shape (systolic peak, early-diastolic
    reversal, small late forward wave)."""
    t = np.linspace(0.0, period, n, endpoint=False)
    p = t / period

    def g(c, w, a):
        d = (p - c + 0.5) % 1.0 - 0.5
        return a * np.exp(-0.5 * (d / w) ** 2)

    f = 0.05 + g(0.14, 0.045, 1.0) - g(0.31, 0.055, 0.42) + g(0.47, 0.06, 0.12)
    return t, f


def rescale(q: np.ndarray, mean: float, peak: float) -> np.ndarray:
    """Affine rescale so that mean(q) = mean and max(q) = peak (reversal allowed)."""
    a = (peak - mean) / (q.max() - q.mean())
    return mean + a * (q - q.mean())


def write_csv(path: str, t: np.ndarray, q_mL_s: np.ndarray) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_s", "flow_mL_s"])
        for a, b in zip(t, q_mL_s):
            w.writerow([f"{a:.5f}", f"{b:.5f}"])


def load_csv(path: str, period: float, mean_mL_s: float | None = None,
             peak_mL_s: float | None = None, samples: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """One cycle, resampled to `samples` equally spaced points in [0, period)."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(r["time_s"]) for r in rows])
    q = np.array([float(r["flow_mL_s"]) for r in rows])
    t = (t - t[0]) * period / (t[-1] - t[0] + (t[1] - t[0]))
    if mean_mL_s is not None and peak_mL_s is not None:
        q = rescale(q, mean_mL_s, peak_mL_s)
    ts = np.linspace(0.0, period, samples, endpoint=False)
    return ts, np.interp(ts, t, q, period=period)


# ---------------------------------------------------------------------------
class Womersley:
    """Analytical Womersley flow in a rigid tube of radius R (m) for a periodic
    flow-rate waveform q(t) (m^3/s) sampled uniformly over one period."""

    def __init__(self, t: np.ndarray, q_m3_s: np.ndarray, period: float, R: float,
                 mu: float = 0.0035, rho: float = 1060.0, n_harm: int = 25):
        self.R, self.mu, self.rho, self.T = R, mu, rho, period
        self.nu = mu / rho
        N = len(q_m3_s)
        c = np.fft.rfft(q_m3_s) / N               # q(t) = Re(sum c_n e^{i w_n t}) with c_n doubled for n>0
        self.n_harm = min(n_harm, len(c) - 1)
        self.c = c[: self.n_harm + 1].copy()
        self.c[1:] *= 2.0
        self.w = 2 * math.pi * np.arange(self.n_harm + 1) / period
        self.alpha = R * np.sqrt(self.w / self.nu)       # Womersley number per harmonic
        self.t0 = t[0]

    @property
    def womersley_number(self) -> float:
        return float(self.alpha[1])

    def _shape(self, n: int, r: np.ndarray):
        """(velocity profile per unit mean velocity, wall shear per unit mean velocity)."""
        R = self.R
        if n == 0:
            return 2.0 * (1.0 - (r / R) ** 2), 4.0 * self.mu / R
        k = (1j ** 1.5) * math.sqrt(self.w[n] / self.nu)
        kR = k * R
        den = 1.0 - 2.0 * jv(1, kR) / (kR * jv(0, kR))
        u = (1.0 - jv(0, k * r) / jv(0, kR)) / den
        tau = -self.mu * (k * jv(1, kR) / jv(0, kR)) / den    # = -mu du/dr at r = R
        return u, tau

    def velocity(self, r: np.ndarray, t: float) -> np.ndarray:
        r = np.asarray(r, float)
        A = math.pi * self.R ** 2
        out = np.zeros_like(r)
        for n in range(self.n_harm + 1):
            u, _ = self._shape(n, r)
            out += np.real(self.c[n] / A * u * np.exp(1j * self.w[n] * (t - self.t0)))
        return out

    def wall_shear(self, t: np.ndarray) -> np.ndarray:
        t = np.atleast_1d(np.asarray(t, float))
        A = math.pi * self.R ** 2
        out = np.zeros_like(t)
        for n in range(self.n_harm + 1):
            _, tau = self._shape(n, np.array([self.R]))
            out += np.real(self.c[n] / A * tau * np.exp(1j * self.w[n] * (t - self.t0)))
        return out

    def flow(self, t: np.ndarray) -> np.ndarray:
        t = np.atleast_1d(np.asarray(t, float))
        return np.real(sum(self.c[n] * np.exp(1j * self.w[n] * (t - self.t0)) for n in range(self.n_harm + 1)))


def reynolds(q_m3_s: float, D: float, mu: float = 0.0035, rho: float = 1060.0) -> float:
    return rho * (q_m3_s / (math.pi * D ** 2 / 4)) * D / mu
