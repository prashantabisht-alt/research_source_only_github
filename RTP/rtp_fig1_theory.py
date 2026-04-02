# rtp_fig1_theory.py
# Exact P(x,t) for 1D RTP on the infinite line (dimensionless: gamma=1, v=1)
# Implements Eq. (10) of Malakar et al. (2018) via trapezoid quadrature.
# Writes: fig1_theory_D0.03.txt, fig1_theory_D0.10.txt, fig1_theory_D0.20.txt

import math
from typing import Iterable

import numpy as np

NBINS = 300
X_MIN = -2.2
X_MAX = 2.2
TGRID = np.round(np.arange(0.2, 2.0 + 1e-12, 0.2), 1)
DLIST = [0.03, 0.10, 0.20]

NY_BASE = 4000
I_TOL = 1e-14
MAX_SER = 2000


def header_line(times: Iterable[float]) -> str:
    return "# x" + "".join([f"   t{t:.1f}" for t in times])


def make_fname(D: float) -> str:
    tag = int(round(100 * D))
    return f"fig1_theory_D0.{tag:02d}.txt"


def I0_series(z: np.ndarray, tol: float = I_TOL, max_terms: int = MAX_SER) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    out = np.ones_like(z)
    term = np.ones_like(z)
    m = 0
    zz = (z * 0.5) ** 2
    while True:
        m += 1
        term *= zz / (m * m)
        out += term
        if (np.abs(term) < tol).all() or m >= max_terms:
            break
    return out


def I1_series(z: np.ndarray, tol: float = I_TOL, max_terms: int = MAX_SER) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    out = 0.5 * z
    term = 0.5 * z
    m = 0
    zz = (z * 0.5) ** 2
    while True:
        m += 1
        term *= zz / (m * (m + 1))
        out += term
        if (np.abs(term) < tol).all() or m >= max_terms:
            break
    return out


def P_theory_xt(xcent: np.ndarray, t: float, D: float, ny: int) -> np.ndarray:
    xcent = xcent.reshape(-1, 1)
    nx = xcent.shape[0]

    term1 = np.cosh(xcent / (2 * D)) / math.sqrt(4 * math.pi * D * t) * np.exp(
        -t - (xcent**2 + t * t) / (4 * D * t)
    )

    if t == 0.0:
        integ = np.zeros((nx, 1))
    else:
        ny = int(ny)
        if ny % 2 == 1:
            ny += 1
        y = np.linspace(-t, t, ny + 1)
        dy = y[1] - y[0]
        z = np.sqrt(np.maximum(0.0, t * t - y * y))
        pref = np.exp(-(xcent - y) ** 2 / (4 * D * t)) / math.sqrt(4 * math.pi * D * t)

        i0z = I0_series(z)
        i1z = I1_series(z)

        add = i0z.copy()
        nz = z > 0
        add[nz] = add[nz] + (t / z[nz]) * i1z[nz]
        add[~nz] = add[~nz] + 0.5 * t

        w = np.ones_like(y)
        w[0] = 0.5
        w[-1] = 0.5
        integ = (pref * add) @ (w * dy)

    p = term1 + 0.5 * np.exp(-t) * integ.reshape(-1, 1)
    return p.ravel()


def main() -> None:
    binw = (X_MAX - X_MIN) / NBINS
    xcent = X_MIN + (np.arange(NBINS) + 0.5) * binw

    for D in DLIST:
        fname = make_fname(D)
        print(f"[theory] D={D:.2f} -> {fname}")
        with open(fname, "w", encoding="utf-8") as f:
            f.write(header_line(TGRID) + "\n")
            for x in xcent:
                row = [f"{x: .12E}"]
                for t in TGRID:
                    ny = max(NY_BASE, int(2000 * t))
                    p = P_theory_xt(np.array([x]), float(t), float(D), ny)[0]
                    row.append(f"{p: .12E}")
                f.write(" ".join(row) + "\n")


if __name__ == "__main__":
    main()
