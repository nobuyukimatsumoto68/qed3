#!/usr/bin/env python3
# diag_fit_claude.py
# Per-diagram const+exp fit  C(dt) = c0 + A exp(-m dt)  of the 10 two-meson diagrams (PS legs),
# translation-averaged, config jackknife.  Reports the mass m per diagram.  CONTACT env (0 / 0.5).
# Fit window [FIT_LO, FIT_HI] (env; default [2,16]).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
from scipy.optimize import curve_fit
import distill_contract_claude as dc
import diag_effmass_claude as de

CONTACT = float(os.environ.get("CONTACT", "0.5"))
FIT_LO = int(os.environ.get("FIT_LO", "2"))
FIT_HI = int(os.environ.get("FIT_HI", "16"))
LABELS = de.LABELS


def cexp(dt, c0, A, m):
    return c0 + A * np.exp(-m * dt)


def fit_one(dts, C):
    A0 = C[0] - C[-1]
    c0 = C[-1]
    m0 = 0.5
    if abs(A0) < 1e-300:
        A0 = C[0] if C[0] != 0 else 1e-12
    try:
        p, _ = curve_fit(cexp, dts, C, p0=[c0, A0, m0], maxfev=20000)
        return p[2]
    except Exception:
        return np.nan


def main():
    de.CONTACT = CONTACT
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  ncfg=%d  CONTACT=%.2f  const+exp per-diagram fit  window dt[%d,%d]"
          % (tag, len(dc.KS), CONTACT, FIT_LO, FIT_HI))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allD = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        D = np.zeros((10, twin))
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += de.diags_pair(Phi, tau, s, s + dt)
            D[:, dt] = (acc / ns).real
        allD.append(D)
    allD = np.array(allD)                              # (ncfg, 10, twin)
    ncfg = allD.shape[0]
    dts = np.arange(FIT_LO, FIT_HI + 1)

    print("\n  diagram         m (const+exp)     |signal|")
    for i in range(10):
        Ci = allD[:, i, :]                             # (ncfg, twin)
        samp = np.array([np.delete(Ci, j, 0).mean(0) for j in range(ncfg)])
        ms = np.array([fit_one(dts, samp[j, dts]) for j in range(ncfg)])
        good = np.isfinite(ms)
        if good.sum() < 2:
            print("  %-14s  fit failed" % LABELS[i])
            continue
        m = ms[good].mean()
        err = np.sqrt((good.sum() - 1) * np.mean((ms[good] - m) ** 2))
        cm = Ci.mean(0)
        sig = abs(cm[FIT_LO] - cm[FIT_HI])
        print("  %-14s  %7.4f(%.4f)   %.3e%s"
              % (LABELS[i], m, err, sig, "" if good.all() else "   [%d/%d fits ok]" % (good.sum(), ncfg)))


if __name__ == "__main__":
    main()
