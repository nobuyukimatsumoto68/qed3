#!/usr/bin/env python3
# f2_dps_gevp_v2_claude.py
#   Coupled 2x2 {D'_S, F^2} GEVP.  D'_S(s)=Tr[Phi tilde_tau Phi tilde_tau] = local connected sigma^2 density
#   (the "ext" 0++ fermionic operator); F^2 = glueball shape op.  TRIPLE SUBTRACTION for the vacuum (per-config
#   t-sum -> per-element plateau), as always.  All three correlators on the SAME window-triangular sampling
#   (s, s+t in [0,twin); cross symmetrized) so Cauchy-Schwarz holds.  Operators normalized at TREF to fix the
#   ~13-orders scale disparity (glue vs fermion loop; GEVP-invariant).  Fixed-t0 GEVP, effmass = log(lam(t)/
#   lam(t+1)) = LATTICE a_t m (overlay units; physical = /a_t).
#   Run: ENS=... NVDIR=distill_Nv24 BINSIZE=10 T0=3 python3 f2_dps_gevp_v2_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import h5py
import distill_contract_claude as dc

Nt = 128
OP_F = int(os.environ.get("OP_F", "0"))
AT = float(os.environ.get("AT", "0.2"))
CONTACT = float(os.environ.get("CONTACT", "0.5"))
TMAX = int(os.environ.get("TMAX", "24"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
T0 = int(os.environ.get("T0", "3"))
TREF = int(os.environ.get("TREF", "2"))
RTOL = float(os.environ.get("RTOL", "1e-9"))


def gluedir():
    return "data_" + dc.ENS


def load_of(k):
    with h5py.File("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k), "r") as g:
        return np.array(g["O"])[OP_F].astype(float)


def dps_series(k, w00):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Iv = np.eye(tau.shape[-1])
    D = np.empty(twin)
    for a in range(twin):
        Phi = (V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T))
        tt = tau[a, a] - CONTACT * Iv
        D[a] = np.trace(Phi @ tt @ Phi @ tt).real
    return D, tsrc0, twin


def inv_sqrt_sym(M, rtol):
    M = 0.5 * (M + M.T)
    w, Vv = np.linalg.eigh(M)
    inv = np.where((w > rtol * w.max()) & (w > 0), 1.0 / np.sqrt(np.abs(w)), 0.0)
    return (Vv * inv) @ Vv.T


def gevp_lat(C, t, t0):
    # fixed-t0 GEVP eigenvalues lambda(t): C(t) v = lambda C(t0) v
    si = inv_sqrt_sym(0.5 * (C[:, :, t0] + C[:, :, t0].T), RTOL)
    M = si @ (0.5 * (C[:, :, t] + C[:, :, t].T)) @ si
    return np.sort(np.linalg.eigvalsh(0.5 * (M + M.T)))[::-1]     # descending -> [0]=lightest


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k))]
    ncfg = len(ks)
    print("# ENS=%s cfg=%d  {D'_S, F^2} GEVP (triple-sub, lattice a_t m)  T0=%d BIN=%d" % (tag, ncfg, T0, BINSIZE))

    # per-config raw correlators (window-triangular, cross symmetrized)
    cff = np.zeros((ncfg, TMAX + 1))
    cfd = np.zeros((ncfg, TMAX + 1))
    cdd = np.zeros((ncfg, TMAX + 1))
    for ic, k in enumerate(ks):
        of = load_of(k)
        D, tsrc0, twin = dps_series(k, w00)
        ofw = of[(tsrc0 + np.arange(twin)) % Nt]
        for t in range(TMAX + 1):
            if t < twin:
                a = np.arange(twin - t)
                cff[ic, t] = np.mean(ofw[a + t] * ofw[a])
                cdd[ic, t] = np.mean(D[a + t] * D[a])
                cfd[ic, t] = 0.5 * (np.mean(ofw[a + t] * D[a]) + np.mean(D[a + t] * ofw[a]))
            else:
                cff[ic, t] = cdd[ic, t] = cfd[ic, t] = np.nan

    # (1) TRIPLE SUB, step 1: per-config t-sum (subtract each config's dt-mean over [1,TMAX))
    def tsum(a):
        return a - a[:, 1:TMAX].mean(1, keepdims=True)
    cff, cfd, cdd = tsum(cff), tsum(cfd), tsum(cdd)

    nb = ncfg // BINSIZE
    def binned(a):
        return np.array([a[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    bff, bfd, bdd = binned(cff), binned(cfd), binned(cdd)

    PLAT = TMAX - 8
    def build_C(sel):
        # step 2: per-element plateau subtraction ; assemble normalized 2x2
        ff = bff[sel].mean(0); ff = ff - np.nanmean(ff[PLAT:TMAX])
        fd = bfd[sel].mean(0); fd = fd - np.nanmean(fd[PLAT:TMAX])
        dd = bdd[sel].mean(0); dd = dd - np.nanmean(dd[PLAT:TMAX])
        C = np.zeros((2, 2, TMAX + 1))
        C[0, 0], C[0, 1], C[1, 0], C[1, 1] = dd, fd, fd, ff     # op0=D'_S, op1=F^2
        return C

    Craw = build_C(np.arange(nb))
    dn = np.sqrt(np.abs([Craw[0, 0, TREF], Craw[1, 1, TREF]]))
    NRM = np.outer(dn, dn)[:, :, None]

    def C_norm(sel):
        return build_C(sel) / NRM

    # diagnostic: rho(t)
    print("\n# DIAG t | C_dd        C_FF         C_Fd        rho")
    for t in range(0, min(TMAX, 12)):
        dd, ff, fd = Craw[0, 0, t], Craw[1, 1, t], Craw[0, 1, t]
        rho = fd / np.sqrt(dd * ff) if (dd > 0 and ff > 0) else np.nan
        print("#  %2d | %11.3e %11.3e %11.3e  %7.3f" % (t, dd, ff, fd, rho))

    def effmass(sel):
        lam = np.full((TMAX + 1, 2), np.nan)
        C = C_norm(sel)
        for t in range(1, TMAX + 1):
            if np.all(np.isfinite(C[:, :, t])) and np.all(np.isfinite(C[:, :, T0])):
                try:
                    lam[t] = gevp_lat(C, t, T0)
                except Exception:
                    pass
        with np.errstate(all="ignore"):
            r = lam[:-1] / lam[1:]
            r[r <= 0] = np.nan
            return np.log(r)

    em_c = effmass(np.arange(nb))
    ems = np.array([effmass(np.delete(np.arange(nb), j)) for j in range(nb)])
    em_e = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, 0))

    print("\n#  t | state0 a_t m(err)   state1 a_t m(err)    [sigma2 gd 0.46, F2 0.576, 2mPS 0.644]")
    for t in range(1, TMAX - 1):
        if np.isfinite(em_c[t, 0]):
            print("#  %2d | %7.4f(%.4f)     %7.4f(%.4f)" % (t, em_c[t, 0], em_e[t, 0], em_c[t, 1], em_e[t, 1]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(TMAX)
    fig, ax = plt.subplots(figsize=(9, 5.8))
    ax.axhline(0.644, color="gray", ls="--", lw=1, alpha=0.6); ax.text(15, 0.65, r"$2m_{PS}=0.644$", fontsize=8, color="gray")
    ax.axhline(0.576, color="black", ls="-.", lw=1, alpha=0.6); ax.text(15, 0.58, r"$F^2=0.576$", fontsize=8, color="black")
    ax.axhline(0.46, color="tab:green", ls=":", lw=1, alpha=0.6); ax.text(15, 0.47, r"$\sigma^2$ gd 0.46", fontsize=8, color="tab:green")
    for n, (c, m, lab) in enumerate([("tab:blue", "o", "state 0"), ("tab:red", "s", "state 1")]):
        g = np.isfinite(em_c[:TMAX, n]) & np.isfinite(em_e[:TMAX, n]) & (em_e[:TMAX, n] < 0.4)
        ax.errorbar(ts[g], em_c[:TMAX][g, n], yerr=em_e[:TMAX][g, n], color=c, marker=m, ms=5, lw=1.1, capsize=2.5, label=lab)
    ax.set_ylim(0.2, 0.95)
    ax.set_xlim(0.5, TMAX - 2)
    ax.set_xlabel(r"$t$ (timeslice)")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$ (lattice)")
    ax.set_title(r"$\{D'_S, F^2\}$ GEVP (triple-sub)  T0=%d  %s L1 %d cfg (bin %d)" % (T0, tag, ncfg, BINSIZE), fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_dps_gevp_v2_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
