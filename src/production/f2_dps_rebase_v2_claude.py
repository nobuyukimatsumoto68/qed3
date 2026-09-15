#!/usr/bin/env python3
# f2_dps_rebase_v2_claude.py
#   {D'_S, F^2} 2-op GEVP, TRIPLE-SUBTRACTED (per-config t-sum -> per-element plateau; vacuum-free, no identity),
#   with T0=0 (clean equal-time metric) + REBASE to NKEEP=2 states at REBT=1 (hankel_rebase machinery,
#   offsets=[0] = pure rebase).  D'_S=Tr[Phi tilde_tau Phi tilde_tau] (local connected sigma^2), F^2 shape op.
#   Window-triangular sampling (CS-consistent); operators normalized at TREF.  Effmass = LATTICE a_t m.
#   Run: ENS=... NVDIR=distill_Nv24 BINSIZE=10 T0=0 REBT=1 NKEEP=2 python3 f2_dps_rebase_v2_claude.py

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
import hankel_rebase_scan_claude as hs

Nt = 128
OP_F = int(os.environ.get("OP_F", "0"))
AT = float(os.environ.get("AT", "0.2"))
CONTACT = float(os.environ.get("CONTACT", "0.5"))
TMAX = int(os.environ.get("TMAX", "24"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
T0 = int(os.environ.get("T0", "0"))
REBT = int(os.environ.get("REBT", "1"))
NKEEP = int(os.environ.get("NKEEP", "2"))
TREF = int(os.environ.get("TREF", "1"))


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


def rebase_em(C2, Vfix):
    Cts = np.transpose(C2, (2, 0, 1))
    Big = hs.hankel_off(Cts, [0])
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k))]
    ncfg = len(ks)
    print("# ENS=%s cfg=%d  {D'_S,F^2} triple-sub  T0=%d rebase %d@%d  BIN=%d" % (tag, ncfg, T0, NKEEP, REBT, BINSIZE))

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

    # TRIPLE SUB step 1: per-config t-sum
    def tsum(a):
        return a - a[:, 1:TMAX].mean(1, keepdims=True)
    cff, cfd, cdd = tsum(cff), tsum(cfd), tsum(cdd)

    nb = ncfg // BINSIZE
    def binned(a):
        return np.array([a[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    bff, bfd, bdd = binned(cff), binned(cfd), binned(cdd)
    PLAT = TMAX - 8

    def build_C(sel):
        ff = bff[sel].mean(0); ff = ff - np.nanmean(ff[PLAT:TMAX])   # step 2: per-element plateau
        fd = bfd[sel].mean(0); fd = fd - np.nanmean(fd[PLAT:TMAX])
        dd = bdd[sel].mean(0); dd = dd - np.nanmean(dd[PLAT:TMAX])
        C = np.zeros((2, 2, TMAX + 1))
        C[0, 0], C[1, 1], C[0, 1], C[1, 0] = dd, ff, fd, fd          # op0=D'_S, op1=F^2
        return C

    Craw = build_C(np.arange(nb))
    dn = np.sqrt(np.abs([Craw[0, 0, TREF], Craw[1, 1, TREF]]))
    NRM = np.outer(dn, dn)[:, :, None]

    def C_norm(sel):
        return build_C(sel) / NRM

    # rho diagnostic
    print("\n# DIAG t | rho=C_Fd/sqrt(C_dd C_FF) (triple-sub, normalized)")
    Cn = C_norm(np.arange(nb))
    for t in range(0, min(TMAX, 10)):
        dd, ff, fd = Cn[0, 0, t], Cn[1, 1, t], Cn[0, 1, t]
        rho = fd / np.sqrt(dd * ff) if (dd > 0 and ff > 0) else np.nan
        print("#  %2d | rho=%7.3f  (C_dd=%.3e C_FF=%.3e)" % (t, rho, dd, ff))

    em_c, Vfix = rebase_em(C_norm(np.arange(nb)), None)
    ems = np.array([rebase_em(C_norm(np.delete(np.arange(nb), j)), Vfix)[0] for j in range(nb)])
    em_e = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, 0))
    tmax = em_c.shape[0]

    print("\n#  t |  " + "  ".join("state%d(err)" % n for n in range(NKEEP)) + "   [sig2 gd 0.46, F2 0.576, 2mPS 0.644]")
    for t in range(0, min(tmax, TMAX - 1)):
        if np.any(np.isfinite(em_c[t])):
            print("#  %2d | " % t + "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(NKEEP)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9, 5.8))
    ax.axhline(0.644, color="gray", ls="--", lw=1, alpha=0.6); ax.text(15, 0.66, r"$2m_{PS}=0.644$", fontsize=8, color="gray")
    ax.axhline(0.576, color="black", ls="-.", lw=1, alpha=0.6); ax.text(15, 0.585, r"$F^2=0.576$", fontsize=8, color="black")
    ax.axhline(0.46, color="tab:green", ls=":", lw=1, alpha=0.6); ax.text(15, 0.475, r"$\sigma^2$ gd 0.46", fontsize=8, color="tab:green")
    cols = ["tab:blue", "tab:red"]
    mks = ["o", "s"]
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.4)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n], marker=mks[n], ms=5, lw=1.1, capsize=2.5,
                    label="state %d" % n)
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.1, 0.95)
    ax.set_xlim(-0.3, TMAX - 4)
    ax.set_xlabel(r"$t$ (timeslice)")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$ (lattice)")
    ax.set_title(r"$\{D'_S,F^2\}$ triple-sub  T0=%d rebase %d@%d  %s L1 %d cfg (bin %d)" % (T0, NKEEP, REBT, tag, ncfg, BINSIZE), fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_dps_rebase_v2_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
