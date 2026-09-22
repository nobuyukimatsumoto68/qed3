#!/usr/bin/env python3
# f2_dps_id_rebase_v2_claude.py
#   {1, D'_S, F^2} GEVP with T0=0 (clean equal-time metric) + REBASE to NKEEP=2 states at REBT=1.
#   The T0=0 metric is positive-definite (all correlators have signal at t=0, unlike large-t where C_FF
#   is noise), and the rebase projects the 3-op basis onto the 2 states dominating at t=1.  Reuses the
#   hankel_rebase machinery (staged_project / rebased_effmass_fixed) with offsets=[0] (no Hankel, pure rebase).
#   Identity absorbs the vacuum (raw correlators).  Effmass = LATTICE a_t m.
#   Run: ENS=... NVDIR=distill_Nv24 BINSIZE=10 T0=0 REBT=1 NKEEP=2 python3 f2_dps_id_rebase_v2_claude.py

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


def rebase_em(C3, Vfix):
    Cts = np.transpose(C3, (2, 0, 1))                  # (T,3,3)
    Big = hs.hankel_off(Cts, [0])                      # pure rebase (no Hankel offsets)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k))]
    ncfg = len(ks)
    print("# ENS=%s cfg=%d  {1,D'_S,F^2} T0=%d rebase %d@%d  BIN=%d" % (tag, ncfg, T0, NKEEP, REBT, BINSIZE))

    cff = np.zeros((ncfg, TMAX + 1))
    cfd = np.zeros((ncfg, TMAX + 1))
    cdd = np.zeros((ncfg, TMAX + 1))
    ofm = np.zeros(ncfg)
    dm = np.zeros(ncfg)
    for ic, k in enumerate(ks):
        of = load_of(k)
        D, tsrc0, twin = dps_series(k, w00)
        ofw = of[(tsrc0 + np.arange(twin)) % Nt]
        ofm[ic] = ofw.mean()
        dm[ic] = D.mean()
        for t in range(TMAX + 1):
            if t < twin:
                a = np.arange(twin - t)
                cff[ic, t] = np.mean(ofw[a + t] * ofw[a])
                cdd[ic, t] = np.mean(D[a + t] * D[a])
                cfd[ic, t] = 0.5 * (np.mean(ofw[a + t] * D[a]) + np.mean(D[a + t] * ofw[a]))
            else:
                cff[ic, t] = cdd[ic, t] = cfd[ic, t] = np.nan

    nb = ncfg // BINSIZE
    def binned(a):
        return np.array([a[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    bff, bfd, bdd, bof, bdm = binned(cff), binned(cfd), binned(cdd), binned(ofm), binned(dm)

    def build_C(sel):
        ff, fd, dd = bff[sel].mean(0), bfd[sel].mean(0), bdd[sel].mean(0)
        of1, dm1 = bof[sel].mean(), bdm[sel].mean()
        C = np.zeros((3, 3, TMAX + 1))
        C[0, 0] = 1.0
        C[0, 1] = C[1, 0] = dm1
        C[0, 2] = C[2, 0] = of1
        C[1, 1], C[2, 2], C[1, 2], C[2, 1] = dd, ff, fd, fd
        return C

    Craw = build_C(np.arange(nb))
    dn = np.sqrt(np.abs([Craw[0, 0, TREF], Craw[1, 1, TREF], Craw[2, 2, TREF]]))
    NRM = np.outer(dn, dn)[:, :, None]

    def C_norm(sel):
        return build_C(sel) / NRM

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
    cols = ["tab:blue", "tab:red", "tab:purple"]
    mks = ["o", "s", "^"]
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.4)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n], marker=mks[n], ms=5, lw=1.1, capsize=2.5,
                    label="state %d" % n)
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.1, 0.95)
    ax.set_xlim(-0.3, TMAX - 4)
    ax.set_xlabel(r"$t$ (timeslice)")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$ (lattice)")
    ax.set_title(r"$\{1,D'_S,F^2\}$ T0=%d rebase %d@%d  %s L1 %d cfg (bin %d)" % (T0, NKEEP, REBT, tag, ncfg, BINSIZE), fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_dps_id_rebase_v2_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
