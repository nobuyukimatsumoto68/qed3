#!/usr/bin/env python3
# sigma2_meson_muweight_gevp_claude.py -- route-1 cross-check for the (2,2): a pure single-meson GEVP over the
#   family O_p = psibar (D~^H D~)^p psi, whose vertex in the distillation basis is DIAGONAL, Phi_p = diag(mu_a^{2p})
#   = diag(evals^p) (evals = mu_a^2 = D~^H D~ eigenvalues, stored per timeslice in /evals).  Different powers p
#   emphasize different mode/energy bands, so the GEVP resolves the single-meson tower: 2E_0 = m_PS (ground),
#   2E_1 = (2,2) (first excited, ell=3/2), ...  NO shell projector, NO m_PS clearing -- the GEVP orthogonalizes.
#   Correlator: C_{pq}(dt) = avg_s -Tr[ Phi_p(t) tau(t,s) Phi_q(s) tau(s,t) ]  (mode space, Nv x Nv, from /peram/tau).
#   Run: ENS=<L2> LREF=2 NVDIR=distill_Nv24_v2 POWERS=0,-1,-2 REBT=4 NKEEP=3 T0=2 BINSIZE=10 python3 sigma2_meson_muweight_gevp_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import h5py
import distill_contract_claude as dc
import hankel_rebase_scan_claude as hs

POWERS = [float(x) for x in os.environ.get("POWERS", "0,-1,-2").split(",")]
DTMAX = int(os.environ.get("DTMAX", "20"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", str(len(POWERS))))
T0 = int(os.environ.get("T0", "2"))
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0").split(",")]
TMAXPLOT = int(os.environ.get("TMAXPLOT", "18"))


def one_config(k):
    V, windows = dc.load_peram_windows(k)
    with h5py.File(dc.PERAM_DIR + "peram.%d.h5" % k, "r") as f:
        evals = f["evals"][:]                     # (Nt, Nv) = mu_a^2 per timeslice (D~^H D~ eigenvalues)
    nop = len(POWERS)
    C = np.zeros((nop, nop, DTMAX))
    nwin = 0
    for (tsrc0, tau, taugw) in windows:
        twin = tau.shape[0]
        Nv = tau.shape[-1]
        # mode-space diagonal vertex diag(mu_a^{2p}) at each ABSOLUTE timeslice, per power
        Phi = {p: [np.diag(evals[tsrc0 + a] ** p) for a in range(twin)] for p in POWERS}
        for dt in range(DTMAX):
            s0s = [s for s in range(twin) if s + dt < twin]
            if not s0s:
                continue
            for ia, pa in enumerate(POWERS):
                for ib, pb in enumerate(POWERS):
                    acc = 0.0
                    for s in s0s:
                        t = s + dt
                        acc += -np.trace(Phi[pa][t] @ tau[t, s] @ Phi[pb][s] @ tau[s, t]).real
                    C[ia, ib, dt] += acc / len(s0s)
        nwin += 1
    return C / nwin


def gevp(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    Vp = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, Vp, T0), Vp


def main():
    tag = dc.ENS.split("nu0")[0]
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    e22 = float(os.environ.get("E22", "0.556"))   # free (2,2)=2E_{3/2} reference
    print("# ENS=%s L=%d  mu-weight single-meson GEVP  POWERS=%s  %d cfg  reb%d@%d T0=%d  m_PS=%.4f 2m_PS=%.4f"
          % (tag, dc.L, POWERS, len(ks), NKEEP, REBT, T0, mps, m2ps))
    allC = np.array([one_config(k) for k in ks])
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = gevp(blk.mean(0), None)
    if nb < 2:                                    # single-config (free testbed): central value only
        em_e = np.zeros_like(em_c)
    else:
        ems = np.array([gevp(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
        em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]
    nk = em_c.shape[1]                            # actual states resolved (GEVP may rank-truncate below NKEEP)

    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(nk)))
    for t in range(T0, min(tmax, 16)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(nk))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1.1, alpha=0.7)
    ax.text(TMAXPLOT * 0.7, m2ps + 0.01, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
    ax.axhline(mps, color="firebrick", ls=":", lw=1.2, alpha=0.8)
    ax.text(TMAXPLOT * 0.7, mps + 0.01, r"$m_{PS}=%.4f$" % mps, fontsize=9, color="firebrick")
    if dc.ENS.startswith("free"):
        ax.axhline(e22, color="tab:purple", ls="-.", lw=1.0, alpha=0.7)
        ax.text(TMAXPLOT * 0.7, e22 + 0.01, r"$(2,2)=%.3f$" % e22, fontsize=9, color="tab:purple")
    cols = ["firebrick", "tab:purple", "tab:blue", "tab:green"]
    mkr = ["o", "D", "s", "^"]
    lab = [r"state 0 ($m_{PS}$)", r"state 1 ($(2,2)$?)", "state 2", "state 3"]
    for n in range(nk):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n % 4], marker=mkr[n % 4], ms=6, lw=1.2,
                    capsize=3, label=lab[n] if n < len(lab) else "state %d" % n)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, TMAXPLOT)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"$\mu$-weight single-meson GEVP  %s L%d %dcfg  powers=%s"
                 % (tag, dc.L, ncfg, POWERS), fontsize=9)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_muweight_gevp_%s_L%d_claude.png" % (tag, dc.L)
    fig.savefig(out, dpi=140)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
