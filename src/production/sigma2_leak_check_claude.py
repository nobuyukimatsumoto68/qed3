#!/usr/bin/env python3
# sigma2_leak_check_claude.py -- direct single-meson LEAK diagnostic for the sigma^2 (P+) channel under basis
#   truncation.  Reuses sigma2_mPS_gevp_claude.one_config (C_11=<sigma_00 sigma_00>, C_12 triangle,
#   C_22=<sigma^2_00 sigma^2_00>, flavfac PP, tau ONLY).  Reports THREE effmasses so the leak is unambiguous:
#     (a) C_11-only  (1x1)  -> the single meson m_PS (sanity: perams reproduce m_PS).
#     (b) C_22-only  (1x1)  -> THE LEAK TEST.  Protected basis -> ground = TWO-MESON (~2m_PS).  Leaky (truncated
#                              one-sided) basis -> ground COLLAPSES to m_PS.  This is NOT masked by any GEVP.
#     (c) {sigma_00, sigma^2} 2x2 GEVP -> state0 = m_PS, state1 = two-meson (the augmented picture).
#   The point: the 2x2 GEVP can RESCUE the two-meson (assign m_PS to sigma_00), so (b) is the clean diagnostic.
#
#   BASIS TRUNCATION via NVKEEP (env, read inside distill_contract loader): keep the lowest NVKEEP modes.
#   A/B at L1 (complete Nv=24):  NVKEEP=12  _sym (protected)  vs  distill_Nv24 one-sided (leaks).
#   Run: ENS=..L1.. LREF=1 NVDIR=distill_Nv24_sym NVKEEP=12 T0=2 REBT=4 BINSIZE=8 python3 sigma2_leak_check_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import hankel_rebase_scan_claude as hs
import sigma2_mPS_gevp_claude as MG

T0 = int(os.environ.get("T0", "2"))
REBT = int(os.environ.get("REBT", "4"))
BINSIZE = int(os.environ.get("BINSIZE", "8"))
NCFG = int(os.environ.get("NCFG", "0"))
NVKEEP = int(os.environ.get("NVKEEP", "0"))


def effmass_1x1(Cdt, V0):
    # 1x1 "GEVP" = plain rebased effmass of a single correlator (Hankel off=[0], rebase at REBT, 1 state).
    Cts = Cdt[:, None, None]
    Big = hs.hankel_off(Cts, [0])
    V = hs.staged_project(Big, [(REBT, 1)], T0) if V0 is None else V0
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = MG.G.antipodal_map()
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    print("# ENS=%s L=%d NVDIR=%s NVKEEP=%d  %d cfg  reb%d T0=%d  m_PS=%.4f 2m_PS=%.4f"
          % (tag, dc.L, os.environ.get("NVDIR"), NVKEEP, len(ks), REBT, T0, mps, m2ps))

    allC = np.array([MG.one_config(k, dualf, wY, Pmap) for k in ks])
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    Cmean = blk.mean(0)

    # jackknife over bins for each of the three estimators
    def run(kind):
        if kind == "C11":
            getC = lambda M: M[0, 0]
        elif kind == "C22":
            getC = lambda M: M[1, 1]
        if kind in ("C11", "C22"):
            emc, V0 = effmass_1x1(getC(Cmean), None)
        else:  # 2x2 GEVP
            emc, V0 = MG.gevp(Cmean, None)
        if nb < 2:                                   # single-config (free testbed): central value only, no jk
            return emc, np.zeros_like(emc)
        if kind in ("C11", "C22"):
            ems = np.array([effmass_1x1(getC(np.delete(blk, i, 0).mean(0)), V0)[0] for i in range(nb)])
        else:
            ems = np.array([MG.gevp(np.delete(blk, i, 0).mean(0), V0)[0] for i in range(nb)])
        eme = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
        return emc, eme

    em11, ee11 = run("C11")
    em22, ee22 = run("C22")
    emg, eeg = run("GEVP")

    print("\n#  t |  C11-only(m_PS?)   C22-only(LEAK TEST)   GEVP s0        GEVP s1")
    tmax = min(em11.shape[0], 16)
    for t in range(T0, tmax):
        def fmt(mc, me, n=0):
            v = mc[t, n] if mc.ndim == 2 else mc[t]
            e = me[t, n] if me.ndim == 2 else me[t]
            return "%7.4f(%.4f)" % (v, e) if np.isfinite(v) and np.isfinite(e) else "    --       "
        print("#  %2d | %s   %s   %s   %s"
              % (t, fmt(em11, ee11), fmt(em22, ee22), fmt(emg, eeg, 0), fmt(emg, eeg, 1)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(em11.shape[0])
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1, alpha=0.7)
    ax.text(ts[-1] * 0.02 + T0, m2ps + 0.008, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
    ax.axhline(mps, color="firebrick", ls=":", lw=1.2, alpha=0.8)
    ax.text(ts[-1] * 0.02 + T0, mps + 0.008, r"$m_{PS}=%.4f$" % mps, fontsize=9, color="firebrick")
    g11 = np.isfinite(em11[:, 0]) & (ee11[:, 0] < 0.3)
    g22 = np.isfinite(em22[:, 0]) & (ee22[:, 0] < 0.3)
    ax.errorbar(ts[g11], em11[g11, 0], yerr=ee11[g11, 0], color="tab:red", marker="o", ms=5, lw=1.1, capsize=2.5,
                label=r"$C_{11}$ only ($\sigma_{00}$, = $m_{PS}$)")
    ax.errorbar(ts[g22] + 0.08, em22[g22, 0], yerr=ee22[g22, 0], color="tab:blue", marker="s", ms=5, lw=1.1,
                capsize=2.5, label=r"$C_{22}$ only ($\sigma^2_{00}$) -- LEAK TEST")
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(em11.shape[0], int(os.environ.get("TMAXPLOT", "14"))))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"sigma^2 leak check  %s L%d  NVDIR=%s NVKEEP=%d  reb%d T0=%d %dcfg"
                 % (tag, dc.L, os.environ.get("NVDIR"), NVKEEP, REBT, T0, ncfg), fontsize=9)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_leak_check_%s_L%d_%s_nv%d_claude.png" % (tag, dc.L, os.environ.get("NVDIR"), NVKEEP)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
