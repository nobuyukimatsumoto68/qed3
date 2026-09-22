#!/usr/bin/env python3
# t00_stress_ham_claude.py  [CHUNK B: fermionic T_00 two-point from the HAMILTONIAN (spatial Wilson) vertex]
#
# Precise energy density (first-order action, NM 2026-09-16): O_T00(t) = eta^H(t) W xi(t) + xi^H(t) W^H eta(t),
# W = spatial Wilson kernel (t00_wilson_kernel_claude.build_W; Eq IV.1 first line, M5=0).  NOT a temporal
# derivative.  Connected two-point (single loop; only cross-propagators nonzero):
#   <xi eta^H> = D_ov^{-1}       = G      = AblkS(a,b)                      (forward)
#   <eta xi^H> = D_ov^{-dag}     = Gtld   = -AblkS(a,b)  (off-diag; = 1/2 I - AblkS at equal time)   (backward)
#   C_T00(t,t0) = -Tr[W G(t,t0) W G(t0,t)] - Tr[W^H Gtld(t,t0) W^H Gtld(t0,t)] , translation-avg over t0.
# Propagator legs verified against the peram code by "Fin: Two-meson" (T2b/T2c); the (backward=-AblkS) identity
# is bar_tau = delta - tau.  See t00_hamiltonian_derivation_claude.md.
#
# Run (free-limit, single config, exact):  ENS=free NVDIR=distill_Nv24 LREF=1 python3 t00_stress_ham_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as fg
import t00_wilson_kernel_claude as wk

NS = dc.NS
DTMAX = int(os.environ.get("DTMAX", "28"))
BINSIZE = int(os.environ.get("BINSIZE", "1"))
NCFG = int(os.environ.get("NCFG", "0"))
GEOM = os.environ.get("GEOM", "../../geometry/data/")
LREF = int(os.environ.get("LREF", "1"))
# ENERGY-DENSITY vertex = the NAIVE e.sigma hop only (R=0): W = 0.5 kappa (e^a sigma_a) Omega.  The Wilson
# -r 1 term (R=1) is an O(a^2), spin-scalar, doubler-removal artifact that overlaps sigma (Delta=2, 0.378) and
# is DROPPED (NM 2026-09-16).  R=0 lands the free plateau on 3/R=0.567 (EOM-identical to psibar sigma_3 D_t psi).
R = float(os.environ.get("R", "0.0"))
ROVER = 1.0 / 0.189


def corr_one_config(k, W):
    AblkS, AblkSt, twin, nsite, U, tau = fg.make_config(k)
    N = NS * nsite
    Iden = np.eye(N)
    Wm = W
    WmH = W.conj().T

    def G(a, b):
        return AblkS(a, b).reshape(N, N)

    def Gtld(a, b):
        g = -AblkS(a, b).reshape(N, N)          # off-diagonal: D_ov^{-dag} = -AblkS
        if a == b:
            g = 0.5 * Iden - AblkS(a, b).reshape(N, N)   # equal-time contact: 1/2 I - AblkS
        return g

    C = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        s_hi = twin - 1 - dt                    # need t0=s and t=s+dt in [0, twin-1]
        if s_hi < 0:
            continue
        acc = 0.0
        cnt = 0
        for s in range(0, s_hi + 1):
            t0 = s
            t = s + dt
            termA = -np.trace(Wm @ G(t, t0) @ Wm @ G(t0, t))
            termB = -np.trace(WmH @ Gtld(t, t0) @ WmH @ Gtld(t0, t))
            acc += termA + termB
            cnt += 1
        C[dt] = (acc / cnt).real
    return C


def main():
    tag = dc.ENS
    W, nns, nsite = wk.build_W(GEOM, LREF, r=R)          # R=0: naive e.sigma energy vertex (Wilson term dropped)
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    if not ks:
        print("# no perambulators found in %s" % dc.PERAM_DIR)
        return
    print("# ENS=%s NVDIR=%s  nsite=%d N=%d  ncfg=%d  DTMAX=%d" % (tag, dc.NVDIR, nsite, NS * nsite, len(ks), DTMAX))

    allC = np.array([corr_one_config(k, W) for k in ks])
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    Cm = blk.mean(0)

    # sign: energy-density connected 2pt should be positive-decaying; report the sign we find
    sgn = np.sign(Cm[2]) if np.isfinite(Cm[2]) else 1.0
    Cpos = sgn * Cm

    def effmass(Cv):
        with np.errstate(all="ignore"):
            return np.log(Cv[:-1] / Cv[1:])

    em_c = effmass(Cpos)
    if nb < 2:
        em_err = np.zeros_like(em_c)
    else:
        ems = np.array([effmass(sgn * np.delete(blk, i, 0).mean(0)) for i in range(nb)])
        em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))

    print("# sign(C)=%+d  (Cpos = sign*C, positive-decaying expected)" % int(sgn))
    print("\n#  dt |    C(dt)         m_eff(err)   [free refs sigma=%.3f T00=%.3f 2m=%.3f]"
          % (2.0 / ROVER, 3.0 / ROVER, 4.0 / ROVER))
    for dt in range(1, DTMAX - 1):
        if not np.isfinite(Cm[dt]):
            continue
        me = em_c[dt] if dt < em_c.shape[0] else np.nan
        ee = em_err[dt] if dt < em_err.shape[0] else np.nan
        print("#  %2d | % .6e   %8.4f(%.4f)" % (dt, Cm[dt], me, ee))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(em_c.shape[0])
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    for val, lab, col in ((2.0 / ROVER, r"$\sigma=2/R$", "tab:green"),
                          (3.0 / ROVER, r"$T_{00}=3/R$", "tab:red"),
                          (4.0 / ROVER, r"$2m=4/R$", "gray")):
        ax.axhline(val, color=col, ls="--", lw=1, alpha=0.7)
        ax.text(1.2, val + 0.006, lab, fontsize=9, color=col)
    g = np.isfinite(em_c) & np.isfinite(em_err) & (em_err < 0.3)
    ax.errorbar(ts[g], em_c[g], yerr=em_err[g], color="tab:red", marker="o", ms=5, lw=1.1,
                capsize=2.5, label=r"$T_{00}$ (Hamiltonian $W$)")
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(1, DTMAX - 1)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"Fermionic $T_{00}$ from the Hamiltonian (spatial Wilson) vertex  %s  %d cfg" % (tag, ncfg),
                 fontsize=10.5)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_stress_ham_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
