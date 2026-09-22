#!/usr/bin/env python3
# t00_stress_claude.py  [CHUNK 1: fermionic stress-tensor T_00 two-point from the distillation perambulator]
#
# Run (free-limit validation, single config, exact):
#   ENS=free NVDIR=distill_Nv24 LREF=1 BINSIZE=1 python3 t00_stress_claude.py
#
# Operator (2-component S^2 x R, ell=0, 0^{++}, \Delta=3):
#   O_T00(t) = sum_x w_x psibar(x,t) sigma_3 (d0^sym psi)(x,t),
#     w_x = dual_areas[x] * Y00 ,  Y00 = 1/sqrt(4 pi) ,  sigma_3 = diag(1,-1) (temporal gamma),
#     d0^sym psi(x,t) = (1/2)[ psi(x,t+1) - psi(x,t-1) ] .
# This is an INTERPOLATOR for the state (energy E), NOT the exactly-conserved lattice T_munu.
#
# Connected two-point (single distillation loop; sink at s+dt, source at s, translation-avg over s):
#   C(dt) = -(1/4) sum_{sx,sy=+-1} sx sy sum_{x,y} w_x w_y
#             Tr_spin[ sigma_3 S(x,s+dt+sx; y,s) sigma_3 S(y,s+sy; x,s+dt) ] ,
#   S(x,ta; y,tb) = AblkS(ta,tb)[x,:,y,:] = V(ta) tau(ta,tb) V(tb)^dag  (fs_gevp_point make_config).
# d0^sym differences the FIRST (out) time index of each propagator block.
# CONTACT: AblkS already subtracts the ultralocal GW contact (A - 1/2 I) at ta==tb, so the dt=1,sx=-1
#   term (which hits AblkS(s,s)) uses the contact-subtracted block -- no special-casing (NM, 2026-09-16).
#
# Vacuum: print <O_T00> (expect ~0 by temporal parity); subtract <O>^2 only if it is nonzero.
# Refs: distillation Peardon 0905.2160; derivative distillation operators (Peardon-Edwards); T_munu \Delta=d=3.
# See t00_stress_impl_plan_claude.md.

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

NS = dc.NS
BINSIZE = int(os.environ.get("BINSIZE", "1"))
DTMAX = int(os.environ.get("DTMAX", "16"))
NCFG = int(os.environ.get("NCFG", "0"))            # 0 = all configs in PERAM_DIR
ROVER = 1.0 / 0.189                                # placeholder 1/R for the free reference lines (a_t m)

SIG3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

# Symmetric first-derivative stencils d0 psi(t) = sum_i c_i psi(t + sh_i), applied on the out-time of each
# propagator leg.  Each is antisymmetric (sh -> -sh flips sign) so it projects the genuine \Delta=3 primary.
#   order 2 : (1/2)[ psi_{t+1} - psi_{t-1} ]                          , error O(a_t^2)   [validated free ->0.567]
#   order 4 : (1/12)[ -psi_{t+2} + 8 psi_{t+1} - 8 psi_{t-1} + psi_{t-2} ] , error O(a_t^4)  [NM 2026-09-16]
# Higher-order = smaller discretization shift of the free plateau off 3/R (NO extra operators, single op).
STENCILS = {
    2: [(1, 0.5), (-1, -0.5)],
    4: [(2, -1.0 / 12.0), (1, 8.0 / 12.0), (-1, -8.0 / 12.0), (-2, 1.0 / 12.0)],
}
STENCIL = int(os.environ.get("STENCIL", "2"))      # 2 = O(a^2) default; 4 = O(a^4) improved


def corr_one_config(k, w, stencil=None):
    # Full connected T_00 two-point C(dt) for one config, plus the one-loop vacuum <O>(t).
    # d0^sym differences the FIRST (out) time index of each block via the chosen finite-difference stencil.
    if stencil is None:
        stencil = STENCILS[STENCIL]
    smax = max(abs(sh) for sh, _ in stencil)
    AblkS, AblkSt, twin, nsite, U, tau = fg.make_config(k)
    C = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        s_lo = max(smax, smax - dt)                 # need s-smax>=0 and s+dt-smax>=0
        s_hi = twin - 1 - dt - smax                 # need s+dt+smax<=twin-1 (also s+smax<=twin-1 for dt>=0)
        if s_hi < s_lo:
            continue
        acc = 0.0
        cnt = 0
        for s in range(s_lo, s_hi + 1):
            val = 0.0
            for sha, ca in stencil:
                B1 = AblkS(s + dt + sha, s)          # [x, a, y, b] = S(x,s+dt+sha ; y,s)
                for shb, cb in stencil:
                    B2 = AblkS(s + shb, s + dt)      # [y, c, x, d] = S(y,s+shb ; x,s+dt)
                    # sum_{x,y} w_x w_y  sigma3_{ab} B1[x,b,y,c] sigma3_{cd} B2[y,d,x,a]
                    t = np.einsum("ab,pbqc,cd,qdpa,p,q->", SIG3, B1, SIG3, B2, w, w, optimize=True)
                    val += (ca * cb) * t
            acc += val
            cnt += 1
        C[dt] = (-acc / cnt).real                    # overall -1 (fermion loop); stencil carries the 1/2, 1/12
    # vacuum one-loop <O>(t) = -sum_x w_x sum_i c_i Tr[ sigma3 S(x,t+sh_i ; x,t) ]
    Ovac = []
    for t in range(smax, twin - smax):
        vv = 0.0
        for sh, c in stencil:
            B = AblkS(t + sh, t)                     # [x, a, x, b]
            vv += c * np.einsum("ab,pbpa->", SIG3, B, optimize=True)
        Ovac.append((-vv).real)
    Ovac = float(np.mean(Ovac))
    return C, Ovac


def main():
    tag = dc.ENS
    dual = dc.dual_areas_from_mesh()
    w = dual * dc.Y00
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    if not ks:
        print("# no perambulators found in %s" % dc.PERAM_DIR)
        return
    print("# ENS=%s NVDIR=%s  nsite=%d  ncfg=%d  DTMAX=%d BINSIZE=%d"
          % (tag, dc.NVDIR, dual.shape[0], len(ks), DTMAX, BINSIZE))

    allC = []
    Ovac_all = []
    for k in ks:
        C, Ovac = corr_one_config(k, w)
        allC.append(C)
        Ovac_all.append(Ovac)
    allC = np.array(allC)                            # (ncfg, DTMAX)
    Ovac_mean = float(np.mean(Ovac_all))
    print("# <O_T00> (one-loop vacuum, config-avg) = % .6e   (expect ~0, temporal parity)" % Ovac_mean)

    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    Cm = blk.mean(0)

    def effmass(Cv):
        with np.errstate(all="ignore"):
            return np.log(Cv[:-1] / Cv[1:])

    em_c = effmass(Cm)
    if nb < 2:
        em_err = np.zeros_like(em_c)                 # single config (free limit): exact, no jackknife
    else:
        ems = np.array([effmass(np.delete(blk, i, 0).mean(0)) for i in range(nb)])
        em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))

    print("\n#  dt |    C(dt)        m_eff(err)      [free refs: sigma=2/R=%.3f  T00=3/R=%.3f  2m=4/R=%.3f]"
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
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    for val, lab, col in ((2.0 / ROVER, r"$\sigma=2/R$", "tab:green"),
                          (3.0 / ROVER, r"$T_{00}=3/R$", "tab:red"),
                          (4.0 / ROVER, r"$2m=4/R$", "gray")):
        ax.axhline(val, color=col, ls="--", lw=1, alpha=0.7)
        ax.text(ts[-1] * 0.02, val + 0.006, lab, fontsize=9, color=col)
    g = np.isfinite(em_c) & np.isfinite(em_err) & (em_err < 0.3)
    ax.errorbar(ts[g], em_c[g], yerr=em_err[g], color="tab:red", marker="o", ms=5, lw=1.1,
                capsize=2.5, label=r"$T_{00}$ effmass")
    ax.set_ylim(0.0, 0.9)
    ax.set_xlim(1, DTMAX - 1)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"Fermionic $T_{00}$ two-point  %s  %d cfg" % (tag, ncfg), fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_stress_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
