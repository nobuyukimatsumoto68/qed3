#!/usr/bin/env python3
# t00_stress_temporal_interacting_claude.py
# INTERACTING time-displaced T_00: O_T(t) = eta^H M D_t^cov xi + xi^H M D_t^cov eta,  M=diag(w) sigma_3,
#   D_t^cov xi(x,s) = (1/2)[ e^{-i theta_tp(s,x)} xi(x,s+1) - e^{+i theta_tp(s-1,x)} xi(x,s-1) ]
# (covariant temporal derivative; phase convention from dirac_ext.h:177-178, naive sigma_3 part).  theta_tp is
# the TEMPORAL block of ckpoint_lat (Nt x n_sites doubles after the spatial block).  signP=signM=1 (interior).
# Contraction (position-space AblkS; both terms; G=AblkS, Gt=-AblkS off-diag / 1/2 I - AblkS eq-time):
#   V_sx(a) = diag_x(w_x e^{i phi(sx,a,x)}) sigma_3 ,  phi(+1,a,x)=-tp[t,x], phi(-1,a,x)=+tp[t-1,x]  (t=tsrc0+a)
#   C(dt) = -(1/4) mean_s sum_{sx,sy} sx sy { Tr[V_sx(t) G(t+sx,t0) V_sy(t0) G(t0+sy,t)]
#                                            + Tr[V_sx(t) Gt(t+sx,t0) V_sy(t0) Gt(t0+sy,t)] } ,  t=t0+dt.
# VALIDATE=1: free peram + tp=0 must reproduce t00_stress_temporal (position-space).
# Run:  ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000 NVDIR=distill_Nv24 \
#       LREF=1 python3 t00_stress_temporal_interacting_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "final/analysis_axial")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as fg
import effmass_axial_tp_l3_perm_hankel_claude as hk

NS = dc.NS
SIG3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
DTMAX = int(os.environ.get("DTMAX", "24"))
BINSIZE = int(os.environ.get("BINSIZE", "80"))
NCFG = int(os.environ.get("NCFG", "0"))
KMIN = int(os.environ.get("KMIN", "20"))
NT = int(os.environ.get("NT", "128"))
NSITES = int(os.environ.get("NSITES", "12"))
NLINKS = int(os.environ.get("NLINKS", "30"))
CONFIGDIR = os.environ.get("CONFIGDIR", dc.ENS)
VALIDATE = int(os.environ.get("VALIDATE", "0"))
OFFS = [int(x) for x in os.environ.get("OFFSETS", "0,4").split(",")]
REBT = int(os.environ.get("REBT", "3"))
NKEEP = int(os.environ.get("NKEEP", "1"))
T0 = int(os.environ.get("T0", "2"))
WLO = int(os.environ.get("WLO", "7"))
WHI = int(os.environ.get("WHI", "10"))
ROVER = 1.0 / 0.189


def read_config_tp(path):
    a = np.fromfile(path, dtype="<f8")
    return a[NT * NLINKS:NT * NLINKS + NT * NSITES].reshape(NT, NSITES)     # tp[t, ix] (temporal block)


def build_V(w, tp, tsrc0, twin):
    # V_plus[a], V_minus[a] (N,N) for window index a (absolute t=tsrc0+a); phase per dirac_ext.h:177-178
    nsite = w.shape[0]
    N = NS * nsite
    Vp = []
    Vm = []
    for a in range(twin):
        t = tsrc0 + a
        php = np.exp(-1j * tp[t])                 # forward: e^{-i tp[t]}
        phm = np.exp(+1j * tp[(t - 1) % NT])      # backward: e^{+i tp[t-1]}
        Ap = np.zeros((N, N), complex)
        Am = np.zeros((N, N), complex)
        for i in range(nsite):
            Ap[NS * i:NS * i + 2, NS * i:NS * i + 2] = w[i] * php[i] * SIG3
            Am[NS * i:NS * i + 2, NS * i:NS * i + 2] = w[i] * phm[i] * SIG3
        Vp.append(Ap)
        Vm.append(Am)
    return Vp, Vm


def corr_one_config(k, w, tp_override=None):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    if tp_override is not None:
        tp = tp_override
    else:
        tp = read_config_tp("%s/ckpoint_lat.%d" % (CONFIGDIR, k))
    AblkS, AblkSt, twin2, nsite, U, tau2 = fg.make_config(k)
    N = NS * nsite
    Iden = np.eye(N)
    Vp, Vm = build_V(w, tp, tsrc0, twin)

    def G(a, b):
        return AblkS(a, b).reshape(N, N)

    def Gt(a, b):
        g = -AblkS(a, b).reshape(N, N)
        if a == b:
            g = 0.5 * Iden - AblkS(a, b).reshape(N, N)
        return g

    def Vpick(a, sgn):
        return Vp[a] if sgn == 1 else Vm[a]

    C = np.full(DTMAX, np.nan)
    signs = (1, -1)
    for dt in range(DTMAX):
        s_lo = 1
        s_hi = twin - 2 - dt
        if s_hi < s_lo:
            continue
        acc = 0.0
        cnt = 0
        for s in range(s_lo, s_hi + 1):
            t0 = s
            t = s + dt
            for sx in signs:
                B1 = G(t + sx, t0)
                Gt1 = Gt(t + sx, t0)
                Vsx = Vpick(t, sx)
                for sy in signs:
                    B2 = G(t0 + sy, t)
                    Gt2 = Gt(t0 + sy, t)
                    Vsy = Vpick(t0, sy)
                    acc += sx * sy * (np.trace(Vsx @ B1 @ Vsy @ B2) + np.trace(Vsx @ Gt1 @ Vsy @ Gt2))
            cnt += 1
        C[dt] = (-0.25 * acc / cnt).real
    return C


def wconst(em, w, lo, hi):
    hi = min(hi, em.shape[0] - 1)
    ts = np.arange(lo, hi + 1)
    if ts.size < 1:
        return np.nan
    e = em[ts]
    wt = w[ts]
    g = np.isfinite(e) & np.isfinite(wt) & (wt > 0)
    if g.sum() < 1:
        return np.nan
    return np.sum(wt[g] * e[g]) / np.sum(wt[g])


def main():
    dual = dc.dual_areas_from_mesh()
    w = dual * dc.Y00
    print("# ENS=%s NVDIR=%s nsite=%d  CONFIGDIR=%s" % (dc.ENS, dc.NVDIR, dual.shape[0], CONFIGDIR))

    if VALIDATE:
        import t00_stress_temporal_claude as tt
        tt.DTMAX = DTMAX
        k = dc.KS[0]
        tp0 = np.zeros((NT, NSITES))
        Cnew = corr_one_config(k, w, tp_override=tp0)
        CA, CB = tt.corr_one_config(k, tt.build_M(dual))
        Cold = CA + CB
        print("# VALIDATE (free, tp=0): interacting-driver vs t00_stress_temporal (ratio should be 1):")
        for dt in (3, 5, 8, 11):
            print("#   dt=%2d  C_new=% .5e  C_old=% .5e  ratio=%.6f" % (dt, Cnew[dt], Cold[dt], Cnew[dt] / Cold[dt]))
        return

    ks = [k for k in dc.KS if k >= KMIN]
    if NCFG:
        ks = ks[:NCFG]
    allC = []
    for i, k in enumerate(ks):
        allC.append(corr_one_config(k, w))
        if (i + 1) % 50 == 0:
            print("#   ... %d/%d" % (i + 1, len(ks)))
    allC = np.array(allC)
    ncfg = allC.shape[0]
    sgn = np.sign(np.nanmean(allC[:, 3]))
    allC = sgn * allC

    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vop = hk.hankel_effmass_scalar(blk.mean(0), OFFS, REBT, NKEEP, T0, 0.2)
    ems = np.array([hk.hankel_effmass_scalar(np.delete(blk, i, 0).mean(0), OFFS, REBT, NKEEP, T0, 0.2, Vop=Vop)[0]
                    for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, axis=0))
    wt = 1.0 / np.where(em_err[:, 0] > 0, em_err[:, 0], np.inf) ** 2
    m0 = wconst(em_c[:, 0], wt, WLO, WHI)
    mjk = np.array([wconst(ems[i, :, 0], wt, WLO, WHI) for i in range(nb)])
    stat = np.sqrt((nb - 1) * np.nanmean((mjk - np.nanmean(mjk)) ** 2))
    sysv = abs(m0 - wconst(em_c[:, 0], wt, WLO + 1, WHI + 1))
    comb = np.sqrt(stat ** 2 + sysv ** 2)

    Cm = blk.mean(0)
    with np.errstate(all="ignore"):
        emp = np.log(Cm[:-1] / Cm[1:])
    print("# ncfg=%d bin%d nb=%d off%s reb%d@%d T0=%d  sign=%+d"
          % (ncfg, BINSIZE, nb, "-".join(map(str, OFFS)), NKEEP, REBT, T0, int(sgn)))
    print("# PLATEAU FIT [%d,%d]: a_t*m = %.4f (stat %.4f)(sys %.4f) -> comb %.4f  [free T00=0.567]"
          % (WLO, WHI, m0, stat, sysv, comb))
    print("\n#  t | Hankel m0(err)   point2point")
    for t in range(em_c.shape[0]):
        if np.isfinite(em_c[t, 0]):
            pp = emp[t] if t < emp.shape[0] else np.nan
            print("#  %2d | %7.4f(%.4f)   %7.4f" % (t, em_c[t, 0], em_err[t, 0], pp))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(em_c.shape[0])
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    ax.axhline(3.0 / ROVER, color="tab:red", ls="--", lw=1, alpha=0.6)
    ax.text(0.2, 3.0 / ROVER + 0.01, r"free $T_{00}=0.567$", fontsize=9, color="tab:red")
    ax_level = os.environ.get("AXIAL", "")
    if ax_level:
        aL = float(ax_level)
        ax.axhline(aL, color="tab:purple", ls=":", lw=1.3, alpha=0.8)
        ax.text(0.2, aL + 0.01, r"axial current $=%.3f$" % aL, fontsize=9, color="tab:purple")
        ax.axhline(1.5 * aL, color="tab:green", ls="-.", lw=1.3, alpha=0.8)
        ax.text(0.2, 1.5 * aL + 0.01, r"$1.5\times$axial $=%.3f$" % (1.5 * aL), fontsize=9, color="tab:green")
    g = np.isfinite(em_c[:, 0]) & (em_err[:, 0] < 0.4)
    ax.errorbar(ts[g], em_c[g, 0], yerr=em_err[g, 0], color="tab:blue", marker="o", ms=5, lw=1.2, capsize=2.5,
                label="block-Hankel ground")
    tp_ = np.arange(emp.shape[0])
    gp = np.isfinite(emp)
    ax.plot(tp_[gp], emp[gp], "-", color="lightgray", lw=1.0, marker=".", ms=4, label="point-to-point")
    ax.fill_between([WLO, WHI], m0 - comb, m0 + comb, color="tab:blue", alpha=0.2)
    ax.plot([WLO, WHI], [m0, m0], color="tab:blue", lw=1.6, label=r"fit [%d,%d]: $%.3f(%.3f)$" % (WLO, WHI, m0, comb))
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(0, em_c.shape[0])
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"Interacting $T_{00}$ time-displaced $O_T$  %s  %d cfg" % (dc.ENS.split("nu0")[0], ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_stress_temporal_interacting_%s_claude.png" % dc.ENS.split("nu0")[0]
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
