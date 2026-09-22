#!/usr/bin/env python3
# t00_stress_ham_interacting_claude.py
# INTERACTING fermionic T_00 (naive e.sigma energy density, r=0) via distillation elementals.
#   O_H(t) = eta^H W(t) xi + xi^H W(t)^H eta,  W(t)_{ij} = 0.5 kappa_ij (e^a sigma_a)(ij) Omega(ij) exp(i u_ij(t)),
#   u_ij(t) = sign_ij * theta_sp[t, il_ij]  (raw HMC spatial link phase; primal_links_n1_claude.dat gives il/sign/kappa).
# Elemental contraction (single loop; G=tau, backward D_ov^{-dag} -> bar_tau=delta-tau; off-diag dt>=1 -> tau):
#   Phi_W(t)  = V(t)^H W(t)  V(t) ,  Phi_WH(t) = V(t)^H W(t)^H V(t)   (Nv x Nv)
#   C(dt) = -(1/N) sum_a { Tr[Phi_W(a+dt) tau(a+dt,a) Phi_W(a) tau(a,a+dt)]
#                        + Tr[Phi_WH(a+dt) tau(a+dt,a) Phi_WH(a) tau(a,a+dt)] }   (a=window index)
# gamma/Omega from omega_n1.dat/alpha_n1.dat (SpinStructureSimp files -> matches D_ov exactly).
# VALIDATE=1: free peram + theta=0 must reproduce t00_stress_ham (position-space) up to the uniform kappa.
# Run interacting:  ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000 \
#                   NVDIR=distill_Nv24 LREF=1 python3 t00_stress_ham_interacting_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "final/analysis_axial")
import numpy as np
import distill_contract_claude as dc
import geom_hopping_claude as gh
import effmass_axial_tp_l3_perm_hankel_claude as hk

NS = dc.NS
S0, S1, S2, S3 = gh.S0, gh.S1, gh.S2, gh.S3
DTMAX = int(os.environ.get("DTMAX", "24"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
KMIN = int(os.environ.get("KMIN", "20"))
NT = int(os.environ.get("NT", "128"))
R = float(os.environ.get("R", "0.0"))
GEOM = os.environ.get("GEOM", "../../geometry/data/")
LREF = int(os.environ.get("LREF", "1"))
CONFIGDIR = os.environ.get("CONFIGDIR", dc.ENS)        # ckpoint_lat.<k> dir (no data_ prefix)
VALIDATE = int(os.environ.get("VALIDATE", "0"))
OFFS = [int(x) for x in os.environ.get("OFFSETS", "0,3").split(",")]
REBT = int(os.environ.get("REBT", "3"))
NKEEP = int(os.environ.get("NKEEP", "1"))
T0 = int(os.environ.get("T0", "2"))
WLO = int(os.environ.get("WLO", "6"))                  # plateau fit window [WLO, WHI]
WHI = int(os.environ.get("WHI", "10"))
ROVER = 1.0 / 0.189


def wconst(em, w, lo, hi):
    # weighted-constant fit (diagonal weights w) over t in [lo,hi]; nan-safe
    ts = np.arange(lo, hi + 1)
    e = em[ts]
    wt = w[ts]
    g = np.isfinite(e) & np.isfinite(wt) & (wt > 0)
    if g.sum() < 1:
        return np.nan
    return np.sum(wt[g] * e[g]) / np.sum(wt[g])


def load_link_table(path):
    tab = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            s = line.split()
            if len(s) < 5:
                continue
            tab[(int(s[0]), int(s[1]))] = (int(s[2]), int(s[3]), float(s[4]))
    return tab


def read_config_sp(path, n_links):
    a = np.fromfile(path, dtype="<f8")
    return a[:NT * n_links].reshape(NT, n_links)          # sp[t, il]  (spatial block first)


def build_W_gauge(om, alpha, nns, nsite, tab, sp_t, r=R):
    N = NS * nsite
    W = np.zeros((N, N), complex)
    for i in range(nsite):
        diag = np.zeros((2, 2), complex)
        for j in nns[i]:
            il, sgn, kap = tab[(i, j)]
            ph = np.exp(1j * sgn * sp_t[il])
            W[NS * i:NS * i + 2, NS * j:NS * j + 2] += 0.5 * kap * ((-r) * S0 + gh.gamma(alpha, i, j)) @ gh.Omega(om, i, j) * ph
            diag += 0.5 * r * kap * S0
        W[NS * i:NS * i + 2, NS * i:NS * i + 2] += diag
    return W


def corr_one_config(k, om, alpha, nns, nsite, tab, sp_override=None):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    if sp_override is not None:
        sp = sp_override
    else:
        sp = read_config_sp("%s/ckpoint_lat.%d" % (CONFIGDIR, k), len(tab) // 2)
    PhiW = []
    PhiWH = []
    for a in range(twin):
        t = tsrc0 + a
        Vt = V[t].T                                        # (2Ns, Nv)
        W = build_W_gauge(om, alpha, nns, nsite, tab, sp[t])
        PhiW.append(Vt.conj().T @ W @ Vt)
        PhiWH.append(Vt.conj().T @ W.conj().T @ Vt)
    C = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        na = twin - dt
        if na <= 0:
            continue
        acc = 0.0
        for a in range(na):
            b = a + dt
            acc += np.trace(PhiW[b] @ tau[b, a] @ PhiW[a] @ tau[a, b])
            acc += np.trace(PhiWH[b] @ tau[b, a] @ PhiWH[a] @ tau[a, b])
        C[dt] = (-acc / na).real
    return C


def main():
    dual = dc.dual_areas_from_mesh()
    om, alpha, nns, nsite = gh.build(GEOM, LREF)
    tab = load_link_table("primal_links_n%d_claude.dat" % LREF)
    n_links = len(tab) // 2
    print("# ENS=%s NVDIR=%s nsite=%d n_links=%d  CONFIGDIR=%s" % (dc.ENS, dc.NVDIR, nsite, n_links, CONFIGDIR))

    if VALIDATE:
        # free peram + theta=0: elemental vs position-space t00_stress_ham (build_W r=0)
        import t00_stress_ham_claude as th
        import t00_wilson_kernel_claude as wk
        th.DTMAX = DTMAX
        k = dc.KS[0]
        sp0 = np.zeros((NT, n_links))
        Cel = corr_one_config(k, om, alpha, nns, nsite, tab, sp_override=sp0)
        Wpos, _, _ = wk.build_W(GEOM, LREF, r=R)
        Cpos = th.corr_one_config(k, Wpos)
        print("# VALIDATE (free, theta=0): elemental vs position-space (ratio should be constant = kappa^2):")
        for dt in (3, 5, 8, 11):
            print("#   dt=%2d  C_el=% .5e  C_pos=% .5e  ratio=%.6f" % (dt, Cel[dt], Cpos[dt], Cel[dt] / Cpos[dt]))
        return

    ks = [k for k in dc.KS if k >= KMIN]
    if NCFG:
        ks = ks[:NCFG]
    tag = dc.ENS.split("nu0")[0]
    cdir = os.environ.get("CACHEDIR", "t00_ham_cache_claude")
    os.makedirs(cdir, exist_ok=True)
    cache = "%s/allC_%s_%s_n%d_R%g_dt%d_km%d_claude.npy" % (cdir, tag.replace(".", "p"), dc.NVDIR, len(ks), R, DTMAX, KMIN)
    if os.path.exists(cache):
        allC = np.load(cache)
        print("# loaded cached correlators <- %s  (%d cfg)" % (cache, allC.shape[0]))
    else:
        allC = []
        for i, k in enumerate(ks):
            allC.append(corr_one_config(k, om, alpha, nns, nsite, tab))
            if (i + 1) % 50 == 0:
                print("#   ... %d/%d configs" % (i + 1, len(ks)))
        allC = np.array(allC)
        np.save(cache, allC)
        print("# cached correlators -> %s" % cache)
    ncfg = allC.shape[0]
    sgn = np.sign(np.nanmean(allC[:, 3]))
    allC = sgn * allC

    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vop = hk.hankel_effmass_scalar(blk.mean(0), OFFS, REBT, NKEEP, T0, 0.2)
    ems = []
    for i in range(nb):
        Ci = np.delete(blk, i, 0).mean(0)
        ems.append(hk.hankel_effmass_scalar(Ci, OFFS, REBT, NKEEP, T0, 0.2, Vop=Vop)[0])
    ems = np.array(ems)
    em_err = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, axis=0))

    # point-to-point too
    Cm = blk.mean(0)
    with np.errstate(all="ignore"):
        emp = np.log(Cm[:-1] / Cm[1:])

    # weighted-constant plateau fit over [WLO,WHI]; stat via jackknife (same central weights), sys via +1 shift
    wt = 1.0 / np.where(em_err[:, 0] > 0, em_err[:, 0], np.inf) ** 2
    m0 = wconst(em_c[:, 0], wt, WLO, WHI)
    mjk = np.array([wconst(ems[i, :, 0], wt, WLO, WHI) for i in range(nb)])
    stat = np.sqrt((nb - 1) * np.nanmean((mjk - np.nanmean(mjk)) ** 2))
    m0_sh = wconst(em_c[:, 0], wt, WLO + 1, WHI + 1)
    sys = abs(m0 - m0_sh)
    comb = np.sqrt(stat ** 2 + sys ** 2)

    print("# ncfg=%d  bin%d nb=%d  off%s reb%d@%d T0=%d" % (ncfg, BINSIZE, nb, "-".join(map(str, OFFS)), NKEEP, REBT, T0))
    print("# PLATEAU FIT [%d,%d] (weighted const): a_t*m = %.4f (stat %.4f)(sys %.4f) -> comb %.4f  [free 0.567]"
          % (WLO, WHI, m0, stat, sys, comb))
    print("\n#  t | Hankel m0(err)   point2point   [T00 free=%.3f]" % (3.0 / ROVER))
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
    ax_level = os.environ.get("AXIAL", "")               # a_t*m level of the axial conserved current
    if ax_level:
        aL = float(ax_level)
        ax.axhline(aL, color="tab:purple", ls=":", lw=1.3, alpha=0.8)
        ax.text(0.2, aL + 0.01, r"axial current $=%.3f$" % aL, fontsize=9, color="tab:purple")
        ax.axhline(1.5 * aL, color="tab:green", ls="-.", lw=1.3, alpha=0.8)
        ax.text(0.2, 1.5 * aL + 0.01, r"$1.5\times$axial $=%.3f$" % (1.5 * aL), fontsize=9, color="tab:green")
    g = np.isfinite(em_c[:, 0]) & (em_err[:, 0] < 0.4)
    ax.errorbar(ts[g], em_c[g, 0], yerr=em_err[g, 0], color="tab:red", marker="o", ms=5, lw=1.2, capsize=2.5,
                label="block-Hankel ground")
    tp = np.arange(emp.shape[0])
    gp = np.isfinite(emp)
    ax.plot(tp[gp], emp[gp], "-", color="lightgray", lw=1.0, marker=".", ms=4, label="point-to-point")
    ax.fill_between([WLO, WHI], m0 - comb, m0 + comb, color="tab:blue", alpha=0.2)
    ax.plot([WLO, WHI], [m0, m0], color="tab:blue", lw=1.6,
            label=r"fit [%d,%d]: $%.3f(%.3f)$" % (WLO, WHI, m0, comb))
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(0, em_c.shape[0])
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"Interacting $T_{00}$ (naive $e\cdot\sigma$)  %s  %d cfg" % (dc.ENS.split("nu0")[0], ncfg),
                 fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_stress_ham_interacting_%s_claude.png" % dc.ENS.split("nu0")[0]
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
