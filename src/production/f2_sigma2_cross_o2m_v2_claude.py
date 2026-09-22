#!/usr/bin/env python3
# f2_sigma2_cross_o2m_v2_claude.py
#   Cross  <F^2(s+dt) O_2m(s)>_c  with the EXISTING antipodal two-meson interpolator
#       O_2m(s) = sum_x A_x sigma(x,s) sigma(P(x),s)   (fs_channels_v2 / fs_gevp_point op 1, EQUAL-TIME).
#   F^2 gluonic -> O_2m's four fermions self-contract into a single connected EQUAL-TIME antipodal loop
#   (the disc tadpole product dies under improvement):
#       L_2m(s) = - sum_x A_x Tr_spin[ P(x,s ; P(x),s) P(P(x),s ; x,s) ]
#   position propagator P = AblkS(s,s) (improved -1/2 on the site-diagonal; antipodal x!=P(x) off-diagonal,
#   so no contact).  Antipodal pairing (x and P(x) maximally separated on S^2) is the best TWO-MESON
#   projector (suppresses the one-meson).  F^2 still cannot reach the 4-fermion two-meson STATE (sea-mediated
#   only); this probes the gluonic admixture of the two-meson-projected single-loop 0++ density.
#
#   Clean t-fold (C(dt)=C(-dt), F^2 global).  Per-config t-sum + jackknife plateau (vacuum) subtraction.
#   Run: ENS=Nf2_gsq1.000000...L1 NVDIR=distill_Nv24 python3 f2_sigma2_cross_o2m_v2_claude.py

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
import fs_gevp_point_claude as G

Nt = 128
OP_F = int(os.environ.get("OP_F", "0"))
DTCALC = int(os.environ.get("DTCALC", "56"))
DTMAX = int(os.environ.get("DTMAX", "40"))
TSUM_LO = int(os.environ.get("TSUM_LO", "1"))
PLAT_LO = int(os.environ.get("PLAT_LO", str(DTCALC - 8)))


def gluedir():
    return "data_" + dc.ENS


def load_of(k):
    with h5py.File("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k), "r") as g:
        return np.array(g["O"])[OP_F].astype(float)     # (Nt,)


def o2m_loops(k, dual, Pmap):
    # per-source connected O_2m self-loop L_2m(s) (equal-time antipodal) ; also tsrc0, twin
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(k)
    idx = np.arange(nsite)
    V, _, _, tsrc0, _ = dc.load_peram(k)
    L = np.full(twin, np.nan)
    for s in range(twin):
        Pf = AblkS(s, s)                                 # (nsite,NS,nsite,NS) equal-time
        B1 = Pf[idx, :, Pmap, :]                          # (nsite,NS,NS): P(x ; P(x))
        B2 = Pf[Pmap, :, idx, :]                          # (nsite,NS,NS): P(P(x) ; x)
        L[s] = -(dual * np.einsum('xab,xba->x', B1, B2)).sum().real
    return L, tsrc0, twin


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh().astype(float)
    Pmap = G.antipodal_map()
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k))]
    print("# ENS=%s matched cfg=%d/%d OP_F=%d  (F^2 x antipodal-equal-time O_2m loop)"
          % (tag, len(ks), len(dc.KS), OP_F))

    allC = []
    for k in ks:
        L, tsrc0, twin = o2m_loops(k, dual, Pmap)
        of = load_of(k)
        svalid = np.arange(twin)
        Cd = np.zeros(DTCALC)
        for dt in range(DTCALC):
            idxf = (tsrc0 + svalid + dt) % Nt
            idxb = (tsrc0 + svalid - dt) % Nt
            ofd = 0.5 * (of[idxf] + of[idxb])
            Cd[dt] = np.mean(ofd * L[svalid])
        allC.append(Cd)
    allC = np.array(allC)                                 # (ncfg, DTCALC)
    ncfg = allC.shape[0]

    ts = allC - allC[:, TSUM_LO:].mean(1, keepdims=True)
    n = ts.shape[0]
    samp = np.array([np.delete(ts, i, 0).mean(0) for i in range(n)])
    samp = samp - samp[:, PLAT_LO:DTCALC].mean(1, keepdims=True)
    cm = samp.mean(0)
    ee = np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    print("\n#  dt |   C(err)             S/N")
    for dt in range(0, DTMAX):
        sn = cm[dt] / ee[dt] if ee[dt] > 0 else 0.0
        print("#  %2d | % .3e(%.1e)  %6.2f" % (dt, cm[dt], ee[dt], sn))
    print("# max|S/N| = %.2f" % np.max(np.abs(cm[:DTMAX] / ee[:DTMAX])))

    # jackknife effmass on the decaying part; lattice a_t m, and physical m = a_t m / a_t (a_t=0.2 -> x5)
    at = float(os.environ.get("AT", "0.2"))
    with np.errstate(all="ignore"):
        em_samp = np.log(samp[:, :-1] / samp[:, 1:])          # valid where same-sign & decaying
    em_c = em_samp.mean(0)
    em_e = np.sqrt((n - 1) * np.mean((em_samp - em_c) ** 2, 0))
    print("\n#  dt | a_t m_eff(err)   ->  m_phys=a_t m/a_t (err)   [F2 glue: a_t m=0.616 lat / 3.08 phys]")
    for dt in range(0, 12):
        if np.isfinite(em_c[dt]) and np.isfinite(em_e[dt]):
            print("#  %2d | %7.4f(%.4f)      ->  %6.3f(%.3f)"
                  % (dt, em_c[dt], em_e[dt], em_c[dt] / at, em_e[dt] / at))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(0, DTMAX)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(dts, np.abs(cm[dts]), yerr=ee[dts], color="tab:blue", marker="s", ms=5, lw=1, capsize=2)
    ax.set_yscale("log")
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$|\langle F^2(s{+}dt)\,O_{2m}(s)\rangle_c|$")
    ax.grid(alpha=0.3)
    ax.set_title(r"$F^2$ x antipodal equal-time $O_{2m}$ (two-meson interpolator), folded, vac-sub  %s L1 %d cfg"
                 % (tag, ncfg), fontsize=11)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_sigma2_cross_o2m_v2_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
