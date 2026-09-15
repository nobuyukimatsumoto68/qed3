#!/usr/bin/env python3
# diag_BE_vs_2mPS_claude.py  [effmass of (B+E) diagrams vs 2 m_PS]
# Run:  ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000 \
#         NVDIR=distill_Nv24 CHANNEL=PS CONTACT=0.5 python3 diag_BE_vs_2mPS_claude.py
#
# B = -T_S (diagram 1), E = C_S^2 (diagram 4) of <PS^2 PS^2>.  Sum B+E, extract effmass (plateau/t-sum
# subtracted, jackknife), and compare to 2 m_PS where m_PS is the single-sigma meson mass from
# C_sigma(dt) = <sigma(t) sigma(0)>_conn = -Tr[Phi(t) tau(t,s) Phi(s) tau(s,t)] (the CS building block).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

PLAT_LO = int(os.environ.get("PLAT_LO", str(de.PLAT_LO)))
DTMAX = int(os.environ.get("DTMAX", "24"))


def main():
    de.CONTACT = float(os.environ.get("CONTACT", "0.5"))
    tag = dc.ENS.split("nu0")[0]
    channel = os.environ.get("CHANNEL", "PS")
    print("# ENS=%s ncfg=%d  (B+E) effmass vs 2 m_PS  CHANNEL=%s CONTACT=%.2f" % (tag, len(dc.KS), channel, de.CONTACT))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allD = []
    allCs = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        leg = tau if channel == "PS" else -taugw
        D = np.zeros((10, twin))
        Cs = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            cs = 0.0
            for s in range(ns):
                t = s + dt
                acc += de.diags_pair(Phi, leg, s, t)
                cs += (-np.trace(Phi[t] @ leg[t, s] @ Phi[s] @ leg[s, t])).real
            D[:, dt] = (acc / ns).real
            Cs[dt] = cs / ns
        allD.append(D)
        allCs.append(Cs)
    allD = np.array(allD)
    allCs = np.array(allCs)
    ncfg, _, twin = allD.shape

    def jk_effmass(C, subtract_plateau):
        if subtract_plateau:
            C = C - C[:, PLAT_LO:].mean(1, keepdims=True)
        n = C.shape[0]
        samp = np.array([np.delete(C, i, 0).mean(0) for i in range(n)])
        with np.errstate(all="ignore"):
            em = np.log(samp[:, :-1] / samp[:, 1:])
        return em.mean(0), np.sqrt((n - 1) * np.mean((em - em.mean(0)) ** 2, 0)), samp.mean(0)

    BE = allD[:, 1, :] + allD[:, 4, :]              # B + E
    emBE, eeBE, cBE = jk_effmass(BE, True)
    emPS, eePS, cPS = jk_effmass(allCs, False)      # single-sigma; no plateau (genuine meson, ->0)

    # pick m_PS plateau (average over a mid window)
    lo, hi = 4, 9
    mPS = np.average(emPS[lo:hi], weights=1.0 / eePS[lo:hi] ** 2)
    mPSe = 1.0 / np.sqrt(np.sum(1.0 / eePS[lo:hi] ** 2))
    print("# m_PS (single-sigma, window [%d,%d)) = %.4f(%.4f)   -> 2 m_PS = %.4f(%.4f)"
          % (lo, hi, mPS, mPSe, 2 * mPS, 2 * mPSe))

    print("\n#  dt |  m_PS(err)      (B+E)(err)     [2 m_PS = %.4f]" % (2 * mPS))
    for dt in range(1, DTMAX):
        sBE = "%6.4f(%.4f)" % (emBE[dt], eeBE[dt]) if np.isfinite(emBE[dt]) else "   ---   "
        sPS = "%6.4f(%.4f)" % (emPS[dt], eePS[dt]) if np.isfinite(emPS[dt]) else "   ---   "
        print("#  %2d |  %s   %s" % (dt, sPS, sBE))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    ax.errorbar(dts, emPS[dts], yerr=eePS[dts], color="tab:blue", marker="s", ms=4, lw=1, capsize=2,
                label=r"$m_{PS}$ (single $\sigma$)")
    ax.errorbar(dts + 0.05, emBE[dts], yerr=eeBE[dts], color="tab:red", marker="o", ms=4, lw=1, capsize=2,
                label=r"$B+E$ effmass")
    # matched-dt curve 2 m_PS(dt) (excited-state contamination partially cancels vs B+E at the same dt)
    ax.errorbar(dts, 2 * emPS[dts], yerr=2 * eePS[dts], color="tab:green", marker="D", ms=3, lw=1,
                ls="--", capsize=2, label=r"$2\, m_{PS}(dt)$ (matched)")
    ax.set_ylim(0.0, 1.3)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("(B+E) vs $2 m_{PS}$   <PS^2 PS^2>  %s L1 %d cfg contact=%.2f" % (tag, ncfg, de.CONTACT))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/diag_BE_vs_2mPS_c%.1f_%s_claude.png" % (de.CONTACT, tag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
