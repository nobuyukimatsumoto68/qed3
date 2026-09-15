#!/usr/bin/env python3
# diag_vacsub_sigma2_claude.py  [vacuum-subtracted <sigma^2 sigma^2> via double plateau subtraction]
# Run:  ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000 \
#         NVDIR=distill_Nv24 CHANNEL=PS CONTACT=0.5 python3 diag_vacsub_sigma2_claude.py
#
# Vacuum subtraction of the full <PS^2 PS^2> = 2 sum_i W10[i] diag_i, done PER JACKKNIFE ENSEMBLE:
#   (1) per-diagram t-sum (plateau) subtraction:  diag_i -> diag_i - <diag_i>_{dt>=PLAT_LO}
#   (2) sum W10-weighted (x2), then a final plateau subtraction on the TOTAL.
# The two-step removes each diagram's own constant and any residual constant in the sum -> the
# vacuum-subtracted sigma^2 two-point function.  Compared to m_PS, 2m_PS, and the connected A..E,G sum.

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
    print("# ENS=%s ncfg=%d  vacuum-subtracted <sigma^2 sigma^2> (double plateau sub, per jk)  CONTACT=%.2f"
          % (tag, len(dc.KS), de.CONTACT))
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
    allD = np.array(allD)                                   # (ncfg,10,twin)
    allCs = np.array(allCs)
    ncfg, _, twin = allD.shape

    def vacsub_total(Djk):
        # (1) per-diagram plateau (t-sum) subtraction ; (2) W10 sum x2 ; (3) total plateau subtraction
        Dsub = Djk - Djk[:, PLAT_LO:].mean(1, keepdims=True)
        tot = 2.0 * (dc.W10[:, None] * Dsub).sum(0)
        tot = tot - tot[PLAT_LO:].mean()
        return tot

    # central + jackknife of the vacuum-subtracted total, and effmass
    Tcen = vacsub_total(allD.mean(0))
    samp = np.array([vacsub_total(np.delete(allD, i, 0).mean(0)) for i in range(ncfg)])
    with np.errstate(all="ignore"):
        emJ = np.log(samp[:, :-1] / samp[:, 1:])
    emT = emJ.mean(0)
    eeT = np.sqrt((ncfg - 1) * np.mean((emJ - emT) ** 2, 0))
    cT = samp.mean(0)
    ceT = np.sqrt((ncfg - 1) * np.mean((samp - cT) ** 2, 0))

    # single-sigma m_PS (no plateau; genuine meson)
    sPS = np.array([np.delete(allCs, i, 0).mean(0) for i in range(ncfg)])
    with np.errstate(all="ignore"):
        emPSj = np.log(sPS[:, :-1] / sPS[:, 1:])
    emPS = emPSj.mean(0)
    eePS = np.sqrt((ncfg - 1) * np.mean((emPSj - emPS) ** 2, 0))

    print("\n#  dt |  C_vacsub(err)      m_eff(err)     m_PS     2m_PS")
    for dt in range(1, DTMAX):
        print("#  %2d |  %10.3e(%.1e)  %6.4f(%.4f)  %6.4f  %6.4f"
              % (dt, cT[dt], ceT[dt], emT[dt], eeT[dt], emPS[dt], 2 * emPS[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.4))
    # left: vacuum-subtracted correlator (log)
    pos = cT[dts] > 0
    a1.errorbar(dts[pos], cT[dts][pos], yerr=ceT[dts][pos], color="tab:red", marker="o", ms=4, lw=1, capsize=2)
    a1.set_yscale("log")
    a1.set_xlabel(r"$dt$")
    a1.set_ylabel(r"$C^{\sigma^2}_\mathrm{vac-sub}(dt)$")
    a1.set_title("vacuum-subtracted $\\langle\\sigma^2\\sigma^2\\rangle$ (log)")
    # right: effmass vs m_PS, 2m_PS
    a2.errorbar(dts, emT[dts], yerr=eeT[dts], color="tab:red", marker="o", ms=4, lw=1, capsize=2,
                label=r"$\sigma^2$ vac-sub effmass")
    a2.errorbar(dts, emPS[dts], yerr=eePS[dts], color="tab:blue", marker="s", ms=3, lw=1, capsize=2,
                label=r"$m_{PS}(dt)$")
    a2.errorbar(dts + 0.05, 2 * emPS[dts], yerr=2 * eePS[dts], color="tab:green", marker="D", ms=3, lw=1,
                ls="--", capsize=2, label=r"$2 m_{PS}(dt)$")
    a2.set_ylim(0.0, 1.3)
    a2.set_xlabel(r"$dt$")
    a2.set_ylabel(r"$a_t m_\mathrm{eff}$")
    a2.set_title("effmass vs $m_{PS}$, $2m_{PS}$")
    a2.legend(fontsize=9)
    fig.suptitle("Vacuum-subtracted $\\langle\\sigma^2\\sigma^2\\rangle$  %s L1 %d cfg contact=%.2f"
                 % (tag, ncfg, de.CONTACT), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs("figs", exist_ok=True)
    out = "figs/diag_vacsub_sigma2_c%.1f_%s_claude.png" % (de.CONTACT, tag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
