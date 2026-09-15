#!/usr/bin/env python3
# diag_sum_vs_mPS_claude.py  [effmass of a W10-weighted sum of chosen diagrams vs m_PS and 2 m_PS]
# Run:  ENS=... NVDIR=distill_Nv24 CHANNEL=PS CONTACT=0.5 DIAGS=0,2,3,6 python3 diag_sum_vs_mPS_claude.py
#   DIAGS = comma-separated diagram indices (0=A,1=B,2=C,3=D,4=E,5=F,6=G,7=H,8=I,9=J).
#   e.g. DIAGS=0,2,3,6 -> A+C+D+G (one-meson-rich); DIAGS=1,4 -> B+E (two-meson).
# Summed diagram = sum_i W10[i]*diag_i (physical weight in <PS^2 PS^2>).  Plateau (t-sum) subtracted,
# jackknife.  m_PS from the single-sigma correlator -Tr[Phi(t) tau(t,s) Phi(s) tau(s,t)].

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
DIAGS = [int(x) for x in os.environ.get("DIAGS", "0,2,3,6").split(",")]


def main():
    de.CONTACT = float(os.environ.get("CONTACT", "0.5"))
    tag = dc.ENS.split("nu0")[0]
    channel = os.environ.get("CHANNEL", "PS")
    lab = "+".join("ABCDEFGHIJ"[i] for i in DIAGS)
    print("# ENS=%s ncfg=%d  effmass of (%s) vs m_PS, 2m_PS  CONTACT=%.2f" % (tag, len(dc.KS), lab, de.CONTACT))
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
    ncfg = allD.shape[0]

    def jk_effmass(C, subtract_plateau):
        if subtract_plateau:
            C = C - C[:, PLAT_LO:].mean(1, keepdims=True)
        n = C.shape[0]
        samp = np.array([np.delete(C, i, 0).mean(0) for i in range(n)])
        with np.errstate(all="ignore"):
            em = np.log(samp[:, :-1] / samp[:, 1:])
        return em.mean(0), np.sqrt((n - 1) * np.mean((em - em.mean(0)) ** 2, 0))

    Wsum = np.zeros_like(allD[:, 0, :])
    for i in DIAGS:
        Wsum += dc.W10[i] * allD[:, i, :]
    emS, eeS = jk_effmass(Wsum, True)
    emPS, eePS = jk_effmass(allCs, False)

    print("\n#  dt |  m_PS(err)     2m_PS       (%s)(err)" % lab)
    for dt in range(1, DTMAX):
        print("#  %2d |  %6.4f(%.4f)  %6.4f    %6.4f(%.4f)"
              % (dt, emPS[dt], eePS[dt], 2 * emPS[dt], emS[dt], eeS[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    ax.errorbar(dts, emS[dts], yerr=eeS[dts], color="tab:red", marker="o", ms=4, lw=1, capsize=2,
                label=r"$(%s)$ effmass" % lab)
    ax.errorbar(dts, emPS[dts], yerr=eePS[dts], color="tab:blue", marker="s", ms=3, lw=1, capsize=2,
                label=r"$m_{PS}(dt)$")
    ax.errorbar(dts + 0.05, 2 * emPS[dts], yerr=2 * eePS[dts], color="tab:green", marker="D", ms=3, lw=1,
                ls="--", capsize=2, label=r"$2 m_{PS}(dt)$")
    ax.set_ylim(0.0, 1.3)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"$(%s)$ vs $m_{PS}$, $2m_{PS}$   <PS^2 PS^2>  %s L1 %d cfg" % (lab, tag, ncfg))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/diag_sum_%s_c%.1f_%s_claude.png" % (lab.replace("+", ""), de.CONTACT, tag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
