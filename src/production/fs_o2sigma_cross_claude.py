#!/usr/bin/env python3
# fs_o2sigma_cross_claude.py  [FS GEVP off-diagonal  <O_sigmasigma^FS(t) sigma_FS^2(s)>, LINEAR]
# Run:  ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000 \
#         NVDIR=distill_Nv24 SPLIT=1 python3 fs_o2sigma_cross_claude.py
#
# O_sigmasigma^FS(t) = sigma_FS,00(t) sigma_FS,00(t+delta)  (time-split two-meson interpolator, delta=SPLIT).
# By the S/tilde-S non-mixing rule (4 sigma_FS vertices) the connected cross splits into all-tau + all-(-tau'):
#   <O_sigmasigma^FS(t) sigma_FS^2(s)>_c = 2 ( M[tau](t,s) M[tau](t+d,s) + M[-tau'](t,s) M[-tau'](t+d,s) ),
#   M[leg](a,b) = -Tr[ Phi(a) leg(a,b) Phi(b) leg(b,a) ]   (single-sigma meson propagator).
# All legs off-diagonal (t,t+d != s for dt>=1) -> NO equal-time contact, automatically vacuum-free.
# The overall tilde-S minus cancels (2 legs per M): M[-tau']=-Tr[Phi tau' Phi tau'], tau'=tau_gw.
# See fs_sigma2_diagram_note_claude.md.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

SPLIT = int(os.environ.get("SPLIT", "1"))
DTMAX = int(os.environ.get("DTMAX", "20"))


def meson_matrix(Phi, leg, twin):
    # M[a,b] = -Tr[ Phi(a) leg(a,b) Phi(b) leg(b,a) ]  (single-sigma meson, sink a, source b)
    M = np.zeros((twin, twin))
    for a in range(twin):
        for b in range(twin):
            M[a, b] = (-np.trace(Phi[a] @ leg[a, b] @ Phi[b] @ leg[b, a])).real
    return M


def per_config(k, w00):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    Mtau = meson_matrix(Phi, tau, twin)
    Mtaup = meson_matrix(Phi, taugw, twin)              # M[-tau'] = M[+tau'] (2 legs, sign cancels)
    d = SPLIT
    S = np.full(DTMAX, np.nan)
    St = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        ss = [s for s in range(twin) if s + dt + d < twin]
        if not ss:
            continue
        S[dt] = np.mean([Mtau[s + dt, s] * Mtau[s + dt + d, s] for s in ss])
        St[dt] = np.mean([Mtaup[s + dt, s] * Mtaup[s + dt + d, s] for s in ss])
    return S, St


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    print("# ENS=%s ncfg=%d  <O_sigmasigma^FS . sigma_FS^2> cross  SPLIT=%d" % (tag, len(dc.KS), SPLIT))
    allS = []
    allSt = []
    for k in dc.KS:
        S, St = per_config(k, w00)
        allS.append(S)
        allSt.append(St)
    allS = np.array(allS)
    allSt = np.array(allSt)
    ncfg = allS.shape[0]

    def jk(C):
        n = C.shape[0]
        samp = np.array([np.delete(C, i, 0).mean(0) for i in range(n)])
        return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    Sm, Se = jk(allS)
    Stm, Ste = jk(allSt)
    totm, tote = jk(2.0 * (allS + allSt))

    print("\n#  dt |  S-part M[tau]^2     Stilde M[tau']^2    TOTAL 2(S+Stilde)")
    for dt in range(0, DTMAX):
        print("#  %2d | %11.4e(%.1e) %11.4e(%.1e) %11.4e(%.1e)"
              % (dt, Sm[dt], Se[dt], Stm[dt], Ste[dt], totm[dt], tote[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(0, DTMAX)
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
    panels = [((Sm, Se), r"$S$-part  $M[\tau]M[\tau]$", "tab:green", "^"),
              ((Stm, Ste), r"$\tilde S$-part  $M[\tau']M[\tau']$", "tab:purple", "D"),
              ((totm, tote), r"TOTAL  $2(S+\tilde S)$", "black", "*")]
    for ax, (cc, lab, col, mk) in zip(axs, panels):
        g = np.isfinite(cc[0])
        ax.errorbar(dts[g], cc[0][g], yerr=cc[1][g], color=col, marker=mk, ms=5, lw=1, capsize=2)
        ax.axhline(0.0, color="gray", lw=0.8, alpha=0.6)
        ax.set_title(lab, fontsize=11)
        ax.set_xlabel(r"$dt$", fontsize=9)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[0].set_ylabel(r"$\langle O_{\sigma\sigma}^{FS}(s{+}dt)\,\sigma_{FS}^2(s)\rangle_c$", fontsize=10)
    fig.suptitle(r"FS cross $\langle O_{\sigma\sigma}\,\sigma_{FS}^2\rangle$ (linear, vacuum-free)  %s L1 %d cfg  $\delta$=%d"
                 % (tag, ncfg, SPLIT), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs("figs", exist_ok=True)
    out = "figs/fs_o2sigma_cross_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)

    # LOG overlay (dt>=1, positive region); distinct color+marker (color-blind safe)
    figL, axL = plt.subplots(figsize=(8, 5.4))
    series = [((Sm, Se), r"$S$-part $M[\tau]^2$", "tab:green", "^"),
              ((Stm, Ste), r"$\tilde S$-part $M[\tau']^2$", "tab:purple", "D"),
              ((totm, tote), r"TOTAL $2(S+\tilde S)$", "black", "*")]
    for (cc, lab, col, mk) in series:
        g = np.isfinite(cc[0]) & (cc[0] > 0)
        g[0] = False
        axL.errorbar(dts[g], cc[0][g], yerr=cc[1][g], color=col, marker=mk, ms=5, lw=1, capsize=2, label=lab)
    axL.set_yscale("log")
    axL.set_xlabel(r"$dt$")
    axL.set_ylabel(r"$\langle O_{\sigma\sigma}^{FS}(s{+}dt)\,\sigma_{FS}^2(s)\rangle_c$")
    axL.set_title(r"FS cross $\langle O_{\sigma\sigma}\,\sigma_{FS}^2\rangle$ (LOG, $dt\geq1$, CONNECTED only)  %s L1 %d cfg  $\delta$=%d"
                  % (tag, ncfg, SPLIT), fontsize=11)
    axL.legend(fontsize=9)
    axL.grid(alpha=0.3, which="both")
    figL.tight_layout()
    outL = "figs/fs_o2sigma_cross_log_%s_claude.png" % tag
    figL.savefig(outL, dpi=130)
    plt.close(figL)
    print("# -> %s" % outL)


if __name__ == "__main__":
    main()
