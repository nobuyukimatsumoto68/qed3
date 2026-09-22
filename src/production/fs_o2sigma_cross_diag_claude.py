#!/usr/bin/env python3
# fs_o2sigma_cross_diag_claude.py  [FULL cross <O_sigmasigma^FS(t,t+d) sigma_FS^2(s)>: connected vs disconnected]
# Run:  ENS=... NVDIR=distill_Nv24 SPLIT=1 python3 fs_o2sigma_cross_diag_claude.py
#
# Brute-force Wick sum over the 4! = 24 pairings of the four sigma vertices p=(t, t+d, s, s):
#   value(perm) = (-1)^{#cycles} prod_cycles Tr[ Phi(p_c0) leg(p_c0,p_c1) Phi(p_c1) ... leg(p_ck,p_c0) ]
#   equal-time contact subtracted on coincident legs:  S-part  leg(a,a)=tau(a,a) - 1/2 I ;
#                                                       Stilde  leg(a,a)= -1/2 ( tau'(a,a)+tau(a,a) ).
#   FS = S-part (leg=tau) + Stilde-part (leg=-tau').  (tau'=tau_gw; the tilde-S minus cancels in even loops.)
# CLASSIFY each perm: CONNECTED iff a cycle contains both a sink vertex {0,1}=(t,t+d) and a source {2,3}=s
#   -> connected = A,B,C,D,E,G analogs (keep) ; disconnected = F,H,I,J (drop).
# Plot connected / disconnected / full totals, LINEAR, t-sum(plateau) subtracted per jackknife.
# See fs_sigma2_diagram_note_claude.md.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
from itertools import permutations
import numpy as np
import distill_contract_claude as dc

SPLIT = int(os.environ.get("SPLIT", "1"))
DTPLOT = int(os.environ.get("DTPLOT", "16"))          # plotted range
DTMAX = int(os.environ.get("DTMAX", "30"))            # computed range (needs a plateau tail)
PLAT_LO = int(os.environ.get("PLAT_LO", "20"))        # plateau (t-sum) window start

# precompute the cycle decomposition + connectivity of all 24 permutations of 4 vertices
PERMS = []
for pi in permutations(range(4)):
    seen = [False] * 4
    cycles = []
    for i in range(4):
        if seen[i]:
            continue
        cyc = []
        j = i
        while not seen[j]:
            seen[j] = True
            cyc.append(j)
            j = pi[j]
        cycles.append(cyc)
    # connected iff some cycle mixes sink {0,1} and source {2,3}
    conn = any(any(v in (0, 1) for v in c) and any(v in (2, 3) for v in c) for c in cycles)
    PERMS.append((cycles, conn))


def per_config(k, w00):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Iv = np.eye(tau.shape[-1])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    def legS(a, b):
        return tau[a, b] - 0.5 * Iv if a == b else tau[a, b]

    def legSt(a, b):
        return -0.5 * (taugw[a, a] + tau[a, a]) if a == b else -taugw[a, b]

    def wick(times, leg):
        # sum over the 24 perms; return (connected, disconnected)
        c = 0.0
        d = 0.0
        for cycles, conn in PERMS:
            val = (-1.0) ** len(cycles)
            for cyc in cycles:
                m = Phi[times[cyc[0]]]
                for idx in range(len(cyc)):
                    a = times[cyc[idx]]
                    b = times[cyc[(idx + 1) % len(cyc)]]
                    m = m @ leg(a, b)
                    if idx + 1 < len(cyc):
                        m = m @ Phi[b]
                val = val * np.trace(m)
            if conn:
                c += val
            else:
                d += val
        return c.real, d.real

    dd = SPLIT
    Cc = np.full(DTMAX, np.nan)
    Cd = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        ss = [s for s in range(twin) if s + dt + dd < twin]
        if not ss:
            continue
        accc = 0.0
        accd = 0.0
        for s in ss:
            t1 = s + dt
            t2 = s + dt + dd
            times = [t1, t2, s, s]                        # sink 0,1 ; source 2,3
            cS, dS = wick(times, legS)
            cSt, dSt = wick(times, legSt)
            accc += cS + cSt
            accd += dS + dSt
        Cc[dt] = accc / len(ss)
        Cd[dt] = accd / len(ss)
    return Cc, Cd


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    print("# ENS=%s ncfg=%d  FULL cross <O_sigmasigma^FS sigma_FS^2>  connected vs disconnected  SPLIT=%d"
          % (tag, len(dc.KS), SPLIT))
    allC = []
    allD = []
    for k in dc.KS:
        Cc, Cd = per_config(k, w00)
        allC.append(Cc)
        allD.append(Cd)
    allC = np.array(allC)
    allD = np.array(allD)
    ncfg = allC.shape[0]

    def jk(C):
        plat = np.nanmean(C[:, PLAT_LO:], axis=1, keepdims=True)
        Cc = C - plat
        n = C.shape[0]
        samp = np.array([np.delete(Cc, i, 0).mean(0) for i in range(n)])
        return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    cm, ce = jk(allC)
    dm, de = jk(allD)
    tm, te = jk(allC + allD)

    print("\n#  dt |  CONNECTED(A,B,C,D,E,G)   DISCONNECTED(F,H,I,J)   FULL")
    for dt in range(1, DTPLOT):
        print("#  %2d | %11.4e(%.1e)  %11.4e(%.1e)  %11.4e(%.1e)"
              % (dt, cm[dt], ce[dt], dm[dt], de[dt], tm[dt], te[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTPLOT)
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
    panels = [((cm, ce), "CONNECTED  (A,B,C,D,E,G) -- keep", "tab:green", "^"),
              ((dm, de), "DISCONNECTED  (F,H,I,J) -- drop", "tab:red", "o"),
              ((tm, te), "FULL (conn+disc)", "black", "*")]
    for ax, (cc, lab, col, mk) in zip(axs, panels):
        ax.errorbar(dts, cc[0][dts], yerr=cc[1][dts], color=col, marker=mk, ms=5, lw=1, capsize=2)
        ax.axhline(0.0, color="gray", lw=0.8, alpha=0.6)
        ax.set_title(lab, fontsize=10)
        ax.set_xlabel(r"$dt$", fontsize=9)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[0].set_ylabel(r"$\langle O_{\sigma\sigma}^{FS}(s{+}dt)\,\sigma_{FS}^2(s)\rangle$", fontsize=10)
    fig.suptitle(r"FS cross $\langle O_{\sigma\sigma}\,\sigma_{FS}^2\rangle$: connected vs disconnected (linear, plateau-sub)  %s L1 %d cfg $\delta$=%d"
                 % (tag, ncfg, SPLIT), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs("figs", exist_ok=True)
    out = "figs/fs_o2sigma_cross_diag_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
