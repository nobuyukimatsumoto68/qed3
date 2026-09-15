#!/usr/bin/env python3
# fs_gevp_connected_claude.py  [FS 2-op connected GEVP {sigma_FS^2, O_sigmasigma^FS}, FULL-connected brute-force]
# Run:  ENS=... NVDIR=distill_Nv24 SPLIT=1 T0=3 BINSIZE=10 python3 fs_gevp_connected_claude.py
#
# Unified: EVERY matrix element uses the SAME full-connected Wick rule (A,B,C,D,E,G kept; F,H,I,J dropped),
# via a brute-force sum over the 4! pairings of the 4 vertices, keeping only sink<->source bridging perms.
#   op0 sigma^2(tau) -> vertices [tau, tau] ;  op1 O_sigmasigma(tau) -> [tau, tau+delta].
#   C_ab(t,s): sink(a,t) ++ source(b,s), sink={0,1}, source={2,3}.
#   FS = S-part(leg=tau, contact tau(a,a)-1/2 I) + Stilde-part(leg=-tau', contact -1/2(tau'+tau)(a,a)).
# Fixes the earlier inconsistency (O_2sigma sector was E-only). See fs_gevp_connected_impl_plan_claude.md.

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
import diag_effmass_claude as de

SPLIT = int(os.environ.get("SPLIT", "1"))
T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
DTMAX = int(os.environ.get("DTMAX", "26"))

# 24 permutations of 4 vertices: cycles + connected(bridging sink{0,1}<->source{2,3})
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
    conn = any(any(v in (0, 1) for v in c) and any(v in (2, 3) for v in c) for c in cycles)
    PERMS.append((cycles, conn))


def wick_connected(times, Phi, legfns):
    # sum of the CONNECTED (bridging) permutations; legfns summed (FS = [legS, legSt])
    tot = 0.0
    for legfn in legfns:
        for cycles, conn in PERMS:
            if not conn:
                continue
            val = (-1.0) ** len(cycles)
            for cyc in cycles:
                m = Phi[times[cyc[0]]]
                for idx in range(len(cyc)):
                    a = times[cyc[idx]]
                    b = times[cyc[(idx + 1) % len(cyc)]]
                    m = m @ legfn(a, b)
                    if idx + 1 < len(cyc):
                        m = m @ Phi[b]
                val = val * np.trace(m)
            tot += val
    return tot.real


def sink_src(op, tau, d):
    return [tau, tau] if op == 0 else [tau, tau + d]


def per_config(k, w00):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Iv = np.eye(tau.shape[-1])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    def legS(a, b):
        return tau[a, b] - 0.5 * Iv if a == b else tau[a, b]

    def legSt(a, b):
        return -0.5 * (taugw[a, a] + tau[a, a]) if a == b else -taugw[a, b]

    legs = [legS, legSt]
    d = SPLIT
    C = np.full((2, 2, DTMAX), np.nan)
    for dt in range(DTMAX):
        acc = np.zeros((2, 2))
        cnt = 0
        for s in range(twin):
            t = s + dt
            # require all vertex times in window
            if max(t + d, s + d) >= twin:
                continue
            for a in range(2):
                for b in range(2):
                    times = sink_src(a, t, d) + sink_src(b, s, d)
                    acc[a, b] += wick_connected(times, Phi, legs)
            cnt += 1
        if cnt > 0:
            C[:, :, dt] = acc / cnt
    return C


def gevp(Ct, C0):
    C0 = 0.5 * (C0 + C0.T)
    w, U = np.linalg.eigh(C0)
    keep = w > 1e-10 * w.max()
    Uk = U[:, keep] / np.sqrt(w[keep])
    M = Uk.T @ (0.5 * (Ct + Ct.T)) @ Uk
    ev = np.sort(np.linalg.eigvals(M).real)[::-1]
    return ev


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    print("# ENS=%s ncfg=%d  FS connected GEVP {sigma_FS^2, O_sigmasigma}  SPLIT=%d T0=%d BINSIZE=%d"
          % (tag, len(dc.KS), SPLIT, T0, BINSIZE))

    # cross-check on first config: C00 S-part only should equal sum_{CONN} W10*diags_pair (FS S-part)
    de.CONTACT = 0.5
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    Phi0 = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    def legS0(a, b):
        return tau[a, b] - 0.5 * Iv if a == b else tau[a, b]

    CONN_IDX = [0, 1, 2, 3, 4, 6]
    dtc = 4
    brute = sum(wick_connected([dtc + s, dtc + s, s, s], Phi0, [legS0]) for s in range(twin - dtc)) / (twin - dtc)
    dpr = np.mean([(dc.W10 * de.diags_pair(Phi0, tau, s, s + dtc))[CONN_IDX].sum() for s in range(twin - dtc)])
    print("# CROSS-CHECK C00 S-part brute=%.6e  vs  sum_CONN W10*diags_pair=%.6e  ratio=%.6f"
          % (brute, dpr, brute / dpr))

    allC = np.array([per_config(k, w00) for k in dc.KS])       # (ncfg, 2, 2, DTMAX)
    ncfg = allC.shape[0]
    # symmetrize off-diagonal
    allC[:, 0, 1, :] = allC[:, 1, 0, :] = 0.5 * (allC[:, 0, 1, :] + allC[:, 1, 0, :])

    # binsize blocking
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])

    def effmass_from(Cmat):
        tmax = Cmat.shape[-1]
        ev = np.full((tmax, 2), np.nan)
        for dt in range(tmax):
            if np.any(~np.isfinite(Cmat[:, :, dt])):
                continue
            try:
                ev[dt] = gevp(Cmat[:, :, dt], Cmat[:, :, T0])
            except Exception:
                pass
        with np.errstate(all="ignore"):
            em = np.log(ev[:-1] / ev[1:])
        return em

    em_c = effmass_from(blk.mean(0))
    ems = np.array([effmass_from(np.delete(blk, i, 0).mean(0)) for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))

    print("\n#  t |   m0(err)        m1(err)")
    for t in range(T0, min(DTMAX - 1, 22)):
        print("#  %2d | %7.4f(%.4f)  %7.4f(%.4f)" % (t, em_c[t, 0], em_err[t, 0], em_c[t, 1], em_err[t, 1]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    M2PS = 0.644
    ts = np.arange(em_c.shape[0])
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.axhline(M2PS, color="gray", ls="--", lw=1, alpha=0.7)
    ax.text(ts[-1] * 0.6, M2PS + 0.012, r"$2m_{PS}\approx0.644$", fontsize=9, color="gray")
    cols = ["tab:green", "tab:red"]
    mkr = ["o", "s"]
    labs = ["state 0", "state 1"]
    for n in range(2):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n], marker=mkr[n], ms=5, lw=1.1,
                    capsize=2.5, label=labs[n])
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(DTMAX - 1, 22))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("FS connected GEVP {$\\sigma_{FS}^2$, $O_{\\sigma\\sigma}$}, full-connected  T0=%d  %s L1 %d cfg"
                 % (T0, tag, ncfg), fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/fs_gevp_connected_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
