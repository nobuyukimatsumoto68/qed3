#!/usr/bin/env python3
# gevp_twosigma_interacting_claude.py  [INTERACTING {1, sigma^2, O_2sigma, O_A} GEVP; config-avg + jackknife]
# Run:  ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000 \
#         NVDIR=distill_Nv24 SPLIT=1 T0=2 python3 gevp_twosigma_interacting_claude.py
#
# O_22 is DROPPED for the interacting theory: the gauge broadens the S^2 shells past their spacing
# (checked: no 4/8/12 clustering at gsq 0.5/1.0/1.5), so the lambda=2 spectral projector Q2 is ill-defined.
# Basis = {1, sigma^2, O_2sigma = sigma_00(t)sigma_00(t+SPLIT), O_A = psibar Phi tilde_tau psi}.
# Single-sigma meson prop M(t,s) = -Tr[Phi(t) tau(t,s) Phi(s) tau(s,t)].
# Full correlators config-averaged; identity carries vacuum; jackknife over configs for errors.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = int(os.environ.get("T0", "2"))
SPLIT = int(os.environ.get("SPLIT", "1"))
NCFG = int(os.environ.get("NCFG", "0"))            # 0 -> all
CONN = int(os.environ.get("CONN", "1"))            # 1 -> connected {sigma^2(ABCDEG), O_A, O_2sigma}, no identity
NOOA = int(os.environ.get("NOOA", "0"))            # 1 -> drop O_A: connected 2x2 {sigma^2, O_2sigma}
BINSCAN = int(os.environ.get("BINSCAN", "0"))      # 1 -> jackknife binsize (blocking) scan for autocorrelation
BINSIZE = int(os.environ.get("BINSIZE", "1"))      # jackknife block size for the main effmass errors (autocorr)
HANKEL = int(os.environ.get("HANKEL", "0"))        # 1 -> block-Hankel proliferation + rebase onto NKEEP states
NSH = int(os.environ.get("NSH", "3"))              # Hankel blocks (shifts 0..NSH-1); Dt=1,2 -> NSH=3, SHIFT=1
SHIFT = int(os.environ.get("SHIFT", "1"))
REBT = int(os.environ.get("REBT", "5"))            # rebase time
NKEEP = int(os.environ.get("NKEEP", "2"))          # rebase onto this many leading states
REBSCAN = int(os.environ.get("REBSCAN", "0"))      # 1 -> scan the rebase time; read masses at RREAD
RREAD = int(os.environ.get("RREAD", "10"))         # readout time for the rebase-point scan
REBT0 = int(os.environ.get("REBT0", "-1"))         # metric time for the rebase eigenvectors; -1 -> T0
CONN_IDX = [0, 1, 2, 3, 4, 6]                       # A,B,C,D,E,G  (drop F,H,I,J vacuum diagrams)


def per_config(k, w00, d):
    de.CONTACT = 0.5
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Nv = V.shape[1]
    Iv = np.eye(Nv)
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]
    Msig = np.zeros((twin, twin))
    for t in range(twin):
        for s in range(twin):
            Msig[t, s] = (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
    hi = twin - 1
    diags = np.full((twin, 10), np.nan)
    CAA = np.full(twin, np.nan)
    Css = np.full(twin, np.nan)
    C2s = np.full(twin, np.nan)
    C2A = np.full(twin, np.nan)
    for dt in range(twin):
        ns = twin - dt
        dacc = np.zeros(10)
        for s in range(ns):
            dacc += (dc.W10 * de.diags_pair(Phi, tau, s, s + dt)).real
        diags[dt] = dacc / ns
        aAA = a2A = ass = a2s = 0.0
        cnt = 0
        for s in range(ns):
            t = s + dt
            if t + d > hi or s + d > hi:
                continue
            aAA += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
            a2A += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PA[t] @ tau[t, s])).real
            ass += Msig[t, s] * Msig[t + d, s + d] + Msig[t, s + d] * Msig[t + d, s]
            a2s += 2.0 * Msig[t, s] * Msig[t + d, s]
            cnt += 1
        if cnt > 0:
            CAA[dt] = aAA / cnt
            Css[dt] = ass / cnt
            C2s[dt] = a2s / cnt
            C2A[dt] = a2A / cnt
    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])
    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(twin)])
    o2s = np.mean([Msig[a, a + d] for a in range(twin - d)]) + \
        np.mean([(-np.trace(Phi[a] @ tt[a, a])).real for a in range(twin)]) ** 2
    return diags, CAA, Css, C2s, C2A, o2, oA, o2s, twin


def build_matrix_conn(D, twin):
    # connected 3x3 {0:sigma^2 (A,B,C,D,E,G only), 1:O_A, 2:O_2sigma}; NO identity (vacuum-free by
    # construction: connected diagrams decay; F,H,I,J dropped).  Energies are normalization-invariant.
    diags = D["diags"]
    if NOOA:
        # connected 2x2 {0:sigma^2 (ABCDEG), 1:O_2sigma}
        Cts = np.full((twin, 2, 2), np.nan)
        for dt in range(twin):
            M = np.zeros((2, 2))
            M[0, 0] = 2.0 * diags[dt][CONN_IDX].sum()
            M[1, 1] = D["Css"][dt]
            M[0, 1] = M[1, 0] = D["C2s"][dt]
            Cts[dt] = M
        return Cts
    Cts = np.full((twin, 3, 3), np.nan)
    for dt in range(twin):
        c22 = 2.0 * diags[dt][CONN_IDX].sum()         # connected <sigma^2 sigma^2> = A+B+C+D+E+G
        M = np.zeros((3, 3))
        M[0, 0] = c22
        M[1, 1] = D["CAA"][dt]                          # <O_A O_A> connected loop
        M[2, 2] = D["Css"][dt]                          # <O_2s O_2s> connected
        M[0, 1] = M[1, 0] = D["C2A"][dt]               # <sigma^2 O_A> = 2*Tri (connected triangle)
        M[0, 2] = M[2, 0] = D["C2s"][dt]               # <sigma^2 O_2s> connected
        M[1, 2] = M[2, 1] = 0.0                         # <O_A O_2s> ~ 0 (2f-4f)
        Cts[dt] = M
    return Cts


def build_matrix(D, twin):
    # D: dict of ensemble-averaged pieces; returns Cts (twin,4,4).  order {0:1,1:sigma^2,2:O_2s,3:O_A}
    Cts = np.full((twin, 4, 4), np.nan)
    diags = D["diags"]
    for dt in range(twin):
        M = np.zeros((4, 4))
        M[0, 0] = 1.0
        M[0, 1] = M[1, 0] = D["o2"]
        M[0, 2] = M[2, 0] = D["o2s"]
        M[0, 3] = M[3, 0] = D["oA"]
        M[1, 1] = diags[dt].sum()
        M[1, 2] = M[2, 1] = D["C2s"][dt] + D["o2_o2s"]
        M[1, 3] = M[3, 1] = D["C2A"][dt] + D["o2_oA"]
        M[2, 2] = D["Css"][dt] + D["o2s_o2s"]
        M[2, 3] = M[3, 2] = 0.0 + D["oA_o2s"]
        M[3, 3] = D["CAA"][dt] + D["oA_oA"]
        Cts[dt] = M
    return Cts


def gevp(Cts, t0, tol=1e-10):
    twin = Cts.shape[0]
    C0 = 0.5 * (Cts[t0] + Cts[t0].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > tol * wv.max()
    Uk = Uv[:, keep] / np.sqrt(wv[keep])
    nlev = int(keep.sum())
    lam = np.full((twin, nlev), np.nan)
    for t in range(twin):
        Mt = Uk.T @ (0.5 * (Cts[t] + Cts[t].T)) @ Uk
        try:
            lam[t] = np.sort(np.linalg.eigvals(Mt).real)[::-1]
        except Exception:
            pass
    return lam, nlev


def zmatrix(Cts, t, t0, tol=1e-10):
    # Z_n^a = <0|O_a|n> = (C(t0) v_n)_a exp(E_n t0/2), v_n normalized v_n^T C(t0) v_n = 1.
    # states sorted by ascending E_n (lightest first).  Returns Z (nops x nlev), E (nlev).
    C0 = 0.5 * (Cts[t0] + Cts[t0].T)
    Ct = 0.5 * (Cts[t] + Cts[t].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > tol * wv.max()
    Uk = Uv[:, keep] / np.sqrt(wv[keep])
    val, vec = np.linalg.eigh(Uk.T @ Ct @ Uk)             # symmetric reduced problem
    order = np.argsort(val)[::-1]                          # descending lambda = ascending E
    val = val[order]
    vec = vec[:, order]
    E = -np.log(np.abs(val)) / (t - t0)
    Vop = Uk @ vec                                        # (nops, nlev), v_n^T C0 v_n = 1
    Z = np.zeros((C0.shape[0], vec.shape[1]))
    for n in range(vec.shape[1]):
        Z[:, n] = (C0 @ Vop[:, n]) * np.exp(E[n] * t0 / 2)
    return Z, E


def hankel(Cts, nsh, shift):
    # block-Hankel: Chat(t)_{(a,i),(b,j)} = Cts[t + (a+b)*shift] ; needs t + 2*(nsh-1)*shift <= twin-1
    twin, N, _ = Cts.shape
    tmax = twin - 2 * (nsh - 1) * shift
    Big = np.full((tmax, nsh * N, nsh * N), np.nan)
    for t in range(tmax):
        for a in range(nsh):
            for b in range(nsh):
                Big[t, a * N:(a + 1) * N, b * N:(b + 1) * N] = Cts[t + (a + b) * shift]
    return Big


def rebase_vectors(Big, t_reb, t0, nkeep, tol=1e-10):
    # leading nkeep GEVP eigenvectors of Big at t_reb (vs t0); columns of Vop (Nbig x nkeep)
    C0 = 0.5 * (Big[t0] + Big[t0].T)
    Ct = 0.5 * (Big[t_reb] + Big[t_reb].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > tol * wv.max()
    Uk = Uv[:, keep] / np.sqrt(wv[keep])
    val, vec = np.linalg.eigh(Uk.T @ Ct @ Uk)
    order = np.argsort(val)[::-1]
    return Uk @ vec[:, order[:nkeep]]


def rebased_effmass(Big, Vop, t0):
    # project Big onto Vop -> nkeep x nkeep, GEVP, effmass
    Cr = np.einsum("ai,tab,bj->tij", Vop, Big, Vop)
    lam, _ = gevp(Cr, t0)
    with np.errstate(all="ignore"):
        return np.log(lam[:-1] / lam[1:])


def ens_avg(store, idx):
    # average per-config arrays/scalars over configs in idx; form one-point products
    D = {}
    for key in ["diags", "CAA", "Css", "C2s", "C2A"]:
        D[key] = np.mean([store[key][i] for i in idx], axis=0)
    o2 = np.array([store["o2"][i] for i in idx])
    oA = np.array([store["oA"][i] for i in idx])
    o2s = np.array([store["o2s"][i] for i in idx])
    D["o2"], D["oA"], D["o2s"] = o2.mean(), oA.mean(), o2s.mean()
    D["o2_o2s"] = np.mean(o2 * o2s)
    D["o2_oA"] = np.mean(o2 * oA)
    D["o2s_o2s"] = np.mean(o2s * o2s)
    D["oA_o2s"] = np.mean(oA * o2s)
    D["oA_oA"] = np.mean(oA * oA)
    return D


def main():
    tag = dc.ENS.split("nu0")[0]
    d = SPLIT
    KS = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    ncfg = len(KS)
    basis = "connected {sigma^2(ABCDEG),O_A,O_2sigma} NO identity" if CONN else "{1,sigma^2,O_2sigma,O_A}"
    print("# ENS=%s  ncfg=%d  %s  split d=%d  T0=%d" % (tag, ncfg, basis, d, T0))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    store = {k: [] for k in ["diags", "CAA", "Css", "C2s", "C2A", "o2", "oA", "o2s"]}
    twin = None
    for j, k in enumerate(KS):
        diags, CAA, Css, C2s, C2A, o2, oA, o2s, tw = per_config(k, w00, d)
        twin = tw
        store["diags"].append(diags)
        store["CAA"].append(CAA)
        store["Css"].append(Css)
        store["C2s"].append(C2s)
        store["C2A"].append(C2A)
        store["o2"].append(o2)
        store["oA"].append(oA)
        store["o2s"].append(o2s)
        if (j + 1) % 100 == 0:
            print("#   ... %d/%d configs" % (j + 1, ncfg))

    bmf = build_matrix_conn if CONN else build_matrix
    allidx = list(range(ncfg))
    Dcen = ens_avg(store, allidx)
    Ccen = bmf(Dcen, twin)

    if BINSCAN:
        keys = ["diags", "CAA", "Css", "C2s", "C2A", "o2", "oA", "o2s"]
        tsel = [int(x) for x in os.environ.get("BSEL", "10,15").split(",")]
        binsizes = [1, 2, 4, 5, 8, 10, 16, 20, 25, 40]
        print("# BINSIZE SCAN  observables: effmass at dt=%s, all levels" % tsel)
        res = {}                                          # (dt,lev) -> list of (b, err)
        for b in binsizes:
            nbin = ncfg // b
            if nbin < 6:
                continue
            bstore = {key: [np.mean([store[key][i * b + r] for r in range(b)], axis=0) for i in range(nbin)]
                      for key in keys}
            ems = []
            for i in range(nbin):
                idx = [j for j in range(nbin) if j != i]
                Ci = bmf(ens_avg(bstore, idx), twin)
                lam_i, _ = gevp(Ci, T0)
                with np.errstate(all="ignore"):
                    ems.append(np.log(lam_i[:-1] / lam_i[1:]))
            ems = np.array(ems)                            # (nbin, twin-1, nlev)
            emm = ems.mean(0)
            eme = np.sqrt((nbin - 1) * np.mean((ems - emm) ** 2, 0))
            for dt in tsel:
                for lev in range(ems.shape[2]):
                    res.setdefault((dt, lev), []).append((b, eme[dt, lev]))
            line = "  ".join("dt%d/l%d=%.4f" % (dt, lev, eme[dt, lev]) for dt in tsel for lev in range(ems.shape[2]))
            print("#  binsize %2d (nbin=%3d):  %s" % (b, nbin, line))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8.5, 5.6))
        for (dt, lev), pts in res.items():
            bs = [p[0] for p in pts]
            er = [p[1] for p in pts]
            ax.plot(bs, er, marker="o", ms=5, lw=1.2, label=r"$dt{=}%d$, level %d" % (dt, lev))
        ax.set_xlabel("bin size (configs/block)")
        ax.set_ylabel(r"jackknife error on $a_t m_\mathrm{eff}$")
        ax.set_title("Binsize (blocking) scan  %s  connected %s GEVP  %d cfg"
                     % (tag, "2op" if NOOA else "3op", ncfg))
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        os.makedirs("figs", exist_ok=True)
        outb = "figs/gevp_twosigma_interacting_binscan_%s_%s_claude.png" % (tag.replace(".", "p"), "2op" if NOOA else "3op")
        fig.savefig(outb, dpi=130)
        plt.close(fig)
        print("\n# -> %s" % outb)
        return

    if HANKEL and REBSCAN:
        keys = ["diags", "CAA", "Css", "C2s", "C2A", "o2", "oA", "o2s"]
        Big = hankel(Ccen, NSH, SHIFT)
        tmax = Big.shape[0]
        rt0 = REBT0 if REBT0 >= 0 else T0                     # metric for the rebase eigenvectors
        rebts = [r for r in range(rt0 + 1, RREAD) if r < tmax]
        if BINSIZE > 1:
            nbin = ncfg // BINSIZE
            jkstore = {key: [np.mean([store[key][i * BINSIZE + rr] for rr in range(BINSIZE)], axis=0)
                             for i in range(nbin)] for key in keys}
            njk = nbin
        else:
            jkstore = store
            njk = ncfg
        Bjk = [hankel(bmf(ens_avg(jkstore, [j for j in range(njk) if j != i]), twin), NSH, SHIFT)
               for i in range(njk)]
        print("# REBASE-POINT SCAN  NSH=%d shift=%d, read masses at t=%d; BINSIZE=%d" % (NSH, SHIFT, RREAD, BINSIZE))
        print("\n# rebase_t | " + "   ".join("m%d(err)" % n for n in range(NKEEP)))
        res = {n: [] for n in range(NKEEP)}
        for rt in rebts:
            Vop = rebase_vectors(Big, rt, rt0, NKEEP)
            em_c = rebased_effmass(Big, Vop, T0)[RREAD]
            emj = np.array([rebased_effmass(Bjk[i], Vop, T0)[RREAD] for i in range(njk)])
            eme = np.sqrt((njk - 1) * np.mean((emj - emj.mean(0)) ** 2, 0))
            print("#    %2d     | %s" % (rt, "   ".join("%6.4f(%.4f)" % (em_c[n], eme[n]) for n in range(NKEEP))))
            for n in range(NKEEP):
                res[n].append((rt, em_c[n], eme[n]))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8.5, 5.6))
        cols = ["tab:green", "tab:red", "tab:blue"]
        for n in range(NKEEP):
            r = [p[0] for p in res[n]]
            m = [p[1] for p in res[n]]
            e = [p[2] for p in res[n]]
            ax.errorbar(r, m, yerr=e, color=cols[n % 3], marker="o", ms=5, lw=1.2, capsize=3, label="state %d" % n)
        ax.set_xlabel("rebase time $t_{\\rm reb}$")
        ax.set_ylabel(r"$a_t m_\mathrm{eff}$ (read at $t=%d$)" % RREAD)
        ax.set_title("Rebase-point scan  Dt=%d Hankel  %s connected %s  %d cfg b=%d"
                     % (SHIFT, tag, "2op" if NOOA else "3op", ncfg, BINSIZE))
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        os.makedirs("figs", exist_ok=True)
        outs = "figs/gevp_twosigma_interacting_rebscan_%s_%s_s%d_claude.png" % (tag.replace(".", "p"), "2op" if NOOA else "3op", SHIFT)
        fig.savefig(outs, dpi=130)
        plt.close(fig)
        print("\n# -> %s" % outs)
        return

    if HANKEL:
        keys = ["diags", "CAA", "Css", "C2s", "C2A", "o2", "oA", "o2s"]
        Big = hankel(Ccen, NSH, SHIFT)
        Vop = rebase_vectors(Big, REBT, T0, NKEEP)     # central rebase basis at t=REBT
        em_c = rebased_effmass(Big, Vop, T0)           # (tmax-1, NKEEP)
        if BINSIZE > 1:
            nbin = ncfg // BINSIZE
            jkstore = {key: [np.mean([store[key][i * BINSIZE + r] for r in range(BINSIZE)], axis=0)
                             for i in range(nbin)] for key in keys}
            njk = nbin
        else:
            jkstore = store
            njk = ncfg
        ems = []
        for i in range(njk):
            idx = [j for j in range(njk) if j != i]
            Bi = hankel(bmf(ens_avg(jkstore, idx), twin), NSH, SHIFT)
            ems.append(rebased_effmass(Bi, Vop, T0))   # FIXED central Vop
        ems = np.array(ems)
        em_err = np.sqrt((njk - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
        tmax = Big.shape[0]
        print("# BLOCK-HANKEL (NSH=%d shift=%d -> %dx%d) + rebase onto %d states at t=%d; BINSIZE=%d"
              % (NSH, SHIFT, NSH * Ccen.shape[1], NSH * Ccen.shape[1], NKEEP, REBT, BINSIZE))
        print("\n#  t | " + "   ".join("m%d(err)" % n for n in range(NKEEP)))
        for t in range(T0 + 1, min(tmax - 1, 25)):
            cells = "   ".join("%6.4f(%.4f)" % (em_c[t, n], em_err[t, n]) if np.isfinite(em_c[t, n]) else "   ---   "
                               for n in range(NKEEP))
            print("#  %2d | %s" % (t, cells))
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ts = np.arange(tmax - 1)
        fig, ax = plt.subplots(figsize=(8.5, 5.6))
        cols = ["tab:green", "tab:red", "tab:blue"]
        for n in range(NKEEP):
            g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n])
            ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n % 3], marker="o", ms=4, lw=1,
                        capsize=2, label="state %d" % n)
        ax.axvline(REBT, color="gray", ls=":", lw=1, alpha=0.7)
        ax.text(REBT + 0.1, 1.15, "rebase t=%d" % REBT, fontsize=8, alpha=0.7)
        ax.set_ylim(0.0, 1.3)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
        ax.set_title("Block-Hankel(Dt=1,2)+rebase %d@t=%d  %s connected %s  %d cfg b=%d"
                     % (NKEEP, REBT, tag, "2op" if NOOA else "3op", ncfg, BINSIZE))
        ax.legend(fontsize=9)
        fig.tight_layout()
        os.makedirs("figs", exist_ok=True)
        outh = "figs/gevp_twosigma_interacting_hankelreb_%s_%s_s%d_claude.png" % (tag.replace(".", "p"), "2op" if NOOA else "3op", SHIFT)
        fig.savefig(outh, dpi=130)
        plt.close(fig)
        print("\n# -> %s" % outh)
        return

    good = np.array([np.all(np.isfinite(Ccen[t])) for t in range(twin)])
    tvec = np.where(good)[0]
    t0i = int(np.argmin(np.abs(tvec - T0)))
    lam_c, nlev = gevp(Ccen[tvec], t0i)
    with np.errstate(all="ignore"):
        em_c = np.log(lam_c[:-1] / lam_c[1:])

    # jackknife over blocks of BINSIZE configs (accounts for autocorrelation)
    if BINSIZE > 1:
        keys = ["diags", "CAA", "Css", "C2s", "C2A", "o2", "oA", "o2s"]
        nbin = ncfg // BINSIZE
        jkstore = {key: [np.mean([store[key][i * BINSIZE + r] for r in range(BINSIZE)], axis=0)
                         for i in range(nbin)] for key in keys}
        njk = nbin
    else:
        jkstore = store
        njk = ncfg
    ems = []
    for i in range(njk):
        idx = [j for j in range(njk) if j != i]
        Ci = bmf(ens_avg(jkstore, idx), twin)
        lam_i, _ = gevp(Ci[tvec], t0i)
        with np.errstate(all="ignore"):
            ems.append(np.log(lam_i[:-1] / lam_i[1:]))
    ems = np.array(ems)                                     # (njk, ntw-1, nlev)
    emmean = ems.mean(0)
    em_err = np.sqrt((njk - 1) * np.mean((ems - emmean) ** 2, axis=0))
    print("# BINSIZE=%d for jackknife errors (nbin=%d)" % (BINSIZE, njk))

    print("# levels kept = %d" % nlev)
    print("\n#  t | " + "   ".join("m%d(err)" % i for i in range(nlev)))
    for a in range(len(tvec) - 2):
        t = tvec[a]
        if t <= T0:
            continue
        cells = []
        for klev in range(nlev):
            ok = np.isfinite(em_c[a, klev]) and (lam_c[a, klev] * lam_c[a + 1, klev] > 0)
            cells.append("%6.3f(%5.3f)" % (em_c[a, klev], em_err[a, klev]) if ok else "    ---      ")
        print("#  %2d | %s" % (t, "  ".join(cells)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = tvec[:-1]
    fig, ax = plt.subplots(figsize=(9, 6))
    cols = ["tab:gray", "tab:green", "tab:red", "tab:orange"]
    for klev in range(nlev):
        g = np.isfinite(em_c[:, klev]) & (lam_c[:-1, klev] * lam_c[1:, klev] > 0)
        ax.errorbar(ts[g], em_c[g, klev], yerr=em_err[g, klev], color=cols[klev % 4],
                    marker="o", ms=3, lw=1, capsize=2, label="level %d" % klev)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    bl = r"connected $\{\sigma^2,O_A,O_{2\sigma}\}$" if CONN else r"$\{1,\sigma^2,O_{2\sigma},O_A\}$"
    ax.set_title(r"%s  interacting %s (ncfg=%d, jk)" % (tag, bl, ncfg))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/gevp_twosigma_interacting_%s_d%d%s_claude.png" % (tag.replace(".", "p"), d, (("_conn2" if NOOA else "_conn") if CONN else "") + ("_b%d" % BINSIZE if BINSIZE > 1 else ""))
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)

    # ---- Z (overlap) matrix heatmaps at t = 5,10,15,20 ----
    if CONN and NOOA:
        oplabels = [r"$\sigma^2$", r"$O_{2\sigma}$"]
    elif CONN:
        oplabels = [r"$\sigma^2$", r"$O_A$", r"$O_{2\sigma}$"]
    else:
        oplabels = [r"$\mathbb{1}$", r"$\sigma^2$", r"$O_{2\sigma}$", r"$O_A$"]
    tshow = [5, 10, 15, 20]
    figz, azs = plt.subplots(1, 4, figsize=(16, 4.4))
    for p, tz in enumerate(tshow):
        Z, E = zmatrix(Ccen, tz, T0)               # Ccen indexed by dt directly; T0 is a dt
        nop, nst = Z.shape
        # normalize each state (column) to unit max |Z| so operator content is comparable
        Zn = Z / (np.abs(Z).max(0, keepdims=True) + 1e-300)
        ax = azs[p]
        im = ax.imshow(Zn, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(nst))
        ax.set_xticklabels([r"$n{=}%d$" "\n" r"$E{=}%.2f$" % (n, E[n]) for n in range(nst)], fontsize=8)
        ax.set_yticks(range(nop))
        ax.set_yticklabels(oplabels[:nop], fontsize=11)
        ax.set_title(r"$t=%d$" % tz, fontsize=11)
        for a in range(nop):
            for n in range(nst):
                ax.text(n, a, "%+.2f" % Zn[a, n], ha="center", va="center", fontsize=8,
                        color="black" if abs(Zn[a, n]) < 0.6 else "white")
    figz.suptitle(r"Overlap $Z_n^a$ (per-state normalized)  %s connected GEVP  T0=%d  %d cfg"
                  % (tag, T0, ncfg), fontsize=12)
    figz.tight_layout(rect=[0, 0, 1, 0.94])
    outz = "figs/gevp_twosigma_interacting_Zheat_%s_d%d%s_claude.png" % (tag.replace(".", "p"), d, "_2op" if NOOA else "")
    figz.savefig(outz, dpi=130)
    plt.close(figz)
    print("# -> %s" % outz)


if __name__ == "__main__":
    main()
