#!/usr/bin/env python3
# fs_gevp_point_perclass_claude.py
#   Diagram-by-diagram (connected classes A,B,C,D,E,G) DIAGONAL correlator for the 3 FS channels
#   {sigma^2_00, O_2m (antipodal), O_1m (coincident split)}, superimposed per class panel.
#   Reuses the validated builder from fs_gevp_point_claude (folded kernels, S+Stilde legs).
#   Run:  ENS=... NVDIR=distill_Nv24 SPLIT=1 NCFG=60 NPROC=4 python3 fs_gevp_point_perclass_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G

SPLIT = G.SPLIT
DTMAX = int(os.environ.get("DTMAX", "18"))
NCFG = int(os.environ.get("NCFG", "60"))
NPROC = int(os.environ.get("NPROC", "4"))
CLASSES = ["A", "B", "C", "D", "E", "G"]


def classify(cycles):
    sizes = sorted(len(c) for c in cycles)
    if sizes == [4]:
        cyc = cycles[0]
        pos = {v: i for i, v in enumerate(cyc)}
        adj = (abs(pos[0] - pos[1]) % 4) in (1, 3)     # sink vertices cyclically adjacent -> A(S_S) ; else B(T_S)
        return "A" if adj else "B"
    if sizes == [1, 3]:
        one = [c for c in cycles if len(c) == 1][0][0]
        return "C" if one in (0, 1) else "D"
    if sizes == [2, 2]:
        return "E"
    return "G"                                          # [1,1,2]


PERM_CLASS = [classify(c) for c in G.PERMS]


CLS_IDX = [CLASSES.index(c) for c in [classify(cy) for cy in G.PERMS]] if False else None


def perclass_one_config(k, KER, OFF, dual):
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(k)
    dualf = dual.astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    nop = len(KER)
    ncl = len(CLASSES)
    cidx = [CLASSES.index(PERM_CLASS[ip]) for ip in range(len(G.PERMS))]
    C = np.full((ncl, nop, nop, DTMAX), np.nan)        # per class, full symmetric matrix
    omax = max(max(o) for o in OFF)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt + omax < twin and s + omax < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        offs = set()
        for a in range(nop):
            for b in range(nop):
                vt = [dt + OFF[a][0], dt + OFF[a][1], OFF[b][0], OFF[b][1]]
                for va in range(4):
                    for vb in range(4):
                        offs.add((vt[va], vt[vb]))
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        bASt = {o: np.array([AblkSt(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        for a in range(nop):
            for b in range(a, nop):
                vt = [dt + OFF[a][0], dt + OFF[a][1], OFF[b][0], OFF[b][1]]
                vspec = G.op_vspec(a, ('i', 'j'), dualf, wY) + G.op_vspec(b, ('k', 'l'), dualf, wY)
                acc = np.zeros(ncl)
                for bA in (bAS, bASt):
                    for ip, cyc in enumerate(G.PERMS):
                        v = G.perm_contrib_folded(cyc, vt, bA, vspec, Pmap)
                        acc[cidx[ip]] += v.real
                C[:, a, b, dt] = acc / len(s0s)
                C[:, b, a, dt] = C[:, a, b, dt]
    return C


_WK = {}


def _init(KER, OFF, dual):
    _WK.update(KER=KER, OFF=OFF, dual=dual)


def _work(k):
    return perclass_one_config(k, _WK["KER"], _WK["OFF"], _WK["dual"])


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    nsite = dual.shape[0]
    Pmap = G.antipodal_map()
    d = SPLIT
    Kw = np.outer(dual * dc.Y00, dual * dc.Y00)
    K2m = np.zeros((nsite, nsite))
    for i in range(nsite):
        K2m[i, Pmap[i]] = dual[i]
    K1m = np.diag(dual.astype(float))
    KER = [Kw, K2m, K1m]
    OFF = [(0, 0), (0, 0), (0, d)]
    ks = dc.KS[:NCFG]
    print("# ENS=%s  per-class diagonal correlators  ncfg=%d NPROC=%d" % (tag, len(ks), NPROC))
    if NPROC > 1:
        import multiprocessing as mp
        with mp.Pool(NPROC, initializer=_init, initargs=(KER, OFF, dual)) as pool:
            res = pool.map(_work, ks)
        allC = np.array(res)
    else:
        allC = np.array([perclass_one_config(k, KER, OFF, dual) for k in ks])   # (ncfg,6,3,DTMAX)
    n = allC.shape[0]
    samp = np.array([np.delete(allC, i, 0).mean(0) for i in range(n)])
    cen = samp.mean(0)
    err = np.sqrt((n - 1) * np.mean((samp - cen) ** 2, 0))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    os.makedirs("figs", exist_ok=True)

    def make(series, fname, title, logabs):
        # series = list of (label, col, mk, ai, bi)  -> cen[ci, ai, bi, :]
        fig, axs = plt.subplots(2, 3, figsize=(14, 8))
        axs = axs.ravel()
        for ci in range(6):
            ax = axs[ci]
            for a, (lab, col, mk, ai, bi) in enumerate(series):
                y = cen[ci, ai, bi, dts]
                yy = np.abs(y) if logabs else y
                ax.errorbar(dts + 0.04 * a, yy, yerr=err[ci, ai, bi, dts], color=col, marker=mk,
                            ms=4, lw=1, capsize=2, label=lab)
            if logabs:
                ax.set_yscale("log")
                ax.grid(alpha=0.3, which="both")
            else:
                ax.axhline(0.0, color="gray", lw=0.8, alpha=0.6)
                ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            ax.set_title("diagram %s" % CLASSES[ci], fontsize=11)
            ax.set_xlabel(r"$dt$", fontsize=9)
            if ci == 0:
                ax.legend(fontsize=8)
        fig.suptitle("%s  %s L1 %d cfg (FS, connected)" % (title, tag, n), fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig("figs/" + fname, dpi=130)
        plt.close(fig)
        print("# -> figs/%s" % fname)

    diag = [("sigma^2_00", "tab:red", "o", 0, 0), ("O_2m (antipodal)", "tab:blue", "s", 1, 1),
            ("O_1m (split)", "tab:green", "^", 2, 2)]
    cross = [("sigma^2_00 x O_2m", "tab:purple", "o", 0, 1), ("sigma^2_00 x O_1m", "tab:orange", "s", 0, 2),
             ("O_2m x O_1m", "tab:brown", "^", 1, 2)]
    make(diag, "fs_gevp_point_perclass_%s_claude.png" % tag, "Per-diagram DIAGONAL correlator, 3 channels superimposed", False)
    make(diag, "fs_gevp_point_perclass_log_%s_claude.png" % tag, "Per-diagram DIAGONAL |correlator| (LOG), 3 channels superimposed", True)
    make(cross, "fs_gevp_point_perclass_cross_%s_claude.png" % tag, "Per-diagram CROSS correlator, 3 crosses superimposed", False)
    make(cross, "fs_gevp_point_perclass_cross_log_%s_claude.png" % tag, "Per-diagram CROSS |correlator| (LOG), 3 crosses superimposed", True)


if __name__ == "__main__":
    main()
