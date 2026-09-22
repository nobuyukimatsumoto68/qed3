#!/usr/bin/env python3
# fs_channels_v2_claude.py  [CORRECTED FS channels: 3 diagonal + 3 cross correlators, per-loop + improved prop]
#   Operators {0:sigma^2_00, 1:O_2m (antipodal), 2:O_1m (coincident split)} -- the SAME three as fs_gevp_point.
#   CORRECTED FS contraction (per fs_diag_corr_v2): per closed loop sum S + Stilde; by GW the Stilde loop
#   collapses to the forward improved loop, so each loop -> 2x (improved), and every TADPOLE loop vanishes
#   (kernel-summed condensate = 0 with the 1/2-improved propagator).  Implementation: forward improved leg
#   only (tau - 1/2 I on equal time) x a per-permutation factor 2^{#cycles}.  (Old fs_gevp_point used the
#   wrong furnished -taugw -> live tadpole / large G.)
#   Vacuum handling: per-config t-sum subtraction + per-channel plateau (double subtraction); LOG|C| plots.
#   Config-parallel (NPROC).  Reuses the validated point-block machinery from fs_gevp_point_claude.
#
# Run (parallel handoff): OMP_NUM_THREADS=1 NPROC=16 SPLIT=1 DTMAX=30 python3 fs_channels_v2_claude.py

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

SPLIT = int(os.environ.get("SPLIT", "1"))
DTMAX = int(os.environ.get("DTMAX", "30"))
NPROC = int(os.environ.get("NPROC", "1"))
NCFG = int(os.environ.get("NCFG", "0"))
VALIDATE = int(os.environ.get("VALIDATE", "0"))
TSUM_LO = int(os.environ.get("TSUM_LO", "1"))
NCYC = np.array([2.0 ** len(c) for c in G.PERMS])       # per-loop FS factor 2^{#cycles}
OPLAB = ["sigma^2_00", "O_2m", "O_1m"]


def matrix_one_config(k, KER, OFF, dual):
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(k)
    dualf = dual.astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    nop = len(KER)
    C = np.full((nop, nop, DTMAX), np.nan)
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
        for a in range(nop):
            for b in range(nop):
                vt = [dt + OFF[a][0], dt + OFF[a][1], OFF[b][0], OFF[b][1]]
                vspec = G.op_vspec(a, ('i', 'j'), dualf, wY) + G.op_vspec(b, ('k', 'l'), dualf, wY)
                v = 0.0
                for ip, cyc in enumerate(G.PERMS):
                    v += NCYC[ip] * G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap)
                C[a, b, dt] = v.real / len(s0s)
    return C


def kernels(dual):
    nsite = dual.shape[0]
    Pmap = G.antipodal_map()
    Kw = np.outer(dual * dc.Y00, dual * dc.Y00)
    K2m = np.zeros((nsite, nsite))
    for i in range(nsite):
        K2m[i, Pmap[i]] = dual[i]
    K1m = np.diag(dual.astype(float))
    return [Kw, K2m, K1m], [(0, 0), (0, 0), (0, SPLIT)]


_WK = {}


def _init(KER, OFF, dual):
    _WK.update(KER=KER, OFF=OFF, dual=dual)


def _work(k):
    return matrix_one_config(k, _WK["KER"], _WK["OFF"], _WK["dual"])


def double_sub(allC_ab, plo, phi):
    # allC_ab: (ncfg, twin) one channel.  per-config t-sum sub -> jackknife + per-channel plateau.
    ts = allC_ab - allC_ab[:, TSUM_LO:].mean(1, keepdims=True)
    n = ts.shape[0]
    samp = np.array([np.delete(ts, k, 0).mean(0) for k in range(n)])
    samp = samp - samp[:, plo:phi].mean(1, keepdims=True)
    return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    KER, OFF = kernels(dual)

    if VALIDATE:
        # tadpole must vanish: C[0,0] (sigma^2_00) connected should be finite & smooth; print a few dt
        C = matrix_one_config(dc.KS[0], KER, OFF, dual)
        print("# VALIDATE corrected FS point channels (1 config)")
        for dt in range(1, 8):
            print("#  dt=%d  C00=% .3e C11=% .3e C22=% .3e C01=% .3e C02=% .3e C12=% .3e"
                  % (dt, C[0, 0, dt], C[1, 1, dt], C[2, 2, dt], C[0, 1, dt], C[0, 2, dt], C[1, 2, dt]))
        return

    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    CACHEDIR = "fs_channels_v2_cache_claude"
    os.makedirs(CACHEDIR, exist_ok=True)
    cache = "%s/fs_channels_v2_%s_%dcfg_d%d_claude.npy" % (CACHEDIR, tag.replace(".", "p"), len(ks), SPLIT)
    if os.path.exists(cache):
        allC = np.load(cache)
        print("# loaded cache <- %s (%d cfg)" % (cache, allC.shape[0]))
    else:
        if NPROC > 1:
            import multiprocessing as mp
            print("# computing %d configs, %d workers ..." % (len(ks), NPROC))
            with mp.Pool(NPROC, initializer=_init, initargs=(KER, OFF, dual)) as pool:
                allC = np.array(pool.map(_work, ks))
        else:
            allC = np.array([matrix_one_config(k, KER, OFF, dual) for k in ks])
        np.save(cache, allC)
        print("# cached -> %s" % cache)
    ncfg = allC.shape[0]
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))          # symmetrize
    twin = allC.shape[-1]
    plo, phi = twin - 8, twin

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    channels = [(0, 0, "sigma^2_00 (diag)"), (1, 1, "O_2m (diag)"), (2, 2, "O_1m (diag)"),
                (0, 1, "sigma^2_00 x O_2m"), (0, 2, "sigma^2_00 x O_1m"), (1, 2, "O_2m x O_1m")]
    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    axs = axs.ravel()
    for p, (a, b, lab) in enumerate(channels):
        cm, ee = double_sub(allC[:, a, b, :], plo, phi)
        axs[p].errorbar(dts, np.abs(cm[dts]), yerr=ee[dts], color="tab:purple", marker="o", ms=4, lw=1, capsize=2)
        axs[p].set_yscale("log")
        axs[p].set_title(lab, fontsize=10)
        axs[p].set_xlabel(r"$dt$", fontsize=9)
        axs[p].grid(alpha=0.3, which="both")
    fig.suptitle("CORRECTED FS channels |correlator| (LOG, double-subtracted)  %s L1 %d cfg" % (tag, ncfg), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = "figs/fs_channels_v2_log_%s_claude.png" % tag
    os.makedirs("figs", exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("# -> %s" % out)


if __name__ == "__main__":
    main()
