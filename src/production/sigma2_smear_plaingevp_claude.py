#!/usr/bin/env python3
# sigma2_smear_plaingevp_claude.py  [chunk 2b: PLAIN pruned fixed-t0 GEVP on the N_v-smearing basis]
#   The smearing tower IS the variational expansion, so use a PLAIN fixed-t0 GEVP (metric pruned to effective
#   rank), NOT Hankel (which double-expands -> over-parametrized garbage).  Reads the 12-op cache, symmetrizes
#   at ensemble level, selects an op subset (OPS), fixed-t0 GEVP with inv_sqrt pruning (rtol), effmass =
#   log(lam(t)/lam(t+1)), jackknife over bins.  Compare state0/state1 to the 3-op baseline 0.46/0.62.
#   Note: mixing the 3 DIFFERENT point-ops at truncated N_v is indefinite (spatial structure washes out);
#   well-conditioned sub-bases are per-operator smearing towers (e.g. sigma^2_00 x {4,8,16,24}).
#   Run: OPS=0,1,2,3 T0=3 RTOL=1e-3 NSTATE=2 BINSIZE=10 python3 sigma2_smear_plaingevp_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import glob
import numpy as np
import distill_contract_claude as dc

T0 = int(os.environ.get("T0", "3"))
RTOL = float(os.environ.get("RTOL", "1e-3"))
NSTATE = int(os.environ.get("NSTATE", "2"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
SPLIT = int(os.environ.get("SPLIT", "1"))
SMEARS = [int(x) for x in os.environ.get("SMEARS", "4,8,16,24").split(",")]
NSM = len(SMEARS)
OPLAB = ["s2_00", "O_2m", "O_1m"]
OPS = os.environ.get("OPS", "0,1,2,3,4,5,6,7,8,9,10,11")
DTMAXP = int(os.environ.get("DTMAXP", "18"))


def gevp(Ct, C0):
    C0 = 0.5 * (C0 + C0.T)
    w, U = np.linalg.eigh(C0)
    keep = w > RTOL * w.max()
    Uk = U[:, keep] / np.sqrt(w[keep])
    M = Uk.T @ (0.5 * (Ct + Ct.T)) @ Uk
    return np.sort(np.linalg.eigvalsh(0.5 * (M + M.T)))[::-1]


def effmass(blkmean, nstate):
    DT = blkmean.shape[-1]
    ev = np.full((DT, nstate), np.nan)
    for t in range(DT):
        try:
            e = gevp(blkmean[:, :, t], blkmean[:, :, T0])
            m = min(nstate, len(e))
            ev[t, :m] = e[:m]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        r = ev[:-1] / ev[1:]
        r[r <= 0] = np.nan
        return np.log(r)


def main():
    tag = dc.ENS.split("nu0")[0]
    cache = max(glob.glob("sigma2_smear_cache_claude/sigma2_smear_%s_*cfg_sm%s_d%d_claude.npy"
                          % (tag.replace(".", "p"), "-".join(map(str, SMEARS)), SPLIT)), key=os.path.getsize)
    allC = np.load(cache)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    sel = [int(x) for x in OPS.split(",")]
    allC = allC[:, sel][:, :, sel]
    ncfg = allC.shape[0]
    labs = ["%s@%d" % (OPLAB[s // NSM], SMEARS[s % NSM]) for s in sel]
    print("# %s  ncfg=%d nop=%d T0=%d RTOL=%g" % (cache.split("/")[-1], ncfg, len(sel), T0, RTOL))
    print("# ops: %s" % ", ".join(labs))

    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c = effmass(blk.mean(0), NSTATE)
    ems = np.array([effmass(np.delete(blk, i, 0).mean(0), NSTATE) for i in range(nb)])
    em_e = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, 0))

    print("\n#  t | " + " ".join("m%d(err)     " % n for n in range(NSTATE)) + "  [ref 0.46, 0.62]")
    for t in range(T0, min(DTMAXP, em_c.shape[0])):
        print("#  %2d | " % t + " ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(NSTATE)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(em_c.shape[0])
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    ax.axhline(0.644, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(DTMAXP * 0.6, 0.65, r"$2m_{PS}=0.644$", fontsize=9, color="gray")
    ax.axhline(0.46, color="tab:green", ls=":", lw=1, alpha=0.5)
    cols = ["tab:green", "tab:red", "tab:blue", "tab:purple"]
    mkr = ["o", "s", "^", "D"]
    for n in range(NSTATE):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.2)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n % 4], marker=mkr[n % 4], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, DTMAXP)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"plain pruned GEVP  nop=%d  T0=%d rtol=%g  %s L1 %dcfg" % (len(sel), T0, RTOL, tag, ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_smear_plaingevp_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
