# vec2_signal_claude.py
# SIGNAL CHECK for the vector^2 (0++) operator O_3(t) = sum_x J_mu(x,t) J^mu(x,t), built from the
# EXISTING 2-hit singlet-vector disc loops h0/disc/ylm/s{a}/l{l}/m{m}/J (a=1,2,3, l=0..3).
# By Parseval the zero-momentum current-square = tower sum  O_3(t) = sum_{a,l,m} J_{a,l,m}(t) conj(J_{a,l,m}(t)).
# UNBIASED estimator: use DISTINCT hits for the two loop factors (h0,h1) -> kills the eta-eta^dag diagonal.
#   O_3(t) = sum_{a,l,m} Re[ J^{h0}_{a,l,m}(t) conj(J^{h1}_{a,l,m}(t)) ] .
# Connected correlator (mirror vec_disc_only_claude.py conventions): source-time-avg two-point, then
#   per-config DC subtraction (remove <O_3>^2 time-mean) + plateau subtraction; config-jackknife SEM.
# Goal: is there a usable plateau / signal-to-noise BEFORE investing in the F^2 O_F(t) dump + GEVP.
# Physics: <O_3 O_3> vanishes only at large-Nf LO (Chester-Pufu Fig.11); finite-Nf it is O(1/Nf), so
#   a weak-but-nonzero signal is EXPECTED, growing toward small Nf.  See sigma_sigma_f2_mixing_impl_plan_claude.md.
import glob, math, re
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
CFG_STRIDE = 1          # subsample configs for speed if needed (1 = all matched)
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0]}
HB = {1: "1.000000", 2: "1.000000"}
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
PWIN = np.arange(16, Nt // 2 + 1)        # plateau region for the DC/plateau subtraction
AVEC = (1, 2, 3)
LMAX = 3
dtp = np.arange(1, Nt // 2)
tt = dtp * at


def esn(nf, L, g):
    return ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, L, HB[L]))


def kof(f):
    return int(re.search(r"corr\.(\d+)\.", f).group(1))


def matched_2hit(nf, L, g):
    # configs with BOTH hits present (h0 and h1), k >= KMIN
    d = esn(nf, L, g) + "corr_ylm_disc_tb2/"
    h0 = {kof(f): f for f in glob.glob(d + "corr.*.h0.h5") if kof(f) >= KMIN}
    h1 = {kof(f): f for f in glob.glob(d + "corr.*.h1.h5") if kof(f) >= KMIN}
    ks = sorted(set(h0) & set(h1))[::CFG_STRIDE]
    return [h0[k] for k in ks], [h1[k] for k in ks]


def load(fs, key):
    o = []
    for fn in fs:
        with h5py.File(fn, 'r') as f:
            o.append(f[key + '/real'][()] + 1j * f[key + '/imag'][()])
    return np.array(o)


def build_O3(h0fs, h1fs):
    # O_3(cfg, t) = sum_{a,l,m} Re[ J^h0 conj(J^h1) ]   (unbiased: distinct hits)
    O3 = None
    for a in AVEC:
        for l in range(LMAX + 1):
            for m in range(-l, l + 1):
                key = "h0/disc/ylm/s%d/l%d/m%d/J" % (a, l, m)
                j0 = load(h0fs, key)
                j1 = load(h1fs, key)
                term = (j0 * np.conj(j1)).real
                O3 = term if O3 is None else O3 + term
    return O3        # (Ncfg, Nt) real


def two_point(O3):
    G = np.zeros(O3.shape)
    for dt in range(Nt):
        G[:, dt] = np.mean(O3 * np.roll(O3, -dt, axis=1), axis=1)
    return G


def connected(O3):
    G = two_point(O3)
    sub = G - G.mean(axis=1, keepdims=True)       # remove <O_3>^2 (time-mean of the 2pt)
    plat = sub[:, PWIN].mean(axis=1)
    return sub - plat[:, None]                     # (Ncfg, Nt)


def jk_sem(samp):
    H = samp.shape[0]
    m = samp.mean(0)
    e = np.sqrt(np.sum((samp - m) ** 2, 0) / (H * (H - 1.0)))
    return m, e


def eff_acosh(C):
    m = np.full(Nt, np.nan)
    for t in range(1, Nt - 1):
        d = 2.0 * C[t]
        if d != 0:
            r = (C[t - 1] + C[t + 1]) / d
            if r > 1.0:
                m[t] = math.acosh(r)
    return m


print("# vector^2 (0++) signal check -- O_3 auto-correlator, unbiased (distinct hits)")
for L in (1, 2):
    fig, axs = plt.subplots(len(NFS), len(GS[L]), figsize=(6.0 * len(GS[L]), 4.2 * len(NFS)),
                            sharex=True, sharey=True, squeeze=False)
    for irow, nf in enumerate(NFS):
        for icol, g in enumerate(GS[L]):
            ax = axs[irow, icol]
            h0fs, h1fs = matched_2hit(nf, L, g)
            n = len(h0fs)
            if n >= 8:
                O3 = build_O3(h0fs, h1fs)
                C = connected(O3)
                m, e = jk_sem(C)
                ax.errorbar(tt, np.abs(m[dtp]), yerr=e[dtp], marker="o", ms=3, capsize=1.5, lw=0.8,
                            color=nfcol[nf], label="n=%d" % n)
                # S:N summary at small dt + crude plateau effmass
                sn = np.abs(m[1:6]) / np.where(e[1:6] > 0, e[1:6], np.nan)
                em = eff_acosh(m) / at
                emwin = em[3:9]
                emwin = emwin[np.isfinite(emwin)]
                emstr = ("%.2f" % np.nanmean(emwin)) if emwin.size else "nan"
                print("  L%d Nf%d g%.1f : n=%3d  S:N(dt1-5)=%s  meff[3:9]~%s" %
                      (L, nf, g, n, np.array2string(sn, precision=1, floatmode="fixed"), emstr))
            else:
                print("  L%d Nf%d g%.1f : n=%d (skip)" % (L, nf, g, n))
            ax.set_yscale("log")
            ax.axvspan(PWIN[0] * at, PWIN[-1] * at, color="gray", alpha=0.08)
            ax.set_title(r"L%d Nf%d g%.1f" % (L, nf, g))
            ax.set_xlim(0, at * Nt / 2)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
            if icol == 0:
                ax.set_ylabel(r"$|\langle O_3(t)O_3(0)\rangle_c|$")
            if irow == len(NFS) - 1:
                ax.set_xlabel(r"$t = dt\,a_t$")
    fig.suptitle(r"VECTOR$^2$ ($0^{++}$) connected auto-correlator $|\langle O_3 O_3\rangle_c|$  L%d  (unbiased, 2-hit)" % L)
    fig.tight_layout()
    fig.savefig("figs/vec2_signal_L%d_claude.png" % L, dpi=150)
    print("# wrote figs/vec2_signal_L%d_claude.png" % L)
