# vec2_signal_fast_claude.py -- same as vec2_signal_claude.py but opens each h5 ONCE (reads all
# (a,l,m) J datasets per file in one pass) instead of reopening per key.  L2 only (active, well-populated).
# O_3(t) = sum_{a,l,m} Re[ J^h0_{a,l,m}(t) conj(J^h1_{a,l,m}(t)) ]  (unbiased vector^2, distinct hits).
import sys, glob, math, re
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0]}
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
PWIN = np.arange(16, Nt // 2 + 1)
AVEC = (1, 2, 3)
LMAX = 3
dtp = np.arange(1, Nt // 2)
tt = dtp * at
KEYS = [(a, l, m) for a in AVEC for l in range(LMAX + 1) for m in range(-l, l + 1)]


def esn(nf, L, g):
    return ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb1.000000"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, L))


def kof(f):
    return int(re.search(r"corr\.(\d+)\.", f).group(1))


def read_all(fn):
    # one open, read every (a,l,m) J -> array (Nkeys, Nt) complex
    out = np.empty((len(KEYS), Nt), dtype=complex)
    with h5py.File(fn, 'r') as f:
        for i, (a, l, m) in enumerate(KEYS):
            g = f["h0/disc/ylm/s%d/l%d/m%d/J" % (a, l, m)]
            out[i] = g['real'][()] + 1j * g['imag'][()]
    return out


def build_O3(nf, L, g):
    d = esn(nf, L, g) + "corr_ylm_disc_tb2/"
    h0 = {kof(f): f for f in glob.glob(d + "corr.*.h0.h5") if kof(f) >= KMIN}
    h1 = {kof(f): f for f in glob.glob(d + "corr.*.h1.h5") if kof(f) >= KMIN}
    ks = sorted(set(h0) & set(h1))
    O3 = np.empty((len(ks), Nt))
    for j, k in enumerate(ks):
        j0 = read_all(h0[k])
        j1 = read_all(h1[k])
        O3[j] = (j0 * np.conj(j1)).real.sum(axis=0)
    return O3, len(ks)


def two_point(O3):
    G = np.zeros(O3.shape)
    for dt in range(Nt):
        G[:, dt] = np.mean(O3 * np.roll(O3, -dt, axis=1), axis=1)
    return G


def connected(O3):
    G = two_point(O3)
    sub = G - G.mean(axis=1, keepdims=True)
    plat = sub[:, PWIN].mean(axis=1)
    return sub - plat[:, None]


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


print("# vector^2 (0++) signal check (fast) -- L%d, unbiased distinct-hit O_3 auto-correlator" % L)
fig, axs = plt.subplots(len(NFS), len(GS[L]), figsize=(6.0 * len(GS[L]), 4.2 * len(NFS)),
                        sharex=True, sharey=True, squeeze=False)
for irow, nf in enumerate(NFS):
    for icol, g in enumerate(GS[L]):
        ax = axs[irow, icol]
        O3, n = build_O3(nf, L, g)
        if n >= 8:
            C = connected(O3)
            m, e = jk_sem(C)
            ax.errorbar(tt, np.abs(m[dtp]), yerr=e[dtp], marker="o", ms=3, capsize=1.5, lw=0.8,
                        color=nfcol[nf], label="n=%d" % n)
            sn = np.abs(m[1:6]) / np.where(e[1:6] > 0, e[1:6], np.nan)
            em = eff_acosh(m) / at
            w = em[3:9]
            w = w[np.isfinite(w)]
            emstr = ("%.2f" % np.nanmean(w)) if w.size else "nan"
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
            ax.set_ylabel(r"$|\langle O_3 O_3\rangle_c|$")
        if irow == len(NFS) - 1:
            ax.set_xlabel(r"$t=dt\,a_t$")
fig.suptitle(r"VECTOR$^2$ ($0^{++}$) connected auto-correlator, L%d (unbiased, 2-hit)" % L)
fig.tight_layout()
fig.savefig("figs/vec2sig_L%d_claude.png" % L, dpi=150)
print("# wrote figs/vec2sig_L%d_claude.png" % L)
