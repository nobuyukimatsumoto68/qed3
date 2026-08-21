# vec_disc_only_claude.py
# The DISCONNECTED contribution to the vector correlator, ALONE (isolated from -C_conn).
# Exactly the disc piece the full uses (effmass_vec_expfit_claude.py vec_full):
#   disc two-point  G(t) = <J(t) J(t+dt)>   (per config, m-averaged, inv4pi)
#   DC-subtract     sub  = G - <G>_t                       (per-config time mean removed)
#   residual plateau plat = mean over dt in [16, Nt/2]
#   disc contribution = sub - plat            <-- this is what is added to -C_conn
# Plotted vs physical time t = dt*a_t, config-jackknife band, l=1, Nf 2/4/6 overlaid.
# L1 (g0.5,1,1.5) and L2 (g1,2,3) -- where disc exists (~398 matched cfg, stride-10).
# For scale, the connected |C_conn| is overlaid faintly (dashed) so the disc size is comparable.
import glob, math, re
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
inv4pi = 1.0 / (4.0 * math.pi)
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0]}
HB = {1: "1.000000", 2: "1.000000"}
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
PWIN = np.arange(16, Nt // 2 + 1)
LVAL = 1
dtp = np.arange(1, Nt // 2)
tt = dtp * at


def esn(nf, L, g):
    return ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, L, HB[L]))


def kof(f):
    return int(re.search(r"corr\.(\d+)\.", f).group(1))


def matched(nf, L, g):
    cfm = {kof(f): f for f in glob.glob(esn(nf, L, g) + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5")
           if kof(f) >= KMIN}
    dfm = {kof(f): f for f in glob.glob(esn(nf, L, g) + "corr_ylm_disc_tb2/corr.*.h0.h5")
           if kof(f) >= KMIN}
    ks = sorted(set(cfm) & set(dfm))
    return [cfm[k] for k in ks], [dfm[k] for k in ks]


def load(fs, key):
    o = []
    for fn in fs:
        with h5py.File(fn, 'r') as f:
            o.append(f[key + '/real'][()] + 1j * f[key + '/imag'][()])
    return np.array(o)


def gl_conn(fs, l):
    return sum(load(fs, "h0/ylm/s3/l%d/m%d/Vpp" % (l, m)) for m in range(-l, l + 1)) / (2 * l + 1)


def two_point_pp(Jc):
    G = np.zeros(Jc.shape, dtype=complex)
    for dt in range(Nt):
        G[:, dt] = np.mean(Jc * np.roll(Jc, -dt, axis=1), axis=1)
    return G


def disc_contrib(dfs, l):
    Jm = [load(dfs, "h0/disc/ylm/s3/l%d/m%d/J" % (l, m)) for m in range(-l, l + 1)]
    samp = (inv4pi * sum(two_point_pp(J) for J in Jm) / (2 * l + 1)).real
    sub = samp - samp.mean(axis=1, keepdims=True)
    plat = sub[:, PWIN].mean(axis=1)
    return sub - plat[:, None]     # (Ncfg, Nt)  disc contribution per config


def jk(samp):
    # samp is RAW per-config data (Ncfg, Nt) -- NOT delete-1 jackknife means.  The plotted quantities
    # (m-averaged disc two-point after per-config DC+plateau subtraction; -C_conn) are LINEAR in the
    # per-config value, so the mean's uncertainty is the standard error of the mean = std/sqrt(H).
    # (Applying the jackknife (H-1)/H factor to raw samples over-inflates the error by a factor H-1.)
    H = samp.shape[0]
    m = samp.mean(0)
    e = np.sqrt(np.sum((samp - m) ** 2, 0) / (H * (H - 1.0)))
    return m, e


# grid: rows = Nf (separated), cols = gsq.  one Nf per panel -> disc (filled circle) + conn (open square).
for L in (1, 2):
    fig, axs = plt.subplots(len(NFS), len(GS[L]), figsize=(6.0 * len(GS[L]), 4.2 * len(NFS)),
                            sharex=True, sharey=True)
    for irow, nf in enumerate(NFS):
        for icol, g in enumerate(GS[L]):
            ax = axs[irow, icol]
            cfm, dfm = matched(nf, L, g)
            if len(dfm) >= 4:
                dc = disc_contrib(dfm, LVAL)
                md, ed = jk(dc)
                # log-y: plot the MAGNITUDE |disc| (the contribution changes sign after the plateau
                # subtraction; near-zero points fall off the bottom of the log axis).
                ax.errorbar(tt, np.abs(md[dtp]), yerr=ed[dtp], marker="o", ms=3, capsize=1.5, lw=0.8,
                            color=nfcol[nf], label="disc (n=%d)" % len(dfm))
                # connected |C_conn| curve (drawn, for comparison to the disc size)
                cc = -gl_conn(cfm, LVAL).real
                mc, ec = jk(cc)
                ax.errorbar(tt, np.abs(mc[dtp]), yerr=ec[dtp], marker="s", ms=2.5, mfc="none",
                            capsize=1.0, lw=0.9, ls="--", color="k", alpha=0.7, label="conn")
            ax.set_yscale("log")
            ax.axvspan(PWIN[0] * at, PWIN[-1] * at, color="gray", alpha=0.08)   # plateau region
            ax.set_title(r"L%d  Nf%d  g%.1f" % (L, nf, g))
            ax.set_xlim(0, at * Nt / 2)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
            if icol == 0:
                ax.set_ylabel(r"$|$disc contribution to $C_V(t)|$")
            if irow == len(NFS) - 1:
                ax.set_xlabel(r"physical time $t = dt\cdot a_t$")
    fig.suptitle(r"VECTOR $\ell=1$ DISCONNECTED contribution alone (circles) vs $|{-}C_\mathrm{conn}|$ (dashed). L%d" % L)
    fig.tight_layout()
    fig.savefig("figs/vec_disc_only_L%d_claude.png" % L, dpi=150)
    print("# wrote vec_disc_only_L%d_claude.png" % L)
