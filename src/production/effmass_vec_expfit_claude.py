# effmass_vec_expfit_claude.py
# Effective-mass fit of the vector-current PHYSICAL full correlator (conn - disc) for l=1,2.
#  full(dt) = -C_conn + (disc_DCsub - plateau[dt 16..Nt/2])   (ipynb scheme, Eq.12 PhysRevD.110.054501)
# EXCITED-STATE effmass fit  m_eff(dt) = m0 + A exp(-c1 dt)  over dt in [8,14]: m0 = plateau mass, the
# A exp(-c1 dt) term follows the residual approach so the fit CURVE bends onto the still-sloping points.
# cosh effmass m_eff(dt) = acosh((C[t-1]+C[t+1])/(2C[t]))/at.  Config-jackknife: fit each delete-1
# sample -> m0^b, jackknife error on m0.  One png per (L,g): l=1,2 effmass (Nf 2/4/6) + fit curve + m0.
import glob, math, re
import numpy as np, h5py
from scipy.optimize import curve_fit
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
LSET = [1]                            # only l=1 is meaningful with current stats
LLIST = [1]                           # L1 only (L2 too noisy for the m0+A exp fit)
PWIN = np.arange(16, Nt // 2 + 1)     # disc plateau window (dt)
FITDT = np.arange(8, 15)             # m0+A exp fit range dt = 8..14
FITW_AX = {1: (4.0, 5.2)}            # axial l=1 conn const-fit plateau (t = a_t n_t)
dtp = np.arange(1, Nt // 2)
tt = dtp * at


def esn(nf, L, g):
    return ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, L, HB[L]))


def kof(fn):
    return int(fn.rsplit('/', 1)[-1].split('.')[1])


def conn_files(nf, L, g):
    fs = sorted(glob.glob(esn(nf, L, g) + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
    return [f for f in fs if kof(f) >= KMIN]


def disc_files(nf, L, g):
    fs = sorted(glob.glob(esn(nf, L, g) + "corr_ylm_disc_tb2_nhits1/corr.*.h0.h5"))
    return [f for f in fs if kof(f) >= KMIN]


def matched(nf, L, g):
    cfm = {kof(f): f for f in conn_files(nf, L, g)}
    dfm = {kof(f): f for f in disc_files(nf, L, g)}
    ks = sorted(set(cfm) & set(dfm))
    return [cfm[k] for k in ks], [dfm[k] for k in ks]


def load(fs, key):
    o = []
    for fn in fs:
        with h5py.File(fn, 'r') as f:
            o.append(f[key + '/real'][()] + 1j * f[key + '/imag'][()])
    return np.array(o)


def gl_conn(fs, grp, l):
    return sum(load(fs, "h0/%s/s3/l%d/m%d/Vpp" % (grp, l, m)) for m in range(-l, l + 1)) / (2 * l + 1)


def two_point_pp(Jc):
    G = np.zeros(Jc.shape, dtype=complex)
    for dt in range(Nt):
        G[:, dt] = np.mean(Jc * np.roll(Jc, -dt, axis=1), axis=1)
    return G


def disc_samples(dfs, l):
    Jm = [load(dfs, "h0/disc/ylm/s3/l%d/m%d/J" % (l, m)) for m in range(-l, l + 1)]
    return inv4pi * sum(two_point_pp(J) for J in Jm) / (2 * l + 1)


def vec_full(nf, L, g, l):
    cfm, dfm = matched(nf, L, g)
    if len(cfm) < 4:
        return None
    conn = -gl_conn(cfm, 'ylm', l).real
    samp = disc_samples(dfm, l).real
    sub = samp - samp.mean(axis=1, keepdims=True)
    plat = sub[:, PWIN].mean(axis=1)
    full = conn + (sub - plat[:, None])
    return full


def meff_acosh(C):
    m = np.full(Nt, np.nan)
    for t in range(1, Nt - 1):
        d = 2.0 * C[t]
        if d == 0:
            continue
        r = (C[t - 1] + C[t + 1]) / d
        if r > 1.0:
            m[t] = np.arccosh(r)
    return m


def effmass_jk(samp):
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at   # (H, Nt) per-sample effmass
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return me, mean, err


def model(dt, m0, A, c1):
    return m0 + A * np.exp(-c1 * dt)


def m0fit(me, mean, err):
    # fit m_eff(dt)=m0+A exp(-c1 dt) over FITDT; central fit + config-jackknife on m0.
    # returns (m0, sig, popt_central) or (None, None, None)
    x = FITDT.astype(float)
    if np.any(~np.isfinite(mean[FITDT])) or np.any(err[FITDT] <= 0):
        return None, None, None
    p0 = [mean[FITDT[-1]], mean[FITDT[0]] - mean[FITDT[-1]], 0.3]
    bnd = ([0.0, -np.inf, 0.0], [np.inf, np.inf, np.inf])
    try:
        popt, _ = curve_fit(model, x, mean[FITDT], sigma=err[FITDT], p0=p0,
                            bounds=bnd, maxfev=20000)
    except (RuntimeError, ValueError):
        return None, None, None
    H = me.shape[0]
    m0b = np.full(H, np.nan)
    for b in range(H):
        yb = me[b, FITDT]
        if np.any(~np.isfinite(yb)):
            continue
        try:
            pb, _ = curve_fit(model, x, yb, sigma=err[FITDT], p0=popt, bounds=bnd, maxfev=20000)
            m0b[b] = pb[0]
        except (RuntimeError, ValueError):
            continue
    good = np.isfinite(m0b)
    if good.sum() < H // 2:
        return None, None, None
    mb = m0b[good]
    n = mb.size
    M = mb.mean()
    sig = math.sqrt((n - 1.0) / n * np.sum((mb - M) ** 2))
    return M, sig, popt


# ---------- axial l=1 conn mass (cosh effmass + const fit) ----------
def C_axial(fs, l):
    return 2.0 * (sum(load(fs, "h0/ylm_axial/s3/l%d/m%d/Vpp" % (l, m)) for m in range(-l, l + 1)) / (2 * l + 1)).real


def axial_l1(nf, L, g):
    fs = conn_files(nf, L, g)
    if len(fs) < 4:
        return None, None
    me, mean, err = effmass_jk(C_axial(fs, 1).astype(float))
    lo, hi = FITW_AX[L]
    idx = dtp[(tt >= lo) & (tt <= hi)]
    good = idx[np.isfinite(mean[idx]) & (err[idx] > 0)]
    if len(good) < 2:
        return None, None
    w = 1.0 / err[good] ** 2
    w /= w.sum()
    Mb = (me[:, good] * w).sum(1)
    M = Mb.mean()
    H = me.shape[0]
    sig = math.sqrt((H - 1.0) / H * np.sum((Mb - M) ** 2))
    return M, sig


rows = []
for L in LLIST:
    for g in GS[L]:
        fig, axs = plt.subplots(1, len(LSET), figsize=(6.5 * len(LSET), 5))
        drew = False
        xc = np.linspace(FITDT[0], FITDT[-1], 60)
        for l, ax in zip(LSET, np.atleast_1d(axs)):
            for nf in NFS:
                full = vec_full(nf, L, g, l)
                if full is None:
                    continue
                full = full.astype(float)
                me, mean, err = effmass_jk(full)
                M, sig, popt = m0fit(me, mean, err)
                ax.errorbar(dtp, mean[dtp], yerr=err[dtp], fmt='o', ms=2.5, capsize=1.5,
                            color=nfcol[nf], label="Nf%d" % nf)
                if M is not None:
                    ax.plot(xc, model(xc, *popt), color=nfcol[nf], lw=1.5)          # m0 + A exp(-c1 dt)
                    ax.axhline(M, color=nfcol[nf], ls="--", lw=0.9, alpha=0.7)       # plateau m0
                    ax.fill_between([FITDT[0], FITDT[-1]], M - sig, M + sig,
                                    color=nfcol[nf], alpha=0.12, linewidth=0)
                    rows.append((L, g, l, nf, M, sig, full.shape[0]))
                drew = True
            ax.axvspan(FITDT[0], FITDT[-1], color="gray", alpha=0.1)
            ax.set_title("vector full l=%d  (m0 = plateau)" % l)
            ax.set_xlim(0, 30)
            ax.set_ylim(0, 6)
            ax.set_xlabel("dt")
            ax.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        fig.suptitle(r"VECTOR full (conn-disc): fit $m_0+A\,e^{-c_1 dt}$ over dt[%d,%d] -- L%d g%.1f"
                     % (FITDT[0], FITDT[-1], L, g))
        fig.tight_layout()
        if drew:
            fig.savefig("expfit_vec_L%d_g%.1f_claude.png" % (L, g), dpi=140)
        plt.close(fig)
        print("L%d g%.1f done" % (L, g))

# ---- table + md ----
lines = ["# Vector-current effmass-fit masses (conn-disc full, m0+A exp(-c1 dt), dt[%d,%d], kmin=20)"
         % (FITDT[0], FITDT[-1]), "",
         "a_t m0 from m_eff(dt)=m0+A exp(-c1 dt), config-jackknife.  m0 = plateau (vector mass).",
         "L1 l=1 only (l=2 / L2 too noisy at current stats).", "",
         "| L | gsq | l | Nf | a_t m0 | ncfg |", "|---|-----|---|----|--------|------|"]
print("\n# L gsq l Nf : a_t m0 (err)  ncfg")
for L, g, l, nf, M, sig, n in rows:
    print("  L%d g%.1f l%d Nf%d : %.4f(%.4f)  %d" % (L, g, l, nf, M, sig, n))
    lines.append("| %d | %.1f | %d | %d | %.4f(%.4f) | %d |" % (L, g, l, nf, M, sig, n))
open("vec_expfit_masses_claude.md", "w").write("\n".join(lines) + "\n")
print("# wrote vec_expfit_masses_claude.md (%d fits)" % len(rows))

# ---- vector l=1 / axial l=1 ratio (L1) vs gsq ----
# vector m0 = m0+A exp fit above; axial l=1 = conn cosh effmass const fit over [4.0,5.2].
# free ref: vector l=1 = sqrt2 (conserved current, 1 photon-like), axial l=1 = 2 -> sqrt2/2 = 1/sqrt2.
vmass = {(L, g, nf): (M, sig) for (L, g, l, nf, M, sig, n) in rows if l == 1}
fig, ax = plt.subplots(figsize=(7, 5))
for nf in NFS:
    xs, ys, es = [], [], []
    for g in GS[1]:
        if (1, g, nf) not in vmass:
            continue
        vm, ve = vmass[(1, g, nf)]
        am, ae = axial_l1(nf, 1, g)
        if am is None:
            continue
        r = vm / am
        er = r * math.sqrt((ve / vm) ** 2 + (ae / am) ** 2)
        xs.append(g)
        ys.append(r)
        es.append(er)
    if xs:
        ax.errorbar(xs, ys, es, marker="o", color=nfcol[nf], capsize=3, ls="-", lw=0.8, label="Nf%d" % nf)
ax.axhline(1.0 / math.sqrt(2.0), color="gray", ls="--", lw=1.0)
ax.text(0.02, 1.0 / math.sqrt(2.0), r" free $1/\sqrt{2}$", color="gray", fontsize=9, va="bottom",
        transform=ax.get_yaxis_transform())
ax.set_xlabel(r"$g_0^2$")
ax.set_ylabel(r"$\Delta_{V,\ell=1}/\Delta_{A,\ell=1}$")
ax.set_title(r"vector $\ell=1$ / axial $\ell=1$ (L1) vs $g_0^2$")
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("ratio_vec_over_axial_l1_L1_claude.png", dpi=150)
print("# wrote ratio_vec_over_axial_l1_L1_claude.png")
