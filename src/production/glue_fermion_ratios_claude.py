# glue_fermion_ratios_claude.py
# Two glueball/fermion mass ratios for the massless redo ensembles, vs gsq and vs 1/L^2:
#   R2 = Delta_{F,l=1} / Delta_{A,l=1}      (glueball F l=1 over axial-current l=1)
#   R3 = Delta_{F^2,0++} / Delta_{F,l=1}    (0++ scalar glueball over F l=1)
# Glueball masses: GEVP per-m + fit_perm variance-avg (F l=1 [0.2,1.0], F^2 0++ [0.2,0.4]).
# Axial l=1: fermionic conn cosh effmass + diagonal const fit (windows L1[4.0,5.2], L2[2.4,4.0],
#   L4[2.4,4.0] provisional), config-jackknife.  kmin=20.
# L4 auto-included per ratio ONLY when the needed data exists (F^2/axial not measured at L4 yet ->
#   L4 skipped today; F l=1 alone is not enough for either ratio).
# Free single-particle refs on S2: F l=1 = sqrt2 (1 photon), F l=2 = sqrt3, F^2 0++ = 2 sqrt2 (2 photons),
#   axial l=1 = 2 (fermion, Delta_l = l+1).  -> R2_free = sqrt2/2 = 1/sqrt2 ; R3_free = 2.
import subprocess, glob, math, re
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
AX_MINCFG = 20
NFS = [2, 4, 6]
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 4: [2.0, 4.0, 6.0]}
HB = {1: "1.000000", 2: "1.000000", 4: "0.400000-1.000000"}
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
Lmk = {1: "o", 2: "s", 4: "^"}
invL2 = {1: 1.0, 2: 0.25, 4: 0.0625}
dxNf = {2: -0.012, 4: 0.0, 6: 0.012}
dt = np.arange(1, Nt // 2)
tt = dt * at
FITW_AX = {1: (4.0, 5.2), 2: (2.4, 4.0), 4: (2.4, 4.0)}


# ---------- glueball mass via fit_perm ----------
def glue_mass(pref, tag, states, tlo, thi):
    jk = "%s_%s_jk_claude.dat" % (pref, tag)
    try:
        open(jk).close()
    except OSError:
        return None, None
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(tlo), str(thi), states],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "variance-avg" in line or (states == "0" and "state 0:" in line):
            seg = line.split("M =")[-1] if "==>" in line else line.split("M=")[-1]
            m = seg.strip().split()[0]
            try:
                v = float(m.split("(")[0])
                e = float(m.split("(")[1].rstrip(")"))
            except (ValueError, IndexError):
                return None, None
            if not np.isfinite(v):
                return None, None
            return v, e
    return None, None


# ---------- axial l=1 mass via conn cosh effmass ----------
def esn(nf, L, g):
    return ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, L, HB[L]))


def afiles(nf, L, g):
    fs = sorted(glob.glob(esn(nf, L, g) + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
    return [f for f in fs if int(re.search(r"corr\.(\d+)\.", f).group(1)) >= KMIN]


def load(fs, key):
    o = []
    for fn in fs:
        with h5py.File(fn, "r") as f:
            o.append(f[key + "/real"][()] + 1j * f[key + "/imag"][()])
    return np.array(o)


def C_axial(fs, l):
    return 2.0 * (sum(load(fs, "h0/ylm_axial/s3/l%d/m%d/Vpp" % (l, m)) for m in range(-l, l + 1)) / (2 * l + 1)).real


def meff_acosh(C):
    m = np.full(Nt, np.nan)
    for t in range(1, Nt - 1):
        den = 2.0 * C[t]
        if den == 0:
            continue
        r = (C[t - 1] + C[t + 1]) / den
        if r > 1.0:
            m[t] = np.arccosh(r)
    return m


def axial_mass(nf, L, g):
    fs = afiles(nf, L, g)
    if len(fs) < AX_MINCFG:
        return None, None
    samp = C_axial(fs, 1).astype(float)
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    lo, hi = FITW_AX[L]
    idx = dt[(tt >= lo) & (tt <= hi)]
    good = idx[np.isfinite(mean[idx]) & (err[idx] > 0)]
    if len(good) < 2:
        return None, None
    w = 1.0 / err[good] ** 2
    w /= w.sum()
    Mb = (me[:, good] * w).sum(1)
    M = Mb.mean()
    sig = math.sqrt((H - 1) * np.mean((Mb - M) ** 2))
    return M, sig


def rerr(m1, e1, m2, e2):
    r = m1 / m2
    return r, r * math.sqrt((e1 / m1) ** 2 + (e2 / m2) ** 2)


# ---------- collect rows: (L, nf, g, ratio, err) ----------
def collect(numfun, denfun):
    rows = []
    for L in [1, 2, 4]:
        for nf in NFS:
            for g in GS[L]:
                tag = "Nf%d_g%.6f_L%d" % (nf, g, L)
                m1, e1 = numfun(nf, L, g, tag)
                m2, e2 = denfun(nf, L, g, tag)
                if m1 is None or m2 is None:
                    continue
                r, er = rerr(m1, e1, m2, e2)
                rows.append((L, nf, g, r, er))
    return rows


F_l1 = lambda nf, L, g, tag: glue_mass("gevp_F", tag, "0,1,2", 0.2, 1.0)
F2_0 = lambda nf, L, g, tag: glue_mass("gevp_f2", tag, "0", 0.2, 0.4)
Axl1 = lambda nf, L, g, tag: axial_mass(nf, L, g)

RATIOS = [
    ("F_over_axial", collect(F_l1, Axl1), 1.0 / math.sqrt(2.0), r"free $1/\sqrt{2}$",
     r"$\Delta_{F,\ell=1}/\Delta_{A,\ell=1}$", "F l=1 / axial l=1"),
    ("F2_over_F", collect(F2_0, F_l1), 2.0, r"free $2$",
     r"$\Delta_{F^2,0^{++}}/\Delta_{F,\ell=1}$", "F^2 0++ / F l=1"),
]

for key, rows, free, freelab, ylab, title in RATIOS:
    print("# %s : %d points" % (key, len(rows)))
    for L, nf, g, r, er in rows:
        print("  L%d Nf%d g%.1f : %.4f(%.4f)" % (L, nf, g, r, er))
    # vs gsq
    fig, ax = plt.subplots(figsize=(7, 5))
    for L in [1, 2, 4]:
        for nf in NFS:
            d = [x for x in rows if x[0] == L and x[1] == nf]
            if not d:
                continue
            ax.errorbar([x[2] for x in d], [x[3] for x in d], [x[4] for x in d],
                        marker=Lmk[L], color=nfcol[nf], capsize=3, ls="-", lw=0.8,
                        mfc="none" if L == 2 else nfcol[nf], label="Nf%d L%d" % (nf, L))
    ax.axhline(free, color="gray", ls="--", lw=1.0)
    ax.text(0.02, free, " " + freelab, color="gray", fontsize=9, va="bottom",
            transform=ax.get_yaxis_transform())
    ax.set_xlabel(r"$g_0^2$")
    ax.set_ylabel(ylab)
    ax.set_title(title + r" vs $g_0^2$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig("ratio_%s_vs_gsq_claude.png" % key, dpi=150)
    print("  wrote ratio_%s_vs_gsq_claude.png" % key)
    # vs 1/L^2
    fig, ax = plt.subplots(figsize=(7, 5))
    for nf in NFS:
        for g in GS[1] + GS[2] + GS[4]:
            d = [x for x in rows if x[1] == nf and abs(x[2] - g) < 1e-9]
            if not d:
                continue
            ax.errorbar([invL2[x[0]] + dxNf[nf] for x in d], [x[3] for x in d], [x[4] for x in d],
                        marker="o", color=nfcol[nf], capsize=3, ls="none",
                        label="Nf%d" % nf if abs(g - GS[1][0]) < 1e-9 else None)
    ax.axhline(free, color="gray", ls="--", lw=1.0)
    ax.text(0.02, free, " " + freelab, color="gray", fontsize=9, va="bottom",
            transform=ax.get_yaxis_transform())
    ax.set_xlabel(r"$1/L^2$")
    ax.set_ylabel(ylab)
    ax.set_title(title + r" vs $1/L^2$")
    ax.set_xlim(-0.05, 1.1)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig("ratio_%s_vs_invL2_claude.png" % key, dpi=150)
    print("  wrote ratio_%s_vs_invL2_claude.png" % key)
