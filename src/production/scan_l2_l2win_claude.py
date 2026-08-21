# scan_l2_l2win_claude.py
# For the axial tp ell=2, L=2 const+exp m0 fit, SCAN the candidate upper-cut windows tried so far and
# report the m0 error per (gsq, Nf) for each window, then pick the MIN-error window per point.
# (The l2/l1 ratio error is dominated by the l2 m0 error -- l1 is well-determined at [10,25] -- so
# minimizing the l2 m0 error minimizes the ratio noise.)  Prints a table + the recommended override
# dict entries keyed (ell=2, L=2, nf, gsq).
import glob, math, re
import numpy as np, h5py
from scipy.optimize import curve_fit

Nt = 128
at = 0.2
KMIN = 20
L = 2
LVAL = 2
GS = [1.0, 2.0, 3.0]
NFS = [2, 4, 6]
HB = "1.000000"
WINS = [(10, 18), (10, 21), (10, 23), (10, 25)]
dtp = np.arange(1, Nt // 2)


def esn(nf, g):
    return ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, L, HB))


def conn_files(nf, g):
    fs = sorted(glob.glob(esn(nf, g) + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
    return [f for f in fs if int(re.search(r"corr\.(\d+)\.", f).group(1)) >= KMIN]


def load(fs, key):
    o = []
    for fn in fs:
        with h5py.File(fn, 'r') as f:
            o.append(f[key + '/real'][()] + 1j * f[key + '/imag'][()])
    return np.array(o)


def C_axial(fs, l):
    return 2.0 * (sum(load(fs, "h0/ylm_axial/s3/l%d/m%d/Vpp" % (l, m)) for m in range(-l, l + 1)) / (2 * l + 1)).real


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
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return me, mean, err


def model(dt, m0, A, c1):
    return m0 + A * np.exp(-c1 * dt)


def m0fit(me, mean, err, D):
    if np.any(~np.isfinite(mean[D])) or np.any(err[D] <= 0):
        return None, None
    x = D.astype(float)
    p0 = [mean[D[-1]], mean[D[0]] - mean[D[-1]], 0.3]
    bnd = ([0.0, -np.inf, 0.0], [np.inf, np.inf, np.inf])
    try:
        popt, _ = curve_fit(model, x, mean[D], sigma=err[D], p0=p0, bounds=bnd, maxfev=40000)
    except (RuntimeError, ValueError):
        return None, None
    H = me.shape[0]
    m0b = np.full(H, np.nan)
    for b in range(H):
        yb = me[b, D]
        if np.any(~np.isfinite(yb)):
            continue
        try:
            pb, _ = curve_fit(model, x, yb, sigma=err[D], p0=popt, bounds=bnd, maxfev=40000)
            m0b[b] = pb[0]
        except (RuntimeError, ValueError):
            continue
    good = np.isfinite(m0b)
    if good.sum() < H // 2:
        return None, None
    mb = m0b[good]
    n = mb.size
    return mb.mean(), math.sqrt((n - 1.0) / n * np.sum((mb - mb.mean()) ** 2))


print("# axial tp ell=2 L=2 m0(err) per window; * = min-error pick")
print("# g   Nf | " + " | ".join("[%d,%d]" % w for w in WINS) + "  -> pick")
best = {}
for g in GS:
    for nf in NFS:
        fs = conn_files(nf, g)
        if len(fs) < 4:
            continue
        me, mean, err = effmass_jk(C_axial(fs, LVAL).astype(float))
        errs = []
        for lo, hi in WINS:
            D = np.arange(lo, hi + 1)
            M, sig = m0fit(me, mean, err, D)
            errs.append(sig)
        vals = [e if e is not None else np.inf for e in errs]
        bi = int(np.argmin(vals))
        best[(nf, g)] = WINS[bi]
        cells = []
        for i, e in enumerate(errs):
            s = "%.4f" % e if e is not None else "  --  "
            cells.append(("*%s*" % s) if i == bi else (" %s " % s))
        print("  %.1f  %d | %s  -> [%d,%d]" % (g, nf, " | ".join(cells), WINS[bi][0], WINS[bi][1]))

print("\n# recommended FITDT_OVR_NF entries (ell=2,L=2,nf,gsq): only where != base [10,25]")
for (nf, g), w in sorted(best.items()):
    if w != (10, 25):
        print("  (2, 2, %d, %.1f): (%d, %d)," % (nf, g, w[0], w[1]))
