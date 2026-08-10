# effmass_plateau_tables_claude.py
# Dump the AXIAL-current cosh effective mass vs t for every massless ensemble (L=1,2,3,4) to a .md,
# so the constant-fit WINDOWS can be chosen from numbers instead of by eye off the png.
# This is the raw m_eff(t) table; the FITTED constants live in effmass_fits_claude.md.
#   m_eff(t) = acosh( (C[t-1]+C[t+1]) / (2 C[t]) ) / a_t ,  config-jackknife error, kmin=20.
#   t = a_t n_t (a_t = 0.2), i.e. the same variable as the FITW windows in effmass_conn_claude.py.
# Nf = 2,4,6 side by side for l=1 and l=2 in one table per (L, gsq).
# Output: effmass_plateau_tables_claude.md (md files stay in the production root; figs/ is for pngs).
import glob, math, re
import numpy as np
import h5py

Nt = 128
at = 0.2
KMIN = 20
NFS = [2, 4, 6]
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
HB = {1: "1.000000", 2: "1.000000", 3: "0.400000-1.000000", 4: "0.400000-1.000000"}
# current windows (effmass_conn_claude.py FITW) -- marked in the table so the tail is visible
FITW = {1: (4.0, 5.2), 2: (2.4, 4.0), 3: (2.4, 4.0), 4: (2.4, 4.0)}
TMIN, TMAX = 1.0, 6.2
dt = np.arange(1, Nt // 2)
tt = dt * at


def esn(nf, L, g):
    return ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, L, HB[L]))


def files(nf, L, g):
    fs = sorted(glob.glob(esn(nf, L, g) + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
    return [f for f in fs if int(re.search(r"corr\.(\d+)\.", f).group(1)) >= KMIN]


def load(fs, key):
    out = []
    for fn in fs:
        with h5py.File(fn, "r") as f:
            out.append(f[key + "/real"][()] + 1j * f[key + "/imag"][()])
    return np.array(out)


def C_axial(fs, l):
    s = sum(load(fs, "h0/ylm_axial/s3/l%d/m%d/Vpp" % (l, m)) for m in range(-l, l + 1))
    return 2.0 * (s / (2 * l + 1)).real


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


def effmass_jk(samp):
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return mean, err


out = ["# Axial-current effective mass $m_\\mathrm{eff}(t)$ tables (redo massless conn)", "",
       "Raw plateau data for CHOOSING the constant-fit windows; the fitted constants are in",
       "`effmass_fits_claude.md`.  $m_\\mathrm{eff}(t)=\\mathrm{acosh}[(C_{t-1}+C_{t+1})/2C_t]/a_t$,",
       "config-jackknife error, kmin=20, $a_t=0.2$, $t=a_t n_t$ (same variable as `FITW`).", "",
       "`*` marks $t$ inside the CURRENT window (L1 [4.0,5.2]; L2/L3/L4 [2.4,4.0], the latter two",
       "PROVISIONAL copies of L2).  Watch the $\\ell=2$ error growth at L3/L4 near the upper end.", ""]

for L in [1, 2, 3, 4]:
    lo, hi = FITW[L]
    for g in GS[L]:
        got = {}
        for nf in NFS:
            fs = files(nf, L, g)
            if len(fs) < 4:
                continue
            got[nf] = (len(fs), effmass_jk(C_axial(fs, 1).astype(float)),
                       effmass_jk(C_axial(fs, 2).astype(float)))
        if not got:
            continue
        ns = ", ".join("Nf%d:%d" % (nf, got[nf][0]) for nf in sorted(got))
        out.append("## L%d  gsq %.1f   (ncfg %s; current window [%.1f, %.1f])" % (L, g, ns, lo, hi))
        out.append("")
        hdr = "| t | " + " | ".join("Nf%d l=1" % nf for nf in sorted(got)) \
              + " | " + " | ".join("Nf%d l=2" % nf for nf in sorted(got)) + " |"
        out.append(hdr)
        out.append("|" + "---|" * (1 + 2 * len(got)))
        for d in dt:
            t = d * at
            if t < TMIN or t > TMAX:
                continue
            mark = "*" if (lo <= t <= hi) else ""
            row = "| %.1f%s " % (t, mark)
            for il in (1, 2):
                for nf in sorted(got):
                    mean, err = got[nf][il]
                    row += "| %.3f(%.3f) " % (mean[d], err[d])
            out.append(row + "|")
        out.append("")
        print("  L%d g%.1f done" % (L, g))

open("effmass_plateau_tables_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote effmass_plateau_tables_claude.md")
