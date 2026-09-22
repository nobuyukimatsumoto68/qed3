# Interacting L1 VSH sp, ELL-general (ell=1 or 2), on the COMPLETE Nv=24 (=2*12) perams = EXACT all-to-all.
# Variant of interacting_vsh_L1_claude.py generalized to ell=2 (leaves the original ell=1 script + its npz
# untouched).  ELL=1 here reproduces the validated ell=1 numbers (sanity check).  Single source window.
# ell=2 at L1 (12 sites) ALIASES hard on the icosahedron (H is 5-dim) -> treat ell=2 L1 as unreliable.
# Same axial contraction f^{ab}_A = -tr[sigma^a G21^dag sigma^b G21].  See vsh_ell2_impl_plan_claude.md.
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
import glob
import re
import numpy as np
ENS = os.environ.get("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
NVDIR = os.environ.get("NVDIR", "distill_Nv24")
os.environ["ENS"] = ENS
os.environ["NVDIR"] = NVDIR
KMIN = int(os.environ.get("KMIN", "20"))
BINCFG = int(os.environ.get("BINCFG", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
ELL = int(os.environ.get("ELL", "1"))
import fs_gevp_point_claude as fg
sys.path.insert(0, "final/analysis_axial")
import effmass_axial_tp_l3_perm_hankel_claude as H

GEO = "/mnt/barracuda22/qed3/qed3/geometry/data"
AT = 0.2
L = 1
NSITE = 12
s1 = np.array([[0, 1], [1, 0]], complex)
s2 = np.array([[0, -1j], [1j, 0]], complex)
s3 = np.array([[1, 0], [0, -1]], complex)
SIG = {1: s1, 2: s2, 3: s3}
CH = ["EA", "MA", "TA"]


def ylm_weights(theta, phi, ell):
    # grad[m]=(df/dtheta, (1/sin t) df/dphi), ysc[m]=Y_lm ; pole-safe ; VSH 1/sqrt(l(l+1)) dropped (cancels in mass)
    st = np.sin(theta)
    ct = np.cos(theta)
    if ell == 1:
        c = np.sqrt(3.0 / (4.0 * np.pi))
        grad = {0: (-c * st, np.zeros_like(st)),
                1: (c * ct * np.cos(phi), -c * np.sin(phi)),
                -1: (c * ct * np.sin(phi), c * np.cos(phi))}
        ysc = {0: c * ct, 1: c * st * np.cos(phi), -1: c * st * np.sin(phi)}
        return grad, ysc, (-1, 0, 1)
    if ell == 2:
        A = 0.25 * np.sqrt(5.0 / np.pi)
        B = 0.5 * np.sqrt(15.0 / np.pi)
        C = 0.25 * np.sqrt(15.0 / np.pi)
        c2t = ct * ct - st * st
        s2t = 2.0 * st * ct
        grad = {0: (-6.0 * A * st * ct, np.zeros_like(st)),
                1: (B * c2t * np.cos(phi), -B * ct * np.sin(phi)),
                -1: (B * c2t * np.sin(phi), B * ct * np.cos(phi)),
                2: (C * s2t * np.cos(2 * phi), -2.0 * C * st * np.sin(2 * phi)),
                -2: (C * s2t * np.sin(2 * phi), 2.0 * C * st * np.cos(2 * phi))}
        ysc = {0: A * (3.0 * ct * ct - 1.0),
               1: B * st * ct * np.cos(phi),
               -1: B * st * ct * np.sin(phi),
               2: C * st * st * np.cos(2 * phi),
               -2: C * st * st * np.sin(2 * phi)}
        return grad, ysc, (-2, -1, 0, 1, 2)
    raise ValueError("ell must be 1 or 2")


def config_corr(k, Wt, Ysc, mlist, w, dmax):
    AblkS, AblkSt, twin, nsite, U, tau = fg.make_config(k)
    g = {c: np.zeros(dmax) for c in CH}
    for dt in range(1, dmax + 1):
        acc = {c: 0.0 for c in CH}
        nc = 0
        for t0 in range(0, twin - dt):
            G21 = AblkS(t0 + dt, t0).transpose(2, 0, 1, 3)
            G21d = np.conj(G21.transpose(0, 1, 3, 2))
            fA = {}
            for a in (1, 2, 3):
                for b in (1, 2, 3):
                    if (a == 3) != (b == 3):
                        continue
                    fA[(a, b)] = -np.einsum('pq,ijqr,rs,ijsp->ij', SIG[a], G21d, SIG[b], G21, optimize=True)
            for m in mlist:
                wt, wp = Wt[m]
                y = Ysc[m]
                WE = {1: wt, 2: wp}
                WM = {1: -wp, 2: wt}
                for a in (1, 2):
                    for b in (1, 2):
                        acc["EA"] += np.sum(np.outer(w * WE[a], w * WE[b]) * fA[(a, b)])
                        acc["MA"] += np.sum(np.outer(w * WM[a], w * WM[b]) * fA[(a, b)])
                acc["TA"] += np.sum(np.outer(w * y, w * y) * fA[(3, 3)])
            nc += 1
        for c in CH:
            g[c][dt - 1] = (acc[c] / nc).real
    return g


def emcurve(C, off):
    sgn = np.sign(C[3])
    em, _ = H.hankel_effmass_scalar((sgn * C).astype(float), off, 4, 1, 1, AT)
    return em[:, 0]


def main():
    ks = sorted(int(re.search(r'peram\.(\d+)\.h5', f).group(1))
                for f in glob.glob("data_%s/%s/peram.*.h5" % (ENS, NVDIR)))
    ks = [k for k in ks if k >= KMIN]
    if NCFG:
        ks = ks[:NCFG]
    dmax = 31
    pts = np.loadtxt("%s/pts_n%d.dat" % (GEO, L))[:NSITE]
    n = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    theta = np.arccos(np.clip(n[:, 2], -1, 1))
    phi = np.arctan2(n[:, 1], n[:, 0])
    w = np.ones(NSITE)
    Wt, Ysc, mlist = ylm_weights(theta, phi, ELL)
    print("# interacting VSH L1 (EXACT Nv=24) %s : %d cfg (kmin=%d), ELL=%d AXIAL"
          % (ENS.split("_hb")[0], len(ks), KMIN, ELL), flush=True)

    gall = {c: [] for c in CH}
    for i, k in enumerate(ks):
        g = config_corr(k, Wt, Ysc, mlist, w, dmax)
        for c in CH:
            gall[c].append(g[c])
        if (i + 1) % 50 == 0:
            print("#  %d/%d cfg" % (i + 1, len(ks)), flush=True)
    nb = len(ks) // BINCFG
    off = [0, 3, 6]
    savez = {}
    out = {}
    for c in CH:
        allc = np.array(gall[c])[:nb * BINCFG]
        binned = allc.reshape(nb, BINCFG, -1).mean(1)
        savez[c + "_binned"] = binned
        cen = binned.mean(0)
        jk = np.array([np.delete(binned, b, 0).mean(0) for b in range(nb)])
        emc = emcurve(cen, off)
        emj = np.array([emcurve(jk[b], off) for b in range(nb)])
        emerr = np.sqrt((nb - 1) * np.nanmean((emj - np.nanmean(emj, 0)) ** 2, 0))
        savez[c + "_em"] = emc
        savez[c + "_err"] = emerr
        out[c] = (emc, emerr)
    tag = "ell%d" % ELL
    np.savez("final/analysis_axial/interacting_vsh_axial_L1_nf2g1_%s_claude.npz" % tag, **savez)

    print("#  dt |  electric Phi      magnetic Psi      tp radial   (ELL=%d)" % ELL, flush=True)
    ncur = len(out["EA"][0])
    for t in range(ncur - 1):
        row = "#  %2d |" % (t + 1)
        for c in ["EA", "MA", "TA"]:
            emc, eme = out[c]
            row += " %7.4f(%.4f)" % (emc[t], eme[t])
        print(row, flush=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"EA": "tab:green", "MA": "tab:red", "TA": "tab:blue"}
    mk = {"EA": "s", "MA": "^", "TA": "o"}
    nm = {"EA": "electric $\\Phi_%d$" % ELL, "MA": "magnetic $\\Psi_%d$" % ELL, "TA": "tp radial ell=%d" % ELL}
    dd = np.arange(1, ncur + 1)
    plt.figure(figsize=(7.6, 5.0))
    for c in ["TA", "EA", "MA"]:
        emc, eme = out[c]
        gmask = np.isfinite(emc) & np.isfinite(eme) & (eme < 0.2)
        plt.errorbar(dd[gmask], emc[gmask], yerr=eme[gmask], marker=mk[c], ms=4, capsize=1.5, lw=0.8, color=col[c], label=nm[c])
    plt.xlabel("dt")
    plt.ylabel("GPOF effmass = $a_t m$")
    plt.title("Interacting L1 Nf2 g1.0 AXIAL sp VSH ell=%d%s" % (ELL, " (ALIASED)" if ELL == 2 else ""))
    plt.legend(fontsize=9)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("final/analysis_axial/interacting_vsh_axial_L1_nf2g1_%s_claude.png" % tag, dpi=150)
    print("# wrote npz + png (%s)" % tag, flush=True)


if __name__ == "__main__":
    main()
