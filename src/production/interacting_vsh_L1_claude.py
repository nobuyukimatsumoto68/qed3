# Interacting L1 VSH sp: reconstruct the exact all-to-all overlap propagator from the COMPLETE-basis
# distillation perambulators (Nv=24=2*12), build f^{ab}(n1,n2), project onto VSH (electric Phi / magnetic
# Psi, ell=1), jackknife over configs, GPOF ground level.  VECTOR current (C_V++): f^{ab}=-tr[s^a A s^b A].
# electric target dim Delta+l-1 (= tp partner); magnetic dim Delta+l.
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
BINCFG = int(os.environ.get("BINCFG", "10"))     # jackknife bin size in CONFIG units
import fs_gevp_point_claude as fg
sys.path.insert(0, "final/analysis_axial")
import effmass_axial_tp_l3_perm_hankel_claude as H

GEO = "/mnt/barracuda22/qed3/qed3/geometry/data"
AT = 0.2
L = 1
s1 = np.array([[0, 1], [1, 0]], complex)
s2 = np.array([[0, -1j], [1j, 0]], complex)
s3 = np.array([[1, 0], [0, -1]], complex)
SIG = {1: s1, 2: s2, 3: s3}
CH = ["EA", "MA", "TA"]   # AXIAL electric/magnetic/tp (t0-leg = forward block, spin-daggered)


def ylm1_weights(theta, phi):
    # tangent-gradient (VSH) weights (pole-safe) AND scalar Y_1m values (for tp)
    c = np.sqrt(3.0 / (4.0 * np.pi))
    st = np.sin(theta)
    ct = np.cos(theta)
    grad = {0: (-c * st, np.zeros_like(st)),
            1: (c * ct * np.cos(phi), -c * np.sin(phi)),
            -1: (c * ct * np.sin(phi), c * np.cos(phi))}
    ysc = {0: c * ct, 1: c * st * np.cos(phi), -1: c * st * np.sin(phi)}
    return grad, ysc


def config_corr(k, Wt, Ysc, w, dmax):
    AblkS, AblkSt, twin, nsite, U, tau = fg.make_config(k)
    g = {c: np.zeros(dmax) for c in CH}
    for dt in range(1, dmax + 1):
        acc = {c: 0.0 for c in CH}
        nc = 0
        for t0 in range(0, twin - dt):
            G21 = AblkS(t0 + dt, t0).transpose(2, 0, 1, 3)          # forward leg R = S(n2@t <- n1@t0)
            # AXIAL (deter C_A+-, gamma5-herm): both legs are the SINGLE forward R, one spin-daggered:
            #   f^{ab}_A = -tr[sigma^a R^dag_spin sigma^b R].  Reproduces d4 axial tp (0.336 vs 0.3346).
            G21d = np.conj(G21.transpose(0, 1, 3, 2))               # spin conj-transpose of forward leg
            fA = {}
            for a in (1, 2, 3):
                for b in (1, 2, 3):
                    if (a == 3) != (b == 3):
                        continue
                    fA[(a, b)] = -np.einsum('pq,ijqr,rs,ijsp->ij', SIG[a], G21d, SIG[b], G21, optimize=True)
            for m in (-1, 0, 1):
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


def plateau(C, off, wlo, whi):
    sgn = np.sign(C[3]) if abs(C[3]) > 0 else 1.0
    em, _ = H.hankel_effmass_scalar((sgn * C).astype(float), off, 4, 1, 1, AT)
    seg = em[wlo:whi + 1, 0]
    seg = seg[np.isfinite(seg)]
    return float(np.mean(seg)) if len(seg) else np.nan


def main():
    ks = sorted(int(re.search(r'peram\.(\d+)\.h5', f).group(1))
                for f in glob.glob("data_%s/%s/peram.*.h5" % (ENS, NVDIR)))
    ks = [k for k in ks if k >= KMIN]
    twin_guess = 32
    dmax = twin_guess - 1
    pts = np.loadtxt("%s/pts_n%d.dat" % (GEO, L))[:12]
    n = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    theta = np.arccos(np.clip(n[:, 2], -1, 1))
    phi = np.arctan2(n[:, 1], n[:, 0])
    w = np.ones(12)
    Wt, Ysc = ylm1_weights(theta, phi)
    print("# interacting VSH L1 %s : %d configs (kmin=%d), ell=1 (V + AXIAL)" % (ENS.split("_hb")[0], len(ks), KMIN), flush=True)

    gall = {c: [] for c in CH}
    for i, k in enumerate(ks):
        g = config_corr(k, Wt, Ysc, w, dmax)
        for c in CH:
            gall[c].append(g[c])
        if (i + 1) % 50 == 0:
            print("#  %d/%d configs" % (i + 1, len(ks)), flush=True)
    nb = len(ks) // BINCFG
    off = [0, 3, 6]
    winmap = {"EA": (10, 15), "TA": (10, 15), "MA": (12, 15)}   # magnetic plateaus later
    lab = {"EA": "electric  (AXIAL)", "MA": "magnetic  (AXIAL)", "TA": "tp radial (AXIAL) [cf d4 tp]"}

    def emcurve(C):
        sgn = np.sign(C[3])
        em, _ = H.hankel_effmass_scalar((sgn * C).astype(float), off, 4, 1, 1, AT)
        return em[:, 0]
    out = {}
    savez = {}
    for c in CH:
        allc = np.array(gall[c])[:nb * BINCFG]
        binned = allc.reshape(nb, BINCFG, -1).mean(1)
        savez[c + "_binned"] = binned
        cen = binned.mean(0)
        jk = np.array([np.delete(binned, b, 0).mean(0) for b in range(nb)])
        emc = emcurve(cen)
        emj = np.array([emcurve(jk[b]) for b in range(nb)])                 # per-dt jk effmass
        emerr = np.sqrt((nb - 1) * np.nanmean((emj - np.nanmean(emj, 0)) ** 2, 0))
        wlo, whi = winmap[c]
        pc = plateau(cen, off, wlo, whi)
        pj = np.array([plateau(jk[b], off, wlo, whi) for b in range(nb)])
        perr = np.sqrt((nb - 1) * np.mean((pj - pj.mean()) ** 2))
        out[c] = (emc, emerr, pc, perr, (wlo, whi))
        savez[c + "_em"] = emc
        savez[c + "_err"] = emerr
        print("  %-30s a_t*m = %.4f +- %.4f  (win dt[%d,%d])" % (lab[c], pc, perr, wlo, whi), flush=True)
    np.savez("final/analysis_axial/interacting_vsh_axial_L1_nf2g1_claude.npz", **savez)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"EA": "tab:green", "MA": "tab:red", "TA": "tab:blue"}
    mk = {"EA": "s", "MA": "^", "TA": "o"}
    nm = {"EA": "electric $\\Phi$  %.4f(%d)", "MA": "magnetic $\\Psi$  %.4f(%d)", "TA": "tp radial  %.4f(%d)"}
    dd = np.arange(1, 16)
    plt.figure(figsize=(7.4, 4.9))
    for c in ["TA", "EA", "MA"]:
        emc, emerr, pc, perr = out[c]
        plt.errorbar(dd, emc[:15], yerr=emerr[:15], marker=mk[c], ms=3.5, capsize=1.5, lw=0.8,
                     color=col[c], label=nm[c] % (pc, round(perr * 1e4)))
    plt.axhline(0.3346, color="gray", ls="--", lw=0.8)
    plt.text(11, 0.3346 + 0.001, "d4 axial tp 0.3346", fontsize=7)
    plt.axhline(0.338, color="silver", ls=":", lw=0.9)
    plt.text(11, 0.338 + 0.001, "scalar-Y sp (aliased)", fontsize=7)
    plt.ylim(0.30, 0.44)
    plt.xlabel("dt")
    plt.ylabel("GPOF effmass = $a_t m$")
    plt.title("Interacting L1 Nf2 g1.0 AXIAL sp (VSH), per-dt jackknife errors")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("final/analysis_axial/interacting_vsh_axial_L1_nf2g1_claude.png", dpi=150)
    print("# wrote plot + npz", flush=True)


if __name__ == "__main__":
    main()
