# MAGNETIC (Psi, ell=1) point-to-point GEVP with a SINGLE fixed point.
#   O_0 = projected magnetic VSH   : sum_n W^a_m(n) j^a(n)         (all 12 sites)
#   O_1 = point magnetic current   : W^a_m(n0) j^a(n0)            (ONE fixed vertex n0, no sum)
# The point op is NOT projected onto a single VSH mode -> overlaps the full tower -> distinct excited
# content from O_0 -> the 2x2 GEVP separates the magnetic ground from contamination.  Picking a SINGLE
# point (not averaging) is what keeps O_1 distinct from O_0 (averaging over sites collapses it back,
# multiplicity=1) and gives a healthy C_11 (real local magnetic density), unlike the isotropic splits.
#
# Master object (per config, per dt), from the AXIAL contraction f^{ab}=-tr[sigma^a G21^dag sigma^b G21]:
#   M_{n1,n2}(t) = sum_m sum_{a,b in 1,2} W^a_m(n1) W^b_m(n2) f^{ab}_A(n1,n2;t)     [12x12]
#   C_00 = sum_{n1,n2} M ;  C_01 = sum_{n1} M[:,n0] ;  C_10 = sum_{n2} M[n0,:] ;  C_11 = M[n0,n0]
# See mag_point2point_gevp_impl_plan_claude.md.
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
T0G = int(os.environ.get("T0G", "3"))
N0 = int(os.environ.get("N0", "0"))       # the single fixed point (vertex index)
import fs_gevp_point_claude as fg
GEO = "/mnt/barracuda22/qed3/qed3/geometry/data"
AT = 0.2
L = 1
s1 = np.array([[0, 1], [1, 0]], complex)
s2 = np.array([[0, -1j], [1j, 0]], complex)
SIG = {1: s1, 2: s2}


def magw(theta, phi):
    # magnetic VSH tangent weight W^a_m(n): from grad Y_1m rotated 90 deg (hat n x grad).
    c = np.sqrt(3.0 / (4.0 * np.pi))
    st = np.sin(theta)
    ct = np.cos(theta)
    grad = {0: (-c * st, np.zeros_like(st)),
            1: (c * ct * np.cos(phi), -c * np.sin(phi)),
            -1: (c * ct * np.sin(phi), c * np.cos(phi))}
    return {m: {1: -grad[m][1], 2: grad[m][0]} for m in (-1, 0, 1)}


def config_M(k, WM, dmax):
    # 12x12 magnetic-directed site-pair matrix M_{n1,n2}(dt), averaged over source time t0.
    AblkS, _, twin, nsite, _, _ = fg.make_config(k)
    M = np.zeros((dmax, nsite, nsite))
    for dt in range(1, dmax + 1):
        acc = np.zeros((nsite, nsite))
        nc = 0
        for t0 in range(0, twin - dt):
            G21 = AblkS(t0 + dt, t0).transpose(2, 0, 1, 3)     # forward leg -> [n1=src, n2=sink, a, b]
            G21d = np.conj(G21.transpose(0, 1, 3, 2))          # spin conj-transpose (leg a<->b, conj)
            fA = {}
            for a in (1, 2):
                for b in (1, 2):
                    fA[(a, b)] = -np.einsum('pq,ijqr,rs,ijsp->ij', SIG[a], G21d, SIG[b], G21, optimize=True)
            for m in (-1, 0, 1):
                for a in (1, 2):
                    for b in (1, 2):
                        acc += np.outer(WM[m][a], WM[m][b]) * fA[(a, b)].real
            nc += 1
        M[dt - 1] = acc / nc
    return M


def gevp2_ground(C2, t0):
    # 2x2 GEVP ground effmass (largest generalized eigenvalue), oriented positive.
    dmax = C2.shape[0]
    lam = np.full(dmax, np.nan)
    C0 = 0.5 * (C2[t0] + C2[t0].T)
    wv, U = np.linalg.eigh(C0)
    if wv.max() <= 0:
        return np.full(dmax - 1, np.nan)
    keep = wv > 1e-10 * wv.max()
    Uk = U[:, keep] / np.sqrt(wv[keep])
    for t in range(dmax):
        Ct = 0.5 * (C2[t] + C2[t].T)
        try:
            ev = np.linalg.eigvalsh(Uk.T @ Ct @ Uk)
            lam[t] = ev[-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        return np.log(lam[:-1] / lam[1:])


def scalar_effmass(C):
    # log-ratio effmass of a single positive correlator (oriented).
    sgn = np.sign(C[T0G]) if abs(C[T0G]) > 0 else 1.0
    Cp = sgn * C
    with np.errstate(all="ignore"):
        return np.log(Cp[:-1] / Cp[1:])


def build_C2(Mbin, n0):
    # from binned 12x12 M -> 2x2 {projected, point-n0} for each bin (or the mean)
    nb = Mbin.shape[0]
    dmax = Mbin.shape[1]
    C2 = np.zeros((nb, dmax, 2, 2))
    for b in range(nb):
        for dt in range(dmax):
            Msym = 0.5 * (Mbin[b, dt] + Mbin[b, dt].T)
            C2[b, dt, 0, 0] = Msym.sum()
            C2[b, dt, 0, 1] = Msym[:, n0].sum()
            C2[b, dt, 1, 0] = Msym[n0, :].sum()
            C2[b, dt, 1, 1] = Msym[n0, n0]
    return C2


def main():
    ks = sorted(int(re.search(r'peram\.(\d+)\.h5', f).group(1))
                for f in glob.glob("data_%s/%s/peram.*.h5" % (ENS, NVDIR)))
    ks = [k for k in ks if k >= KMIN]
    pts = np.loadtxt("%s/pts_n%d.dat" % (GEO, L))[:12]
    n = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    theta = np.arccos(np.clip(n[:, 2], -1, 1))
    phi = np.arctan2(n[:, 1], n[:, 0])
    WM = magw(theta, phi)
    dmax = 28
    print("# magnetic point-to-point GEVP {projected, point n0=%d}  %s  n0-site test" % (N0, ENS.split("_hb")[0]), flush=True)
    Ms = []
    for i, k in enumerate(ks):
        Ms.append(config_M(k, WM, dmax))
        if (i + 1) % 50 == 0:
            print("#  %d/%d" % (i + 1, len(ks)), flush=True)
    Ms = np.array(Ms)
    nb = len(ks) // BINCFG
    Mbin = Ms[:nb * BINCFG].reshape(nb, BINCFG, dmax, 12, 12).mean(1)
    print("# %d cfg -> %d bins (bin%d), kmin=%d" % (len(ks), nb, BINCFG, KMIN), flush=True)

    C2 = build_C2(Mbin, N0)
    C2m = C2.mean(0)
    print("# C2(t0=%d) =\n%s" % (T0G, np.array2string(C2m[T0G], precision=5)), flush=True)
    print("#   C_11/C_00 at t0 = %.4f (partner overlap fraction; ~0 = decoupled)" % (C2m[T0G, 1, 1] / C2m[T0G, 0, 0]), flush=True)

    # central + jackknife effmasses: GEVP ground, projected-only C00, point-to-point C11
    emg_c = gevp2_ground(C2m, T0G)
    emg_j = np.array([gevp2_ground(np.delete(C2, b, 0).mean(0), T0G) for b in range(nb)])
    emg_e = np.sqrt((nb - 1) * np.nanmean((emg_j - np.nanmean(emg_j, 0)) ** 2, 0))

    C00 = C2m[:, 0, 0]
    C00j = np.array([np.delete(C2, b, 0).mean(0)[:, 0, 0] for b in range(nb)])
    em00_c = scalar_effmass(C00)
    em00_j = np.array([scalar_effmass(C00j[b]) for b in range(nb)])
    em00_e = np.sqrt((nb - 1) * np.nanmean((em00_j - np.nanmean(em00_j, 0)) ** 2, 0))

    C11 = C2m[:, 1, 1]
    C11j = np.array([np.delete(C2, b, 0).mean(0)[:, 1, 1] for b in range(nb)])
    em11_c = scalar_effmass(C11)
    em11_j = np.array([scalar_effmass(C11j[b]) for b in range(nb)])
    em11_e = np.sqrt((nb - 1) * np.nanmean((em11_j - np.nanmean(em11_j, 0)) ** 2, 0))

    print("#  dt |  GEVP ground        C00 projected       C11 point-to-point", flush=True)
    for t in range(1, dmax - 1):
        print("#  %2d | %7.4f(%.4f)   %7.4f(%.4f)   %7.4f(%.4f)"
              % (t + 1, emg_c[t], emg_e[t], em00_c[t], em00_e[t], em11_c[t], em11_e[t]), flush=True)

    np.savez("final/analysis_axial/mag_point2point_gevp_claude.npz",
             Mbin=Mbin, C2=C2, emg=emg_c, emg_err=emg_e, em00=em00_c, em00_err=em00_e,
             em11=em11_c, em11_err=em11_e, n0=N0)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dd = np.arange(2, dmax)
    plt.figure(figsize=(7.6, 5.0))
    sl = slice(1, dmax - 1)
    plt.errorbar(dd, emg_c[sl], yerr=emg_e[sl], marker="o", ms=4, capsize=2, lw=0.9, color="tab:red", label="GEVP ground {proj, point}")
    plt.errorbar(dd, em00_c[sl], yerr=em00_e[sl], marker="s", ms=3.5, capsize=1.5, lw=0.7, color="gray", alpha=0.8, label="C00 projected (single op)")
    plt.errorbar(dd, em11_c[sl], yerr=em11_e[sl], marker="^", ms=3.5, capsize=1.5, lw=0.7, color="tab:blue", alpha=0.7, label="C11 point-to-point n0")
    plt.axhline(0.3606, color="tab:red", ls=":", lw=0.8)
    plt.text(dmax - 8, 0.3606 + 0.003, "single-op magnetic 0.361", fontsize=7, color="tab:red")
    plt.ylim(0.25, 0.55)
    plt.xlabel("dt")
    plt.ylabel("$a_t m$")
    plt.title("Magnetic point-to-point GEVP (single point n0=%d)  L1 Nf2 g1.0" % N0)
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("final/analysis_axial/mag_point2point_gevp_claude.png", dpi=150)
    print("# wrote npz + png", flush=True)


if __name__ == "__main__":
    main()
