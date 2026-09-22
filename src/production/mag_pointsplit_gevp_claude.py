# MAGNETIC (Psi, ell=1) variational GEVP with a POINT-SPLIT interpolator basis:
#   O_0 = psibar(n) sigma psi(n)                 (local)
#   O_1 = sum_{n' in NN(n)} psibar(n) sigma psi(n')   (nearest-neighbor spatial split)
# These are DIFFERENT single-particle operators (site-diagonal vs site-off-diagonal) -> independent
# magnetic-T1 copies -> a genuine 2x2 GEVP on L1.  Reconstruct all-to-all overlap prop from perams.
# Correlator (vector; V=A degeneracy):  smear the psi (sink-site) leg of each propagator by K^i = {I, A}:
#   f_ij^{ab}(n1,n2) = -tr[sigma^a (K^i P1)[n2,n1] sigma^b (K^j P2)[n1,n2]],  P1=AblkS(t,t0), P2=AblkS(t0,t)
# then project both VERTEX sites (n1,n2) with the magnetic VSH weight.  GEVP with magnetic-conj sign fix.
import os
os.environ.setdefault("OMP_NUM_THREADS", "4"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys, glob, re
import numpy as np
ENS = os.environ.get("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
NVDIR = os.environ.get("NVDIR", "distill_Nv24"); os.environ["ENS"] = ENS; os.environ["NVDIR"] = NVDIR
NCFG = int(os.environ.get("NCFG", "0")); BINCFG = int(os.environ.get("BINCFG", "10")); T0G = int(os.environ.get("T0G", "3"))
import fs_gevp_point_claude as fg
GEO = "/mnt/barracuda22/qed3/qed3/geometry/data"; AT = 0.2
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]], complex); SIG = {1: s1, 2: s2}


def magw(theta, phi):
    c = np.sqrt(3.0 / (4.0 * np.pi)); st = np.sin(theta); ct = np.cos(theta)
    grad = {0: (-c * st, 0 * st), 1: (c * ct * np.cos(phi), -c * np.sin(phi)), -1: (c * ct * np.sin(phi), c * np.cos(phi))}
    return {m: {1: -grad[m][1], 2: grad[m][0]} for m in (-1, 0, 1)}   # magnetic (sigma1,sigma2) weights


def adjacency(nvec):
    ns = nvec.shape[0]; dot = nvec @ nvec.T; A = np.zeros((ns, ns))
    for i in range(ns):
        nb = np.argsort(dot[i])[::-1][1:6]     # 5 nearest neighbours (icosahedron)
        A[i, nb] = 1.0
    return A


def config_matrix(k, WM, w, Ks, dmax):
    AblkS, _, twin, ns, _, _ = fg.make_config(k)
    nb = len(Ks); C = np.zeros((dmax, nb, nb))
    for dt in range(1, dmax + 1):
        M = np.zeros((nb, nb)); ncnt = 0
        for t0 in range(0, twin - dt):
            t = t0 + dt
            P1 = AblkS(t, t0)      # [nsink,alpha,nsrc,beta]  S(nsink@t <- nsrc@t0)
            P2 = AblkS(t0, t)      # [nsink,alpha,nsrc,beta]  S(nsink@t0 <- nsrc@t)
            for i in range(nb):
                P1w = np.einsum('mn,najb->majb', Ks[i], P1)     # smear sink of P1 by K^i
                for j in range(nb):
                    P2w = np.einsum('mn,najb->majb', Ks[j], P2)
                    # f_ij[n1,n2] = -tr[sig^a P1w[n2,n1] sig^b P2w[n1,n2]]
                    P1arr = np.transpose(P1w, (2, 0, 1, 3))    # [nsrc=n1, nsink=n2, a, b]
                    P2arr = np.transpose(P2w, (0, 2, 1, 3))    # [nsink=n1, nsrc=n2, a, b]
                    val = 0.0
                    for a in (1, 2):
                        for b in (1, 2):
                            f = -np.einsum('pq,ijqr,rs,ijsp->ij', SIG[a], P1arr, SIG[b], P2arr, optimize=True)
                            for m in (-1, 0, 1):
                                val += np.sum(np.outer(w * WM[m][a], w * WM[m][b]) * f)
                    M[i, j] += val.real
            ncnt += 1
        C[dt - 1] = M / ncnt
    return C


def gevp_ground(C, t0):
    tmax = C.shape[0]; lam = np.full(tmax, np.nan)
    C0 = 0.5 * (C[t0] + C[t0].T); wv, U = np.linalg.eigh(C0)
    if wv.max() <= 0: return np.full(tmax - 1, np.nan)
    keep = wv > 1e-10 * wv.max(); Uk = U[:, keep] / np.sqrt(np.abs(wv[keep]))
    for t in range(tmax):
        try: lam[t] = np.sort(np.linalg.eigvalsh(Uk.T @ (0.5 * (C[t] + C[t].T)) @ Uk))[-1]
        except Exception: pass
    with np.errstate(all="ignore"): return np.log(lam[:-1] / lam[1:])


def main():
    ks = sorted(int(re.search(r'peram\.(\d+)\.h5', f).group(1)) for f in glob.glob("data_%s/%s/peram.*.h5" % (ENS, NVDIR)))
    ks = [k for k in ks if k >= 20]
    if NCFG: ks = ks[:NCFG]
    pts = np.loadtxt("%s/pts_n1.dat" % GEO)[:12]; nvec = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    th = np.arccos(np.clip(nvec[:, 2], -1, 1)); ph = np.arctan2(nvec[:, 1], nvec[:, 0]); w = np.ones(12)
    WM = magw(th, ph); A = adjacency(nvec); Ks = [np.eye(12), A]; dmax = 22
    print("# magnetic POINT-SPLIT GEVP {local, NN-split}  %s  %d cfg  t0=%d" % (ENS.split("_hb")[0], len(ks), T0G), flush=True)
    Cs = []
    for idx, k in enumerate(ks):
        Cs.append(config_matrix(k, WM, w, Ks, dmax))
        if (idx + 1) % 50 == 0: print("#  %d/%d" % (idx + 1, len(ks)), flush=True)
    Cs = np.array(Cs); nbb = len(ks) // BINCFG
    binned = Cs[:nbb * BINCFG].reshape(nbb, BINCFG, dmax, 2, 2).mean(1)
    cen = binned.mean(0)
    print("# raw C(t0=%d) =\n%s  eig %s" % (T0G, np.array2string(cen[T0G], precision=4), np.linalg.eigvalsh(0.5 * (cen[T0G] + cen[T0G].T))), flush=True)
    # orient: overall sign so ground positive; if diagonals opposite-sign, apply magnetic-conj minus to split row/col
    s = np.array([1.0, np.sign(cen[T0G, 0, 0] * cen[T0G, 1, 1])])   # +1 if same-sign diag, else flip split
    B = binned.copy(); B[:, :, 1, :] *= s[1]
    B = np.sign(cen[T0G, 0, 0]) * 0.5 * (B + np.transpose(B, (0, 1, 3, 2)))   # symmetrize + orient ground +
    cen2 = B.mean(0)
    print("# oriented C(t0) eig %s (split-sign=%.0f)" % (np.linalg.eigvalsh(cen2[T0G]), s[1]), flush=True)
    emc = gevp_ground(cen2, T0G)
    emj = np.array([gevp_ground(np.delete(B, b, 0).mean(0), T0G) for b in range(nbb)])
    err = np.sqrt((nbb - 1) * np.nanmean((emj - np.nanmean(emj, 0)) ** 2, 0))
    print("# magnetic point-split GEVP-ground (a_t*m):", flush=True)
    for t in range(3, 18):
        if np.isfinite(emc[t]): print("  dt=%2d  %.4f +- %.4f" % (t + 1, emc[t], err[t]))
    np.savez("final/analysis_axial/mag_pointsplit_gevp_claude.npz", em=emc, err=err, binned=binned)


if __name__ == "__main__":
    main()
