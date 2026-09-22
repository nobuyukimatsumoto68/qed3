# MAGNETIC (Psi, ell=1) variational GEVP with a TIME-SPLIT interpolator basis:
#   O_0 = psibar_t  sigma^a psi_t      (local)
#   O_1 = psibar_t  sigma^a psi_{t+1}  (temporally point-split)
# 2x2 correlator matrix C_ij(dt) = <O_i(sink) O_j(src)>, magnetic-VSH projected, per config, jackknife.
# GEVP C(dt) v = lam C(t0) v -> ground effmass.  Vector contraction (V=A degeneracy).
# Sink O_i field-times (a_i=psibar, b_i=psi): O_0->(t,t), O_1->(t,t+1).  Source O_j (c_j=psibar,d_j=psi):
# O_0->(t0,t0), O_1->(t0,t0+1).  f_ij(n1,n2)=-tr[sigma^a P1 sigma^b P2],
#   P1=S(n2@b_i <- n1@c_j)=AblkS(b_i,c_j).transpose(2,0,1,3);  P2=S(n1@d_j <- n2@a_i)=AblkS(d_j,a_i).transpose(0,2,1,3).
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys, glob, re
import numpy as np
ENS = os.environ.get("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
NVDIR = os.environ.get("NVDIR", "distill_Nv24")
os.environ["ENS"] = ENS; os.environ["NVDIR"] = NVDIR
NCFG = int(os.environ.get("NCFG", "0"))
BINCFG = int(os.environ.get("BINCFG", "10"))
T0G = int(os.environ.get("T0G", "3"))
import fs_gevp_point_claude as fg
GEO = "/mnt/barracuda22/qed3/qed3/geometry/data"; AT = 0.2
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]], complex)
SIG = {1: s1, 2: s2}


def magw(theta, phi):
    c = np.sqrt(3.0 / (4.0 * np.pi)); st = np.sin(theta); ct = np.cos(theta)
    grad = {0: (-c * st, np.zeros_like(st)), 1: (c * ct * np.cos(phi), -c * np.sin(phi)),
            -1: (c * ct * np.sin(phi), c * np.cos(phi))}
    WM = {}
    for m in (-1, 0, 1):
        wt, wp = grad[m]
        WM[m] = {1: -wp, 2: wt}     # magnetic: (sigma1,sigma2) weights = (-dphiY/sin, dthetaY)
    return WM


def config_matrix(k, WM, w, dmax):
    AblkS, _, twin, ns, _, _ = fg.make_config(k)
    C = np.zeros((dmax, 2, 2))
    SNK = [(0, 0), (1, 0)]   # sink = O_i^dag: local (t,t); split^dag = psibar_{t+1} psi_t -> (a=t+1,b=t)
    SRC = [(0, 0), (0, 1)]   # source = O_j: local (t0,t0); split = psibar_{t0} psi_{t0+1} -> (c=t0,d=t0+1)
    for dt in range(1, dmax + 1):
        M = np.zeros((2, 2)); nt = 0
        for t0 in range(0, twin - dt - 1):
            t = t0 + dt
            for i, (ao, bo) in enumerate(SNK):
                for j, (co, do) in enumerate(SRC):
                    P1 = AblkS(t + bo, t0 + co).transpose(2, 0, 1, 3)   # [n1,n2] S(n2@b_i <- n1@c_j)
                    P2 = AblkS(t0 + do, t + ao).transpose(0, 2, 1, 3)   # [n1,n2] S(n1@d_j <- n2@a_i)
                    fab = {(a, b): -np.einsum('pq,ijqr,rs,ijsp->ij', SIG[a], P1, SIG[b], P2, optimize=True)
                           for a in (1, 2) for b in (1, 2)}
                    val = 0.0
                    for m in (-1, 0, 1):
                        Wm = WM[m]
                        for a in (1, 2):
                            for b in (1, 2):
                                val += np.sum(np.outer(w * Wm[a], w * Wm[b]) * fab[(a, b)])
                    M[i, j] += val.real
            nt += 1
        C[dt - 1] = M / nt
    return C


def gevp_ground(C, t0):
    tmax = C.shape[0]; lam = np.full(tmax, np.nan)
    C0 = 0.5 * (C[t0] + C[t0].T)
    wv, U = np.linalg.eigh(C0)
    if wv.max() <= 0:
        return np.full(tmax - 1, np.nan)
    keep = wv > 1e-10 * wv.max(); Uk = U[:, keep] / np.sqrt(np.abs(wv[keep]))
    for t in range(tmax):
        Ct = 0.5 * (C[t] + C[t].T)
        try:
            lam[t] = np.sort(np.linalg.eigvalsh(Uk.T @ Ct @ Uk))[-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        return np.log(lam[:-1] / lam[1:])


def main():
    ks = sorted(int(re.search(r'peram\.(\d+)\.h5', f).group(1)) for f in glob.glob("data_%s/%s/peram.*.h5" % (ENS, NVDIR)))
    ks = [k for k in ks if k >= 20]
    if NCFG: ks = ks[:NCFG]
    pts = np.loadtxt("%s/pts_n1.dat" % GEO)[:12]; n = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    th = np.arccos(np.clip(n[:, 2], -1, 1)); ph = np.arctan2(n[:, 1], n[:, 0])
    w = np.ones(12); WM = magw(th, ph); dmax = 26
    print("# magnetic time-split GEVP {O_t,t ; O_t,t+1}  ENS=%s  %d cfg  t0=%d" % (ENS.split("_hb")[0], len(ks), T0G), flush=True)
    Cs = []
    for idx, k in enumerate(ks):
        Cs.append(config_matrix(k, WM, w, dmax))
        if (idx + 1) % 50 == 0: print("#  %d/%d" % (idx + 1, len(ks)), flush=True)
    Cs = np.array(Cs)
    nb = len(ks) // BINCFG
    binned = Cs[:nb * BINCFG].reshape(nb, BINCFG, dmax, 2, 2).mean(1)
    # MAGNETIC split conjugation: minus on the split-sink row (i=1) -> C(t0) positive-definite (NM).
    binned[:, :, 1, :] *= -1.0
    binned = -0.5 * (binned + np.transpose(binned, (0, 1, 3, 2)))   # symmetrize + orient ground positive
    cen = binned.mean(0)
    emc = gevp_ground(cen, T0G)
    emj = np.array([gevp_ground(np.delete(binned, b, 0).mean(0), T0G) for b in range(nb)])
    err = np.sqrt((nb - 1) * np.nanmean((emj - np.nanmean(emj, 0)) ** 2, 0))
    print("# magnetic GEVP-ground effmass (a_t*m):", flush=True)
    for t in range(3, 20):
        if np.isfinite(emc[t]): print("  dt=%2d  %.4f +- %.4f" % (t + 1, emc[t], err[t]))
    np.savez("final/analysis_axial/mag_split_gevp_claude.npz", em=emc, err=err, binned=binned)


if __name__ == "__main__":
    main()
