# Time-displaced GEVP (Prony) for the MAGNETIC (Psi) channel, on the saved binned correlator.
# Build the Hankel matrix H(t)_{ij} = g(t + off_i + off_j) from time-displaced copies of the single
# magnetic operator, solve the GEVP H(t) v = lam H(t0) v, take the GROUND generalized eigenvalue,
# effmass = log(lam0(t)/lam0(t+1)).  Jackknife over config bins.  Scan basis size / t0 to reduce the
# excited-state contamination the single-operator GPOF leaves.
import numpy as np

AT = 0.2
NPZ = "final/analysis_axial/interacting_vsh_axial_L1_nf2g1_claude.npz"


def hankel_gevp_ground(g, offs, t0):
    # g: 1D correlator (dt index 0 = dt1). offs: displacement list. returns ground effmass curve.
    N = len(offs)
    tmax = len(g)
    lam0 = np.full(tmax, np.nan)

    def Hmat(t):
        M = np.empty((N, N))
        for i in range(N):
            for j in range(N):
                idx = t + offs[i] + offs[j]
                M[i, j] = g[idx] if 0 <= idx < tmax else np.nan
        return 0.5 * (M + M.T)
    C0 = Hmat(t0)
    if not np.all(np.isfinite(C0)):
        return lam0
    # symmetric GEVP via Cholesky of a lightly-regularized C0
    w, U = np.linalg.eigh(C0)
    keep = w > 1e-10 * w.max()
    Uk = U[:, keep] / np.sqrt(w[keep])
    for t in range(tmax):
        Ct = Hmat(t)
        if not np.all(np.isfinite(Ct)):
            continue
        try:
            ev = np.linalg.eigvalsh(Uk.T @ Ct @ Uk)
            lam0[t] = np.sort(ev)[-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        em = np.log(lam0[:-1] / lam0[1:])
    return em


def jk_curve(binned, offs, t0):
    nb = binned.shape[0]
    cen = binned.mean(0)
    sgn = np.sign(cen[3])
    emc = hankel_gevp_ground(sgn * cen, offs, t0)
    emj = np.array([hankel_gevp_ground(sgn * np.delete(binned, b, 0).mean(0), offs, t0) for b in range(nb)])
    err = np.sqrt((nb - 1) * np.nanmean((emj - np.nanmean(emj, 0)) ** 2, 0))
    return emc, err


def main():
    d = np.load(NPZ)
    gb = d["MA_binned"]     # (nb, ndt) magnetic
    print("# magnetic time-displaced GEVP scan (a_t*m); single-op GPOF gave ~0.361(5)")
    for offs, t0 in [([0, 1], 2), ([0, 2], 2), ([0, 1, 2], 2), ([0, 2, 4], 3), ([0, 1, 2, 3], 2), ([0, 3, 6], 4)]:
        em, err = jk_curve(gb, offs, t0)
        seg = "  ".join("%.3f(%d)" % (em[t], round(err[t] * 1e3)) for t in range(3, 15) if np.isfinite(em[t]))
        print("off=%-12s t0=%d:  dt4..: %s" % (str(offs), t0, seg))


if __name__ == "__main__":
    main()
