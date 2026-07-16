#!/usr/bin/env python3
# Per-m correlated constant fit + variance-weighted average, on a per-m GEVP jackknife dump.
#   jk dump: line1 "# nbins ntpts nstates at"; then nbins rows of ntpts*nstates (col = ti*nstates + s).
# For each fitted state s (an m-block ground), a CORRELATED constant fit over the t-range uses the
# jackknife covariance across t: fixed fit weights w = C^{-1}1/(1^T C^{-1}1), M^b = w.d^b (linear ->
# per-sample), giving M_s, sigma_s, chi2/dof. Then the states are combined by INVERSE-VARIANCE weights
# u_s = (1/sigma_s^2)/sum, applied PER SAMPLE (M^b = sum_s u_s M_s^b) and jackknifed -> correct error
# through the m-m correlations.
#   usage: fit_perm_claude.py <jkfile> <tlo> <thi> <states e.g. 0,1,2 | single s>
import sys
import numpy as np

def load(fn):
    with open(fn) as f:
        h = f.readline().split()
    nbins, ntpts, nstates, at = int(h[1]), int(h[2]), int(h[3]), float(h[4])
    a = np.loadtxt(fn).reshape(nbins, ntpts, nstates)
    return a, at

def fit_state(jk, ti_list, s):
    # DIAGONAL chi^2 const fit of state s over ti_list: keep only the diagonal (variances) in the chi^2
    # -> inverse-variance weights w_t = (1/var_t)/sum (center passes THROUGH the points; the strong
    # t-t correlations no longer drive the center outside the data range). The ERROR is still the
    # JACKKNIFE of the linear estimator M^b = w.d^b, so it remains correlation-aware. Returns
    # (Mb per-sample, M, sigma, chi2_diag, dof).
    nbins = jk.shape[0]
    d = jk[:, ti_list, s]                 # (nbins, npts)
    npts = d.shape[1]
    mean = d.mean(0)
    dev = d - mean
    var = (nbins-1.0)/nbins * np.sum(dev*dev, axis=0)   # diagonal jackknife variance per t
    w = (1.0/var); w /= w.sum()           # inverse-variance (diagonal) weights
    Mb = d @ w
    M = Mb.mean()
    sig = np.sqrt((nbins-1.0)/nbins * np.sum((Mb-M)**2))   # jackknife error (correlation-aware)
    r = mean - M
    chi2 = float(np.sum(r*r/var))         # DIAGONAL chi^2 (off-diagonal dropped)
    return Mb, M, sig, chi2, npts-1

if __name__ == "__main__":
    fn = sys.argv[1]
    tlo, thi = float(sys.argv[2]), float(sys.argv[3])
    states = [int(x) for x in sys.argv[4].split(",")]
    jk, at = load(fn)
    ntpts = jk.shape[1]
    ti_list = [ti for ti in range(ntpts) if tlo-1e-9 <= (ti+1)*at <= thi+1e-9]
    nbins = jk.shape[0]
    print(f"# {fn}: bins={nbins} fit t in [{tlo},{thi}] ({len(ti_list)} pts) states={states}")
    Mbs, sigs = [], []
    for s in states:
        Mb, M, sig, chi2, dof = fit_state(jk, ti_list, s)
        Mbs.append(Mb); sigs.append(sig)
        print(f"  state {s}: M={M:.4f}({sig:.4f})  chi2/dof={chi2:.1f}/{dof}={chi2/dof:.2f}")
    if len(states) == 1:
        sys.exit(0)
    # inverse-variance combine across states (per-sample, fixed weights), jackknifed
    Mbs = np.array(Mbs)                    # (nstates, nbins)
    u = np.array([1.0/sig**2 for sig in sigs]); u /= u.sum()
    Mc = (u[:,None] * Mbs).sum(0)          # (nbins,)
    M = Mc.mean()
    err = np.sqrt((nbins-1.0)/nbins * np.sum((Mc-M)**2))
    print(f"  ==> variance-avg over states {states}: M = {M:.4f}({err:.4f})   weights={np.round(u,3).tolist()}")
