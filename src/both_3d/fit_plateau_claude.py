#!/usr/bin/env python3
# Correlated constant (plateau) fit to the GEVP effective masses, using the jackknife covariance.
#
# The three l=1 states are EXACTLY degenerate (T_1 irrep of I_h) -> there is no meaningful continuity
# / level-ordering in t, so we fit ALL three levels AND the first nfit t-points SIMULTANEOUSLY to a
# single constant M. Points p index (ti, s); the covariance C_pq is the jackknife covariance over bins.
#
# Fit (fixed-covariance correlated constant): with 1 = ones vector,
#   w = C^{-1} 1 / (1^T C^{-1} 1)   (fit weights)
#   M = w . mean                    (best-fit constant)
#   M_b = w . d_b  (per jk sample)  -> jackknife error on M  (= 1/sqrt(1^T C^{-1} 1))
#   chi2 = (mean - M 1)^T C^{-1} (mean - M 1) ,  dof = npts - 1
#
# jk dump format (from glue_gevp_analysis_claude.cu, arg22): line1 "# nbins ntpts nstates at";
# then nbins rows of ntpts*nstates values, column = ti*nstates + s.
import sys
import numpy as np

def load_jk(fn):
    with open(fn) as f:
        hdr = f.readline().split()
    nbins = int(hdr[1])
    ntpts = int(hdr[2])
    nstates = int(hdr[3])
    at = float(hdr[4])
    a = np.loadtxt(fn)                       # (nbins, ntpts*nstates)
    a = a.reshape(nbins, ntpts, nstates)
    return a, at

def fit_const(jk, ti_list, s_list):
    # jk : (nbins, ntpts, nstates) leave-one-out jackknife samples of the effmass
    nbins = jk.shape[0]
    d = jk[:, ti_list][:, :, s_list].reshape(nbins, -1)   # (nbins, npts)
    npts = d.shape[1]
    mean = d.mean(axis=0)
    dev = d - mean
    # jackknife covariance of the mean: (nbins-1)/nbins * sum_b dev_b dev_b^T
    C = (nbins - 1.0) / nbins * (dev.T @ dev)
    if nbins <= npts + 2:
        return dict(M=np.nan, err=np.nan, chi2=np.nan, dof=npts-1, npts=npts, nbins=nbins, singular=True)
    Cinv = np.linalg.inv(C)
    one = np.ones(npts)
    denom = one @ Cinv @ one
    w = (Cinv @ one) / denom
    M = w @ mean
    Mb = d @ w
    err = np.sqrt((nbins - 1.0) / nbins * np.sum((Mb - Mb.mean())**2))
    r = mean - M * one
    chi2 = float(r @ Cinv @ r)
    dof = npts - 1
    return dict(M=float(M), err=float(err), chi2=chi2, dof=dof, npts=npts, nbins=nbins, singular=False)

if __name__ == "__main__":
    # usage: fit_plateau_claude.py <jkfile> <nfit> <states e.g. 0,1,2>
    fn = sys.argv[1]
    nfit = int(sys.argv[2])
    s_list = [int(x) for x in sys.argv[3].split(",")]
    jk, at = load_jk(fn)
    ntpts = jk.shape[1]
    ti_list = list(range(min(nfit, ntpts)))
    r = fit_const(jk, ti_list, s_list)
    tlo = (ti_list[0]+1)*at
    thi = (ti_list[-1]+1)*at
    if r["singular"]:
        print(f"{fn}: SINGULAR (nbins={r['nbins']} <= npts+2={r['npts']+2}); need more bins")
    else:
        print(f"{fn}: M={r['M']:.4f} +/- {r['err']:.4f}  chi2/dof={r['chi2']:.1f}/{r['dof']}={r['chi2']/r['dof']:.2f}  "
              f"(t=[{tlo:.1f},{thi:.1f}], states={s_list}, nbins={r['nbins']})")
