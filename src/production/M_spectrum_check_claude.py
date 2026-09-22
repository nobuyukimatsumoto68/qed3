#!/usr/bin/env python3
# M_spectrum_check_claude.py -- why does Im Tr(M^3) vanish in FREE but not INTERACTING?  M = D_ov^{-1} - 1/2 is
#   anti-hermitian (GW) so eig(M) = i mu_k, mu_k real, and Re Tr(M^3)=0 ALWAYS.  Im Tr(M^3) = -Tr(H^3) = -sum mu_k^3
#   (H = -i M Hermitian).  This vanishes ONLY if the spectrum {mu_k} is symmetric under mu -> -mu (an extra
#   spectral-reflection / sigma3-hermiticity Gamma M Gamma^{-1} = -M).  FREE has it, the gauge interaction BREAKS it.
#   Uses the equal-time block tt = tau(0,0) - 1/2 I_Nv (anti-hermitian) as the finite illustration.
#   Reports: sum mu, sum mu^3, and the reflection-asymmetry of the sorted spectrum |mu_k + mu_{Nv-1-k}|.
#   Run: ENS=.. LREF=.. NVDIR=.. python3 M_spectrum_check_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

KCFG = int(os.environ.get("KCFG", str(dc.KS[0])))


def main():
    V, tau, taugw, tsrc0, twin = dc.load_peram(KCFG)
    nvk = tau.shape[-1]
    tt = tau[0, 0] - 0.5 * np.eye(nvk)
    ah = np.linalg.norm(tt + tt.conj().T) / np.linalg.norm(tt)
    ev = np.linalg.eigvals(tt)                 # ~ i mu_k
    mu = np.sort(ev.imag)                        # real spectrum of H = -i M
    remax = np.max(np.abs(ev.real))             # should be ~0 (anti-herm -> purely imaginary eig)
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s L=%d k=%d Nv=%d  M=tau(0,0)-1/2 (anti-herm: ||M+M^dag||/||M||=%.2e, max|Re eig|=%.2e)"
          % (tag, dc.L, KCFG, nvk, ah, remax))
    s1 = mu.sum()
    s3 = (mu ** 3).sum()
    n3 = np.abs(mu) ** 3
    print("#   sum mu   = %+.6e   (0 if reflection-symmetric)" % s1)
    print("#   sum mu^3 = %+.6e   (= -Im Tr(M^3); 0 <=> Im Tr(M^3)=0)   |sum mu^3|/sum|mu|^3 = %.3e"
          % (s3, np.abs(s3) / n3.sum()))
    # reflection asymmetry: pair k-th smallest with k-th largest; |mu_k + mu_{N-1-k}| ~ 0 if mu -> -mu symmetric
    asym = np.abs(mu + mu[::-1])
    print("#   reflection mu -> -mu:  max|mu_k + mu_{N-1-k}| = %.3e   mean = %.3e   (0 => symmetric spectrum)"
          % (asym.max(), asym.mean()))
    print("#   Tr(M^3) = %s  (Re should be ~0; Im = -sum mu^3)" % np.array2string(np.trace(tt @ tt @ tt), precision=4))


if __name__ == "__main__":
    main()
