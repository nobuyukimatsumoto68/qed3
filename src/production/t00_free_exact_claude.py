#!/usr/bin/env python3
# t00_free_exact_claude.py
# EXACT free-field overlap propagator on S^2 x R (no GPU, no distillation) for any refinement L
# (plan: t00_free_exact_impl_plan_claude.md).  Conventions copied from includes/dirac_ext.h + overlap.h:
#   D_W[(s,ix),(s,iy)] = 0.5 kappa_il (-r + gamma(ix,iy)) Omega(ix,iy)            (spatial hop; free phase=1)
#   D_W[(s+1,ix),(s,ix)] = 0.5 signP kappa_t (-r - sigma3),  D_W[(s-1,ix),(s,ix)] = 0.5 signM kappa_t (-r + sigma3)
#   diag = sum_nns 0.5 r kappa + r kappa_t + M5,  kappa_t = dual_area/mean_ell/at,  antiperiodic in t (signP/M=-1 at edge)
#   D_ov = 1 + D_W (D_W^dag D_W)^{-1/2}   (massless; Zolotarev -> exact sign here)
# Temporal momentum space: D(p) = diag + spatial + tmpP e^{-ip} + tmpM e^{+ip}, p=(2n+1)pi/Nt ;
#   G(s,s') = (1/Nt) sum_p e^{ip(s-s')} D_ov(p)^{-1}.
# Observables: C_sigma(dt) = -Tr[S G(dt) S G(-dt)] (S=diag(A_x Y00)) ; C_ell(dt) = -2 Re Tr[W_m G(dt) W_m G(-dt)]_m-avg,
#   W_m = W^{r=0} * Y_lm(link midpoint)  (same as t00_ell2_claude.py) .
# Run: LREF=3 python3 t00_free_exact_claude.py     (VALIDATE=1 compares G with the peram at L1/L2)

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
LREF = int(os.environ.get("LREF", "1"))
os.environ["LREF"] = str(LREF)
os.environ.setdefault("NVDIR", "distill_Nv24" if LREF == 1 else "distill_Nv84")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import geom_hopping_claude as gh
import t00_stress_ham_interacting_claude as ti
import t00_ell2_claude as e2

NS = 2
NT = int(os.environ.get("NT", "128"))
AT = float(os.environ.get("AT", "0.2"))
M5 = float(os.environ.get("M5", "-1.0"))
RW = 1.0
DTMAX = int(os.environ.get("DTMAX", "40"))
ELLS = [int(x) for x in os.environ.get("ELLS", "0,2").split(",")]
VALIDATE = int(os.environ.get("VALIDATE", "0"))
GEOM = "../../geometry/data/"
ROVER = 1.0 / 0.189
S0 = np.eye(2)
S3 = np.diag([1.0, -1.0])


def load_sites(path):
    mean_ell = None
    areas = {}
    for ln in open(path):
        if ln.startswith("# mean_ell="):
            mean_ell = float(ln.split("=")[1])
        elif not ln.startswith("#"):
            a, b = ln.split()
            areas[int(a)] = float(b)
    return np.array([areas[i] for i in range(len(areas))]), mean_ell


def build_Gt(om, alpha, nns, nsite, tab, dual, mean_ell, shifts):
    N = NS * nsite
    n_links = len(tab) // 2
    Wsp = ti.build_W_gauge(om, alpha, nns, nsite, tab, np.zeros(n_links), r=RW)   # hop + 0.5 r sum kappa diag
    kt = dual / mean_ell / AT
    Kt = np.kron(np.diag(kt), S0)
    diag = RW * Kt + M5 * np.eye(N)
    P = 0.5 * np.kron(np.diag(kt), (-RW * S0 - S3))       # tmpP (signs handled by antiperiodic momenta)
    Mm = 0.5 * np.kron(np.diag(kt), (-RW * S0 + S3))      # tmpM
    D0 = Wsp + diag
    ps = (2.0 * np.arange(NT) + 1.0) * np.pi / NT
    G = {d: np.zeros((N, N), complex) for d in shifts}
    for p in ps:
        Dp = D0 + P * np.exp(-1j * p) + Mm * np.exp(1j * p)
        H = Dp.conj().T @ Dp
        w, U = np.linalg.eigh(H)
        Hinvhalf = (U * (1.0 / np.sqrt(w))) @ U.conj().T
        Dov = np.eye(N) + Dp @ Hinvhalf
        Ginv = np.linalg.inv(Dov)
        for d in shifts:
            G[d] += np.exp(1j * p * d) * Ginv / NT
    return G


def main():
    om, alpha, nns, nsite = gh.build(GEOM, LREF)
    tab = ti.load_link_table("primal_links_n%d_claude.dat" % LREF)
    dual, mean_ell = load_sites("primal_sites_n%d_claude.dat" % LREF)
    pts = dc.load_vec3(GEOM + "pts_n%d.dat" % LREF)
    N = NS * nsite
    shifts = list(range(-DTMAX, DTMAX + 1))
    cache = "t00_ham_cache_claude/free_exact_G_n%d_nt%d_at%g_claude.npz" % (LREF, NT, AT)
    if os.path.exists(cache):
        z = np.load(cache)
        G = {int(d): z["G%d" % d] for d in shifts}
    else:
        G = build_Gt(om, alpha, nns, nsite, tab, dual, mean_ell, shifts)
        np.savez(cache, **{"G%d" % d: G[d] for d in shifts})
    print("# L=%d nsite=%d N=%d Nt=%d at=%g M5=%g mean_ell=%.6f  GW check |G0+G0^H-1|=%.2e"
          % (LREF, nsite, N, NT, AT, M5, mean_ell, np.abs(G[0] + G[0].conj().T - np.eye(N)).max()))

    if VALIDATE:
        V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
        worst = 0.0
        for (a, b) in [(0, 0), (5, 0), (0, 5), (12, 3), (20, 20), (31, 0)]:
            Gp = V[tsrc0 + a].T @ tau[a, b] @ V[tsrc0 + b].conj()
            dif = np.abs(Gp - G[a - b]).max()
            worst = max(worst, dif)
            print("#   validate G(%2d,%2d): max|peram - exact| = %.2e  (|G|max=%.2e)" % (a, b, dif, np.abs(Gp).max()))
        print("# VALIDATE worst = %.2e" % worst)

    Sg = np.kron(np.diag(dc.dual_areas_from_mesh() * dc.Y00), S0)
    C = {}
    Cs = np.zeros(DTMAX + 1)
    for dt in range(1, DTMAX + 1):
        Cs[dt] = -np.trace(Sg @ G[dt] @ Sg @ G[-dt]).real
    C["sigma"] = Cs
    W0 = ti.build_W_gauge(om, alpha, nns, nsite, tab, np.zeros(len(tab) // 2), r=0.0)
    for ell in ELLS:
        Ym = e2.link_weight_matrices(pts, nns, nsite, ell)
        Wm = [W0 * np.kron(Ym[m], np.ones((NS, NS))) for m in range(Ym.shape[0])]
        Cl = np.zeros(DTMAX + 1)
        for dt in range(1, DTMAX + 1):
            acc = 0.0
            for m in range(len(Wm)):
                acc += -2.0 * np.trace(Wm[m] @ G[dt] @ Wm[m] @ G[-dt]).real
            Cl[dt] = acc / len(Wm)
        C["ell%d" % ell] = Cl
    keys = ["sigma"] + ["ell%d" % ell for ell in ELLS]
    with np.errstate(all="ignore"):
        em = {k: np.log(C[k][1:-1] / C[k][2:]) for k in keys}
    print("# refs 2/R=%.4f 3/R=%.4f 4/R=%.4f" % (2 / ROVER, 3 / ROVER, 4 / ROVER))
    print("#  dt | " + " | ".join("%-8s m_eff" % k for k in keys))
    for dt in range(1, DTMAX - 1):
        print("#  %2d | " % dt + " | ".join("%14.5f" % em[k][dt - 1] for k in keys))
    np.savez("t00_ham_cache_claude/free_exact_corr_n%d_nt%d_at%g_claude.npz" % (LREF, NT, AT), **C)


if __name__ == "__main__":
    main()
