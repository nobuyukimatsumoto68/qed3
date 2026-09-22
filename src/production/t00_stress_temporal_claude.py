#!/usr/bin/env python3
# t00_stress_temporal_claude.py
# Direct time-displacement T_00 in the first-order eta/xi formalism, Y00-projected, BOTH terms:
#   O(t) = eta^H(t) M D_t xi(t)  +  xi^H(t) M D_t eta(t),   M = diag(w) sigma_3,  w_x = A_x Y00,
#   D_t f(t) = (1/2)[ f(t+1) - f(t-1) ]  (symmetric temporal difference on the propagator leg).
# The OLD t00_stress_claude.py had only the FIRST term (eta^H..xi, D_ov^{-1}).  Here we add the SECOND term
# (xi^H..eta) with the DH propagator <eta xi^H> = D_ov^{-dag} = -AblkS (off-diag; 1/2 I - AblkS equal-time).
# Connected (single loop; cross terms vanish; G=AblkS, Gt=-AblkS):
#   C(dt) = -(1/4) sum_{sx,sy=+-1} sx sy { Tr[M G(t+sx,t0) M G(t0+sy,t)] + Tr[M Gt(t+sx,t0) M Gt(t0+sy,t)] }.
# Prints term A / term B separately.  See t00_hamiltonian_derivation_claude.md.
# Run:  ENS=free NVDIR=distill_Nv24 LREF=1 python3 t00_stress_temporal_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as fg

NS = dc.NS
DTMAX = int(os.environ.get("DTMAX", "28"))
BINSIZE = int(os.environ.get("BINSIZE", "1"))
NCFG = int(os.environ.get("NCFG", "0"))
ROVER = 1.0 / 0.189
SIG3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


def build_M(dual):
    # M = diag(w) sigma_3 (N,N), w_x = A_x Y00  (Y00 projection; block-diagonal in site, sigma_3 in spin)
    nsite = dual.shape[0]
    N = NS * nsite
    M = np.zeros((N, N), complex)
    w = dual * dc.Y00
    for i in range(nsite):
        M[NS * i:NS * i + 2, NS * i:NS * i + 2] = w[i] * SIG3
    return M


def corr_one_config(k, M):
    AblkS, AblkSt, twin, nsite, U, tau = fg.make_config(k)
    N = NS * nsite
    Iden = np.eye(N)

    def G(a, b):
        return AblkS(a, b).reshape(N, N)

    def Gt(a, b):
        g = -AblkS(a, b).reshape(N, N)
        if a == b:
            g = 0.5 * Iden - AblkS(a, b).reshape(N, N)
        return g

    CA = np.full(DTMAX, np.nan)
    CB = np.full(DTMAX, np.nan)
    signs = (1, -1)
    for dt in range(DTMAX):
        s_lo = 1
        s_hi = twin - 2 - dt                        # need s-1>=0 and s+dt+1<=twin-1
        if s_hi < s_lo:
            continue
        accA = 0.0
        accB = 0.0
        cnt = 0
        for s in range(s_lo, s_hi + 1):
            t0 = s
            t = s + dt
            a = 0.0
            b = 0.0
            for sx in signs:
                for sy in signs:
                    a += sx * sy * np.trace(M @ G(t + sx, t0) @ M @ G(t0 + sy, t))
                    b += sx * sy * np.trace(M @ Gt(t + sx, t0) @ M @ Gt(t0 + sy, t))
            accA += -0.25 * a
            accB += -0.25 * b
            cnt += 1
        CA[dt] = (accA / cnt).real
        CB[dt] = (accB / cnt).real
    return CA, CB


def main():
    tag = dc.ENS
    dual = dc.dual_areas_from_mesh()
    M = build_M(dual)
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    if not ks:
        print("# no perambulators in %s" % dc.PERAM_DIR)
        return
    print("# ENS=%s NVDIR=%s nsite=%d ncfg=%d DTMAX=%d" % (tag, dc.NVDIR, dual.shape[0], len(ks), DTMAX))

    resA = []
    resB = []
    for k in ks:
        CA, CB = corr_one_config(k, M)
        resA.append(CA)
        resB.append(CB)
    CA = np.mean(resA, 0)
    CB = np.mean(resB, 0)
    Cfull = CA + CB

    sgn = np.sign(Cfull[2]) if np.isfinite(Cfull[2]) else 1.0
    with np.errstate(all="ignore"):
        em = np.log((sgn * Cfull)[:-1] / (sgn * Cfull)[1:])

    print("\n#  dt |   termA(D^-1)    termB(D^-H)    B/A     full C        m_eff  [T00=%.3f sigma=%.3f]"
          % (3.0 / ROVER, 2.0 / ROVER))
    for dt in range(1, DTMAX - 1):
        if not np.isfinite(Cfull[dt]):
            continue
        ratio = CB[dt] / CA[dt] if CA[dt] != 0 else np.nan
        m = em[dt] if dt < em.shape[0] else np.nan
        print("#  %2d | % .5e  % .5e  %6.3f  % .5e  %7.4f" % (dt, CA[dt], CB[dt], ratio, Cfull[dt], m))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(em.shape[0])
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    for val, lab, col in ((2.0 / ROVER, r"$\sigma=2/R$", "tab:green"),
                          (3.0 / ROVER, r"$T_{00}=3/R=0.567$", "tab:red"),
                          (4.0 / ROVER, r"$2m=4/R$", "gray")):
        ax.axhline(val, color=col, ls="--", lw=1, alpha=0.7)
        ax.text(1.2, val + 0.006, lab, fontsize=9, color=col)
    g = np.isfinite(em) & (np.abs(em) < 3)
    ax.plot(ts[g], em[g], "o-", color="tab:red", ms=5, lw=1.1, label=r"$T_{00}$ temporal (both terms, $Y_{00}$)")
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(1, DTMAX - 1)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"Time-displacement $T_{00}$ (both terms, $Y_{00}$-projected)  %s  %d cfg" % (tag, len(ks)),
                 fontsize=10.5)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_stress_temporal_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
