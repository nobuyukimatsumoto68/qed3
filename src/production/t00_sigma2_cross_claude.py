#!/usr/bin/env python3
# t00_sigma2_cross_claude.py
#   <T_00(t) sigma^2_P+(0)> cross-correlator: does the stress tensor (Delta=3) mix with our
#   parity-even two-meson sigma^2 channel (Delta=4)?  RESULT: = 0 EXACTLY (machine precision, free AND
#   interacting per-config).  MECHANISM (corrected -- NOT sigma3-herm; 2+1D 2-comp = no chirality): the r=0
#   e.sigma kernel W is ANTI-HERMITIAN, so O_H = eta^dag W xi + xi^dag W^dag eta has phi_WH = -phi_W; the
#   triangle is LINEAR in the sink vertex, so the two h.c. halves cancel (Tri[phi_W]+Tri[-phi_W]=0).  In
#   <T_00 T_00> the vertex is QUADRATIC -> they ADD.  See gw_antiherm_exclusion_mechanism_claude.md.
#
#   TOPOLOGY: sink O_H = one bilinear (2 fields) ; source sigma^2 = two bilinears (4 fields) => 3 props
#   = the triangle of diag_OA_sig2_linear_claude.py, with sink vertex Phi_W (+ Phi_WH) instead of Phi*tt.
#
#   SINK (from Fin: Stress tensor):  O_H = eta^dag W xi + xi^dag W^dag eta  (r=0 naive e.sigma hop)
#     Phi_W(a)  = V[tsrc0+a]^dag W  V[tsrc0+a]     ties in with tau   (eta^dag W xi term)
#     Phi_WH(a) = V[tsrc0+a]^dag W^dag V[tsrc0+a]   ties in with tau^dag = delta - tau  (h.c. term)
#   Phi_W/Phi_WH read from the npz Fin dumps (aligned to my free-L1 window).
#   LEG (T2b-verified, NOT tau_gw):  <xi eta^dag>=tau ,  <eta xi^dag>=delta-tau .
#
#   SANITY: <T_00 T_00> = -Tr[Phi_W(b) tau(b,a) Phi_W(a) tau(a,b)] - (WH,tau) must plateau at 3/R=0.567.
#   Refs: single-meson (2,2)=0.556 ; two-meson 2m_PS=0.756 ; stress tensor Delta=3 -> 0.567.
#   Run: ENS=free NVDIR=distill_Nv24 python3 t00_sigma2_cross_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

PHIW_NPZ = os.environ.get("PHIW_NPZ", "t00_cross_claude/phiW_free_L1_claude.npz")
PLAT_LO = int(os.environ.get("PLAT_LO", "16"))
DTMAX = int(os.environ.get("DTMAX", "16"))


def effmass_plat(c):
    plat = c[PLAT_LO:].mean()
    cc = c - plat
    with np.errstate(all="ignore"):
        return np.log(cc[:-1] / cc[1:])


def triangle(PWa_list, Phi, tau, dtau, leg_sink, twin):
    # <O_H(t) sigma^2(0)> for ONE sink term with sink vertex list PWa_list[a] and sink->source leg leg_sink
    # (tau for the W term, dtau for the WH term).  Source sigma^2: vertex Phi, internal tau + contact.
    # 4 Wick classes (mult): Tri(x2), Semi(x2), sinkTad(x1), Disc(x1).  Sink at t, source at s (s<t via dt).
    Iv = np.eye(tau.shape[-1])
    tt = [tau[a, a] - 0.5 * Iv for a in range(twin)]           # source equal-time contact
    out = np.zeros((4, twin))
    for dt in range(twin):
        ns = twin - dt
        acc = np.zeros(4)
        for s in range(ns):
            t = s + dt
            PW = PWa_list[t]                                   # sink vertex at t
            Ph = Phi[s]
            tss = tt[s]
            Ls_ts = leg_sink[t, s]                            # sink(t) <- source(s)
            Ls_st = leg_sink[s, t]                            # source(s) <- sink(t)
            oaloop = np.trace(PW @ Ls_ts @ Ph @ Ls_st).real
            dS = np.trace(Ph @ tss).real
            dpS = np.trace(Ph @ tss @ Ph @ tss).real
            dW = np.trace(PW @ leg_sink[t, t]).real           # sink tadpole <O_H>(t)
            tri = -np.trace(PW @ Ls_ts @ Ph @ tss @ Ph @ Ls_st).real
            semi = oaloop * dS
            sinktad = dW * dpS
            disc = -dW * dS * dS
            acc += np.array([tri, semi, sinktad, disc])
        out[:, dt] = acc / ns
    return out


def main():
    tag = dc.ENS
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    k = dc.KS[0]
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Iv = np.eye(tau.shape[-1])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    # tau^dag = delta - tau  (GW; complete distillation basis)
    dtau = -tau.copy()
    for a in range(twin):
        dtau[a, a] = Iv - tau[a, a]

    d = np.load(PHIW_NPZ)
    PhiW = [d["phiW"][a] for a in range(twin)]
    PhiWH = [d["phiWH"][a] for a in range(twin)]
    print("# ENS=%s twin=%d Nv=%d  loaded %s" % (tag, twin, tau.shape[-1], PHIW_NPZ))

    # ---- SANITY: <T_00 T_00> (both terms tau; even-degree so tau ok) -> must plateau ~0.567 ----
    Chh = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        acc = 0.0
        for s in range(ns):
            t = s + dt
            acc += -np.trace(PhiW[t] @ tau[t, s] @ PhiW[s] @ tau[s, t]).real
            acc += -np.trace(PhiWH[t] @ tau[t, s] @ PhiWH[s] @ tau[s, t]).real
        Chh[dt] = acc / ns
    emhh = effmass_plat(Chh)
    print("\n# SANITY <T_00 T_00> effmass (must -> 3/R = 0.567):")
    for dt in range(1, 12):
        print("#  dt=%2d  m=%7.4f" % (dt, emhh[dt]))

    # ---- CROSS <T_00 sigma^2>: BOTH sink terms tie with plain tau (Fin-confirmed; h.c. baked into PhiWH,
    #      and the connected triangle's two sink legs make tau^dag=-tau signs cancel anyway) ----
    triW = triangle(PhiW, Phi, tau, dtau, tau, twin)          # W term, sink leg tau
    triWH = triangle(PhiWH, Phi, tau, dtau, tau, twin)        # WH term, sink leg tau
    mult = np.array([2.0, 2.0, 1.0, 1.0])
    totW = (mult[:, None] * triW).sum(0)
    totWH = (mult[:, None] * triWH).sum(0)
    cross = totW + totWH
    # cross-check: WH with tau^dag on the sink legs -- should agree on the connected Tri (2 legs, +1)
    triWH_dag = triangle(PhiWH, Phi, tau, dtau, dtau, twin)
    cross_dagWH = totW + (mult[:, None] * triWH_dag).sum(0)

    # connected triangle only (Tri class), summed W+WH(correct legs)
    triW_conn = 2.0 * triW[0]
    triWH_conn = 2.0 * triWH[0]
    conn = triW_conn + triWH_conn

    emc = effmass_plat(cross)
    emconn = effmass_plat(conn)
    # scale reference: |cross| vs the sigma^2 diagonal norm (two-meson channel) at each dt
    print("\n# <T_00 sigma^2> cross  (refs: single (2,2)=0.556, two-meson 2m_PS=0.756)")
    print("#  dt |   cross(full)   conn-Tri    |  m_eff(full)  m_eff(conn) |  cross(WH=tau^dag xcheck)")
    for dt in range(1, DTMAX):
        print("#  %2d | %11.3e  %11.3e |   %7.4f     %7.4f   |  %11.3e"
              % (dt, cross[dt], conn[dt],
                 emc[dt] if np.isfinite(emc[dt]) else np.nan,
                 emconn[dt] if np.isfinite(emconn[dt]) else np.nan,
                 cross_dagWH[dt]))

    # magnitude check: is the cross negligible vs a same-normalization sigma^2 two-meson amplitude?
    s2diag = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a = 0.0
        for s in range(ns):
            t = s + dt
            a += -np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t]).real   # single-meson-scale ref
        s2diag[dt] = a / ns
    print("\n# |cross|/|C_S(sigma 2pt)| at dt=5,8,10 = %.3e %.3e %.3e  (small => decoupled)"
          % (abs(cross[5]) / abs(s2diag[5]), abs(cross[8]) / abs(s2diag[8]), abs(cross[10]) / abs(s2diag[10])))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for y, lab, col in [(0.556, "single (2,2)=0.556", "tab:orange"),
                        (0.756, "two-meson 2m_PS=0.756", "tab:blue"),
                        (0.567, "stress T Delta=3 (0.567)", "tab:green")]:
        ax.axhline(y, color=col, ls="--", lw=1, alpha=0.6)
        ax.text(DTMAX * 0.45, y + 0.008, lab, color=col, fontsize=8)
    ax.plot(dts, emc[dts], color="tab:red", marker="o", ms=5, lw=1.1, label=r"$\langle T_{00}\,\sigma^2\rangle$ full")
    ax.plot(dts, emconn[dts], color="tab:purple", marker="^", ms=5, lw=1.0, ls="--", label="connected triangle")
    ax.set_ylim(0.3, 1.0)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L1  $\langle T_{00}\,\sigma^2_{P+}\rangle$ cross (mixing check)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_sigma2_cross_free_L1_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
