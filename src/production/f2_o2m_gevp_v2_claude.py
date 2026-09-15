#!/usr/bin/env python3
# f2_o2m_gevp_v2_claude.py
#   Coupled 2x2 {F^2, O_2m} GEVP -- the Chester-Pufu 0++ mixing of the glueball F^2 with the antipodal
#   two-meson sigma^2 interpolator O_2m(s)=sum_x A_x sigma(x,s)sigma(P(x),s).  Same conventions as the
#   glueball GEVP glue_gevp_analysis_claude.cu: vacuum subtraction C_ij - <O_i><O_j>, moving inv_sqrt_sym
#   metric, Delta_eff = -log(lambda)/(dt*at) (PHYSICAL, at=0.2).  F^2 diagonal reproduces a_t m_F2=2.88.
#
#   Correlators (all translation-averaged; F^2 = full Nt, O_2m = window twin=32):
#     C_FF(t) = <O_F(s+t) O_F(s)>_s            (F^2 shape op OP_F, folded by Nt periodicity)
#     C_Fs(t) = <O_F(s+/-t) L_2m(s)>_s          (folded; F^2 global)
#     C_ss(t) = <L_2m(s+t) L_2m(s)>_s           (O_2m antipodal self-loop; window-bound, t<=twin-1)
#   L_2m(s) = - sum_x A_x Tr[P(x,s;P(x),s) P(P(x),s;x,s)]  (equal-time antipodal loop, improved P=AblkS).
#
#   Run: ENS=Nf2_gsq1.000000...L1 NVDIR=distill_Nv24 BINSIZE=10 T0=... python3 f2_o2m_gevp_v2_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import h5py
import distill_contract_claude as dc
import fs_gevp_point_claude as G

Nt = 128
OP_F = int(os.environ.get("OP_F", "0"))
AT = float(os.environ.get("AT", "0.2"))
TMAX = int(os.environ.get("TMAX", "24"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
DT = int(os.environ.get("DT", "1"))              # effmass step (timeslices)
RTOL = float(os.environ.get("RTOL", "1e-8"))
EVEC_T = int(os.environ.get("EVEC_T", "3"))      # t at which to report the mixing eigenvectors


def gluedir():
    return "data_" + dc.ENS


def load_of(k):
    with h5py.File("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k), "r") as g:
        return np.array(g["O"])[OP_F].astype(float)


def o2m_loop(k, dual, Pmap):
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(k)
    idx = np.arange(nsite)
    V, _, _, tsrc0, _ = dc.load_peram(k)
    L = np.full(twin, np.nan)
    for s in range(twin):
        Pf = AblkS(s, s)
        B1 = Pf[idx, :, Pmap, :]
        B2 = Pf[Pmap, :, idx, :]
        L[s] = -(dual * np.einsum('xab,xba->x', B1, B2)).sum().real
    return L, tsrc0, twin


def inv_sqrt_sym(M, rtol):
    M = 0.5 * (M + M.T)
    w, V = np.linalg.eigh(M)
    wmax = w.max()
    inv = np.where((w > rtol * wmax) & (w > 0), 1.0 / np.sqrt(np.abs(w)), 0.0)
    return (V * inv) @ V.T


def gevp_states(C, t, dt, at):
    # C: (2,2,T) ; moving-metric GEVP effmass -log(lambda)/(dt*at), states sorted lightest->heaviest
    A = 0.5 * (C[:, :, t] + C[:, :, t].T)
    B = 0.5 * (C[:, :, t + dt] + C[:, :, t + dt].T)
    si = inv_sqrt_sym(A, RTOL)
    M = si @ B @ si
    w, vv = np.linalg.eigh(0.5 * (M + M.T))
    order = np.argsort(w)[::-1]          # descending lambda -> ascending effmass
    lam = w[order]
    em = -np.log(np.abs(lam)) / (dt * at)
    # eigenvectors back in the operator basis (v = si @ eigvec)
    evec = si @ vv[:, order]
    return em, lam, evec


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh().astype(float)
    Pmap = G.antipodal_map()
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k))]
    ncfg = len(ks)
    print("# ENS=%s cfg=%d OP_F=%d BINSIZE=%d  coupled {F^2,O_2m} GEVP" % (tag, ncfg, OP_F, BINSIZE))

    # per-config raw correlators + one-points
    cff = np.zeros((ncfg, TMAX + DT))
    cfs = np.zeros((ncfg, TMAX + DT))
    css = np.zeros((ncfg, TMAX + DT))
    ofm = np.zeros(ncfg)
    lm = np.zeros(ncfg)
    for ic, k in enumerate(ks):
        of = load_of(k)
        L, tsrc0, twin = o2m_loop(k, dual, Pmap)
        ofw = of[(tsrc0 + np.arange(twin)) % Nt]        # F^2 on the SAME window as O_2m (consistent sampling)
        ofm[ic] = ofw.mean()
        lm[ic] = L.mean()
        # ALL correlators on the window-triangular sampling s, s+t in [0,twin); cross symmetrized
        for t in range(TMAX + DT):
            if t < twin:
                a = np.arange(twin - t)
                cff[ic, t] = np.mean(ofw[a + t] * ofw[a])                           # F^2 auto (window)
                css[ic, t] = np.mean(L[a + t] * L[a])                               # O_2m auto (window)
                cfs[ic, t] = 0.5 * (np.mean(ofw[a + t] * L[a]) + np.mean(L[a + t] * ofw[a]))  # symmetrized cross
            else:
                cff[ic, t] = css[ic, t] = cfs[ic, t] = np.nan

    # bin then jackknife
    nb = ncfg // BINSIZE
    def binned(a):
        return np.array([a[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    bff, bfs, bss, bof, blm = binned(cff), binned(cfs), binned(css), binned(ofm), binned(lm)

    TREF = int(os.environ.get("TREF", "2"))          # timeslice to normalize the operators (fix scale disparity)

    def build_C_raw(sel):
        # vacuum-subtracted 2x2 correlator over selected bins (RAW operator normalization)
        ff = bff[sel].mean(0) - bof[sel].mean() ** 2
        fs = bfs[sel].mean(0) - bof[sel].mean() * blm[sel].mean()
        ss = bss[sel].mean(0) - blm[sel].mean() ** 2
        C = np.zeros((2, 2, TMAX + DT))
        C[0, 0], C[0, 1], C[1, 0], C[1, 1] = ff, fs, fs, ss
        return C

    allsel = np.arange(nb)
    Craw = build_C_raw(allsel)
    # fixed operator normalization from the FULL sample diagonal at TREF (GEVP-invariant; fixes conditioning)
    dnorm = np.sqrt(np.abs([Craw[0, 0, TREF], Craw[1, 1, TREF]]))
    NRM = np.outer(dnorm, dnorm)

    def build_C(sel):
        C = build_C_raw(sel)
        return C / NRM[:, :, None]

    # DIAGNOSTIC: raw diagonals + correlation coefficient rho(t) = C_Fs/sqrt(C_FF C_ss)
    print("\n# DIAG  t(phys) |   C_FF        C_ss         C_Fs        rho=C_Fs/sqrt(C_FF C_ss)")
    for t in range(0, min(TMAX, 12)):
        ff, ss, fs = Craw[0, 0, t], Craw[1, 1, t], Craw[0, 1, t]
        rho = fs / np.sqrt(ff * ss) if (ff > 0 and ss > 0) else np.nan
        print("#  %5.2f  | %11.3e %11.3e %11.3e   %7.3f" % (t * AT, ff, ss, fs, rho))

    Cc = build_C(allsel)
    # central + jackknife effmass per state
    em_c = np.full((TMAX, 2), np.nan)
    for t in range(1, TMAX):
        try:
            em_c[t], _, _ = gevp_states(Cc, t, DT, AT)
        except Exception:
            pass
    ems = []
    for j in range(nb):
        sel = np.delete(allsel, j)
        Cj = build_C(sel)
        e = np.full((TMAX, 2), np.nan)
        for t in range(1, TMAX):
            try:
                e[t], _, _ = gevp_states(Cj, t, DT, AT)
            except Exception:
                pass
        ems.append(e)
    ems = np.array(ems)
    em_e = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, 0))

    print("\n#  t(phys) |  light state (err)   heavy state (err)     [F2 glue=2.88, m_PS=1.6, 2m_PS=3.2]")
    for t in range(1, TMAX):
        if np.isfinite(em_c[t, 0]):
            print("#  %5.2f  |  %7.3f(%.3f)      %7.3f(%.3f)"
                  % (t * AT, em_c[t, 0], em_e[t, 0], em_c[t, 1], em_e[t, 1]))

    # mixing eigenvectors at EVEC_T (unit-diagonal-normalized so the angle is meaningful)
    d = np.sqrt(np.abs([Cc[0, 0, EVEC_T], Cc[1, 1, EVEC_T]]))
    Cn = Cc.copy()
    for a in range(2):
        for b in range(2):
            Cn[a, b] = Cc[a, b] / (d[a] * d[b])
    em_n, lam_n, evec_n = gevp_states(Cn, EVEC_T, DT, AT)
    print("\n# mixing at t=%.2f (normalized ops):  light evec=[F2 %.3f, O2m %.3f]  heavy evec=[F2 %.3f, O2m %.3f]"
          % (EVEC_T * AT, evec_n[0, 0], evec_n[1, 0], evec_n[0, 1], evec_n[1, 1]))
    ang = np.degrees(np.arctan2(evec_n[1, 0], evec_n[0, 0]))
    print("# light-state mixing angle (F2 vs O2m) = %.1f deg  |  normalized off-diag C_Fs/sqrt(C_FF C_ss) @t=%.2f = %.3f"
          % (ang, EVEC_T * AT, Cc[0, 1, EVEC_T] / (d[0] * d[1])))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(TMAX) * AT
    fig, ax = plt.subplots(figsize=(9, 5.6))
    for line, lab, c in [(2.88, r"$F^2$ glueball 2.88", "gray"),
                         (1.6, r"$m_{PS}$ 1.6", "tab:orange"),
                         (2.15, r"$\sigma^2$ ground 2.15", "tab:green"),
                         (3.2, r"$2m_{PS}$ 3.2", "tab:brown")]:
        ax.axhline(line, ls="--", lw=0.9, color=c, alpha=0.6)
        ax.text(ts[-1] * 0.72, line + 0.03, lab, fontsize=8, color=c)
    for n, (col, mk, lab) in enumerate([("tab:blue", "o", "light state"), ("tab:red", "s", "heavy state")]):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 1.5)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=col, marker=mk, ms=5, lw=1.1, capsize=2.5, label=lab)
    ax.set_ylim(0.5, 4.0)
    ax.set_xlabel(r"$t$ (physical, $=t_{\rm lat}a_t$)")
    ax.set_ylabel(r"$\Delta_{\rm eff}=-\log\lambda/(dt\,a_t)$")
    ax.set_title(r"Coupled $\{F^2,O_{2m}\}$ GEVP  %s L1 %d cfg (bin %d)" % (tag, ncfg, BINSIZE), fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_o2m_gevp_v2_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
