#!/usr/bin/env python3
# f2_dps_id_gevp_v2_claude.py
#   EXPERIMENT: add the identity operator to the {D'_S, F^2} GEVP -> {1, D'_S, F^2} 3x3.
#   The identity absorbs the vacuum variationally, so we use RAW (un-subtracted) correlators:
#     C_11 = <1 1> = 1 ;  C_1j(t) = <O_j> (constant one-point) ;  C_ij = <O_i(t) O_j(0)> raw.
#   op0 = 1 (identity) , op1 = D'_S = Tr[Phi tilde_tau Phi tilde_tau] , op2 = F^2 shape op.
#   Same window-triangular sampling (Cauchy-Schwarz-consistent), operators normalized at TREF.
#   Fixed-t0 GEVP; effmass = log(lam(t)/lam(t+1)) = LATTICE a_t m.  Vacuum -> level 0 (Delta~0).
#   Run: ENS=... NVDIR=distill_Nv24 BINSIZE=10 T0=3 python3 f2_dps_id_gevp_v2_claude.py

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

Nt = 128
OP_F = int(os.environ.get("OP_F", "0"))
AT = float(os.environ.get("AT", "0.2"))
CONTACT = float(os.environ.get("CONTACT", "0.5"))
TMAX = int(os.environ.get("TMAX", "24"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
T0 = int(os.environ.get("T0", "3"))
TREF = int(os.environ.get("TREF", "2"))
RTOL = float(os.environ.get("RTOL", "1e-9"))


def gluedir():
    return "data_" + dc.ENS


def load_of(k):
    with h5py.File("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k), "r") as g:
        return np.array(g["O"])[OP_F].astype(float)


def dps_series(k, w00):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Iv = np.eye(tau.shape[-1])
    D = np.empty(twin)
    for a in range(twin):
        Phi = (V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T))
        tt = tau[a, a] - CONTACT * Iv
        D[a] = np.trace(Phi @ tt @ Phi @ tt).real
    return D, tsrc0, twin


def inv_sqrt_sym(M, rtol):
    M = 0.5 * (M + M.T)
    w, Vv = np.linalg.eigh(M)
    inv = np.where((w > rtol * w.max()) & (w > 0), 1.0 / np.sqrt(np.abs(w)), 0.0)
    return (Vv * inv) @ Vv.T


def gevp_lat(C, t, t0):
    si = inv_sqrt_sym(0.5 * (C[:, :, t0] + C[:, :, t0].T), RTOL)
    M = si @ (0.5 * (C[:, :, t] + C[:, :, t].T)) @ si
    return np.sort(np.linalg.eigvalsh(0.5 * (M + M.T)))[::-1]


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k))]
    ncfg = len(ks)
    print("# ENS=%s cfg=%d  {1, D'_S, F^2} GEVP (identity, RAW corr)  T0=%d BIN=%d" % (tag, ncfg, T0, BINSIZE))

    cff = np.zeros((ncfg, TMAX + 1))
    cfd = np.zeros((ncfg, TMAX + 1))
    cdd = np.zeros((ncfg, TMAX + 1))
    ofm = np.zeros(ncfg)
    dm = np.zeros(ncfg)
    for ic, k in enumerate(ks):
        of = load_of(k)
        D, tsrc0, twin = dps_series(k, w00)
        ofw = of[(tsrc0 + np.arange(twin)) % Nt]
        ofm[ic] = ofw.mean()
        dm[ic] = D.mean()
        for t in range(TMAX + 1):
            if t < twin:
                a = np.arange(twin - t)
                cff[ic, t] = np.mean(ofw[a + t] * ofw[a])
                cdd[ic, t] = np.mean(D[a + t] * D[a])
                cfd[ic, t] = 0.5 * (np.mean(ofw[a + t] * D[a]) + np.mean(D[a + t] * ofw[a]))
            else:
                cff[ic, t] = cdd[ic, t] = cfd[ic, t] = np.nan

    nb = ncfg // BINSIZE
    def binned(a):
        return np.array([a[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    bff, bfd, bdd = binned(cff), binned(cfd), binned(cdd)
    bof, bdm = binned(ofm), binned(dm)

    def build_C(sel):
        # RAW 3x3 with identity op0 ; C_1j = one-point (constant in t)
        ff = bff[sel].mean(0)
        fd = bfd[sel].mean(0)
        dd = bdd[sel].mean(0)
        of1 = bof[sel].mean()
        dm1 = bdm[sel].mean()
        C = np.zeros((3, 3, TMAX + 1))
        C[0, 0] = 1.0
        C[0, 1] = C[1, 0] = dm1
        C[0, 2] = C[2, 0] = of1
        C[1, 1] = dd
        C[2, 2] = ff
        C[1, 2] = C[2, 1] = fd
        return C

    allsel = np.arange(nb)
    Craw = build_C(allsel)
    dn = np.sqrt(np.abs([Craw[0, 0, TREF], Craw[1, 1, TREF], Craw[2, 2, TREF]]))
    NRM = np.outer(dn, dn)[:, :, None]

    def C_norm(sel):
        return build_C(sel) / NRM

    def effmass(sel):
        lam = np.full((TMAX + 1, 3), np.nan)
        C = C_norm(sel)
        for t in range(1, TMAX + 1):
            if np.all(np.isfinite(C[:, :, t])) and np.all(np.isfinite(C[:, :, T0])):
                try:
                    e = gevp_lat(C, t, T0)
                    lam[t, :len(e)] = e[:3]
                except Exception:
                    pass
        with np.errstate(all="ignore"):
            r = lam[:-1] / lam[1:]
            r[r <= 0] = np.nan
            return np.log(r)

    em_c = effmass(allsel)
    ems = np.array([effmass(np.delete(allsel, j)) for j in range(nb)])
    em_e = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, 0))

    print("\n#  t | lvl0 (vac?)        lvl1               lvl2         [sig2 gd 0.46, F2 0.576, 2mPS 0.644]")
    for t in range(1, TMAX - 1):
        if np.any(np.isfinite(em_c[t])):
            print("#  %2d | %7.4f(%.4f)   %7.4f(%.4f)   %7.4f(%.4f)"
                  % (t, em_c[t, 0], em_e[t, 0], em_c[t, 1], em_e[t, 1], em_c[t, 2], em_e[t, 2]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(TMAX)
    fig, ax = plt.subplots(figsize=(9, 5.8))
    ax.axhline(0.644, color="gray", ls="--", lw=1, alpha=0.6); ax.text(15, 0.66, r"$2m_{PS}=0.644$", fontsize=8, color="gray")
    ax.axhline(0.576, color="black", ls="-.", lw=1, alpha=0.6); ax.text(15, 0.585, r"$F^2=0.576$", fontsize=8, color="black")
    ax.axhline(0.46, color="tab:green", ls=":", lw=1, alpha=0.6); ax.text(15, 0.475, r"$\sigma^2$ gd 0.46", fontsize=8, color="tab:green")
    for n, (c, m, lab) in enumerate([("tab:gray", "x", "level 0 (vac)"), ("tab:blue", "o", "level 1"), ("tab:red", "s", "level 2")]):
        g = np.isfinite(em_c[:TMAX, n]) & np.isfinite(em_e[:TMAX, n]) & (em_e[:TMAX, n] < 0.4)
        ax.errorbar(ts[g], em_c[:TMAX][g, n], yerr=em_e[:TMAX][g, n], color=c, marker=m, ms=5, lw=1.1, capsize=2.5, label=lab)
    ax.set_ylim(-0.1, 0.95)
    ax.set_xlim(0.5, TMAX - 2)
    ax.set_xlabel(r"$t$ (timeslice)")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$ (lattice)")
    ax.set_title(r"$\{1, D'_S, F^2\}$ GEVP (identity, raw)  T0=%d  %s L1 %d cfg (bin %d)" % (T0, tag, ncfg, BINSIZE), fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_dps_id_gevp_v2_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
