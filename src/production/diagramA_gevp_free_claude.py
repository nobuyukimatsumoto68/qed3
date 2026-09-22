#!/usr/bin/env python3
# diagramA_gevp_free_claude.py
#   3x3 {sigma^2_00, O_2m, O_1m} GEVP built from DIAGRAM A ONLY (class-A bridging perms), S-leg
#   contact-subtracted, single free config.  Isolates the single-meson / (2,2) sector variationally.
#   Plain GEVP at T0 (full generalized eigenvalues; no Hankel) unless NOHANKEL=0 with a rebase.
#   Run: ENS=free LREF=2 NVDIR=distill_Nv84 LTAG=L2 M_PS=0.393 STATE_A=0.69 TWO_MESON=0.786 T0=1 python3 diagramA_gevp_free_claude.py

import os
os.environ.setdefault("ENS", "free")
os.environ.setdefault("LREF", "2")
os.environ.setdefault("NVDIR", "distill_Nv84")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import fs_gevp_point_perclass_claude as pc
import hankel_rebase_scan_claude as hs

DTMAX = int(os.environ.get("DTMAX", "24"))
LTAG = os.environ.get("LTAG", "L2")
T0 = int(os.environ.get("T0", "1"))
SPLIT = int(os.environ.get("SPLIT", "1"))
M_PS = float(os.environ.get("M_PS", "0.393"))
STATE_A = float(os.environ.get("STATE_A", "0.69"))
TWO_MESON = float(os.environ.get("TWO_MESON", "0.786"))
REB = int(os.environ.get("REB", "0"))                    # 1 = rebase (staged_project) instead of plain GEVP
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "1"))
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0").split(",")]
OFF_GEOM = [(0, 0), (0, 0), (0, SPLIT)]           # sigma^2_00, O_2m, O_1m time offsets (as in the 9-op builder)


def build_diagA_matrix():
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(dc.KS[0])
    dual = dc.dual_areas_from_mesh().astype(float)
    wY = dual * dc.Y00
    Pmap = G.antipodal_map()
    cls_of = [pc.classify(cy) for cy in G.PERMS]
    aperms = [ip for ip in range(len(G.PERMS)) if cls_of[ip] == "A"]
    C = np.full((3, 3, DTMAX), np.nan)
    omax = max(max(o) for o in OFF_GEOM)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt + omax < twin and s + omax < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        offs = set()
        for ga in range(3):
            for gb in range(3):
                vt = [dt + OFF_GEOM[ga][0], dt + OFF_GEOM[ga][1], OFF_GEOM[gb][0], OFF_GEOM[gb][1]]
                offs.update((vt[va], vt[vb]) for va in range(4) for vb in range(4))
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        for ga in range(3):
            for gb in range(3):
                vt = [dt + OFF_GEOM[ga][0], dt + OFF_GEOM[ga][1], OFF_GEOM[gb][0], OFF_GEOM[gb][1]]
                vspec = G.op_vspec(ga, ('i', 'j'), dual, wY) + G.op_vspec(gb, ('k', 'l'), dual, wY)
                acc = 0.0
                for ip in aperms:
                    acc += G.perm_contrib_folded(G.PERMS[ip], vt, bAS, vspec, Pmap).real
                C[ga, gb, dt] = acc / len(s0s)
    return C


def plain_gevp(Cmat):
    C0 = 0.5 * (Cmat[:, :, T0] + Cmat[:, :, T0].T)
    nop = Cmat.shape[0]
    lam = np.full((DTMAX, nop), np.nan)
    for dt in range(DTMAX):
        Ct = 0.5 * (Cmat[:, :, dt] + Cmat[:, :, dt].T)
        try:
            lam[dt] = np.sort(np.linalg.eigvals(np.linalg.solve(C0, Ct)).real)[::-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        return np.log(lam[:-1] / lam[1:])


def hankel_reb(Cmat):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0)
    return hs.rebased_effmass_fixed(Big, V, T0)


def main():
    C = build_diagA_matrix()
    C = 0.5 * (C + np.swapaxes(C, 0, 1))
    if REB:
        em = hankel_reb(C)
        method = "rebase Dt=%s reb%d@%d T0=%d" % (OFFSETS, NKEEP, REBT, T0)
    else:
        em = plain_gevp(C)
        method = "plain T0=%d" % T0
    tmax, nstates = em.shape

    print("# FREE %s diagram-A-ONLY 3x3 {sigma^2_00,O_2m,O_1m} GEVP (%s)" % (LTAG, method))
    print("# refs: m_PS=%.3f  (diagram A)=%.3f  2m_PS=%.3f" % (M_PS, STATE_A, TWO_MESON))
    print("#  t |   m0        m1        m2")
    for t in range(T0 + 1, min(tmax, 20)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f" % em[t, n] if np.isfinite(em[t, n]) else "  ---  " for n in range(nstates))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    for y, lab, col in [(M_PS, r"$m_{PS}=%.3f$" % M_PS, "tab:gray"),
                        (STATE_A, r"(diagram A)$=%.3f$" % STATE_A, "tab:orange"),
                        (TWO_MESON, r"two-meson$=%.3f$" % TWO_MESON, "tab:blue")]:
        ax.axhline(y, color=col, ls="--", lw=1, alpha=0.6)
        ax.text(tmax * 0.55, y + 0.008, lab, fontsize=9, color=col)
    cols = ["tab:red", "tab:purple", "tab:brown"]
    mkr = ["o", "s", "^"]
    for n in range(nstates):
        g = np.isfinite(em[:, n])
        ax.plot(ts[g], em[g, n], color=cols[n % 3], marker=mkr[n % 3], ms=5, lw=1.1, label="state %d" % n)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(tmax, 18))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE %s diagram-A-only 3x3 $\{\sigma^2_{00},O_{2m},O_{1m}\}$ GEVP (T0=%d)" % (LTAG, T0))
    ax.legend(fontsize=9, loc="lower left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/diagramA_gevp_free_%s_T0%d_claude.png" % (LTAG, T0)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
