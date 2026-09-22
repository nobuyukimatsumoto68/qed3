#!/usr/bin/env python3
# sigma2_3x3_free_gevp_plot_claude.py
#   FREE-limit 3x3 geometry GEVP {sigma^2_00, O_2m, O_1m} (PS block = FS by GW), block-Hankel, with the
#   CORRECT FREE reference lines.  Variant of sigma2_Peven_6x6_gevp_claude.py whose hardcoded 0.644/0.46/0.92
#   are the INTERACTING values and are misleading in the free case.  Free references:
#     m_PS = m_sigma = 2E0 = 0.378   (one-meson, Delta=2 ; contact-subtracted out of clean sigma^2)
#     {1,1,1,1} = (2,2) = 2E1 = 0.556 (Delta=4 single-meson radial excitation)
#     two-meson = 2 m_PS = 0.756      (Delta=4 threshold)
#   Reads the production free cache; PS block = ops 0,1,2.
#   Run: OFFSETS=0,4 REBT=4 NKEEP=2 T0=3 python3 sigma2_3x3_free_gevp_plot_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import hankel_rebase_scan_claude as hs

T0 = int(os.environ.get("T0", "3"))
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "2"))
NOHANKEL = int(os.environ.get("NOHANKEL", "0"))          # 1 = plain 3x3 GEVP at T0, no block-Hankel
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0,4").split(",")]
OPS = [int(x) for x in os.environ.get("OPS", "0,1,2").split(",")]   # PS geometry: 0=sigma^2_00 1=O_2m 2=O_1m
LTAG = os.environ.get("LTAG", "L1")              # "L1" or "L2" (for the cache name + output name)
CACHE = os.environ.get("CACHE",
    "sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_free_%s_1cfg_d1_claude.npy"
    % ("" if LTAG == "L1" else "L2_"))
CACHE = CACHE.replace("_free__", "_free_")       # L1 uses the bare 'free' name

M_PS = float(os.environ.get("M_PS", "0.378"))            # 2E0: L1 0.378, L2 0.393
STATE_1111 = float(os.environ.get("STATE_A", "0.556"))   # (diagram A)=2E1: L1 0.556, L2 ~0.69
TWO_MESON = float(os.environ.get("TWO_MESON", "0.756"))  # 2 m_PS: L1 0.756, L2 0.786

# basis label derived from OPS (cache layout: flavor*3+geom, flavor 0=PP 1=FF 2=FP)
NOP = len(OPS)
FLAV = ["PP", "FF", "FP"]
GEO_TX = ["sigma^2_00", "O_2m", "O_1m"]
GEO_LX = [r"\sigma^2_{00}", "O_{2m}", "O_{1m}"]
_fl = sorted(set(o // 3 for o in OPS))
_ge = sorted(set(o % 3 for o in OPS))
_factor = (OPS == [f * 3 + g for f in _fl for g in _ge])
if _factor and len(_fl) > 1:
    BASIS_TX = "{%s}x{%s}" % (",".join(FLAV[f] for f in _fl), ",".join(GEO_TX[g] for g in _ge))
    BASIS_LX = r"\{%s\}\times\{%s\}" % (",".join(FLAV[f] for f in _fl), ",".join(GEO_LX[g] for g in _ge))
elif _factor:
    BASIS_TX = "{%s}" % ",".join(GEO_TX[g] for g in _ge)
    BASIS_LX = r"\{%s\}" % ",".join(GEO_LX[g] for g in _ge)
else:
    BASIS_TX = "{%s}" % ",".join("%s-%s" % (FLAV[o // 3], GEO_TX[o % 3]) for o in OPS)
    BASIS_LX = r"\{%s\}" % ",".join(r"%s\text{-}%s" % (FLAV[o // 3], GEO_LX[o % 3]) for o in OPS)


def hankel_reb(Cmat):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0)
    return hs.rebased_effmass_fixed(Big, V, T0)


def plain_gevp(Cmat):
    # no Hankel: full generalized eigenvalues eig(C0^{-1} C(t)) at metric point T0 (C0 may be indefinite)
    DT = Cmat.shape[-1]
    C0 = 0.5 * (Cmat[:, :, T0] + Cmat[:, :, T0].T)
    nop = Cmat.shape[0]
    lam = np.full((DT, nop), np.nan)
    for dt in range(DT):
        Ct = 0.5 * (Cmat[:, :, dt] + Cmat[:, :, dt].T)
        try:
            lam[dt] = np.sort(np.linalg.eigvals(np.linalg.solve(C0, Ct)).real)[::-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        return np.log(lam[:-1] / lam[1:])


def main():
    allC = np.load(CACHE)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    allC = allC[:, OPS][:, :, OPS]
    Cmat = allC.mean(0)                            # single free config
    if NOHANKEL:
        em = plain_gevp(Cmat)
        method = "plain GEVP (NO Hankel) T0=%d" % T0
    else:
        em = hankel_reb(Cmat)
        method = "Hankel Dt=%s reb%d@%d T0=%d" % (OFFSETS, NKEEP, REBT, T0)
    tmax, nstates = em.shape

    print("# FREE %s %dx%d %s  %s" % (LTAG, NOP, NOP, BASIS_TX, method))
    print("# refs: m_PS=%.3f  (diagram A)=%.3f  2m_PS=%.3f" % (M_PS, STATE_1111, TWO_MESON))
    print("#  t | " + "  ".join("m%d" % n for n in range(nstates)))
    for t in range(T0, min(tmax, 18)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f" % em[t, n] if np.isfinite(em[t, n]) else "  ---  " for n in range(nstates))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    for y, lab, col in [(M_PS, r"$m_{PS}=m_\sigma=0.378$ (one-meson, $\Delta=2$)", "tab:green"),
                        (STATE_1111, r"(diagram A) $=0.556$", "tab:orange"),
                        (TWO_MESON, r"$2m_{PS}=0.756$ (two-meson)", "tab:blue")]:
        ax.axhline(y, color=col, ls="--", lw=1, alpha=0.6)
        ax.text(tmax * 0.50, y + 0.008, lab, fontsize=9, color=col)
    cols = ["tab:red", "tab:purple", "tab:brown"]
    mkr = ["o", "s", "^"]
    for n in range(nstates):
        g = np.isfinite(em[:, n])
        ax.plot(ts[g], em[g, n], color=cols[n % 3], marker=mkr[n % 3], ms=5, lw=1.1, label="state %d" % n)
    if not NOHANKEL:
        ax.axvline(REBT, color="k", ls=":", lw=0.8, alpha=0.3)
    ax.set_ylim(0.3, 1.0)
    ax.set_xlim(T0, min(tmax, 16))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    mtag = "no Hankel T0=%d" % T0 if NOHANKEL else "Hankel Dt=%s reb%d@%d T0=%d" % (OFFSETS, NKEEP, REBT, T0)
    ax.set_title(r"FREE %s %dx%d $%s$ GEVP  (%s)" % (LTAG, NOP, NOP, BASIS_LX, mtag))
    ax.legend(fontsize=9, loc="lower left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    ftag = "noHankel_T0%d" % T0 if NOHANKEL else "off%s_reb%d_T0%d" % ("".join(str(o) for o in OFFSETS), NKEEP, T0)
    out = "figs/sigma2_%dx%d_free_%s_gevp_%s_claude.png" % (NOP, NOP, LTAG, ftag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
