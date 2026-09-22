# Hankel/GPOF tuning for the PROJECTED magnetic operator (C00 = VSH-projected Psi_ell=1).
# Loads the saved binned magnetic correlator MA_binned (from interacting_vsh_L1) -- NO recompute.
# Exposes the block-Hankel params via env so they are one-line tweaks:
#   OFF (offsets, comma list), REBT, NKEEP, T0H, WLO, WHI.  Prints plateau + writes effmass plot.
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
import numpy as np
sys.path.insert(0, "final/analysis_axial")
import effmass_axial_tp_l3_perm_hankel_claude as H

AT = 0.2
NPZ = "final/analysis_axial/interacting_vsh_axial_L1_nf2g1_claude.npz"
OFF = [int(x) for x in os.environ.get("OFF", "0,3,6").split(",")]
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "1"))
T0H = int(os.environ.get("T0H", "1"))
WLO = int(os.environ.get("WLO", "12"))
WHI = int(os.environ.get("WHI", "15"))
NLEV = int(os.environ.get("NLEV", "1"))     # how many GPOF levels to plot (<= NKEEP)


def emcurve(C, lev):
    sgn = np.sign(C[3]) if abs(C[3]) > 0 else 1.0
    em, _ = H.hankel_effmass_scalar((sgn * C).astype(float), OFF, REBT, NKEEP, T0H, AT)
    return em[:, lev]


def main():
    d = np.load(NPZ)
    binned = d["MA_binned"]         # (nb, dmax)
    nb = binned.shape[0]
    dmax = binned.shape[1]
    cen = binned.mean(0)
    jk = np.array([np.delete(binned, b, 0).mean(0) for b in range(nb)])
    print("# projected magnetic Hankel: OFF=%s REBT=%d NKEEP=%d T0H=%d win[%d,%d]  (%d bins)"
          % (OFF, REBT, NKEEP, T0H, WLO, WHI, nb), flush=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = ["tab:red", "tab:orange", "tab:brown"]
    mks = ["o", "s", "^"]
    plt.figure(figsize=(7.6, 5.0))
    plat = {}
    for lev in range(NLEV):
        emc = emcurve(cen, lev)
        emj = np.array([emcurve(jk[b], lev) for b in range(nb)])
        eme = np.sqrt((nb - 1) * np.nanmean((emj - np.nanmean(emj, 0)) ** 2, 0))
        # WLO,WHI are dt values; effmass array index i sits at dt=i+1, so the dt window maps to indices [WLO-1,WHI-1]
        seg = emc[WLO - 1:WHI]
        segj = emj[:, WLO - 1:WHI]
        pv = np.nanmean(seg)
        pj = np.nanmean(segj, 1)
        pe = np.sqrt((nb - 1) * np.nanmean((pj - np.nanmean(pj)) ** 2))
        plat[lev] = (pv, pe)
        dd = np.arange(1, len(emc) + 1)
        g = np.isfinite(emc) & np.isfinite(eme) & (eme < 0.2)
        plt.errorbar(dd[g], emc[g], yerr=eme[g], marker=mks[lev], ms=4, capsize=2, lw=0.9,
                     color=cols[lev], label="GPOF level %d = %.4f(%d)" % (lev, pv, round(pe * 1e4)))
        # fitted plateau: central line + \pm 1 sigma error band across the fit window (dt WLO..WHI)
        plt.hlines(pv, WLO, WHI, color=cols[lev], ls="-", lw=1.6, zorder=5)
        plt.fill_between([WLO, WHI], pv - pe, pv + pe, color=cols[lev], alpha=0.25, lw=0, zorder=4)
        print("#   level %d plateau a_t*m = %.4f +- %.4f  (dt window [%d,%d])" % (lev, pv, pe, WLO, WHI), flush=True)
    plt.axvspan(WLO, WHI, color="gray", alpha=0.12)
    plt.axhline(0.3606, color="silver", ls=":", lw=0.9)
    plt.text(2, 0.3606 + 0.002, "prior magnetic 0.3606", fontsize=7, color="gray")
    plt.ylim(0.30, 0.50)
    plt.xlim(1, min(dmax, 24))
    plt.xlabel("dt")
    plt.ylabel("GPOF effmass = $a_t m$")
    plt.title("Projected magnetic $\\Psi_{\\ell=1}$ Hankel  OFF=%s REBT=%d NKEEP=%d T0=%d  L1 Nf2 g1.0"
              % (OFF, REBT, NKEEP, T0H))
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    out = "final/analysis_axial/hankel_mag_proj_claude.png"
    plt.savefig(out, dpi=150)
    print("# wrote %s" % out, flush=True)


if __name__ == "__main__":
    main()
